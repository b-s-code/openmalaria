/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "PkPd/Drug/LSTMDrugThreeComp.h"
#include "Host/WithinHost/Infection/CommonInfection.h"
#include "util/errors.h"
#include "util/StreamValidator.h"

#include <gsl/gsl_integration.h>
#include <limits>

using namespace std;

namespace OM { namespace PkPd {

LSTMDrugThreeComp::LSTMDrugThreeComp(const LSTMDrugType& type, LocalRng& rng) :
    LSTMDrug(type.sampleVd(rng)),
    m_typeData(type),
    m_concA(0.0), m_concB(0.0), m_concC(0.0), m_concAbc(0.0),
    m_elimSample(type.sampleElimRate(rng)),
    m_k12(type.sampleK12(rng)),
    m_k21(type.sampleK21(rng)),
    m_k13(type.sampleK13(rng)),
    m_k31(type.sampleK31(rng)),
    m_nka(-type.sampleKa(rng)),
    m_lastBm(numeric_limits<double>::quiet_NaN()),
    m_na(numeric_limits<double>::quiet_NaN()),
    m_nb(numeric_limits<double>::quiet_NaN()),
    m_ng(numeric_limits<double>::quiet_NaN()),
    A(numeric_limits<double>::quiet_NaN()),
    B(numeric_limits<double>::quiet_NaN()),
    C(numeric_limits<double>::quiet_NaN())
{
    // These are from the Monolix article, pp38-39.
}

size_t LSTMDrugThreeComp::getIndex() const {
    return m_typeData.getIndex();
}
double LSTMDrugThreeComp::getConcentration(size_t index) const {
    if( index == m_typeData.getIndex() ) return conc();
    else return 0.0;
}

void LSTMDrugThreeComp::updateCached(double bodyMass) const{
    if( m_lastBm == bodyMass ) return;
    
    double k = m_elimSample * pow(bodyMass, m_typeData.negMExponent());
    
    const double a0 = k*m_k21*m_k31;
    const double a1 = k*m_k31 + m_k21*m_k31 + m_k21*m_k13 + k*m_k21 + m_k31*m_k12;
    const double a2 = k + m_k12 + m_k13 + m_k21 + m_k31;
    
    const double third = 1.0 / 3.0;
    const double p = a1 - third * a2 * a2;
    const double q = (2.0 / 27.0) * a2 * a2 * a2 - third * a1 * a2 + a0;
    
    const double at = third * a2;
    const double rt = -third * p;
    const double r2 = 2.0 * sqrt(rt);
    
    const double phi = acos(-q / (rt * r2)) * third;
    const double pi23 = (2.0 / 3.0) * PI;
    
    // negative alpha, beta, gamma:
    m_na = r2 * cos(phi) - at;
    m_nb = r2 * cos(phi + pi23) - at;
    m_ng = r2 * cos(phi + 2.0 * pi23) -at;
    
    // A, B, C from Monolix 1.3.3 (p44):
    const double kaOv = -m_nka / (m_volDist * bodyMass);
    A = kaOv * (m_k21 + m_na) * (m_k31 + m_na) / ((m_na - m_nka) * (m_nb - m_na) * (m_ng - m_na));
    B = kaOv * (m_k21 + m_nb) * (m_k31 + m_nb) / ((m_nb - m_nka) * (m_na - m_nb) * (m_ng - m_nb));
    C = kaOv * (m_k21 + m_ng) * (m_k31 + m_ng) / ((m_ng - m_nka) * (m_nb - m_ng) * (m_na - m_ng));
    
    m_lastBm = bodyMass;
}

/// Parameters for func_fC
struct ParamsFC {
    double cA, cB, cC, cABC;    // concentration parameters
    double na, nb, ng, nka;     // decay parameters
    double n;       // slope: unitless
    double V;       // max killing rate: unitless
    double Kn;      // IC50^n: (mg/kg) ^ n
};
/** Function for calculating concentration and then killing function at time t
 * 
 * @param t The variable being integrated over (in this case, time since start
 *      of day or last dose, units days)
 * @param pp Pointer to a Params_fC struct
 * @return killing rate (unitless)
 */
double funcFC( double t, void* pp ){
    const ParamsFC& p = *static_cast<const ParamsFC*>( pp );
    
    // exponential decay of drug concentration:
    const double concA = p.cA * exp(p.na * t);
    const double concB = p.cB * exp(p.nb * t);
    const double concC = p.cC * exp(p.ng * t);
    const double concABC = p.cABC * exp(p.nka * t);
    const double conc = concA + concB + concC - concABC;      // mg/l
    
    const double cn = pow(conc, p.n);        // (mg/l) ^ n
    const double fC = p.V * cn / (cn + p.Kn);       // unitless
    return fC;
}
const size_t GslIntgMaxIter = 1000;     // 10 seems enough, but no harm in using a higher value
gsl_integration_workspace *gslIntgrWksp = gsl_integration_workspace_alloc (GslIntgMaxIter);
//NOTE: we "should" free, but mem-leaks at end of program aren't really important
// gsl_integration_workspace_free (gsl_intgr_wksp);

double LSTMDrugThreeComp::calculateFactor(const ParamsFC& p, double duration) const{
    gsl_function F;
    F.function = &funcFC;
    // gsl_function doesn't accept const; we re-apply const later
    F.params = static_cast<void*>(const_cast<ParamsFC*>(&p));
    
    // NOTE: tolerances are arbitrary, but seem to be sufficient
    const double absEps = 1e-4, relEps = 1e-4;
    // NOTE: 1 through 6 are different algorithms of increasing complexity
    const int qagRule = 1;     // alg 1 seems to be good enough
    double intfC, errEps;
    
    int r = gsl_integration_qag (&F, 0.0, duration, absEps, relEps,
                                 GslIntgMaxIter, qagRule, gslIntgrWksp, &intfC, &errEps);
    if( r != 0 ){
        throw TRACED_EXCEPTION( "calculateFactor: error from gsl_integration_qag",util::Error::GSL );
    }
    if( errEps > 5e-2 ){
        // This could be a warning, except that warnings tend to be ignored.
        ostringstream msg;
        msg << "calculateFactor: error epsilon is large: "<<errEps<<" (integral is "<<intfC<<")";
        throw TRACED_EXCEPTION( msg.str(), util::Error::GSL );
    }
    
    return 1.0 / exp( intfC );  // drug factor
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrugThreeComp::calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double bodyMass) const {
    if( conc() == 0.0 && m_doses.size() == 0 ) return 1.0; // nothing to do
    updateCached(bodyMass);
    
    ParamsFC p;
    p.cA = m_concA;       p.cB = m_concB;   p.cC = m_concC;   p.cABC = m_concAbc;
    p.na = m_na;  p.nb = m_nb;      p.ng = m_ng;      p.nka = m_nka;
    const LSTMDrugPD& pd = m_typeData.getPD(inf->genotype());
    p.n = pd.slope();   p.V = pd.maxKillingRate();
    p.Kn = pd.IC50_pow_slope(rng, m_typeData.getIndex(), inf);
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    for( const TimeConc& timeConc : m_doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( timeConc.first < 1.0 /*i.e. today*/ ){
            if( time < timeConc.first ){
                double duration = timeConc.first - time;
                totalFactor *= calculateFactor(p, duration);
                p.cA *= exp(p.na * duration);
                p.cB *= exp(p.nb * duration);
                p.cC *= exp(p.ng * duration);
                p.cABC *= exp(p.nka * duration);
                time = timeConc.first;
            }else{ assert( time == timeConc.first ); }
            // add dose:
            p.cA += A * timeConc.second;
            p.cB += B * timeConc.second;
            p.cC += C * timeConc.second;
            p.cABC += (A + B + C) * timeConc.second;
        }else /*i.e. tomorrow or later*/{
            // ignore
        }
    }
    if( time < 1.0 ){
        totalFactor *= calculateFactor(p, 1.0 - time);
    }
    
    return totalFactor;
}

void LSTMDrugThreeComp::updateConcentration (double bodyMass) {
    if( conc() == 0.0 && m_doses.size() == 0 ) return;     // nothing to do
    updateCached(bodyMass);
    
    // exponential decay of existing quantities:
    //TODO(performance): is it faster to pre-calculate these and either store extra
    // parameters or adapt uses of alpha, beta, gamma, etc. below?
    m_concA *= exp(m_na);
    m_concB *= exp(m_nb);
    m_concC *= exp(m_ng);
    m_concAbc *= exp(m_nka);
    
    size_t dosesTaken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& timeConc : m_doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( timeConc.first < 1.0 /*i.e. today*/ ){
            // add dose:
            const double qty = timeConc.second;
            m_concA += A * qty * exp(m_na * (1.0 - timeConc.first));
            m_concB += B * qty * exp(m_nb * (1.0 - timeConc.first));
            m_concC += C * qty * exp(m_ng * (1.0 - timeConc.first));
            m_concAbc += (A + B + C) * qty * exp(m_nka * (1.0 - timeConc.first));
            dosesTaken += 1;
        }else /*i.e. tomorrow or later*/{
            timeConc.first -= 1.0;
        }
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    m_doses.erase(m_doses.begin(), m_doses.begin() + dosesTaken);
    
    util::streamValidate( conc() );
    if( conc() < m_typeData.getNegligibleConcentration() ){
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        m_concA = 0.0;
        m_concB = 0.0;
        m_concC = 0.0;
        m_concAbc = 0.0;
    }
}

void LSTMDrugThreeComp::checkpoint(ostream& stream){
    m_concA & stream;
    m_concB & stream;
    m_concC & stream;
    m_concAbc & stream;
    m_elimSample & stream;
    m_k12 & stream;
    m_k21 & stream;
    m_k13 & stream;
    m_k31 & stream;
    m_nka & stream;
    m_lastBm & stream;
    m_na & stream;
    m_nb & stream;
    m_ng & stream;
    A & stream;
    B & stream;
    C & stream;
}
void LSTMDrugThreeComp::checkpoint(istream& stream){
    m_concA & stream;
    m_concB & stream;
    m_concC & stream;
    m_concAbc & stream;
    m_elimSample & stream;
    m_k12 & stream;
    m_k21 & stream;
    m_k13 & stream;
    m_k31 & stream;
    m_nka & stream;
    m_lastBm & stream;
    m_na & stream;
    m_nb & stream;
    m_ng & stream;
    A & stream;
    B & stream;
    C & stream;
}

}
}
