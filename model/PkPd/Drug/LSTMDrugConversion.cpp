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

#include "PkPd/Drug/LSTMDrugConversion.h"
#include "Host/WithinHost/Infection/CommonInfection.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <gsl/gsl_integration.h>
#include <limits>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrugConversion::LSTMDrugConversion(const LSTMDrugType& parent,
                                       const LSTMDrugType& metabolite,
                                       LocalRng& rng) :
    LSTMDrug(parent.sampleVd(rng)),
    m_parentType(parent),
    m_metaboliteType(metabolite),
    m_qtyG(0.0), m_qtyP(0.0), m_qtyM(0.0),
    m_nka(-parent.sampleKa(rng)),
    m_nkPSample(-parent.sampleElimRate(rng)),
    m_nconvSample(-parent.sampleConvRate(rng)),
    m_nkMSample(-metabolite.sampleElimRate(rng)),
    m_volDistMetabolite(metabolite.sampleVd(rng))
{}

//LSTMDrugConversion::~LSTMDrugConversion(){}

size_t LSTMDrugConversion::getIndex() const {
    return m_parentType.getIndex();       // parent drug's index should work
}
double LSTMDrugConversion::getConcentration(size_t index) const {
    if( index == m_parentType.getIndex() ){
        return getParentConcentration();
    }else if( index == m_metaboliteType.getIndex() ){
        return getMetaboliteConcentration();
    }else return 0.0;
}
double LSTMDrugConversion::getParentConcentration() const {
    return m_qtyP / (m_volDist * m_lastBm);
}
double LSTMDrugConversion::getMetaboliteConcentration() const {
    return m_qtyM / (m_volDistMetabolite * m_lastBm);
}


/// Parameters for func_convFactor and LSTMDrugConversion::calculateFactor
struct ParamsConvFactor {
    // Quantities of parent in gut, parent in circulation and metabolite in
    // circulation, in mg (A, B, C in paper):
    double qtyG, qtyP, qtyM;
    // nka = -x (negative of absorption rate)
    // nkM = -k
    // nl = -(y + z)
    double nka, nkM, nl;
    // terms involving x, y, z and the molecular weight ratio:
    double f, g, h, i, j;
    
    double invVdP, invVdM;     // 1.0 / (Vd * body mass) for parent and metabolite; units: 1/l
    double nP, nM;       // slope: unitless
    double VP, VM;       // max killing rate: unitless
    double KnP, KnM;      // IC50^n: (mg/kg) ^ n
};

double calculateParentQuantity( struct ParamsConvFactor p, double expAbsorb, double expPLoss ) {
    return p.f * p.qtyG * expAbsorb
        + (p.qtyP - p.f * p.qtyG) * expPLoss;
}

double calculateParentDrugFactor( struct ParamsConvFactor p, double expAbsorb, double expPLoss ) {
    const double qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
    const double cP = qtyP * p.invVdP;                  // concentrations; mg/l*/
    const double cnP = pow(cP, p.nP);                   // (mg/l) ^ n
    const double fCP = p.VP * cnP / (cnP + p.KnP);      // unitless
    return fCP;
}

double calculateMetaboliteQuantity(struct ParamsConvFactor p, double expAbsorb, double expPLoss, double t) {
    return p.g * p.qtyG * expAbsorb
        + (p.h * p.qtyG - p.i * p.qtyP) * expPLoss
        + (p.j * p.qtyG + p.i * p.qtyP + p.qtyM) * exp(p.nkM * t);
}

double calculateMetaboliteDrugFactor( struct ParamsConvFactor p, double expAbsorb, double expPLoss, double t ) {
    const double qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, t);
    const double cM = qtyM * p.invVdM;              // concentrations; mg/l
    const double cnM = pow(cM, p.nM);               // (mg/l) ^ n
    const double fCM = p.VM * cnM / (cnM + p.KnM);  // unitless
    return fCM;
}

// size_t intg_steps = 0;
/** Function for calculating concentration and then killing function at time t
 * 
 * @param t The variable being integrated over (in this case, time since start
 *      of day or last dose, units days)
 * @param pp Pointer to a Params_convFactor struct
 * @return killing rate (unitless)
 */
double funcConvFactor( double t, void* pp ){
//     intg_steps += 1;
    const ParamsConvFactor& p = *static_cast<const ParamsConvFactor*>( pp );
    const double expAbsorb = exp(p.nka * t), expPLoss = exp(p.nl * t);
    const double fCP = calculateParentDrugFactor( p, expAbsorb, expPLoss );
    const double fCM = calculateMetaboliteDrugFactor( p, expAbsorb, expPLoss, t );
    // use the most effective killing factor (from area under the drug kill curve), which is the one with the bigger number
    return max(fCP,fCM);
}

const size_t GslIntgConvMaxIter = 1000;     // 10 seems enough, but no harm in using a higher value
gsl_integration_workspace *gslIntgrConvWksp = gsl_integration_workspace_alloc (GslIntgConvMaxIter);
//NOTE: we "should" free, but mem-leaks at end of program aren't really important
// gsl_integration_workspace_free (gsl_intgr_conv_wksp);
double LSTMDrugConversion::calculateFactor(const ParamsConvFactor& p, double duration) const{
    gsl_function F;
    F.function = &funcConvFactor;
    // gsl_function doesn't accept const; we re-apply const later
    F.params = static_cast<void*>(const_cast<ParamsConvFactor*>(&p));
    
    // We use exp(-result), so small absolute differences can matter (but also
    // using smaller abs_eps is cheap). We likely don't need high rel precision.
    const double absEps = 1e-5, relEps = 1e-2;
    // NOTE: 1 through 6 are different algorithms of increasing complexity
    const int qagRule = 1;     // alg 1 seems to be good enough
    double intfC, errEps;      // intfC will carry our result; err_eps is a measure of accuracy of the result
    
//     intg_steps = 0;
    int r = gsl_integration_qag (&F, 0.0, duration, absEps, relEps,
                                 GslIntgConvMaxIter, qagRule, gslIntgrConvWksp, &intfC, &errEps);
    if( r != 0 ){
        throw TRACED_EXCEPTION( "calculateFactor: error from gsl_integration_qag",util::Error::GSL );
    }
    // Testing err_eps is redundant with GSL's built-in tests
//     cout << "integration steps: " << intg_steps << endl;
//     cout << "duration: " << duration << ", AUC: " << intfC << endl;
    return exp( -intfC );  // drug factor
}

void LSTMDrugConversion::setConversionParameters(ParamsConvFactor& p,
        double bodyMass) const
{
    p.qtyG = m_qtyG; p.qtyP = m_qtyP; p.qtyM = m_qtyM;
    
    // decay "constants" (dependent on body mass):
    const double nkP = m_nkPSample * pow(bodyMass, m_parentType.negMExponent());      // -y
    const double nconv = m_nconvSample * pow(bodyMass, m_parentType.negMExponent());  // -z
    p.nka = m_nka;        // -x
    p.nkM = m_nkMSample * pow(bodyMass, m_metaboliteType.negMExponent());  // -k
    p.nl = nkP + nconv;    // -(y + z)
    
    p.f = m_nka / (p.nl - m_nka);                           // -x / (x-y-z) = x / (y+z-x)
    // Let u = parentType.molecular_weight_ratio() = M_C / M_B:
    const double rz = m_parentType.molecularWeightRatio() * nconv; // -u*z
    p.g = rz * m_nka / ((m_nka - p.nl) * (m_nka - p.nkM));    // u*z*x / ((y+z-x) * (k-x))
    p.h = rz * m_nka / ((m_nka - p.nl) * (p.nkM - p.nl));   // u*z*x / ((y+z-x) * (y+z-k))
    p.i = rz / (p.nl - p.nkM);                          // -u*z / (k - y - z) = u*z / (y+z-k)
    p.j = rz * m_nka / ((p.nkM - p.nl) * (p.nkM - m_nka));  // u*z*x / ((y+z-k) * (x-k))
    
    p.invVdP = 1.0 / (m_volDist * bodyMass); p.invVdM = 1.0 / (m_volDistMetabolite * bodyMass);
}

void LSTMDrugConversion::setKillingParameters(LocalRng& rng, ParamsConvFactor& p,
        WithinHost::CommonInfection *inf) const
{
    uint32_t genotype = inf->genotype();
    const LSTMDrugPD& pdP = m_parentType.getPD(genotype), &pdM =
            m_metaboliteType.getPD(genotype);
    p.nP = pdP.slope();
    p.VP = pdP.maxKillingRate();
    p.nM = pdM.slope();
    p.VM = pdM.maxKillingRate();
    
    // Use custom code here because we need to handle covariance
    auto pIndex = m_parentType.getIndex();
    if (inf->Kn.count(pIndex) != 0) {
        // Read cached values: IC50 ^ n
        p.KnP = inf->Kn.at(pIndex);
        p.KnM = inf->Kn.at(m_metaboliteType.getIndex());
    } else {
        // First usage for this infection / treatment: sample, optionally with correlation.
        auto zscore = NormalSample::generate(rng);
        p.KnP = pdP.IC50_pow_slope(zscore);
        inf->Kn[pIndex] = p.KnP;
        
        auto metabZscore = m_parentType.iC50CorrelatedSample(zscore, rng);
        p.KnM = pdM.IC50_pow_slope(metabZscore);
        inf->Kn[m_metaboliteType.getIndex()] = p.KnM;
    }
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrugConversion::calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double bodyMass) const {
    if( m_qtyG == 0.0 && m_qtyP == 0.0 && m_qtyM == 0.0 && m_doses.size() == 0 ){
        return 1.0; // nothing to do
    }
    
    ParamsConvFactor p;
    setConversionParameters(p, bodyMass);
    setKillingParameters(rng, p, inf);
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    for( const TimeConc& timeConc : m_doses ){
        // we iterate through doses in time order (since doses are sorted)
        if( timeConc.first < 1.0 /*i.e. today*/ ){
            if( time < timeConc.first ){
                double duration = timeConc.first - time;
                totalFactor *= calculateFactor(p, duration);
                const double expAbsorb = exp(m_nka * duration), expPLoss = exp(p.nl * duration);
                p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
                p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
                p.qtyG *= expAbsorb;
                time = timeConc.first;
            }else{ assert( time == timeConc.first ); }
            // add to quantity of drug in gut:
            p.qtyG += timeConc.second;   // units: mg
        }else /*i.e. tomorrow or later*/{
            // ignore
        }
    }
    if( time < 1.0 ){
        totalFactor *= calculateFactor(p, 1.0 - time);
    }
    
    return totalFactor;
}

void LSTMDrugConversion::updateConcentration( double bodyMass ){
    if( m_qtyG == 0.0 && m_qtyP == 0.0 && m_qtyM == 0.0 && m_doses.size() == 0 ){
        return; // nothing to do
    }
    m_lastBm = bodyMass;
    
    ParamsConvFactor p;
    setConversionParameters(p, bodyMass);
    
    double time = 0.0, duration;
    size_t dosesTaken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& timeConc : m_doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( timeConc.first < 1.0 /*i.e. today*/ ){
            if( (duration = timeConc.first - time) > 0.0 ){
                const double expAbsorb = exp(m_nka * duration), expPLoss = exp(p.nl * duration);
                p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
                p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
                p.qtyG *= expAbsorb;
                time = timeConc.first;
            }else{ assert( time == timeConc.first ); }
            // add to quantity of drug in gut:
            p.qtyG += timeConc.second;   // units: mg
            dosesTaken += 1;
        }else /*i.e. tomorrow or later*/{
            timeConc.first -= 1.0;
        }
    }
    if( time < 1.0 ){
        duration = 1.0 - time;
        const double expAbsorb = exp(m_nka * duration), expPLoss = exp(p.nl * duration);
        p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
        p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
        p.qtyG *= expAbsorb;
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    m_doses.erase(m_doses.begin(), m_doses.begin() + dosesTaken);
    
    m_qtyM = p.qtyM;
    m_qtyP = p.qtyP;
    m_qtyG = p.qtyG;
    
    util::streamValidate( m_qtyM );
    if( m_qtyP < m_parentType.getNegligibleConcentration() * bodyMass * m_volDist &&
            m_qtyM < m_metaboliteType.getNegligibleConcentration() * bodyMass * m_volDistMetabolite )
    {
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        m_qtyG = m_qtyP = m_qtyM = 0.0;
    }
}

void LSTMDrugConversion::checkpoint(ostream& stream){
    m_qtyG & stream;
    m_qtyP & stream;
    m_qtyM & stream;
    m_nka & stream;
    m_nkPSample & stream;
    m_nconvSample & stream;
    m_nkMSample & stream;
}
void LSTMDrugConversion::checkpoint(istream& stream){
    m_qtyG & stream;
    m_qtyP & stream;
    m_qtyM & stream;
    m_nka & stream;
    m_nkPSample & stream;
    m_nconvSample & stream;
    m_nkMSample & stream;
}

}
}
