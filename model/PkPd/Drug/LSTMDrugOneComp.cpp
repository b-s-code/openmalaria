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

#include "PkPd/Drug/LSTMDrugOneComp.h"
#include "Host/WithinHost/Infection/CommonInfection.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrugOneComp::LSTMDrugOneComp(const LSTMDrugType& type, LocalRng& rng) :
    LSTMDrug(type.sampleVd(rng)),
    m_typeData(type),
    m_concentration (0.0),
    m_negElimSample(-type.sampleElimRate(rng))
{}

//LSTMDrugOneComp::~LSTMDrugOneComp(){}

size_t LSTMDrugOneComp::getIndex() const {
    return m_typeData.getIndex();
}
double LSTMDrugOneComp::getConcentration(size_t index) const {
    if( index == m_typeData.getIndex() ) return m_concentration;
    else return 0.0;
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrugOneComp::calculateDrugFactor(LocalRng& rng, WithinHost::CommonInfection *inf, double bodyMass) const {
    if( m_concentration == 0.0 && m_doses.size() == 0 ) return 1.0; // nothing to do
    
    /* Survival factor of the parasite (this multiplies the parasite density).
    Calculated below for each time interval. */
    double totalFactor = 1.0;
    
    // Make a copy of concetration and use that over today. Don't adjust concentration because this
    // function may be called multiple times (or not at all) in a day.
    double concentrationToday = m_concentration; // mg / l
    double negElimRate = m_negElimSample * pow(bodyMass, m_typeData.negMExponent());
    
    const LSTMDrugPD& drugPD = m_typeData.getPD(inf->genotype());
    const double Kn = drugPD.IC50_pow_slope(rng, m_typeData.getIndex(), inf);
    
    double time = 0.0;
    typedef pair<double,double> TimeConc;
    for( TimeConc timeConc : m_doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( timeConc.first < 1.0 /*i.e. today*/ ){
            if( time < timeConc.first ){
                totalFactor *= drugPD.calcFactor( Kn, negElimRate, &concentrationToday, timeConc.first - time );
                time = timeConc.first;
            }else{ assert( time == timeConc.first ); }
            // add dose (instantaneous absorption):
            concentrationToday += timeConc.second / (m_volDist * bodyMass);
        }else/*i.e. tomorrow or later*/{
            break;
        }
    }
    if( time < 1.0 ){
        totalFactor *= drugPD.calcFactor( Kn, negElimRate, &concentrationToday, 1.0 - time );
    }
    
    return totalFactor; // Drug effect per day per drug per parasite
}

void LSTMDrugOneComp::updateConcentration( double bodyMass ){
    if( m_concentration == 0.0 && m_doses.size() == 0 ) return;     // nothing to do
    
    // exponential decay of drug concentration (portion without new doses):
    //TODO: is it faster to pre-calculate this and either store an extra
    // parameter or adapt uses of neg_elim_rate, etc. below?
    double negElimRate = m_negElimSample * pow(bodyMass, m_typeData.negMExponent());
    m_concentration *= exp(negElimRate);
    size_t dosesTaken = 0;
    typedef pair<double,double> TimeConc;
    for( TimeConc& timeConc : m_doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( timeConc.first < 1.0 /*i.e. today*/ ){
            // calculate decayed dose and add:
            m_concentration += timeConc.second / (m_volDist * bodyMass) * exp(negElimRate * (1.0 - timeConc.first));
            dosesTaken += 1;
        }else /*i.e. tomorrow or later*/{
            timeConc.first -= 1.0;
        }
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    m_doses.erase(m_doses.begin(), m_doses.begin() + dosesTaken);
    
    util::streamValidate( m_concentration );
    if( m_concentration < m_typeData.getNegligibleConcentration() ){
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        m_concentration = 0.0;
    }
}

void LSTMDrugOneComp::checkpoint(ostream& stream){
    m_concentration & stream;
    m_negElimSample & stream;
}
void LSTMDrugOneComp::checkpoint(istream& stream){
    m_concentration & stream;
    m_negElimSample & stream;
}

}
}
