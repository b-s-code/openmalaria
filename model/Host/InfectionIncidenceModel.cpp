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

#include "Host/InfectionIncidenceModel.h"
#include "Host/Human.h"
#include "Transmission/PerHost.h"
#include "Host/WithinHost/WHInterface.h"
#include "Parameters.h"
#include "mon/Continuous.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/errors.h"

#include <stdexcept>
#include <cmath>

namespace OM { namespace Host {
    using namespace OM::util;

// ———  model constants set by init()  ———

/* Shape constant of (Gamma) distribution of availability
real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
double baselineAvailShapeParam, baselineAvailOffset;

//VARIABLES INCLUDED IN CORE GETs of number of infections 
//! Describes the shape of the Infectionrate distribution, related to the baseline availabilty distr.
double infRateShapeParam, infRateOffset;

/** @brief Parameters used in the "expected number of infections" model */
//@{
//!Steepness of relationship between success of inoculation and Xp in Phase A model 
double gammaP;
//!Lower limit of success probability of inoculations at high exposure in Phase A model 
double sinf;
//!Lower limit of success probability of inoculations in immune individuals in Phase A model 
double simm;
//!1 over the critical value of cumulative number of entomologic inoculations in Phase A model 
double Xstar_pInv;
//!1 over the critical value of EIR in Phase A pre-erythrocytic model 
double estarInv;
//@}

// model options:
bool optNegBinMassAction = false, optLognormalMassAction = false,
        optNoPreErythrocytic = false, optAnyHet = false, optVaccineGenotype = false;

// ———  variables  ———
int InfectionIncidenceModel::s_ctsNewInfections = 0;

// -----  static initialisation  -----

void InfectionIncidenceModel::init ( const Parameters& parameters ) {
    baselineAvailShapeParam = parameters[Parameter::BASELINE_AVAILABILITY_SHAPE];
    baselineAvailOffset = 0.5 * pow(baselineAvailShapeParam, 2);
    
    gammaP=parameters[Parameter::GAMMA_P];
    sinf=1-exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_SINF]);
    simm=parameters[Parameter::SIMM];
    estarInv = 1.0/parameters[Parameter::E_STAR];
    Xstar_pInv = 1.0/parameters[Parameter::X_STAR_P];
    
    //! constant defining the constraint for the Gamma shape parameters
    /// Used for the case where availability is assumed gamma distributed
    double r_square_Gamma;
    /*
    //! Expected number of inoculations
    /// Product of measured EIR, susceptibility and length of time Global::interval
    double gsi = 1.0;
    
    r_square_Gamma=(totalInfectionrateVariance**2-gsi*BaselineAvailabilityMean)/(gsi*BaselineAvailabilityMean)**2
    r_square_Gamma must be greater than zero, so r_square_LogNormal is also. 
    */
    r_square_Gamma=0.649; //such that r_square_LogNormal =0.5
    
    optNoPreErythrocytic = util::ModelOptions::option (util::NO_PRE_ERYTHROCYTIC);
    optNegBinMassAction = util::ModelOptions::option (util::NEGATIVE_BINOMIAL_MASS_ACTION);
    optVaccineGenotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    if (optNegBinMassAction) {
        infRateShapeParam = (baselineAvailShapeParam+1.0) / (r_square_Gamma*baselineAvailShapeParam - 1.0);
        infRateShapeParam=std::max(infRateShapeParam, 0.0);
    } else if (util::ModelOptions::option (util::LOGNORMAL_MASS_ACTION)) {
        optLognormalMassAction = true;
        
        //! constant defining the constraint for the log Normal variance
        /// Used for the case where availability is assumed log Normally distributed
        double r_square_LogNormal = log(1.0+r_square_Gamma);
        
        infRateShapeParam = sqrt(r_square_LogNormal - 1.86 * pow(baselineAvailShapeParam, 2));
        infRateShapeParam=std::max(infRateShapeParam, 0.0);
        infRateOffset = 0.5 * pow(infRateShapeParam, 2);
        if( (std::isnan)(infRateShapeParam) ){
            throw util::XmlScenarioError( "bad parameter 16 (BASELINE_AVAILABILITY_SHAPE)" );
        }
    }else{
        if( util::ModelOptions::option( util::TRANS_HET ) ||
            util::ModelOptions::option( util::COMORB_TRANS_HET ) ||
            util::ModelOptions::option( util::TRANS_TREAT_HET ) ||
            util::ModelOptions::option( util::TRIPLE_HET ) )
        {
            optAnyHet = true;
            cerr << "Warning: will use heterogeneity workaround." << endl;
        }
    }
    
    mon::Continuous::registerCallback( "new infections", "\tnew infections", &InfectionIncidenceModel::ctsReportNewInfections );
}


// -----  other static methods  -----

void InfectionIncidenceModel::ctsReportNewInfections (const vector<Host::Human> &, ostream& stream){
    stream << '\t' << s_ctsNewInfections;
    s_ctsNewInfections = 0;
}


// -----  non-static non-checkpointing constructors  -----

InfectionIncidenceModel* InfectionIncidenceModel::createModel () {
    if (optNegBinMassAction) return new NegBinomMAII ();
    if(optLognormalMassAction) return new LogNormalMAII ();
    if (optAnyHet) return new HeterogeneityWorkaroundII();
    return new InfectionIncidenceModel();
}

InfectionIncidenceModel::InfectionIncidenceModel () :
  m_pInfected(0.0), m_cumulativeEIRa(0.0)
{}


// -----  non-static methods  -----

double InfectionIncidenceModel::getAvailabilityFactor(LocalRng& rng, double baseAvailability) {
  return baseAvailability;
}
double NegBinomMAII::getAvailabilityFactor(LocalRng& rng, double baseAvailability) {
    // Gamma sample with k=BaselineAvailabilityShapeParam; mean is baseAvailability
  return rng.gamma(baselineAvailShapeParam,
		 baseAvailability/baselineAvailShapeParam);
}
double LogNormalMAII::getAvailabilityFactor(LocalRng& rng, double baseAvailability) {
    // given BaselineAvailabilityShapeParam = sqrt (log (1 + variance/mean²))
    // and baseAvailability = mean, this is a draw from the log-normal distribution.
    if( baseAvailability != 1.0 ){
        // NOTE: shouldn't the normal_mean parameter be adjusted when baseAvailability != 1.0?
        throw TRACED_EXCEPTION_DEFAULT("LogNormalMAII::getAvailabilityFactor");
    }
  return rng.logNormal (log(baseAvailability) - baselineAvailOffset,
		     baselineAvailShapeParam);
}

void InfectionIncidenceModel::summarize (const Host::Human& human) {
    mon::recordStat(mon::measure("nExpectd"), human, m_pInfected);
}


double InfectionIncidenceModel::getModelExpectedInfections (LocalRng& rng, double effectiveEIR, const Transmission::PerHost& phTrans) {
  // First two lines are availability adjustment: S_1(i,t) from AJTMH 75 (suppl 2) p12 eqn. (5)
  // Note that NegBinomMAII and LogNormalMAII supercede this model; see below
  return (sinf+(1-sinf) / 
    (1 + effectiveEIR/Sim::oneTS()*estarInv)) *
    susceptibility() * effectiveEIR;
}
double HeterogeneityWorkaroundII::getModelExpectedInfections (LocalRng& rng, double effectiveEIR, const Transmission::PerHost& phTrans) {
  return (sinf+(1-sinf) / 
    (1 + effectiveEIR/(Sim::oneTS()*phTrans.relativeAvailabilityHet)*estarInv)) *
    susceptibility() * effectiveEIR;
}
double NegBinomMAII::getModelExpectedInfections (LocalRng& rng, double effectiveEIR, const Transmission::PerHost&) {
  // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pmed.1001157.s009
  return rng.gamma(infRateShapeParam,
      effectiveEIR * susceptibility() / infRateShapeParam);
}
double LogNormalMAII::getModelExpectedInfections (LocalRng& rng, double effectiveEIR, const Transmission::PerHost&) {
    // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pmed.1001157.s009
    if (effectiveEIR < 1e-10)
        return 0.0;
    double meanlog = log(effectiveEIR * susceptibility()) - infRateOffset;
    return rng.logNormal(meanlog, infRateShapeParam);
}

double InfectionIncidenceModel::susceptibility () {
  if (optNoPreErythrocytic) {
    //! The average proportion of bites from sporozoite positive mosquitoes resulting in infection. 
    /*! 
    This is computed as 0.19 (the value S from a neg bin mass action model fitted 
    to Saradidi data, divided by 0.302 (the ratio of body surface area in a 
    0.5-6 year old child (as per Saradidi) to adult) 
    \sa getExpectedNumberOfInfections() 
    */ 
    return 0.702;
  } else {
    // S_2(i,t) from AJTMH 75 (suppl 2) p12 eqn. (7)
    return simm + (1.0-simm) /
      (1.0 + pow(m_cumulativeEIRa*Xstar_pInv, gammaP));
  }
}

double InfectionIncidenceModel::expectedNumNewInfections(OM::Host::Human& human, double effectiveEIR)
{
    // error check (should be OK if kappa is checked, for nonVector model):
    if( !(std::isfinite)(effectiveEIR) ){
        ostringstream out;
        out << "effectiveEIR is not finite: " << effectiveEIR << endl;
        throw TRACED_EXCEPTION (out.str(), util::Error::EffectiveEIR);
    }

    double expectedNumInfections = getModelExpectedInfections (human.rng, effectiveEIR, human.perHostTransmission);    
    
    // Only to be consistent with old simulation runs when set to false
    // Setting this option to true will only affect reporting
    if(optVaccineGenotype == false)
        //Introduce the effect of vaccination. Note that this does not affect cumEIR.
        expectedNumInfections *= human.vaccine.getFactor( interventions::Vaccine::PEV );
  
    //Update pre-erythrocytic immunity
    m_cumulativeEIRa += effectiveEIR;

    m_pInfected = 1.0 - exp(-expectedNumInfections) * (1.0-m_pInfected);
    if (m_pInfected < 0.0)
        m_pInfected = 0.0;
    else if (m_pInfected > 1.0)
        m_pInfected = 1.0;
    
    return expectedNumInfections;
}

int InfectionIncidenceModel::numNewInfections (Human& human, double expectedNumInfections)
{ 
  if (expectedNumInfections > 0.0000001){
    int n = human.rng.poisson(expectedNumInfections);

    if( n > WithinHost::WHInterface::MaxInfections )
        n = WithinHost::WHInterface::MaxInfections;
    return n;
  }
  
  if ( (std::isnan)(expectedNumInfections) ){	// check for not-a-number
      // bad Params::BASELINE_AVAILABILITY_SHAPE ?
      throw TRACED_EXCEPTION( "numNewInfections: NaN", util::Error::NumNewInfections );
  }
  return 0;
}

void InfectionIncidenceModel::reportNumNewInfections(Human& human, int newNumInfections_i, int newNumInfections_l)
{
    mon::recordEvent(mon::measure("nNewInfections"), human, newNumInfections_i+newNumInfections_l);
    mon::recordEvent(mon::measure("nNewInfections_Introduced"), human, newNumInfections_i);
    mon::recordEvent(mon::measure("nNewInfections_Indigenous"), human, newNumInfections_l);
    s_ctsNewInfections += newNumInfections_i+newNumInfections_l;
}

} }
