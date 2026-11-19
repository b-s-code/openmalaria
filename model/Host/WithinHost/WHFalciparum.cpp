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

#include "Host/WithinHost/WHFalciparum.h"
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "Host/WithinHost/Infection/DummyInfection.h"
#include "Host/WithinHost/Infection/EmpiricalInfection.h"
#include "Host/WithinHost/Infection/MolineauxInfection.h"
#include "Host/WithinHost/Infection/PennyInfection.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Treatments.h"
#include "Host/WithinHost/Genotypes.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/checkpoint_containers.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"

#include <cmath>
#include <gsl/gsl_cdf.h>


namespace OM {
namespace WithinHost {

using namespace OM::util;

// -----  private model parameters  -----

//Standard dev innate immunity for densities
static double sigmaI;

// contribution of parasite densities to acquired immunity in the presence of fever
static double immPenalty_22;

// Remaining immunity against asexual parasites(after time step, each of 2 components y and h)
// This variable decays the effectors m_cumulative_h and m_cumulative_Y in a way that their
// effects on densities (1-Dh and 1-Dy) decay exponentially.
static double asexImmRemain;

// Remaining immunity against asexual parasites(after each time step, each of 2 components y and h)
// This variable decays the effectors m_cumulative_h and m_cumulative_Y exponentially.
static double immEffectorRemain;

/// Critical value for immunity trigger (cumulative densities)
static double invCumulativeYstar;
/// Critical value for immunity trigger (cumulative inoculations)
static double invCumulativeHstar;

/// Maternal protection at birth
static double alphaM;

/// More or less (up to 0.693) inverse quantity of alphaMStar (AJTM p. 9 eq. 12),
/// decay rate of maternal protection in years^(-1).
static double decayM;

SimTime Infection::s_latentP = Sim::never();
int WHFalciparum::s_yLagLen = 0;



// -----  static functions  -----

void WHFalciparum::init( const OM::Parameters& parameters, const scnXml::Model& model ) {
    sigmaI=sqrt(parameters[Parameter::SIGMA_I_SQ]);
    immPenalty_22=1-exp(parameters[Parameter::IMMUNITY_PENALTY]);
    immEffectorRemain=exp(-parameters[Parameter::IMMUNE_EFFECTOR_DECAY]);
    asexImmRemain=exp(-parameters[Parameter::ASEXUAL_IMMUNITY_DECAY]);
    
    // calculate inverses here, so we can use multiplication later (faster):
    invCumulativeYstar = 1.0 / parameters[Parameter::CUMULATIVE_Y_STAR];
    invCumulativeHstar = 1.0 / parameters[Parameter::CUMULATIVE_H_STAR];
    alphaM = 1.0 - exp(-parameters[Parameter::NEG_LOG_ONE_MINUS_ALPHA_M]);
    decayM = parameters[Parameter::DECAY_M];
    
    s_yLagLen = Sim::inSteps(Sim::fromDays(20)) + 1;
    
    //NOTE: should also call cleanup() on the PathogenesisModel, but it only frees memory which the OS does anyway
    Pathogenesis::PathogenesisModel::init( parameters, model.getClinical(), false );
    
    // Default value to use if no <parameters> element is present in XML.
    SimTime latentP = UnitParse::readShortDuration("15d", UnitParse::DAYS);
    try
    {
        if (model.getParameters().present())
        {
            latentP = UnitParse::readShortDuration(
                model.getParameters().get().getLatentp(), UnitParse::NONE);
        }
    }
    catch( const util::FormatError& e )
    {
            throw util::XmlScenarioError( string("model/parameters/latentP: ").append(e.message()) );
    }

    Infection::init( latentP );
}

void WHFalciparum::setParams(double cumYStar, double cumHStar, double aM, double dM) {
    invCumulativeYstar = cumYStar;
    invCumulativeHstar = cumHStar;
    alphaM = aM;
    decayM = dM;
}


// -----  Non-static  -----

WHFalciparum::WHFalciparum( LocalRng& rng, double comorbidityFactor ):
    WHInterface(),
    m_cumulativeH(0.0), m_cumulativeY(0.0), m_cumulativeYLag(0.0),
    m_totalDensity(0.0), m_hrp2Density(0.0), m_timeStepMaxDensity(0.0),
    m_pathogenesisModel( Pathogenesis::PathogenesisModel::createPathogenesisModel( comorbidityFactor ) ),
    m_treatExpiryLiver(0), m_treatExpiryBlood(0)
{
    // NOTE: negating a Gaussian sample with mean 0 is pointless — except that
    // the individual samples change. In any case the overhead is negligible.
    //FIXME: Should this be allowed to be greater than 1?
    // Oldest code on GoogleCode: _innateImmunity=(double)(W_GAUSS((0), (sigma_i)));
    m_innateImmSurvFact = exp(-rng.gauss(0.0, sigmaI));
    
    m_yLagI.resize(s_yLagLen * Genotypes::N());
    m_yLagL.resize(s_yLagLen * Genotypes::N());
}

WHFalciparum::~WHFalciparum()
{
}

double WHFalciparum::immunitySurvivalFactor (double ageInYears, double cumulativeExposureJ) {
    if (std::isnan(ageInYears) || std::isnan(cumulativeExposureJ) ||
        std::isnan(m_cumulativeH) || std::isnan(m_cumulativeY)) {
        throw BaseException("nan in immunitySurvivalFactor");
    }
    
    // Documentation: AJTMH pp22-23
    
    // Effect of cumulative Parasite density (named Dy in AJTM)
    double dY = 1.0;
    // Effect of number of infections experienced since birth (named Dh in AJTM)
    double dH = 1.0;
    if (m_cumulativeH > 1.0) {
        dH = 1.0 / (1.0 + (m_cumulativeH - 1.0) * invCumulativeHstar);
        dY = 1.0 / (1.0 + (m_cumulativeY - cumulativeExposureJ) * invCumulativeYstar);
    }
    
    // Effect of age-dependent maternal immunity (named Dm in AJTM)
    double dA = 1.0 - alphaM * exp(-decayM * ageInYears);
    
    return util::streamValidate( std::min(dY*dH*dA, 1.0) );
}

// Infectiousness parameters: see AJTMH p.33; tau=1/sigmag**2 
const double PTM_beta1=1.0;
const double PTM_beta2=0.46;
const double PTM_beta3=0.17;
const double PTM_tau= 0.066;
const double PTM_tau_prime = 1.0 / sqrt(1.0 / PTM_tau);
const double PTM_mu= -8.1;

inline double computeYLag(const vector<double> &yLag, int yLagLen, vector<double> &yLagG)
{
    const size_t n = Genotypes::N();

    // Take weighted sum of total asexual blood stage density 10, 15 and 20 days
    // before. Add y_lag_len to index to ensure positive.
    size_t d10 = modNn(yLagLen + Sim::inSteps(Sim::ts1() - Sim::fromDays(10)), yLagLen);
    size_t d15 = modNn(yLagLen + Sim::inSteps(Sim::ts1() - Sim::fromDays(15)), yLagLen);
    size_t d20 = modNn(yLagLen + Sim::inSteps(Sim::ts1() - Sim::fromDays(20)), yLagLen);

    double y10Sum = 0.0, y15Sum = 0.0, y20Sum = 0.0;
    for( size_t g = 0; g < n; ++g ){
        double y10 = yLag[d10 * n + g];
        double y15 = yLag[d15 * n + g];
        double y20 = yLag[d20 * n + g];
        y10Sum += y10;
        y15Sum += y15;
        y20Sum += y20;
        yLagG[g] = PTM_beta1 * y10 + PTM_beta2 * y15 + PTM_beta3 * y20;
    }
    
    // Weighted sum:
    return PTM_beta1 * y10Sum + PTM_beta2 * y15Sum + PTM_beta3 * y20Sum;
}

double WHFalciparum::probTransmissionToMosquito(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l) const{
    // This model (often referred to as the gametocyte model) was designed for
    // 5-day time steps. We use the same model (sampling 10, 15 and 20 days
    // ago) for 1-day time steps to avoid having to design and analyse a new
    // model. Description: AJTMH pp.32-33 and p9.
    
    // Note: we don't allow for gametocydal treatments (e.g. Primaquine).
    const size_t n = Genotypes::N();

    // Sum lagged densities across genotypes:
    vector<double> yLagGI(n), yLagGL(n);
    const double yLagSumI = computeYLag(m_yLagI, s_yLagLen, yLagGI);
    const double yLagSumL = computeYLag(m_yLagL, s_yLagLen, yLagGL);
    const double yLagSum = yLagSumI + yLagSumL;

    if( yLagSum < 0.001 ) return 0.0; // cut off for uninfectious humans
    
    // Get a zval, convert to equivalent Normal sample:
    const double zval = (log(yLagSum) + PTM_mu) * PTM_tau_prime;
    const double pone = gsl_cdf_ugaussian_P(zval);
    double pTransmit = pone*pone;

    // pTransmit has to be between 0 and 1:
    pTransmit=std::max(pTransmit, 0.0);
    pTransmit=std::min(pTransmit, 1.0);
    
    if(pTransmit <= 0.0)
        return pTransmit;

    for( size_t g = 0; g < n; ++g )
    {
        probTransGenotype_i[g] = pTransmit * yLagGI[g] / yLagSum;
        probTransGenotype_l[g] = pTransmit * yLagGL[g] / yLagSum;
    }

    // Include here the effect of transmission-blocking vaccination:
    util::streamValidate( pTransmit );
    return pTransmit;
}

bool WHFalciparum::diagnosticResult( LocalRng& rng, const Diagnostic& diagnostic ) const{
    return diagnostic.isPositive( rng, m_totalDensity, m_hrp2Density );
}

void WHFalciparum::treatment( Host::Human& human, TreatmentId treatId ){
    const Treatments& treat = Treatments::select( treatId );
    treatSimple( human, treat.liverEffect(), treat.bloodEffect() );
    
    // triggered intervention deployments:
    treat.deploy( human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
}
bool WHFalciparum::treatSimple( Host::Human& human, SimTime timeLiver, SimTime timeBlood ){
    if( timeLiver != Sim::zero() ){
        if( timeLiver < Sim::zero() )
            clearInfections( Treatments::LIVER );
        else
            m_treatExpiryLiver = max( int(m_treatExpiryLiver), int(Sim::nowOrTs1()) + timeLiver );
        mon::recordEvent(mon::measure("nLiverStageTreatments"), human);
    }
    if( timeBlood != Sim::zero() ){
        if( timeBlood < Sim::zero() )
            clearInfections( Treatments::BLOOD );
        else
            m_treatExpiryBlood = max( int(m_treatExpiryBlood), int(Sim::nowOrTs1()) + timeBlood );
        return true;    // blood stage treatment
    }
    return false;    // no blood stage treatment
}

Pathogenesis::StatePair WHFalciparum::determineMorbidity( Host::Human& human,
        double ageYears, bool isDoomed )
{
    Pathogenesis::StatePair result = m_pathogenesisModel->determineState( human,
            ageYears, m_timeStepMaxDensity, m_totalDensity, isDoomed );
    
    /* Note: this model can easily be re-enabled, but is not used and not considered to be a good model.
    if( (result.state & Pathogenesis::MALARIA) && util::ModelOptions::option( util::PENALISATION_EPISODES ) ){
        // This does immunity penalisation:
        m_cumulative_Y = m_cumulative_Y_lag - immPenalty_22*(m_cumulative_Y-m_cumulative_Y_lag);
        if (m_cumulative_Y < 0) {
            m_cumulative_Y=0.0;
        }
    }*/
    
    return result;
}


// -----  immunity  -----

void WHFalciparum::updateImmuneStatus() {
    if (immEffectorRemain < 1) {
        m_cumulativeH *= immEffectorRemain;
        m_cumulativeY *= immEffectorRemain;
    }
    if (asexImmRemain < 1) {
        m_cumulativeH *= asexImmRemain /
                      (1+(m_cumulativeH*(1-asexImmRemain) * invCumulativeHstar));
        m_cumulativeY *= asexImmRemain /
                      (1+(m_cumulativeY*(1-asexImmRemain) * invCumulativeYstar));
    }
    m_cumulativeYLag = m_cumulativeY;
}


// -----  Checkpointing  -----

void WHFalciparum::checkpoint (istream& stream) {
    WHInterface::checkpoint( stream );
    m_innateImmSurvFact & stream;
    m_cumulativeH & stream;
    m_cumulativeY & stream;
    m_cumulativeYLag & stream;
    m_totalDensity & stream;
    m_hrp2Density & stream;
    m_timeStepMaxDensity & stream;
    m_yLagI & stream;
    m_yLagL & stream;
    (*m_pathogenesisModel) & stream;
    m_treatExpiryLiver & stream;
    m_treatExpiryBlood & stream;
}
void WHFalciparum::checkpoint (ostream& stream) {
    WHInterface::checkpoint( stream );
    m_innateImmSurvFact & stream;
    m_cumulativeH & stream;
    m_cumulativeY & stream;
    m_cumulativeYLag & stream;
    m_totalDensity & stream;
    m_hrp2Density & stream;
    m_timeStepMaxDensity & stream;
    m_yLagI & stream;
    m_yLagL & stream;
    (*m_pathogenesisModel) & stream;
    m_treatExpiryLiver & stream;
    m_treatExpiryBlood & stream;
}

}
}
