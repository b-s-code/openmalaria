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

#include "Host/WithinHost/Infection/MolineauxInfection.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/checkpoint_containers.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

// ———  model constants set from input parameters  ———

// depends on lognormal/gamma distribution for first_local_max
static bool firstLocalMaximumGamma = false;
// mean/shape and sd/scale of the first local maximum density for lognormal/gamma distribution
static double meanShapeFirstLocalMax = std::numeric_limits<double>::signaling_NaN();
static double sdScaleFirstLocalMax = std::numeric_limits<double>::signaling_NaN();

// depends on between lognormal/gamma distribution for mean_duration diff_pos_days
static bool meanDurationGamma = false;
// mean/shape and sd/scale of the difference between last positive and first positive days for lognormal/gamma distribution
static double meanShapeDiffPosDays = std::numeric_limits<double>::signaling_NaN();
static double sdScaleDiffPosDays = std::numeric_limits<double>::signaling_NaN();

// boolean choosing between gamma and lognormal distribution for equation 11
static bool multiFactorGamma = false;

// pairwise sample of case-specific P* parameters
static bool pairwisePStarSample = false;

// q^(i+1) array
// All the values of q^1... q^v are stored in this array.
// This avoids the recalculation of those values every second time step. */
static double qPow[MolineauxInfection::V];

// ———  hard-coded model constants  ———

// sigma, rho: decay parameters, per day, of the acquired variant-specific and variant-transcending immune responses
const double Sigma = 0.02;
const double SigmaDecay = exp(-2.0*Sigma);
// const double rho = 0.0;      // NOTE: rho is not used; if setting this, also change the code below (search rho)
// beta: Minimum value of the probability that a parasite escapes control by the acquired and variant-transcending immune response
const double Beta = 0.01;
// s: the fraction of parasites switching among variants per two-day cycle
const double S = 0.02;
// q: Parameter of the geometric distribution of switching probabilities
const double Q = 0.3;

// mu_m, sigma_m: Mean and standard deviation to use for the normal distribution setting the variant specific multiplication factor.
const double MuM = 16.0;
const double SigmaM = 10.4;
// constants for sampling multiplication factors from a Gamma distribution instead of normal (above)
const double ShapeM = 2.4;
const double ScaleM = 6.8;

// k_c,k_m: constants allowing calculation of Pc_star and Pm_star from host-specific data
const double KC = 0.2;
const double KM = 0.04;
// Pv_star: critical density of a variant, common to all variants
const double inv_Pv_star = 1.0 / 30.0;
// kappa_c, kappa_m, kappa_v: Stiffness parameters for saturation of immune responses
const int KappaC = 3;
const int KappaM = 1;
const int KappaV = 3;
// C: Maximum daily antigenic stimulus, per mul, of the acquired variant-transcending immune response
const double C = 1.0;

// Blood volume per kg can be found, for example, here:
// http://hypertextbook.com/facts/1998/LanNaLee.shtml
// Futher searches give values between around 65 and 85 ml/kg.
const double BloodVolPerKg = 7e4;    // μl/kg; (equals 70 ml/kg)
// PRBC = parasited red blood cell
const double InitialDens = 0.1;        // initial parasite density in PRBC/μl
const double ElimParasites = 50;       // minimum number PRBC to avoid elimination of infection and variant

/* Case-specific parameters.
 * 
 * For each of the 35 Malaria Therapy patients, this is:
 *  0) the duration of the infection as last day with positive density minus first positive
 *  1) the parasite density in parasites/microlitre at the first local maximum
 * 
 * Indexing is 2*patient + param where 0 <= patient < 35 and 0 <= param < 2.
 */
const double CaseSpecificData[] = {
        216, 18600, // G131
        198, 13080, // G140
        206, 45720, // G141
        366, 23760, // G142
        230, 60840, // G143
        172, 6000, // G146
        100, 2340, // G147
        236, 31440, // G159
        236, 453600, // G161
        120, 4240, // G173
        176, 195840, // G174
        178, 60120, // G178
        36, 8720, // G184
        44, 8000, // G193
        242, 395280, // G200
        70, 28320, // G210
        292, 200160, // G212
        248, 59320, // G217
        98, 66480, // G23
        176, 61200, // G230
        234, 169920, // G240
        226, 46800, // G264
        270, 19260, // G265
        278, 86040, // G279
        212, 110160, // G290
        264, 43200, // G385
        364, 133920, // G407
        184, 222480, // G408
        160, 21420, // G414
        220, 74160, // G416
        132, 210960, // G423
        176, 89280, // G439
        208, 105840, // G445
        330, 21600, // G457
        404, 156240 // G48
};

// ———  static (non-class-member) code  ———

CommonInfection* createMolineauxInfection (LocalRng& rng, uint32_t protID, int origin) {
    return new MolineauxInfection (rng, protID, origin);
}

CommonInfection* checkpointedMolineauxInfection (istream& stream) {
    return new MolineauxInfection (stream);
}

void MolineauxInfection::init( const Parameters& parameters ){
    CommonWithinHost::s_createInfection = &createMolineauxInfection;
    CommonWithinHost::s_checkpointedInfection = &checkpointedMolineauxInfection;
    
    multiFactorGamma = util::ModelOptions::option( util::PARASITE_REPLICATION_GAMMA );
    
    pairwisePStarSample = util::ModelOptions::option( util::MOLINEAUX_PAIRWISE_SAMPLE );
    if( !pairwisePStarSample ){
        meanShapeFirstLocalMax = parameters[Parameter::MEAN_LOCAL_MAX_DENSITY];
        sdScaleFirstLocalMax = parameters[Parameter::SD_LOCAL_MAX_DENSITY];
        
        meanShapeDiffPosDays = parameters[Parameter::MEAN_DIFF_POS_DAYS];
        sdScaleDiffPosDays = parameters[Parameter::SD_DIFF_POS_DAYS];
        
        
        // with gamma distribution shape and scale parameters has to be recalculated 
        if (util::ModelOptions::option (util::FIRST_LOCAL_MAXIMUM_GAMMA)) {
            //NOTE: variant does not appear to work well (according to unit test stats)
            firstLocalMaximumGamma = true;
            meanShapeFirstLocalMax = pow(meanShapeFirstLocalMax,2) / pow(sdScaleFirstLocalMax,2);
            sdScaleFirstLocalMax = pow(sdScaleFirstLocalMax,2) / meanShapeFirstLocalMax;
        } else {
            firstLocalMaximumGamma = false;
        }
        
        // with gamma distribution shape and scale parameters has to be recalculated
        if(util::ModelOptions::option (util::MEAN_DURATION_GAMMA)) {
            //NOTE: variant does not appear to work well (according to unit test stats)
            meanDurationGamma = true;
            meanShapeDiffPosDays = pow(meanShapeDiffPosDays,2) / pow(sdScaleDiffPosDays,2);
            sdScaleDiffPosDays = pow(sdScaleDiffPosDays,2) / meanShapeDiffPosDays;
        } else {
            meanDurationGamma = false;
        }
    }
    
   for( size_t i = 0; i < V; i++ ){
       qPow[i] = pow(Q, static_cast<double>(i+1));
   }
}

// ———  MolineauxInfection: initialisation  ———

MolineauxInfection::MolineauxInfection(LocalRng& rng, uint32_t genotype, int origin):
        CommonInfection(genotype, origin)
{
    for( size_t i = 0; i < V; i++ ){
        // Molineaux paper, equation 11
        if( multiFactorGamma ){
            do{
                m_mi[i] = static_cast<float>(rng.gamma(ShapeM,ScaleM));
            }while( m_mi[i]<1.0 );
        }else{
            do{
                m_mi[i] = static_cast<float>(rng.gauss(MuM, SigmaM));
            }while( m_mi[i]<1.0 );
        }
    }
    
    for( size_t tau=0; tau<Taus; tau++ ){
        lagged_Pc[tau] = 0.0;
    }
    
    Sm_summation = 0.0;
    
    if( pairwisePStarSample ){
        int patient = rng.uniform( 35 );
        Pc_star = static_cast<float>( KC * CaseSpecificData[2 * patient + 1] );
        Pm_star = static_cast<float>( KM * CaseSpecificData[2 * patient + 0] );
    }else{
        // Sampling of the first local maxima:
        if( firstLocalMaximumGamma ){
            Pc_star = static_cast<float>( KC * pow(10.0, rng.gamma(
                meanShapeFirstLocalMax, sdScaleFirstLocalMax)) );
        } else {
            Pc_star = static_cast<float>( KC * pow(10.0, rng.gauss(
                meanShapeFirstLocalMax, sdScaleFirstLocalMax)) );
        }
        
        // Sampling of duration:
        if( meanDurationGamma ) {
            Pm_star = static_cast<float>( KM * pow(10.0, rng.gamma(
                meanShapeDiffPosDays,sdScaleDiffPosDays)) );
        } else {
            Pm_star = static_cast<float>( KM * pow(10.0, rng.gauss(
                meanShapeDiffPosDays,sdScaleDiffPosDays)) );
        }
    }
}

MolineauxInfection::Variant::Variant () :
        Pi1(0.0), Pi2(0.0), Si_summation(0.0)
{
    for(size_t tau=0; tau<Taus; tau++){
        lagged_Pi[tau] =  0.0;
    }
}

// ———  MolineauxInfection: density updates  ———

bool MolineauxInfection::updateDensity( LocalRng&, double survivalFactor, SimTime age_BS, double bodyMass ){
    // bsAge : age of blood stage; 0 implies initial density (0.1; time t=0 in
    // paper, t=1 in MP's Matlab code), age 2 days is after first update step
    // survivalFactor : probabilty of merozoites surviving drugs, inter-infection immunity and vaccines
    
    double bloodVolume = BloodVolPerKg * bodyMass;
    double elimDens = ElimParasites / bloodVolume;   // 50 / 5e6 = 5e-5
    
    // ———  1. Update m_density (Pc), Pi and related  ———
    double Pi[V] = { 0.0f };
    
    if (age_BS == Sim::zero()){
        // The first variant starts with a pre-set density (regardless of blood
        // volume; this is an assumption by DH; paper assumes fixed volume)
        m_variants.resize(1);
        Pi[0] = InitialDens;
        m_density = InitialDens;
    }else{
        double sum = 0.0;
        for( size_t i = 0; i < m_variants.size(); i++ ){
            double newP = survivalFactor * m_variants[i].Pi1;
            Pi[i] = newP;
            m_variants[i].Pi1 = static_cast<float>(survivalFactor * m_variants[i].Pi2);
            sum += newP;
        }
        m_density = sum;
    }
    
    // Integration with inter-infection immunity model
    m_cumulativeExposureJ += m_density;
    
    if( m_density <= elimDens ){
        return true;    // infection goes extinct
    }
    
    // Every other step (if bsAge is even) we calculate new densities (the
    // Molineaux model uses a 2-day step; we have simply added interpolation
    // for the interleaving densities). In other cases we are finished now.
    if( modNn(int(age_BS), 2) != 0 ){
        return false;   // end of update, not extinct
    }
    
    
    // ———  2-4: Calculate infection-specific immune responses  ———
    // Sc (probability that a parasite escapes control by innate and variant-transcending immune response)
    // Sm (                        "                      acquired and variant-transcending immune response)
    // S[i] (                      "                      acquired and variant-specific immune response)
    
    // ———  2. innate immunity (equation 5)  ———
    static_assert( KappaC == 3, "kappa_c == 3" );        // allows us to optimise pow(..., kappa_c) to multiplication
    const double base = m_density/Pc_star;
    const double Sc = 1.0 / (1.0 + base*base*base);
    
    // ———  3. variant-transcending immune response (equations 7, 8)  ———
    // 3.a) Update the sum in (7), based on previous result
    const size_t tau = modNn(age_BS / 2, Taus);    // 8 days ago has same index as today
    // NOTE: rho == 0 so we can optimise this code:
    //Sm_summation = Sm_summation * exp(-2.0*rho) + laggedPc[index];
    Sm_summation = Sm_summation + lagged_Pc[tau];
    
    // 3.b) Update laggedPc (equation 8, but only stored for the last tau steps)
    //We could use min here, but it seems that min has problems with static const double C
    lagged_Pc[tau] = static_cast<float>(m_density < C ? m_density : C);
    
    // 3.c) calculate S_m(t) (equation 7)
    static_assert( KappaM == 1, "kappa_m == 1" );        // allows us to optimise out pow(..., kappa_m):
    //double Sm = (1.0 - beta) / (1.0 + pow(Sm_summation / Pm_star, kappa_m)) + beta
    const double Sm = (1.0 - Beta) / (1.0 + Sm_summation / Pm_star) + Beta;
    
    // ———  4. variant-specific immune response (equation 6)  ———
    double Si[V];       // calculate value for each variant
    double sum_qj_Sj=0.0;       // simultaneously calculation summataion in equation 4
    
    for(size_t i = 0; i < V; i++){
        if( i < m_variants.size() ){
            // 4.a) Update the sum in (6) based on the last step's value
            //note: sigma_decay = exp(-2*sigma)
            m_variants[i].Si_summation = static_cast<float>(
                m_variants[i].Si_summation * SigmaDecay + m_variants[i].lagged_Pi[tau]);
            // 4.b) update history of density (P_i(t))
            m_variants[i].lagged_Pi[tau] = static_cast<float>(Pi[i]);
            
            // 4.c) calculate S_i(t) (equation 6)
            static_assert( KappaV == 3, "kappa_v == 3" );        // again, optimise pow to multiplication
            const double base = m_variants[i].Si_summation * inv_Pv_star;
            Si[i] = 1.0 / (1.0 + base*base*base);        // eqn 6, given κ_v = 3
        }else{
            Si[i] = 1.0; // eqn 6 for the case when P_i(τ) = 0 for τ ≤ t - δ_m
        }
        sum_qj_Sj += qPow[i] * Si[i];    // sum in eqn 4
    }
    
    // ———  5. Variant densities, equations 1, 2 and 4  ———
    for(size_t i = 0; i < V; i++ ){
        // 4.a) Calculate p_i, variant selection probability (eqn 4)
        double pI = 0.0;
        if( Si[i] >= 0.1 ){
            //note: qPow[i] = pow(q, i+1)
            pI = qPow[i] * Si[i] / sum_qj_Sj;
        }
        
        // 4.b) calculate P_i'(t+2) [eqn 1] then P_i(t+2) [eqn 2]
        // This is the growth rate after taking immune effect into account:
        double growthFactor = m_mi[i] * Si[i] * Sc * Sm;   // part of eqn 1
        if( i < m_variants.size() ){
            // Pi_prime: the variant's density at time t+2 (eqn 1)
            double Pi_prime = ( (1.0 - S) * Pi[i] + S * pI * m_density ) * growthFactor;
            
            if( Pi_prime < elimDens ) Pi_prime = 0.0;    // eqn 2
            
            m_variants[i].Pi1 = static_cast<float>(sqrt(Pi[i] * Pi_prime));
            m_variants[i].Pi2 = static_cast<float>(Pi_prime);
        }else{
            // In this case P_i(τ) = 0 for all τ ≤ t, and we haven't allocated storage.
            
            // Pi_prime: the variant's density at time t+2 (eqn 1 in paper)
            double Pi_prime = ( S * pI * m_density ) * growthFactor;
            
            // Molineaux paper equation 2
            if( Pi_prime >= elimDens ){    // [if not, P_i(t+2) = 0]
                // express a new variant at time t+2:
                m_variants.resize( i+1 ); // allocate (potentially for multiple variants)
                m_variants[i].Pi2 = static_cast<float>(Pi_prime);
            }
        }
    }
    
    return false;       // end of update, not extinct
}

// ———  MolineauxInfection: checkpointing  ———

MolineauxInfection::MolineauxInfection (istream& stream) :
        CommonInfection(stream)
{
    Sm_summation & stream;
    for(size_t i=0;i<V;i++) {
        m_mi[i] & stream;
    }
    m_variants & stream;
    for(size_t j=0;j<Taus;j++){
        lagged_Pc[j] & stream;
    }
    Pc_star & stream;
    Pm_star & stream;
}

void MolineauxInfection::checkpoint (ostream& stream) {
    CommonInfection::checkpoint (stream);

    Sm_summation & stream;
    for(size_t i=0;i<V;i++) {
        m_mi[i] & stream;
    }
    m_variants & stream;
    for(size_t j=0;j<Taus;j++){
        lagged_Pc[j] & stream;
    }
    Pc_star & stream;
    Pm_star & stream;
}

void MolineauxInfection::Variant::operator& (istream& stream) {
    bool nonZero;
    nonZero & stream;
    if( nonZero ){
        Pi1 & stream;
        Pi2 & stream;
        Si_summation & stream;
        for(size_t i = 0; i < Taus; ++i){
            lagged_Pi[i] & stream;
        }
    }
    // else: all members are zero-initialised by the constructor, so don't do anything
}

void MolineauxInfection::Variant::operator& (ostream& stream) {
    bool nonZero =
            Pi1 != 0.0 ||
            Pi2 != 0.0 ||
            Si_summation != 0.0;

    nonZero & stream;
    if( nonZero ){
        Pi1 & stream;
        Pi2 & stream;
        Si_summation & stream;
        for(size_t i = 0; i < Taus; ++i){
            lagged_Pi[i] & stream;
        }
    }
}

}
}
