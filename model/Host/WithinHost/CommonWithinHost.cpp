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

#include "Host/Human.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/AgeGroupInterpolation.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "schema/scenario.h"

using namespace std;

namespace OM {
namespace WithinHost {

CommonInfection* (* CommonWithinHost::s_createInfection) (LocalRng& rng, uint32_t protID, int origin);
CommonInfection* (* CommonWithinHost::s_checkpointedInfection) (istream& stream);

double hetMassMultStdDev = std::numeric_limits<double>::signaling_NaN();
double minHetMassMult = std::numeric_limits<double>::signaling_NaN();
util::AgeGroupInterpolator massByAge;

bool reportInfectedOrPatentInfected = false;

// -----  Initialization  -----

void CommonWithinHost::init( const scnXml::Scenario& scenario ){
    const scnXml::Human human = scenario.getModel().getHuman();
    if( !human.getWeight().present() ){
        // Technically this is needed by the PK/PD and Molineaux models
        throw util::XmlScenarioError( "model->human->weight element required by certain models" );
    }
    massByAge.set( human.getWeight().get(), "weight" );
    hetMassMultStdDev = human.getWeight().get().getMultStdDev();
    // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
    minHetMassMult = 0.5 / massByAge.eval( 0.0 );
    
    reportInfectedOrPatentInfected = mon::isUsed(mon::measure("totalInfs")) ||
        mon::isUsed(mon::measure("totalPatentInf"));
    
    PkPd::LSTMModel::init( scenario );
}

CommonWithinHost::CommonWithinHost( LocalRng& rng, double comorbidityFactor ) :
        WHFalciparum( rng, comorbidityFactor )
{
    assert( Sim::oneTS() == Sim::fromDays(1) || Sim::oneTS() == Sim::fromDays(5) );
    
    m_optVaccineGenotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    // Sample a weight heterogeneity factor
#ifndef NDEBUG
    int counter = 0;
#endif
    do {
        m_hetMassMultiplier = rng.gauss( 1.0, hetMassMultStdDev );
#ifndef NDEBUG
        assert( counter < 100 );        // too many resamples: resamples should rarely be needed...
        ++counter;
#endif
    } while( m_hetMassMultiplier < minHetMassMult );
}

CommonWithinHost::~CommonWithinHost() {
    for( auto inf = m_infections.begin(); inf != m_infections.end(); ++inf ){
        delete *inf;
    }
    m_infections.clear();
}

// -----  Simple infection adders/removers  -----

void CommonWithinHost::clearInfections( Treatments::Stages stage ){
    for(auto inf = m_infections.begin(); inf != m_infections.end();) {
        if( stage == Treatments::BOTH ||
            (stage == Treatments::LIVER && !(*inf)->bloodStage()) ||
            (stage == Treatments::BLOOD && (*inf)->bloodStage())
        ){
            delete *inf;
            inf = m_infections.erase( inf );
        }else{
            ++inf;
        }
    }
    m_numInfs = m_infections.size();
}

// -----  interventions -----

void CommonWithinHost::treatPkPd(size_t schedule, size_t dosage, double age, double delayD){
    double mass = massByAge.eval( age ) * m_hetMassMultiplier;
    m_pkpdModel.prescribe( schedule, dosage, age, mass, delayD );
}
void CommonWithinHost::clearImmunity() {
    for(auto inf = m_infections.begin(); inf != m_infections.end(); ++inf) {
        (*inf)->clearImmunity();
    }
    m_cumulativeH = 0.0;
    m_cumulativeYLag = 0.0;
}
void CommonWithinHost::importInfection(LocalRng& rng){
    if( m_numInfs < MaxInfections ){
        m_cumulativeH += 1;
        m_numInfs += 1;
        // This is a hook, used by interventions. The newly imported infections
        // should use initial frequencies to select genotypes.
        vector<double> weights( 0 );        // zero length: signal to use initial frequencies
        uint32_t genotype = Genotypes::sampleGenotype(rng, weights);
        m_infections.push_back(s_createInfection(rng, genotype, WithinHost::InfectionOrigin::Imported));
    }
    assert( m_numInfs == static_cast<int>(m_infections.size()) );
}


// -----  Density calculations  -----

void CommonWithinHost::update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotypeWeightsI, vector<double>& genotypeWeightsL, double ageInYears)
{
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    nNewInfs_l = min(nNewInfs_l,MaxInfections-m_numInfs);
    nNewInfs_i = min(nNewInfs_i,MaxInfections-m_numInfs-nNewInfs_l);

    int nNewInfsIgnored = nNewInfs_i + nNewInfs_l - (MaxInfections-m_numInfs);

    assert( m_numInfs>=0 && m_numInfs<=MaxInfections );

    int nNewInfsDiscarded = 0;
    for( int i=0; i<nNewInfs_i; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotypeWeightsI);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( m_optVaccineGenotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                m_infections.push_back(s_createInfection (rng, genotype, InfectionOrigin::Introduced));
            else
                nNewInfsDiscarded++;
        }
        else if (m_optVaccineGenotype == false)
            m_infections.push_back(s_createInfection (rng, genotype, InfectionOrigin::Introduced));
    }
    // Update nNewInfs, this is the number that will be reported in Human
    nNewInfs_i -= nNewInfsDiscarded;
    m_numInfs += nNewInfs_i;

    nNewInfsDiscarded = 0;
    for( int i=0; i<nNewInfs_l; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotypeWeightsL);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( m_optVaccineGenotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                m_infections.push_back(s_createInfection (rng, genotype, InfectionOrigin::Indigenous));
            else
                nNewInfsDiscarded++;
        }
        else if (m_optVaccineGenotype == false)
            m_infections.push_back(s_createInfection (rng, genotype, InfectionOrigin::Indigenous));
    }
    // Update nNewInfs, this is the number that will be reported in Human
    nNewInfs_l -= nNewInfsDiscarded;
    m_numInfs += nNewInfs_l;

    assert( m_numInfs == static_cast<int>(m_infections.size()) );
    
    updateImmuneStatus ();

    m_totalDensity = 0.0;
    m_hrp2Density = 0.0;
    m_timeStepMaxDensity = 0.0;
    
    bool treatmentLiver = m_treatExpiryLiver > Sim::ts0();
    bool treatmentBlood = m_treatExpiryBlood > Sim::ts0();
    
    double bodyMass = massByAge.eval( ageInYears ) * m_hetMassMultiplier;
    
    for( SimTime now = Sim::ts0(), end = Sim::ts0() + Sim::oneTS(); now < end; now = now + Sim::oneDay() ){
        // every day, medicate drugs, update each infection, then decay drugs
        m_pkpdModel.medicate(rng);
        
        double sumLogDens = 0.0;
        
        for(auto inf = m_infections.begin(); inf != m_infections.end();) {
            // Note: this is only one treatment model; there is also the PK/PD model
            bool expires = ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver);
            
            if( !expires ){     /* no expiry due to simple treatment model; do update */
                const double drugFactor = m_pkpdModel.getDrugFactor(rng, *inf, bodyMass);
                const double immFactor = immunitySurvivalFactor(ageInYears, (*inf)->cumulativeExposureJ());
                const double bsvFactor = human.vaccine.getFactor(interventions::Vaccine::BSV, m_optVaccineGenotype? (*inf)->genotype() : 0);
                const double survivalFactor = bsvFactor * m_innateImmSurvFact * immFactor * drugFactor;
                // update, may result in termination of infection:
                expires = (*inf)->update(rng, survivalFactor, now, bodyMass);
            }
            
            if( expires ){
                delete *inf;
                inf = m_infections.erase(inf);        // inf points to next infection now so don't increment with ++inf
                --m_numInfs;
            } else {
                double density = (*inf)->getDensity();
                m_totalDensity += density;
                if( !(*inf)->isHrp2Deficient() ){
                    m_hrp2Density += density;
                }
                m_timeStepMaxDensity = max(m_timeStepMaxDensity, density);
                if( density > 0 ){
                    // Base 10 logarithms are usually used; +1 because it avoids negatives in output while having very little affect on high densities
                    sumLogDens += log10(1.0 + density);
                }
                ++inf;
            }
        }
        m_pkpdModel.decayDrugs (bodyMass);
    }
    
    // As in AJTMH p22, cumulative_h (X_h + 1) doesn't include infections added
    // this time-step and cumulative_Y only includes past densities, thus we
    // increment these after the update.
    m_cumulativeH += nNewInfs_i + nNewInfs_l;
    m_cumulativeY += m_totalDensity;
    
    util::streamValidate(m_totalDensity);
    util::streamValidate(m_hrp2Density);
    assert( (std::isfinite)(m_totalDensity) );        // inf probably wouldn't be a problem but NaN would be
    
    // Cache total density for infectiousness calculations
    int yLagI = Sim::moduloSteps(Sim::ts1(), s_yLagLen);
    for( size_t g = 0; g < Genotypes::N(); ++g )
    {
        m_yLagI[yLagI * Genotypes::N() + g] = 0.0;
        m_yLagL[yLagI * Genotypes::N() + g] = 0.0;
    }

    for( auto inf = m_infections.begin(); inf != m_infections.end(); ++inf )
    {
        if((*inf)->origin() == InfectionOrigin::Imported)
            m_yLagI[yLagI * Genotypes::N() + (*inf)->genotype()] += (*inf)->getDensity();
        else
            m_yLagL[yLagI * Genotypes::N() + (*inf)->genotype()] += (*inf)->getDensity();
    }

    // This is a bug, we keep it this way to be consistent with old simulations
    if(nNewInfsIgnored > 0)
        nNewInfs_l += nNewInfsIgnored;
}

void CommonWithinHost::addProphylacticEffects(const vector<double>& pClearanceByTime) {
    // this should actually be easy; it just isn't needed yet
    throw util::UnimplementedException( "prophylactic effects on 1-day time step" );
}

InfectionOrigin CommonWithinHost::getInfectionOrigin()const
{
    return WithinHost::getInfectionOrigin(m_infections);
}

// -----  Summarize  -----

// Used in summarizeInfs.
vector<CommonInfection*> sortedInfs;
struct InfGenotypeSorter {
    bool operator() (CommonInfection* i, CommonInfection* j){
        return i->genotype() < j->genotype();
    }
} infGenotypeSorter;

bool CommonWithinHost::summarize( Host::Human& human )const{
    m_pathogenesisModel->summarize( human );
    m_pkpdModel.summarize( human );
    
    InfectionOrigin infectionType = WithinHost::getInfectionOrigin(m_infections);
    
    // If the number of infections is 0 and parasite density is positive we default to Indigenous
    if( m_infections.size() > 0 ){
        mon::recordStat(mon::measure("nInfect"), human);
        if(infectionType == InfectionOrigin::Indigenous)
            mon::recordStat(mon::measure("nInfect_Indigenous"), human);
        else if(infectionType == InfectionOrigin::Introduced)
            mon::recordStat(mon::measure("nInfect_Introduced"), human);
        else
            mon::recordStat(mon::measure("nInfect_Imported"), human);

        if( reportInfectedOrPatentInfected ){
            for(auto inf = m_infections.begin(); inf != m_infections.end(); ++inf)
            {
                uint32_t genotype = (*inf)->genotype();
                mon::recordStat(mon::measure("totalInfs"), human, 1, 0, genotype, 0);
                if((*inf)->origin() == InfectionOrigin::Indigenous)
                    mon::recordStat(mon::measure("totalInfs_Indigenous"), human, 1, 0, genotype, 0);
                else if((*inf)->origin() == InfectionOrigin::Introduced)
                    mon::recordStat(mon::measure("totalInfs_Introduced"), human, 1, 0, genotype, 0);
                else
                    mon::recordStat(mon::measure("totalInfs_Imported"), human, 1, 0, genotype, 0);

                if( Diagnostics::monitoringDiagnostic().isPositive( human.rng, (*inf)->getDensity(), std::numeric_limits<double>::quiet_NaN() ) ){
                    mon::recordStat(mon::measure("totalPatentInf"), human, 1, 0, genotype, 0);
                    if((*inf)->origin() == InfectionOrigin::Indigenous)
                        mon::recordStat(mon::measure("totalPatentInf_Indigenous"), human, 1, 0, genotype, 0);
                    else if((*inf)->origin() == InfectionOrigin::Introduced)
                        mon::recordStat(mon::measure("totalPatentInf_Introduced"), human, 1, 0, genotype, 0);
                    else
                        mon::recordStat(mon::measure("totalPatentInf_Imported"), human, 1, 0, genotype, 0);
                }
            }
        }
        
        if( reportInfectionsByGenotype ){
            // Instead of storing nInfs and total density by genotype we sort
            // infections by genotype and report each in sequence.
            // We don't sort in place since that would affect random number sampling
            // order when updating, and the monitoring system should not in my
            // opinion affect outputs (since it would make testing harder).
            sortedInfs.assign( m_infections.begin(), m_infections.end() );
            sort( sortedInfs.begin(), sortedInfs.end(), infGenotypeSorter );
            auto inf = sortedInfs.begin();
            while( inf != sortedInfs.end() ){
                uint32_t genotype = (*inf)->genotype();
                double dens = 0.0;
                do{     // at start: genotype is that of the current infection, dens is 0
                    dens += (*inf)->getDensity();
                    ++inf;
                }while( inf != sortedInfs.end() && (*inf)->genotype() == genotype );
                // we had at least one infection of this genotype
                mon::recordStat(mon::measure("nInfectByGenotype"), human, 1, 0, genotype, 0);
                if( Diagnostics::monitoringDiagnostic().isPositive(human.rng, dens, std::numeric_limits<double>::quiet_NaN()) ){
                    mon::recordStat(mon::measure("nPatentByGenotype"), human, 1, 0, genotype, 0);
                    mon::recordStat(mon::measure("logDensByGenotype"), human, log(dens), 0, genotype, 0);
                }
            }
        }
    }
    
    // Some treatments (simpleTreat with steps=-1) clear infections immediately
    // (and are applied after update()), thus infections.size() may be 0 while
    // totalDensity > 0. Here we report the last calculated density.
    if( Diagnostics::monitoringDiagnostic().isPositive(human.rng, m_totalDensity, std::numeric_limits<double>::quiet_NaN()) ){
        mon::recordStat(mon::measure("nPatent"), human);
        if(infectionType == InfectionOrigin::Imported)
            mon::recordStat(mon::measure("nPatent_Imported"), human);
        else if(infectionType == InfectionOrigin::Introduced)
            mon::recordStat(mon::measure("nPatent_Introduced"), human);
        else if(infectionType == InfectionOrigin::Indigenous)
            mon::recordStat(mon::measure("nPatent_Indigenous"), human);

        if(m_totalDensity == 0.0)
            mon::recordStat(mon::measure("sumlogDens"), human, 0.0);
        else
            mon::recordStat(mon::measure("sumlogDens"), human, log(m_totalDensity));
        return true;    // patent
    }
    return false;       // not patent
}


void CommonWithinHost::checkpoint (istream& stream) {
    WHFalciparum::checkpoint (stream);
    m_hetMassMultiplier & stream;
    m_pkpdModel & stream;
    for(int i = 0; i < m_numInfs; ++i) {
        m_infections.push_back (s_checkpointedInfection (stream));
    }
    assert( m_numInfs == static_cast<int>(m_infections.size()) );
}

void CommonWithinHost::checkpoint (ostream& stream) {
    WHFalciparum::checkpoint (stream);
    m_hetMassMultiplier & stream;
    m_pkpdModel & stream;
    for(auto inf = m_infections.begin(); inf != m_infections.end(); ++inf) {
        (**inf) & stream;
    }
}
}
}
