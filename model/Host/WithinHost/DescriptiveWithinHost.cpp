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

#include "Global.h"
#include "Host/Human.h"
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "util/errors.h"

#include <cassert>

using namespace std;

namespace OM {
namespace WithinHost {

extern bool bugfixMaxDens;    // DescriptiveInfection.cpp
bool reportPatentInfected = false;
// -----  Initialization  -----

void DescriptiveWithinHostModel::initDescriptive(){
    reportPatentInfected = (
                                mon::isUsed(mon::measure("totalPatentInf"))
                                || mon::isUsed(mon::measure("totalPatentInf_Imported"))
                                || mon::isUsed(mon::measure("totalPatentInf_Introduced"))
                                || mon::isUsed(mon::measure("totalPatentInf_Indigenous"))
                            );
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel( LocalRng& rng, double comorbidityFactor ) :
        WHFalciparum( rng, comorbidityFactor )
{
    assert( Sim::oneTS() == Sim::fromDays(5) );
    m_optVaccineGenotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);
}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {
    for( auto inf = m_infections.begin(); inf != m_infections.end(); ++inf ){
        delete *inf;
    }
    m_infections.clear();
}


// -----  Simple infection adders/removers  -----

void DescriptiveWithinHostModel::loadInfection(istream& stream) {
    m_infections.push_back(new DescriptiveInfection(stream));
}

void DescriptiveWithinHostModel::clearInfections( Treatments::Stages stage ){
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

// -----  Interventions  -----

void DescriptiveWithinHostModel::clearImmunity() {
    for(auto inf = m_infections.begin(); inf != m_infections.end(); ++inf) {
        (*inf)->clearImmunity();
    }
    m_cumulativeH = 0.0;
    m_cumulativeYLag = 0.0;
}

void DescriptiveWithinHostModel::importInfection(LocalRng& rng){
    if( m_numInfs < MaxInfections ){
        m_cumulativeH += 1;
        m_numInfs += 1;
        // This is a hook, used by interventions. The newly imported infections
        // should use initial frequencies to select genotypes.
        vector<double> weights( 0 );        // zero length: signal to use initial frequencies
        uint32_t genotype = Genotypes::sampleGenotype(rng, weights);
        m_infections.push_back(new DescriptiveInfection(rng, genotype, WithinHost::InfectionOrigin::Imported));
    }
    assert( m_numInfs == static_cast<int>(m_infections.size()) );
}


// -----  Density calculations  -----

void DescriptiveWithinHostModel::update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotypeWeightsI, vector<double>& genotypeWeightsL, double ageInYears)
{
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    int nNewInfsToBeCreated_i = nNewInfs_i;
    int nNewInfsToBeCreated_l = nNewInfs_l;

    nNewInfs_l = min(nNewInfs_l,MaxInfections-m_numInfs);
    nNewInfs_i = min(nNewInfs_i,MaxInfections-m_numInfs-nNewInfs_l);

    m_numInfs += nNewInfs_i;
    assert( m_numInfs>=0 && m_numInfs<=MaxInfections );
    for( int i=0; i<nNewInfs_i; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotypeWeightsI);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( m_optVaccineGenotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                m_infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Introduced));
        }
        else if (m_optVaccineGenotype == false)
            m_infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Introduced));
    }
    assert( m_numInfs == static_cast<int>(m_infections.size()) );

    m_numInfs += nNewInfs_l;
    assert( m_numInfs>=0 && m_numInfs<=MaxInfections );
    for( int i=0; i<nNewInfs_l; ++i ) {
        uint32_t genotype = Genotypes::sampleGenotype(rng, genotypeWeightsL);

        // If opt_vaccine_genotype is true the infection is discarded with probability 1-vaccineFactor
        if( m_optVaccineGenotype )
        {
            double vaccineFactor = human.vaccine.getFactor( interventions::Vaccine::PEV, genotype );
            if(vaccineFactor == 1.0 || human.rng.bernoulli(vaccineFactor))
                m_infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Indigenous));
        }
        else if (m_optVaccineGenotype == false)
            m_infections.push_back(new DescriptiveInfection (rng, genotype, InfectionOrigin::Indigenous));
    }
    assert( m_numInfs == static_cast<int>(m_infections.size()) );

    updateImmuneStatus ();

    m_totalDensity = 0.0;
    m_hrp2Density = 0.0;
    m_timeStepMaxDensity = 0.0;

    bool treatmentLiver = m_treatExpiryLiver > Sim::ts0();
    bool treatmentBlood = m_treatExpiryBlood > Sim::ts0();

    for(auto inf = m_infections.begin(); inf != m_infections.end();) {
        //NOTE: it would be nice to combine this code with that in
        // CommonWithinHost.cpp, but a few changes would be needed:
        // INNATE_MAX_DENS and MAX_DENS_CORRECTION would need to be required
        // (couldn't support old parameterisations using buggy versions of code
        // any more).
        // SP drug action and the PK/PD model would need to be abstracted
        // behind a common interface.
        if ( (*inf)->expired() /* infection has self-terminated */ ||
            ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver) )
        {
            delete *inf;
            inf=m_infections.erase(inf);
            m_numInfs--;
            continue;
        }
        
        // Should be: infStepMaxDens = 0.0, but has some history.
        // See MAX_DENS_CORRECTION in DescriptiveInfection.cpp.
        double infStepMaxDens = m_timeStepMaxDensity;
        double immSurvFact = immunitySurvivalFactor(ageInYears, (*inf)->cumulativeExposureJ());
        double bsvFactor = human.vaccine.getFactor(interventions::Vaccine::BSV, m_optVaccineGenotype? (*inf)->genotype() : 0);

        (*inf)->determineDensities(rng, m_cumulativeH, infStepMaxDens, immSurvFact, m_innateImmSurvFact, bsvFactor);

        if (bugfixMaxDens)
            infStepMaxDens = std::max(infStepMaxDens, m_timeStepMaxDensity);
        m_timeStepMaxDensity = infStepMaxDens;

        double density = (*inf)->getDensity();
        m_totalDensity += density;
        if( !(*inf)->isHrp2Deficient() ){
            m_hrp2Density += density;
        }

        ++inf;
    }
    
    // As in AJTMH p22, cumulative_h (X_h + 1) doesn't include infections added
    // this time-step and cumulative_Y only includes past densities.
    m_cumulativeH += nNewInfs_i + nNewInfs_l;
    m_cumulativeY += Sim::oneTS() * m_totalDensity;
    
    util::streamValidate( m_totalDensity );
    util::streamValidate( m_hrp2Density );
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
    if(m_optVaccineGenotype == false)
    {
        nNewInfs_i = nNewInfsToBeCreated_i;
        nNewInfs_l = nNewInfsToBeCreated_l;
    }
}

InfectionOrigin DescriptiveWithinHostModel::getInfectionOrigin()const
{
    return WithinHost::getInfectionOrigin(m_infections);
}

// -----  Summarize  -----

bool DescriptiveWithinHostModel::summarize( Host::Human& human )const{
    m_pathogenesisModel->summarize( human );
    
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

        int nImported = 0, nIntroduced = 0, nIndigenous = 0;
        for( auto inf = m_infections.begin(); inf != m_infections.end(); ++inf )
        {
            if((*inf)->origin() == InfectionOrigin::Indigenous) nIndigenous++;
            else if((*inf)->origin() == InfectionOrigin::Introduced) nIntroduced++;
            else nImported++;
        }

        // (patent) infections are reported by genotype, even though we don't have
        // genotype in this model
        mon::recordStat(mon::measure("totalInfs"), human, static_cast<int>(m_infections.size()));
        mon::recordStat(mon::measure("totalInfs_Imported"), human, nImported);
        mon::recordStat(mon::measure("totalInfs_Introduced"), human, nIntroduced);
        mon::recordStat(mon::measure("totalInfs_Indigenous"), human, nIndigenous);

        if( reportPatentInfected ){
            for(auto inf = m_infections.begin(); inf != m_infections.end(); ++inf)
            {
                if( Diagnostics::monitoringDiagnostic().isPositive( human.rng, (*inf)->getDensity(), std::numeric_limits<double>::quiet_NaN() ) )
                {
                    mon::recordStat(mon::measure("totalPatentInf"), human);
                    if((*inf)->origin() == InfectionOrigin::Indigenous)
                        mon::recordStat(mon::measure("totalPatentInf_Indigenous"), human);
                    else if((*inf)->origin() == InfectionOrigin::Introduced)
                        mon::recordStat(mon::measure("totalPatentInf_Introduced"), human);
                    else
                        mon::recordStat(mon::measure("totalPatentInf_Imported"), human);
                }
            }
        }
        if( reportInfectionsByGenotype ){
            // accumulate total density by genotype
            map<uint32_t, double> densByGtype;
            for(auto inf = m_infections.begin(); inf != m_infections.end(); ++inf)
                densByGtype[(*inf)->genotype()] += (*inf)->getDensity();
            
            for( auto gtype: densByGtype ){
                // we had at least one infection of this genotype
                mon::recordStat(mon::measure("nInfectByGenotype"), human, 1, 0, gtype.first, 0);
                if( Diagnostics::monitoringDiagnostic().isPositive(human.rng, gtype.second, std::numeric_limits<double>::quiet_NaN()) ){
                    mon::recordStat(mon::measure("nPatentByGenotype"), human, 1, 0, gtype.first, 0);
                    mon::recordStat(mon::measure("logDensByGenotype"), human, log(gtype.second), 0, gtype.first, 0);
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

        if(m_totalDensity > 1e-10)
            mon::recordStat(mon::measure("sumlogDens"), human, log(m_totalDensity));
        return true;    // patent
    }
    return false;       // not patent
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::checkpoint (istream& stream) {
    WHFalciparum::checkpoint (stream);
    for(int i=0; i<m_numInfs; ++i) {
        loadInfection(stream);  // create infections using a virtual function call
    }
    assert( m_numInfs == static_cast<int>(m_infections.size()) );
}
void DescriptiveWithinHostModel::checkpoint (ostream& stream) {
    WHFalciparum::checkpoint (stream);
    for(DescriptiveInfection* inf : m_infections)
        *inf & stream;
}

char const*const NotImpl = "feature not available with the \"descriptive\" within-host model";
void DescriptiveWithinHostModel::treatPkPd(size_t schedule, size_t dosages, double age, double delayD){
    throw TRACED_EXCEPTION( NotImpl, util::Error::WHFeatures ); }

}
}
