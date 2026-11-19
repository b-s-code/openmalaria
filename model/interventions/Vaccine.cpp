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

#include "interventions/Vaccine.h"
#include "interventions/HumanComponents.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "schema/interventions.h"
#include "util/StreamValidator.h"
#include "Host/WithinHost/Genotypes.h"

#include <limits>
#include <cmath>

namespace OM {
namespace interventions {
using namespace OM::util;

vector<VaccineComponent*> VaccineComponent::s_params;
ComponentId VaccineComponent::s_reportComponent = ComponentId::wholePop();

VaccineComponent::VaccineComponent( ComponentId component, const scnXml::VaccineDescription& vd, Vaccine::Types type ) :
        HumanInterventionComponent(component),
        m_type(type),
        m_decayFunc(DecayFunction::makeObject( vd.getDecay(), "decay" )),
        m_efficacyB(vd.getEfficacyB().getValue())
{
    if( type == Vaccine::BSV && ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
        throw util::UnimplementedException( "blood stage vaccines (BSV) cannot be used with vivax model" );
    
    m_optVaccineGenotype = util::ModelOptions::option (util::VACCINE_GENOTYPE);

    if( s_reportComponent == ComponentId::wholePop() /* the initial value */ ){
        // set to the first type described
        s_reportComponent = component;
    }

    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd.getInitialEfficacy();
    const scnXml::VaccineDescription::PhenotypeSequence ps = vd.getPhenotype();

    //cout << ies.size() << " " << ps.size() << " " << WithinHost::Genotypes::N() << endl;

    if(m_optVaccineGenotype == false && ies.size() == 0)
        throw util::XmlScenarioError( "Vaccine: initialEfficacy is required but is missing" );
    
    if(m_optVaccineGenotype == true && ies.size() != 0)
        throw util::XmlScenarioError( "Vaccine: initialEfficacy outside of phenotype definition with VACCINE_GENOTYPE option enabled is not allowed" );
    
/*    if(opt_vaccine_genotype == true && ies.size() == 0 && ps.size() != WithinHost::Genotypes::N())
        throw util::XmlScenarioError( "Vaccine: Phenotype sequence must reference all genotypes" );
    */
    if(m_optVaccineGenotype == true)
    {
        m_perGenotypeInitialMeanEfficacy.resize(WithinHost::Genotypes::N());

        for(size_t i = 0; i < ps.size(); ++i)
        {
            //cout << "Phenotype: " << ps[i].getName() << endl;
            
            const scnXml::VaccinePhenotype::RestrictionSequence restrictions = ps[i].getRestriction();
            const scnXml::VaccinePhenotype::InitialEfficacySequence vpies = ps[i].getInitialEfficacy();

            for(size_t j = 0; j < restrictions.size(); ++j)
            {
                //cout << "\tonLocus: " << restrictions[j].getOnLocus() << ", " << "toAllelle: " << restrictions[j].getToAllele() << endl;
                //cout << "\tmatching allelle code: " << WithinHost::Genotypes::findAlleleCode(restrictions[j].getOnLocus(), restrictions[j].getToAllele()) << endl;

                uint32_t alleleCode = WithinHost::Genotypes::findAlleleCode(restrictions[j].getOnLocus(), restrictions[j].getToAllele());

                if( alleleCode == numeric_limits<uint32_t>::max() ){
                    throw util::XmlScenarioError("phenotype has "
                    "restriction on locus " + string(restrictions[j].getOnLocus())
                    + " allele " + string(restrictions[j].getToAllele()) + " but this locus/allele "
                    "has not been defined in parasiteGenetics section");
                }

                if (alleleCode >= WithinHost::Genotypes::N())
                    throw util::XmlScenarioError( string("Vaccine phenotype: Invalid Locus ") + restrictions[j].getOnLocus() + string(" and allele ") + restrictions[j].getToAllele());

                if (m_perGenotypeInitialMeanEfficacy[alleleCode].size() != 0)
                    throw util::XmlScenarioError( string("Vaccine phenotype: there are multiple restrictions locus " + string(restrictions[j].getOnLocus() + string(" and allele ") + restrictions[j].getToAllele())));

                if(vpies.size() == 0)
                    throw util::XmlScenarioError( string("Vaccine phenotype: initialEfficacy is required but is missing for phenotype: " + string(restrictions[j].getOnLocus() + string(" ") + restrictions[j].getToAllele())));

                m_perGenotypeInitialMeanEfficacy[alleleCode].resize(vpies.size());

                for(size_t k = 0; k < vpies.size(); ++k)
                    m_perGenotypeInitialMeanEfficacy[alleleCode][k] = vpies[k].getValue();
            }
        }

        // Check that all genotypes have been adressed
        for(size_t i = 0; i < WithinHost::Genotypes::N(); ++i)
        {
            if(m_perGenotypeInitialMeanEfficacy[i].size() == 0)
                throw util::XmlScenarioError( "Vaccine phenotype: missing restriction for at least one locus and allele pair" );
        }
    }
    else
    {
        m_perGenotypeInitialMeanEfficacy.resize(1);
        m_perGenotypeInitialMeanEfficacy[0].resize(ies.size());
        for(size_t i = 0; i < ies.size(); ++i)
            m_perGenotypeInitialMeanEfficacy[0][i] = ies[i].getValue();
    }

   /* initialMeanEfficacy.resize (ies.size());
    for(size_t i = 0; i < initialMeanEfficacy.size(); ++i)
        initialMeanEfficacy[i] = ies[i].getValue();*/
    
    if( s_params.size() <= component.id ) s_params.resize( component.id + 1 );
    assert( s_params[component.id] == 0 );
    s_params[component.id] = this;
}

void VaccineComponent::deploy(Host::Human& human, mon::Deploy::Method method, VaccineLimits vaccLimits) const
{
    bool administered = human.vaccine.possiblyVaccinate( human, id(), vaccLimits );
    if( administered && VaccineComponent::s_reportComponent == id() ){
        mon::recordDeploy(mon::measure("nMassVaccinations"), mon::measure("nEPIVaccinations"), human, method);
    }
    assert(m_type == Vaccine::PEV || m_type == Vaccine::BSV || m_type == Vaccine::TBV);
    if (method == mon::Deploy::TREAT)
        mon::recordEvent(mon::measure("nTreatDeployments"), human);
}

Component::Type VaccineComponent::componentType() const
{
    if( m_type == Vaccine::PEV ) return Component::PEV;
    else if( m_type == Vaccine::BSV ) return Component::BSV;
    else if( m_type == Vaccine::TBV ) return Component::TBV;
    else throw SWITCH_DEFAULT_EXCEPTION;
}

void VaccineComponent::printDetails(ostream& out) const
{
        out << id().id << "\t" << (m_type == Vaccine::PEV ? "PEV" : (m_type == Vaccine::BSV ? "BSV" : "TBV"));
}

double VaccineComponent::getInitialEfficacy (LocalRng& rng, size_t numPrevDoses, uint32_t genotype) const
{
    /* If initialMeanEfficacy.size or more doses have already been given, use
     * the last efficacy. */
    if (numPrevDoses >= m_perGenotypeInitialMeanEfficacy[genotype].size())
        numPrevDoses = m_perGenotypeInitialMeanEfficacy[genotype].size() - 1;
    double ime = m_perGenotypeInitialMeanEfficacy[genotype][numPrevDoses];
    //NOTE(validation): With extra valiadation in random, the first difference is noticed here:
    util::streamValidate(ime);
    util::streamValidate(m_efficacyB);
    if (ime == 0.0){
        return 0.0;
    } else if (ime < 1.0) {
        double result = rng.betaWithMean (ime, m_efficacyB);
        //NOTE(validation):: Without extra validation in random, the first difference is noticed here:
        util::streamValidate(result);
        //TODO(validation):: Why the difference? Bug in/limitation of StreamValidatior?
        // Extra memory allocation due to the extra logging causes some bad
        // memory usage to manifest differently? Am I forgetting to initialise
        // some variable?
        return result;
    } else {
        return 1.0;
    }
}

#if 0
/* R_0: this used to use _the_ vaccine configuration; this is just a check.
 * If reimplementing, use a new vaccine instance; implement some hook to ensure
 * that TBV is given to everyone after enabling the intervention, etc. */
void Vaccine::verifyEnabledForR_0 (){
    if( !types[PEV].active || !types[TBV].active )
        throw util::XmlScenarioError("PEV and TBV vaccines must have a "
                "description to use the insertR_0Case intervention");
}
#endif

// this is only used for checkpointing; loaded values are unimportant
PerEffectPerHumanVaccine::PerEffectPerHumanVaccine() :
    m_component(ComponentId::wholePop()),
    m_numDosesAdministered(0)
{
}

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine( LocalRng& rng, ComponentId id, const VaccineComponent& params ) :
    m_component( id ), m_numDosesAdministered(0)
{
    m_hetSample = params.m_decayFunc->hetSample(rng);
}

void PerEffectPerHumanVaccine::checkpointHet(ostream& stream) {
    VaccineComponent::getParams(m_component).m_decayFunc->checkpointSample(m_hetSample, stream);
}

void PerEffectPerHumanVaccine::checkpointHet(istream& stream) {
    VaccineComponent::getParams(m_component).m_decayFunc->checkpointSample(m_hetSample, stream);
}

double PerHumanVaccine::getFactor( Vaccine::Types type, uint32_t genotype) const {
    double factor = 1.0;
    for( EffectList::const_iterator effect = m_effects.begin(); effect != m_effects.end(); ++effect ){
        if( VaccineComponent::getParams(effect->m_component).m_type == type ){
            SimTime age = Sim::ts1() - effect->m_timeLastDeployment;  // implies age 1 TS on first use
            double decayFactor = effect->m_hetSample->eval( age );
            factor *= 1.0 - effect->m_perGenotypeInitialEfficacy[genotype] * decayFactor;
        }
    }
    return factor;
}

bool PerHumanVaccine::possiblyVaccinate( Host::Human& human,
        ComponentId componentId, interventions::VaccineLimits vaccLimits )
{
    PerEffectPerHumanVaccine* effect = 0;
    for( EffectList::iterator it = m_effects.begin(); it != m_effects.end(); ++it ){
        if( it->m_component == componentId ){
            effect = &*it;
            break;
        }
    }
    
    uint32_t numDosesAdministered = (effect == 0) ? 0 : effect->m_numDosesAdministered;
    if( numDosesAdministered < vaccLimits.minPrevDoses ||
        numDosesAdministered >= vaccLimits.maxCumDoses )
        return false;   // no vaccination (this replaces the old schedule for continuous doses)
    
    const VaccineComponent& params = VaccineComponent::getParams( componentId );
    
    if( effect == 0 ){
        m_effects.push_back( PerEffectPerHumanVaccine( human.rng, componentId, params ) );
        effect = &m_effects.back();
    }
    
    effect->m_perGenotypeInitialEfficacy.resize(WithinHost::Genotypes::N());

    for(size_t g=0; g<params.m_perGenotypeInitialMeanEfficacy.size(); g++)
        effect->m_perGenotypeInitialEfficacy[g] = params.getInitialEfficacy(human.rng, numDosesAdministered, g);

    util::streamValidate(effect->m_perGenotypeInitialEfficacy);
    
    effect->m_numDosesAdministered = numDosesAdministered + 1;
    effect->m_timeLastDeployment = Sim::nowOrTs1();
    
    return true;
}

}
}
