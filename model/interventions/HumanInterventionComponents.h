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

// Note: this file is included by exactly one source (Interventions.cpp)
// and contains definitions as well as declarations.

// The includes here are more for documentation than required.
#include "interventions/Interfaces.h"
#include "interventions/InterventionManager.h"
#include "util/ModelOptions.h"
#include "Clinical/ESCaseManagement.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Host/Human.h"
#include "Transmission/TransmissionModel.h"
#include "Host/WithinHost/Diagnostic.h"
#include "PkPd/LSTMTreatments.h"
#include "mon/Monitoring.h"
#include "util/random.h"
#include <schema/healthSystem.h>
#include <schema/interventions.h>

namespace OM { namespace interventions {
using Host::Human;
using WithinHost::Diagnostic;
using WithinHost::Diagnostics;

// ———  HumanIntervention  ———

bool componentCmp(const HumanInterventionComponent *a, const HumanInterventionComponent *b){
    return a->componentType() < b->componentType();
}
HumanIntervention::HumanIntervention(
    const xsd::cxx::tree::sequence<scnXml::Component>& componentList )
{
    m_components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it ) {
        ComponentId id = InterventionManager::getComponentId( it->getId() );
        const HumanInterventionComponent* component = &InterventionManager::getComponent( id );
        m_components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( m_components.begin(), m_components.end(), componentCmp );
}
HumanIntervention::HumanIntervention(
    const xsd::cxx::tree::sequence<scnXml::Component>& componentList,
    const xsd::cxx::tree::sequence<scnXml::Condition>& conditionList )
{
    m_components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it ) {
        ComponentId id = InterventionManager::getComponentId( it->getId() );
        const HumanInterventionComponent* component = &InterventionManager::getComponent( id );
        m_components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( m_components.begin(), m_components.end(), componentCmp );
    
    for( auto it = conditionList.begin(), end = conditionList.end(); it != end; ++it ) {
        double min = it->getMinValue().present() ? it->getMinValue().get() : -numeric_limits<double>::infinity();
        double max = it->getMaxValue().present() ? it->getMaxValue().get() : numeric_limits<double>::infinity();
        m_conditions.push_back(mon::setupCondition( it->getMeasure(), min, max, it->getInitialState() ));
    }
}
HumanIntervention::HumanIntervention( const xsd::cxx::tree::sequence<scnXml::DTDeploy>& componentList ){
    m_components.reserve( componentList.size() );
    for( auto it = componentList.begin(), end = componentList.end(); it != end; ++it ) {
        ComponentId id = InterventionManager::getComponentId( it->getComponent() );
        const HumanInterventionComponent* component = &InterventionManager::getComponent( id );
        m_components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( m_components.begin(), m_components.end(), componentCmp );
}

void HumanIntervention::deploy( Human& human, mon::Deploy::Method method,
    VaccineLimits vaccLimits ) const
{
    for( size_t cond : m_conditions ){
        // Abort deployment if any condition is false.
        if( !mon::checkCondition(cond) ) return;
    }
    
    for( auto it = m_components.begin(); it != m_components.end(); ++it ) {
        const interventions::HumanInterventionComponent& component = **it;
        // we must report first, since it can change cohort and sub-population
        // which may affect what deployment does (at least in the case of reporting deployments)
        human.addToCohort(component.id(), component.duration());
        component.deploy( human, method, vaccLimits );
    }
}

void HumanIntervention::printDetails( std::ostream& out )const{
    out << "human:";
    for( auto it = m_components.begin(); it != m_components.end(); ++it ){
        out << '\t' << (*it)->id().id;
    }
}

// ———  Utilities  ———

TriggeredDeployments::TriggeredDeployments(const scnXml::TriggeredDeployments& elt){
    m_lists.reserve( elt.getDeploy().size() );
    for( auto it = elt.getDeploy().begin(), end = elt.getDeploy().end(); it != end; ++it ) {
        m_lists.push_back( SubList( *it ) );
    }
}
void TriggeredDeployments::deploy(Human& human,
        mon::Deploy::Method method, VaccineLimits vaccLimits) const
{
    for( auto it = m_lists.begin(); it != m_lists.end(); ++it ){
        it->deploy( human, method, vaccLimits );
    }
}

TriggeredDeployments::SubList::SubList( const scnXml::TriggeredDeployments::DeployType& elt ) :
        HumanIntervention( elt.getComponent() ),
        minAge( Sim::fromYearsN( elt.getMinAge() ) ),
        maxAge( Sim::future() ),
        coverage( elt.getP() )
{
    if( elt.getMaxAge().present() )
        maxAge = Sim::fromYearsN( elt.getMaxAge().get() );
    
    if( minAge < Sim::zero() || maxAge < minAge ){
        throw util::XmlScenarioError("triggered intervention must have 0 <= minAge <= maxAge");
    }
    
    if( coverage < 0.0 || coverage > 1.0 ){
        throw util::XmlScenarioError( "triggered intervention must have 0 <= coverage <= 1" );
    }
    
    // Zero coverage: optimise
    if( coverage <= 0.0 || minAge >= maxAge ) m_components.clear();
}
void TriggeredDeployments::SubList::deploy( Host::Human& human,
        mon::Deploy::Method method, VaccineLimits vaccLimits )const
{
    // This may be used within a time step; in that case we use human age at
    // the beginning of the step since this gives the expected result with age
    // limits of 0 to 1 time step.
    //TODO: why does the above contradict what we do? Is this due to a date being shifted later?
    SimTime age = human.age(Sim::nowOrTs1());
    if( age >= minAge && age < maxAge ){
        if( coverage >= 1.0 || human.rng.bernoulli( coverage ) ){
            HumanIntervention::deploy( human, method, vaccLimits );
        }
    }
}


// ———  Derivatives of HumanInterventionComponent  ———

class RecruitmentOnlyComponent : public HumanInterventionComponent {
public:
    RecruitmentOnlyComponent( ComponentId id ) :
        HumanInterventionComponent(id)
    {}
    virtual ~RecruitmentOnlyComponent() {}
    
    /// Reports to monitoring, nothing else
    virtual void deploy( Host::Human& human, mon::Deploy::Method method, VaccineLimits ) const{
        mon::recordDeploy(mon::measure("nMassRecruitOnly"), mon::measure("nCtsRecruitOnly"), human, method);
    }
    
    virtual Component::Type componentType() const{
        return Component::RECRUIT_ONLY;
    }
    
    virtual void printDetails( std::ostream& out )const{
        out << id().id << "\tRecruit only";
    }
};

class ScreenComponent : public HumanInterventionComponent {
public:
    ScreenComponent( ComponentId id, const scnXml::Screen& elt ) :
        HumanInterventionComponent(id),
        m_diagnostic(Diagnostics::get(elt.getDiagnostic())),
        m_positive( elt.getPositive() ),
        m_negative( elt.getNegative() )
    {}
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits vaccLimits ) const{
        mon::recordDeploy(mon::measure("nMassScreenings"), mon::measure("nCtsScreenings"), human, method);
        if( human.withinHostModel->diagnosticResult(human.rng, m_diagnostic) ){
            m_positive.deploy( human, method, vaccLimits );
        }else{
            m_negative.deploy( human, method, vaccLimits );
        }
    }
    
    virtual Component::Type componentType() const{ return Component::SCREEN; }
    
    virtual void printDetails( std::ostream& out )const{
        out << id().id << "\tScreen";
    }
    
private:
    const Diagnostic& m_diagnostic;
    HumanIntervention m_positive, m_negative;
};

/// Simple treatment: no PK/PD, just remove parasites
class TreatSimpleComponent : public HumanInterventionComponent {
public:
    TreatSimpleComponent( ComponentId id, const scnXml::DTTreatSimple& elt ) :
        HumanInterventionComponent(id)
    {
        //NOTE: this code is currently identical to that in CMDTTreatSimple
        try{
            SimTime durL = UnitParse::readShortDuration( elt.getDurationLiver(), UnitParse::NONE ),
                durB = UnitParse::readShortDuration( elt.getDurationBlood(), UnitParse::NONE );
            SimTime neg1 = -Sim::oneTS();
            if( durL < neg1 || durB < neg1 ){
                throw util::XmlScenarioError( "treatSimple: cannot have durationBlood or durationLiver less than -1" );
            }
            if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
                if( durL != Sim::zero() || durB != neg1 )
                    throw util::UnimplementedException( "vivax model only supports timestepsLiver=0, timestepsBlood=-1" );
                // Actually, the model ignores these parameters; we just don't want somebody thinking it doesn't.
            }
            m_timeLiver = durL;
            m_timeBlood = durB;
        }catch( const util::FormatError& e ){
            throw util::XmlScenarioError( string("treatSimple: ").append(e.message()) );
        }
    }
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits ) const{
        mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        human.withinHostModel->treatSimple( human, m_timeLiver, m_timeBlood );
    }
    
    virtual Component::Type componentType() const{ return Component::TREAT_SIMPLE; }
    
    virtual void printDetails( std::ostream& out )const{
        out << id().id << "\ttreatSimple";
    }
    
private:
    SimTime m_timeLiver, m_timeBlood;
};

//NOTE: this is similar to CMDTTreatPKPD, but supporting only a single "treatment"
class TreatPKPDComponent : public HumanInterventionComponent {
public:
    TreatPKPDComponent( ComponentId id, const scnXml::DTTreatPKPD& elt ) :
        HumanInterventionComponent(id),
            m_schedule(PkPd::LSTMTreatments::findSchedule(elt.getSchedule())),
            m_dosage(PkPd::LSTMTreatments::findDosages(elt.getDosage())),
            m_delayD(elt.getDelay_h() / 24.0 /* convert to days */)
    {
    }
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits ) const{
        mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        double age = Sim::inYears(human.age(Sim::nowOrTs1()));
        human.withinHostModel->treatPkPd( m_schedule, m_dosage, age, m_delayD );
    }
    
    virtual Component::Type componentType() const{ return Component::TREAT_PKPD; }
    
    virtual void printDetails( std::ostream& out )const{
        out << id().id << "\ttreatPKPD";
    }
    
private:
    size_t m_schedule;        // index of the schedule
    size_t m_dosage;          // name of the dosage table
    double m_delayD;         // delay in days
};

class DecisionTreeComponent : public HumanInterventionComponent {
public:
    DecisionTreeComponent( ComponentId id, const scnXml::DecisionTree& elt ) :
        HumanInterventionComponent(id),
        m_tree(Clinical::CMDecisionTree::create(elt, false))
    {}
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits )const{
        Clinical::CMDTOut out = m_tree.exec( Clinical::CMHostData(human,
                Sim::inYears(human.age(Sim::nowOrTs1())),
                human.clinicalModel->getLatestState()) );
        if( out.treated ){
            mon::recordDeploy(mon::measure("nMDAs"), mon::measure("nCtsMDA"), human, method);
        }
        if( out.screened ){
            mon::recordDeploy(mon::measure("nMassScreenings"), mon::measure("nCtsScreenings"), human, method);
        }
    }
    
    virtual Component::Type componentType() const{ return Component::CM_DT; }
    
    virtual void printDetails( std::ostream& out )const{
        out << id().id << "\tDecision tree";
    }
    
private:
    const Clinical::CMDecisionTree& m_tree;
};

class ClearImmunityComponent : public HumanInterventionComponent {
public:
    ClearImmunityComponent( ComponentId id ) :
        HumanInterventionComponent(id) {}
    
    void deploy( Human& human, mon::Deploy::Method method, VaccineLimits )const{
        human.withinHostModel->clearImmunity();
        //NOTE: not reported
    }
    
    virtual Component::Type componentType() const{ return Component::CLEAR_IMMUNITY; }
    
    virtual void printDetails( std::ostream& out )const{
        out << id().id << "\tclear immunity";
    }
};

} }
