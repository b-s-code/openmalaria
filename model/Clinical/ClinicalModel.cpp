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

#include "Clinical/ClinicalModel.h"

#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/DecisionTree5Day.h"
#include "Host/NeonatalMortality.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"

namespace OM { namespace Clinical {

bool optEventScheduler = false;
bool optImmOutcomes = false;

bool ClinicalModel::s_indirectMortBugfix;
SimTime ClinicalModel::s_healthSystemMemory{ Sim::never() };

//log odds ratio of case-fatality in community compared to hospital
double oddsRatioThreshold;

util::AgeGroupInterpolator ClinicalModel::s_caseFatalityRate;
util::AgeGroupInterpolator ClinicalModel::s_pSequelaeInpatient;

// -----  static methods  -----

void ClinicalModel::init( const Parameters& parameters, const scnXml::Scenario& scenario ) {
    const scnXml::Clinical& clinical = scenario.getModel().getClinical();
    try{
        s_indirectMortBugfix = util::ModelOptions::option (util::INDIRECT_MORTALITY_FIX);
        //NOTE: if changing XSD, this should not have a default unit:
        s_healthSystemMemory = UnitParse::readShortDuration(clinical.getHealthSystemMemory(), UnitParse::STEPS);
        oddsRatioThreshold = exp( parameters[Parameter::LOG_ODDS_RATIO_CF_COMMUNITY] );
        InfantMortality::init( parameters );
    }catch( const util::FormatError& e ){
        throw util::XmlScenarioError( string("model/clinical/healthSystemMemory: ").append(e.message()) );
    }
    
    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER)){
        optEventScheduler = true;
        ClinicalEventScheduler::init( parameters, clinical );
    }else{
        if( scenario.getHealthSystem().getImmediateOutcomes().present() ){
            optImmOutcomes = true;
            
            if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                cerr << "Deprecation warning: healthSystem: use of ImmediateOutcomes can be replaced by the more flexible DecisionTree5Day (optional)" << endl;
            }
        }
        // else: decision tree 5 day
        
        CM5DayCommon::init();
    }
}

void ClinicalModel::setHS( const scnXml::HealthSystem& healthSystem ){
    s_caseFatalityRate.set( healthSystem.getCFR(), "CFR" );
    s_pSequelaeInpatient.set( healthSystem.getPSequelaeInpatient(), "pSequelaeInpatient" );
    if( optEventScheduler ){
        if( !healthSystem.getEventScheduler().present() ){
            throw util::XmlScenarioError ("Expected EventScheduler "
                "section in healthSystem data (initial or intervention)");
        }
        ESCaseManagement::setHealthSystem(healthSystem.getEventScheduler().get());
    }else if( optImmOutcomes ){
        if ( !healthSystem.getImmediateOutcomes().present() ){
            throw util::XmlScenarioError ("Expected ImmediateOutcomes "
                "section in healthSystem data (initial or intervention)");
        }
        ImmediateOutcomes::setHealthSystem(healthSystem.getImmediateOutcomes().get());
    }else{
        if( !healthSystem.getDecisionTree5Day().present() ){
            throw util::XmlScenarioError ("Expected DecisionTree5Day "
                "section in healthSystem data (initial or intervention)");
        }
        DecisionTree5Day::setHealthSystem(healthSystem.getDecisionTree5Day().get());
    }
}

unique_ptr<ClinicalModel> ClinicalModel::createClinicalModel (double tSF) {
    if (optEventScheduler){
        return unique_ptr<ClinicalModel>(new ClinicalEventScheduler (tSF));
    }else if( optImmOutcomes ){
        return unique_ptr<ClinicalModel>(new ImmediateOutcomes (tSF));
    }else{
        return unique_ptr<ClinicalModel>(new DecisionTree5Day( tSF ));
    }
}

double ClinicalModel::getCommunityCFR (double caseFatalityRatio){
    double x = caseFatalityRatio * oddsRatioThreshold;
    return x / (1 - caseFatalityRatio + x);
}



// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel () :
    m_doomed(NOT_DOOMED)
{}
ClinicalModel::~ClinicalModel () {
  // latestReport is reported, if any, by destructor
}


// -----  other non-static methods  -----

bool ClinicalModel::isDead( SimTime age ){
    if( age >= Sim::maxHumanAge())       // too old (reached age limit)
        m_doomed = DOOMED_TOO_OLD;
    return m_doomed > NOT_DOOMED;	// killed by some means
}

void ClinicalModel::update (Human& human, double ageYears, bool newBorn) {
    if (m_doomed < NOT_DOOMED)	// Countdown to indirect mortality
        m_doomed -= Sim::oneTS();
    
    //indirect death: if this human's about to die, don't worry about further episodes:
    if (m_doomed <= DOOMED_EXPIRED) {	//clinical bout 6 intervals before
        mon::recordEvent(mon::measure("nIndDeaths"), human);
        m_doomed = DOOMED_INDIRECT;
        return;
    }
    if(newBorn /* i.e. first update since birth */) {
        // Chance of neonatal mortality:
        if (Host::NeonatalMortality::eventNeonatalMortality(human.rng)) {
            mon::recordEvent(mon::measure("nIndDeaths"), human);
            m_doomed = DOOMED_NEONATAL;
            return;
        }
    }
    
    const bool isDoomed = m_doomed != NOT_DOOMED;
    WithinHost::Pathogenesis::StatePair pg = human.withinHostModel->determineMorbidity( human, ageYears, isDoomed );
    m_latestState = static_cast<Episode::State>(pg.state);

    doClinicalUpdate (human, ageYears, pg);
}

void ClinicalModel::updateInfantDeaths( SimTime age ){
    // update array for the infant death rates
    if (age < Sim::oneYear()){
        size_t index = age / Sim::oneTS();
        
        // Testing doomed == DOOMED_NEXT_TS gives very slightly different results than
        // testing doomed == DOOMED_INDIRECT (due to above if(..))
        bool isDoomed = m_doomed == DOOMED_COMPLICATED
                || m_doomed == DOOMED_NEXT_TS
                || m_doomed == DOOMED_NEONATAL;
        InfantMortality::reportRisk( index, isDoomed );
    }
}


void ClinicalModel::checkpoint (istream& stream) {
    m_latestReport & stream;
    m_doomed & stream;
}
void ClinicalModel::checkpoint (ostream& stream) {
    m_latestReport & stream;
    m_doomed & stream;
}



// ——— infant mortality ———

// Infant death summaries (checkpointed).
vector<int> infantDeaths;
vector<int> infantIntervalsAtRisk;

/// Non-malaria mortality in under 1year olds.
/// Set by init ()
double nonMalariaMortality;

void InfantMortality::init( const OM::Parameters& parameters ){
    infantDeaths.resize(Sim::stepsPerYear());
    infantIntervalsAtRisk.resize(Sim::stepsPerYear());
    nonMalariaMortality=parameters[Parameter::NON_MALARIA_INFANT_MORTALITY];
}

void InfantMortality::preMainSimInit() {
    std::fill(infantDeaths.begin(), infantDeaths.end(), 0);
    std::fill(infantIntervalsAtRisk.begin(), infantIntervalsAtRisk.end(), 0);
}

void InfantMortality::staticCheckpoint(istream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}
void InfantMortality::staticCheckpoint (ostream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}

void InfantMortality::reportRisk(size_t index, bool isDoomed) {
    infantIntervalsAtRisk[index] += 1;     // baseline
    if (isDoomed)
        infantDeaths[index] += 1;  // deaths
}

double InfantMortality::allCause(){
    double infantPropSurviving=1.0;       // use to calculate proportion surviving
    for( size_t i = 0; i < Sim::stepsPerYear(); i += 1 ){
        // multiply by proportion of infants surviving at each interval
        infantPropSurviving *= double(infantIntervalsAtRisk[i] - infantDeaths[i])
            / double(infantIntervalsAtRisk[i]);
    }
    // Child deaths due to malaria (per 1000), plus non-malaria child deaths. Deaths per 1000 births is the return unit.
    return (1.0 - infantPropSurviving) * 1000.0 + nonMalariaMortality;
}


} }
