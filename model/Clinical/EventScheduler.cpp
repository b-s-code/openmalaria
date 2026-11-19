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

#include "Clinical/EventScheduler.h"
#include "util/random.h"
#include "Host/WithinHost/WHInterface.h"
#include "mon/Monitoring.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include <schema/scenario.h>

#include <limits>

namespace OM { namespace Clinical {
    using namespace OM::util;

SimTime ClinicalEventScheduler::s_maxUcSeekingMemory(Sim::never());
SimTime ClinicalEventScheduler::s_uncomplicatedCaseDuration(Sim::never());
SimTime ClinicalEventScheduler::s_complicatedCaseDuration(Sim::never());
SimTime ClinicalEventScheduler::s_extraDaysAtRisk(Sim::never());
vector<double> ClinicalEventScheduler::s_cumDailyPrImmUcts;
double ClinicalEventScheduler::s_negV;
double ClinicalEventScheduler::s_alpha;

double ClinicalEventScheduler::s_logOddsAbBase = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::s_logOddsAbNegTest = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::s_logOddsAbPosTest = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::s_logOddsAbNeed = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::s_logOddsAbInformal = std::numeric_limits<double>::signaling_NaN();
double ClinicalEventScheduler::s_oneMinusEfficacyAb = std::numeric_limits<double>::signaling_NaN();
AgeGroupInterpolator ClinicalEventScheduler::s_severeNmfMortality;

bool optNonMalariaFevers = false;

AgeGroupInterpolator ClinicalEventScheduler::s_nmfNeedAntibiotic;
AgeGroupInterpolator ClinicalEventScheduler::s_mfNeedAntibiotic;

// -----  static init  -----

void ClinicalEventScheduler::init( const Parameters& parameters, const scnXml::Clinical& clinical )
{
    if( Sim::oneTS() != Sim::oneDay() )
        throw util::XmlScenarioError ("ClinicalEventScheduler is only designed for a 1-day time step.");
    
    optNonMalariaFevers = util::ModelOptions::option( util::NON_MALARIA_FEVERS );
    
    s_alpha = exp( -parameters[Parameter::CFR_NEG_LOG_ALPHA] );
    if( !(0.0<=s_alpha && s_alpha<=1.0) ){
        throw util::XmlScenarioError(
            "Clinical outcomes: propDeathsFirstDay should be within range [0,1]");
    }
    s_negV = -parameters[Parameter::CFR_SCALE_FACTOR];
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::XmlScenarioError("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc2 =
            clinical.getNonMalariaFevers().get();
        if( !nmfDesc2.getPrNeedTreatmentMF().present() ||
            !nmfDesc2.getPrNeedTreatmentNMF().present() ){
            throw util::XmlScenarioError( "prNeedTreatmentMF and "
                "prNeedTreatmentNMF elements required in model->clinical->"
                "NonMalariaFevers" );
        }
        s_nmfNeedAntibiotic.set( nmfDesc2.getPrNeedTreatmentNMF().get(), "prNeedTreatmentNMF" );
        s_mfNeedAntibiotic.set( nmfDesc2.getPrNeedTreatmentMF().get(), "prNeedTreatmentMF" );
    }
}

void ClinicalEventScheduler::setParameters(const scnXml::HSEventScheduler& esData) {
    const scnXml::ClinicalOutcomes& coData = esData.getClinicalOutcomes();
    
    s_maxUcSeekingMemory = Sim::fromDays(coData.getMaxUCSeekingMemory());
    s_uncomplicatedCaseDuration = Sim::fromDays(coData.getUncomplicatedCaseDuration());
    s_complicatedCaseDuration = Sim::fromDays(coData.getComplicatedCaseDuration());
    s_extraDaysAtRisk = Sim::fromDays(coData.getComplicatedRiskDuration()) - s_complicatedCaseDuration;
    if( s_uncomplicatedCaseDuration < Sim::fromDays(1)
	|| s_complicatedCaseDuration < Sim::fromDays(1)
	|| s_maxUcSeekingMemory < Sim::zero()
	|| s_extraDaysAtRisk + s_complicatedCaseDuration < Sim::fromDays(1) // at risk at least 1 day
	|| s_extraDaysAtRisk > Sim::zero()        // at risk longer than case duration
    ){
	throw util::XmlScenarioError(
            "Clinical outcomes: constraints on case/risk/memory duration not met (see documentation)");
    }
    
    s_cumDailyPrImmUcts.reserve( coData.getDailyPrImmUCTS().size() );
    double cumP = 0.0;
    for( auto it = coData.getDailyPrImmUCTS().begin(); it != coData.getDailyPrImmUCTS().end(); ++it ){
        cumP += *it;
        s_cumDailyPrImmUcts.push_back( cumP );
    }
    if( cumP < 0.999 || cumP > 1.001 ){
        throw util::XmlScenarioError( "Event scheduler: dailyPrImmUCTS seq must add up to 1" );
    }
    s_cumDailyPrImmUcts.back() = 1.0;
    
    s_caseFatalityRate.scale( s_alpha );
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !esData.getNonMalariaFevers().present() ){
            throw util::XmlScenarioError("NonMalariaFevers element of healthSystem->EventScheduler required");
        }
        const scnXml::HSESNMF& nmfDesc = esData.getNonMalariaFevers().get();
        
        double pT = nmfDesc.getPrTreatment();
        s_logOddsAbBase = log( pT / (1.0 - pT) );
        s_logOddsAbNegTest = log(nmfDesc.getEffectNegativeTest());
        s_logOddsAbPosTest = log(nmfDesc.getEffectPositiveTest());
        s_logOddsAbNeed = log(nmfDesc.getEffectNeed());
        s_logOddsAbInformal = log(nmfDesc.getEffectInformal());
        s_oneMinusEfficacyAb = 1.0 - nmfDesc.getTreatmentEfficacy();
        s_severeNmfMortality.set( nmfDesc.getCFR(), "CFR" );
    }
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double tSF) :
        m_pgState (Episode::NONE),
        m_caseStartTime (Sim::never()),
        m_timeOfRecovery (Sim::never()),
        m_timeLastTreatment (Sim::never()),
        m_previousDensity (numeric_limits<double>::quiet_NaN())
{
    if( tSF != 1.0 ){
	// p(treatment seeking) is part of tree & handled generically, so we
	// don't have a way of modifying it.
	throw XmlScenarioError("treatment seeking heterogeneity not supported");
    }
}


// -----  other methods  -----

bool ClinicalEventScheduler::isExistingCase() {
        // If treated in the recent past:
        return Sim::now() > m_timeLastTreatment && Sim::now() <= m_timeLastTreatment + s_healthSystemMemory;
}

void ClinicalEventScheduler::doClinicalUpdate (Human& human, double ageYears, WithinHost::Pathogenesis::StatePair &pg){
    WHInterface& withinHostModel = *human.withinHostModel;
    Episode::State newState = static_cast<Episode::State>( pg.state );
    util::streamValidate( (newState << 16) & m_pgState );
    
    if ( Sim::ts0() == m_timeOfRecovery ) {
	if( m_pgState & Episode::DIRECT_DEATH ){
	    // Human dies this time step (last day of risk of death)
	    m_doomed = DOOMED_COMPLICATED;
	    
	    m_latestReport.update (human, m_pgState);
        } else {
	    if ( m_pgState & Episode::COMPLICATED ) {
                const double pSequelae = s_pSequelaeInpatient.eval( ageYears );
                mon::recordStat(mon::measure("expectedSequelae"), human, pSequelae);
		if( human.rng.uniform01() < pSequelae ){
		    m_pgState = Episode::State (m_pgState | Episode::SEQUELAE);
                }else{
		    m_pgState = Episode::State (m_pgState | Episode::RECOVERY);
                }
	    } else {
		m_pgState = Episode::State (m_pgState | Episode::RECOVERY);
            }
	    // report bout, at conclusion of episode:
	    m_latestReport.update (human, m_pgState);
	    
	    // Individual recovers (and is immediately susceptible to new cases)
	    m_pgState = Episode::NONE;	// recovery (reset to healthy state)
	    
	    // And returns to transmission (if was removed)
	    human.perHostTransmission.outsideTransmission = false;
	}
    }
    
    // The bug fixed here is that indirectMortality was not set if no event occurred
    bool indirectMortality = false;
    if( s_indirectMortBugfix ) indirectMortality = pg.indirectMortality;
    
    if ( newState & Episode::SICK ){
        // we have some new case: is it severe/complicated?
        if ( newState & Episode::COMPLICATED ){
            if ( m_pgState & Episode::COMPLICATED ) {
                // previously severe: no events happen for course of medication
            } else {
                // previously healthy or UC: progress to severe
                m_pgState = Episode::State (m_pgState | newState | Episode::RUN_CM_TREE);
                indirectMortality = pg.indirectMortality;
                m_caseStartTime = Sim::ts0();
            }
        } else {
            // uncomplicated case (UC/UC2/NMF): is it new?
            if (m_pgState & Episode::SICK) {
                // previously UC; nothing to do
            }else{
                if( m_caseStartTime < Sim::ts0() ) {
                    // new UC case
                    m_pgState = Episode::State (m_pgState | newState | Episode::RUN_CM_TREE);
                    indirectMortality = pg.indirectMortality;
                    
                    double uVariate = human.rng.uniform01();
                    size_t i = 0;       // units: days
                    for(; i < s_cumDailyPrImmUcts.size(); ++i){
                        if( uVariate < s_cumDailyPrImmUcts[i] ){
                            goto gotDelay;
                        }
                    }
                    assert(false);      // should have uVariate < 1 = cumDailyPrImmUCTS[len-1]
                    gotDelay:
                    // set start time: current time plus length of delay (days)
                    m_caseStartTime = Sim::ts0() + Sim::fromDays(i);
                }
            }
        }
        
        if (indirectMortality && m_doomed == NOT_DOOMED)
            m_doomed = -Sim::oneTS(); // start indirect mortality countdown
    }
    
    if( m_caseStartTime == Sim::ts0() && (m_pgState & Episode::RUN_CM_TREE) ){
        // OK, we're about to run the CM tree
        m_pgState = Episode::State (m_pgState & ~Episode::RUN_CM_TREE);
        
        // If last treatment prescribed was in recent memory, consider second line.
        if( m_timeLastTreatment + s_healthSystemMemory > Sim::ts0() ){
            m_pgState = Episode::State (m_pgState | Episode::SECOND_CASE);
        }
	
	CMDTOut auxOut = ESCaseManagement::execute(
	    CMHostData( human, ageYears, m_pgState ) );
	
        if( auxOut.treated ){	// I.E. some treatment was given
            m_timeLastTreatment = Sim::ts0();
            if( m_pgState & Episode::COMPLICATED ){
                mon::recordEvent(mon::measure("nTreatments3"), human);
            }else{
                if( m_pgState & Episode::SECOND_CASE ){
                    mon::recordEvent(mon::measure("nTreatments2"), human);
                }else{
                    mon::recordEvent(mon::measure("nTreatments1"), human);
                }
            }
        }
        if( auxOut.screened ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
        }
	
	if ( true /*FIXME auxOut.hospitalisation != CMAuxOutput::NONE*/ ) {	// in hospital
	    m_pgState = Episode::State (m_pgState | Episode::EVENT_IN_HOSPITAL);
	    
	    /*FIXME if( auxOut.hospitalisation == CMAuxOutput::DELAYED )
		++caseStartTime;*/
	}
	
	// Case fatality rate (first day of illness)
	// P(death) is some fixed input scaled by age-specific CFR.
	if( (m_pgState & Episode::COMPLICATED)
	    && !(m_pgState & Episode::DIRECT_DEATH)
	) {
            const bool inHospital = true /*FIXME auxOut.hospitalisation == CMAuxOutput::IMMEDIATE */;
	    double pDeath = s_caseFatalityRate.eval( ageYears );
	    // community fatality rate when not in hospital or delayed hospital entry
            if( !inHospital )
                pDeath = getCommunityCFR( pDeath );
            mon::recordStat(mon::measure("expectedDirectDeaths"), human, pDeath);
            if( inHospital )
                mon::recordStat(mon::measure("expectedHospitalDeaths"), human, pDeath);
	    if (human.rng.uniform01() < pDeath) {
		m_pgState = Episode::State (m_pgState | Episode::DIRECT_DEATH | Episode::EVENT_FIRST_DAY);
		// Human is killed at end of time at risk
		//timeOfRecovery += extraDaysAtRisk;	(no point updating; will be set later: ATORWD)
	    }
	    m_previousDensity = withinHostModel.getTotalDensity();
	}
	
	if( optNonMalariaFevers ){
            if( (m_pgState & Episode::SICK) && !(m_pgState & Episode::COMPLICATED) ){
                // Have a NMF or UC malaria case
                
                /** Given a non-malaria fever, return the probability of it requiring
                 * treatment. */
                bool isMalarial = m_pgState & Episode::MALARIA;
                double pNeedTreat = isMalarial ?
                        s_mfNeedAntibiotic.eval( ageYears ) :
                        s_nmfNeedAntibiotic.eval( ageYears );
                bool needTreat = human.rng.bernoulli(pNeedTreat);
                
                // Calculate chance of antibiotic administration:
                double pTreatment = 0.0;
                /*FIXME
                if( auxOut.AB_provider != CMAuxOutput::NO_AB ){
                    // Treatment seeking; may be at a health facility or from
                    // the informal sector.
                    double logOddsAb = logOddsAbBase - logOddsAbNeed * pNeedTreat;
                    if( auxOut.AB_provider == CMAuxOutput::FACILITY ){
                        if( auxOut.diagnostic == CMAuxOutput::NEGATIVE ){
                            logOddsAb += logOddsAbNegTest;
                        }else if( auxOut.diagnostic == CMAuxOutput::POSITIVE ){
                            logOddsAb += logOddsAbPosTest;
                        }
                        if( needTreat ){
                            logOddsAb += logOddsAbNeed;
                        }
                    }else{
                        assert( auxOut.AB_provider == CMAuxOutput::INFORMAL );
                        logOddsAb += logOddsAbInformal;
                    }
                    
                    double oddsTreatment = exp( logOddsAb );
                    pTreatment = oddsTreatment / (1.0 + oddsTreatment);
                }*/
                
                double treatmentEffectMult = 1.0;
                
                if( human.rng.uniform01() < pTreatment ){
                    /*FIXME: impossible due to above; NMF output removed
                    Survey::current().addInt( Report::MI_NMF_TREATMENTS, human, 1 );
                    treatmentEffectMult = oneMinusEfficacyAb;
                    */
                }
                
                // In a severe NMF case (only when not malarial), there is a
                // chance of death:
                if( needTreat ){
                    double pDeath = s_severeNmfMortality.eval( ageYears ) * treatmentEffectMult;
                    if( human.rng.uniform01() < pDeath ){
                        m_pgState = Episode::State (m_pgState | Episode::DIRECT_DEATH);
                    }
                }
            }
        }
    } else {
	// No new event (haven't changed state this time step).
	
	// Case fatality rate (subsequent days)
	// Complicated case & at risk of death (note: extraDaysAtRisk <= 0)
	if( (m_pgState & Episode::COMPLICATED)
	    && !(m_pgState & Episode::DIRECT_DEATH)
	    && (Sim::ts0() < m_timeOfRecovery + s_extraDaysAtRisk)
	) {
	    // In complicated episodes, S(t), the probability of survival on
	    // subsequent days t, is described by log(S(t)) = -v(Y(t)/Y(t-1)),
	    // for parasite density Y(t). v_neg below is -v.
            // TODO: this model should be revisited, and if possible placed
            // within the WHFalciparum class (or a subclass).
	    if( withinHostModel.getTotalDensity() > 0.0 ){	// avoid division by zero
		double parasiteReductionEffect = withinHostModel.getTotalDensity() / m_previousDensity;
		double pDeath = 1.0 - exp( s_negV * parasiteReductionEffect );
		// community fatality rate when not in hospital
		if( !(m_pgState & Episode::EVENT_IN_HOSPITAL) )
		    pDeath = getCommunityCFR( pDeath );
                mon::recordStat(mon::measure("expectedDirectDeaths"), human, pDeath);
                if( m_pgState & Episode::EVENT_IN_HOSPITAL )
                    mon::recordStat(mon::measure("expectedHospitalDeaths"), human, pDeath);
		if (human.rng.uniform01() < pDeath) {
		    m_pgState = Episode::State (m_pgState | Episode::DIRECT_DEATH);
		    // Human is killed at end of time at risk
		    m_timeOfRecovery = m_timeOfRecovery + s_extraDaysAtRisk;	// may be re-set later (see ATORWD)
		}
	    }
	    m_previousDensity = withinHostModel.getTotalDensity();
	}
    }
    
    // Start of case. Not necessarily start of sickness due to treatment-seeking
    // delays and travel time.
    if( m_caseStartTime == Sim::ts0() ){
	// Patients in hospital are removed from the transmission cycle.
	// This should have an effect from the start of the next time step.
	// NOTE: This is not very accurate, but considered of little importance.
	if (m_pgState & Episode::EVENT_IN_HOSPITAL)
	    human.perHostTransmission.outsideTransmission = true;
	
	if (m_pgState & Episode::COMPLICATED) {
	    // complicatedCaseDuration should to some respects be associated
	    // with medication duration, however ongoing medications after
	    // exiting hospital are OK and medications terminating before the
	    // end of hospitalisation shouldn't matter too much if the person
	    // can't recieve new infections due to zero transmission in hospital.
	    m_timeOfRecovery = Sim::ts0() + s_complicatedCaseDuration;
	    // Time should be adjusted to end of at-risk period when patient dies:
	    if( m_pgState & Episode::DIRECT_DEATH )	// death may already have been determined
		m_timeOfRecovery = m_timeOfRecovery + s_extraDaysAtRisk;	// ATORWD (search keyword)
	} else {
	    m_timeOfRecovery = Sim::ts0() + s_uncomplicatedCaseDuration;
	}
    }
    
    // Remove on first models...
    if( m_timeLastTreatment == Sim::ts0() )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_TREATMENT);
    
    if( m_pgState & Episode::SICK )
        human.removeFirstEvent(interventions::SubPopRemove::ON_FIRST_BOUT);
}


void ClinicalEventScheduler::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    int s;
    s & stream;
    m_pgState = Episode::State(s);
    m_caseStartTime & stream;
    m_timeOfRecovery & stream;
    m_timeLastTreatment & stream;
    m_previousDensity & stream;
}
void ClinicalEventScheduler::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    m_pgState & stream;
    m_caseStartTime & stream;
    m_timeOfRecovery & stream;
    m_timeLastTreatment & stream;
    m_previousDensity & stream;
}

} }
