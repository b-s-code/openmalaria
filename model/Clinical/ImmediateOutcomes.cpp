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

#include "Clinical/ImmediateOutcomes.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/random.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;

// These parameters are set by setHealthSystem() and do not need checkpointing.
double ImmediateOutcomes::s_cureRateUcOfficial[ImmediateOutcomes::NumCaseTypes];
double ImmediateOutcomes::s_cureRateUcSelfTreat[NumCaseTypes];
WithinHost::TreatmentId ImmediateOutcomes::s_treatmentUc[NumCaseTypes];
bool ImmediateOutcomes::s_useDiagnosticUc = false;


// ———  static, utility functions  ———

double getHealthSystemACRByName (const scnXml::TreatmentDetails& td, const string& drug)
{
    double val;
    if (drug == "CQ")
        val = td.getCQ().present() ? td.getCQ().get().getValue() : 0;
    else if (drug == "SP")
        val = td.getSP().present() ? td.getSP().get().getValue() : 0;
    else if (drug == "AQ")
        val = td.getAQ().present() ? td.getAQ().get().getValue() : 0;
    else if (drug == "SPAQ")
        val = td.getSPAQ().present() ? td.getSPAQ().get().getValue() : 0;
    else if (drug == "ACT")
        val = td.getACT().present() ? td.getACT().get().getValue() : 0;
    else if (drug == "QN")
        val = td.getQN().present() ? td.getQN().get().getValue() : 0;
    else if (drug == "selfTreatment")
        val = td.getSelfTreatment().getValue();
    else {
        throw util::XmlScenarioError ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
    }
    
    if( !(val >= 0.0 && val <= 1.0) ){
        throw util::XmlScenarioError (string("healthSystem initialACR/compliance/nonCompliersEffective: ").append(drug).append(" value must be in range [0,1]"));
    }
    return val;
}

// note: "CQOptional" is the same type as for other drugs
const scnXml::TreatmentActions::CQOptional& getHealthSystemTreatmentOptional(
    const scnXml::TreatmentActions& treatActs, const string& drug )
{
    if( drug == "CQ" ) return treatActs.getCQ();
    if( drug == "SP" ) return treatActs.getSP();
    if( drug == "AQ" ) return treatActs.getAQ();
    if( drug == "SPAQ" ) return treatActs.getSPAQ();
    if( drug == "ACT" ) return treatActs.getACT();
    if( drug == "QN" ) return treatActs.getQN();
    throw util::XmlScenarioError ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
}
WithinHost::TreatmentId getHealthSystemTreatmentByName( const scnXml::TreatmentActions& treatActs, const string& drug )
{
    const scnXml::TreatmentActions::CQOptional& elt = getHealthSystemTreatmentOptional( treatActs, drug );
    if( !elt.present() )
        throw util::XmlScenarioError( string("healthSystem.treatmentActions: description required for ").append(drug) );
    return WithinHost::WHInterface::addTreatment( elt.get() );
}


// ———  static, init  ———

void ImmediateOutcomes::setHealthSystem (const scnXml::HSImmediateOutcomes& hsDescription) {
    const string &firstLine = hsDescription.getDrugRegimen().getFirstLine(),
        &secondLine = hsDescription.getDrugRegimen().getSecondLine(),
        &inpatient = hsDescription.getDrugRegimen().getInpatient();
    
    // ———  calculate probability of getting treatment  ———

    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    s_accessUcSelfTreat[FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    s_accessUcSelfTreat[SecondLine] = 0.0;
    s_accessUcAny[FirstLine] = accessUCOfficial1 + s_accessUcSelfTreat[FirstLine];
    s_accessUcAny[SecondLine] = accessUCOfficial1 + s_accessUcSelfTreat[SecondLine];
    s_accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();
    
    if( accessUCOfficial1 < 0.0 || s_accessUcSelfTreat[FirstLine] < 0.0 || s_accessUcAny[FirstLine] > 1.0 ||
        accessUCOfficial2 < 0.0 || s_accessUcSelfTreat[SecondLine] < 0.0 || s_accessUcAny[SecondLine] > 1.0 ||
        s_accessSevere < 0.0 || s_accessSevere > 1.0 )
    {
        throw util::XmlScenarioError ("healthSystem: "
            "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
            " at least 0 and their sum must be at most 1, and"
            "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
            "be in range [0,1]");
    }
    
    // ———  calculate probability of clearing parasites  ———
    
    const double complianceFirstLine = getHealthSystemACRByName (hsDescription.getCompliance(), firstLine);
    const double complianceSecondLine = getHealthSystemACRByName (hsDescription.getCompliance(), secondLine);

    const double cureRateFirstLine = getHealthSystemACRByName (hsDescription.getInitialACR(), firstLine);
    const double cureRateSecondLine = getHealthSystemACRByName (hsDescription.getInitialACR(), secondLine);

    const double nonCompliersEffectiveFirstLine = getHealthSystemACRByName (hsDescription.getNonCompliersEffective(), firstLine);
    const double nonCompliersEffectiveSecondLine = getHealthSystemACRByName (hsDescription.getNonCompliersEffective(), secondLine);

    const double complianceSelfTreatment = hsDescription.getCompliance().getSelfTreatment().getValue();
    const double cureRateSelfTreatment = hsDescription.getInitialACR().getSelfTreatment().getValue();
    if( !(
        complianceSelfTreatment >= 0.0 && complianceSelfTreatment <= 1.0
        && cureRateSelfTreatment >= 0.0 && cureRateSelfTreatment <= 1.0
    ) ){
        throw util::XmlScenarioError ("healthSystem initialACR/compliance/nonCompliersEffective: pSelfTreatment must be in range [0,1]");
    }
    
    
    s_cureRateUcOfficial[FirstLine] = complianceFirstLine * cureRateFirstLine
            + (1 - complianceFirstLine) * nonCompliersEffectiveFirstLine;
    s_cureRateUcOfficial[SecondLine] = complianceSecondLine * cureRateSecondLine
            + (1 - complianceSecondLine) * nonCompliersEffectiveSecondLine;
    s_cureRateUcSelfTreat[FirstLine] = complianceSelfTreatment * cureRateSelfTreatment
            + (1 - complianceSelfTreatment) * nonCompliersEffectiveFirstLine;
    s_cureRateUcSelfTreat[SecondLine] = 0.0;      // as with access, this is a historical asymmetry
    
    s_cureRateSevere = getHealthSystemACRByName (hsDescription.getInitialACR(), inpatient);
    
    s_treatmentUc[FirstLine] = getHealthSystemTreatmentByName(hsDescription.getTreatmentActions(), firstLine);
    
    if( secondLine == firstLine ) s_treatmentUc[SecondLine] = s_treatmentUc[FirstLine];
    else s_treatmentUc[SecondLine] = getHealthSystemTreatmentByName(hsDescription.getTreatmentActions(), secondLine);
    
    if( inpatient == firstLine ) s_treatmentSevere = s_treatmentUc[FirstLine];
    else if( inpatient == secondLine ) s_treatmentSevere = s_treatmentUc[SecondLine];
    else s_treatmentSevere = getHealthSystemTreatmentByName(hsDescription.getTreatmentActions(), inpatient);
    
    s_useDiagnosticUc = hsDescription.getUseDiagnosticUC();
    
    if( ModelOptions::option(util::VIVAX_SIMPLE_MODEL) ){
        WithinHost::WHVivax::setHSParameters(
            hsDescription.getLiverStageDrug().present() ?
            &hsDescription.getLiverStageDrug().get() : 0 );
    }else if( hsDescription.getLiverStageDrug().present() ){
        throw util::XmlScenarioError( "health-system's liverStageDrug element only supported by vivax" );
    }
}


// ———  per-human, update  ———

void ImmediateOutcomes::uncomplicatedEvent (
    Human& human,
    Episode::State pgState
) {
    // If last treatment prescribed was in recent memory, consider second line.
    CaseType regimen = FirstLine;
    if (tLastTreatment + s_healthSystemMemory > Sim::ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = SecondLine;
    }

    m_latestReport.update (human, Episode::State( pgState ) );

    double x = human.rng.uniform01();
    if( x < s_accessUcAny[regimen] * m_treatmentSeekingFactor ){
        // UC1: official care OR self treatment
        // UC2: official care only
        
        if( s_useDiagnosticUc ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
            auto diag = WithinHost::Diagnostics::monitoringDiagnostic();
            if( !human.withinHostModel->diagnosticResult(human.rng, diag) )
                return; // negative outcome: no treatment
        }
        
        tLastTreatment = Sim::ts0();
        mon::recordEvent(s_measures[regimen], human);
        
        double p = ( x < s_accessUcSelfTreat[regimen] * m_treatmentSeekingFactor ) ?
            s_cureRateUcSelfTreat[regimen] : s_cureRateUcOfficial[regimen];
        if( human.rng.bernoulli(p) ){
            // Could report Episode::RECOVERY to latestReport,
            // but we don't report out-of-hospital recoveries anyway.
            human.withinHostModel->treatment( human, s_treatmentUc[regimen] );
        } else {
            // No change in parasitological status: treated outside of hospital
        }
        
        human.withinHostModel->optionalPqTreatment(human);
    } else {
        // No change in parasitological status: non-treated
    }
}

}
}
