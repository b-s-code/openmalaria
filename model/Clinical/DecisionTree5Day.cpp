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

#include "Clinical/DecisionTree5Day.h"
#include "Clinical/CMDecisionTree.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/errors.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;

enum CaseType { FirstLine, SecondLine, NumCaseTypes };

// These parameters are set by setHealthSystem() and do not need checkpointing.
const CMDecisionTree *treeUCOfficial = 0, *treeUCSelfTreat = 0;

// ———  static, init  ———

void DecisionTree5Day::setHealthSystem(const scnXml::HSDT5Day& hsDescription){
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
    
    treeUCOfficial = &CMDecisionTree::create( hsDescription.getTreeUCOfficial(), true );
    treeUCSelfTreat = &CMDecisionTree::create( hsDescription.getTreeUCSelfTreat(), true );
    
    s_cureRateSevere = hsDescription.getCureRateSevere().getValue();
    s_treatmentSevere = WHInterface::addTreatment( hsDescription.getTreatmentSevere() );
    
    if( ModelOptions::option(util::VIVAX_SIMPLE_MODEL) ){
        WithinHost::WHVivax::setHSParameters(
            hsDescription.getLiverStageDrug().present() ?
            &hsDescription.getLiverStageDrug().get() : 0 );
    }else if( hsDescription.getLiverStageDrug().present() ){
        throw util::XmlScenarioError( "health-system's liverStageDrug element only supported by vivax" );
    }
}


// ———  per-human, update  ———

void DecisionTree5Day::uncomplicatedEvent ( Human& human, Episode::State pgState ){    
    // If last treatment prescribed was in recent memory, consider second line.
    CaseType regimen = FirstLine;
    if (tLastTreatment + s_healthSystemMemory > Sim::ts0()){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = SecondLine;
    }

    m_latestReport.update (human, Episode::State( pgState ) );

    double x = human.rng.uniform01();
    if( x < s_accessUcAny[regimen] * m_treatmentSeekingFactor ){
        CMHostData hostData( human, Sim::inYears(human.age(Sim::ts0())), pgState );
        
        // Run tree (which may deploy treatment)
        CMDTOut output = ( x < s_accessUcSelfTreat[regimen] * m_treatmentSeekingFactor ) ?
            treeUCSelfTreat->exec( hostData ) :
            treeUCOfficial->exec( hostData );
        
        if( output.treated ){   // if any treatment or intervention deployed
            tLastTreatment = Sim::ts0();
            mon::recordEvent(s_measures[regimen], human);
        }
        if( output.screened ){
            mon::recordEvent(mon::measure("nTreatDiagnostics"), human);
        }
        
        human.withinHostModel->optionalPqTreatment(human);
    } else {
        // No care sought
    }
}

}
}
