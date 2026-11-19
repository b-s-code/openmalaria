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

#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/WHVivax.h"
#include "Host/WithinHost/DescriptiveWithinHost.h"
#include "Host/WithinHost/CommonWithinHost.h"
#include "Host/WithinHost/Infection/DummyInfection.h"
#include "Host/WithinHost/Infection/EmpiricalInfection.h"
#include "Host/WithinHost/Infection/MolineauxInfection.h"
#include "Host/WithinHost/Infection/PennyInfection.h"
#include "Host/WithinHost/Treatments.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "schema/scenario.h"
//using namespace std;

#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

bool reportInfectionsByGenotype = false;

bool optVivaxSimple = false,
        optDummyWhm = false, optEmpiricalWhm = false,
        optMolineauxWhm = false, optPennyWhm = false,
        optCommonWhm = false;

// -----  static functions  -----

void WHInterface::init( const OM::Parameters& parameters, const scnXml::Scenario& scenario ) {
    reportInfectionsByGenotype = mon::isUsed(mon::measure("nInfectByGenotype")) ||
        mon::isUsed(mon::measure("nPatentByGenotype")) ||
        mon::isUsed(mon::measure("logDensByGenotype"));
    
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        optVivaxSimple = true;
        WHVivax::init( parameters, scenario.getModel() );
    }else{
        WHFalciparum::init( parameters, scenario.getModel() );
        
        if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL)) {
            optDummyWhm = true;
            DummyInfection::init ();
        } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
            optEmpiricalWhm = true;
            EmpiricalInfection::init();    // 1-day time step check
        } else if (util::ModelOptions::option (util::MOLINEAUX_WITHIN_HOST_MODEL)) {
            optMolineauxWhm = true;
            MolineauxInfection::init( parameters );
        } else if (util::ModelOptions::option (util::PENNY_WITHIN_HOST_MODEL)) {
            optPennyWhm = true;
            PennyInfection::init();
        } else {
            DescriptiveInfection::init( parameters );      // 5-day time step check
        }
        optCommonWhm = optDummyWhm || optEmpiricalWhm
                || optMolineauxWhm || optPennyWhm;
        
        if( optCommonWhm ){
            CommonWithinHost::init( scenario );
        }else{
            DescriptiveWithinHostModel::initDescriptive();
        }
    }
}

TreatmentId WHInterface::addTreatment(const scnXml::TreatmentOption& desc){
    return Treatments::addTreatment( desc );
}

unique_ptr<WHInterface> WHInterface::createWithinHostModel(LocalRng& rng, double comorbidityFactor) {
    if( optVivaxSimple ) {
        return unique_ptr<WHInterface>(new WHVivax( rng, comorbidityFactor ));
    } else if( optCommonWhm ) {
        return unique_ptr<WHInterface>(new CommonWithinHost( rng, comorbidityFactor ));
    } else {
        return unique_ptr<WHInterface>(new DescriptiveWithinHostModel( rng, comorbidityFactor ));
    }
}


// -----  Non-static  -----


void WHInterface::checkpoint (istream& stream) {
    m_numInfs & stream;

    if (m_numInfs > MaxInfections)
        throw util::CheckpointError("numInfs: " + to_string(m_numInfs));
}
void WHInterface::checkpoint (ostream& stream) {
    m_numInfs & stream;
}

}
}
