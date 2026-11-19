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
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"
#include "mon/init.h"
#include "mon/Monitoring.h"

namespace OM {

void Sim::init( const scnXml::Scenario& scenario, util::ModelNameProvider mnp ){
    initInterval(scenario, mnp);
    SimData::s_stepsPerYear = Sim::inSteps(Sim::oneYear());
    SimData::s_yearsPerStep = 1.0 / SimData::s_stepsPerYear;
    Sim::s_maxHumanAge = Sim::fromYearsD( scenario.getDemography().getMaximumAgeYrs() );
    
    Sim::s_start = Sim::origin();
    auto mon = scenario.getMonitoring();
    if( mon.getStartDate().present() ){
        try{
            // on failure, this throws or returns Sim::never()
            Sim::s_start = UnitParse::parseDate( mon.getStartDate().get() );
            if( Sim::s_start == Sim::never() ){
                throw util::FormatError( "invalid format (expected YYYY-MM-DD)" );
            }
        }catch( const util::FormatError& e ){
            throw util::XmlScenarioError( string("monitoring/startDate: ").append(e.message()) );
        }
    }
    
    Sim::s_interv = Sim::never(); // large negative number
    Sim::s_end = mon::readSurveyDates( mon );
}

void Sim::initInterval(const scnXml::Scenario& scenario, util::ModelNameProvider mnp)
{
    // If scenario states to use a named model, set the relevant time step duration for
    // that named model here.  It may possibly get overwritten later in this method if
    // the scenario also states an explicit interval value - a behaviour that exists so
    // that users have the ability to override the time step duration if they need to.
    const util::ModelNames namedModelToUse = mnp.getModelName();
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel)
    {
        if (namedModelToUse == util::ModelNames::base)
        {
            // This is where we define the length of the default time step, in days,
            // for the base model.
            SimData::s_interval = 5;
        }
        else
        {
            throw util::XmlScenarioError(
                "No pre-set interval value is available for the specified model name.");
        }
    }

    // Apply the interval duration explicitly written in the scenario, if one exists.
    const bool parameterBlockExistsInScenario = scenario.getModel().getParameters().present();
    if (parameterBlockExistsInScenario)
    {
        // Schema dictates that interval attribute is required if parameter element exists.
        // So no need to check for existsence of interval attribute in this case.
        SimData::s_interval = scenario.getModel().getParameters().get().getInterval();
    }
}

}
