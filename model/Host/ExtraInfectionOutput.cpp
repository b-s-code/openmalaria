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

#include "Host/ExtraInfectionOutput.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "Clinical/ClinicalModel.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "util/errors.h"
#include "schema/scenario.h"

#include <fstream>
#include <vector>
#include <iostream>

namespace OM { namespace Host { namespace ExtraInfectionOutput {

using namespace std;

namespace {
    ofstream outFile;
    bool enabled = false;

    /// True when a startDate and/or endDate was given, so calendar-date
    /// filtering is applied.
    bool dateFilter = false;
    /// Inclusive calendar-date bounds (absolute SimTime, comparable to
    /// sim::intervDate()). Only meaningful when dateFilter is true.
    SimTime startBound = sim::never();
    SimTime endBound = sim::future();

    /// Map the human's current clinical (pathogenesis) state to a small integer
    /// code for output: 0 = STATE_NONE, 1 = STATE_MALARIA, 2 = STATE_SEVERE.
    int clinicalStatusCode( Clinical::Episode::State state ){
        if( state & Clinical::Episode::SEVERE ) return 2;   // STATE_SEVERE
        if( state & Clinical::Episode::MALARIA ) return 1;  // STATE_MALARIA
        return 0;                                           // STATE_NONE
    }

    /// Derive the CSV file name from the scenario file name: strip any trailing
    /// ".xml" and append ".csv".
    string deriveOutputName( const string& scenarioFileName ){
        string base = scenarioFileName;
        const string xmlExt = ".xml";
        if( base.size() >= xmlExt.size() &&
            base.compare( base.size() - xmlExt.size(), xmlExt.size(), xmlExt ) == 0 ){
            base.erase( base.size() - xmlExt.size() );
        }
        return base + ".csv";
    }

    /// Find the EXTRA_INFECTION_OUTPUT <option> element, if present, so its
    /// optional startDate/endDate attributes can be read.
    const scnXml::Option* findExtraInfectionOption( const scnXml::Model& model ){
        if( !model.getModelOptions().present() ) return nullptr;
        const scnXml::OptionSet::OptionSequence& optSeq =
            model.getModelOptions().get().getOption();
        for( const scnXml::Option& opt : optSeq ){
            if( opt.getName() == "EXTRA_INFECTION_OUTPUT" ) return &opt;
        }
        return nullptr;
    }

    /// Parse a calendar date (YYYY-MM-DD) into an absolute SimTime, throwing a
    /// scenario error if it is not a valid date.
    SimTime parseBound( const string& attr, const string& what ){
        SimTime date = UnitParse::parseDate( attr );
        if( date == sim::never() ){
            throw util::xml_scenario_error(
                "EXTRA_INFECTION_OUTPUT/" + what + ": expected a date (YYYY-MM-DD) but got \"" + attr + "\"" );
        }
        return date;
    }
}

void init( const string& scenarioFileName, const scnXml::Model& model ){
    enabled = util::ModelOptions::option( util::EXTRA_INFECTION_OUTPUT );

    const scnXml::Option* opt = findExtraInfectionOption( model );

    // Read the optional date range from the EXTRA_INFECTION_OUTPUT option.
    dateFilter = false;
    startBound = sim::never();
    endBound = sim::future();
    if( opt != nullptr && (opt->getStartDate().present() || opt->getEndDate().present()) ){
        if( !enabled ){
            cerr << "Warning: startDate/endDate given on the EXTRA_INFECTION_OUTPUT "
                    "option but the option is disabled; date range ignored." << endl;
        } else {
            dateFilter = true;
            // Omitted start => from the beginning of the monitored period (the
            // first step with a valid calendar date). Omitted end => to the end.
            startBound = opt->getStartDate().present()
                ? parseBound( opt->getStartDate().get(), "startDate" )
                : sim::startDate();
            endBound = opt->getEndDate().present()
                ? parseBound( opt->getEndDate().get(), "endDate" )
                : sim::future();
            if( startBound > endBound ){
                throw util::xml_scenario_error(
                    "EXTRA_INFECTION_OUTPUT: startDate is after endDate" );
            }
        }
    }

    if( !enabled ) return;

    const string fileName = deriveOutputName( scenarioFileName );
    outFile.open( fileName.c_str() );
    if( !outFile.good() ){
        throw util::base_exception( "EXTRA_INFECTION_OUTPUT: unable to open " + fileName + " for writing" );
    }

    // Header. Columns shared by both row types are listed first; "human" rows
    // leave the infection columns blank and "infection" rows leave the
    // human-immunity columns blank.
    outFile << "row_type,time_days,human_id,human_age_years,"
            << "m_cumulative_h,m_cumulative_Y,clinical_status,"
            << "infection_id,infection_age,parasite_density,m_cumulativeExposureJ\n";
}

bool isEnabled(){
    return enabled;
}

void report( const Human& human, SimTime time ){
    if( !enabled ) return;

    // Only emit during the intervention (monitored) period, to match the rest of
    // OpenMalaria's monitoring: warm-up and EIR-calibration steps produce no
    // regular monitoring output, so they produce no rows here either.
    // sim::intervTime() is negative during those phases and >= 0 once the
    // intervention period begins.
    if( sim::intervTime() < sim::zero() ) return;

    // Optional calendar-date range within the monitored period.
    if( dateFilter ){
        const SimTime date = sim::intervDate();
        if( date < startBound || date > endBound ) return;
    }

    // Time column: days since the start of the intervention (monitored) period,
    // as used by the rest of OpenMalaria's monitoring (first monitored step is
    // 5). Note human age is computed from absolute time (`time`) since birth
    // dates are absolute.
    const SimTime timeDays = sim::intervTime();
    const uint32_t humanID = human.getId();
    const double ageYears = sim::inYears( human.age( time ) );

    // Human row: row attributes + human attributes, infection columns blank.
    outFile << "human," << timeDays << ',' << humanID << ',' << ageYears << ','
            << human.withinHostModel->getCumulative_h() << ','
            << human.withinHostModel->getCumulative_Y() << ','
            << clinicalStatusCode( human.clinicalModel->getLatestState() )
            << ",,,,\n";     // infection_id, infection_age, parasite_density, m_cumulativeExposureJ (blank)

    // Infection rows: one per infection object, human-only columns blank.
    vector<WithinHost::ReportedInfectionData> infections;
    human.withinHostModel->getInfectionData( infections );
    for( const auto& inf : infections ){
        outFile << "infection," << timeDays << ',' << humanID << ',' << ageYears << ','
                << ",,,"     // m_cumulative_h, m_cumulative_Y, clinical_status (blank)
                << inf.id << ',' << inf.ageDays << ',' << inf.density << ','
                << inf.cumulativeExposureJ << '\n';
    }
}

void close(){
    if( outFile.is_open() ) outFile.close();
    enabled = false;
    dateFilter = false;
}

} } }
