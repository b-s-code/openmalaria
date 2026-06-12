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
#include "util/ModelOptions.h"
#include "util/errors.h"

#include <fstream>
#include <vector>

namespace OM { namespace Host { namespace ExtraInfectionOutput {

using namespace std;

namespace {
    ofstream outFile;
    bool enabled = false;

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
}

void init( const string& scenarioFileName ){
    enabled = util::ModelOptions::option( util::EXTRA_INFECTION_OUTPUT );
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
            << "m_cumulative_h,m_cumulative_Y,"
            << "infection_id,infection_age,parasite_density,m_cumulativeExposureJ\n";
}

bool isEnabled(){
    return enabled;
}

void report( const Human& human, SimTime time ){
    if( !enabled ) return;

    const uint32_t humanID = human.getId();
    const double ageYears = sim::inYears( human.age( time ) );

    // Human row: row attributes + human attributes, infection columns blank.
    outFile << "human," << time << ',' << humanID << ',' << ageYears << ','
            << human.withinHostModel->getCumulative_h() << ','
            << human.withinHostModel->getCumulative_Y() << ','
            << ",,,\n";

    // Infection rows: one per infection object, human-immunity columns blank.
    vector<WithinHost::ReportedInfectionData> infections;
    human.withinHostModel->getInfectionData( infections );
    for( const auto& inf : infections ){
        outFile << "infection," << time << ',' << humanID << ',' << ageYears << ','
                << ",,"      // m_cumulative_h, m_cumulative_Y (blank)
                << inf.id << ',' << inf.ageDays << ',' << inf.density << ','
                << inf.cumulativeExposureJ << '\n';
    }
}

void close(){
    if( outFile.is_open() ) outFile.close();
    enabled = false;
}

} } }
