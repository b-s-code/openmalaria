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

#ifndef Hmod_ExtraInfectionOutput
#define Hmod_ExtraInfectionOutput

#include "Global.h"
#include <string>

namespace scnXml { class Model; }

namespace OM { namespace Host {

class Human;

/** Optional per-timestep CSV output enabled by the EXTRA_INFECTION_OUTPUT model
 * option.
 *
 * When enabled, a single CSV file (named after the scenario XML file) is written
 * containing, for every time step:
 *  - one "human" row per human (with that human's immunity summary), and
 *  - one "infection" row per Infection object belonging to each human.
 *
 * This is a debugging/analysis output and can produce very large files. */
namespace ExtraInfectionOutput {

    /** If the EXTRA_INFECTION_OUTPUT option is active, open the output CSV file
     * (derived from `scenarioFileName`) and write its header row. Otherwise do
     * nothing. Must be called after model options and sim time have been
     * initialised.
     *
     * `model` is inspected for optional startDate/endDate attributes on the
     * EXTRA_INFECTION_OUTPUT `<option>` element, which restrict output to a
     * calendar-date range (see report()). */
    void init( const std::string& scenarioFileName, const scnXml::Model& model );

    /** True if output is enabled and the file is open. */
    bool isEnabled();

    /** Write the human row and one infection row per infection for `human` at
     * the given simulation time (in days). No-op when output is disabled or when
     * the current step's calendar date falls outside the configured date range.
     *
     * `time` is the internal simulation time in days (as written to the
     * time_days column); date-range filtering is applied separately using the
     * current calendar date (sim::intervDate()). */
    void report( const Human& human, SimTime time );

    /** Flush and close the output file. Safe to call when disabled. */
    void close();

} } }
#endif
