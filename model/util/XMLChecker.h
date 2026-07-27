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

#ifndef Hmod_util_XMLChecker
#define Hmod_util_XMLChecker

#include "Global.h"
#include <schema/scenario.h>

namespace OM
{
    namespace util
    {
        class XMLChecker
        {
        public:
            XMLChecker(){}
            ~XMLChecker(){}

            /*
            * Should be called after the input XML has been loaded and before all subsequent
            * initialisation. I.e. well before the simulation starts running.
            *
            * Performs custom checks on the input XML which are not enforced in the schema itself.
            *
            * It would be better to handle as many such issues as possible in the schema itself
            * However, some required checks are not possible to do using the version of the XSD
            * spec supported by the XML validation library used by OpenMalaria.
            *
            * Throws iff a check fails.
            *
            * The purpose of these XML checks to:
            *  - identify certain problems in the input XML as early as possible, and
            *  - enable the user to obtain a more informative error message.
            */
            void PerformPostValidationChecks(const scnXml::Scenario& scenario)
            {
                CheckModelOptionsAndParams(scenario);
                CheckVaxEfficacyVsCumulativeExposure(scenario);
            }
        private:

            /*
            * Verifies that, if no model name is written in input XML, then both parameters and
            * model options and written explicitly.
            */
            void CheckModelOptionsAndParams(const scnXml::Scenario& scenario)
            {
                // For each relevant element, determine whether it's specified explicitly in XML.
                const bool modelName = scenario.getModel().getModelName().present();
                const bool parameters = scenario.getModel().getParameters().present();
                const bool modelOptions = scenario.getModel().getModelOptions().present();
                if (!modelName && (!parameters || !modelOptions))
                {
                    throw util::xml_scenario_error(
                        "If a model name is not specified then both <ModelOptions> and <parameters> must be specified");
                }
            }

            /*
            * If a vaccine-efficacy-vs-cumulative-exposure model option is enabled in the
            * XML, then every vaccine description (PEV / BSV / TBV) must specify the optional
            * coefficient attribute corresponding to that option:
            *   - VAX_EFFICACY_VS_CUMULATIVE_INFS    requires cumulativeInfsCoeff
            *   - VAX_EFFICACY_VS_CUMULATIVE_DENSITY requires cumulativeDensityCoeff
            * Throw otherwise.
            *
            * Both options may be enabled together in a single scenario, allowing
            * different vaccine descriptions to use different variants. When both are
            * enabled, each vaccine description need only specify at least one of the two
            * coefficients; an absent coefficient is treated as zero (a no-op, exp(0)=1).
            * However, a single vaccine description may not use both variants at once: it
            * is an error for both cumulativeInfsCoeff and cumulativeDensityCoeff to be
            * present and non-zero on the same description. Setting one coefficient to
            * zero (or omitting it) is how a description selects which single variant it
            * uses (or opts out of the feature entirely).
            */
            void CheckVaxEfficacyVsCumulativeExposure(const scnXml::Scenario& scenario)
            {
                const bool needInfs = IsModelOptionEnabled(scenario, "VAX_EFFICACY_VS_CUMULATIVE_INFS");
                const bool needDensity = IsModelOptionEnabled(scenario, "VAX_EFFICACY_VS_CUMULATIVE_DENSITY");
                if (!needInfs && !needDensity)
                    return;

                if (!scenario.getInterventions().getHuman().present())
                    return;

                const scnXml::HumanInterventions& human = scenario.getInterventions().getHuman().get();
                for (auto it = human.getComponent().begin(), end = human.getComponent().end(); it != end; ++it)
                {
                    const scnXml::HumanInterventionComponent& component = *it;
                    if (component.getPEV().present())
                        RequireCoeffs(component.getPEV().get(), "PEV", component.getId(), needInfs, needDensity);
                    if (component.getBSV().present())
                        RequireCoeffs(component.getBSV().get(), "BSV", component.getId(), needInfs, needDensity);
                    if (component.getTBV().present())
                        RequireCoeffs(component.getTBV().get(), "TBV", component.getId(), needInfs, needDensity);
                }
            }

            /*
            * Returns true if the named option appears in the scenario's <ModelOptions>
            * with an explicit or defaulted value of true. Returns false if the option
            * is absent or explicitly disabled.
            */
            bool IsModelOptionEnabled(const scnXml::Scenario& scenario, const std::string& optionName)
            {
                if (!scenario.getModel().getModelOptions().present())
                    return false;

                const scnXml::OptionSet::OptionSequence& optSeq =
                    scenario.getModel().getModelOptions().get().getOption();
                for (auto it = optSeq.begin(); it != optSeq.end(); ++it)
                {
                    if (it->getName() == optionName)
                        return it->getValue();
                }
                return false;
            }

            void RequireCoeffs(const scnXml::VaccineDescription& vd,
                               const std::string& vaccineType,
                               const std::string& componentId,
                               bool needInfs,
                               bool needDensity)
            {
                const bool infsPresent = vd.getCumulativeInfsCoeff().present();
                const bool densityPresent = vd.getCumulativeDensityCoeff().present();

                if (needInfs && needDensity)
                {
                    // Both options enabled: at least one coefficient must be given; the
                    // other is assumed zero (a no-op).
                    if (!infsPresent && !densityPresent)
                    {
                        throw util::xml_scenario_error(
                            "Vaccine component \"" + componentId + "\" (" + vaccineType +
                            "): at least one of cumulativeInfsCoeff / cumulativeDensityCoeff "
                            "is required when both the VAX_EFFICACY_VS_CUMULATIVE_INFS and "
                            "VAX_EFFICACY_VS_CUMULATIVE_DENSITY model options are enabled");
                    }
                }
                else if (needInfs && !infsPresent)
                {
                    throw util::xml_scenario_error(
                        "Vaccine component \"" + componentId + "\" (" + vaccineType +
                        "): cumulativeInfsCoeff attribute is required when the "
                        "VAX_EFFICACY_VS_CUMULATIVE_INFS model option is enabled");
                }
                else if (needDensity && !densityPresent)
                {
                    throw util::xml_scenario_error(
                        "Vaccine component \"" + componentId + "\" (" + vaccineType +
                        "): cumulativeDensityCoeff attribute is required when the "
                        "VAX_EFFICACY_VS_CUMULATIVE_DENSITY model option is enabled");
                }

                // A single vaccine description may use at most one variant: reject the
                // combined form where both coefficients are present and non-zero. This
                // can only arise when both options are enabled.
                if (needInfs && needDensity &&
                    infsPresent && vd.getCumulativeInfsCoeff().get() != 0.0 &&
                    densityPresent && vd.getCumulativeDensityCoeff().get() != 0.0)
                {
                    throw util::xml_scenario_error(
                        "Vaccine component \"" + componentId + "\" (" + vaccineType +
                        "): a single vaccine description may not use both the "
                        "VAX_EFFICACY_VS_CUMULATIVE_INFS and VAX_EFFICACY_VS_CUMULATIVE_DENSITY "
                        "variants at once; at most one of cumulativeInfsCoeff / "
                        "cumulativeDensityCoeff may be non-zero. Set the unused coefficient to 0");
                }
            }
        };
    }
}

#endif
