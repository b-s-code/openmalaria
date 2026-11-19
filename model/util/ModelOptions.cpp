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

#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "schema/util.h"

#include <sstream>
#include <iostream>

namespace OM { namespace util {
    std::bitset<NUM_OPTIONS> ModelOptions::s_options;
    
    // Utility: converts option strings to codes and back
    class OptionCodeMap {
	// Lookup table to translate the strings used in the XML file to the internal enumerated values:
	map<string,OptionCodes> m_codeMap;
        set<string> m_ignoreOptions;
	
    public:
	OptionCodeMap () {
// 	    codeMap["PENALISATION_EPISODES"] = PENALISATION_EPISODES;
	    m_codeMap["NEGATIVE_BINOMIAL_MASS_ACTION"] = NEGATIVE_BINOMIAL_MASS_ACTION;
	    // codeMap["ATTENUATION_ASEXUAL_DENSITY"] = ATTENUATION_ASEXUAL_DENSITY;
	    m_codeMap["LOGNORMAL_MASS_ACTION"] = LOGNORMAL_MASS_ACTION;
	    m_codeMap["NO_PRE_ERYTHROCYTIC"] = NO_PRE_ERYTHROCYTIC;
	    m_codeMap["MAX_DENS_CORRECTION"] = MAX_DENS_CORRECTION;
	    m_codeMap["INNATE_MAX_DENS"] = INNATE_MAX_DENS;
	    // 	codeMap["MAX_DENS_RESET"] = MAX_DENS_RESET;
	    m_codeMap["DUMMY_WITHIN_HOST_MODEL"] = DUMMY_WITHIN_HOST_MODEL;
	    m_codeMap["PREDETERMINED_EPISODES"] = PREDETERMINED_EPISODES;
	    m_codeMap["NON_MALARIA_FEVERS"] = NON_MALARIA_FEVERS;
        m_ignoreOptions.insert("INCLUDES_PK_PD");
	    m_codeMap["CLINICAL_EVENT_SCHEDULER"] = CLINICAL_EVENT_SCHEDULER;
	    m_codeMap["MUELLER_PRESENTATION_MODEL"] = MUELLER_PRESENTATION_MODEL;
	    m_codeMap["TRANS_HET"] = TRANS_HET;
	    m_codeMap["COMORB_HET"] = COMORB_HET;
	    m_codeMap["TREAT_HET"] = TREAT_HET;
	    m_codeMap["COMORB_TRANS_HET"] = COMORB_TRANS_HET;
	    m_codeMap["TRANS_TREAT_HET"] = TRANS_TREAT_HET;
	    m_codeMap["COMORB_TREAT_HET"] = COMORB_TREAT_HET;
	    m_codeMap["TRIPLE_HET"] = TRIPLE_HET;
	    m_codeMap["EMPIRICAL_WITHIN_HOST_MODEL"] = EMPIRICAL_WITHIN_HOST_MODEL;
        m_codeMap["MOLINEAUX_WITHIN_HOST_MODEL"] = MOLINEAUX_WITHIN_HOST_MODEL;
        m_codeMap["PENNY_WITHIN_HOST_MODEL"] = PENNY_WITHIN_HOST_MODEL;
        m_codeMap["GARKI_DENSITY_BIAS"] = GARKI_DENSITY_BIAS;
//      codeMap["IPTI_SP_MODEL"] = IPTI_SP_MODEL;
//      codeMap["REPORT_ONLY_AT_RISK"] = REPORT_ONLY_AT_RISK; 
	    m_codeMap["MEAN_DURATION_GAMMA"] = MEAN_DURATION_GAMMA;
	    m_codeMap["FIRST_LOCAL_MAXIMUM_GAMMA"] = FIRST_LOCAL_MAXIMUM_GAMMA;
	    m_codeMap["IMMUNE_THRESHOLD_GAMMA"] = IMMUNE_THRESHOLD_GAMMA;
	    m_codeMap["UPDATE_DENSITY_GAMMA"] = UPDATE_DENSITY_GAMMA;
	    m_codeMap["PARASITE_REPLICATION_GAMMA"] = PARASITE_REPLICATION_GAMMA;
        m_codeMap["VECTOR_LIFE_CYCLE_MODEL"] = VECTOR_LIFE_CYCLE_MODEL;
        m_codeMap["VECTOR_SIMPLE_MPD_MODEL"] = VECTOR_SIMPLE_MPD_MODEL;
        m_codeMap["MOLINEAUX_PAIRWISE_SAMPLE"] = MOLINEAUX_PAIRWISE_SAMPLE;
        m_ignoreOptions.insert("PROPHYLACTIC_DRUG_ACTION_MODEL");
        m_codeMap["VIVAX_SIMPLE_MODEL"] = VIVAX_SIMPLE_MODEL;
        m_codeMap["INDIRECT_MORTALITY_FIX"] = INDIRECT_MORTALITY_FIX;
        m_codeMap["VACCINE_GENOTYPE"] = VACCINE_GENOTYPE;
        m_codeMap["CFR_PF_USE_HOSPITAL"] = CFR_PF_USE_HOSPITAL;
        m_codeMap["HEALTH_SYSTEM_MEMORY_FIX"] = HEALTH_SYSTEM_MEMORY_FIX;
        m_codeMap["USE_EXACT_NV0_SOLVER"] = USE_EXACT_NV0_SOLVER;
	}
	
	OptionCodes operator[] (const string s) {
	    auto codeIt = m_codeMap.find (s);
	    if( codeIt != m_codeMap.end() ) return codeIt->second;
            if( m_ignoreOptions.count(s) ){
                if( CommandLine::option(CommandLine::DEPRECATION_WARNINGS) ){
                    cerr << "Deprecation warning: model option " << s << " is no longer used" << endl;
                }
                return IGNORE;
            }
            
            ostringstream msg;
            if( s == "PENALISATION_EPISODES" || s == "ATTENUATION_ASEXUAL_DENSITY" ){
                msg << "Please use schema 31 or earlier to use option "
                    << s << "; it is not available in later versions.";
            }else if( s == "IPTI_SP_MODEL" ){
                msg << "The IPT model is no longer available. Use MDA instead.";
            }else if( s == "REPORT_ONLY_AT_RISK" ){
                msg << "Option " << s << " has been replaced by <SurveyOptions onlyNewEpisode=\"true\">.";
            }else{
                msg << "Unrecognised model option: " << s;
            }
            throw XmlScenarioError(msg.str());
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const OptionCodes code) {
	    for(auto codeIt = m_codeMap.begin(); codeIt != m_codeMap.end(); ++codeIt) {
		if (codeIt->second == code)
		    return codeIt->first;
	    }
	    throw TRACED_EXCEPTION_DEFAULT ("toString called with unknown code");	// this is a code error
	}
    };

    std::bitset<NUM_OPTIONS> ModelOptions::getBaseModelOptions()
    {
        // The "base" model is very simple
        bitset<NUM_OPTIONS> allOptionsOff;

        return allOptionsOff;
    }

    std::bitset<NUM_OPTIONS> ModelOptions::getLegacyDefaultModelOptions()
    {
        bitset<NUM_OPTIONS> defaultOptSet;

        // Kept for legacy reasons.  That is, when an input XML does not specify a model name
        // and does not set these model options to false, they are turned on *for backwards
        // compatibility*.
        defaultOptSet.set (MAX_DENS_CORRECTION);
        defaultOptSet.set (INNATE_MAX_DENS);
        defaultOptSet.set (INDIRECT_MORTALITY_FIX);
        defaultOptSet.set (HEALTH_SYSTEM_MEMORY_FIX);

        return defaultOptSet;
    }

    void ModelOptions::init(scnXml::Model::ModelOptionsOptional& optionsElt,
                            util::ModelNameProvider mnp)
    {
    OptionCodeMap codeMap;

	// Set options to legacy defaults, then override any given in the XML file:
	s_options = getLegacyDefaultModelOptions();

    const util::ModelNames namedModelToUse = mnp.getModelName();
    const bool useNamedModel = namedModelToUse != util::ModelNames::none;
    if (useNamedModel)
    {
        if (namedModelToUse == util::ModelNames::base)
        {
            // This completely discards any legacy model options set above.
            s_options = getBaseModelOptions();
        }
        else
        {
            throw util::XmlScenarioError(
                "No pre-set model options are available for the specified model name.");
        }
    }

    // Apply any user-specified overrides.
    if (optionsElt.present())
    {
	    const scnXml::OptionSet::OptionSequence& optSeq = optionsElt.get().getOption();
	    for(auto it = optSeq.begin(); it != optSeq.end(); ++it) {
                OptionCodes opt = codeMap[it->getName()];
                if( opt != IGNORE ) s_options[opt] = it->getValue();
	    }
    }
	
	// Print non-default model options:
	if (CommandLine::option (CommandLine::PRINT_MODEL_OPTIONS)) {
        const std::bitset<NUM_OPTIONS> defaultOptSet = getLegacyDefaultModelOptions();
	    cout << "Non-default model options:";
	    for(int i = 0; i < NUM_OPTIONS; ++i) {
		if (s_options[i] != defaultOptSet[i])
		    cout << "\t" << codeMap.toString(OptionCodes(i)) << "=" << s_options[i];
	    }
	    cout << endl;
	}
	
	// Test for incompatible options
	
	// An incompatibility triangle, listing options incompatible with each option
	// Doesn't check required versions
	// Originally from "description of variables for interface" excel sheet
	bitset<NUM_OPTIONS> incompatibilities[NUM_OPTIONS];	// all default to 0
	
	incompatibilities[NEGATIVE_BINOMIAL_MASS_ACTION]
	    .set(LOGNORMAL_MASS_ACTION)
	    .set(TRANS_HET)	.set(COMORB_TRANS_HET)
	    .set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[LOGNORMAL_MASS_ACTION]
	    .set(TRANS_HET)	.set(COMORB_TRANS_HET)
	    .set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	
	// Note: MAX_DENS_CORRECTION is irrelevant when using new
	// within-host models, but we don't mark it incompatible so that we can
	// leave MAX_DENS_CORRECTION on by default.
	incompatibilities[DUMMY_WITHIN_HOST_MODEL]
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(EMPIRICAL_WITHIN_HOST_MODEL)
	    .set(MOLINEAUX_WITHIN_HOST_MODEL)
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);

	incompatibilities[EMPIRICAL_WITHIN_HOST_MODEL]
            .set(MOLINEAUX_WITHIN_HOST_MODEL)
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);
	    
	    
	incompatibilities[MOLINEAUX_WITHIN_HOST_MODEL]
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);
	    
	incompatibilities[PENNY_WITHIN_HOST_MODEL]
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA);
	    
	incompatibilities[NON_MALARIA_FEVERS]
	    .set(MUELLER_PRESENTATION_MODEL);
	
	incompatibilities[TRANS_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TREAT_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TRANS_HET]
	    .set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TRANS_TREAT_HET]
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TREAT_HET]
	    .set(TRIPLE_HET);
        
        incompatibilities[MOLINEAUX_PAIRWISE_SAMPLE]
            .set(FIRST_LOCAL_MAXIMUM_GAMMA)
            .set(MEAN_DURATION_GAMMA);
        
        
        incompatibilities[VIVAX_SIMPLE_MODEL]
            .set( DUMMY_WITHIN_HOST_MODEL )
            .set( EMPIRICAL_WITHIN_HOST_MODEL )
            .set( MOLINEAUX_WITHIN_HOST_MODEL )
            .set( PENNY_WITHIN_HOST_MODEL );
        
	for(size_t i = 0; i < NUM_OPTIONS; ++i) {
	    if (s_options [i] && (s_options & incompatibilities[i]).any()) {
		ostringstream msg;
		msg << "Incompatible model options: " << codeMap.toString(OptionCodes(i)) << "=" << s_options[i]
			<< " is incompatible with flags:";
		bitset<NUM_OPTIONS> incompat = (s_options & incompatibilities[i]);
		for(int j = 0; j < NUM_OPTIONS; ++j) {
		    if (incompat[j])
			msg << "\t" << codeMap.toString(OptionCodes(j)) << "=" << s_options[j];
		}
		throw XmlScenarioError (msg.str());
	    }
	}
	
        // Required options (above table can't check these):
        if (s_options[INNATE_MAX_DENS] && !s_options[MAX_DENS_CORRECTION])
            throw XmlScenarioError ("INNATE_MAX_DENS requires MAX_DENS_CORRECTION");
        if( !s_options[MOLINEAUX_WITHIN_HOST_MODEL] && (
            s_options[FIRST_LOCAL_MAXIMUM_GAMMA] ||
            s_options[MEAN_DURATION_GAMMA] ||
            s_options[PARASITE_REPLICATION_GAMMA] ) )
            throw XmlScenarioError( "Molineaux model option used without MOLINEAUX_WITHIN_HOST_MODEL option" );
        if( !s_options[PENNY_WITHIN_HOST_MODEL] && (
            s_options[IMMUNE_THRESHOLD_GAMMA] ||
            s_options[UPDATE_DENSITY_GAMMA] ) )
            throw XmlScenarioError( "Penny model option used without PENNY_WITHIN_HOST_MODEL option" );
        
        if( Sim::oneTS() == Sim::fromDays(5) ){
            // 5 day TS is okay; some tests specific to this TS:
            bitset<NUM_OPTIONS> require1DayTS;
            require1DayTS
                .set( CLINICAL_EVENT_SCHEDULER );
            
            for(size_t i = 0; i < NUM_OPTIONS; ++i) {
                if (s_options [i] && require1DayTS[i]) {
                    ostringstream msg;
                    msg << "Model option " << codeMap.toString(OptionCodes(i)) << " is only compatible with a 1-day time step.";
                    throw XmlScenarioError (msg.str());
                }
            }
        }else if( Sim::oneTS() == Sim::fromDays(1) ){
            // 1 day TS is also okay
        }else{
            ostringstream msg;
            msg << "Time step set to " << Sim::oneTS() << " days but only 1 and 5 days are supported.";
            throw XmlScenarioError (msg.str());
        }
    }
    
} }
