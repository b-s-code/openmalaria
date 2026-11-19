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

#include "Host/WithinHost/Diagnostic.h"
#include "Parameters.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "schema/scenario.h"
#include <limits>

namespace OM { namespace WithinHost {

// ———  Diagnostic (non-static)  ———

Diagnostic::Diagnostic( const Parameters& parameters, const scnXml::Diagnostic& elt ){
    m_densLim = numeric_limits<double>::quiet_NaN();
    m_specificity = numeric_limits<double>::quiet_NaN();
    
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        return; // both parameters NaN (set above)
    }
    
    if( elt.getDeterministic().present() ){
        m_densLim = elt.getDeterministic().get().getMinDensity();
    }else if( elt.getStochastic().present() ){
        m_densLim = elt.getStochastic().get().getDens_50();
        if( m_densLim == 0.0 ){
            // The equation used for stochastic diagnostics breaks down when
            // dens=dens_lim=0 and for other cases the deterministic model is
            // the same when dens_lim=0.
        }else{
            m_specificity = elt.getStochastic().get().getSpecificity();
            if( m_specificity < 0.0 || m_specificity > 1.0 ){
                throw util::XmlScenarioError(
                    string("diagnostics/diagnostic(").append(elt.getName())
                    .append("): specificity must be in range [0,1]") );
            }
        }
    }else{
        // This should be impossible since according to schema one of these
        // elements must be present.
        throw SWITCH_DEFAULT_EXCEPTION;
    }
    if( m_densLim < 0.0 ){
        throw util::XmlScenarioError(
            string("diagnostics/diagnostic(").append(elt.getName())
            .append("): must have density ≥ 0") );
    }
    
    // We use a bias factor to adjust the "units" used to specify the density
    // of this diagnostic, since estimates from Garki and the standard
    // non-Garki sources are not equivalent to those from the Malariatherapy
    // data (which is used internally).
    if( !elt.getUnits().present() ){
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            // User must be explicit in this case, because presumably the Garki
            // bias is to be used for some diagnostics but likely not all
            // (e.g. neonatal mortality).
            throw util::XmlScenarioError( "diagnostics/diagnostic(*)/units: must specify this attribute when GARKI_DENSITY_BIAS is set" );
        }
        // otherwise we assume "Other"
        m_densLim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Other" ){
        m_densLim *= parameters[Parameter::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Garki" ){
        m_densLim *= parameters[Parameter::DENSITY_BIAS_GARKI];
    }else{
        assert( elt.getUnits().get() == "Malariatherapy" );
        // in this case we don't need to use a bias factor
    }
    
    m_usesHrp2 = elt.getMechanism() == "HRP2";
}

Diagnostic::Diagnostic(double minDens){
    m_specificity = numeric_limits<double>::quiet_NaN();
    m_densLim = minDens;
    m_usesHrp2 = false;
}

bool Diagnostic::isPositive( LocalRng& rng, double dens, double densHRP2 ) const {
    if( m_usesHrp2 ){
        assert( densHRP2 == densHRP2 ); // monitoring diagnostic passes NaN; use of HRP2 is not supported
        dens = densHRP2;
    }
    if( (std::isnan)(m_specificity) ){
        // use deterministic test
        return dens >= m_densLim;
    }else{
        // dens_lim is dens_50 in this case
        double pPositive = 1.0 + m_specificity * (dens / (dens + m_densLim) - 1.0);
//         double pPositive = (dens + dens_lim - dens_lim * specificity) / (dens + dens_lim);       // equivalent
        return rng.bernoulli( pPositive );
    }
}

bool Diagnostic::allowsFalsePositives() const{
    if( (std::isnan)(m_specificity) ) return m_densLim <= 0.0;
    return m_specificity < 1.0;
}


// ———  diagnostics (static)  ———

map<string, Diagnostic> diagnosticSet;
const Diagnostic* Diagnostics::s_monitoringDiagnostic = 0;

void Diagnostics::clear(){
    diagnosticSet.clear();
}

void Diagnostics::init( const Parameters& parameters, const scnXml::Scenario& scenario ){
    diagnosticSet.clear(); // compatibility with unit tests
    
    if(scenario.getDiagnostics().present()){
        for( const scnXml::Diagnostic& diagnostic : scenario.getDiagnostics().get().getDiagnostic() ){
            string name = diagnostic.getName();     // conversion fails without this extra line
            bool inserted = diagnosticSet.insert( make_pair(std::move(name), Diagnostic(parameters, diagnostic)) ).second;
            if( !inserted ){
                throw util::XmlScenarioError( string("diagnostic with this name already set: ").append(diagnostic.getName()) );
            }
        }
    }

    const scnXml::Surveys& surveys = scenario.getMonitoring().getSurveys();
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        // So far the implemented Vivax code does not produce parasite
        // densities, thus this diagnostic model cannot be used.
        Diagnostics::s_monitoringDiagnostic = &Diagnostics::makeDeterministic( numeric_limits<double>::quiet_NaN() );
    }else if( surveys.getDetectionLimit().present() ){
        if( surveys.getDiagnostic().present() ){
            throw util::XmlScenarioError( "monitoring/surveys: do not "
                "specify both detectionLimit and diagnostic" );
        }
        if( util::CommandLine::option( util::CommandLine::DEPRECATION_WARNINGS ) ){
            std::cerr << "Deprecation warning: monitoring/surveys: "
                "specification of \"diagnostic\" is suggested over \"detectionLimit\"" << std::endl;
        }

        // This controls whether the detection limit is specified relative to
        // the Garki or other methods.
        double densitybias = numeric_limits<double>::quiet_NaN();
        if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
            densitybias = parameters[Parameter::DENSITY_BIAS_GARKI];
        } else {
            if( scenario.getAnalysisNo().present() ){
                int analysisNo = scenario.getAnalysisNo().get();
                if ((analysisNo >= 22) && (analysisNo <= 30)) {
                    cerr << "Warning: these analysis numbers used to mean "
                        "use Garki density bias. If you do want to use this, "
                        "specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
                }
            }
            densitybias = parameters[Parameter::DENSITY_BIAS_NON_GARKI];
        }
        double detectionLimit = surveys.getDetectionLimit().get() * densitybias;
        Diagnostics::s_monitoringDiagnostic = &Diagnostics::makeDeterministic( detectionLimit );
    }else{
        if( !surveys.getDiagnostic().present() ){
            throw util::XmlScenarioError( "monitoring/surveys: require "
                "either detectionLimit or diagnostic" );
        }
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            throw util::XmlScenarioError( "Use of GARKI_DENSITY_BIAS is not "
                "appropriate when monitoring/surveys/diagnostic is used." );
        }
        Diagnostics::s_monitoringDiagnostic = &Diagnostics::get( surveys.getDiagnostic().get() );
    }
    
    if( Diagnostics::s_monitoringDiagnostic->m_usesHrp2 ){
        throw util::XmlScenarioError( "the diagnostic used for "
            "monitoring may not use HRP2 as its mechanism" );
    }
}

const Diagnostic& Diagnostics::get( const string& name ){
    auto it = diagnosticSet.find(name);
    if( it == diagnosticSet.end() ){
        throw util::XmlScenarioError( string("diagnostic not found: ").append(name) );
    }
    return it->second;
}

const Diagnostic& Diagnostics::makeDeterministic(double minDens){
    string name = "";   // anything not matching an existing name is fine
    assert( diagnosticSet.count(name) == 0 );  // we shouldn't need to call make_deterministic more than once, so I think this will do
    auto ret = diagnosticSet.insert( make_pair(std::move(name), Diagnostic( minDens )) );
    assert(ret.second);
    return ret.first->second;
}

} }
