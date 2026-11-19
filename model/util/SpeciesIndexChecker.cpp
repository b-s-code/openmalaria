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

#include "util/SpeciesIndexChecker.h"
#include "util/errors.h"
#include <sstream>

namespace OM { namespace util {
using std::ostringstream;

size_t SpeciesIndexChecker::getIndex( string species ){
    if( m_found.count( species ) > 0 ){
        ostringstream msg;
        msg << "Intervention \"" << m_intervName
            << "\" has multiple descriptions for vector species \""
            << species << "\"";
        throw util::XmlScenarioError( msg.str() );
    }
    auto sIndex = m_indices.find( species );
    if( sIndex == m_indices.end() ){
        ostringstream msg;
        msg << "Intervention \"" << m_intervName
            << "\" has a description for vector species \"" << species
            << "\", but this species is not mentioned in the entomology section";
        throw util::XmlScenarioError( msg.str() );
    }
    m_found.insert( species );
    return sIndex->second;
}

void SpeciesIndexChecker::checkNoneMissed() const{
    for( auto it = m_indices.begin(); it != m_indices.end(); ++it ){
        if( m_found.count( it->first ) == 0 ){
            ostringstream msg;
            msg << "Intervention \"" << m_intervName
                << "\" has a no description for vector species \""
                << it->first << "\"";
            throw util::XmlScenarioError( msg.str() );
        }
    }
}

 } }
 
