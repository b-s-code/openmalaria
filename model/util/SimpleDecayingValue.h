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

#ifndef Hmod_SimpleDecayingValue
#define Hmod_SimpleDecayingValue

#include "Global.h"
#include "util/DecayFunction.h"

namespace OM {
namespace util {

/** Parameters to for a decaying intervention whose effect can be described
 * by a single scalar output.
 * 
 * This class outputs a scalar which defaults to and decays to 0, but
 * immediately after deployment will be the set _initial value_.
 *****************************************************************************/
class SimpleDecayingValue
{
public:
    /** Default construction: always return 0. */
    SimpleDecayingValue() :
            m_initial(0.0),
            m_deployT(Sim::never()) {}
    
    /** Configure from an XML element. */
    inline void set (double initialValue, const scnXml::DecayFunction& elt, const char* name){
        m_decay = DecayFunction::makeObject(elt, name);
        m_initial = initialValue;
    }

    /** Configure from an XML element. */
    inline void set (double initialValue){
        m_initial = initialValue;
    }
    
    /** Trigger a deployment, if decay function was set. */
    inline void deploy (LocalRng& rng, SimTime time){
        if( m_decay.get() == 0 ) return;  // cannot deploy
        m_deployT = time;
        m_het = m_decay->hetSample(rng);
    }
    
    /** Get the value (0 if before any deployment or after complete decay,
     * also 0 if no decay function or initial value was set,
     * otherwise between zero and the inital value). */
    inline double currentValue (SimTime time) const{
        if( m_decay.get() == 0 || m_het == nullptr) return 0.0;  // decay wasn't set: the always return 0
        else return m_initial * m_het->eval( time - m_deployT );
    }
    
    /** Checkpointing: only checkpoint parameters which change after initial
     * set-up. */
    template<class S>
    void operator& (S& stream) {
        if (m_decay) m_decay->checkpointSample(m_het, stream);
        m_deployT & stream;
    }
    
private:
    /** Description of decay of effects on emergence. */
    unique_ptr<util::DecayFunction> m_decay;
    
    /** Initial value. Is initialised to 0. */
    double m_initial;
    
    /** Description of larviciding decay. */
    unique_ptr<util::DecayFunction> m_het;
    
    /** Time of larviciding deployment. */
    SimTime m_deployT = Sim::never();
};

} }
#endif
