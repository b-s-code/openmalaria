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

#ifndef OM_IMPORTED_INFECTIONS
#define OM_IMPORTED_INFECTIONS

#include "Global.h"
#include "schema/interventions.h"
#include "util/errors.h"
#include "util/UnitParse.h"
#include "Host/Human.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Infection/Infection.h"

namespace OM {
    class Population;

namespace Host {

    class ImportedInfections {
    public:
        ImportedInfections() : m_period(Sim::zero()), m_lastIndex(0) {}
        
        /** Initialise, passing intervention description
         * 
         * @param iiElt The scenario element with description of rates
         * @returns True if any infections are imported
         */
        void init( const scnXml::ImportedInfections& iiElt )
        {
            const scnXml::ImportedInfections::TimedType& tElt = iiElt.getTimed();
            try{
                //NOTE: if changing XSD, this should not have a default unit:
                m_period = UnitParse::readDuration( tElt.getPeriod(), UnitParse::STEPS );
                if( m_period < Sim::zero() ){
                    throw util::FormatError( "cannot be negative" );
                }
            }catch( const util::FormatError& e ){
                throw util::XmlScenarioError( string("interventions/importedInfections/timed/period: ").append(e.message()) );
            }
            m_rate.reserve( tElt.getRate().size() );
            try{
                for( auto it = tElt.getRate().begin(); it != tElt.getRate().end(); ++it ){
                    SimTime date = UnitParse::readDate( it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/ );
                    // convert to per-time-step, per-person
                    double rateVal = it->getValue() * Sim::yearsPerStep() * (1.0 / 1000.0);
                    m_rate.push_back( Rate( date - Sim::startDate(), rateVal ) );
                }
                sort( m_rate.begin(), m_rate.end() );
                if( m_rate.size() > 0 ){
                    if( m_period != Sim::zero() && m_rate[0].time != Sim::zero() ){
                        throw util::XmlScenarioError( "must specify rate at time zero when period is not zero" );
                    }
                    // remove useless repeated entries from list
                    double lastRate = m_rate[0].value;
                    for( size_t i = 1; i < m_rate.size(); ){
                        if( m_rate[i].value == lastRate ){
                            m_rate.erase( m_rate.begin() + i );
                        }else{
                            lastRate = m_rate[i].value;
                            ++i;
                        }
                    }
                }
            }catch( const util::FormatError& e ){
                throw util::XmlScenarioError( string("interventions/importedInfections/timed/time: ").append(e.message()) );
            }
        }
        
        /** Import this time-step's imported infections, according to initialised rates
         * 
         *  The probability of an host to import an infection is calculated
         *  from the importedInfectionsPerThousandHosts. The bernoulli distribution
         *  is then used to predict if an human has imported the infection in the
         *  population or not. A maximum of one infection can be imported per
         *  person.
         * 
         * @param pop The Population class encapsulating all humans */
        void import(vector<Human> &population)
        {
            if( m_rate.size() == 0 ) return;      // no imported infections
            SimTime now = Sim::intervTime();
            assert( now >= Sim::zero() );
            if( m_period > Sim::zero() ){
                now = modNn(now, m_period);
            }
            if( m_rate[m_lastIndex].time > now ){
                m_lastIndex = 0;  // gone round in a loop: back to start of period
            }
            if( now > m_rate[m_lastIndex].time && m_rate.size() > m_lastIndex+1 ){
                if( m_rate[m_lastIndex+1].time <= now ){
                    ++m_lastIndex;
                    // won't be necessary to do this more than once, since import() is called every time step
                }
            }
            
            double rateNow = m_rate[m_lastIndex].value;
            if( rateNow > 0.0 ){
                for(Human& human : population){
                    if(human.rng.bernoulli( rateNow )){
                        human.withinHostModel->importInfection(human.rng);
                        mon::recordEvent(mon::measure("nNewInfections"), human);
                        mon::recordEvent(mon::measure("nNewInfections_Imported"), human);
                    }
                }
            }
        }
        
        /// Checkpointing
        template<class S>
        void operator& (S& stream) {
            using namespace OM::util::checkpoint;
            // period and rate are set from XML and not changed
            m_lastIndex & stream;
        }
        
    private:
        SimTime m_period = Sim::never();
        uint32_t m_lastIndex;
        struct Rate {
            Rate( SimTime t, double v ): time(t), value(v) {}
            SimTime time = Sim::never();
            double value;
            inline bool operator< (const Rate& that) const{
                return time < that.time;
            }
        };
        vector<Rate> m_rate;
    };
} }

#endif
