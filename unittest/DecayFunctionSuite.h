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

#ifndef Hmod_DecayFunctionSuite
#define Hmod_DecayFunctionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "util/DecayFunction.h"

using ::OM::util::DecayFunction;

class DecayFunctionSuite : public CxxTest::TestSuite
{
public:
    DecayFunctionSuite () :
        m_rng(0, 0), m_dfElt( "" )
    {
        m_dfElt.setL( "10y" );
        m_dfElt.setK( 1.6 );
    }
    
    void setUp() {
        m_rng.seed(0, 721347520444481703);
        UnittestUtil::initTime(5);
    }
    
    void testBad() {
        m_dfElt.setFunction( "unknown" );
        string errMsg = "decay function DecayFunctionSuite of type unknown unrecognized";
        TS_ASSERT_THROWS_EQUALS( m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" ),
                          const std::runtime_error &e, e.what(), errMsg );
    }
    
    void testConstant () {
        m_dfElt.setFunction( "constant" );
        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        // Second test: with an appropriately sampled helper value, we should get the results we want
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ) , 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(10 ) ), 1.0 );
        // Third test: time of decay (age plus now) should always be in the future
        TS_ASSERT_EQUALS( m_df->sampleAgeOfDecay(m_rng), Sim::future() );
    }
    
    void testStep () {
        m_dfElt.setFunction( "step" );
        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(6 ) ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(20 ) ), 0.0 );
        TS_ASSERT_EQUALS( m_df->sampleAgeOfDecay(m_rng), Sim::fromYearsI(10) );
    }
    
    void testLinear () {
        m_dfElt.setFunction( "linear" );
        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(6 ) ), 0.4 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(20 ) ), 0.0 );
    }
    
    void testExponential () {
        m_dfElt.setFunction( "exponential" );
        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(6 ) ), 0.65975394736842108 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(20 ) ), 0.25 );
    }
    
    void testWeibull () {
        m_dfElt.setFunction( "weibull" );
        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(6 ) ), 0.73631084210526321 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(20 ) ), 0.122306 );
    }
    
    void testHill () {
        m_dfElt.setFunction( "hill" );
        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(6 ) ), 0.6936673684210527 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(20 ) ), 0.24805074736842106 );
    }
    
    void testSmoothCompact () {
        m_dfElt.setFunction( "smooth-compact" );
        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ), 1.0 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(6 ) ), 0.40656965789473687 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(20 ) ), 0.0 );
    }

    void testEmax () {
        m_dfElt.setFunction( "emax" );
        m_dfElt.setEmax( 0.8 );
        m_dfElt.setIC50( 1.0 );
        m_dfElt.setSlope( 1.0 );
        m_dfElt.setInitialConcentration( 2.0 );
        m_dfElt.getDecay().clear();
        scnXml::DecayFunction inner( "exponential" );
        inner.setL( "10y" );
        m_dfElt.getDecay().push_back( inner );

        m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" );
        unique_ptr<DecayFunction> dHet = m_df->hetSample(m_rng);
        TS_ASSERT_APPROX( dHet->eval( Sim::zero() ), 0.5333333333333333 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(10) ), 0.4 );
        TS_ASSERT_APPROX( dHet->eval( Sim::fromYearsI(20) ), 0.26666666666666666 );
    }

    void testEmaxBadIC50 () {
        m_dfElt.setFunction( "emax" );
        m_dfElt.setEmax( 0.8 );
        m_dfElt.setIC50( -1.0 );
        m_dfElt.setSlope( 0.5 );
        m_dfElt.setInitialConcentration( 2.0 );
        m_dfElt.getDecay().clear();
        scnXml::DecayFunction inner( "constant" );
        m_dfElt.getDecay().push_back( inner );

        string errMsg = "Emax function: IC50 must be > 0, got -1.000000";
        TS_ASSERT_THROWS_EQUALS( m_df = DecayFunction::makeObject( m_dfElt, "DecayFunctionSuite" ),
                          const std::runtime_error &e, e.what(), errMsg );
    }
    
private:
    LocalRng m_rng;
    scnXml::DecayFunction m_dfElt;
    unique_ptr<DecayFunction> m_df;
};

#endif
