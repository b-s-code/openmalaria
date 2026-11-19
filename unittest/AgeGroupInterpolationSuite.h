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

#ifndef Hmod_AgeGroupInterpolationSuite
#define Hmod_AgeGroupInterpolationSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "util/AgeGroupInterpolation.h"

using util::AgeGroupInterpolator;

class AgeGroupInterpolationSuite : public CxxTest::TestSuite
{
public:
    AgeGroupInterpolationSuite () {
        // Set maximum age to 90 years:
        UnittestUtil::initTime( 5 );
        m_agvElt = unique_ptr<scnXml::AgeGroupValues>(new scnXml::AgeGroupValues());
        scnXml::AgeGroupValues::GroupSequence& seq = m_agvElt->getGroup();
        seq.resize( DataLen, scnXml::Group(0.0,0.0) );
        for( size_t i = 0; i < DataLen; ++i ){
            seq[ i ].setLowerbound( StdLbounds[ i ] );
            seq[ i ].setValue( StdValues[ i ] );
        }
    }
    ~AgeGroupInterpolationSuite () {
    }
    
    void setUp () {
    }
    void tearDown () {
    }
    
    void testDummy() {
        AgeGroupInterpolator o;
        TS_ASSERT_THROWS( o.eval(5.7), const std::exception& );
    }
    
    void testPiecewiseConst () {
        m_agvElt->setInterpolation( "none" );
        AgeGroupInterpolator o;
        o.set( *m_agvElt, "testPiecewiseConst" );
        for( size_t i = 0; i < TestLen; ++i ){
            TS_ASSERT_APPROX( o.eval( TestAges[ i ] ), PiecewiseConstValues[ i ] );
        }
    }
    
    void testLinearInterp () {
        m_agvElt->setInterpolation( "linear" );
        AgeGroupInterpolator o;
        o.set( *m_agvElt, "testLinearInterp" );
        for( size_t i = 0; i < TestLen; ++i ){
            TS_ASSERT_APPROX( o.eval( TestAges[ i ] ), LinearInterpValues[ i ] );
        }
    }
    
private:
    static const size_t DataLen = 5;
    static const size_t TestLen = 8;
    static const double StdLbounds[DataLen];
    static const double StdValues[DataLen];
    static const double TestAges[TestLen];
    static const double PiecewiseConstValues[TestLen];
    static const double LinearInterpValues[TestLen];
    unique_ptr<scnXml::AgeGroupValues> m_agvElt;
};

const double AgeGroupInterpolationSuite::StdLbounds[] = {
    0,5,10,15,60
};
const double AgeGroupInterpolationSuite::StdValues[] = {
    6.08,3.81,2.62,4.05,5.41
};
const double AgeGroupInterpolationSuite::TestAges[] = {
    // various ages, designed to test limits, boundary points and interpolation
    15.2,18.09,7.0,2.5,
    0.0,20.0,900.0,62.0
};
const double AgeGroupInterpolationSuite::PiecewiseConstValues[] = {
    4.0499999999999998, 4.0499999999999998, 3.8100000000000001, 6.0800000000000001,
    6.0800000000000001, 4.0499999999999998, 5.4100000000000001, 5.4100000000000001
};
const double AgeGroupInterpolationSuite::LinearInterpValues[] = {
    2.7744400000000002, 2.9397479999999998, 4.0369999999999999, 6.0800000000000001,
    6.0800000000000001, 3.0489999999999999, 5.4100000000000001, 4.938533333333333
};

#endif
