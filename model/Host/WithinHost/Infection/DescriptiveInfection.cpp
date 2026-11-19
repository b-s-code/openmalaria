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

#include "Host/WithinHost/Infection/DescriptiveInfection.h"
#include "util/random.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"

#include <sstream>
#include <string>
#include <cmath>
#include <fstream>

namespace OM {
namespace WithinHost {
using namespace util;

// static class variables (see description in header file):
double DescriptiveInfection::s_meanLogParasiteCount[NumDurations][NumDurations];
double DescriptiveInfection::s_sigma0sq;
double DescriptiveInfection::s_xNuStar;

bool bugfixMaxDens = true, bugfixInnateMaxDens = true;


// ———  static init/clear ———

void DescriptiveInfection::init (const Parameters& parameters) {
    // Error checks
    if( Sim::oneTS() != 5 ){
        // To support non-5-day time-step models, either different data would
        // be needed or times need to be adjusted when accessing
        // meanLogParasiteCount. Probably the rest would be fine.
        throw util::XmlScenarioError ("DescriptiveInfection only supports using an interval of 5");
    }
    // Bug fixes: these are enabled by default but may be off in old parameterisations
    bugfixInnateMaxDens = util::ModelOptions::option (util::INNATE_MAX_DENS);
    // Warning: if MAX_DENS_CORRECTION is off, infections not yet at the blood
    // stage could result in BSVEfficacy and potentially innateImmSurvFact
    // being applied to timeStepMaxDensity more than once in some cases.
    bugfixMaxDens = util::ModelOptions::option (util::MAX_DENS_CORRECTION);
    
    // Read parameters
    s_sigma0sq=parameters[Parameter::SIGMA0_SQ];
    s_xNuStar=parameters[Parameter::X_NU_STAR];
    
    // Read file empirical parasite densities
    string densitiesFilename = util::CommandLine::lookupResource ("densities.csv");
    ifstream fMTherapyDensities( densitiesFilename.c_str() );
    if( !fMTherapyDensities.good() ){
        throw util::BaseException( string("Cannot read ").append(densitiesFilename), util::Error::FileIO );
    }
    
    //read header of file (unused)
    string csvLine;
    getline(fMTherapyDensities,csvLine);

    // read every line from the stream
    while (getline(fMTherapyDensities, csvLine)) {
        //Empirical description of single Malaria infections in naive individuals
        //counter variables, i stands for 5 day time interval, j for duration of infection
        int i;
        int j;
        double meanlogdens;

        std::istringstream csvStream(csvLine);

        string csvField1, csvField2, csvField3;

        // read every element from the line that is seperated by commas
        getline(csvStream, csvField1, ',');
        getline(csvStream, csvField2, ',');
        getline(csvStream, csvField3, ',');

        istringstream csvNum1(csvField1), csvNum2(csvField2), csvNum3(csvField3);

        csvNum1 >> i;
        csvNum2 >> j;
        csvNum3 >> meanlogdens;
        
        if( i < 1 || i > NumDurations || j < 1 || j > NumDurations ){
            throw util::BaseException( string("error in data file: ")
                .append(densitiesFilename), Error::InputResource );
        }

        //fill initial matrix
        s_meanLogParasiteCount[i-1][j-1]=meanlogdens;
        //fill also the triangle that will not be used (to ensure everything is initialised)
        if (j!=i) {
            s_meanLogParasiteCount[j-1][i-1]=0.0;
        }

    }
}


// ———  non-static init/destruction  ———

DescriptiveInfection::DescriptiveInfection (LocalRng& rng, uint32_t genotype, int origin) :
	Infection(genotype, origin),
        m_duration(infectionDuration(rng)),
        m_notPrintedMdWarning(true)
{
    assert( Sim::oneTS() == 5 );
}

SimTime DescriptiveInfection::infectionDuration(LocalRng& rng) {
    // Forgive the excess precision; it just avoids having to update all expected results
    double durMean = 5.1300001144409179688;
    double durSigma = 0.80000001192092895508;
    double dur=rng.logNormal(durMean, durSigma);
    
    //TODO:
    // Model did say infection is cleared on day dur+1 converted to a time-step
    // let interval = Sim::oneTS() in:
    // ((1+floor(dur))/interval); now it says the last interval is:
    // floor((1+dur)/interval)-1 = floor((dur+1-interval)/interval)
    // Is this reasonable, or should we change?
    return Sim::fromDays( static_cast<int>(1.0 + dur) ) - Sim::oneTS();
}


// ———  time-step updates  ———
void DescriptiveInfection::determineDensities(
        LocalRng& rng,
        double cumulativeh,
        double &timeStepMaxDensity,
        double immSurvFact,
        double innateImmSurvFact,
        double bsvFactor)
{
    // Age of patent blood stage infection. Note: liver stage is fixed at one
    // 5-day time step and prepatent blood stage is latentp - 1 time steps.
    SimTime infage = Sim::ts0() - m_startDate - s_latentP;
    if ( infage < Sim::zero()) {
        m_density = 0.0;
        if (bugfixMaxDens) timeStepMaxDensity = 0.0;
    }else{
        timeStepMaxDensity = 0.0;
        
        int32_t infAge = min( Sim::inSteps(infage), MaxDurationTs );
        int32_t infDur = min( Sim::inSteps(m_duration), MaxDurationTs );
        m_density=max (s_meanLogParasiteCount[infAge][infDur], 0.0);
        
        // The expected parasite density in the non naive host (AJTM p.9 eq. 9)
        // Note that in published and current implementations Dx is zero.
        m_density = m_density * immSurvFact;
        
        //Perturb m_density using a lognormal
        double varlog = s_sigma0sq / (1.0 + (cumulativeh / s_xNuStar));
        double stdlog = sqrt(varlog);
        
        /*
        This code samples from a log normal distribution with mean equal to the predicted density
        n.b. AJTM p.9 eq 9 implies that we sample the log of the density from a normal with mean equal to
        the log of the predicted density.  If we really did the latter then this bias correction is not needed.
        */
        if (stdlog > 0.0000001) {
            // Calculate the expected density on the day of sampling:
            double meanlog = m_density - stdlog*stdlog / 2.0;
            m_density = rng.logNormal(meanlog, stdlog);
            // Calculate additional samples for T-1 days (T=days per step):
            if( true /*Sim::oneTS() > 1, always true for this model*/ ){
                timeStepMaxDensity = rng.maxMultiLogNormal (m_density,
                        Sim::oneTS() - 1, meanlog, stdlog);
            } else {
                timeStepMaxDensity = m_density;
            }
        }
        if (timeStepMaxDensity > MaxDens && m_notPrintedMdWarning){
            cerr << "TSMD hit limit:\t" << m_density << ",\t" << timeStepMaxDensity << endl;
            m_notPrintedMdWarning = false;
        }
        m_density = min(m_density, MaxDens);
        timeStepMaxDensity = min(timeStepMaxDensity, MaxDens);
        
        //Compute the proportion of parasites remaining after innate blood stage effect
        m_density *= innateImmSurvFact;
        
        //Include here the effect of blood stage vaccination
        m_density *= bsvFactor;
        
        m_cumulativeExposureJ += Sim::oneTS() * m_density;
    }
    
    if (bugfixInnateMaxDens) timeStepMaxDensity *= innateImmSurvFact;
    timeStepMaxDensity *= bsvFactor;
}


// ———  checkpointing  ———

DescriptiveInfection::DescriptiveInfection (istream& stream) :
        Infection(stream), m_duration(Sim::never())
{
    m_duration & stream;
    m_notPrintedMdWarning & stream;
}
void DescriptiveInfection::checkpoint (ostream& stream) {
    Infection::checkpoint (stream);
    m_duration & stream;
    m_notPrintedMdWarning & stream;
}

}
}
