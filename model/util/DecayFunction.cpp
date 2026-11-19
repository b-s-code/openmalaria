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

#include <algorithm>
#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "util/UnitParse.h"
#include "util/DecayFunction.h"

namespace OM {
namespace util {
    
inline double readLToDays( const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw XmlScenarioError( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(elt.getL().get(), UnitParse::YEARS);
}

class ConstantDecayFunction : public DecayFunction {
public:
    ConstantDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV())
    {}
    
    double compute(double effectiveAge) const {
        // Note: we now require all decay functions to return 0 when time > 0
        // and the DecayFunction is default-constructed. So const *after deployment*.
        if( effectiveAge * m_hetFactor == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    }
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return Sim::future();        // decay occurs "in the future" (don't use Sim::never() because that is interpreted as being in the past)
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<ConstantDecayFunction> copy = std::make_unique<ConstantDecayFunction>(*this);
        return std::move(copy);
    }
};

class StepDecayFunction : public DecayFunction {
public:
    StepDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        m_invL( 1.0 / readLToDays(elt) )
    {}
    
    double compute(double effectiveAge) const{
        if( effectiveAge * m_invL * m_hetFactor < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return Sim::roundToTSFromDays( 1.0 / m_invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<StepDecayFunction> copy = make_unique<StepDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double m_invL;
};

class LinearDecayFunction : public DecayFunction {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        m_invL( 1.0 / readLToDays(elt) )
    {}
    
    double compute(double effectiveAge) const{
        if( effectiveAge * m_invL * m_hetFactor < 1.0 )
            return 1.0 - effectiveAge * m_invL * m_hetFactor;
        else
            return 0.0;
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // Note: rounds to nearest. Object may decay instantly or at time L.
        return Sim::roundToTSFromDays(rng.uniform01() / m_invL);
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<LinearDecayFunction> copy = make_unique<LinearDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double m_invL;
};

class ExponentialDecayFunction : public DecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        m_invL( log(2.0) / readLToDays(elt) )
    {}
    
    double compute(double effectiveAge) const{
        return exp( -effectiveAge * m_invL * m_hetFactor);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return Sim::roundToTSFromDays( -log(rng.uniform01()) / m_invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<ExponentialDecayFunction> copy = make_unique<ExponentialDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double m_invL;
};

class WeibullDecayFunction : public DecayFunction {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        m_constOverLambda( pow(log(2.0),1.0/elt.getK()) / readLToDays(elt) ),
        m_k( elt.getK() )
    {}
    
    double compute(double effectiveAge) const{
        double p = -pow(effectiveAge * m_constOverLambda * m_hetFactor, m_k);
        if(p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return Sim::roundToTSFromDays( pow( -log(rng.uniform01()), 1.0/m_k ) / m_constOverLambda );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<WeibullDecayFunction> copy = make_unique<WeibullDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double m_constOverLambda, m_k;
};

class HillDecayFunction : public DecayFunction {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        m_invL( 1.0 / readLToDays(elt) ),
        m_k( elt.getK() )
    {}
    
    double compute(double effectiveAge) const{
        return 1.0 / (1.0 + pow(effectiveAge * m_invL * m_hetFactor, m_k));
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return Sim::roundToTSFromDays( pow( 1.0 / rng.uniform01() - 1.0, 1.0/m_k ) / m_invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<HillDecayFunction> copy = make_unique<HillDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double m_invL, m_k;
};

class SmoothCompactDecayFunction : public DecayFunction {
public:
    SmoothCompactDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        m_invL( 1.0 / readLToDays(elt) ),
        m_k( elt.getK() )
    {}

    double compute(double effectiveAge) const{
        if( effectiveAge * m_invL * m_hetFactor < 1.0 )
            return exp( m_k - m_k / (1.0 - pow(effectiveAge * m_invL * m_hetFactor, 2.0)) );
        else
            return 0.0;
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return Sim::roundToTSFromDays( sqrt( 1.0 - m_k / (m_k - log( rng.uniform01() )) ) / m_invL );
    }

    unique_ptr<DecayFunction> makeHetSample(double) const {
        unique_ptr<SmoothCompactDecayFunction> copy = make_unique<SmoothCompactDecayFunction>(*this);
        return std::move(copy);
    }
    
private:
    double m_invL, m_k;
};

template<class T>
class OperatorDecayFunction : public DecayFunction {
public:
    OperatorDecayFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 2)
            throw XmlScenarioError("Operator decay function expects two decay functions, " + to_string(decaySequence.size()) +"  were given.");

        m_f1 = DecayFunction::makeObject(decaySequence[0], "Operator::f1");
        m_f2 = DecayFunction::makeObject(decaySequence[1], "Operator::f2");
    }

    OperatorDecayFunction(const OperatorDecayFunction &copy, unique_ptr<DecayFunction> f1, unique_ptr<DecayFunction> f2) : 
        DecayFunction(copy),
        m_f1(std::move(f1)), m_f2(std::move(f2)) {}

    double compute(double effectiveAge) const {
        return max(min(m_op(m_f1->eval(effectiveAge), m_f2->eval(effectiveAge)), 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return Sim::roundToTSFromDays( max(m_f1->sampleAgeOfDecay(rng), m_f2->sampleAgeOfDecay(rng)) );
    }

    unique_ptr<DecayFunction> makeHetSample(double hetFactor) const {
        unique_ptr<DecayFunction> f1hetSample = m_f1->hetSample(hetFactor);
        unique_ptr<DecayFunction> f2hetSample = m_f2->hetSample(hetFactor);
        unique_ptr<OperatorDecayFunction> copy = make_unique<OperatorDecayFunction>(*this, std::move(f1hetSample), std::move(f2hetSample));
        return std::move(copy);
    }

private:
    unique_ptr<DecayFunction> m_f1, m_f2;
    T m_op;
};

class EmaxFunction : public DecayFunction {
public:
    EmaxFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 1)
            throw util::XmlScenarioError("The Emax function expects exactly one user function, but " + to_string(decaySequence.size()) +"  were given.");
        if(!elt.getEmax().present())
            throw util::XmlScenarioError("Emax function: the Emax parameter is not declared");
        if(!elt.getIC50().present())
            throw util::XmlScenarioError("Emax function: the IC50 parameter is not declared");
        if(!elt.getSlope().present())
            throw util::XmlScenarioError("Emax function: the Slope parameter is not declared");
        if(!elt.getInitialConcentration().present())
            throw util::XmlScenarioError("Emax function: the InitialConcentration parameter is not declared");
        
        m_f = makeObject(decaySequence[0], "Emax::f");

        m_emax = elt.getEmax().get();
        const double rawIC50 = elt.getIC50().get();
        m_slope = elt.getSlope().get();
        m_initialConcentration = elt.getInitialConcentration().get();

        if (m_emax < 0.0 || m_emax > 1.0) {
            throw util::XmlScenarioError(
                "Emax function: Emax must be in [0, 1], got " + std::to_string(m_emax));
        }

        if (rawIC50 <= 0.0) {
            throw util::XmlScenarioError(
                "Emax function: IC50 must be > 0, got " + std::to_string(rawIC50));
        }

        if (m_slope <= 0.0) {
            throw util::XmlScenarioError(
                "Emax function: Slope must be > 0, got " + std::to_string(m_slope));
        }

        if (m_initialConcentration < 0.0) {
            throw util::XmlScenarioError(
                "Emax function: InitialConcentration must be >= 0, got " + std::to_string(m_initialConcentration));
        }

        // Store IC50^slope once since compute() evaluates this repeatedly.
        IC50 = std::pow(rawIC50, m_slope);
    }

    EmaxFunction(const EmaxFunction &copy, std::unique_ptr<DecayFunction> f) :
        DecayFunction(copy),
        m_f(std::move(f)),
        m_emax(copy.m_emax),
        IC50(copy.IC50),
        m_slope(copy.m_slope),
        m_initialConcentration(copy.m_initialConcentration) {}

    double compute(double effectiveAge) const {
        // Emax model: E(C)=Emax*C^slope/(C^slope + IC50^slope), C=C0*f(t).
        const double concentrationPow =
            std::pow(m_initialConcentration * m_f->eval(effectiveAge), m_slope);
        const double efficacy = m_emax * concentrationPow / (concentrationPow + IC50);
        return std::max(std::min(efficacy, 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return Sim::roundToTSFromDays( m_f->sampleAgeOfDecay(rng) );
    }

    std::unique_ptr<DecayFunction> makeHetSample(double hetFactor) const {
        std::unique_ptr<DecayFunction> fhetSample = m_f->hetSample(hetFactor);
        std::unique_ptr<EmaxFunction> copy =
            std::make_unique<EmaxFunction>(*this, std::move(fhetSample));
        return std::move(copy);
    }

private:
    std::unique_ptr<DecayFunction> m_f;
    double m_emax, IC50, m_slope, m_initialConcentration;
};

// -----  interface / static functions  -----
void DecayFunction::checkpointSample(
        const unique_ptr<DecayFunction>& sample, ostream& stream
) const {
    bool present = bool(sample);
    present & stream;
    if (present) sample->m_hetFactor & stream;
}

void DecayFunction::checkpointSample(
        unique_ptr<DecayFunction>& sample, istream& stream
) const {
    bool present;
    present & stream;
    if (present) {
        double factor;
        factor & stream;
        sample = hetSample(factor);
    } else {
        sample.reset();
    }
}

unique_ptr<DecayFunction> DecayFunction::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    if( func == "constant" )
        return make_unique<ConstantDecayFunction>( elt );
    else if( func == "step" )
        return make_unique<StepDecayFunction>( elt );
    else if( func == "linear" )
        return make_unique<LinearDecayFunction>( elt );
    else if( func == "exponential" )
        return make_unique<ExponentialDecayFunction>( elt );
    else if( func == "weibull" )
        return make_unique<WeibullDecayFunction>( elt );
    else if( func == "hill" )
        return make_unique<HillDecayFunction>( elt );
    else if( func == "smooth-compact" )
        return make_unique<SmoothCompactDecayFunction>( elt );
    else if( func == "plus" )
        return make_unique<OperatorDecayFunction<std::plus<double>>>( elt );
    else if( func == "minus" )
        return make_unique<OperatorDecayFunction<std::minus<double>>>( elt );
    else if( func == "divides" )
        return make_unique<OperatorDecayFunction<std::divides<double>>>( elt );
    else if( func == "multiplies" )
        return make_unique<OperatorDecayFunction<std::multiplies<double>>>( elt );
    else if( func == "emax" )
        return make_unique<EmaxFunction>( elt );
    else
        throw XmlScenarioError("decay function " + string(eltName) + " of type " + string(func) + " unrecognized");
}

} }
