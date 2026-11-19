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

#include "util/sampler.h"
#include "util/errors.h"
#include "util/random.h"
#include <cmath>

namespace OM { namespace util {

double NormalSample::asNormal( double mu, double sigma )const{
    return sigma*m_x + mu;
}

double NormalSample::asLognormal( double mu, double sigma )const{
    return exp( sigma*m_x + mu );
}

NormalSample NormalSample::generate(LocalRng& rng) {
    return NormalSample( rng.gauss(0.0, 1.0) );
}

NormalSample NormalSample::generateCorrelated(NormalSample base, double correlation, double factor, LocalRng& rng) {
    if( correlation == 1.0 ) { return base; }
    
    double e = rng.gauss(0.0, factor);
    return NormalSample( base.m_x * correlation + e );
}

void NormalSampler::setParams( double m, double s ){
    m_mu = m;
    m_sigma = s;
}

void NormalSampler::setParams(const scnXml::SampledValueN& elt){
    m_mu = elt.getMean();
    if( elt.getDistr() == "const" ){
        if( elt.getSD().present() && elt.getSD().get() != 0.0) {
            throw util::XmlScenarioError( "attribute SD must be zero or omitted when distr=\"const\" or is omitted" );
        }
        m_sigma = 0.0;
        return;
    }
    
    if( elt.getDistr() != "normal" ){
        throw util::XmlScenarioError( "expected distr to be one of \"const\", \"lognormal\" (note: not all distributions are supported here)" );
    }
    if( !elt.getSD().present() ){
        throw util::XmlScenarioError( "attribute \"SD\" required for sampled normal value when distr is not const" );
    }
    m_sigma = elt.getSD().get();
}

double NormalSampler::sample(LocalRng& rng) const{
    if( m_sigma == 0.0 ){
        return m_mu;
    }
    return rng.gauss( m_mu, m_sigma );
}

unique_ptr<util::LognormalSampler> LognormalSampler::fromMeanCV( double mean, double CV, std::optional<double> truncate )
{
    unique_ptr<LognormalSampler> sampler = std::make_unique<LognormalSampler>();
    sampler->m_truncate = truncate;

    // The distribution is "const"
    if( CV == 0.0 )
    {
        sampler->m_sigma = 0.0;
        // as a special case, we can support mean == CV == 0
        if( mean == 0.0 )
            sampler->m_mu = -numeric_limits<double>::infinity();
        else if( mean > 0 )
            sampler->m_mu = log(mean);
        else
            throw util::XmlScenarioError( "const: required mean > 0" );
        return sampler;
    }

    // The distribution is not "const"
    if( mean <= 0 )
        throw util::XmlScenarioError( "log-normal: required mean > 0" );
    if( CV < 0 )
        throw util::XmlScenarioError( "log-normal: required CV >= 0" );
    
    // The following formulae are derived from:
    // X ~ lognormal(μ, σ)
    // E(X) = exp(μ + σ² / 2)
    // Var(X) = exp(σ² + 2μ)(exp(σ²) - 1)
    // Var = (CV * mean)²
    const double a = 1 + (CV * CV);

    sampler->CV = CV;
    sampler->m_mu = log( mean / sqrt(a) );
    sampler->m_sigma = sqrt(log(a));

    return sampler;
}

unique_ptr<util::LognormalSampler> LognormalSampler::fromMeanVariance( double mean, double variance, std::optional<double> truncate ){
    unique_ptr<LognormalSampler> sampler = std::make_unique<LognormalSampler>();
    sampler->m_truncate = truncate;

    // The distribution is "const"
    if( variance == 0.0)
    {
        sampler->m_sigma = 0.0;
        if(mean > 0)
            sampler->m_mu = log(mean);
        else
            throw util::XmlScenarioError( "const: required mean > 0" );
        return sampler;
    }

    // The distribution is not "const"
    if( mean <= 0 )
        throw util::XmlScenarioError( "log-normal: required mean > 0" );
    if( variance < 0 )
        throw util::XmlScenarioError( "log-normal: required variance >= 0" );

    // The following formulae are derived from:
    // X ~ lognormal(μ, σ)
    // E(X) = exp(μ + σ² / 2)
    // Var(X) = exp(σ² + 2μ)(exp(σ²) - 1)
    // Var = (CV * mean)²
    const double CV = variance / mean;
    const double a = 1 + (CV * CV);
    sampler->m_mu = log( mean / sqrt(a) );
    sampler->m_sigma = sqrt(log(a));

    return sampler;
}

void LognormalSampler::scaleMean(double scalar){
    m_mu += log(scalar);
}

double LognormalSampler::mean() const{
    return exp(m_mu + 0.5*m_sigma*m_sigma);
}

double LognormalSampler::sample(LocalRng& rng) const{
     double value;

    if( m_sigma == 0.0 )
        value = exp( m_mu );
    else
        value = rng.logNormal( m_mu, m_sigma );

    if (m_truncate && value > *m_truncate)
        return *m_truncate;

    return value;
}

double LognormalSampler::cdf(double x) const {
    // 1) anything at or below 0 has CDF = 0
    if (x <= 0.0)
        return 0.0;

    // 2) Degenerate case
    if (m_sigma == 0.0)
    {
        if(log(x) >= m_mu) return 1.0;
        else return 0.0;
    }
    return gsl_cdf_lognormal_P(x, m_mu, m_sigma);
}

unique_ptr<util::GammaSampler> GammaSampler::fromMeanCV( double mean, double CV, std::optional<double> truncate )
{
    if( mean <= 0 )
        throw util::XmlScenarioError( "gamma: required mean > 0" );
    if( CV < 0 )
        throw util::XmlScenarioError( "gamma: required CV >= 0" );

    unique_ptr<GammaSampler> sampler = std::make_unique<GammaSampler>();
    sampler->m_truncate = truncate;
    sampler->m_mu = mean;
    sampler->CV = CV;

    if( CV == 0.0 )
        return sampler;

    // 1 / sqrt(k) = CV
    // sqrt(k) = 1/CV
    // k = 1 / CV^2
    sampler->m_k = 1.0/(CV*CV);
    // k * theta = mean
    // theta = mean / k
    sampler->m_theta = sampler->m_mu / sampler->m_k;

    return sampler;
}

unique_ptr<util::GammaSampler> GammaSampler::fromMeanVariance( double mean, double variance, std::optional<double> truncate )
{
    if( mean <= 0 )
        throw util::XmlScenarioError( "gamma: required mean > 0" );
    if( variance < 0 )
        throw util::XmlScenarioError( "gamma: required variance >= 0" );

    unique_ptr<GammaSampler> sampler = std::make_unique<GammaSampler>();
    sampler->m_truncate = truncate;
    sampler->m_mu = mean;
    sampler->m_variance = variance;

    if( variance == 0.0 )
        return sampler;

    // sigma / mu = 1 / sqrt(k)
    // sqrt(k) = mu / sigma
    // k = mu^2 / variance
    sampler->m_k = (sampler->m_mu*sampler->m_mu)/sampler->m_variance;
    // k * theta = mean
    // theta = mean / k
    sampler->m_theta = sampler->m_mu / sampler->m_k;

    return sampler;
}

double GammaSampler::mean() const {
    return m_mu;
}

double GammaSampler::sample(LocalRng& rng) const{
    double value;

    if (std::isnan(m_theta)) // CV=0
        value = m_mu;
    else
        value = rng.gamma(m_k, m_theta);

    if (m_truncate && value > *m_truncate)
        return *m_truncate;

    return value;
}

double GammaSampler::cdf(double x) const {
    // 1) anything at or below 0 has CDF = 0
    if (x <= 0.0)
        return 0.0;

    // 2) Degenerate case
    if(isnan(m_theta)) // CV=0 or variance=0
    {
        if(x >= m_mu) return 1.0;
        else return 0.0;
    }

    return gsl_cdf_gamma_P(x, m_k, m_theta);
}
        
void BetaSampler::setParamsMV( double mean, double variance ){
    if( variance > 0.0 ){
        // double c = mean / (1.0 - mean);
        // double cp1 = c + 1.0;
        // double s = cp1*cp1*variance;
        // b = (s - c) / (s*cp1);
        // a = c*b;
        
        if (variance >= mean * (1.0 - mean))
            throw util::XmlScenarioError("BetaSampler::setParamsMV: require variance < mean(1.0 - mean)");

        double c = mean * ((1.0 - mean) / variance) - 1.0;
        m_a = mean * c;
        m_b = c - m_a;
    }else if (variance == 0.0)
    {
        // Not using the Beta distirbution, we simply return the mean when sampling
        if( variance != 0.0 ){
            throw util::XmlScenarioError("BetaSampler::setParamsMV: require variance ≥ 0");
        }
        m_a = mean;
        m_b = 0.0;
    }
    else
        throw util::XmlScenarioError("BetaSampler::setParamsMV: require variance ≥ 0");
}

double BetaSampler::sample(LocalRng& rng) const{
    if( m_b == 0.0 ){
        return m_a;
    }else{
        assert( m_a>0.0 && m_b>0.0 );
        return rng.beta( m_a, m_b );
    }
}


} }
