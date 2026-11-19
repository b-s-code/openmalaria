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

// Note: this file is included by exactly one source (InterventionManager.cpp)
// and contains definitions as well as declarations.

// The includes here are more for documentation than required.
#include "Global.h"
#include "util/errors.h"
#include "Clinical/ClinicalModel.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/VectorModel.h"
#include "Transmission/PerHost.h"
#include "util/random.h"
#include "util/errors.h"
#include <schema/interventions.h>

#include "util/SpeciesIndexChecker.h"
#include "Host/Human.h"

namespace OM { namespace interventions {

struct ByDeployTime;    // forward decl for friend

// ———  TimedDeployment and derivatives  ———

/** Interface for timed deployment of an intervention. */
class TimedDeployment {
public:
    /// Create, passing time of deployment
    explicit TimedDeployment(SimTime deployDate) :
            date( deployDate )
    {
        if( deployDate < Sim::startDate() ){
            throw util::XmlScenarioError("timed intervention deployment: may not be negative");
        }else if( deployDate >= Sim::endDate() ){
            cerr << "Warning: timed intervention deployment at date "<< deployDate;
            cerr << " happens after last survey" << endl;
        }
    }
    virtual ~TimedDeployment() {}
    
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) =0;
    
    virtual void printDetails( std::ostream& out )const =0;
    
    // Read access required in this file; don't really need protection:
    SimTime date = Sim::never();
};

class DummyTimedDeployment : public TimedDeployment {
public:
    DummyTimedDeployment() :
        TimedDeployment( Sim::startDate() /* hack */ )
    {
        // TimedDeployment's ctor checks that the deployment time-step is
        // within the intervention period. We want this time to be after the
        // last time-step, so set the time here after TimedDeployment's ctor
        // check has been done (hacky).
        date = Sim::future();
    }
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {}
    virtual void printDetails( std::ostream& out )const{
        out << "Dummy";
    }
};

class TimedChangeHSDeployment : public TimedDeployment {
public:
    TimedChangeHSDeployment( SimTime date, const scnXml::ChangeHS::TimedDeploymentType& hs ) :
        TimedDeployment( date ),
        m_newHs( hs._clone() )
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        Clinical::ClinicalModel::setHS( *m_newHs );
        delete m_newHs;
        m_newHs = 0;
    }
    virtual void printDetails( std::ostream& out )const{
        out << "ChangeHS";
    }
    
private:
    scnXml::HealthSystem *m_newHs;
};

class TimedChangeEIRDeployment : public TimedDeployment {
public:
    TimedChangeEIRDeployment( SimTime date, const scnXml::ChangeEIR::TimedDeploymentType& nv ) :
        TimedDeployment( date ),
        m_newEir( nv._clone() )
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        transmission.changeEIRIntervention( *m_newEir );
        delete m_newEir;
        m_newEir = 0;
    }
    virtual void printDetails( std::ostream& out )const{
        out << "ChangeEIR";
    }
    
private:
    scnXml::NonVector *m_newEir;
};

class TimedUninfectVectorsDeployment : public TimedDeployment {
public:
    TimedUninfectVectorsDeployment( SimTime date ) :
        TimedDeployment( date )
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        transmission.uninfectVectors();
    }
    virtual void printDetails( std::ostream& out )const{
        out << "UninfectVectors";
    }
};

void VaccineLimits::set( const scnXml::DeploymentBase& deploy ){
    if( deploy.getVaccMinPrevDoses().present() ){
        int v = deploy.getVaccMinPrevDoses().get();
        if( v < 0 ) throw util::XmlScenarioError( "vaccMinPrevDoses: min value is 0" );
        minPrevDoses = v;
    }
    if( deploy.getVaccMaxCumDoses().present() ){
        int v = deploy.getVaccMaxCumDoses().get();
        if( v < 0 ) throw util::XmlScenarioError( "vaccMaxCumDoses: min value is 0" );
        maxCumDoses = v;
    }
}

/// Base class for TimedHumanDeployment and ContinuousHumanDeployment
class HumanDeploymentBase {
protected:
    /**
     * @param deploy XML element describing deployment
     * @param intervention The intervention to deploy (list of components)
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     * @param complement Whether to take the complement of the sub-population
     *  to which deployment will be restricted
     */
    HumanDeploymentBase( const scnXml::DeploymentBase& deploy,
                         shared_ptr<const HumanIntervention> intervention,
                         ComponentId subPop, bool complement ) :
            m_coverage( deploy.getCoverage() ),
            m_subPop( subPop ),
            m_complement( complement ),
            m_intervention( intervention )
    {
        if( !(m_coverage >= 0.0 && m_coverage <= 1.0) ){
            throw util::XmlScenarioError("intervention deployment coverage must be in range [0,1]");
        }
        // If coverageCorr or coverageVar is present, then the user is asking for the copula
        if(deploy.getCoverageCorr().present() || deploy.getCoverageVar().present())
        {
            if(!(deploy.getCoverageCorr().present() && deploy.getCoverageVar().present()))
                throw std::runtime_error("[Gaussian copula]: both coverageCorr and coverageVar must be present");

            if(deploy.getCoverageVar().get() < 0)
                throw std::runtime_error("[Gaussian copula]: negative variance is not allowed");
            
            // if variance is 0, then no copula and fall back to probability = coverage
            if(deploy.getCoverageVar().get() > 0)
            {
                m_coverageCorr = deploy.getCoverageCorr().get();
                m_coverageVar = deploy.getCoverageVar().get();
                m_copula = true;
            }
        }

        m_vaccLimits.set( deploy );
    }
    
    inline void deployToHuman( Host::Human& human, mon::Deploy::Method method ) const{
        m_intervention->deploy( human, method, m_vaccLimits );
    }
    
    bool m_copula = false;
    double m_coverage;    // proportion coverage within group meeting above restrictions
    double m_coverageCorr = 0.0, m_coverageVar = 0.0;
    VaccineLimits m_vaccLimits;
    ComponentId m_subPop;      // ComponentId::wholePop() if deployment is not restricted to a sub-population
    bool m_complement;
    shared_ptr<const HumanIntervention> m_intervention;
};


/// Timed deployment of human-specific interventions
class TimedHumanDeployment : public TimedDeployment, protected HumanDeploymentBase {
public:
    /** 
     * @param date Time of deployment
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     */
    TimedHumanDeployment( SimTime date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement ) :
        TimedDeployment( date ),
        HumanDeploymentBase( mass, intervention, subPop, complement ),
        m_minAge( Sim::fromYearsN( mass.getMinAge() ) ),
        m_maxAge( Sim::future() ),
        m_minAvailability( 0 ),
        m_maxAvailability( 100 )
    {
        if( mass.getMaxAge().present() )
            m_maxAge = Sim::fromYearsN( mass.getMaxAge().get() );
            
        if( m_minAge < Sim::zero() || m_maxAge < m_minAge )
            throw util::XmlScenarioError("timed intervention must have 0 <= minAvailability <= maxAvailability <= 100");

        if( mass.getMaxAvailability().present() )
            m_maxAvailability = mass.getMaxAvailability().get();

        if( mass.getMinAvailability().present() )
            m_minAvailability = mass.getMinAvailability().get();
            
        if( m_minAvailability < Sim::zero() || m_maxAvailability < m_minAvailability )
            throw util::XmlScenarioError("timed intervention must have 0 <= minAvailability <= maxAvailability <= 100");
    }
    
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        for(Human& human : population) {
            const SimTime age = human.age(Sim::now());
            if( age >= m_minAge && age < m_maxAge ){
                double availability = 0.0;
                for(size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(); ++i)
                    availability += human.perHostTransmission.anophEntoAvailability[i];

                double probability = m_coverage;

                if(m_copula && m_coverage > 0)
                {
                    try 
                    {
                        // Beta distribution for intervention (Beta distributed)
                        const double betaMean = m_coverage;  // Assuming this value is given
                        const double betaVar = m_coverageVar;       // Given as well
                        const double alpha = ((1.0 - betaMean) / betaVar - 1.0 / betaMean) * (betaMean * betaMean);
                        const double beta = alpha * (1.0 / betaMean - 1.0);

                        if(Transmission::PerHostAnophParams::numSpecies() > 1)
                            throw std::runtime_error("only supports one mosquito species");

                        const double ux = Transmission::PerHostAnophParams::get(0).entoAvailability->cdf(human.perHostTransmission.anophEntoAvailabilityRaw[0]);
                        if(ux < 1)
                        {
                            double xx = gsl_cdf_ugaussian_Pinv(ux);                      // Unit interval to Normal

                            // Apply the Gaussian copula transformation for correlation
                            const ComponentId cid = m_subPop;

                            /** human.rng.gauss(0.0, 1.0) is calculated once per component and per human. Re-deployment of
                             * the same component on the same human must use the same gaussian sample. Therefore we store
                             * this value in the human the first time it is calculated. */
                            double g = 0.0;
                            const auto it = human.perHostTransmission.copulaGaussianSamples.find(cid);
                            if (it != human.perHostTransmission.copulaGaussianSamples.end())
                                g = it->second;
                            else
                            {
                                g = human.rng.gauss(0.0, 1.0);
                                human.perHostTransmission.copulaGaussianSamples[cid] = g;
                            }
                            const double yy = m_coverageCorr * xx + g * sqrt(1 - m_coverageCorr * m_coverageCorr);

                            // Transform back to unit interval
                            const double uy = gsl_cdf_ugaussian_P(yy);

                            if (alpha < 0 || beta < 0)
                                throw std::runtime_error("resulting alpha and beta parameters must be positive");

                            // Get the final probability from Beta distribution
                            probability = gsl_cdf_beta_Pinv(uy, alpha, beta);
                        }
                        else // Gamma with CV=0
                            probability = human.rng.beta(alpha, beta);
                    }
                    catch (const std::exception& e) {
                        std::ostringstream oss;
                        oss << "[Gaussian copula]: computing correlated intervention probability using Gaussian copula and Beta distribution: "
                            << e.what()
                            << ". Possible reason: coverage might be too high and/or population size might be too small.";
                        throw std::runtime_error(oss.str());
                    }
                }

                if( availability >= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(m_minAvailability) && availability <= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(m_maxAvailability) ) {
                    if( m_subPop == ComponentId::wholePop() || (human.isInSubPop( m_subPop ) != m_complement) ){
                        if( human.rng.bernoulli( probability ) ){
                            deployToHuman( human, mon::Deploy::TIMED );
                        }
                    }
                }
            }
        }
    }
    
    virtual void printDetails( std::ostream& out )const{
        out << "Human:" << endl;
        out << "\tage: " << Sim::inYears(m_minAge) << "y\tmax age: " << Sim::inYears(m_maxAge) << "y\tsubpop: ";
        if( m_subPop == ComponentId::wholePop() ) out << "whole";
        else out << m_subPop.id;
        out << "\tcomplement: " << m_complement << "\tcoverage: " << m_coverage << "\t";
        m_intervention->printDetails( out );
    }
    
protected:
    // restrictions on deployment
    SimTime m_minAge, m_maxAge;
    SimTime m_minAvailability, m_maxAvailability;
};

/// Timed deployment of human-specific interventions in cumulative mode
class TimedCumulativeHumanDeployment : public TimedHumanDeployment {
public:
    /** 
     * @param date Time of deployment
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param subPop Either ComponentId::wholePop() or a sub-population to which deployment is restricted
     * @param cumCuvId Id of component to test coverage for
     */
    TimedCumulativeHumanDeployment( SimTime date,
                           const scnXml::MassDeployment& mass,
                           shared_ptr<const HumanIntervention> intervention,
                           ComponentId subPop, bool complement,
                           ComponentId cumCuvId ) :
        TimedHumanDeployment( date, mass, intervention, subPop, complement ),
        m_cumCovInd( cumCuvId )
    {
    }
    
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        // Cumulative case: bring target group's coverage up to target coverage
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally subPop
        for(Host::Human &human : population) {
            SimTime age = human.age(Sim::now());
            if( age >= m_minAge && age < m_maxAge ){
                double availability = 0.0;
                for(size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(); ++i)
                    availability += human.perHostTransmission.anophEntoAvailability[i];
                if( availability >= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(m_minAvailability) && availability <= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(m_maxAvailability) ) {
                    if( m_subPop == ComponentId::wholePop() || (human.isInSubPop( m_subPop ) != m_complement) ){
                        total+=1;
                        if( !human.isInSubPop(m_cumCovInd) )
                            unprotected.push_back( &human );
                    }
                }
            }
        }
        
        if( total == 0 ) return;        // no humans to deploy to; avoid divide by zero
        double propProtected = static_cast<double>( total - unprotected.size() ) / static_cast<double>( total );
        if( propProtected < m_coverage ){
            // Proportion propProtected are already covered, so need to
            // additionally cover the proportion (coverage - propProtected),
            // selected from the list unprotected.
            double additionalCoverage = (m_coverage - propProtected) / (1.0 - propProtected);
            cerr << "cum deployment: prop protected " << propProtected << "; additionalCoverage " << additionalCoverage << "; total " << total << endl;
            for(Human* human : unprotected) {
                if( human->rng.uniform01() < additionalCoverage ){
                    deployToHuman( *human, mon::Deploy::TIMED );
                }
            }
        }
    }
    
protected:
    ComponentId m_cumCovInd;
};

class TimedVectorDeployment : public TimedDeployment {
public:
    TimedVectorDeployment( SimTime date, size_t instance ) :
        TimedDeployment( date ),
        m_inst(instance)
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
            vectorModel->deployVectorPopInterv(m_inst);
    }
    virtual void printDetails( std::ostream& out )const{
        out << "Vector";
    }
    
private:
    size_t m_inst;
};

class TimedTrapDeployment : public TimedDeployment {
public:
    TimedTrapDeployment( SimTime date, size_t instance, double ratio, SimTime lifespan ) :
        TimedDeployment(date), m_inst(instance), m_ratio(ratio), m_lifespan(lifespan)
    {}
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission) {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            double number = population.size() * m_ratio;
            vectorModel->deployVectorTrap(m_inst, number, m_lifespan);
        }
    }
    virtual void printDetails( std::ostream& out )const{
        out << "VectorTrap:" << endl;
        out << "\tratio: " << m_ratio << "\tlifespan: " << m_lifespan; 
    }
    
private:
    size_t m_inst;
    double m_ratio;
    SimTime m_lifespan = Sim::never();
};

class TimedNonHumanHostsDeployment : public TimedDeployment {
public:
    TimedNonHumanHostsDeployment( SimTime date, size_t instance, string intervName, const scnXml::Description2::AnophelesSequence list, const scnXml::DecayFunction &decay, Transmission::TransmissionModel& transmission) :
        TimedDeployment( date ),
        m_instance(instance),
        m_intervName(intervName)
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            util::SpeciesIndexChecker checker(intervName, vectorModel->speciesIndex);
            for (const scnXml::NonHumanHostsSpeciesIntervention &anoph : list)
            {
                const string &mosqName = anoph.getMosquito();
                size_t i = checker.getIndex(mosqName);

                Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();

                if (anophModel->reduceNhhAvailability[intervName].size() <= instance) anophModel->reduceNhhAvailability[intervName].resize(instance + 1);
                if (anophModel->reduceP_B_I[intervName].size() <= instance) anophModel->reduceP_B_I[intervName].resize(instance + 1);
                if (anophModel->reduceP_C_I[intervName].size() <= instance) anophModel->reduceP_C_I[intervName].resize(instance + 1);
                if (anophModel->reduceP_D_I[intervName].size() <= instance) anophModel->reduceP_D_I[intervName].resize(instance + 1);
                if (anophModel->reduceFecundity[intervName].size() <= instance) anophModel->reduceFecundity[intervName].resize(instance + 1);

                if (anoph.getAvailabilityReduction().present())
                {
                    const scnXml::AvailabilityReduction &decayFunc = anoph.getAvailabilityReduction().get();
                    if (decayFunc.getInitial() > 1.0) throw util::XmlScenarioError("availabilityReduction intervention: initial effect must be <= 1");
                    anophModel->reduceNhhAvailability[intervName][instance].set(decayFunc.getInitial(), decay, "availabilityReduction");
                }
                if (anoph.getPreprandialKillingEffect().present())
                {
                    const scnXml::PreprandialKillingEffect &decayFunc = anoph.getPreprandialKillingEffect().get();
                    if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                        throw util::XmlScenarioError("PreprandialKillingEffect intervention: initial effect must be between 0 and 1");
                    anophModel->reduceP_B_I[intervName][instance].set(decayFunc.getInitial(), decay, "reduceP_B_I");
                }
                if (anoph.getPostprandialKillingEffect().present())
                {
                    const scnXml::PostprandialKillingEffect &decayFunc = anoph.getPostprandialKillingEffect().get();
                    if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                        throw util::XmlScenarioError("PostprandialKillingEffect intervention: initial effect must be between 0 and 1");
                    anophModel->reduceP_C_I[intervName][instance].set(decayFunc.getInitial(), decay, "reduceP_C_I");
                }
                if (anoph.getRestingKillingEffect().present())
                {
                    const scnXml::RestingKillingEffect &decayFunc = anoph.getRestingKillingEffect().get();
                    if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                        throw util::XmlScenarioError("RestingKillingEffect intervention: initial effect must be be between 0 and 1");
                    anophModel->reduceP_D_I[intervName][instance].set(decayFunc.getInitial(), decay, "reduceP_D_I");
                }
                if (anoph.getFecundityReduction().present())
                {
                    const scnXml::FecundityReduction &decayFunc = anoph.getFecundityReduction().get();
                    if (decayFunc.getInitial() < 0 || decayFunc.getInitial() > 1)
                        throw util::XmlScenarioError("FecundityReduction intervention: initial effect must be be between 0 and 1");
                    anophModel->reduceFecundity[intervName][instance].set(decayFunc.getInitial(), decay, "reduceFecundity");
                }
            }
            checker.checkNoneMissed();
        }
    }
    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission)
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i)
            {
                Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();

                if (anophModel->nhhInstances.count(m_intervName) == 0)
                    throw util::XmlScenarioError("non human hosts type " + m_intervName + " not deployed during non human hosts intervention deployment");

                anophModel->reduceNhhAvailability[m_intervName][m_instance].deploy(vectorModel->rng, Sim::now());
                anophModel->reduceP_B_I[m_intervName][m_instance].deploy(vectorModel->rng, Sim::now());
                anophModel->reduceP_C_I[m_intervName][m_instance].deploy(vectorModel->rng, Sim::now());
                anophModel->reduceP_D_I[m_intervName][m_instance].deploy(vectorModel->rng, Sim::now());
                anophModel->reduceFecundity[m_intervName][m_instance].deploy(vectorModel->rng, Sim::now());
            }
        }
    }
    virtual void printDetails( std::ostream& out )const{
        out << "NonHumanHosts(" << m_intervName << ")";
    }
    
private:
    size_t m_instance;
    string m_intervName;
};

const string VecModeErr = "vector interventions can only be used in dynamic transmission mode (mode=\"dynamic\")";

class TimedAddNonHumanHostsDeployment : public TimedDeployment {
public:
    TimedAddNonHumanHostsDeployment( SimTime date, const string &intervName, SimTime lifespan, const scnXml::Description3::AnophelesSequence list, Transmission::TransmissionModel& transmission) :
        TimedDeployment( date ), m_intervName(intervName), m_lifespan(lifespan)
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            util::SpeciesIndexChecker checker(intervName, vectorModel->speciesIndex);
            for (const scnXml::NonHumanHostsVectorSpecies &anoph : list)
            {
                const string &mosqName = anoph.getMosquito();
                checker.getIndex(mosqName);

                NhhParamsInterv nhh;
                nhh.mosqRelativeAvailabilityHuman = anoph.getMosqRelativeAvailabilityHuman().getValue();
                nhh.mosqProbBiting = anoph.getMosqProbBiting().getValue();
                nhh.mosqProbFindingRestSite = anoph.getMosqProbFindRestSite().getValue();
                nhh.mosqProbResting = anoph.getMosqProbResting().getValue();
                nhh.hostFecundityFactor = anoph.getHostFecundityFactor().getValue();
                m_nhhParams[mosqName] = nhh;
            }
            checker.checkNoneMissed();
        }
    }

    virtual void deploy (vector<Host::Human> &population, Transmission::TransmissionModel& transmission)
    {
        Transmission::VectorModel *vectorModel = dynamic_cast<Transmission::VectorModel *>(&transmission);
        if(vectorModel)
        {
            double popSize = population.size();
            if (vectorModel->interventionMode != Transmission::SimulationMode::dynamicEIR) { throw util::XmlScenarioError(VecModeErr); }
            for (size_t i = 0; i < vectorModel->speciesIndex.size(); ++i)
            {
                Transmission::Anopheles::AnophelesModel *anophModel = vectorModel->species[i].get();
        
                if (anophModel->nhhInstances.count(m_intervName) != 0)
                    throw util::XmlScenarioError("non human hosts type " + m_intervName + " already deployed during non human hosts deployment");

                const NhhParamsInterv &p = m_nhhParams[anophModel->mosq.name];
                const double adultAvail = Transmission::PerHostAnophParams::get(i).entoAvailability->mean();
                const double availI = Transmission::PerHostAnophParams::get(i).entoAvailabilityFactor * popSize * adultAvail * p.mosqRelativeAvailabilityHuman;
                Transmission::Anopheles::Nhh nhh;

                nhh.availI = availI;
                nhh.pBI = p.mosqProbBiting;
                nhh.pCI = p.mosqProbFindingRestSite;
                nhh.pDI = p.mosqProbResting;
                nhh.relFecundity = p.hostFecundityFactor;
                nhh.expiry = Sim::now() + m_lifespan;

                // add the nhh to the active nhh instances
                anophModel->nhhInstances[m_intervName] = nhh;
            }
        }
    }
    virtual void printDetails( std::ostream& out )const{
        out << "AddNonHumanHosts(" << m_intervName << "):" << endl;
        out << "\tlifespan: " << m_lifespan;
    }
    
private:
    struct NhhParamsInterv {
        double mosqRelativeAvailabilityHuman;
        double mosqProbBiting;
        double mosqProbFindingRestSite;
        double mosqProbResting;
        double hostFecundityFactor;
    };
    string m_intervName;
    SimTime m_lifespan = Sim::never();
    std::map<string, NhhParamsInterv> m_nhhParams; 
};

// ———  ContinuousHumanDeployment  ———

/** Interface for continuous deployment of an intervention. */
class ContinuousHumanDeployment : protected HumanDeploymentBase {
public:
    /// Create, passing deployment age
    ContinuousHumanDeployment( SimTime begin, SimTime end,
                                 const ::scnXml::ContinuousDeployment& elt,
                                 shared_ptr<const HumanIntervention> intervention,
                                 ComponentId subPop, bool complement ) :
            HumanDeploymentBase( elt, intervention, subPop, complement ),
            m_begin( begin ), m_end( end ),
            m_deployAge( Sim::fromYearsN( elt.getTargetAgeYrs() ) ),
            m_minAvailability( 0 ),
            m_maxAvailability( 100 )
    {
        if( begin < Sim::startDate() || end < begin ){
            throw util::XmlScenarioError("continuous intervention must have startDate <= begin <= end");
        }
        if( m_deployAge <= Sim::zero() ){
            ostringstream msg;
            msg << "continuous intervention with target age "<<elt.getTargetAgeYrs();
            msg << " years corresponds to time step " << Sim::inSteps(m_deployAge);
            msg << "; must be at least 1.";
            throw util::XmlScenarioError( msg.str() );
        }
        if( m_deployAge > Sim::maxHumanAge() ){
            ostringstream msg;
            msg << "continuous intervention must have target age no greater than ";
            msg << Sim::inYears(Sim::maxHumanAge());
            throw util::XmlScenarioError( msg.str() );
        }
        if( elt.getMaxAvailability().present() )
            m_maxAvailability = elt.getMaxAvailability().get();

        if( elt.getMinAvailability().present() )
            m_minAvailability = elt.getMinAvailability().get();
            
        if( m_minAvailability < Sim::zero() || m_maxAvailability < m_minAvailability )
            throw util::XmlScenarioError("timed intervention must have 0 <= minAvailability <= maxAvailability <= 100");
    }
    
    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy( Host::Human& human ) const{
        SimTime age = human.age(Sim::now());
        if( m_deployAge > age ){
            // stop processing continuous deployments for this
            // human for now because remaining ones happen in the future
            return false;
        }else if( m_deployAge == age ){
            double availability = 0.0;
            for(size_t i = 0; i < Transmission::PerHostAnophParams::numSpecies(); ++i)
                availability += human.perHostTransmission.anophEntoAvailabilityRaw[i];

            if( availability >= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(m_minAvailability) && availability <= Transmission::PerHostAnophParams::getEntoAvailabilityPercentile(m_maxAvailability) )
            {
                auto now = Sim::intervDate();
                if( m_begin <= now && now < m_end &&
                    ( m_subPop == ComponentId::wholePop() ||
                        (human.isInSubPop( m_subPop ) != m_complement)
                    ) &&
                    human.rng.uniform01() < m_coverage )     // RNG call should be last test
                {
                    deployToHuman( human, mon::Deploy::CTS );
                }
            }
        }//else: for some reason, a deployment age was missed; ignore it
        return true;
    }
    
    inline void printDetails( std::ostream& out )const{
        out << "Human:" << endl;
        out << "\tbegin: " << m_begin << "\tend: ";
        if( m_end == Sim::future() ) out << "(none)";
        else out << m_end;
        out << "\tdeployAge: " << Sim::inYears(m_deployAge) << "y\tsubpop: ";
        if( m_subPop == ComponentId::wholePop() ) out << "whole";
        else out << m_subPop.id;
        out << "\tcomplement: " << m_complement << "\tcoverage: " << m_coverage << "\t";
        m_intervention->printDetails( out );
    }
    
protected:
    SimTime m_begin, m_end;    // first time step active and first time step no-longer active
    SimTime m_deployAge = Sim::never();
    SimTime m_minAvailability, m_maxAvailability;

    friend ByDeployTime;
};

struct ByDeployTime {
    bool operator() (const unique_ptr<TimedDeployment>& a, const unique_ptr<TimedDeployment>& b) const{
        return a->date < b->date;
    }
    
    bool operator()(const ContinuousHumanDeployment& a, const ContinuousHumanDeployment& b) const{
        return a.m_deployAge < b.m_deployAge;
    }
} byDeployTime;

} }
