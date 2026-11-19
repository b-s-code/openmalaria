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

#ifndef Hmod_NonVectorModel
#define Hmod_NonVectorModel

#include "TransmissionModel.h"
#include "util/vectors.h"
#include "Host/Human.h"
#include "mon/Monitoring.h"

namespace scnXml
{
class NonVector;
}

namespace OM
{
namespace Transmission
{

//! Base transmission model, as used in Phase A
class NonVectorModel : public TransmissionModel
{
public:
    NonVectorModel(vector<double> initEIR, int interventionMode, const scnXml::Entomology &entoData, const scnXml::NonVector &nonVectorData, int eipDuration)
        : TransmissionModel(std::move(initEIR), interventionMode, 1 /*this model doesn't support multiple genotypes*/)
        , m_nSpore(Sim::fromDays(eipDuration))
    {
        laggedKappa.resize(Sim::inSteps(m_nSpore) + 1, 0.0);

        vector<int> nDays(Sim::stepsPerYear(), 0);
        // The minimum EIR allowed in the array. The product of the average EIR and a constant.
        double minEIR = min_EIR_mult * averageEIR(nonVectorData);

        const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily();
        if (daily.size() < static_cast<size_t>(Sim::oneYear()))
            throw util::XmlScenarioError("insufficient EIRDaily data for a year");

        for (SimTime mpcday = Sim::zero(), endDay = Sim::fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + Sim::oneDay())
        {
            double eiRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

            // Index 0 of initialisationEIR refers to the EIR affecting the
            // first day(s) of the year. Correspondingly, the first 1 or 5 values
            // of EIRDaily affect this (1- or 5-day) time-step.
            size_t i = modNn(Sim::inSteps(mpcday), Sim::stepsPerYear());

            nDays[i] += 1;
            initialisationEIR[i] += eiRdaily;
        }

        // Calculate total annual EIR
        // divide by number of records assigned to each interval (usually one per day)
        for (size_t indTS = 0; indTS < Sim::stepsPerYear(); indTS += 1)
        {
            initialisationEIR[indTS] *= Sim::oneTS() / static_cast<double>(nDays[indTS]);
            annualEIR += initialisationEIR[indTS];
        }

        m_initialKappa.assign(Sim::inSteps(Sim::fromYearsI(m_nYearsWarmupData)), 0.0);
    }

    virtual ~NonVectorModel() {}

    virtual void init2(const vector<Host::Human> &population) { simulationMode = forcedEIR; } // no set-up needed; just indicate we're ready to roll

    virtual void scaleEIR(double factor)
    {
        util::vectors::scale(initialisationEIR, factor);
        annualEIR = util::vectors::sum(initialisationEIR);
    }

    virtual SimTime initIterate()
    {
        simulationMode = interventionMode;
        if (simulationMode != dynamicEIR) { return Sim::zero(); }

        // initialKappa is used in calculateEIR
        assert(m_initialKappa.size() >= Sim::stepsPerYear());
        assert(modNn(m_initialKappa.size(), Sim::stepsPerYear()) == 0);
        for (size_t i = Sim::stepsPerYear(); i < m_initialKappa.size(); ++i)
        {
            m_initialKappa[modNn(i, Sim::stepsPerYear())] += m_initialKappa[i];
        }
        double factor = static_cast<double>(Sim::stepsPerYear()) / static_cast<double>(m_initialKappa.size());
        m_initialKappa.resize(Sim::stepsPerYear());
        for (size_t i = 0; i < m_initialKappa.size(); ++i)
        {
            m_initialKappa[i] *= factor;
            // error check:
            if (!(m_initialKappa[i] > 0.0)) // if not positive
                throw TRACED_EXCEPTION("initialKappa is invalid", util::Error::InitialKappa);
        }

        return Sim::zero(); // nothing to do
    }

    virtual void changeEIRIntervention(const scnXml::NonVector &nonVectorData)
    {
        // Note: requires Sim::intervTime() >= Sim::zero(), but this can only be
        // called in intervention period anyway.
        simulationMode = transientEIRknown;

        if (m_nSpore != Sim::fromDays(nonVectorData.getEipDuration()))
            throw util::XmlScenarioError("change-of-EIR intervention cannot change EIP duration");

        const scnXml::NonVector::EIRDailySequence &daily = nonVectorData.getEIRDaily();
        vector<int> nDays(Sim::inSteps(Sim::fromDays(daily.size() - 1)) + 1, 0);
        m_interventionEir.assign(nDays.size(), 0.0);
        size_t requiredDays = static_cast<size_t>((Sim::endDate() - Sim::startDate()) + 1);
        if (daily.size() < requiredDays)
        {
            cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size() << "\nRequired: " << requiredDays << endl;
            throw util::XmlScenarioError("Insufficient intervention phase EIR values provided");
        }
        // The minimum EIR allowed in the array. The product of the average EIR and a constant.
        double minEIR = min_EIR_mult * averageEIR(nonVectorData);
        for (SimTime mpcday = Sim::zero(), endDay = Sim::fromDays(daily.size()); mpcday < endDay; mpcday = mpcday + Sim::oneDay())
        {
            double eiRdaily = std::max(static_cast<double>(daily[mpcday]), minEIR);

            // istep is the time period to which the day is assigned.
            size_t istep = Sim::inSteps(mpcday);
            nDays[istep]++;
            m_interventionEir[istep] += eiRdaily;
        }
        // divide by number of records assigned to each interval (usually one per day)
        for (size_t i = 0; i < m_interventionEir.size(); ++i)
        {
            m_interventionEir[i] *= Sim::oneTS() / static_cast<double>(nDays[i]);
        }

        // I've no idea what this should be, so until someone asks it can be NaN.
        // It was -9.99 and later 0.0. It could of course be recalculated from interventionEIR.
        annualEIR = numeric_limits<double>::quiet_NaN();
    }

    virtual void uninfectVectors()
    {
        if (simulationMode != dynamicEIR) cerr << "Warning: uninfectVectors is not efficacious with forced EIR" << endl;
        // reset history of human infectivity, which scales dynamic EIR:
        laggedKappa.assign(laggedKappa.size(), 0.0);
    }

    virtual double updateKappa(const vector<Host::Human> &population)
    {
        double currentKappa = TransmissionModel::updateKappa(population);
        if (simulationMode == forcedEIR) { m_initialKappa[Sim::moduloSteps(Sim::ts1(), m_initialKappa.size())] = currentKappa; }
        return currentKappa;
    }

    virtual void calculateEIR(Host::Human &human, double ageYears, vector<double> &EIR_i, vector<double> &EIR_l) const
    {
        EIR_i.resize(1);
        EIR_l.resize(1); // no support for per-genotype tracking in this model (possible, but we're lazy)
        // where the full model, with estimates of human mosquito transmission is in use, use this:
        if (simulationMode == forcedEIR) { EIR_l[0] = initialisationEIR[Sim::moduloYearSteps(Sim::ts0())]; }
        else if (simulationMode == transientEIRknown)
        {
            // where the EIR for the intervention phase is known, obtain this from
            // the interventionEIR array
            EIR_l[0] = m_interventionEir[Sim::inSteps(Sim::intervTime())];
        }
        else if (simulationMode == dynamicEIR)
        {
            EIR_l[0] = initialisationEIR[Sim::moduloYearSteps(Sim::ts0())];
            if (Sim::intervTime() >= Sim::zero())
            {
                // we modulate the initialization based on the human infectiousness time steps ago in the
                // simulation relative to infectiousness at the same time-of-year, pre-intervention.
                // nspore gives the sporozoite development delay.
                size_t t = Sim::inSteps(Sim::ts1() - m_nSpore);
                EIR_l[0] *= laggedKappa[modNn(t, laggedKappa.size())] / m_initialKappa[modNn(t, Sim::stepsPerYear())];
            }
        }
        else
        {
            throw util::XmlScenarioError("Invalid simulation mode");
        }
#ifndef NDEBUG
        if (!(std::isfinite)(EIR_l[0]))
        {
            size_t t = Sim::inSteps(Sim::ts1() - m_nSpore);
            ostringstream msg;
            msg << "Error: non-vect eir is: " << EIR_l[0] << "\nlaggedKappa:\t" << laggedKappa[modNn(t, laggedKappa.size())]
                << "\ninitialKappa:\t" << m_initialKappa[modNn(t, Sim::stepsPerYear())] << endl;
            throw TRACED_EXCEPTION(msg.str(), util::Error::InitialKappa);
        }
#endif
        EIR_l[0] *= human.perHostTransmission.relativeAvailabilityHetAge(ageYears);

        mon::recordStat(mon::measure("innoculationsPerAgeGroup"), human, EIR_l[0]);
        mon::recordStat(mon::measure("innoculationsPerVector"), human, EIR_l[0]);
    }

private:
    virtual void checkpoint(istream &stream)
    {
        TransmissionModel::checkpoint(stream);
        m_nSpore &stream;
        m_interventionEir &stream;
        m_initialKappa &stream;
    }

    virtual void checkpoint(ostream &stream)
    {
        TransmissionModel::checkpoint(stream);
        m_nSpore &stream;
        m_interventionEir &stream;
        m_initialKappa &stream;
    }

    //! multiplier used to calculate a positive EIR value where the measured value is zero
    /*
      0.01 was old pv(30) Now a constant. min_EIR_mult multiplies the average EIR to obtain a value used for the EIR during periods when
      it is too low to be measureable. The value of 0.01 was old pv(30) Now a constant. 0.01 was old pv(30) Now a constant.
    */
    const double min_EIR_mult = 0.01;

    /** @brief Variables set by constructor.
     *
     * There shouldn't be any need to checkpoint these, at least before
     * interventions take effect. */
    //@{
    //! Variance of Infection Rate according to fielddata
    const double m_totalInfectionrateVariance = 1.0;

    const int m_nYearsWarmupData = 5;

    //! The duration of sporogony in time steps
    // doesn't need checkpointing
    SimTime m_nSpore = Sim::never();
    //@}

    /** EIR per time interval during the intervention period. Value at index
     * Sim::intervTime().inSteps() used each time-step.
     *
     * Units: inoculations per adult per time step */
    vector<double> m_interventionEir;

    /** When simulationMode == dynamicEIR, this is the annual cycle of kappa
     * from the warmup phase and has length 1 year (in time steps).
     *
     * When simulationMode == equilibriumMode, this may be multiple years long and
     * is used to collect values of kappa (human infectiousness).
     *
     * In either case, Sim::ts0().moduloSteps(initialKappa.size()) is the index
     * for the current infectiousness during updates. */
    vector<double> m_initialKappa;
};
} // namespace Transmission
} // namespace OM
#endif
