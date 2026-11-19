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

#include "PkPd/Drug/LSTMDrugType.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Infection/CommonInfection.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "PkPd/Drug/LSTMDrugOneComp.h"
#include "PkPd/Drug/LSTMDrugThreeComp.h"
#include "PkPd/Drug/LSTMDrugConversion.h"
#include <schema/pharmacology.h>

#include <cmath>
#include <iomanip>

using namespace std;

namespace OM {
namespace PkPd {
using namespace OM::util;
using WithinHost::Genotypes;

// -----  Static variables and functions  -----

// The list of drugTypes drugs. Not checkpointed.
vector<LSTMDrugType> drugTypes;
map<string,size_t> drugTypeNames;
// List of all indices of drugs being used
vector<size_t> drugsInUse;


LSTMDrugPD::LSTMDrugPD( const scnXml::Phenotype& phenotype ){
    m_n = phenotype.getSlope ();
    IC50 = util::createSampler( phenotype.getIC50() );
    V = phenotype.getMax_killing_rate ();  
}

double LSTMDrugPD::calcFactor( double Kn, double negElimRate, double* conc, double duration ) const{
    const double C0 = *conc;
    const double C1 = C0 * exp(negElimRate * duration);
    
    // From Hastings & Winter 2011 paper
    // Note: these look a little different from original equations because Kn
    // is calculated when parameters are read from the scenario document instead of now.
    const double numerator = Kn + pow(C1, m_n);
    const double denominator = Kn + pow(C0, m_n);
    const double power = V / (-negElimRate * m_n);
    
    *conc = C1;    // conc is an in/out parameter
    return pow( numerator / denominator, power );       // unitless
}

double LSTMDrugPD::IC50_pow_slope(LocalRng& rng, size_t index, WithinHost::CommonInfection *inf) const{
    double Kn;  // gets sampled once per infection
    if (inf->Kn.count(index) != 0) {
        Kn = inf->Kn.at(index);
    } else {
        Kn = pow(IC50->sample(rng), m_n);
        inf->Kn[index] = Kn;
    }
    return Kn;
}


void LSTMDrugType::init (const scnXml::Pharmacology::DrugsType& drugData) {
    drugTypes.reserve(drugData.getDrug().size());
    for( const scnXml::PKPDDrug& drug : drugData.getDrug() ){
        const string& abbrev = drug.getAbbrev();
        // Check drug doesn't already exist
        if (drugTypeNames.find (abbrev) != drugTypeNames.end())
            throw TRACED_EXCEPTION_DEFAULT (string ("Drug added twice: ").append(abbrev));
        
        size_t i = drugTypes.size();
        drugTypes.emplace_back( i, drug );
        drugTypeNames[abbrev] = i;
    }
}
void LSTMDrugType::clear()
{
    drugTypes.clear();
    drugTypeNames.clear();
}

size_t LSTMDrugType::numDrugTypes(){
    return drugTypes.size();
}

// Add index to drugsInUse if not already present. Not fast but doesn't need to be.
void drugIsUsed(size_t index){
    for( size_t i : drugsInUse ){
        if( i == index ) return;        // already in list
    }
    drugsInUse.push_back(index);
}
size_t LSTMDrugType::findDrug(string abbreviation) {
    auto iter = drugTypeNames.find (abbreviation);
    if (iter == drugTypeNames.end())
        throw util::XmlScenarioError (string ("attempt to use drug without description: ").append(abbreviation));
    size_t index = iter->second;
    
    // We assume that drugs are used when and only when findDrug returns their
    // index or they are a metabolite of a drug returned here.
    drugIsUsed(index);
    if( drugTypes[index].m_conversionRate /*conversion model used*/ ){
        drugIsUsed(drugTypes[index].m_metabolite);
    }
    
    return index;
}

LSTMDrugType& LSTMDrugType::get(size_t index) { //static
    return drugTypes.at(index);
}

const vector< size_t >& LSTMDrugType::getDrugsInUse(){
    return drugsInUse;
}

unique_ptr<LSTMDrug> LSTMDrugType::createInstance(LocalRng& rng, size_t index) {
    LSTMDrugType& typeData = drugTypes[index];
    if( typeData.m_conversionRate ){
        LSTMDrugType& metaboliteData = drugTypes[typeData.m_metabolite];
        return unique_ptr<LSTMDrug>(new LSTMDrugConversion( typeData, metaboliteData, rng ));
    }else if( typeData.m_k12 ){
        // k21 is set when k12 is set; k13 and k31 may be set
        return unique_ptr<LSTMDrug>(new LSTMDrugThreeComp( typeData, rng ));
    }else{
        // none of k12/k21/k13/k31 should be set in this case
        return unique_ptr<LSTMDrug>(new LSTMDrugOneComp( typeData, rng ));
    }
}



// -----  Non-static LSTMDrugType functions  -----

LSTMDrugType::LSTMDrugType (size_t index, const scnXml::PKPDDrug& drugData) :
        m_index (index), m_metabolite(0),
        m_negligibleConcentration(numeric_limits<double>::quiet_NaN()),
        m_negMExp(numeric_limits<double>::quiet_NaN()),
        m_mwr(numeric_limits<double>::quiet_NaN()),
        m_ic50LogCorr(numeric_limits<double>::quiet_NaN()),
        m_ic50CorrFactor(numeric_limits<double>::quiet_NaN())
{
    // ———  PK parameters  ———
    const scnXml::PK& pk = drugData.getPK();
    m_negligibleConcentration = pk.getNegligible_concentration();
    if( pk.getHalf_life().present() ){
        if( pk.getK().present() || pk.getM_exponent().present() ){
            throw util::XmlScenarioError( "PK data must specify one of half_life or (k, m_exponent); it specifies both" );
        }
        m_eliminationRate = LognormalSampler::fromMeanCV( log(2.0) / pk.getHalf_life().get(), 0.0 );
        m_negMExp = 0.0; // no dependence on body mass
    }else{
        if( !(pk.getK().present() && pk.getM_exponent().present()) ){
            throw util::XmlScenarioError( "PK data must include either half_life or (k and m_exponent)" );
        }
        m_eliminationRate = util::createSampler( pk.getK().get() );
        m_negMExp = -pk.getM_exponent().get();
    }
    m_volDist = util::createSampler( pk.getVol_dist() );
    if( pk.getCompartment2().present() ){
        m_k12 = util::createSampler( pk.getCompartment2().get().getK12() );
        m_k21 = util::createSampler( pk.getCompartment2().get().getK21() );
        if( pk.getCompartment3().present() ){
            m_k13 = util::createSampler( pk.getCompartment3().get().getK13() );
            m_k31 = util::createSampler( pk.getCompartment3().get().getK31() );
        }else{
            // 2-compartment model: use 3-compartment code with these parameters set to zero
            // LognormalSampler supports 0,0 as a special case
            m_k13 = LognormalSampler::fromMeanCV(0.0, 0.0);
            m_k31 = LognormalSampler::fromMeanCV(0.0, 0.0);
        }
    }else if( pk.getCompartment3().present() ){
        throw util::XmlScenarioError( "PK model specifies parameters for "
                "compartment3 without compartment2" );
    }
    if( pk.getConversion().present() ){
        if( m_k12 ){
            throw util::XmlScenarioError( "PK conversion model is incompatible with 2/3-compartment model" );
        }
        const scnXml::Conversion& conv = pk.getConversion().get();
        try{
            m_metabolite = findDrug(conv.getMetabolite());
        }catch( util::XmlScenarioError& e ){
            throw util::XmlScenarioError( "PK: metabolite drug not found; metabolite must be defined *before* parent drug!" );
        }
        m_conversionRate = util::createSampler( conv.getRate() );
        m_mwr = conv.getMolRatio();
        double correlation = conv.getIC50_log_correlation();
        if( correlation < 0.0 || correlation > 1.0 ){
            throw util::XmlScenarioError( "PK conversion model: IC50 correlation must be between 0 and 1" );
        }
        m_ic50LogCorr = correlation;
        m_ic50CorrFactor = sqrt(1.0 - correlation * correlation);
    }
    if( pk.getK_a().present() ){
        if( !m_k12 && !m_conversionRate ){
            throw util::XmlScenarioError( "PK models only allow an "
                "absorption rate parameter (k_a) when compartment2 or "
                " conversion parameters are present" );
        }
        m_absorptionRate = util::createSampler(pk.getK_a().get());
    }else{
        if( m_k12 ){
            throw util::XmlScenarioError( "PK models require an absorption "
                "rate parameter (k_a) when compartment2 is present" );
        }
    }
    
    // ———  PD parameters  ———
    const scnXml::PD::PhenotypeSequence& pd = drugData.getPD().getPhenotype();
    assert( pd.size() > 0 );  // required by XSD
    
    set<string> lociPerPhenotype;
    string firstPhenotypeName;
    size_t nPhenotypes = pd.size();
    PD.reserve (nPhenotypes);
    // per phenotype (first index), per locus-of-restriction (second list), a list of alleles
    vector<vector<vector<uint32_t> > > phenotypeRestrictions;
    phenotypeRestrictions.reserve( nPhenotypes );
    for(size_t i = 0; i < nPhenotypes; ++i) {
        const scnXml::Phenotype& phenotype = pd[i];
        if( i == 0 ){
            if( phenotype.getName().present() )
                firstPhenotypeName = phenotype.getName().get();
            else
                firstPhenotypeName = to_string(i);
            for( const scnXml::Restriction& restriction : phenotype.getRestriction() ){
                lociPerPhenotype.insert( restriction.getOnLocus() );
            }
        }else{
            set<string> expectedLoci = lociPerPhenotype;
            for( const scnXml::Restriction& restriction : phenotype.getRestriction() ){
                size_t n = expectedLoci.erase( restriction.getOnLocus() );
                if( n == 0 ){
                    string thisPhenotypeName;
                    if( phenotype.getName().present() )
                        thisPhenotypeName = phenotype.getName().get();
                    else
                        thisPhenotypeName = to_string(i);
                    throw util::XmlScenarioError(
                        "pharmacology/drugs/drug/pd/phenotype/restriction/onLocus:"
                        "locus " + string(restriction.getOnLocus()) + " included in restriction of phenotype "
                        + string(thisPhenotypeName) + " but not " + string(firstPhenotypeName));
                }
            }
            if( !expectedLoci.empty() ){
                string thisPhenotypeName;
                if( phenotype.getName().present() )
                    thisPhenotypeName = phenotype.getName().get();
                else
                    thisPhenotypeName = to_string(i);
                throw util::XmlScenarioError(
                    "pharmacology/drugs/drug/pd/phenotype/restriction/onLocus:"
                    "locus " + string(*expectedLoci.begin()) + " included in restriction of phenotype "
                    + string(firstPhenotypeName) + " but not " + string(thisPhenotypeName));
            }
        }
        map<string,size_t> loci;
        vector<vector<uint32_t> > locAlleles;
        for( const scnXml::Restriction& restriction : phenotype.getRestriction() ){
            uint32_t allele = Genotypes::findAlleleCode( restriction.getOnLocus(), restriction.getToAllele() );
            if( allele == numeric_limits<uint32_t>::max() ){
                throw util::XmlScenarioError("phenotype has "
                    "restriction on locus " + string(restriction.getOnLocus())
                    + " allele " + string(restriction.getToAllele()) + " but this locus/allele "
                    "has not been defined in parasiteGenetics section");
            }
            auto iter = loci.find(restriction.getOnLocus());
            if( iter == loci.end() ){
                loci[restriction.getOnLocus()] = locAlleles.size();
                locAlleles.push_back( vector<uint32_t>(1,allele) );
            }else{
                locAlleles[iter->second].push_back( allele );
            }
        }
        phenotypeRestrictions.push_back( locAlleles );
        PD.push_back( LSTMDrugPD( pd[i] ) );
    }
    
    if( lociPerPhenotype.size() == 0 ){
        if( pd.size() > 1 ){
            throw util::XmlScenarioError( "pharmacology/drugs/drug/pd/phenotype:"
                " restrictions required when num. phenotypes > 1" );
        }else{
            // All genotypes map to phenotype 0
            m_genotypeMapping.assign( Genotypes::getGenotypes().size(), 0 );
        }
    }else{
        const vector<Genotypes::Genotype>& genotypes = Genotypes::getGenotypes();
        m_genotypeMapping.assign( genotypes.size(), 0 );
        for( size_t j = 0; j < genotypes.size(); ++j ){
            uint32_t phenotype = numeric_limits<uint32_t>::max();
            for( size_t i = 0; i < nPhenotypes; ++i ){
                // genotype matches this phenotype when, for every locus among
                // the phenotype restriction rules, one allele matches a
                // genotype allele
                bool match = true;
                for( const vector<uint32_t>& locAlleles : phenotypeRestrictions[i] ){
                    bool matchForLocus = false;
                    for( uint32_t restrictAllele : locAlleles ){
                        if( genotypes[j].alleles.count( restrictAllele ) > 0 ){
                            assert( !matchForLocus ); // shouldn't be two matches
                            matchForLocus = true;
                        }
                    }
                    if( !matchForLocus ){
                        match = false;
                        break;  // can skip rest of rules
                    }
                }
                if( match ){
                    if( phenotype == numeric_limits<uint32_t>::max() ){
                        phenotype = i;
                    }else{
                        // We could try to convert genotype's allele codes to
                        // locus/allele names, but requires several lookups
                        throw util::XmlScenarioError("phenotype restrictions "
                            "not restrictive enough: multiple phenotypes match "
                            "some genotype");
                    }
                }
            }
            if( phenotype == numeric_limits<uint32_t>::max() ){
                // We could try to convert genotype's allele codes to
                // locus/allele names, but requires several lookups
                throw util::XmlScenarioError("phenotype restrictions too"
                    "restrictive: no phenotype matching some genotype");
            }
            m_genotypeMapping[j] = phenotype;
        }
        
        if( util::CommandLine::option( util::CommandLine::PRINT_GENOTYPES ) ){
            // Debug output to see which genotypes correspond to phenotypes
            cout << endl;
            cout << "Phenotype mapping: " << endl << "----------" << endl;
            cout << "|" << setw(8) << "genotype" << "|" << setw(14) << "phenotype" << "|" << endl;
            cout << "|" << setw(8) << "--------" << "|" << setw(14) << "-------------" << "|" << endl;
            
            uint32_t genotype = 0;
            stringstream phenotypeName;
            for( auto phen = m_genotypeMapping.begin(); phen !=m_genotypeMapping.end(); ++phen ){
                phenotypeName.str("");
                if(pd[*phen].getName().present()){
                    phenotypeName << pd[*phen].getName().get();
                } else {
                    phenotypeName << "no: " << *phen;
                }
                cout << "|" << setw(8) << (genotype*100000) << "|" << setw(14) << phenotypeName.str() << "|" << endl;
                genotype += 1;
            }
        }
    }
}
LSTMDrugType::~LSTMDrugType () {
}

const LSTMDrugPD& LSTMDrugType::getPD( uint32_t genotype ) const {
    return PD[m_genotypeMapping[genotype]];
}

}
}
