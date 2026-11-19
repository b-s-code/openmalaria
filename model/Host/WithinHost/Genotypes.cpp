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

#include "Host/WithinHost/Genotypes.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/vectors.h"
#include "util/CommandLine.h"
#include "schema/scenario.h"

#include <iomanip>

namespace OM {
namespace WithinHost {

namespace GT /*for genotype impl details*/{
// ———  Model constants (after init)  ———
// keys are cumulative probabilities; last entry should equal 1; values are genotype codes
map<double,uint32_t> cumInitialFreqs;

// we give each allele of each loci a unique code
map<string, map<string, uint32_t> > alleleCodes;
uint32_t nextAlleleCode = 0;

// Represent a set of loci: all possible combinations of alleles
// This is actually just machinery to calculate the list of all genotypes
struct LocusSet{
    LocusSet( const scnXml::ParasiteLocus& eltL ){
        alleles.reserve( eltL.getAllele().size() );
        double cumP = 0.0;
        for( const scnXml::ParasiteAllele& eltA : eltL.getAllele() ){
            uint32_t alleleCode = (nextAlleleCode++);
            alleleCodes[eltL.getName()][eltA.getName()] = alleleCode;
            cumP += eltA.getInitialFrequency();
            alleles.push_back( Genotypes::Genotype( alleleCode, eltA.getInitialFrequency(), eltA.getFitness(), eltA.getHrp2_deletion() ) );
        }
        if( cumP < 0.999 || cumP > 1.001 )
            throw util::XmlScenarioError("expected sum of initial probabilities of alleles to be 1, but for the " 
                + to_string(alleles.size()) + " alleles under locus " + string(eltL.getName()) + " this is " + to_string(cumP));

        alleles[0].initFreq += 1.0 - cumP;     // account for any small errors by adjusting the first frequency
    }
    void include( LocusSet& that ){
        vector<Genotypes::Genotype> newAlleles;
        newAlleles.reserve( alleles.size() * that.alleles.size() );
        for( size_t i = 0; i < alleles.size(); i += 1 ){
            for( size_t j = 0; j < that.alleles.size(); j += 1 ){
                // note that we generate new codes each time; we just waste the old ones
                newAlleles.push_back( alleles[i].cross(that.alleles[j]) );
            }
        }
        alleles.swap(newAlleles);
    }
    
    vector<Genotypes::Genotype> alleles;
};

vector<Genotypes::Genotype> genotypes;
// ———  Model variables  ———
enum SampleMode{
    SAMPLE_FIRST,      // always choose first genotype (essentially the off switch)
    SAMPLE_INITIAL,    // sample from initial probabilities
    SAMPLE_TRACKING    // sample from tracked success at genotype level (no recombination)
};
// Mode to use now (until switched) and from the start of the intervention period.
SampleMode currentMode = SAMPLE_FIRST, intervMode = SAMPLE_FIRST;
}
size_t Genotypes::s_nGenotypes = 1;

Genotypes::Genotype Genotypes::Genotype::cross( const Genotype& that )const{
    Genotype result(*that.alleles.begin(),
                    initFreq*that.initFreq,
                    fitness*that.fitness,
                    hrp2Deficient || that.hrp2Deficient);
    result.alleles.insert( alleles.begin(), alleles.end() );
    result.alleles.insert( that.alleles.begin(), that.alleles.end() );
    return result;
}

void Genotypes::initSingle()
{
    // no specification implies there is a single genotype
    GT::genotypes.assign( 1, Genotypes::Genotype(
        0 /*allele code*/, 1.0/*frequency*/, 1.0/*fitness*/, false /*hrp2 deficiency*/) );
    s_nGenotypes = 1;
}

// utility function: does a vector contain an element?
template<class T>
bool contains(const vector<T>& vec, const T& x){
    for( auto it = vec.begin(), end = vec.end(); it != end; ++it ){
        if( *it == x ) return true;
    }
    return false;
}

void Genotypes::init( const scnXml::Scenario& scenario ){
    if( scenario.getParasiteGenetics().present() ){
        const scnXml::ParasiteGenetics& genetics =
            scenario.getParasiteGenetics().get();
        
        GT::currentMode = GT::SAMPLE_INITIAL;      // turn on sampling
        if( genetics.getSamplingMode() == "initial" ){
            GT::intervMode = GT::SAMPLE_INITIAL;
        }else if( genetics.getSamplingMode() == "tracking" ){
            GT::intervMode = GT::SAMPLE_TRACKING;
            if( !scenario.getEntomology().getVector().present() ){
                throw util::XmlScenarioError( "incompatibility; either use "
                    "entomology/vector (not nonVector) or set "
                    "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")" );
            }
            if( scenario.getEntomology().getMode() != "dynamic" ){
                throw util::XmlScenarioError( "incompatibility; either set "
                    "entomology/mode to \"dynamic\" (not \"forced\") or set "
                    "parasiteGenetics/samplingMode to \"initial\" (not \"tracking\")" );
            }
        }else{
            throw util::XmlScenarioError( "parasiteGenetics/samplingMode: expected \"initial\" or \"tracking\"" );
        }
        
        // Build the list of all allele combinations by iterating over loci (plural of locus):
        GT::LocusSet loci( genetics.getLocus()[0] );
        for( size_t i = 1; i < genetics.getLocus().size(); i += 1 ){
            GT::LocusSet newLocus(genetics.getLocus()[i]);
            loci.include( newLocus );
        }
        GT::genotypes.swap( loci.alleles );
        s_nGenotypes = GT::genotypes.size();
        
        double cumP = 0.0;
        for( size_t i = 0; i < GT::genotypes.size(); ++i ){
            cumP += GT::genotypes[i].initFreq;
            uint32_t genotypeId = static_cast<uint32_t>(i);
            assert( genotypeId == i );
            GT::cumInitialFreqs.insert( make_pair(cumP, genotypeId) );
        }
        
        // Test cum_p is approx. 1.0 in case the input tree is wrong.
        if (cumP < 0.999 || cumP > 1.001)
            throw util::XmlScenarioError ("decision tree (random node): expected probability sum to be 1.0 but found " + to_string(cumP));

        // last cum_p might be slightless less than 1 due to arithmetic errors; add a failsafe:
        GT::cumInitialFreqs[1.0] = GT::genotypes.size() - 1;
    }else{
        initSingle();
        GT::cumInitialFreqs[1.0] = GT::genotypes.size() - 1;
    }
    
    if( util::CommandLine::option( util::CommandLine::PRINT_GENOTYPES ) ){
        // reorganise GT::alleleCodes so that we can look up codes, not names
        vector<pair<string,string> > alleleCodes( GT::cumInitialFreqs.size() );
        for( auto i = GT::alleleCodes.begin(), iend = GT::alleleCodes.end(); i != iend; ++i ) {
            const string locus = i->first;
            for( auto j = i->second.begin(),
                jend = i->second.end(); j != jend; ++j )
            {
                uint32_t code = j->second;
                const string allele = j->first;
                assert( code < alleleCodes.size() );
                alleleCodes[code] = make_pair( locus, allele );
            }
        }
        
        // determine our columns
        map<string,uint32_t> longest;      // longest name in column; key is locus
        for( auto i = GT::genotypes.begin(),
            iend = GT::genotypes.end(); i != iend; ++i )
        {
            if( longest.size() == 0 ){
                for( auto j = i->alleles.begin();
                    j != i->alleles.end(); ++j ){
                    const string& locus = alleleCodes[*j].first;

                    longest[locus] = locus.length();  // locus name is included in column
                }
            }else assert( longest.size() == i->alleles.size() );
            
            for( auto j = i->alleles.begin();
                j != i->alleles.end(); ++j ){
                auto it = longest.find(alleleCodes[*j].first);
                assert( it != longest.end() );
                uint32_t lenAllele = alleleCodes[*j].second.length();
                if( lenAllele > it->second ) it->second = lenAllele;
            }
        }
        
        // find original loci order
        vector<string> loci;
        loci.reserve( longest.size() );
        for( size_t i = 0; loci.size() < longest.size(); ++i ){
            assert( i < alleleCodes.size() );
            const string& locus = alleleCodes[i].first;
            if( !contains(loci, locus) ){
                loci.push_back( locus );
            }
        }
        
//         string fm = "|%1$-12d|%|14t|%2$-12d|";
        cout << endl;
        stringstream fmt;
        fmt << "|%8d|";
        for( auto it = loci.begin(); it !=loci.end(); ++it )
            fmt << "%" << longest[*it] << "s|";
        fmt << "%9.3f|%7.3f|";
        
        // Table header:
        cout << "|" << std::setw(8) << "Genotype"; 
        for( auto it = loci.begin(); it !=loci.end(); ++it )
            cout << "|" << std::setw(longest[*it]) << *it;
        cout << "|" << std::setw(9) << "init freq" << "|" << std::setw(7) << "fitness" << "|" << endl;
        
        cout << "|" << std::setw(8) << "--------"; 
        for( auto it = loci.begin(); it !=loci.end(); ++it )
        {
            if(longest[*it] > 0)
                cout << "|" << "---------";
            else
                cout << "|";
        }

        cout << "|" << std::setw(9) << "---------" << "|" << std::setw(7) << "-------" << "|" << endl;

        cout << std::fixed << setprecision(3);
        for( size_t i = 0; i < GT::genotypes.size(); ++i ){
            const Genotype& genotype = GT::genotypes[i];
            map<string,string> locusAllele;    // to find allele for each locus
            for( auto a = genotype.alleles.begin(); a != genotype.alleles.end(); ++a ){
                assert( *a < alleleCodes.size() );
                locusAllele[alleleCodes[*a].first] = alleleCodes[*a].second;
            }

            cout << "|" << std::setw(8) << (i*100000);
            for( auto it = loci.begin(); it !=loci.end(); ++it ){
                auto la = locusAllele.find( *it );
                assert( la != locusAllele.end() );
                cout << "|" << std::setw(longest[*it]) << la->second << "|";

            }
            cout << std::setw(9) << genotype.initFreq << "|" << std::setw(7) << genotype.fitness << "|" << endl;
        }
    }
}

void Genotypes::preMainSimInit(){
    GT::currentMode = GT::intervMode;
}

uint32_t Genotypes::findAlleleCode(const string& locus, const string& allele){
    auto it = GT::alleleCodes.find(locus);
    if( it == GT::alleleCodes.end() ) return numeric_limits<uint32_t>::max();
    auto it2 = it->second.find( allele );
    if( it2 == it->second.end() ) return numeric_limits<uint32_t>::max();
    return it2->second;
}

const vector< Genotypes::Genotype >& Genotypes::getGenotypes(){
    return GT::genotypes;
}

uint32_t Genotypes::sampleGenotype( LocalRng& rng, vector<double>& genotypeWeights ){
    if( GT::currentMode == GT::SAMPLE_FIRST ){
        return 0;       // always the first genotype code
    }else if( GT::currentMode == GT::SAMPLE_INITIAL
            || genotypeWeights.size() == 0 )
    {
        double sample = rng.uniform01();
        auto it = GT::cumInitialFreqs.upper_bound( sample );
        assert( it != GT::cumInitialFreqs.end() );
        return it->second;
    }else{
        assert( GT::currentMode == GT::SAMPLE_TRACKING );
        assert( genotypeWeights.size() == s_nGenotypes );
        double weightSum = util::vectors::sum( genotypeWeights );
        assert( weightSum >= 0.0 && weightSum < 1e5 );        // possible loss of precision or other error
        double sample = rng.uniform01() * weightSum;
        double cum = 0.0;
        for( size_t g = 0; g < s_nGenotypes; ++g ){
            cum += genotypeWeights[g];
            if( sample < cum ) return g;
        }
        return 0;       // just to be safe (could happen if weight_sum == 0.0)
    }
}

double Genotypes::initialFreq( size_t genotype ){
    if( GT::genotypes.size() == 0 ){
        assert( genotype == 0 );
        assert( s_nGenotypes == 1 );
        return 1.0;     // only 1 genotype, thus initial frequency is 100%
    }else{
        return GT::genotypes[genotype].initFreq;
    }
}


// ———  checkpointing  ———

void Genotypes::staticCheckpoint( ostream& stream ){
    int t = GT::currentMode;
    t & stream;
}
void Genotypes::staticCheckpoint( istream& stream ){
    int t;
    t & stream;
    GT::currentMode = static_cast<GT::SampleMode>(t);
    assert( GT::currentMode == GT::SAMPLE_FIRST ||
        GT::currentMode == GT::SAMPLE_INITIAL ||
        GT::currentMode == GT::SAMPLE_TRACKING );
}

}
}
