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

#include "Global.h"
#include "util/errors.h"
#include "interventions/Interfaces.h"

#include <limits>
#include <sstream>
#include <assert.h>
using namespace std;

namespace OM { namespace util { namespace checkpoint {
    
    // Header test constants
    const unsigned int h_BOM = 0x50434D4F;      // "OMCP" in little-endian: OpenMalaria CheckPoint
    const bool HB = true;
    const unsigned char HC = 0xA5;     // binary: 10100101; don't care if char is signed as long as first bit is read and written correctly
    const double HN0 = -0.0;
    const double HNan = numeric_limits<double>::quiet_NaN();
    
    
    void validateListSize (long length, long max) {
        if (length < 0 || length > max) {
            ostringstream s;
            s << "List length out of range: " << length;
            throw CheckpointError(s.str());
        }
    }
    
    /** @brief Binary checkpointing
     *
     * Writes value in direct binary (not suitible for pointers). Currently
     * machine dependent (don't transfer a checkpoint to another machine).
     */
    //@{
    template<class T>
    void binaryWrite (T x, ostream& stream) {
        stream.write (reinterpret_cast<char*>(&x), sizeof(x));
    }
    template<class T>
    void binaryRead (T& x, istream& stream) {
        stream.read (reinterpret_cast<char*>(&x), sizeof(x));
        if (!stream || stream.gcount() != sizeof(x))
            throw CheckpointError ("stream read error binary");
    }
    //@}
    
    void staticChecks () {
        static_assert (sizeof(char) == 1, "sizeof(char) == 1");
        static_assert (sizeof(short int) == 2, "sizeof(short int) == 2");
        static_assert (sizeof(int) == 4, "sizeof(int) == 4");
        static_assert (sizeof(float) == 4, "sizeof(float) == 4");
        static_assert (sizeof(double) == 8, "sizeof(double) == 8");
    }
    
    void header (ostream& stream) {
        staticChecks();
        
        binaryWrite (h_BOM, stream);
        binaryWrite (HB, stream);
        binaryWrite (HC, stream);
        binaryWrite (HN0, stream);
        binaryWrite (HNan, stream);
    }
    void header (istream& stream) {
        staticChecks ();
        unsigned int BOM;
        bool b;
        unsigned char c;
        double n0;
        double nan;
        
        binaryRead (BOM, stream);
        binaryRead (b, stream);
        binaryRead (c, stream);
        binaryRead (n0, stream);
        binaryRead (nan, stream);
        
        // Check. Use binary check for doubles since it's not the same as numeric ==
        if (BOM != h_BOM ||
            b != HB ||
            c != HC ||
            memcmp (&n0, &HN0, sizeof(double)) ||
            memcmp (&nan, &HNan, sizeof(double))
            )
            throw CheckpointError ("invalid header");
    }
    
    ///@brief Operator& for simple data-types
    //@{
    void operator& (bool x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (bool& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (signed char x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (signed char& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (short x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (short& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (int x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (int& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (long x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (long& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (long long x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (long long& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (unsigned char x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (unsigned char& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (unsigned short x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (unsigned short& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (unsigned int x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (unsigned int& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (unsigned long x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (unsigned long& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (unsigned long long x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (unsigned long long& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (float x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (float& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (double x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (double& x, istream& stream) {
        binaryRead (x, stream);
    }
    
    void operator& (long double x, ostream& stream) {
        binaryWrite (x, stream);
    }
    void operator& (long double& x, istream& stream) {
        binaryRead (x, stream);
    }
    //@}
    
    // string
    void operator& (string x, ostream& stream) {
        x.length() & stream;
        stream.write (x.c_str(), x.length());
    }
    void operator& (string& x, istream& stream) {
        size_t len;
        len & stream;
        validateListSize (len);
        x.resize (len);
        stream.read (&x[0], x.length());
        if (!stream || stream.gcount() != streamsize(len))
            throw CheckpointError ("stream read error string");
    }
    
    void operator& (const set<interventions::ComponentId>&x, ostream& stream){
        x.size() & stream;
        for( auto it = x.begin(); it != x.end(); ++it )
            *it & stream;
    }
    void operator& (set<interventions::ComponentId>& x, istream& stream){
        size_t len;
        len & stream;
        validateListSize(len);
        for( size_t i = 0; i < len; ++i ){
            interventions::ComponentId id( stream );
            x.insert( id );
        }
    }
    
    // map<S,T> — template doesn't work on gcc
    void operator& (const map<string,double>& x, ostream& stream) {
        x.size() & stream;
        for(auto pos = x.begin (); pos != x.end() ; ++pos) {
            pos->first & stream;
            pos->second & stream;
        }
    }
    void operator& (map<string,double>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.clear ();
        auto pos = x.begin ();
        for(size_t i = 0; i < l; ++i) {
            string s;
            double t;
            s & stream;
            t & stream;
            pos = x.insert (pos, make_pair (s,t));
        }
    }

    void operator& (const map<double,double>& x, ostream& stream) {
        x.size() & stream;
        for(auto pos = x.begin (); pos != x.end() ; ++pos) {
            pos->first & stream;
            pos->second & stream;
        }
    }
    void operator& (map<double,double>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.clear ();
        auto pos = x.begin ();
        for(size_t i = 0; i < l; ++i) {
            double s, t;
            s & stream;
            t & stream;
            pos = x.insert (pos, make_pair (s,t));
        }
    }

    void operator& (const map<interventions::ComponentId,SimTime>& x, ostream& stream) {
        x.size() & stream;
        for(auto pos = x.begin (); pos != x.end() ; ++pos) {
            pos->first & stream;
            SimTime t = pos->second;    // use a copy: avoid const restriction
            t & stream;
        }
    }
    void operator& (map<interventions::ComponentId,SimTime>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.clear ();
        auto pos = x.begin ();
        for(size_t i = 0; i < l; ++i) {
            interventions::ComponentId s( stream );
            SimTime t = Sim::never();
            t & stream;
            pos = x.insert (pos, make_pair (s,t));
        }
    }

    void operator& (const multimap<double,double>& x, ostream& stream) {
        x.size() & stream;
        for(auto pos = x.begin (); pos != x.end() ; ++pos) {
            pos->first & stream;
            pos->second & stream;
        }
    }
    void operator& (multimap<double,double>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.clear ();
        auto pos = x.begin ();
        for(size_t i = 0; i < l; ++i) {
            double s, t;
            s & stream;
            t & stream;
            pos = x.insert (pos, make_pair (s,t));
        }
        assert( x.size() == l );
    }
    
} } }
