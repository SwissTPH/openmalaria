/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#ifndef OM_util_checkpoint_containers
#define OM_util_checkpoint_containers

#ifndef Hmod_Global
#error "Please include Global.h first."
// otherwise "using ..." declaration in Global.h won't work
#endif

/** Provides some extra functions. See checkpoint.h. */
namespace OM {
namespace util {
namespace checkpoint {

    ///@brief Operator& for stl containers
    //@{
    template<class U, class V, class S>
    void operator& (pair<U,V> x, S& stream) {
        x.first & stream;
        x.second & stream;
    }
    
    template<class T>
    void operator& (vector<T> x, ostream& stream) {
        x.size() & stream;
        foreach (T& y, x) {
            y & stream;
        }
    }
    template<class T>
    void operator& (vector<T>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l);
        foreach (T& y, x) {
            y & stream;
        }
    }
    /// Version of above taking an element to initialize each element from.
    template<class T>
    void checkpoint (vector<T>& x, istream& stream, T templateInstance) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l, templateInstance);
        foreach (T& y, x) {
            y & stream;
        }
    }
    
    template<class T>
    void operator& (list<T> x, ostream& stream) {
        x.size() & stream;
        foreach (T& y, x) {
            y & stream;
        }
    }
    template<class T>
    void operator& (list<T>& x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.resize (l);
        foreach (T& y, x) {
            y & stream;
        }
    }
    
    /* The following templates compile on gcc but don't appear to work when S is a string and T a double.
    // Template templates:
    template<class S, class T>
    void operator& (map<S,T> x, ostream& stream) {
        x.size() & stream;
        for( typename map<S,T>::const_iterator it = x.begin(); it != x.end(); ++it ){
            it->first & stream;
            it->second & stream;
        }
        cerr<<"operator&(map<S,T> x, ostream&) where S="<<typeid(S).name()<<", T="<<typeid(T).name()<<", x.size()="<<x.size()<<endl;
    }
    template<class S, class T>
    void operator& (map<S,T> x, istream& stream) {
        size_t l;
        l & stream;
        validateListSize (l);
        x.clear ();
        typename map<S,T>::iterator pos = x.begin ();
        for(size_t i = 0; i < l; ++i) {
            S s;
            T t;
            s & stream;
            t & stream;
            pos = x.insert (pos, make_pair (s,t));
        }
    }
    */
    
    void operator& (const set<interventions::ComponentId>&x, ostream& stream);
    void operator& (set<interventions::ComponentId>& x, istream& stream);
    
    void operator& (const map<string,double>& x, ostream& stream);
    void operator& (map<string, double >& x, istream& stream);
    
    void operator& (const map<double,double>& x, ostream& stream);
    void operator& (map<double, double>& x, istream& stream);
    
    void operator& (const map<interventions::ComponentId,SimTime>& x, ostream& stream);
    void operator& (map<interventions::ComponentId,SimTime>& x, istream& stream);
    
    void operator& (const multimap<double,double>& x, ostream& stream);
    void operator& (multimap<double, double>& x, istream& stream);
    //@}
    
} } }   // end of namespaces
#endif
