//
// VoidPtrStore.h
//
// This file contains the BTL_PtrSet BTL_IIand BTL_CII classes. These
// classes were developed specifically for form part of the implementation
// of the BTL_Graph and related classes. Outside of this context, the classes
// within this file should be used WITH CAUTION.
//
// This code is part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997,1998 Birkbeck College, Malet Street, London WC1E 7HX, 
// U.K. (classlib@mail.cryst.bbk.ac.uk)
//
// This library is free software; you can
// redistribute it and/or modify it under the terms of
// the GNU Library General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////

#if !defined(BTL_VOIDPTRSTORE_H)
#define BTL_VOIDPTRSTORE_H

#include <set>
#include "BTL.h"

using namespace std;


_BTL_BEGIN_NAMESPACE

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

// This class replaces the STL set container when a pointer type is to be
// stored. It implemented in such a way that only one set type is instantiated.
// This reduces code bloat and, more importantly for the graph classes, it
// reduces template function signatures that had been reaching compiler limits.
//
// N.B. this case should only be used for pointer types that can be safely
// converted to and from void*. The iterators once dereferenced need to be caste
// (T*) or alteratively converted to a BTL IndirectIterator (BTL_II or BTL_CII).
// These iterators do the type casting for you.

                /**#: [Hidden] */

template <class T>
class BTL_PtrSet
{
public:

    typedef set<void*, less<void*> > 	    PointerStore;
    typedef PointerStore::iterator  	    iterator;
    typedef PointerStore::const_iterator    const_iterator;
    typedef PointerStore::size_type 	    size_type;
//    typedef T*      	    	    	    value_type;

private:
    
    PointerStore rep;
    
public:

    iterator
    insert(iterator position, T* t) { return rep.insert(position,t); } 

    iterator 
    find(T* t) const { return rep.find(t); }
    
    size_type
    size() const { return rep.size(); }
    
    void 
    erase(iterator i) { rep.erase(i); }
    
    iterator 
    begin() { return rep.begin(); }

    const_iterator 
    begin() const { return rep.begin(); }
    
    iterator 
    end() { return rep.end(); }

    const_iterator 
    end() const { return rep.end(); }

};

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


//     This is an iterator class to be used for interating through containers of
//     pointers. The derefencing operator deferences twice to give the object
//     referenced by the pointer at the current position within the container.
//     This is a bidirectional iterator. N.B. this iterator works equally well
//     with a container of void*. Containers of this sort can be used to reduce
//     code bloat etc.."] [Summary = "generic iterator for use with containers
//     of pointers"]  [Authors = "W.R.Pitt"] [Files = "<A
//     HREF=./btl/II.h>II.h</A>"] [Friends="None"] [Dependencies="None"]
//

                /**#: [Hidden] */

template<class T>
class BTL_CII;

                /**#: [Hidden] */

template<class T>
class BTL_II
{
    typedef BTL_PtrSet<T>   	    	    PointerStore;
    typedef BTL_TYPENAME PointerStore::iterator iterator;
    typedef T*  		    	    Pointer;

    friend class BTL_CII<T>;
    
private:
    
    iterator current;
    
public:

	    	/**#: [Description="Construct new floating II"] */

    BTL_II() {}

	    	/**#: [Description="Construct new II from a
	    	       PointerStore iterator."] */

    BTL_II(const iterator& i) : current(i) {}
    
	    	/**#: [Description="Copy operation"] */

    BTL_II&
    operator=(const BTL_II& i) { current=i.current;return *this; }

	    	/**#: [Description="Assign II from a
	    	       PointerStore iterator."] */

    BTL_II&
    operator=(const iterator& i) { current = i; return *this; }

	    	/**#: [Description="Prefix ++ operator."] */

    BTL_II&
    operator++() { current++; return *this;}
    
	    	/**#: [Description="Postfix ++ operator."] */

    BTL_II
    operator++(int) 
    { 
    	BTL_II tmp = *this;  
    	current++;
    	return tmp;
    }
    
	    	/**#: [Description="Prefix -- operator."] */

    BTL_II&
    operator--() { current--; return *this;}

	    	/**#: [Description="Postfix -- operator."] */

    BTL_II
    operator--(int) 
    {
    	BTL_II tmp = *this;
    	current--; 
    	return tmp;
    }

	    	/**#: [Description="Equality operator."] */

    bool		
    operator==(BTL_II i) const { return current==i.current; }

	    	/**#: [Description="Deference twice, via pointer stored in the
	    	       container. "] */

    T&
    operator*() const 
    {   
    	// A type caste is used here because current references a void*
    	return * ( (T*) (*current));
    }
    
	    	/**#: [Description="Dereference, but only as far as pointer
	    	       stored in the container."] */

    	// A type caste is used here because current references a void*

    Pointer
    GetPointer() const { return (T*) (*current); }
};

//------------------------------------------------------------------------------
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//------------------------------------------------------------------------------


//     [Description ="The const version of BTL_II."]
//     [Summary = "const version of BTL_II"] 
//     [Authors = "W.R.Pitt"]
//     [Files = "<A HREF=./btl/II.h>II.h</A>"]
//     [Friends="None"]
//     [Dependencies="None"]

                /**#: [Hidden] */

template<class T>
class BTL_CII
{

    typedef BTL_PtrSet<T>   	    	    	    PointerStore;
    typedef BTL_TYPENAME PointerStore::iterator  	    iterator;
    typedef BTL_TYPENAME PointerStore::const_iterator   const_iterator;
    typedef T*  			    	    Pointer;
    typedef BTL_II<T>	    	    	    	    BTL_II;

private:
    
    const_iterator current;
    
public:

	    	/**#: [Description="Construct new floating 
	    	       CII"] */

    BTL_CII() {}

	    	/**#: [Description="Construct new CII from a
	    	       PointerStore iterator."] */

    BTL_CII(const iterator& i) : current(i) {}

// 	    	/**#: [Description="Construct new CII from a
// 	    	       PointerStore const_iterator."] */
// 
//    BTL_CII(const const_iterator& i) : current(i) {}
    
	    	/**#: [Description="Copy operation"] */

    BTL_CII&
    operator=(const BTL_CII& i) 
    	    	    	    	    { current=i.current;return *this; }

	    	/**#: [Description="Copy from BTL_II"] */

    BTL_CII&
    operator=(const BTL_II& i) { current=i.current;return *this; }

	    	/**#: [Description="Assign CII from a
	    	       PointerStore iterator."] */

    BTL_CII&
    operator=(const iterator& i) { current = i; return *this; }

	    	/**#: [Description="Prefix ++ operator."] */

    BTL_CII&
    operator++() { current++; return *this;}
    
	    	/**#: [Description="Postfix ++ operator."] */

    BTL_CII
    operator++(int) 
    { 
    	BTL_CII tmp = *this;  
    	current++;
    	return tmp;
    }
    
	    	/**#: [Description="Prefix -- operator."] */

    BTL_CII&
    operator--() { current--; return *this;}

	    	/**#: [Description="Postfix -- operator."] */

    BTL_CII
    operator--(int) 
    {
    	BTL_CII tmp = *this;
    	current--; 
    	return tmp;
    }

	    	/**#: [Description="Equality operator."] */

    bool		
    operator==(BTL_CII i) const { return current==i.current; }

	    	/**#: [Description="Deference twice, via pointer stored in the
	    	       container. "] */

    T
    operator*() const 
    { 
    	// A type caste is used here because current references a void*
    	return *((T*) (*current));
    }
    
	    	/**#: [Description="Dereference, but only as far as pointer
	    	       stored in the container."] */

    	// A type caste is used here because current references a void*
    Pointer
    GetPointer() const { return (T*) (*current); }
};

_BTL_END_NAMESPACE

#endif
