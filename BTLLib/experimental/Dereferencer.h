// Dereferencer.h
//
// This file contains the BTL_InD and BTL_ItD classes. These
// classes were developed specifically for form part of the implementation
// of the BTL_Graph and related classes. Outside of this context, the classes
// within this file should be used WITH CAUTION.
//
// This code is part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997, 1998 Birkbeck College, Malet Street, London WC1E 7HX,
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


#if !defined(BTL_DEREFERENCER_H)
#define BTL_DEREFERENCER_H

_BTL_BEGIN_NAMESPACE


//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


//  [Description ="This class, together with BTL_ItD provide
//   a uniform interface STL container classes. These classes are necessary if
//   persistant (i.e. for the life time of the container object) references to
//   individual items with the container are required. This class can be used
//   when individual items within a container are to be referenced using an
//   index number. This is required for vector and deque as the iterators for
//   these containers can become invalid after insertions and erasures. All
//   other STL containers should be used with the BTL_ItD."]
//  [Summary = "Allows stable dereferencing to the vector and deque containers
//   via index numbers."]
//  [Usage = "This is template class the container type to be used must be
//   provided as the template parameter. It will not work for containers that
//   do not provide an operator[int] type function."]
//  [Authors = "W.R.Pitt"]
//  [Files = "<A HREF=./btl/Derefer.h>Derefer.h</A>"]
//  [Friends="None"]
//  [Dependencies="None"]


                /**#: [Hidden] */

template <class Container>
class BTL_InD
{
public:

    typedef Container container_type;	    	    	    // The container
    typedef BTL_TYPENAME Container::value_type value_type;	    	    // Its contents
    typedef BTL_TYPENAME Container::size_type reference_type;    	    // The index type
    typedef BTL_TYPENAME Container::iterator iterator;   	    	    // Its iterator
    typedef BTL_TYPENAME Container::const_iterator const_iterator;	    // Const iterator
    typedef BTL_InD<Container> Derefer;  // This class

private:

    // Pointer to container which is to be referenced
    //
    Container *data;

public:

    /**#: [Description="Set pointer to the container which is to be 
           referenced"] */

    void
    SetContainer(Container *d) { data = d; }

    /**#: [Description="Insert an item into the container and return an index
    	   reference to it. With this dereferencer new items are always inserted
    	   at the end of the container"] */
    
    reference_type
    Insert(const value_type& t)
    { return data->insert(data->end(),t) - data->begin(); }

    /**#: [Description="Remove an item from the container. Assumes that i is
     	   greater or equal to  zero and less than the size of the container."]
	   */
    void
    Remove(reference_type i) { data->erase(data->begin()+i); }

    /**#: [Description="When the data item referenced by i2 is removed, then i1
           may need updating"] */
    void
    Update(reference_type& toBeUpdated, const reference_type& toBeErased) const
    {
    	if (toBeUpdated > toBeErased) toBeUpdated--;
    }

    /**#: [Description="Convert an index reference into an iterator. (There is a
           const version of this function as well)."] */
    iterator
    GetIterator(reference_type i) { return data->begin() + i; }

    /**#: [Hidden] */
    const_iterator
    GetIterator(reference_type i) const { return data->begin() + i; }
};

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


//  [Description ="This class, together with BTL_InD provide a
//   uniform interface for the BTL_Graph class to the STL containers. This
//   class is used when individual items within a container are to be
//   referenced using an iterator. This is required for all STL containers
//   except the vector and deque as the iterators for these containers can
//   become invalid after insertions and erasures. "]
//  [Summary = "Allows stable dereferencing to the containers
//   via iterators. SHOULD NOT BE USED FOR vector AND deque LIKE CONTAINERS."]
//  [Usage = "This is template class the container type to be used must be
//   provided as the template parameter. SHOULD NOT BE USED FOR vector AND
//   deque LIKE CONTAINERS."]
//  [Authors = "W.R.Pitt"]
//  [Files = "<A HREF=./btl/Derefer.h>Derefer.h</A>"]
//  [Friends="None"]
//  [Dependencies="None"]



                /**#: [Hidden] */

template <class Container>
class BTL_ItD
{
public:

    typedef Container container_type;
    typedef BTL_TYPENAME Container::value_type value_type;
    typedef BTL_TYPENAME Container::iterator reference_type;
    typedef BTL_TYPENAME Container::iterator iterator;
    typedef BTL_TYPENAME Container::const_iterator const_iterator;
    typedef BTL_ItD<Container> Derefer;

private:

    // Pointer to container which is to be referenced
    //
    Container *data;

public:

    /**#: [Description="Set pointer to the container which is to be
           referenced"] */
    void
    SetContainer(Container *d) { data = d; }

    /**#: [Description="Insert an item into the container and return an iterator
	   reference to it"] */
    reference_type
    Insert(const value_type& t) { return data->insert(data->end(),t); }

    /**#: [Description="Remove an item from the container. Assumes that the
           input iterator is a valid reference to the container"] */
    void
    Remove(reference_type i) { data->erase(i); }

    /**#: [Description="The next 2 methods do nothing but are necessary to
 	   provide the same interface as the InD class. Providing
 	   the same interface enables interoperability."] */

    void
    Update(const reference_type& i1, const reference_type& i2) const
    { i1 == i2; }

    /**#: [Description="See Update description."] */
    iterator
    GetIterator(reference_type i) { return i; }
    
    /**#: [Hidden] */
    const_iterator
    GetIterator(reference_type i) const { return i; }
};

_BTL_END_NAMESPACE

#endif // Dereferencer.h
