//
// DFSIterator.h
//
// This file contains the BTL_DFSIterator class.
// This code is part of the Bioinformatics Template Library (BTL).
//
/// Copyright (C) 1997,1998 Birkbeck College, Malet Street, London WC1E 7HX,U.K.
// (classlib@mail.cryst.bbk.ac.uk)
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

// What if user wants to create new iterator that visits vertices not visited
// before (not connected to previous searches) ?

#if !defined(BTL_DFSITERATOR_H)
#define BTL_DFSITERATOR_H

#if !defined(SGI_CC)
#include <set>
#include <stack>
#else
#include <set.h>
#include <stack.h>
#endif

using namespace std;

_BTL_BEGIN_NAMESPACE


/**#: [Description ="
    This is a an iterator for iterating over a connected sub-graph in depth
    first order. It is designed to work with the BTL_Graph and
    BTL_GraphWithEdges but takes that type of graph as its template parameter
    and could be used on any graph and vertex type if they have a given set of
    functions."]
    [Summary = "an iterator for iterating over a connected sub-graph in depth
    first order."] 
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/DFSIterator.h>DFSIterator.h</A>"]
    [Dependencies="None"]

*/

template<class Graph>
class BTL_DFSIterator
{
    typedef BTL_TYPENAME Graph::VertexIterator	iterator;
    typedef BTL_TYPENAME Graph::Vertex   	    	Vertex;
    typedef BTL_TYPENAME Graph::value_type      	value_type;
    typedef Vertex*  	    	    	    	VertexPtr;

private:

    set<void*, less<void*> >	visited;

#if defined(GCC_272)
    stack< vector<void*> >   	notVisited;
#else
#if defined(GCC_280)
    stack<void*>    	    	notVisited;
#else
    stack< void*, vector<void*> > notVisited;
#endif
#endif

    VertexPtr 	    	    	current;

public:

	    	/**#: [Description="Default constructor."] */

    BTL_DFSIterator() { current = NULL; }

	    	/**#: [Description="Construct with iterator. This is essential
		       if operator++ etc. is to be used."] */

    BTL_DFSIterator(const iterator& i)
    {
    	// Initialise pointer to Vertex
    	//
    	current = i.GetPointer();
    }

	    	/**#: [Description="Move iterator forward one in depth first
		       search order."] */

    BTL_DFSIterator&
    operator++()
    {
    	if (current == NULL) return *this;

    	// Add current to visited set
    	//
    	visited.insert(visited.end(),current);

    	// Add any neighbours of current to unvisited stack
    	//
    	for (BTL_TYPENAME Vertex::VertexIterator v=current->begin();
	                                     v!=current->end(); v++)
    	    notVisited.push(v.GetPointer());

    	// Try and find the next Vertex to visit
    	//
    	bool nextFound=false;
        while (!notVisited.empty() && !nextFound)
        {
            // Get next vertex iterator off stack
            //
            current = (VertexPtr) notVisited.top();
            notVisited.pop();

            // Check that this vertex has not already been visited
            //
            nextFound = visited.find(current)==visited.end();
    	}

	// If an unvisited vertex is not found set current to NULL to indicate
	// termination
	//
	if (!nextFound)
	    current = NULL;

    	return *this;
    }

	    	/**#: [Description="Returns false when current position is one
		past the end of the last element in the container."] */

    // Return false if current == NULL - indicating termination, true otherwise
    //
    bool
    operator!() const { return current != NULL; }

	    	/**#: [Description="Move iterator forward one in depth first
		       search order."] */

    BTL_DFSIterator
    operator++(int) { BTL_DFSIterator tmp = *this; operator++(); return tmp; }

	    	/**#: [Description="Equality operator"] */

    bool
    operator==(BTL_DFSIterator i) const { return current==i.current; }

	    	/**#: [Description="Dereference operator"] */

    value_type&
    operator*() { return *(current->GetIterator()); }
};

/**#: [Description ="The const version of the BTL_DFSIterator."]
    [Summary = "the const version of the BTL_DFSIterator."]
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/DFSIterator.h>DFSIterator.h</A>"]
    [Dependencies="None"]

*/

template<class Graph>
class BTL_ConstDFSIterator
{
    typedef BTL_TYPENAME Graph::VertexIterator	iterator;
    typedef BTL_TYPENAME Graph::Vertex   	    	Vertex;
    typedef BTL_TYPENAME Graph::value_type      	value_type;
    typedef Vertex*  	    	    	    	VertexPtr;

private:

    set<void*, less<void*> >	visited;

#if defined(GCC_272)
    stack< vector<void*> >   	notVisited;
#else
#if defined(GCC_280)
    stack<void*>    	    	notVisited;
#else
    stack< void*, vector<void*> > notVisited;
#endif
#endif

    VertexPtr 	    current;

public:

	    	/**#: [Description="Default constructor."] */

    BTL_ConstDFSIterator() { current = NULL; }

	    	/**#: [Description="Construct with iterator. This is essential
		       if operator++ etc. is to be used."] */

    BTL_ConstDFSIterator(const iterator& i) 
    {
    	// Initialise pointer to Vertex  
    	//
    	current = i.GetPointer();
    }

	    	/**#: [Description="Move iterator forward one in depth first
		       search order."] */

    BTL_ConstDFSIterator&
    operator++() 
    {
    	if (current == NULL) return *this;
     
    	// Add current to visited set
    	//
    	visited.insert(visited.end(),current);
 
    	// Add any neighbours of current to unvisited stack
    	//
    	for (BTL_TYPENAME Vertex::VertexIterator v=current->begin(); 
	                                     v!=current->end(); v++)
    	    notVisited.push(v.GetPointer());

    	// Try and find the next Vertex to visit
    	//
    	bool nextFound=false;
        while (!notVisited.empty() && !nextFound)
        {
            // Get next vertex iterator off stack
            //
            current = (VertexPtr) notVisited.top();
            notVisited.pop();
 
            // Check that this vertex has not already been visited
            //
            nextFound = visited.find(current)==visited.end();
    	}
 
	// If an unvisited vertex is not found set current to NULL to indicate
	// termination
	//
	if (!nextFound)
	    current = NULL;

    	return *this;
    }
    
	    	/**#: [Description="Returns false when current position is one
		past the end of the last element in the container."] */

    // Return false if current == NULL - indicating termination, true otherwise
    //
    bool
    operator!() const { return current != NULL; }

	    	/**#: [Description="Move iterator forward one in depth first
		       search order."] */

    BTL_ConstDFSIterator
    operator++(int) 
    { 
    	BTL_ConstDFSIterator tmp = *this;
	operator++();
	return tmp; 
    }

	    	/**#: [Description="Equality operator"] */

    bool		
    operator==(BTL_ConstDFSIterator i) const { return current==i.current; }

	    	/**#: [Description="Dereference operator"] */

    value_type
    operator*() const { return *(current->GetIterator()); }
};

_BTL_BEGIN_NAMESPACE

#endif // BTL_DFSITERATOR_H
