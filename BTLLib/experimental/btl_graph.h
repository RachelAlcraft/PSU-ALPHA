
//
// btl_graph.h
//
// This file contains the initiation of the static member of BTL_VertexBase.
// This code is part of the Bioinformatics Template Library (BTL).
//
// Copyright (C) 1997, 1998 Birkbeck College, Malet Street, London WC1E 7HX,
// U.K. (classlib@mail.cryst.bbk.ac.uk)
//
// This library is free software; you can redistribute it and/or modify 
// it under the terms of the GNU Library General Public License as published 
// by the Free Software Foundation; either version 2 of the License, or 
// (at your option) any later version.  This library is distributed in the hope
// that it will be useful, but WITHOUT ANY WARRANTY; without even the
// implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU Library General Public License for more details.
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
///////////////////////////////////////////////////////////////////////////


#if !defined(BTL_GRAPH_H)
#define BTL_GRAPH_H


#include "BTL.h"
#include "VoidPtrStore.h"
#include "Dereferencer.h"

using namespace std;

_BTL_BEGIN_NAMESPACE

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

// Compare 'object' used for sorting vertices and edges by their id.
// Should be a member function but template member
// functions are not supported by many compilers..

                /**#: [Hidden] */

template <class T>
struct BTL_ById : binary_function<T, T, bool> {
    bool operator()(T x, T y) const { return x->GetId() < y->GetId(); }
};


//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

//     [Description ="This abstract base class is used to generate an individual integer
//       identifier for BTL_Vertex objects."]
//      [Summary = "base class for BTL_Vertex. Provides access to each object's
//       id."]
//      [Authors = "W.R.Pitt"]
//      [Files = "<A HREF=./btl/btl_graph.h>btl_graph.h</A>"]
//      [Friends="operator&lt&lt and operator=="]
//      [Dependencies="None."]

	    	/**#: [Hidden] */

class BTL_VertexBase
{
private:

    static unsigned int nextId;     // The id number for the next vertex
    unsigned int    	id; 	    // The id number for this vertex

	    	/**#: [Hidden] */
    virtual void
    Print(ostream& os) const = 0;

public:

	    	/**#: [Description="Construct vertex with a new id."] */

    BTL_VertexBase() { id = nextId++; }

	    	/**#: [Description="Copy constructor"] */

    BTL_VertexBase(const BTL_VertexBase& v) { id = v.id; }

	    	/**#: [Description="Read the id of this vertex."] */

    unsigned int
    GetId() const { return id; }

	    	/**#: [Description="Equality operator."] */

    friend bool
    operator==(const BTL_VertexBase& v1, const BTL_VertexBase& v2) 
    { return v1.GetId() == v2.GetId(); }

	    	/**#: [Description="Output operator."] */

    friend ostream& operator<<(ostream& os, const BTL_VertexBase& v)
    	{ v.Print(os); return os; }

};


//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

/**#: [Description ="
    This class encapsulates a graph vertex and is designed for use with the
    graph class. Each vertex contains a reference to a data object and
    pointers to other vertices. It is a template class the single template
    parameter is a type of dereferencer."]
    [Summary = "a template graph vertex for use with graph"] 
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/btl_graph.h>btl_graph.h</A>"]
    [Dependencies="<A HREF=#BTL_II>BTL_II</A>, and 
                   <A HREF=#BTL_II>BTL_IteratorDerefencer</A>, or 
		   <A HREF=#BTL_II>BTL_InDencer</A>"]

*/

template <class T> class BTL_Graph;

template <class Dereferencer> 
class BTL_Vertex : public BTL_VertexBase
{
public:

#if !defined(SGI_CC)
    typedef BTL_TYPENAME Dereferencer::container_type   container_type;
    typedef BTL_TYPENAME Dereferencer::value_type 	    value_type;
    typedef BTL_TYPENAME Dereferencer::reference_type   reference_type;
    typedef BTL_TYPENAME Dereferencer::iterator  	    iterator;
    typedef BTL_TYPENAME Dereferencer::const_iterator   const_iterator; 
#else
    typedef Dereferencer::container_type   	    container_type;
    typedef Dereferencer::value_type	   	    value_type;
    typedef Dereferencer::reference_type   	    reference_type;
    typedef Dereferencer::iterator	   	    iterator;
    typedef Dereferencer::const_iterator    	    const_iterator;
#endif
    typedef BTL_Vertex<Dereferencer>   	    	    Vertex;
    typedef Vertex* 	    	    	    	    VertexPtr;
    typedef BTL_PtrSet<Vertex>	    	    	    VertexPtrStore;
    typedef BTL_II<Vertex>     	    	    	    VertexIterator;
    typedef BTL_CII<Vertex>  	    	    	    ConstVertexIterator;

    friend class BTL_Graph<Dereferencer>;
    
private:

    Dereferencer    	*deref;     // Pointer to deferencer object stored in
    	    	    	    	    // graph object that this vertex is a part
    reference_type   	ref;	    // Reference to a data item - the contents
    	    	    	    	    // of this vertex
    VertexPtrStore    	links;	    // container of pointers to vertices which
    	    	    	    	    // are linked to this vertex via edges of
				    // a graph
    bool    	    	visited;    // visited flag used by depth first search
    	    	    	    	    // algorithms and the like
        
protected:
		// function	    	
	    	// Add pointer/link to another vertex. Return false if link
	    	// already there or if it points to this  vertex. Takes O(logN)

	    	/**#: [Hidden] */
    bool
    AddLink(VertexPtr p) 
    {
	if (p == this) return false; 	// Avoids self links.
	return links.insert(links.end(), p) != links.end(); 
    }
    	    	// Remove all links to this vertex present in other vertices
		// N.B. This should be used for undirected graphs only.
		
	    	/**#: [Hidden] */
    void
    RemoveLinks()
    {
    	for (VertexIterator i = begin_vertex(); i!= end_vertex(); i++)
    	    (*i).RemoveLink(this);
    }
	    	// Remove pointer/link to other vertex. Return false if link 
	    	// not found. Takes O(logN). 
	    
	    	/**#: [Hidden] */
    bool
    RemoveLink(VertexPtr p) 
    {
	VertexPtrStore::iterator search = links.find(p);
	if (search == links.end()) return false;
	links.erase(search);
	return true;
    }
    	    	// Update reference to vertex data item. This is only necessary
		// for vertices using an InDencer
    
	    	/**#: [Hidden] */
    void
    UpdateReference(VertexPtr p) 
    {
	deref->Update(ref,p->ref);
    }
    
    	    	// Return the reference to the data item associated with this
		// vertex. The type of reference could be an iterator or an
		// index number.
		    
	    	/**#: [Hidden] */
    reference_type
    GetRef()
    {
    	return ref;
    }
    	    	// Print vertex contents and links
    	    	//
		
	    	/**#: [Hidden] */
    virtual void
    Print(ostream& os) const
    {
    	os << "\nVertex " << GetId()
    	   << "\nContents: " << *GetIterator()
    	   << "\nConnected to vertices: ";

    	typedef set<unsigned int, less<unsigned int> > Ids;
    	Ids ids;
    	for (ConstVertexIterator i=begin_vertex(); i!= end_vertex(); i++)
    	    ids.insert(ids.end(),(*i).GetId());
	    
    	for (Ids::iterator i=ids.begin(); i!=ids.end(); i++)
    	    os << *i << " ";
    	os << '\n';
    }

	    	// Construct Vertex belonging particular graph and associated
		// with a particular data item stored within this graph

	    	/**#: [Hidden] */
	    	
    BTL_Vertex(Dereferencer& d, reference_type r) : BTL_VertexBase() 
    { 
    	visited = false;
    	deref = &d; 
	ref = r; 
    }
    

public:

	    	/**#: [Description="Construct new empty Vertex"] */

    BTL_Vertex() : BTL_VertexBase() { visited = false; }


	    	/**#: [Description="Copy constructor."] */

    BTL_Vertex(const BTL_Vertex& v) : BTL_VertexBase(v) 
    {  
    	deref = v.deref;
    	ref = v.ref;
    	links = v.links;
    	visited = v.visited;
    }
    	    	    	       
    // N.B. The Graph deletes the data item referenced by this vertex. 

	    	/**#: [Description="Destructor"] */
    virtual
    ~BTL_Vertex() {}
    
    
	    	/**#: [Description="Get an iterator that references the data
		item associated with this vertex. With graphs that use STL
		an vector or deque to store data items this iterator may become
		invalid after insertions or deletions of vertices. (There is also a
      const version of this function.)"] */
    iterator
    GetIterator() { return deref->GetIterator(ref); }

    const_iterator
    GetIterator() const { return deref->GetIterator(ref); }

	    	/**#: [Description="Read/write access to visited flag. This
	    	       can be used in Graph searching routines."] */
    bool&
    Visited() { return visited; }

	    	/**#: [Description="Return iterator that references first
	    	       vertex connected to this vertex.(There is a
                       const version of this function as well)"] */
    virtual VertexIterator
    begin_vertex() { return links.begin(); }

    virtual ConstVertexIterator
    begin_vertex() const { return links.begin(); }

	    	/**#: [Description="Return iterator that references the
	    	       position in memory after the last vertex connected to
		       this vertex.(There is a const version of this function
		       as well)"] */
    virtual VertexIterator
    end_vertex() { return links.end(); }

    virtual ConstVertexIterator
    end_vertex() const { return links.end(); }

	    	/**#: [Description="Return iterator that references first
	    	       vertex connected to this vertex. (There is a
                       const version of this function as well)"] */
    VertexIterator
    begin() { return links.begin(); }

    ConstVertexIterator
    begin() const { return links.begin(); }

	    	/**#: [Description="Return iterator that references the
	    	       position in memory after the last vertex connected to
		       this vertex.(There is a const version of this function
		       as well)"] */
    VertexIterator
    end() { return links.end(); }

    ConstVertexIterator
    end() const { return links.end(); }
};

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


/**#: [Description ="
    Graph container with any type of data at vertices of the graph and no data
    associated with edges. The default behaviour is have undirected edges.
    This can be changed to directed if specified in the constructor. No
    self links (edges from one vertex to itself) are allowed."]
    [Summary = "A template graph with edges as simple pointers to vertices."]
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/Graph.h>Graph.h</A>"]
    [Friends="operator&lt&lt"]
    [Dependencies="<A HREF=#BTL_Vertex>BTL_VertexBase and its subclasses.</A>"]

*/

template<class Dereferencer>
class BTL_Graph
{
public:

    typedef BTL_TYPENAME Dereferencer::container_type   container_type;
    typedef BTL_TYPENAME Dereferencer::value_type 	    value_type;
    typedef BTL_TYPENAME Dereferencer::reference_type   reference_type;
    typedef BTL_TYPENAME Dereferencer::iterator  	    iterator;
    typedef BTL_TYPENAME Dereferencer::const_iterator   const_iterator;
    typedef BTL_Vertex<Dereferencer>  	    	    Vertex;
    typedef BTL_VertexBase*  	    	    	    VertexBasePtr;
    typedef Vertex* 	    	    	    	    VertexPtr;
    typedef BTL_PtrSet<Vertex>      	    	    VertexPtrStore;
    typedef BTL_II<Vertex>   	    		    VertexIterator;
    typedef BTL_CII<Vertex>  	    		    ConstVertexIterator;
    typedef VertexPtrStore::size_type	    	    size_type;

    enum Directionality { directed, undirected };   // The two types
    	    	    	    	    	    	    // graph allowed
private:

    Dereferencer	deref;	    	// Object used to reference vertex data
    container_type  	vertexData;	// Vertex data
    Directionality      directionality; // Type of edges
    VertexPtrStore  	vertices;       // Vertices within this graph

protected:

    // Check that a pointer points to a vertex in this graph O(logN)
    //
	    	/**#: [Hidden] */
    bool
    ValidVertexPointer(VertexPtr v)
    {
    	return vertices.find(v) != vertices.end();
    }
    
    // Remove link between two vertices. Return false if vertices are not part
    // of this graph or if no edge exists between them.
    //
	    	/**#: [Hidden] */
    bool
    RemoveLink(VertexPtr p1, VertexPtr p2)
    {
    	if (p1 == p2) return false;
   	if (!ValidVertexPointer(p1) || !ValidVertexPointer(p2)) return false;
    	if (!p1->RemoveLink(p2)) return false; 
    	if (directionality == undirected) p2->RemoveLink(p1);
    	return true;
    }
    
    // Remove vertex from this graph. Return false if vertex is not part of this
    // graph.
	    	/**#: [Hidden] */
    bool
    RemoveVertex(VertexPtr p)
    {
    	// Check that vertex is in this graph
    	//
    	VertexPtrStore::iterator target = vertices.find(p);
    	if (target == vertices.end()) return false;
    	
    	// Remove all links to the vertex
    	//	
    	if (directionality==undirected)
    	    p->RemoveLinks();
    	else
    	{
    	    for (VertexIterator i=begin_vertex(); i!=end_vertex(); i++)
    	    	(*i).RemoveLink(p);
    	}    	
    	
	// Update references from vertices to data. This is needed only when
	// indices are used instead of iterators. Does nothing (over and over
	// again when iterator references are used)
    	
    	for (VertexIterator i=begin_vertex(); i!=end_vertex(); i++)
    	    (*i).UpdateReference(p);
	    
    	// Remove data item
    	//
   	deref.Remove(p->GetRef());
    	
    	// Delete vertex object
    	//    
    	delete p;
	
	// Remove vertex from list of vertices within this graph
	//    
    	vertices.erase(target);
    	
    	return true;
    }	
     
	    	/**#: [Hidden] */
    reference_type
    InsertVertexItem(const value_type& v)
    {
   	 return deref.Insert(v);
    }

	    	/**#: [Hidden] */
    Dereferencer&
    GetDerefObj()
    {
    	return deref;
    }

	    	/**#: [Hidden] */
    VertexIterator
    InsertVertex(VertexPtr p)
    {
    	return vertices.insert(vertices.end(),p);
    }

	    	// Virtual method from actually doing the printing
	    	//
	    	/**#: [Hidden] */
    virtual void
    Print(ostream& os) const
    {
    	os << "\nVERTICES:\n";

    	// Print each vertex in the order they were created
    	//
	typedef set< VertexBasePtr, BTL_ById<VertexBasePtr> > SortedSet;
	SortedSet sortedSet;

    	for (ConstVertexIterator i = begin_vertex(); i!=end_vertex();i++)
	    sortedSet.insert(sortedSet.end(),i.GetPointer());

	for (SortedSet::iterator i=sortedSet.begin(); i!=sortedSet.end(); i++)
    	    os << **i;

    	os << '\n';
    }

public:

	    	/**#: [Description="Construct new undirected graph"] */

    BTL_Graph() { deref.SetContainer(&vertexData); directionality = undirected; }

	    	/**#: [Description="Construct new undirected or 
	    	       directed graph i.e. BTL_Graph g1(directed); or BTL_Graph
	    	       g2(undirected); "] */

    BTL_Graph(const Directionality& e) : directionality(e) 
    { 
    	deref.SetContainer(&vertexData); 
    }

	    	/**#: [Description="Delete Graph and all the vertices in it"] */
    virtual 
    ~BTL_Graph() 
    { 
    	VertexPtrStore::iterator i;
    	for (i=vertices.begin(); i != vertices.end(); i++)
    	    delete *i;
    }

	    	/**#: [Description="Get the edge directionality of this
		       graph."] */
    Directionality
    GetDirectionality() const { return directionality; }

	    	/**#: [Description="Add new vertex to graph."] */
    virtual VertexIterator
    AddVertex(const value_type& v)
    {
     	// Create new Vertex object with no links.
    	//
	reference_type ref = InsertVertexItem(v);
    	VertexPtr p = new Vertex(deref,ref);
    	if (p == NULL) FATAL_ERROR("No more memory.");
    	return vertices.insert(vertices.end(),p);
    }

	    	/**#: [Description="Remove a vertex (and all pointers to it)
	    	       from the graph. Return false if vertex was not found
	    	       in the graph."] */
    virtual bool
    RemoveVertex(VertexIterator in)
    {

    	VertexPtr p = in.GetPointer();
	return RemoveVertex(p);
    }

	    	/**#: [Description="Add a pointer between two vertices.
	    	       Return false if either vertex is not in this graph or
	    	       is link already present"] */
    bool
    AddLink(VertexIterator i1, VertexIterator i2)
    {
    	VertexPtr p1 = i1.GetPointer(), p2 = i2.GetPointer();
    	if (!ValidVertexPointer(p1) || !ValidVertexPointer(p2)) return false;

    	// Try add link(s) return false if it is already present or if
    	// p1 == p2
	//
    	if (!p1->AddLink(p2)) return false;
    	if (directionality == undirected) p2->AddLink(p1);
    	return true;
    }
	    	/**#: [Description="Remove pointer/link between two vertices.
	    	       Return false if link not there"] */
    bool
    RemoveLink(VertexIterator i1, VertexIterator i2)
    {
     	VertexPtr p1 = i1.GetPointer(), p2 = i2.GetPointer();
	return RemoveLink(p1,p2);
    }

	    	/**#: [Description="Set the visited field of all vertices in
	    	       this graph."] */
    virtual void
    SetAllVisited(bool truth)
    {
    	for (VertexIterator i=begin_vertex(); i!=end_vertex(); i++)
    	    (*i).Visited() = truth;
    }
	    	/**#: [Description="Find a vertex within this Graph.(There is a
		       const version of this function as well)"] */
    virtual VertexIterator
    find(Vertex& v)
    {
    	return vertices.find(&v);
    }

    virtual ConstVertexIterator
    find(Vertex& v) const
    {
    	return vertices.find(&v);
    }
	    	/**#: [Description="Return iterator that references first
	    	       data item.(There is a const version of this function
		       as well)"] */
    iterator
    begin() { return vertexData.begin(); }

    const_iterator
    begin() const { return vertexData.begin(); }

	    	/**#: [Description="Return iterator that references the
	    	       position in memory after the last data item. (There is a
		       const version of this function as well) "] */
    iterator
    end() { return vertexData.end(); }

    const_iterator
    end() const { return vertexData.end(); }

	    	/**#: [Description="Return iterator that references the first
	    	       Vertex. (There is a const version of this function
		       as well)"] */
    virtual VertexIterator
    begin_vertex() { return vertices.begin(); }

    virtual ConstVertexIterator
    begin_vertex() const { return vertices.begin(); }

	    	/**#: [Description="Return iterator that references the
	    	       position in memory after the last Vertex. (There is a
		       const version of this function as well)"] */
    virtual VertexIterator
    end_vertex() { return vertices.end(); }

    virtual ConstVertexIterator
    end_vertex() const { return vertices.end(); }

	    	/**#: [Description="Give the number of vertices within the
	    	       Graph"] */
    size_type
    size() const { return vertices.size(); }

	    	/**#: [Description="Print list of vertices in this graph and
		       their connections."] */
    friend ostream&
    operator<<(ostream& os, const BTL_Graph& g) { g.Print(os); return os; }

};
	    	/**#: [Hidden] */

template <class Container>
class BTL_Graph1 : public
BTL_Graph< BTL_ItD<Container> >
{};
	    	/**#: [Hidden] */

template <class Container>
class BTL_Graph2 : public
BTL_Graph< BTL_InD<Container> >
{};

_BTL_END_NAMESPACE

#endif // Graph.h

unsigned int BTL_VertexBase::nextId = 1;
