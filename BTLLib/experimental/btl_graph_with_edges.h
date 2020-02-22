
//
// GraphWithEdges.h
//
// This file contains the BTL_GraphWithEdges, BTL_VertexWithEdges, and BTL_Edge
// classes.
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


#if !defined(BTL_GRAPHWITHEDGES_H)
#define BTL_GRAPHWITHEDGES_H

#include "btl_graph.h"

_BTL_BEGIN_NAMESPACE

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

/**#: [Description ="This class represents an vertex that has links to edge
       objects. It is to be used in combination with class BTL_GraphWithEdges."]
      [Summary = "template graph vertex for use with BTL_GraphWithEdges"] 
      [Authors = "W.R.Pitt"]
      [Files = "<A HREF=./btl/btl_graph_with_edges.h>btl_graph_with_edges.h</A>"]
      [Friends="friend class BTL_GraphWithEdges&ltT1,T2&gt>; <P>
       friend ostream& operator&lt&lt(ostream &os, 
       const BTL_VertexWithEdges &m);"]
      [Dependencies="None, except by inheritance."]

*/


template <class T1, class T2> class BTL_GraphWithEdges;
template <class T1, class T2> class BTL_Edge;

template <class VertexDereferencer, class EdgeDereferencer>
class BTL_VertexWithEdges : public BTL_Vertex<VertexDereferencer>
{
public:

    typedef BTL_TYPENAME VertexDereferencer::reference_type     reference_type;

    typedef BTL_Vertex<VertexDereferencer>    	    	    	    	Base;
    typedef BTL_VertexWithEdges<VertexDereferencer,EdgeDereferencer>	Vertex;
    typedef BTL_Edge<VertexDereferencer,EdgeDereferencer>	    	Edge;
    typedef Edge*   	    	EdgePtr;
    typedef Vertex*   	    	VertexPtr;
    typedef BTL_PtrSet<Edge>    EdgePtrStore;
    typedef BTL_II<Edge>        EdgeIterator;
    typedef BTL_CII<Edge>       ConstEdgeIterator;

    friend class BTL_GraphWithEdges<VertexDereferencer, EdgeDereferencer>;

private:


    EdgePtrStore    connectedEdges;
    
	    	// Connect this vertex with an edge.
	    	//
	    	/**#: [Hidden] */
    bool
    AddEdge(EdgePtr edge) 
    { 
    	return connectedEdges.insert(connectedEdges.end(),edge) !=
	       connectedEdges.end(); 
    }

	    	
	    	// Disconnect this vertex from an edge.
	    	//
	    	/**#: [Hidden] */
    bool
    RemoveEdge(EdgePtr edge)
    {
	EdgePtrStore::iterator search = connectedEdges.find(edge);
	if (search == connectedEdges.end()) return false;
	connectedEdges.erase(search);
	return true;
    }
	    	// Construct Vertex belonging particular graph

	    	/**#: [Hidden] */
	    	
    BTL_VertexWithEdges(VertexDereferencer& d, reference_type r) : Base(d,r) {} 

public:

	    	/**#: [Description="Constructor empty vertex."] */
	    	
    BTL_VertexWithEdges() : Base() {}

	    	/**#: [Description="Delete vertex and its contents."] */
	    	
    virtual
    ~BTL_VertexWithEdges() {}


	    	/**#: [Description="Return iterator that references first
	    	       edge connected to Vertex. (There is also a const version
                of this function.)"] */
    EdgeIterator
    begin_edge() { return connectedEdges.begin(); }

    ConstEdgeIterator
    begin_edge() const { return connectedEdges.begin(); }

	    	/**#: [Description="Return iterator that references the
	    	       position in memory after the last edge connected to
	    	       Vertex. (There is also a const version
                of this function.)"] */
    EdgeIterator
    end_edge() { return connectedEdges.end(); }

    ConstEdgeIterator
    end_edge() const { return connectedEdges.end(); }

	    	// Virtual method to print vertex and its connections.
	    	//
	    	/**#: [Hidden] */
    virtual void
    Print(ostream& os) const
    {
    	Base::Print(os);

    	typedef set<unsigned int, less<unsigned int> > Ids;
    	Ids ids;

    	for(ConstEdgeIterator i=begin_edge(); i!=end_edge(); i++)
    	    ids.insert(ids.end(),(*i).GetId());

    	os  << "Connected to edges: ";

    	for (Ids::iterator i=ids.begin(); i!=ids.end(); i++)
    	    os << *i << " ";

    	os << '\n';
    }

}; // class BTL_VertexWithEdges

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

// [Description ="This class is used to generate an individual integer
//       identifier for BTL_Edge objects. It is necessary as static member
//       templates are not implemented in gcc 2.7.2, but should become obsolete
//       before very long ."]
//      [Summary = "base class for BTL_Edge providing access to an object's id."]
//      [Authors = "W.R.Pitt"]
//      [Files = "<A HREF=./btl/btl_graph_with_edges.h>btl_graph_with_edges.h</A>"]
//      [Dependencies="None."]

	    	/**#: [Hidden] */

class BTL_EdgeBase
{
private:

    static unsigned int nextId;
    unsigned int    	id;

	    	/**#: [Hidden] */
    virtual void
    Print(ostream& os) const = 0;

public:

	    	/**#: [Description="Construct vertex with a new id."] */

    BTL_EdgeBase() { id = nextId++; }

	    	/**#: [Description="Read the id of this vertex."] */

    unsigned int
    GetId() const { return id; }

	    	/**#: [Description="Equality operator."] */

    friend bool
    operator==(const BTL_EdgeBase& e1, const BTL_EdgeBase& e2)
    { return e1.GetId() == e2.GetId(); }

	    	/**#: [Description="Output operator."] */

    friend ostream& operator<<(ostream& os, const BTL_EdgeBase& e)
    	{ e.Print(os); return os; }

};

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

/**#: [Description ="This class represents an edge of a graph. It is to be used
       in combination with class BTL_GraphWithEdges."]
      [Summary = "a template graph edge for use with BTL_GraphWithEdges"]
      [Authors = "W.R.Pitt"]
      [Files = "<A HREF=./btl/btl_graph_with_edges.h>btl_graph_with_edges.h</A>"]
      [Friends="friend class BTL_GraphWithEdges&ltT1,T2&gt>; <P>
       friend ostream& operator&lt&lt(ostream &os,
       const BTL_Edge &m);"]
      [Dependencies="<A HREF=#BTL_VertexWithEdges>BTL_VertexWithEdges</A>"]

*/

template <class VertexDereferencer, class EdgeDereferencer>
class BTL_Edge : public BTL_EdgeBase
{
public:

    typedef EdgeDereferencer 	    	    	    	dereferencer_type;
    typedef BTL_TYPENAME EdgeDereferencer::container_type 	container_type;
    typedef BTL_TYPENAME EdgeDereferencer::value_type    	value_type;
    typedef BTL_TYPENAME EdgeDereferencer::reference_type 	reference_type;
    typedef BTL_TYPENAME EdgeDereferencer::iterator      	iterator;
    typedef BTL_TYPENAME EdgeDereferencer::const_iterator 	const_iterator;

    typedef BTL_VertexWithEdges<VertexDereferencer,EdgeDereferencer> VertexWE;
    typedef BTL_Edge<VertexDereferencer,EdgeDereferencer>   	     Edge;
    typedef VertexWE* 	    	    VertexWEPtr;
    typedef Edge* 	    	    EdgePtr;

    friend class BTL_GraphWithEdges<VertexDereferencer,EdgeDereferencer>;

private:

    EdgeDereferencer    *deref;     // Pointer to deferencer object stored in
    	    	    	    	    // graph object that this vertex is a part
    reference_type   	ref;	    // Reference to a data item - the contents
    	    	    	    	    // of this vertex
    VertexWEPtr	    	leftVertex;
    VertexWEPtr	    	rightVertex;


	    	/**#: [Hidden] */
    void
    UpdateReference(EdgePtr i)
    {
	// Update reference to data if necessary
	//
	deref->Update(ref,i->ref);
    }

	    	/**#: [Hidden] */
    virtual void
    Print(ostream& os) const
    {
    	os  << "\nEdge " << GetId()
    	    << "\nContents: " << *GetIterator()
    	    << "\nConnects vertices: "
    	    << leftVertex->GetId()
    	    << " "
    	    << rightVertex->GetId()
    	    << '\n';
    }
	    	// Construct Edge belonging particular graph

	    	/**#: [Hidden] */

    BTL_Edge(EdgeDereferencer& d, reference_type r,
             VertexWEPtr p1, VertexWEPtr p2) : BTL_EdgeBase()
    {
    	deref = &d;
	ref = r;
	leftVertex = p1;
	rightVertex = p2;
    }

public:

	    	/**#: [Description="Construct empty edge"] */

    BTL_Edge() : BTL_EdgeBase() {}

	    	/**#: [Description="Copy constructor"] */

    BTL_Edge(const BTL_Edge& e) : BTL_EdgeBase(e)
    {
    	deref = e.deref;
    	ref = e.ref;
	leftVertex = e.leftVertex;
	rightVertex = e.rightVertex;
    }

    // N.B. The Graph deletes the data item referenced by this edge.

	    	/**#: [Description="Destructor"] */
    virtual
    ~BTL_Edge() {}

    	    	/**#: [Description="Returns an iterator that can be deferenced
		       to obtain the data object associated with this edge
		       (There is a const version of this function as well)"] */
    iterator
    GetIterator() { return deref->GetIterator(ref); }

    const_iterator
    GetIterator() const { return deref->GetIterator(ref); }

	    	/**#: [Description="Returns pointer to the first vertex
	    	       connected to this edge."] */
    VertexWEPtr
    LeftVertex() const { return leftVertex; }

	    	/**#: [Description="Returns pointer to the second vertex
	    	       connected to this edge."] */
    VertexWEPtr
    RightVertex() const { return rightVertex; }

};

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


/**#: [Description ="
    Graph container with any type of data at the vertices of the graph and any
    type of data at the edges. The default behaviour is have undirected
    edges. This can be changed to directed if specified in the
    constructor. No self links (edges from one vertex to itself) are allowed."]
    [Summary = "a template graph with two template parameters determining the
     type of vertices and edges"]
    [Authors = "W.R.Pitt"]
    [Files = "<A HREF=./btl/btl_graph_with_edges.h>btl_graph_with_edges.h</A>"]
    [Friends="friend ostream&
              operator&lt&lt<(ostream &os, const BTL_GraphWithEdges &m);"]
    [Dependencies="<A HREF=#BTL_Vertex>BTL_Vertex</A>,
                   <A HREF=#BTL_II>BTL_II</A>"]

*/

template <class VertexDereferencer, class EdgeDereferencer>
class BTL_GraphWithEdges : public BTL_Graph<VertexDereferencer>
{
public:

    typedef BTL_Edge<VertexDereferencer,EdgeDereferencer>   Edge;
    typedef BTL_TYPENAME Edge::container_type    	EdgeContainer;
    typedef BTL_TYPENAME Edge::iterator  	    	EdgeDataIterator;
    typedef BTL_TYPENAME Edge::const_iterator  	    	EdgeDataConstIterator;
    typedef BTL_TYPENAME Edge::value_type 	    	EdgeData;
    typedef BTL_TYPENAME Edge::reference_type    	EdgeReference;
    typedef BTL_TYPENAME VertexDereferencer::value_type value_type;
    typedef BTL_Graph<VertexDereferencer>	Base;
    typedef Edge*   	    	    	    	EdgePtr;
    typedef BTL_PtrSet<Edge>        	    	EdgePtrStore;
    typedef BTL_II<Edge>       			EdgeIterator;
    typedef BTL_CII<Edge>      			ConstEdgeIterator;

    typedef BTL_VertexWithEdges<VertexDereferencer,EdgeDereferencer> VertexWE;
    typedef VertexWE* VertexWEPtr;

private:

    EdgeDereferencer	deref;
    EdgeContainer  	data;
    EdgePtrStore	edges;

    // Delete pointer to edge from edge store (but don't delete object).
    //
	    	/**#: [Hidden] */
    bool
    RemoveEdgePtrFromStore(EdgePtr edge)
    {
    	EdgePtrStore::iterator i = edges.find(edge);
    	if (i == edges.end()) return false;
    	edges.erase(i);
    	return true;
    }


	    	/**#: [Hidden] */
    void
    ConnectNewEdge(EdgePtr e, VertexWEPtr left, VertexWEPtr right)
    {
    	// Connect the two vertices to the new edge
    	//
    	left->AddEdge(e);
    	right->AddEdge(e);
    }

    // Delete all edges connected to a particular vertex
    //
	    	/**#: [Hidden] */
    void
    DeleteVertexEdges(VertexWEPtr v)
    {
    	typedef BTL_TYPENAME VertexWE::EdgeIterator EIterator;
	EdgePtrStore vertex_edges;
	
	// Store pointers to edges connected to the vertex. Edges are deleted
	// from the vertex's own container by DeleteEdge
	//
    	for (EIterator i = v->begin_edge(); i!= v->end_edge(); i++)
    	{
    	    EdgePtr edgePtr = i.GetPointer();
	    vertex_edges.insert(vertex_edges.end(),edgePtr);
    	}
	for (EdgePtrStore::iterator i=vertex_edges.begin();
	                           i!=vertex_edges.end(); i++)
    	    DeleteEdge((EdgePtr) *i);
    }

    // Delete Edge
    //
	    	/**#: [Hidden] */
    void
    DeleteEdge(EdgePtr edgePtr)
    {
    	// Remove pointers to this edge from connected vertices
	//
        edgePtr->leftVertex->RemoveEdge(edgePtr);
        edgePtr->rightVertex->RemoveEdge(edgePtr);
	
	// Remove this edge from the list of edges in this graph
	//
    	RemoveEdgePtrFromStore(edgePtr);

	// Update references from edges to data. This is needed only when
	// indices are used instead of iterators. Does nothing (over and over
	// again when iterator references are used)
    	
    	for (EdgeIterator i=begin_edge(); i!=end_edge(); i++)
    	    (*i).UpdateReference(edgePtr);
	    
    	// Remove edge data item
    	//
   	deref.Remove(edgePtr->ref);

    	// Delete edge object
	//
    	delete edgePtr;
    }


    // Make link functions private to enforce use of edge functions

	    	/**#: [Hidden] */
    bool
    AddLink(VertexIterator i1, VertexIterator i2) 
    { 
    	return Base::AddLink(i1,i2);
    }

	    	/**#: [Hidden] */
    bool
    RemoveLink(VertexIterator i1, VertexIterator i2)
    { 
    	return Base::RemoveLink(i1,i2);
    }

public:


	    	/**#: [Description="Construct empty undirected graph."] */

    BTL_GraphWithEdges() : Base() { deref.SetContainer(&data); }

	    	/**#: [Description="Construct new undirected or
	    	       directed graph i.e.
	    	       BTL_GraphWithEdges g1(directed); or BTL_Graph
	    	       g2(undirected); "] */

    BTL_GraphWithEdges(const Directionality& e) : Base(e) 
    {
    	deref.SetContainer(&data);
    }


	    	/**#: [Description="Delete graph and all vertices and edges in
	    	       it."] */
    virtual
    ~BTL_GraphWithEdges()
    {
    	// Delete all the Edges before the Graph itself is deleted.
    	//
    	EdgePtrStore::iterator i;
    	for (i=edges.begin(); i!=edges.end(); i++)
    	    delete *i;
    }

	    	/**#: [Description="Add new vertex to graph."] */
    virtual VertexIterator
    AddVertex(const value_type& v)
    {
    	// Insert data item into vertex data container
	//
	reference_type ref = InsertVertexItem(v);

     	// Create new VertexWithEdges object with no links.
    	//
    	VertexWEPtr p = new VertexWE(GetDerefObj(),ref);
    	if (p == NULL) FATAL_ERROR("No more memory.");

	// Insert new VertexWithEdges into vertex store (implicit up cast from
	// VertexWithEdge* to Vertex*).
	//
    	return InsertVertex(p);
    }

	    	/**#: [Description="Delete vertex from graph. Return false
	    	       if vertex not found in this graph."] */
    virtual bool
    RemoveVertex(VertexIterator in)
    {
    	VertexPtr v = in.GetPointer();

    	// Check that input is valid
    	//
    	if (!ValidVertexPointer(v)) return false;

    	// Delete connected edges
    	//
#if defined(GCC_272)
	VertexWEPtr vwe = (VertexWEPtr) v;
#else
	VertexWEPtr vwe = dynamic_cast<VertexWEPtr>(v);
#endif
    	DeleteVertexEdges(vwe);

    	// Delete vertex.
    	//
    	Base::RemoveVertex(v);

    	return true;
    }

	    	/**#: [Description="Insert edge into this graph. Return
	    	       endEdge() if input is invalid."] */
    EdgeIterator
    AddEdge(VertexIterator i1, VertexIterator i2, const EdgeData& d)
    {
    	// Add link(s) between vertices. This returns false if both iterators
	// reference the same vertex.
    	//
    	if (!AddLink(i1,i2)) return edges.end();

    	VertexPtr v1 = i1.GetPointer(), v2 = i2.GetPointer();

#if defined(GCC_272)
	VertexWEPtr vwe1 = (VertexWEPtr) v1;
	VertexWEPtr vwe2 = (VertexWEPtr) v2;
#else
	VertexWEPtr vwe1 = dynamic_cast<VertexWEPtr>(v1);
	VertexWEPtr vwe2 = dynamic_cast<VertexWEPtr>(v2);
#endif

     	// Create new Edge object.
    	//
   	EdgeReference ref = deref.Insert(d);
    	EdgePtr e = new Edge(deref,ref,vwe1,vwe2);
    	if (e == NULL) FATAL_ERROR("No more memory.");

    	// Make connection from connected vertices to this edge
	//
     	ConnectNewEdge(e, vwe1, vwe2);

    	// Store pointer to edge
    	//
    	return edges.insert(edges.end(),e);
    }


	    	/**#: [Description="Remove edge from this graph. Return
	    	       false if input is invalid."] */
    bool
    RemoveEdge(EdgeIterator edgeIter)
    {
    	// Convert iterator to pointer
    	//
    	EdgePtr edgePtr = edgeIter.GetPointer();

    	// Check edge belongs to this graph
	//
	if (edges.find(edgePtr) == edges.end()) return false;

    	// Find the two vertices connected by edge
    	//
    	VertexWEPtr v1 = edgePtr->leftVertex,
    	            v2 = edgePtr->rightVertex;

    	// Remove link in BTL_Graph (implicit up cast of vertex pointers)
    	//
    	if (!Base::RemoveLink(v1,v2))
    	    FATAL_ERROR("Unexpected occurrence.");

    	// Delete edge
    	//
      	DeleteEdge(edgePtr);

    	return true;
    }

	    	/**#: [Description="Return iterator that references first
	    	       edge data item. (There is a const version of this
		       function as well)"] */
    EdgeDataIterator
    BeginEdgeData() { return data.begin(); }

    EdgeDataConstIterator
    BeginEdgeData() const { return data.begin(); }

	    	/**#: [Description="Return iterator that references the
	    	       position in memory after the last edge data item.
	    	       (There is a const version of this function as well)"] */
    EdgeDataIterator
    EndEdgeData() { return data.end(); }

    EdgeDataConstIterator
    EndEdgeData() const { return data.end(); }

	    	/**#: [Description="Return EdgeIterator that references first
	    	       edge in graph. (There is a const version of this function
		       as well)"] */
    EdgeIterator
    begin_edge() { return edges.begin(); }

    ConstEdgeIterator
    begin_edge() const { return edges.begin(); }

	    	/**#: [Description="Return EdgeIterator that references the
	    	       position in memory after the last edge in graph.(There
		       is a const version of this function as well)"] */
    EdgeIterator
    end_edge() { return edges.end(); }

    ConstEdgeIterator
    end_edge() const { return edges.end(); }


    	    	// Print details of all vertices and edges in graph and the
    	    	// connections betwen them.
    	    	//
	    	/**#: [Hidden] */
    virtual void
    Print(ostream& os) const
    {
    	Base::Print(os);

    	os << "EDGES:\n";
    	for (ConstEdgeIterator i=begin_edge(); i!=end_edge(); i++)
    	    os << *i;
    	os << '\n';
    }

}; // class BTL_GraphWithEdges

// Graph where vertex data and edge data are be stored in the same type
// of container which is neither a vector or a deque

	    	/**#: [Hidden] */

template <class Container>
class BTL_Graph3 : public
BTL_GraphWithEdges< BTL_ItD<Container>, BTL_ItD<Container> >
{};

// Graph where vertex data and edge data are be stored in the same type
// of container which is either a vector or a deque

	    	/**#: [Hidden] */

template <class Container>
class BTL_Graph4 : public
BTL_GraphWithEdges< BTL_InD<Container>, BTL_InD<Container> >
{};

// Graph where vertices and edges can be stored in different types of containers
// as long as NEITHER are STL vectors or the like.

	    	/**#: [Hidden] */

template <class Container1, class Container2>
class BTL_Graph5 : public
BTL_GraphWithEdges< BTL_ItD<Container1>, BTL_ItD<Container2> >
{};

// Graph where vertices and edges can be stored in different types of containers
// as long as BOTH are STL vectors or the like.

	    	/**#: [Hidden] */

template <class Container1, class Container2>
class BTL_Graph6 : public
BTL_GraphWithEdges< BTL_InD<Container1>, BTL_InD<Container2> >
{};

_BTL_END_NAMESPACE

#endif // GraphWithEdges.h

unsigned int BTL_EdgeBase::nextId = 1;
