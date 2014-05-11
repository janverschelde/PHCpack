with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Brackets;                          use Brackets;

package Localization_Posets is

-- DESCRIPTION :
--   This package provides an abstraction of a poset of localization patterns.
--   A localization pattern is characterized by top and bottom pivots,
--   represented by a pair of brackets.  In constructing the poset we can
--   decrement the top pivots or increment the bottom pivots, when moving
--   to lower-dimensional spaces that contain the special p-planes.
--   We denote the Grassmannian of all p-planes in n-space by G(p,m+p),
--   where n = m+p.

-- DATASTRUCTURES :

  type Node_Type is (top,bottom,mixed);

  type Node;
  type Link_to_Node is access Node;

  type Array_of_Nodes is array ( integer32 range <> ) of Link_to_Node;
  type Link_to_Array_of_Nodes is access Array_of_Nodes;
  type Array_of_Array_of_Nodes is 
    array ( integer32 range <> ) of Link_to_Array_of_Nodes;

  type Matrix_of_Nodes is
    array ( integer32 range <>, integer32 range <> ) of Link_to_Node;

  type Node ( p : integer32 ) is record
    tp : Node_Type;                            -- type of node
    level,label : integer32;                   -- coordinates in poset
    roco : integer32;                          -- root count, also as flag
    top,bottom : Bracket(1..p);                -- top(i) <= bottom(i)
    prev_sibling : Link_to_Node;               -- previous node at same level
    next_sibling : Link_to_Node;               -- next node at same level
    children : Matrix_of_Nodes(0..p,0..p);
    child_labels : Link_to_Vector;
  end record;
  -- The level of a node is a natural number between m*p+q*(m+p) and 0,
  -- children nodes have a level number that is one less.
  -- At a leaf in the poset, level = 0 and the p-plane is completely specified.
  -- At the root of the poset, level = m*p and any p-plane will fit it.
  -- For the quantum case, the bottom level is m*p+q*(m+p).
  -- In going up in the poset, the top pivots are incremented and
  -- the bottom pivots are decremented, keeping top(i) <= bottom(i).
  -- The field roco (root count) is determined by the root counting procedure.
  -- At a leaf, roco = 1.  At a node, roco = sum of roco's of its children.
  -- The labels of the children are determined in a indexed poset.

-- CREATORS :

  function Trivial_Root ( m,p : natural32 ) return Node;

  -- DESCRIPTION :
  --   Returns the root of the localization poset, with top and bottom pivots
  --   for the trivial localization pattern for any p-plane in (m+p)-space.
  --   The level of the trivial root equals m+p.

  function Trivial_Root ( m,p,q : natural32 ) return Node;

  -- DESCRIPTION :
  --   Returns the root for the poset to count all curves of degree q.
  --   The case q = 0 is the default, as implemented by the routine above.

  procedure Top_Create ( root : in Link_to_Node; n : in natural32 );

  -- DESCRIPTION :
  --   Creates a poset consistently incrementing top pivots,
  --   where n is the maximal entry in a top pivot.

  procedure Q_Top_Create ( root : in out Link_to_Node; n,lag : in natural32 );

  -- DESCRIPTION :
  --   Creates a poset consistently incrementing top pivots,
  --   where n is the maximal entry in a top pivot and lag the maximal
  --   space between first and last entry, typically, lag = m+p.
  --   For root = Trivial_Root(m,p,q) this poset will serve to count all
  --   maps of degree q in G(p,m+p) that meet m*p+q*(m+p) general m-planes.

  procedure Top_Create ( root : in Link_to_Node;
                         k : in Bracket; n : in natural32 );

  -- DESCRIPTION :
  --   Creates the poset by incrementing top pivots, with k = co-dimensions.
  --   By default, compared to the other Top_Create, all k's are one.
  --   So all nodes created are of type "top".

  procedure Q_Top_Create ( root : in Link_to_Node;
                           k : in Bracket; n,lag : in natural32 );

  -- DESCRIPTION :
  --   Quantum analogue for creation of poset incrementing top pivots
  --   for given co-dimension conditions.

  procedure Bottom_Create ( root : in Link_to_Node );

  -- DESCRIPTION :
  --   Creates a poset consistently decrementing bottom pivots.
  --   So all nodes created are of type "bottom".

  procedure Q_Bottom_Create ( root : in Link_to_Node; lag : in natural32 );

  -- DESCRIPTION :
  --   Creates a poset consistently decrementing bottom pivots for counting
  --   all maps of degree q with root = Trivial_Root(m,p,q).
  --   The parameter lag > max space between first and last entry in brackets.

  procedure Bottom_Create ( root : in Link_to_Node; k : in Bracket );

  -- DESCRIPTION :
  --   Creates the poset by decrementing bottom pivots, k = co-dimensions.
  --   By default, compared to the other Bottom_Create, all k's are one.

  procedure Q_Bottom_Create ( root : in Link_to_Node; k : in Bracket;
                              lag : in natural32 );

  -- DESCRIPTION :
  --   Quantum analogue for creation of poset decrementing bottom pivots
  --   for given co-dimension conditions.

  procedure Top_Bottom_Create ( root : in Link_to_Node; n : in natural32 );

  -- DESCRIPTION :
  --    Creates a poset by incrementing top and decrementing bottom pivots.

  procedure Q_Top_Bottom_Create ( root : in Link_to_Node;
                                  n,lag : in natural32 );

  -- DESCRIPTION :
  --    Creates a poset by incrementing top and decrementing bottom pivots
  --    to count all maps of degree q in G(p,m+p) that meet m*p+q*(m+p)
  --    general m-planes.

  procedure Top_Bottom_Create ( root : in Link_to_Node;
                                k : in Bracket; n : in natural32 );

  -- DESCRIPTION :
  --   Creates a poset by incrementing top and decrementing bottom pivots.
  --   The vector k contains co-dimensions.  All nodes are of type "mixed".

  procedure Q_Top_Bottom_Create ( root : in Link_to_Node;
                                  k : in Bracket; n,lag : in natural32 );

  -- DESCRIPTION :
  --   Quantum analogue for creation of poset in mixed fashion,
  --   for given co-dimension conditions.

  function Create_Leveled_Poset ( root : Link_to_Node ) return Array_of_Nodes;

  -- DESCRIPTION :
  --   The array on return contains points to the first node at each level.

  function Create_Indexed_Poset ( poset : Array_of_Nodes )
                                return Array_of_Array_of_Nodes;

  -- DESCRIPTION :
  --   Returns the poset structure where every node has a unique label.
  --   Every node contains a vector with labels to the children.
  --   The coordinates of every node in the poset on return is uniquely
  --   determined by level and label, as poset(node.level)(node.label).

-- SELECTORS :

  function Equal ( nd1,nd2 : Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the levels, top and bottom pivots are equal
  --   for both nodes.  False is returned otherwise.

  function Is_Leaf ( nd : Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the node has all its children empty.

  function Find_Node ( root : Link_to_Node; lvl : integer32 )
                     return Link_to_Node;

  -- DESCRIPTION :
  --   Finds the first node at the given level.

  function Number_of_Siblings ( nd : Link_to_Node ) return natural32;

  -- DESCRIPTION :
  --   If nd = null, then 0 is returned else the number of return
  --   equals the number of nonzero siblings plus one.

  function Number_of_Children ( nd : Node ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of children of the node.

-- ITERATORS :

  generic
    with procedure Report ( lvlnd : in Node; continue : out boolean );
  procedure Enumerate_Siblings ( nd : in Node );

  -- DESCRIPTION :
  --   Visits all siblings of the given node and calls Report on each of them.

  generic
    with procedure Report ( lnd : in Link_to_Node; continue : out boolean );
  procedure Enumerate_Grand_Children ( nd : in Node; k : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all grandchildren of the node nd, k generations deep.
  --   For k = 1, these are just the children of nd.

  generic
    with procedure Modify ( lvlnd : in out Node; continue : out boolean );
  procedure Modify_Siblings ( nd : in out Node );

  -- DESCRIPTION :
  --   Runs through all siblings of the node and calls Modify on each of them.

-- COMBINATORIAL ROOT COUNTING :

  procedure Count_Roots ( poset : in out Array_of_Nodes );

  -- DESCRIPTION :
  --   Counts the roots by assigning roco fields in the leveled poset.
  --   The total number of roots equals poset(poset'last).roco.

  function Row_Root_Count_Sum
              ( poset : Array_of_Nodes; i : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Sums up all the root counts at the i-th level in the poset.

  function Root_Count_Sum ( poset : Array_of_Nodes ) return natural32;

  -- DESCRIPTION :
  --   Returns the sum of all root counts over levels 1 and higher.
  --   This amount equals the number of paths that need to be traced.

-- DESTRUCTORS :

  procedure Clear ( nd : in out Node );

  -- DESCRIPTION :
  --   Deallocation of the siblings.

  procedure Clear ( lnd : in out Link_to_Node );

  -- DESCRIPTION :
  --   Releases the memory.

  procedure Clear ( arrnd : in out Array_of_Nodes );

  -- DESCRIPTION :
  --   Deallocation of the memory sibling after sibling.

  procedure Clear ( arrnd : in out Link_to_Array_of_Nodes );

  -- DESCRIPTION :
  --   Release of the pointers and deallocation of memory.

  procedure Clear ( arrnd : in out Array_of_Array_of_Nodes );

  -- DESCRIPTION :
  --   Deallocation of all nodes in all arrays.

  procedure Clear ( matnd : in out Matrix_of_Nodes );

  -- DESCRIPTION :
  --   Applies the nodes destructor to every element in the matrix.

end Localization_Posets;
