with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Brackets;                          use Brackets;

package Pieri_Trees is

-- DESCRIPTION :
--   A Pieri tree is a combinatorial model for the Pieri homotopy algorithm.
--   The nodes contain brackets.  The root of a Pieri tree is the bracket
--   [1 2 .. d], which is the bracket with the lowest sum of entries.
--   Nodes situated at the same height in the tree all have the same sum
--   of their entries.  If this sum is S, then nodes at the next level
--   have sum S+1, nodes at the previous level have sum S-1.
--   Given a vector r(r(0),r(1),..,r(r'last)), the Pieri tree t(r) collects
--   all chains of length r(0)+r(1)+..+r(r'last) or equivalently of height
--   r(0)+r(1)+..+r(r'last)+1, that increase everywhere, except possibly
--   at levels r(0)+1, r(0)+r(1)+1, .., and r(0)+r(1)+..+r(r'last-1)+1.

-- DATA STRUCTURES :

  type Pieri_Node;
  type Link_to_Pieri_Node is access Pieri_Node;
  type Array_of_Pieri_Nodes is
    array ( integer32 range <> ) of Link_to_Pieri_Node;

  type Pieri_Node ( d : integer32 ) is record
    c,i,h : natural32; 
    node : Bracket(1..d);
    children : Array_of_Pieri_Nodes(1..d);
    ancestor : Link_to_Pieri_Node;
  end record;
  -- The "c" is the number of the plane that is folded in at the closest
  -- node down that is a jumping-decreasing node.
  -- A jumping-decreasing node is a node where decreasing may occur.
  -- The "i" indicates how close the node is to the next jumping-decreasing
  -- node down in the Pieri tree.  At such a node, i = 0.
  -- The "h" is the height of the node in the Pieri tree, at the root h = 0.
  -- When children(i) /= null, the ith index has increased by one.

  type Pieri_Tree ( d,a : integer32 ) is record
    branches : Vector(0..a);          -- # branching-with-decrease levels
    root : Link_to_Pieri_Node;
  end record;

  type Bracket_Array is array ( integer32 range <> ) of Link_to_Bracket;

-- CREATOR :
  
  function Create ( n,d : natural32; r : Vector ) return Pieri_Tree;

  -- DESCRIPTION :
  --   Returns a Pieri tree T(r(0),r(1),..,r(r'last)), where the nodes
  --   contain brackets with d entries chosen from {1,2,..,n}.
  --   The Pieri tree contains r(0)+r(1)+..+r(r'last)+1 levels of nodes.
  --   Decreasing may occur at levels r(0)+1,r(0)+r(1)+1,.., up to
  --   level r(0)+r(1)+..+r(r'last-1)+1.

-- SELECTORS :

  function Height ( t : Pieri_Tree ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of levels in the Pieri tree t.

  function Is_Leaf ( nd : Pieri_Node ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the node has no children.

  function Jump ( b1,b2 : Bracket ) return integer32;

  -- DESCRIPTION :
  --   Returns the largest index j such that b1(j) < b2(j),
  --   or zero when no such index exists.

  function Jump ( nd : Pieri_Node ) return integer32;

  -- DESCRIPTION :
  --   Returns the largest index of increase with the ancestor nodes,
  --   or zero when there is no ancestor.

  function Lower_Jump_Decrease ( nd : Pieri_Node ) return Bracket;

  -- DESCRIPTION :
  --   Returns a bracket at the node where decreasing may occur.
  --   If the current node is such a jumping-decreasing node,
  --   then the bracket at the node nd is returned, 
  --   otherwise, the bracket at the next jumping-decreasing node
  --   down the tree is returned.
  --   In case nd.c = 0, or nd has no ancestors, then also nd.node
  --   is returned.

  function Lowest_Jump_Decrease ( nd : Pieri_Node ) return Bracket;

  -- DESCRIPTION :
  --    Returns the bracket at the lowest node where decreasing may occur.

  function Upper_Jump_Decrease ( nd : Pieri_Node ) return Bracket;

  -- DESCRIPTION :
  --   Returns a bracket at the node where decreasing may occur,
  --   similarly as the function above, but now up in the tree,
  --   each time taking the last index in the brackets to jump.

  generic
    with procedure Visit_Node ( lnd : in Link_to_Pieri_Node;
                                continue : out boolean );
  procedure Enumerate_Nodes ( t : in Pieri_Tree; level : in natural32 );

  -- DESCRIPTION :
  --   Enumerates all nodes at the same level.  The root is at level 1.

  generic
    with procedure Visit_Chain
                     ( b : in Bracket_Array; continue : out boolean );
  procedure Enumerate_Chains ( t : in Pieri_Tree );

  -- DESCRIPTION :
  --   Enumerates all chains in the Pieri tree.

  generic
    with procedure Visit_Paired_Chain ( b1,b2 : in Bracket_Array;
                                        continue : out boolean );
  procedure Enumerate_Paired_Chains ( t1,t2 : in Pieri_Tree );

  -- DESCRIPTION :
  --   Enumerates all pairs of chains.

  function Pieri_Condition ( n : natural32; b1,b2 : Bracket ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the pair of brackets satisfies Pieri's condition.
  --   The number n equals m+p, and is the range of numbers to choose from.

  -- REQUIRED : b1'range = b2'range and b1(i) <= n >= b2(i).

-- DESTRUCTORS :

  procedure Clear ( nd : in out Link_to_Pieri_Node );

  -- DESCRIPTION :
  --   Deallocation of the memory space occupied by the node.

  procedure Clear ( t : in out Pieri_Tree );

  -- DESCRIPTION :
  --   Deallocation of the memory space occupied by the Pieri tree t.

end Pieri_Trees;
