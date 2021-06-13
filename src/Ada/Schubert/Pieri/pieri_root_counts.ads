with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Generic_Lists;
with Pieri_Trees;                        use Pieri_Trees;

package Pieri_Root_Counts is

-- DESCRIPTION :
--   This package executes the combinatorial root counting method based on
--   the construction of Pieri trees and provides data structures to perform
--   the Pieri deformations.  These deformations start at the leaves that
--   satisfy Pieri's condition for which a triple intersection of planes
--   leads to a unique solution.  They continue following the chains that
--   end at those leaves.  These chains are indexed by the pairs of leaves.
--   So a list of paired nodes suffices to run the Pieri homotopy algorithm.
--   For efficiency one wants to avoid the repeated generation of the same
--   equations at every node.  Therefore a tree of paired nodes is needed.
--   In the creation of this tree, a chain is traversed in increasing order.

-- DATA STRUCTURES :

  type Paired_Nodes is record
    left,right : Link_to_Pieri_Node;
  end record;

  package Lists_of_Paired_Nodes is new Generic_Lists(Paired_Nodes);
  type List_of_Paired_Nodes is new Lists_of_Paired_Nodes.List;

  type Paired_Chain is array ( integer32 range <> ) of Paired_Nodes;

  type Nodal_Pair;
  type Link_to_Nodal_Pair is access Nodal_Pair;
  type Nodal_Matrix is 
     array (integer32 range <>,integer32 range <>) of Link_to_Nodal_Pair;

  type Nodal_Pair ( d : integer32 ) is record
    pnd : Paired_Nodes;
    sols : natural32;                     -- #paths arriving at the node
    children : Nodal_Matrix(0..d,0..d);   -- indices in children are jumps
    ancestor : Link_to_Nodal_Pair;
  end record;

-- CREATORS :

  function Create ( n,d : integer32; t1,t2 : Pieri_Tree )
                  return List_of_Paired_Nodes;

  -- DESCRIPTION :
  --   Creates a list of paired nodes that satisfy Pieri's condition.

  function Create ( pnd : Paired_Nodes ) return Paired_Chain;

  -- DESCRIPTION :
  --   Returns the chain of nodes that start at the first branch point down
  --   and end at the given pair of leaves pnd.
  --   This is an intermediate data structure in the Create operation below.

  function Create
             ( d : integer32; lp : List_of_Paired_Nodes ) return Nodal_Pair;

  -- DESCRIPTION :
  --   Creates a tree of nodal pairs and returns its root.
  --   The tree structure allows repeated nodes.
  --   A node that occured already once has sols equal to zero.

-- SELECTORS :

  function Height ( pnd : Paired_Nodes ) return natural32;

  -- DESCRIPTION :
  --   Returns Max(pnd.left.h,pnd.right.h).

  function Equal ( pnd1,pnd2 : Paired_Nodes ) return boolean;

  -- DESCRIPTION :
  --   Returns true when both nodes are the same.

  function At_First_Branch_Point ( pnd : Paired_Nodes ) return boolean;

  -- DESCRIPTION :
  --   Returns true when the current pair of nodes is at the first
  --   branch point in the Pieri trees.

  function At_Leaves ( pnd : Paired_Nodes ) return boolean;

  -- DESCRIPTION :
  --   Returns true when the current pair of nodes consists of leaves.

  function Ancestor ( pnd : Paired_Nodes ) return Paired_Nodes;

  -- DESCRIPTION :
  --   Returns the pair of ancestor nodes to the given pair of nodes.

  function First_Branch_Point ( pnd : Paired_Nodes ) return Paired_Nodes;

  -- DESCRIPTION :
  --   Returns the pair of nodes down in the chains of pnd where the
  --   first branch point occurs.

  function Height ( np : Nodal_Pair ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of nodes above the current nodal pair.

  function Is_In ( root : Nodal_Pair; pnd : Paired_Nodes ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the paired nodes occur somewhere in the tree of
  --   nodal pairs starting at the root.

  function Number_of_Paths ( root : Nodal_Pair ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of paths that need to be followed.

-- FORMATTED OUTPUT :

  procedure Write ( file : in file_type; chn : in Paired_Chain );

  -- DESCRIPTION :
  --   Writes the chain of paired nodes in a formatted way.

  procedure Write ( file : in file_type; root : in Nodal_Pair );

  -- DESCRIPTION :
  --   Writes the tree of paired nodes on file in a formatted way.

end Pieri_Root_Counts;
