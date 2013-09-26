with Generic_Lists;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;

package Trees_of_Vectors is

-- DESCRIPTION :
--   This package provides a data abstraction for working
--   with trees of vectors.

-- THE DATA STRUCTURES :

  type Tree_of_Vectors;
  type Link_to_Tree_of_Vectors is access Tree_of_Vectors;

  type node is record
    d : Link_to_Vector;
    ltv : Link_to_Tree_of_Vectors;
  end record;

  package Link_to_Vector_Trees is new Generic_Lists(node);
  type Tree_of_Vectors is new Link_to_Vector_Trees.List;

-- SELECTORS :

  function Is_In ( tv : Tree_of_Vectors; v : Vector ) return boolean;
  function Is_In ( tv : Tree_of_Vectors; v : Link_to_Vector ) return boolean;

  -- DESCRIPTION :
  --   returns true if v belongs to the top level of the tree.

  generic

    with procedure Process ( nd : in node; continue : out boolean );

  procedure Iterator ( tv : in Tree_of_Vectors );

  -- DESCRIPTION :
  --   A walk through the tree will be made, left branch first.

-- DESTRUCTORS :

  procedure Clear ( nd : in out node );
  procedure Clear ( tv : in out Tree_of_Vectors );
  procedure Clear ( ltv : in out Link_to_Tree_of_Vectors );

  -- DESCRIPTION :
  --   All allocated memory space will be freed.

end Trees_of_Vectors;
