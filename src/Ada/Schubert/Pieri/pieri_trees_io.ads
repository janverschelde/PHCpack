with text_io;                         use text_io;
with Standard_Natural_Numbers;        use Standard_Natural_Numbers;
with Pieri_Trees;                     use Pieri_Trees;

package Pieri_Trees_io is

-- DESCRIPTION :
--   This package provides output operations for Pieri trees.

  procedure put ( nd : in Pieri_Node );
  procedure put ( file : in file_type; nd : in Pieri_Node );

  -- DESCRIPTION :
  --   Writes the bracket and the triplet (c,i,h).

  procedure put ( lnd : in Link_to_Pieri_Node );
  procedure put ( file : in file_type; lnd : in Link_to_Pieri_Node );

  -- DESCRIPTION :
  --   Writes the Pieri node and the whole chain down in the tree.

  procedure put ( t : in Pieri_Tree );
  procedure put ( file : in file_type; t : in Pieri_Tree );

  -- DESCRIPTION :
  --   Writes all nodes of the Pieri tree, level by level.

  procedure put ( t : in Pieri_Tree; level : in natural32 );
  procedure put ( file : in file_type; t : in Pieri_Tree;
                  level : in natural32 );

  -- DESCRIPTION :
  --   Writes all nodes of the Pieri tree at the indicated level.

  procedure Write_Tree ( t : in Pieri_Tree );
  procedure Write_Tree ( file : in file_type; t : in Pieri_Tree );

  -- DESCRIPTION :
  --   Writes the Pieri tree in a formatted way.

end Pieri_Trees_io;
