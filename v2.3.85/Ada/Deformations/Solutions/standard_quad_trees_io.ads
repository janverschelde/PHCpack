with text_io;                            use text_io;
with Standard_Quad_Trees;                use Standard_Quad_Trees;

package Standard_Quad_Trees_io is

-- DESCRIPTION :
--   This package provides output routines for the cardinalities
--   of quad trees, useful to see the distribution.

  procedure put ( root : Link_to_Quad_Node );
  procedure put ( file : file_type; root : Link_to_Quad_Node );

  -- DESCRIPTION :
  --   Writes the cardinalities of all the nodes in the quad tree.

end Standard_Quad_Trees_io;
