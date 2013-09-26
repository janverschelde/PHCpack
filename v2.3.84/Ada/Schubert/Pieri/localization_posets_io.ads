with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Brackets;                          use Brackets;
with Localization_Posets;               use Localization_Posets;  

package Localization_Posets_io is

-- DESCRIPTION :
--   Provides output routines for localization posets.

  procedure put ( top,bottom : in Bracket; roco : in natural32 );
  procedure put ( file : in file_type;
                  top,bottom : in Bracket; roco : in natural32 );

  -- DESCRIPTION :
  --   Writes the top and bottom pivots as a couple of brackets,
  --   added with the root count, without taking a new line.

  procedure put ( root : in Node; lvl : in integer32 );
  procedure put ( file : in file_type; root : in Node; lvl : in integer32 );

  -- DESCRIPTION :
  --   Writes all nodes in the poset at the given level.

  procedure put ( poset : in Node );
  procedure put ( file : in file_type; poset : in Node );

  -- DESCRIPTION :
  --   Writes all nodes in the poset, level by level.

  procedure put ( poset : in Array_of_Nodes );
  procedure put ( file : in file_type; poset : in Array_of_Nodes );

  -- DESCRIPTION :
  --   Writes the leveled poset on standard output or on file.

  procedure put ( poset : in Array_of_Array_of_Nodes );
  procedure put ( file : in file_type; poset : in Array_of_Array_of_Nodes );

  -- DESCRIPTION :
  --   Writes the indexed poset on standard output or on file.
  --   The labels of the children to each node are in curled braces.

  procedure put_roco ( poset : in Array_of_Array_of_Nodes );
  procedure put_roco
              ( file : in file_type; poset : in Array_of_Array_of_Nodes );

  -- DESCRIPTION :
  --    Writes at each level the root count only.

end Localization_Posets_io;
