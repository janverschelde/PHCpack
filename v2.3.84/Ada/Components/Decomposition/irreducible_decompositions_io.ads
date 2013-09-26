with text_io;                           use text_io;
with Irreducible_Decompositions;        use Irreducible_Decompositions;

package Irreducible_Decompositions_io is

-- DESCRIPTION :
--   This package offers input-output facilities for irreducible
--   decompositions of solution sets of polynomial systems.

-- REQUIRED : 
--   The symbol table is initialized with the symbols for the variables
--   of the original system, and the zz-symbols of the slack variables.

  procedure get ( dc : in out Standard_Irreducible_Decomposition );
  procedure get ( file : in file_type;
                  dc : in out Standard_Irreducible_Decomposition );

  -- DESCRIPTION :
  --   Reads an irreducible decomposition from standard input or
  --   from file.  The input format matches the output format.

  procedure put ( dc : in Standard_Irreducible_Decomposition );
  procedure put ( file : in file_type;
                  dc : in Standard_Irreducible_Decomposition );

  -- DESCRIPTION :
  --   Writes an irreducible decomposition on standard output
  --   or on file.  The output format consists of the top dimension
  --   as first natural number, followed by the embedded systems
  --   with the generic points in the format of a solution list,
  --   for all dimensions, ranging from 0 to the top dimension.

end Irreducible_Decompositions_io;
