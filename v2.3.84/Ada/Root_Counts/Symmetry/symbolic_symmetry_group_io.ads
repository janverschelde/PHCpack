with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Permutations,Symmetry_Group;        use Permutations,Symmetry_Group;

package Symbolic_Symmetry_Group_io is

-- DESCRIPTION :
--   This package contains routines for input and output of permutations
--   and groups of permutations, with the use of the symbol table.
--   This means that symbols rather than positions are read and written.

-- REQUIRED : The symbol table may not be empty!

  procedure get ( p : in out Permutation );
  procedure get ( file : in file_type; p : in out Permutation );

  -- DESCRIPTION :
  --   Reads a permutation from standard input or from file.
  --   When read from standard input, the get will ask the user to retry
  --   until a valid permutation is entered.
   
  procedure get ( L : in out List_of_Permutations; n,nb : in natural32 );
  procedure get ( file : in file_type;
		  L : in out List_of_Permutations; n,nb : in natural32 );

  -- DESCRIPTION :
  --   Reads a list of permutations from standard in put or from file.
  --   When read from file, nonvalid permutations are ignored.
  --   When read from standard input, the get askes the user to retry
  --   until a valid permutation is entered.

  -- ON ENTRY :
  --   n          the number of elements in the permutations;
  --   nb         the total number of permutations that must be read.
      
  procedure put ( p : in Permutation );
  procedure put ( file : in file_type; p : in Permutation );

  -- DESCRIPTION :
  --   Writes a permutation on standard output or on file.

  procedure put ( L : in List_of_Permutations );
  procedure put ( file : in file_type; L : in List_of_Permutations );

  -- DESCRIPTION :
  --   Writes a list of permutations on standard output or on file.

end Symbolic_Symmetry_Group_io;
