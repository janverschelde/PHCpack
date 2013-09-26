with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Permutations,Symmetry_Group;        use Permutations,Symmetry_Group;

package Symmetry_Group_io is

-- DESCRIPTION :
--   This package contains routines for input and output of 
--   permutations and groups of permutations.

  procedure get ( p : out Permutation );
  procedure get ( file : in file_type; p : out Permutation );

  -- DESCRIPTION :
  --   Reads a permutation from standard input or from file.
   
  procedure get ( L : in out List_of_Permutations; n,nb : in integer32 );
  procedure get ( file : in file_type;
		  L : in out List_of_Permutations; n,nb : in integer32 );

  -- DESCRIPTION :
  --   Reads a list of permutations from standard in put or from file.
  --   Vectors that do not represent permutations are ignored.

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

end Symmetry_Group_io;
