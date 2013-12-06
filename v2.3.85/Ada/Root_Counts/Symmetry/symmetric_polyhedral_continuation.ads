with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Symmetry_Group;                     use Symmetry_Group;

package Symmetric_Polyhedral_Continuation is

-- DESCRIPTION :
--   Polyhedral continuation based on symmetric mixed subdivision.

  function Symmetric_Mixed_Solve
               ( file : file_type; grp : List_of_Permutations; sign : boolean;
    	         p : Laur_Sys; mixsub : Mixed_Subdivision;
                 n : integer32; mix : Vector ) return Solution_List;

  -- DESCRIPTION :
  --   This function computes the generating solutions of a given
  --   Laurent polynomial system, by making use of its mixed subdivision.

  -- ON ENTRY :
  --   file      a file to write intermediate results on;
  --   grp       representations of the symmetry group;
  --   sign      if true, then there is sign symmetry;
  --   p         a lifted Laurent polynomial system;
  --   mixsub    the mixed subdivision of the supports of p;
  --   n         the number of polynomials in p;
  --   mix(k)    indicates the number of occurencies of the kth support.

  -- REQUIRED :
  --   The polynomials in p should be ordered according to the
  --   information in the vector `mixed_type'!

  function Symmetric_Mixed_Solve
               ( file : file_type; sign : boolean; p : Laur_Sys; 
                 mixsub : Mixed_Subdivision;
                 n : integer32; mix : Vector ) return Solution_List;

  -- DESCRIPTION :
  --   Here the general permutation group is assumed.

end Symmetric_Polyhedral_Continuation;
