with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Permutations,Symmetry_Group;        use Permutations,Symmetry_Group;

package Drivers_for_Orbits_of_Solutions is

-- DESCRIPTION :
--   This package provides two drivers for reporting on the
--   computation of the orbits of a given list of solutions.

  procedure Driver_for_Orbits_of_Solutions
                  ( file : in file_type; sols : in out Solution_List;
                    v : in List_of_Permutations; allperms,signsym : in boolean;
                    tol : in double_float );

  procedure Driver_for_Orbits_of_Solutions
                  ( file : in file_type; sols : in out Solution_List;
                    v : in List_of_Permutations; allperms,signsym : in boolean;
                    tol : in double_float; orbi : out Permutation );

  -- DESCRIPTION :
  --   Computes the orbits of the given list of solutions, creates a
  --   list with only the generating solutions and reports on file.

  -- ON ENTRY :
  --   file         to write the results on, must be opened for output;
  --   sols         a solution list;
  --   v            list of permutations;
  --   allperms     when true, then v is the full permutation group;
  --   signsym      when true, there is additional sign symmetry;
  --   tol          tolerance for comparing the solution vectors.

  -- ON RETURN :
  --   sols         generating list of solutions;
  --   orbi         permutation vector, indicating the orbits, if provided
  --                as output parameter.

end Drivers_for_Orbits_of_Solutions;
