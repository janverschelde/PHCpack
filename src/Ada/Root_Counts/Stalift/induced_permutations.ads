with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;

package Induced_Permutations is

-- DESCRIPTION :
--   This package provides tools to deal with supports, 
--   permuted by MixedVol and needed for processing with
--   the black box polyhedral homotopies and solver.

  function Remove_Artificial_Origin
             ( L : Lists_of_Floating_Vectors.List;
               b : double_float )
             return Lists_of_Floating_Vectors.List;
  function Remove_Artificial_Origin
             ( L : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               b : double_float )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Returns copies of the same points in the lists L,
  --   with the exception of those points which have
  --   the artificial origin with lifting value equal to b.

  procedure Remove_Artificial_Origin
              ( L : in out Lists_of_Floating_Vectors.List;
                b : in double_float );
  procedure Remove_Artificial_Origin
              ( L : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                b : in double_float );

  -- DESCRIPTION :
  --   Removes the artificial origin from the lists L.
  --   The artificial origin has lifting value equal to b.

  function Is_Subset
             ( lifted,original : Lists_of_Floating_Vectors.List )
             return boolean;

  -- DESCRIPTION :
  --   Returns true if every point in lifted without its lifting
  --   belongs to the original list.

  function Permutation
             ( s,ls : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               mix : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Given lifted supports in ls ordered along type of mixture in mix,
  --   and the original supports in s, the vector p on return has s'range
  --   and p(k) is the new position of the k-th list in s for the 
  --   corresponding lifted support in ls to be a subset.
  --   In case of multiple correspondences, the smallest enclosing list
  --   of the lifted support is selected.

  function Shift_Indices 
             ( p : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   If p'first = 1, then p is returned, otherwise p'first = 0
  --   is assumed and the vector on return has the same contents as p,
  --   but the indexing starts at one.

  function Relabel_for_Zero
             ( p : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   If p contains a zero, then the vector on return contains
  --   all entries of p increased by one, otherwise p is returned.
  --   This is a patch for the permutation returned by MixedVol,
  --   which starts its indexing apparently at zero.
  --   The Shift_Indices of above is also applied to the vector on return.
  --   The relabeling happens automatically in the Permute operations below.

  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out Standard_Complex_Laur_Systems.Laur_Sys );
  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out DoblDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out DoblDobl_Complex_Laur_Systems.Laur_Sys );
  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out QuadDobl_Complex_Poly_Systems.Poly_Sys );
  procedure Permute ( p : in Standard_Integer_Vectors.Vector;
                      f : in out QuadDobl_Complex_Laur_Systems.Laur_Sys );

  -- DESCRIPTION :
  --   Applies the permutation p to the system f.

end Induced_Permutations;
