with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Permutations;                       use Permutations;

package Shuffle_Polynomials is

-- DESCRIPTION :
--   This package provides some operations to rearrange the order of
--   equations in a polynomial system to investigate the influence of
--   the order on the performance of the equation-by-equation solver.

  function Increasing_Degree_Sort ( f : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   If q = Increasing_Degree_Sort(f), Degree(q(i)) > Degree(q(j)), i<j.

  function Decreasing_Degree_Sort ( f : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   If q = Decreasing_Degree_Sort(f), Degree(q(i)) > Degree(q(j)), i<j.

  function Permute ( f : Poly_Sys; p : Permutation ) return Poly_Sys;

  -- DESCRIPTION :
  --   If q = Permute(f,p), then we have q(i) = f(p(i)).

  function Random_Permutation ( n : integer32 ) return Permutation;

  -- DESCRIPTION :
  --   Returns a random permutation to permute n polynomials.

  function Random_Permute ( f : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Applies a random permutation to the polynomials in f.

end Shuffle_Polynomials;
