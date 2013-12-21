with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer64_Vectors;         use Standard_Integer64_Vectors;
with Standard_Integer64_Simplices;       use Standard_Integer64_Simplices;
with Standard_Integer64_Triangulations;  use Standard_Integer64_Triangulations;

package Dynamic64_Lifting_Functions is

-- DESCRIPTION :
--   This package contains operations to control the lifting function
--   in order to obtain triangulations either by placing or by pulling,
--   using 64-bit integers.

  function Lift_to_Place ( s : Simplex; x : Vector ) return integer64;
  function Lift_to_Place ( t : Triangulation; x : Vector ) return integer64;

  -- DESCRIPTION :
  --   Determines the next value of the lifting function for the given
  --   point x w.r.t. each simplex s, such that the product of the inner
  --   normal of s with x is higher than the product with any point in s.

  -- REQUIRED :
  --   x(x'last) provides a lower bound for the lifting function.
  --   To obtain an optimal lifting, x(x'last) can be set equal to 1.

  function Lift_to_Pull ( s : Simplex; x : Vector ) return integer64;
  function Lift_to_Pull ( t : Triangulation; x : Vector ) return integer64;

  -- DESCRIPTION :
  --   Determines the lifting value of x such that the product of the inner
  --   normal of each simplex s with x is strictly lower than the product
  --   with any point of s.
  --   The lifting value such that x can be pulled in s will be returned.

  -- REQUIRED :
  --   x(x'last) contains already a lower bound on the lifting for x.
  --   To obtain an optimal lifting, x(x'last) can be set equal to -1.

  function Degenerate ( t : Triangulation; x : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there exists a cell in t, such that the lifted point
  --   x lies in the same hyperplane as the corresponding simplex on the
  --   lower hull.

  function Lift_to_Pull
             ( t1,t2 : Triangulation; x : Vector ) return integer64;

  -- DESCRIPTION :
  --   Determines the lifting value of x such that for all cells in t1,
  --   the inner product of the inner normal with x is strictly lower
  --   than the inner product with any point of that cell.
  --   The triangulation t2 is used to be guarantee that the lifting does
  --   not induce a degeneracy in the lower hull.
  --   The lifting value such that x can be pulled in t will be returned.

  -- REQUIRED :
  --   x(x'last) contains already a lower bound on the lifting for x.
  --   To obtain an optimal lifting, x(x'last) can be set equal to -1.

end Dynamic64_Lifting_Functions;
