with Standard_Integer_Numbers;          use standard_Integer_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;

package Standard_Affine_Solutions is

-- DESCRIPTION :
--   An affine solution is a solution of a polynomial system in affine space,
--   i.e.: we can write the solution as a sum of a basis point b and
--   a linear combination of a set v of linearly independent vectors.
--   This package provides these rewriting operations.

-- NOTICE : this is a remainder, see Standard_Intrinsic_Solutions.

  function Rewrite ( sol : Solution; n : integer32; b : Vector; v,w : VecVec )
                   return Solution;

  -- DESCRIPTION :
  --   The solution is represented by coefficients of the linear combination
  --   with respect to the offset vector b and directions v.  The same affine
  --   plane is given by b and w and this routine returns the new coefficients
  --   of the linear combination of the directions in w.

  function Rewrite ( sols : Solution_List; n : integer32;
                     b : Vector; v,w : VecVec ) return Solution_List;
  function Rewrite ( sols : Array_of_Solution_Lists; n : integer32;
                     b : Vector; v,w : Array_of_VecVecs )
                   return Array_of_Solution_Lists;

  -- DESCRIPTION :
  --   Applies the rewrite for every solution in the lists.

end Standard_Affine_Solutions;
