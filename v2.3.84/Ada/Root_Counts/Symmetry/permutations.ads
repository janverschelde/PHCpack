with Standard_Integer_Vectors;

package Permutations is

-- DESCRIPTION :
--   This package defines the type Permutation.

  type Permutation is new Standard_Integer_Vectors.Vector;
  type Link_to_Permutation is new Standard_Integer_Vectors.Link_to_Vector;

  -- DESCRIPTION :
  --   A permutation p defines the image of i -> p(i).
  --   As also negative entries are alowed, sign permutations
  --   will be modelled as follows: 
  --     let perm = (1 3 -2), applied to F=(f1,f2,f3):
  --     perm*F = (f1,f3,-f2).

  function Is_Permutation ( p : Permutation ) return boolean;

  -- DESCRIPTION :
  --   Checks whether the vector p models a permutation:
  --    p(i) /= p(j) and p(i) /= -p(j), for all i /= j and
  --    -p'last <= p(i) <= p'last.

  function Equal ( p1,p2 : Permutation ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both permutations are equal.

  function "*" ( p1,p2 : Permutation ) return Permutation;

  -- DESCRIPTION :
  --   returns p1 `after' p2
  --   (1 3 -2) * (2 -1 3) = (3 -1 -2)
  -- REQUIRED :
  --   p1'range = p2'range

  function inv ( p : Permutation ) return Permutation;

  -- DESCRIPTION :
  --   inv(p)*p = p*inv(p) = the identical transformation

end Permutations;
