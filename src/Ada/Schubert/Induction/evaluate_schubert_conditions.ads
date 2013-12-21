with Standard_Complex_Vectors;
with Standard_Complex_Matrices;

package Evaluate_Schubert_Conditions is

-- DESCRIPTION :
--   This package provides functions to evaluate the Schubert
--   intersection conditions.

  function Eval ( X,P,B : Standard_Complex_Matrices.Matrix;
                  h,m : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the result of the evaluation of 
  --
  --         (  [ X | P ]*B*m )
  --     y = (                )
  --         (  < h , m > - 1 )
  -- 
  --   to express that Rank([ X | P ]) = r.
  --
  -- REQUIRED :
  --   X'range(1) = P'range(1) = 1..n, n = dimension of ambient space;
  --   if X'range(2) = 1..k, P'range(2) = 1..l, then B'range(1) = 1..l;
  --   B'range(2) = h'range = m'range = 1..r+1.
  --
  -- ON ENTRY : 
  --   X        is a k-plane in n-space, represented by an n-by-k matrix,
  --            with in its columns the generators for the k-plane X;
  --   P        is an l-plane in n-space, represented by an n-by-l matrix,
  --            with in its columns the generators for the l-plane P;
  --   B        is a random (k+l)-by-(r+1)-matrix;
  --   m        is a vector of r+1 multipliers for the colums of [X|P]*B;
  --   h        is a random (r+1)-vector to scale the multipliers with.
  --
  -- ON RETURN :
  --   if y is zero, then X meets P in the proper way.

end Evaluate_Schubert_Conditions;
