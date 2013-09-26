with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;

package Standard_Jacobian_Evaluations is

-- DESCRIPTION :
--   This package offers evaluation of a polynomial system and its matrix
--   of all partial derivatives using the Speelpenning example.
--   The system is assumed to be sparse but all polynomials may share
--   the same support, or share some monomials.

  function Integer_to_Natural
              ( v : Standard_Integer_VecVecs.VecVec )
              return Standard_Natural_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Converts the type of the vector of vectors from integer to natural.

  procedure Standard_Jacobian_Evaluation
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   The common factors of the exponent vectors are in f
  --   and the Speelpenning products are defined by the vectors in b.

  procedure Standard_Jacobian_Evaluation
              ( v : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given the flattened representation (v,c,k) of a system and x,
  --   returns in z the value of the polynomial system at x
  --   and in A the matrix of all partial derivatives evaluated at x. 

  procedure Standard_Jacobian_Evaluation
              ( v : in Standard_Integer_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given the flattened representation (v,c,k) of a system and x,
  --   with integer vectors as exponents, first converts to naturals,
  --   returns in z the value of the polynomial system at x
  --   and in A the matrix of all partial derivatives evaluated at x. 

end Standard_Jacobian_Evaluations;
