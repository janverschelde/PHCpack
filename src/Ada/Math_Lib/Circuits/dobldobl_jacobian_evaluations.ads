with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;

package DoblDobl_Jacobian_Evaluations is

-- DESCRIPTION :
--   This package offers evaluation of a polynomial system and its matrix
--   of all partial derivatives using the Speelpenning example.
--   The system is assumed to be sparse but all polynomials may share
--   the same support, or share some monomials.

  procedure DoblDobl_Jacobian_Evaluation
              ( v : in Standard_Natural_VecVecs.VecVec;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given the flattened representation (v,c,k) of a system and x,
  --   returns in z the value of the polynomial system at x
  --   and in A the matrix of all partial derivatives evaluated at x. 

  procedure DoblDobl_Jacobian_Evaluation
              ( v : in Standard_Integer_VecVecs.VecVec;
                c : in DoblDobl_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in DoblDobl_Complex_Vectors.Vector;
                z : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Given the flattened representation (v,c,k) of a system and x,
  --   returns in z the value of the polynomial system at x
  --   and in A the matrix of all partial derivatives evaluated at x. 

end DoblDobl_Jacobian_Evaluations;
