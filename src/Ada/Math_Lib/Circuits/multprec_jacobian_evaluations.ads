with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;

package Multprec_Jacobian_Evaluations is

-- DESCRIPTION :
--   This package offers evaluation of a polynomial system and its matrix
--   of all partial derivatives using the Speelpenning example,
--   in arbitrary multiprecision arithmetic.
--   The system is assumed to be sparse but all polynomials may share
--   the same support, or share some monomials.

  procedure EvalDiff
              ( b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : out Multprec_Complex_Matrices.Matrix );
  procedure EvalDiff
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : out Multprec_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the system defined by (f,b,c).
  --   The common factors of the exponent vectors are in f
  --   and the Speelpenning products are defined by the vectors in b.
  --   In this version, the matrix A is the two dimensional structure.

  -- ON ENTRY :
  --   f        common factors in the exponent vectors,
  --            if omitted, then there are no common factors;
  --   b        bit vectors defining the positions in the products;
  --   c        coefficients of the terms in the system;
  --   k        indices in the flattened representation;
  --   x        point where to evaluate and differentiate the system;
  --   y        allocated work space for the evaluated monomials.

  -- ON RETURN :
  --   z        evaluated system at x;
  --   y        work space;
  --   A        evaluated Jacobian matrix.

  procedure EvalDiff
              ( b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : in Multprec_Complex_VecVecs.VecVec );
  procedure EvalDiff
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Multprec_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Multprec_Complex_Vectors.Vector;
                z : out Multprec_Complex_Vectors.Vector;
                y : out Multprec_Complex_VecVecs.VecVec;
                A : in Multprec_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates and differentiates the system defined by (f,b,c).
  --   The common factors of the exponent vectors are in f
  --   and the Speelpenning products are defined by the vectors in b.
  --   In this version, the matrix is stored as a vector of columns.

  -- ON ENTRY :
  --   f        common factors in the exponent vectors,
  --            if omitted, then there are no common factors;
  --   b        bit vectors defining the positions in the products;
  --   c        coefficients of the terms in the system;
  --   k        indices in the flattened representation;
  --   x        point where to evaluate and differentiate the system;
  --   y        allocated work space for the evaluated monomials;
  --   A        allocated array of columns of the Jacobian matrix.

  -- ON RETURN :
  --   z        evaluated system at x;
  --   y        work space;
  --   A        evaluated Jacobian matrix.

end Multprec_Jacobian_Evaluations;
