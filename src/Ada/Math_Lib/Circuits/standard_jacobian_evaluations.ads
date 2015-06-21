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

  procedure EvalDiff
              ( b : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix );
  procedure EvalDiff
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix );

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
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_VecVecs.VecVec;
                A : in Standard_Complex_VecVecs.VecVec );
  procedure EvalDiff
              ( f,b : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_VecVecs.VecVec;
                A : in Standard_Complex_VecVecs.VecVec );

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

-- WRAPPERS :

  procedure Standard_Jacobian_Evaluation
              ( v : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix );
  procedure Standard_Jacobian_Evaluation
              ( v : in Standard_Natural_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_VecVecs.VecVec;
                A : in Standard_Complex_VecVecs.VecVec );

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
                y : out Standard_Complex_VecVecs.VecVec;
                A : out Standard_Complex_Matrices.Matrix );
  procedure Standard_Jacobian_Evaluation
              ( v : in Standard_Integer_VecVecs.VecVec;
                c : in Standard_Complex_VecVecs.VecVec;
                k : in Standard_Natural_VecVecs.VecVec;
                x : in Standard_Complex_Vectors.Vector;
                z : out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_VecVecs.VecVec;
                A : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Given the flattened representation (v,c,k) of a system and x,
  --   with integer vectors as exponents, first converts to naturals,
  --   returns in z the value of the polynomial system at x
  --   and in A the matrix of all partial derivatives evaluated at x. 

end Standard_Jacobian_Evaluations;
