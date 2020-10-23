with DecaDobl_Complex_Numbers;            use DecaDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Matrices;

package DecaDobl_Hessian_Updaters is

-- DESCRIPTION :
--   To update Hessian matrices for monomials where only 1, 2, 3, or 4
--   variables are present, call Speel1, Speel2, Speel3, or Speel4,
--   respectively.  For any number of variables, call SpeelN.
--   All computations occur in deca double precision.

  procedure Speel1 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of one variable raised to
  --   some power higher than 1.

  -- ON ENTRY :
  --   H        the current Hessian matrix,
  --            initialized with zero if called for the first time.
  --   c        coefficient of the term in the circuit;
  --   xps      exponents of all variables in the monomial;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   x        values for all variables;
  --   pwt      values of higher powers of x to evaluate the common factor.

  -- ON RETURN :
  --   H        updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last = fac'last = 1.

  procedure Speel2 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of two variables,
  --   where at least one variable is raised to a power higher than 1.

  -- ON ENTRY :
  --   H        the current Hessian matrix,
  --            initialized with zero if called for the first time.
  --   c        coefficient of the term in the circuit;
  --   xps      exponents of all variables in the monomial;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   x        values for all variables;
  --   pwt      values of higher powers of x to evaluate the common factor.

  -- ON RETURN :
  --   H        updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last = 2 >= fac'last >= 1.

  procedure Speel3 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of three variables,
  --   where at least one variable is raised to a power higher than 1.

  -- ON ENTRY :
  --   H        the current Hessian matrix,
  --            initialized with zero if called for the first time.
  --   c        coefficient of the term in the circuit;
  --   xps      exponents of all variables in the monomial;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   x        values for all variables;
  --   pwt      values of higher powers of x to evaluate the common factor.

  -- ON RETURN :
  --   H        updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last = 3 >= fac'last >= 1.

  procedure Speel4 ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of four variables,
  --   where at least one variable is raised to a power higher than 1.

  -- ON ENTRY :
  --   H        the current Hessian matrix,
  --            initialized with zero if called for the first time.
  --   c        coefficient of the term in the circuit;
  --   xps      exponents of all variables in the monomial;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   x        values for all variables;
  --   pwt      values of higher powers of x to evaluate the common factor.

  -- ON RETURN :
  --   H        updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last = 4 >= fac'last >= 1.

  procedure SpeelN ( H : in out DecaDobl_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in DecaDobl_Complex_Vectors.Vector;
                     fwd : in DecaDobl_Complex_Vectors.Link_to_Vector;
                     bck : in DecaDobl_Complex_Vectors.Link_to_Vector;
                     pwt : in DecaDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Updates the Hessian H for any number of variables,
  --   for higher powers and with computed forward and backward
  --   products of the values.

  -- ON ENTRY :
  --   H        the current Hessian matrix,
  --            initialized with zero if called for the first time.
  --   c        coefficient of the term in the circuit;
  --   xps      exponents of all variables in the monomial;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   x        values for all variables;
  --   fwd      forward products of the variables;
  --   bck      backward products of the variables;
  --   pwt      values of higher powers of x to evaluate the common factor.

  -- ON RETURN :
  --   H        updated Hessian matrix, only for upper triangular part.

  -- REQUIRED : idx'last >= 5 >= fac'last >= 1.

end DecaDobl_Hessian_Updaters;
