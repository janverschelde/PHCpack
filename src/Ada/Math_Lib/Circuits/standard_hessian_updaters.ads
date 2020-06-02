with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;

package Standard_Hessian_Updaters is

-- DESCRIPTION :
--   To update Hessian matrices for monomials where only 1, 2, 3, or 4
--   variables are present, call Speel1, Speel2, Speel3, or Speel4,
--   respectively.  For any number of variables, call SpeelN.

  procedure Speel1 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of one variable raised to
  --   some power higher than one.

  -- REQUIRED : idx'last = fac'last = 1.

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

  procedure Speel1 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Mirrors the other Speel1, but then for coefficient circuits,
  --   to deal with the special case of one variable raised to some
  --   power higher than one.

  -- REQUIRED : idx'last = fac'last = 1.

  -- ON ENTRY :
  --   hrp      real parts of the Hessian matrix;
  --   hip      imaginary parts of the Hessian matrix;
  --   rcff     real part of the coefficient of a term in the circuit;
  --   icff     imaginary part of the coefficient of a term in the circuit;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   xr       real parts of the values for all variables;
  --   xi       imaginary parts of the values for all variables;
  --   rpwt     real parts of the values of higher powers of x 
  --            to evaluate the common factor;
  --   ipwt     imaginary parts of the values of higher powers of x 
  --            to evaluate the common factor.

  -- ON RETURN :
  --   hrp      updated real parts of the Hessian matrix,
  --            only for upper triangular part;
  --   hip      updated imaginary parts of the Hessian matrix,
  --            only for upper triangular part.

  procedure Speel2 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of two variables,
  --   where at least one variable is raised to a power higher than one.

  -- REQUIRED : idx'last = 2 >= fac'last >= 1.

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

  procedure Speel2 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Mirrors the other Speel1, but then for coefficient circuits,
  --   to deal with the special case of two variables, where at
  --   least one variable is raised to a power higher than one.

  -- REQUIRED : idx'last = 2 >= fac'last >= 1.

  -- ON ENTRY :
  --   hrp      real parts of the Hessian matrix;
  --   hip      imaginary parts of the Hessian matrix;
  --   rcff     real part of the coefficient of a term in the circuit;
  --   icff     imaginary part of the coefficient of a term in the circuit;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   xr       real parts of the values for all variables;
  --   xi       imaginary parts of the values for all variables;
  --   rpwt     real parts of the values of higher powers of x 
  --            to evaluate the common factor;
  --   ipwt     imaginary parts of the values of higher powers of x 
  --            to evaluate the common factor.

  -- ON RETURN :
  --   hrp      updated real parts of the Hessian matrix,
  --            only for upper triangular part;
  --   hip      updated imaginary parts of the Hessian matrix,
  --            only for upper triangular part.

  procedure Speel3 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of three variables,
  --   where at least one variable is raised to a power higher than 1.

  -- REQUIRED : idx'last = 3 >= fac'last >= 1.

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

  procedure Speel3 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Mirrors the other Speel1, but then for coefficient circuits,
  --   to deal with the special case of two variables, where at
  --   least one variable is raised to a power higher than one.

  -- REQUIRED : idx'last = 3 >= fac'last >= 1.

  -- ON ENTRY :
  --   hrp      real parts of the Hessian matrix;
  --   hip      imaginary parts of the Hessian matrix;
  --   rcff     real part of the coefficient of a term in the circuit;
  --   icff     imaginary part of the coefficient of a term in the circuit;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   xr       real parts of the values for all variables;
  --   xi       imaginary parts of the values for all variables;
  --   rpwt     real parts of the values of higher powers of x 
  --            to evaluate the common factor;
  --   ipwt     imaginary parts of the values of higher powers of x 
  --            to evaluate the common factor.

  -- ON RETURN :
  --   hrp      updated real parts of the Hessian matrix,
  --            only for upper triangular part;
  --   hip      updated imaginary parts of the Hessian matrix,
  --            only for upper triangular part.

  procedure Speel4 ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Deals with the special case of four variables,
  --   where at least one variable is raised to a power higher than one.

  -- REQUIRED : idx'last = 4 >= fac'last >= 1.

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

  procedure Speel4 ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Mirrors the previous Speel4, but then for coefficient circuits,
  --   to deal with the special case of four variables,
  --   where at least one variable is raised to a power higher than one.

  -- REQUIRED : idx'last = 4 >= fac'last >= 1.

  -- ON ENTRY :
  --   hrp      real parts of the Hessian matrix;
  --   hip      imaginary parts of the Hessian matrix;
  --   rcff     real part of the coefficient of a term in the circuit;
  --   icff     imaginary part of the coefficient of a term in the circuit;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   xr       real parts of the values for all variables;
  --   xi       imaginary parts of the values for all variables;
  --   rpwt     real parts of the values of higher powers of x 
  --            to evaluate the common factor;
  --   ipwt     imaginary parts of the values of higher powers of x 
  --            to evaluate the common factor.

  -- ON RETURN :
  --   hrp      updated real parts of the Hessian matrix,
  --            only for upper triangular part;
  --   hip      updated imaginary parts of the Hessian matrix,
  --            only for upper triangular part.

  procedure SpeelN ( H : in out Standard_Complex_Matrices.Matrix;
                     c : in Complex_Number;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     x : in Standard_Complex_Vectors.Vector;
                     fwd : in Standard_Complex_Vectors.Link_to_Vector;
                     bck : in Standard_Complex_Vectors.Link_to_Vector;
                     pwt : in Standard_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Updates the Hessian H for any number of variables,
  --   for higher powers and with computed forward and backward
  --   products of the values.

  -- REQUIRED : idx'last >= 5 >= fac'last >= 1.

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

  procedure SpeelN ( hrp : in Standard_Floating_VecVecs.VecVec;
                     hip : in Standard_Floating_VecVecs.VecVec;
                     rcff,icff : in double_float;
                     xps : in Standard_Integer_Vectors.Vector;
                     idx : in Standard_Integer_Vectors.Vector;
                     fac : in Standard_Integer_Vectors.Vector;
                     xr : in Standard_Floating_Vectors.Vector;
                     xi : in Standard_Floating_Vectors.Vector;
                     rfwd : in Standard_Floating_Vectors.Link_to_Vector;
                     ifwd : in Standard_Floating_Vectors.Link_to_Vector;
                     rbck : in Standard_Floating_Vectors.Link_to_Vector;
                     ibck : in Standard_Floating_Vectors.Link_to_Vector;
                     rpwt : in Standard_Floating_VecVecs.VecVec;
                     ipwt : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Mirrors the previous SpeelN, but then for coefficient circuits.
  --   Updates the Hessian H for any number of variables,
  --   for higher powers and with computed forward and backward
  --   products of the values.

  -- REQUIRED : idx'last >= 5 >= fac'last >= 1.

  -- ON ENTRY :
  --   hrp      real parts of the Hessian matrix;
  --   hip      imaginary parts of the Hessian matrix;
  --   rcff     real part of the coefficient of a term in the circuit;
  --   icff     imaginary part of the coefficient of a term in the circuit;
  --   idx      index of the participating variables;
  --   fac      indices to the variables in the common factor;
  --   xr       real parts of the values for all variables;
  --   xi       imaginary parts of the values for all variables;
  --   rfwd     real parts of the forward products of the variables;
  --   ifwd     imaginary parts of the forward products of the variables;
  --   rbck     real parts of the backward products of the variables;
  --   ibck     imaginary parts of the backward products of the variables;
  --   rpwt     real parts of the values of higher powers of x 
  --            to evaluate the common factor;
  --   ipwt     imaginary parts of the values of higher powers of x 
  --            to evaluate the common factor.

  -- ON RETURN :
  --   hrp      updated real parts of the Hessian matrix,
  --            only for upper triangular part;
  --   hip      updated imaginary parts of the Hessian matrix,
  --            only for upper triangular part.

end Standard_Hessian_Updaters;
