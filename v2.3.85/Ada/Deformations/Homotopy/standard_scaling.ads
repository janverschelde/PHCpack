with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_Scaling is

-- DESCRIPTION :
--   This package provides routines for scaling polynomial systems.

  procedure Scale ( p : in out Poly );
  procedure Scale ( s : in out Poly_Sys );

  -- DESCRIPTION :
  --   In each polynomial the coefficients are divided by the average
  --   coefficient of that polynomial.

  procedure Scale ( s : in out Poly_Sys; bas : in natural32 := 10;
                    diff : in boolean; cond : out double_float;
                    sccff : out vector );

  -- DESCRIPTION :
  --   Equation and variable scaling of a polynomial system.

  -- ON ENTRY :
  --   s        a polynomial system;
  --   bas      must be 2 or 10;
  --   diff     true to reduce the difference between the coefficients,
  --            false, otherwise.
 
  -- ON RETURN :
  --   s        is the scaled polynomial system by centering the coefficients
  --            close to units and eventually by reducing the difference 
  --            between the coefficients;
  --   cond     is an estimate for condition number of the linear system
  --            that had to be solved for scaling the polynomial system;
  --            it is an indication how good or bad the polynomial system 
  --            was scaled;
  --   sccff    is a vector containing the the exponents (w.r.t. basis 10) of
  --            the factors to scale the solutions back to the original
  --            coordinates.

  procedure Scale ( basis : in natural32; sccff : in Vector;
                    s : in out Solution );
  procedure Scale ( basis : in natural32; sccff : in Vector;
                    sols : in out Solution_List );

  -- DESCRIPTION :
  --   The solution(s) is (are) scaled backwards to the original coordinates,
  --   the vector sccff has been constructed by the procedure above.

end Standard_Scaling;
