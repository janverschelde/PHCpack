with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with TripDobl_Complex_Numbers;           use TripDobl_Complex_Numbers;
with TripDobl_Complex_Vectors;           use TripDobl_Complex_Vectors;
with TripDobl_Complex_Matrices;          use TripDobl_Complex_Matrices;
with TripDobl_Complex_Poly_Systems;      use TripDobl_Complex_Poly_Systems;

package TripDobl_Homotopy is

-- DESCRIPTION :
--   Defines polynomial homotopies with triple double complex coefficients.

-- CONSTRUCTORS :

  procedure Create ( p,q : in Poly_Sys; k : in natural32; 
                     a : in Complex_Number );

  -- DESCRIPTION :
  --   The following artificial-parameter homotopy is build :
  --     H(x,t) = a * ((1 - t)^k) * q + (t^k) * p.

  procedure Create ( p,q : in Poly_Sys; k : in natural32; a : in Vector );

  -- DESCRIPTION :
  --   This operation is similar to the Create above, except for the fact
  --   that every equation can have a different random constant.

  procedure Create ( p,q : in Poly_Sys; k : in natural32; a,b : in Vector;
                     linear : in boolean );

  -- DESCRIPTION :
  --   If linear, then the following artificial-parameter homotopy is build :
  --     H(i)(x,t) = a(i) * ((1 - t)^k) * q(i) + b(i)(t^k) * p(i),
  --   for i in p'range.  
  --   Otherwise, if not linear, then
  --     H(i)(x,t)
  --       = (1-[t-t*(1-t)*a(i)])^k * q(i) + (t - t*(1-t)*b(i))^k * p(i),
  --   for i in p'range.

  procedure Create ( p : in Poly_Sys; k : in integer32 );

  -- DESCRIPTION :
  --   Given a polynomial system p with dimension n*(n+1) and
  --   k an index, then t = x_k as continuation parameter.

-- SELECTORS :

  function Relaxation_Power return natural32;

  -- DESCRIPTION :
  --   Returns the value of the relaxation power k given at creation
  --   of an artificial-parameter homotopy.  If the homotopy was not
  --   created, or if the homotopy is a natural-parameter homotopy,
  --   then zero is returned.

  function Accessibility_Constant return Complex_Number;

  -- DESCRIPTION :
  --   Returns the value of the accessibility (or gamma) constant a
  --   given at the creation of an artificial-parameter homotopy.
  --   If the homotopy is empty or of natural-parameter type,
  --   then zero is returned.

  function Dimension return integer32;

  -- DESCRIPTION :
  --   Returns the number of equations in the homotopy or zero
  --   if there is no homotopy defined.

  function Target_System return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the target system in an artificial-parameter homotopy.

  -- REQUIRED : Dimension > 0 and the homotopy is artificial parameter.

  function Start_System return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the start system in an artificial-parameter homotopy.

  -- REQUIRED : Dimension > 0 and the homotopy is artificial parameter.

  function Homotopy_System return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the homotopy system in the unknowns (x,t).

-- SYMBOLIC ROUTINES :

  function Eval ( t : Complex_Number ) return Poly_Sys;

  -- DESCRIPTION :
  --   The homotopy is evaluated in t and a polynomial system is returned

  function Diff ( t : Complex_Number ) return Poly_Sys;

  -- DESCRIPTION :
  --   The homotopy is symbolically differentiated w.r.t. t.

-- NUMERIC ROUTINES :

  function Eval ( x : Vector; t : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   The homotopy is evaluated in x and t and a vector is returned.

  function Diff ( x : Vector; t : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   The homotopy is differentiated wr.t. t and is evaluated in (x,t).

  function Diff ( x : Vector; t : Complex_Number ) return Matrix;

  -- DESCRIPTION :
  --   The homotopy is differentiated to x and the Jacobi matrix
  --   of H(x,t) is returned.

  function Diff ( x : Vector; t : Complex_Number; k : integer32 )
                return Vector;

  -- DESCRIPTION :
  --   The returning vector contains all derivatives from the homotopy
  --   to the unknown x_k; note that t = x_n+1.

  function Diff ( x : Vector; t : Complex_Number; k : integer32 )
                return Matrix;

  -- DESCRIPTION :
  --   The Jacobi matrix of the homotopy is returned where the kth
  --   column has been deleted; note that Diff(x,t,n+1) = Diff(x,t).

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   The homotopy is cleared.

end TripDobl_Homotopy;
