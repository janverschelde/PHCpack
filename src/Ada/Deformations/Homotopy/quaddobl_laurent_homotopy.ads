with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Laur_Systems;      use QuadDobl_Complex_Laur_Systems;

package QuadDobl_Laurent_Homotopy is

-- DESCRIPTION :
--   This package provides implementation of cheater's homotopy for
--   Laurent polynomial systems with coefficients in double double precision.

-- CONSTRUCTORS :

  procedure Create ( p,q : in Laur_Sys; k : in natural32; 
                     a : in Complex_Number );

  -- DESCRIPTION :
  --   The following artificial-parameter homotopy is build :
  --     H(x,t) = a * ((1 - t)^k) * q + (t^k) * p.

  procedure Create ( p,q : in Laur_Sys; k : in natural32; a : in Vector );

  -- DESCRIPTION :
  --   This operation is similar to the Create above, except for the fact
  --   that every equation can have a different random constant.

  procedure Create ( p,q : in Laur_Sys; k : in natural32; a,b : in Vector;
                     linear : in boolean );

  -- DESCRIPTION :
  --   If linear, then the following artificial-parameter homotopy is build :
  --     H(i)(x,t) = a(i) * ((1 - t)^k) * q(i) + b(i)(t^k) * p(i),
  --   for i in p'range.  
  --   Otherwise, if not linear, then
  --     H(i)(x,t)
  --       = (1-[t-t*(1-t)*a(i)])^k * q(i) + (t - t*(1-t)*b(i))^k * p(i),
  --   for i in p'range.

  procedure Create ( p : in Laur_Sys; k : in integer32 );

  -- DESCRIPTION :
  --   Given a polynomial system p with dimension n*(n+1) and
  --   k an index, then t = x_k as continuation parameter.

-- SELECTOR :

  function Homotopy_System return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the homotopy system in the unknowns (x,t).

-- SYMBOLIC ROUTINES :

  function Eval ( t : Complex_Number ) return Laur_Sys;

  -- DESCRIPTION :
  --   The homotopy is evaluated in t and a polynomial system is returned

  function Diff ( t : Complex_Number ) return Laur_Sys;

  -- DESCRIPTION :
  --   The homotopy is symbolically differentiated w.r.t. t.

-- NUMERIC ROUTINES :

  function Eval ( x : Vector; t : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   The homotopy is evaluated in x and t and a vector is returned.

  function Diff ( x : Vector; t : Complex_Number ) return Vector;

  -- DESCRIPTION :
  --   The homotopy is differentiated wr.t. t and is evaluated in (x,t).

  function Diff ( x : Vector; t : Complex_Number ) return matrix;

  -- DESCRIPTION :
  --   The homotopy is differentiated to x and the Jacobi matrix
  --   of H(x,t) is returned.

  function Diff ( x : Vector; t : Complex_Number; k : integer32 ) return Vector;

  -- DESCRIPTION :
  --   The returning vector contains all derivatives from the homotopy
  --   to the unknown x_k; note that t = x_n+1.

  function Diff ( x : Vector; t : Complex_Number; k : integer32 ) return matrix;

  -- DESCRIPTION :
  --   The Jacobi matrix of the homotopy is returned where the kth
  --   column has been deleted; note that Diff(x,t,n+1) = Diff(x,t).

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   The homotopy is cleared.

end QuadDobl_Laurent_Homotopy;
