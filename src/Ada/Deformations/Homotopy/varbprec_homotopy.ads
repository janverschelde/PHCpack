with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;

package Varbprec_Homotopy is

-- DESCRIPTION :
--   A variable precision homotopy can be evaluated and differentiated
--   at several levels of precision.

-- CONSTRUCTORS :

  procedure Create ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                     a : in Standard_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Stores the data for the artificial parameter homotopy
  --   H(x,t) = a * ((1-t)^k) * q + (t^k) * p, so q is the start system
  --   and p is the target system. 
  --   In storing the string representations of the polynomial systems,
  --   deep copies are made.

-- SELECTORS :

  function Standard_Homotopy_System
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Evaluates the homotopy stored as string representation by Create
  --   and returns the evaluated system in standard double precision.

  function Eval ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the homotopy at x and t in standard double precision.

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the homotopy differentiated with respect to t
  --   and evaluated at x and t.

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the evaluated Jacobian matrix of the homotopy,
  --   differentiated with respect to x, and evaluated at x and t
  --   in standard double precision.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   The homotopy is cleared.

end Varbprec_Homotopy;
