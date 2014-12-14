with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;

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

  procedure Create ( h : in Link_to_Array_of_Strings; txk : in integer32 );

  -- DESCRIPTION :
  --   Stores the data for a natural parameter homotopy,
  --   where the index of the continuation parameter is the variable
  --   with index txk in the polynomials of h.

-- SELECTORS :

  function Standard_Homotopy_System
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Evaluates the homotopy stored as string representation by Create
  --   and returns the evaluated system in standard double precision.

  function DoblDobl_Homotopy_System
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Evaluates the homotopy stored as string representation by Create
  --   and returns the evaluated system in double double precision.

  function QuadDobl_Homotopy_System
             return QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Evaluates the homotopy stored as string representation by Create
  --   and returns the evaluated system in quad double precision.

  function Multprec_Homotopy_System ( deci : natural32 )
             return Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  -- DESCRIPTION :
  --   Evaluates the homotopy stored as string representation by Create
  --   and returns the evaluated system in arbitrary multiprecision,
  --   with as many decimal places as the value of deci.

  function Eval ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the homotopy at x and t in standard double precision.

  function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the homotopy at x and t in double double precision.

  function Eval ( x : QuadDobl_Complex_Vectors.Vector;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the homotopy at x and t in quad double precision.

  function Eval ( x : Multprec_Complex_Vectors.Vector;
                  t : Multprec_Complex_Numbers.Complex_Number;
                  d : natural32 )
                return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the homotopy at x and t in arbitrary multiprecision,
  --   with as many decimal places as d.

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the homotopy differentiated with respect to t
  --   and evaluated at x and t, evaluated in standard double precision.

  function Diff ( x : DoblDobl_Complex_Vectors.Vector;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the homotopy differentiated with respect to t
  --   and evaluated at x and t, evaluated in double double precision.

  function Diff ( x : QuadDobl_Complex_Vectors.Vector;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the homotopy differentiated with respect to t
  --   and evaluated at x and t, evaluated in quad double precision.

  function Diff ( x : Multprec_Complex_Vectors.Vector;
                  t : Multprec_Complex_Numbers.Complex_Number;
                  d : natural32 )
                return Multprec_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the homotopy differentiated with respect to t
  --   and evaluated at x and t, evaluated in arbitrary multiprecision,
  --   using as many decimal places as the value of d.

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the evaluated Jacobian matrix of the homotopy,
  --   differentiated with respect to x, and evaluated at x and t
  --   in standard double precision.

  function Diff ( x : DoblDobl_Complex_Vectors.Vector;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the evaluated Jacobian matrix of the homotopy,
  --   differentiated with respect to x, and evaluated at x and t
  --   in double double precision.

  function Diff ( x : QuadDobl_Complex_Vectors.Vector;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the evaluated Jacobian matrix of the homotopy,
  --   differentiated with respect to x, and evaluated at x and t
  --   in quad double precision.

  function Diff ( x : Multprec_Complex_Vectors.Vector;
                  t : Multprec_Complex_Numbers.Complex_Number;
                  d : natural32 )
                return Multprec_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the evaluated Jacobian matrix of the homotopy,
  --   differentiated with respect to x, and evaluated at x and t
  --   in arbitrary multiprecision using d decimal places.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   All homotopies are cleared.

end Varbprec_Homotopy;
