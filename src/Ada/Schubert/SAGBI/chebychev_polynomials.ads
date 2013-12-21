with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;

package Chebychev_Polynomials is

-- DESCRIPTION :
--   This is a simple implementation to generate, differentiate and
--   evaluate Chebychev polynomials.  The polynomials are represented
--   as vectors of range 0..d, with d their degree.

  function Create ( k : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Creates the kth Chebychev polynomial.

  function Eval ( k : natural32; x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Evaluates the kth Chebychev polynomial at x.

  -- REQUIRED : x lies in [-1,+1].

  function Eval ( p : Vector; x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Evaluates the polynomial at x.

  function Diff ( p : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the 1st derivative of p.

  function Diff ( p : Vector; k : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the kth derivative of the polynomial p.

  function Int ( p : Vector ) return Vector;

  -- DESCRIPTION :
  --   Returns the antiderivative of p.

  function Int ( p : Vector; k : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the kth antiderivative of p.

end Chebychev_Polynomials;
