with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Series;           use QuadDobl_Complex_Series;

package QuadDobl_Complex_Series_Norms is

-- DESCRIPTION :
--   The norm of a series is the square root of the series
--   multiplied by its complex conjugate,
--   computed in quad double precision.

  function Conjugate ( s : Series ) return Series;

  -- DESCRIPTION :
  --   The complex conjugate of a series s has as coefficients
  --   the complex conjugates of the coefficients of s.

  function Conjugate ( s : Link_to_Series ) return Link_to_Series;

  -- DESCRIPTION :
  --   Returns s if s is null, otherwise returns the complex
  --   conjugate of all coefficients in s.

  function Norm ( s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the square root of s multiplied by Conjugate(s).

  procedure Normalize ( s : in out Series );

  -- DESCRIPTION :
  --   Divides the series by its norm.

  function Normalize ( s : Series ) return Series;

  -- DESCRIPTION :
  --   Returns the normalization of the series s.

  function Max_Norm ( s : Series ) return quad_double;

  -- DESCRIPTION :
  --   The max norm of a series is the absolute value
  --   of the largest coefficient.

  function Two_Norm ( s : Series ) return quad_double;

  -- DESCRIPTION :
  --   The two norm of a series is the 2-norm of Norm(s).

end QuadDobl_Complex_Series_Norms;
