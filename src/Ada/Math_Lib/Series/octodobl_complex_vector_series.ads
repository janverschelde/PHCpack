with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with OctoDobl_Complex_Numbers;          use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_VecVecs;
with OctoDobl_Complex_Series_Vectors;

package OctoDobl_Complex_Vector_Series is

-- DESCRIPTION :
--   A series vector is a vector of truncated power series, with complex
--   numbers as coefficients.  A vector series is a truncated power series,
--   where the coefficients are vectors.  As a data type, a vector series
--   is represented as a vector of vectors, of a certain degree.
--   Coefficients are in octo double precision.

  type Vector ( deg : integer32 ) is record
     -- the last power in the series -- the error is of order deg+1
    cff : OctoDobl_Complex_VecVecs.VecVec(0..deg);
  end record;
  type Link_to_Vector is access Vector;

-- CONSTRUCTORS :

  function Create ( v : OctoDobl_Complex_Series_Vectors.Vector )
                  return OctoDobl_Complex_Vector_Series.Vector;

  -- DESCRIPTION :
  --   A vector which has as entries truncated power series
  --   is converted into a series with vectors as coefficients.
  --   The range of the vector on entry is 1..n, where n is the
  --   number of series in as entries in v.
  --   The range of the vector on return is 1..d, where d is the
  --   degree of each series in v.

  -- REQUIRED :
  --   The vector v is not empty.
  --   The assumption is that every series in v has the same degree d.

  function Create ( v : OctoDobl_Complex_Vector_Series.Vector )
                  return OctoDobl_Complex_Series_Vectors.Vector;

  -- DESCRIPTION :
  --   A truncated power series with vectors as coefficients
  --   is converted into a vector of truncated power series.
  --   The vector on entry has range 1..d, where d is v.deg.
  --   The vector on return has range 1..n, where n = v(i)'last.

  -- REQUIRED :
  --   The degree of v must be at least zero, at least one of
  --   the coefficients in v must be defined.

-- EVALUATORS :

  function Eval ( v : OctoDobl_Complex_Vector_Series.Vector;
                  t : octo_double )
                return OctoDobl_Complex_Vectors.Vector;
  function Eval ( v : OctoDobl_Complex_Vector_Series.Vector;
                  t : Complex_Number )
                return OctoDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the values of all series in v at t.

  -- REQUIRED : v.deg >= 0 and v.cff(0) is defined.

-- DESTRUCTORS :

  procedure Clear ( v : in out OctoDobl_Complex_Vector_Series.Vector );
  procedure Clear ( v : in out OctoDobl_Complex_Vector_Series.Link_to_Vector );

  -- DESCRIPTION :
  --   Deallocates all coefficients in the series.

end OctoDobl_Complex_Vector_Series;
