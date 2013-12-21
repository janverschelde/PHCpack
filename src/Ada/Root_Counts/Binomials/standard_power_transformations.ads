with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer64_Vectors;
with Standard_Integer64_Matrices;

package Standard_Power_Transformations is

-- DESCRIPTION :
--   This package provides functions to form the unimodular matrices 
--   to represent the power transformations.  Power transformations 
--   operate in the space of the powers of a Laurent polynomial system.
--   The routines use standard 32-bit or 64-bit integer arithmetic.

  function Identity_Matrix
              ( n : natural32 ) return Standard_Integer_Matrices.Matrix;
  function Identity_Matrix
              ( n : natural32 ) return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the n-by-n identity matrix.

  function Pivot ( v : Standard_Integer_Vectors.Vector ) return integer32;
  function Pivot ( v : Standard_Integer64_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns v'first - 1 if the vector is zero,
  --   otherwise the first nonzero element in v is returned.

  function Pivot ( v : Standard_Integer_Vectors.Vector;
                   i : integer32 ) return integer32;
  function Pivot ( v : Standard_Integer64_Vectors.Vector;
                   i : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the first nonzero v(k) for k in i..v'last.
  --   If there is no such k, then v'last + 1 is returned.

  function Rotate ( v : Standard_Integer_Vectors.Vector;
                    i : integer32 ) return Standard_Integer_Matrices.Matrix;
  function Rotate ( v : Standard_Integer64_Vectors.Vector;
                    i : integer32 ) return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a unimodular transformation to reduce v to d times
  --   the i-th standard unit vector, where d is gcd(v).
  --   This coordinate transformation rotates a vector into
  --   a direction parallel to a coordinate axis.

  -- REQUIRED : v(i) /= 0, ensured by i = Pivot(v).

  function Eliminate ( v : Standard_Integer_Vectors.Vector; i : integer32 )
                     return Standard_Integer_Matrices.Matrix;
  function Eliminate ( v : Standard_Integer64_Vectors.Vector; i : integer32 )
                     return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a unimodular coordinate transformation to eliminate
  --   the i-th coordinate using a nonzero vector v.
  --   Applying the transformation to exponents will result
  --   in the same degree for the i-th coordinate for all points
  --   that make the same inner product with the vector v.
  --   Because of this property, we may divide out the i-th variable
  --   from all points that span the face with inner normal v.

  -- REQUIRED : v(i) /= 0, ensured by i = Pivot(v).

end Standard_Power_Transformations;
