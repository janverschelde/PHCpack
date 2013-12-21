with text_io;                           use text_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;

package Standard_Flag_Representations is

-- DESCRIPTION :
--   A flag is a sequence of linear spaces, nested so that any lower
--   dimensional subspace belongs to any larger linear space.
--   The extrinsic representation of a flag is a sequence of hyperplanes,
--   whose coefficients are stored in the rows of a matrix.
--   The intrinsic representation of a flag is a sequence of generators.

  function Create_Intrinsic ( f : Matrix ) return Matrix;

  -- DESCRIPTION :
  --   Given an extrinsic description of a flag,
  --   the matrix on return has in its column the corresponding
  --   intrinsic description of the same flag.

  procedure Test_Flag ( file : in file_type;
                        extrinsic_flag,intrinsic_flag : in Matrix );

  -- DESCRIPTION :
  --   As we take points from larger spaces in the flag, chosen at
  --   random from the generators in the intrinsic flag, we will see
  --   that fewer and fewer hyperplanes in the extrinsic flag are zero.
  --   More precisely, a point from a d-space will only satisfy the
  --   first n-d equations in extrinsic_flag.
  --   This procedure generates random points and prints their value
  --   at the extrinsic hyperplanes to file.

end Standard_Flag_Representations;
