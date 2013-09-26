with text_io;                           use text_io;
with Complex_Polynomial_Matrices;       use Complex_Polynomial_Matrices;

package Complex_Polynomial_Matrices_io is

-- DESCRIPTION :
--   This package offers input/output routines for matrices of polynomials
--   in one variable with complex floating-point coefficients.

  procedure Interactive_get ( lpm : out Link_to_Polynomial_Matrix );

  -- DESCRIPTION :
  --   Prompts the user for the dimensions, the degrees, and the coefficients
  --   of each polynomial in the matrix.

  procedure Interactive_get ( pm : in out Polynomial_Matrix );

  -- DESCRIPTION :
  --   Interactive reading of a polynomial matrix of fixed dimenion.

  procedure get ( lpm : out Link_to_Polynomial_Matrix );
  procedure get ( file : file_type; lpm : out Link_to_Polynomial_Matrix );

  -- DESCRIPTION :
  --   Reads from standard input or from file first the dimensions of the
  --   matrix and then, for every polynomial, first its degree, followed
  --   by as many coefficients as the degree of the polynomial plus one.
  --   The empty polynomial has degree minus one.

  procedure get ( pm : in out Polynomial_Matrix );
  procedure get ( file : file_type; pm : in out Polynomial_Matrix );

  -- DESCRIPTION :
  --   Reads a polynomial matrix of fixed dimension.

  procedure put ( pm : in Polynomial_Matrix );
  procedure put ( file : in file_type; pm : in Polynomial_Matrix );

  -- DESCRIPTION :
  --   Writes the polynomial matrix pm on standard output or on file.

  procedure put ( apm : in Array_of_Polynomial_Matrices );
  procedure put ( file : in file_type; apm : in Array_of_Polynomial_Matrices );

  -- DESCRIPTION :
  --   Writes an array of polynomial matrices on standard output or on file.

end Complex_Polynomial_Matrices_io;
