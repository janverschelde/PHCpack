with text_io;                           use text_io;
with Standard_Dense_Matrix_Series2;     use Standard_Dense_Matrix_Series2;

package Standard_Dense_Matrix_Series2_io is

-- DESCRPTION :
--   Provides very basic output for series with matrix coefficients.

  procedure put ( A : in Matrix );
  procedure put ( file : in file_type; A : in Matrix );

  -- DESCRIPTION :
  --   Writes the coefficients of the Matrix series A to standard output,
  --   or to file.  The matrix coefficients are separated by a new line.

end Standard_Dense_Matrix_Series2_io;
