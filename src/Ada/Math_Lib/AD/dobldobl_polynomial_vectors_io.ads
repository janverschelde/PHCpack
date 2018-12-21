with text_io;                           use text_io;
with DoblDobl_Polynomial_Vectors;       use DoblDobl_Polynomial_Vectors;

package DoblDobl_Polynomial_Vectors_io is

-- DESCRIPTION :
--   A basic output procedure is defined by this package.

  procedure put ( v : in Polynomial_Vector );
  procedure put ( file : in file_type; v : in Polynomial_Vector );
  procedure put ( v : in Link_to_Polynomial_Vector );
  procedure put ( file : in file_type; v : in Link_to_Polynomial_Vector );

  -- DESCRIPTION :
  --   Writes the polynomial vector to standard output or to file.

  procedure put ( s : in System );
  procedure put ( file : in file_type; s : in System );
  procedure put ( s : in Link_to_System );
  procedure put ( file : in file_type; s : in Link_to_System );

  -- DESCRIPTION :
  --   Writes the contents of the system to standard output or to file.

end DoblDobl_Polynomial_Vectors_io;
