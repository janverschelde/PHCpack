with text_io;                           use text_io;
with Standard_Polynomial_Vectors;       use Standard_Polynomial_Vectors;

package Standard_Polynomial_Vectors_io is

-- DESCRIPTION :
--   Basic output procedures are defined by this package.

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

end Standard_Polynomial_Vectors_io;
