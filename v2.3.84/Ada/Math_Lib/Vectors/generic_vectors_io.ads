with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Abstract_Ring_io;
with Generic_Vectors;

generic

  with package Ring_io is new Abstract_Ring_io(<>);
  with package Vectors is new Generic_Vectors(Ring_io.Ring);

package Generic_Vectors_io is

-- DESCRIPTION :
--   Provides input/output routines for vectors with any entries.

  use Vectors;

  procedure get ( v : in out Vector );
  procedure get ( file : in file_type; v : in out Vector );

  -- DESCRIPTION :
  --   Numbers will be read from standard input or from file,
  --   until all entries of v are filled.
  --   The numbers must be separated by spaces.

  procedure get ( n : in natural32; v : in out Link_to_Vector );
  procedure get ( file : in file_type; n : in natural32;
                  v : in out Link_to_Vector );

  -- DESCRIPTION :
  --   The range of v on return will be 1..n, and it will be filled up
  --   with n numbers that are read from standard input or from file.
  --   The numbers must be separated by spaces or line breaks.

  procedure put ( v : in Vector );
  procedure put ( file : in file_type; v : in Vector );
  procedure put ( v : in Link_to_Vector );
  procedure put ( file : in file_type; v : in Link_to_Vector );

  -- DESCRIPTION :
  --   The vector v is written on standard output or on file.
  --   The elements appear on the same line and are separated by a space.

  procedure put_line ( v : in Vector );
  procedure put_line ( file : in file_type; v : in Vector );
  procedure put_line ( v : in Link_to_Vector );
  procedure put_line ( file : in file_type; v : in Link_to_Vector );

  -- DESCRIPTION :
  --   The vector v is written on standard output or on file.
  --   The elements are written on separate lines.

  procedure put ( v : in Vector; dp : in natural32 );
  procedure put ( file : in file_type; v : in Vector; dp : in natural32 );
  procedure put ( v : in Link_to_Vector; dp : in natural32 );
  procedure put ( file : in file_type;
                  v : in Link_to_Vector; dp : in natural32 );

  -- DESCRIPTION :
  --   The vector v is written on standard output or on file,
  --   with dp decimal places. 
  --   The elements appear on the same line and are separated by a space.

  procedure put_line ( v : in Vector; dp : in natural32 );
  procedure put_line ( file : in file_type; v : in Vector; dp : in natural32 );
  procedure put_line ( v : in Link_to_Vector; dp : in natural32 );
  procedure put_line ( file : in file_type;
                       v : in Link_to_Vector; dp : in natural32 );

  -- DESCRIPTION :
  --   The vector v is written on standard output or on file,
  --   with dp decimal places.
  --   The elements are written on separate lines.

end Generic_Vectors_io;
