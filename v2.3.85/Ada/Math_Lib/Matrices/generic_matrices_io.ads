with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Abstract_Ring_io;
with Generic_Vectors;
with Generic_Matrices;

generic

  with package Ring_io is new Abstract_Ring_io(<>);  use Ring_io.Ring;
  with package Vectors is new Generic_Vectors(Ring_io.Ring);
  with package Matrices is new Generic_Matrices(Ring_io.Ring,Vectors);

package Generic_Matrices_io is

-- DESCRIPTION :
--   Provides input/output routines for matrices with any entries.

  use Matrices;

  procedure get ( m : in out Matrix );
  procedure get ( file : in file_type; m : in out Matrix );

  procedure get ( m : in out Matrix; rw1,rw2 : in integer32 );
  procedure get ( file : in file_type;
                  m : in out Matrix; rw1,rw2 : in integer32 );

  -- DESCRIPTION :
  --   Reads an integer matrix m or m(rw1..rw2,m'range(2))
  --   from standard input or on file.

  procedure put ( m : in Matrix );
  procedure put ( file : in file_type; m : in Matrix );

  procedure put ( m : in Matrix; rw1,rw2 : in integer32 );
  procedure put ( file : in file_type;
                  m : in Matrix; rw1,rw2 : in integer32 );

  -- DESCRIPTION :
  --   Writes a matrix m or submatrix m(rw1..rw2,m'range(2))
  --   on standard output or on file.

  procedure put ( m : in Matrix; dp : in natural32 );
  procedure put ( file : in file_type; m : in Matrix; dp : in natural32 );

  procedure put ( m : in Matrix; rw1,rw2 : in integer32; dp : in natural32 );
  procedure put ( file : in file_type;
                  m : in Matrix; rw1,rw2 : in integer32; dp : in natural32 );

  -- DESCRIPTION :
  --   Writes a matrix m or submatrix m(rw1..rw2,m'range(2))
  --   on standard output or on file, with dp decimal places.

end Generic_Matrices_io;
