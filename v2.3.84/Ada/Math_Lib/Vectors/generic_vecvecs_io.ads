with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Abstract_Ring_io;
with Generic_Vectors;
with Generic_Vectors_io;
with Generic_VecVecs;

generic

  with package Ring_io is new Abstract_Ring_io(<>);
  with package Vectors is new Generic_Vectors(Ring_io.Ring);
  with package Vectors_io is new Generic_Vectors_io(Ring_io,Vectors);
  with package VecVecs is new Generic_VecVecs(Ring_io.Ring,Vectors);

package Generic_VecVecs_io is

-- DESCRIPTION :
--   Provides input/output routines for vectors with any entries.

  use VecVecs;

  procedure get ( n : in natural32; v : in out VecVec );
  procedure get ( file : in file_type; n : in natural32; v : in out VecVec );

  -- DESCRIPTION :
  --   Numbers will be read from standard input or from file,
  --   until all entries of v are filled with vectors of range 1..n.
  --   The numbers must be separated by spaces or line breaks.

  procedure get ( n1,n2 : in natural32; v : in out Link_to_VecVec );
  procedure get ( file : in file_type; n1,n2 : in natural32;
                  v : in out Link_to_VecVec );

  -- DESCRIPTION :
  --   The vector on return will be of range 1..n1 and will be filled
  --   with vectors of range 1..n2, with numbers read from standard
  --   input or from file.
  --   The numbers must be separated by spaces or line breaks.

  procedure put ( v : in VecVec );
  procedure put ( file : in file_type; v : in VecVec );
  procedure put ( v : in Link_to_VecVec );
  procedure put ( file : in file_type; v : in Link_to_VecVec );

  -- DESCRIPTION :
  --   The vector of vectors v is written on standard output or on file.
  --   The elements of v are written on separate lines.
  --   The elements of the elements of v appear on the same line and are 
  --   separated by a space.

  procedure put_line ( v : in VecVec );
  procedure put_line ( file : in file_type; v : in VecVec );
  procedure put_line ( v : in Link_to_VecVec );
  procedure put_line ( file : in file_type; v : in Link_to_VecVec );

  -- DESCRIPTION :
  --   The vector of vectors v is written on standard output or on file.
  --   The elements are written on separate lines.

  procedure put ( v : in VecVec; dp : in natural32 );
  procedure put ( file : in file_type; v : in VecVec; dp : in natural32 );
  procedure put ( v : in Link_to_VecVec; dp : in natural32 );
  procedure put ( file : in file_type;
                  v : in Link_to_VecVec; dp : in natural32 );

  -- DESCRIPTION :
  --   The vector of vectors v is written on standard output or on file.
  --   The elements of v are written on separate lines with dp decimal places.
  --   The elements of the elements of v appear on the same line and are 
  --   separated by a space.

  procedure put_line ( v : in VecVec; dp : in natural32 );
  procedure put_line ( file : in file_type; v : in VecVec; dp : in natural32 );
  procedure put_line ( v : in Link_to_VecVec; dp : in natural32 );
  procedure put_line ( file : in file_type;
                       v : in Link_to_VecVec; dp : in natural32 );

  -- DESCRIPTION :
  --   The vector of vectors v is written on standard output or on file.
  --   The elements are written on separate lines with dp decimal places.

end Generic_VecVecs_io;
