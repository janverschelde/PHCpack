with text_io;                            use text_io;
with Abstract_Ring_io;
with Generic_Vectors;
with Generic_Vectors_io;
with Generic_NesVecs;

generic

  with package Ring_io is new Abstract_Ring_io(<>);
  with package Vectors is new Generic_Vectors(Ring_io.Ring);
  with package Vectors_io is new Generic_Vectors_io(Ring_io,Vectors);
  with package NesVecs is new Generic_NesVecs(Ring_io.Ring,Vectors);

package Generic_NesVecs_io is

-- DESCRIPTION :
--   Provides input/output routines for nested vectors.

  use NesVecs;

  procedure get ( v : in out NesVec );
  procedure get ( file : in file_type; v : in out NesVec );
  procedure get ( v : in out Link_to_NesVec );
  procedure get ( file : in file_type; v : in out Link_to_NesVec );
  procedure get ( v : in out Array_of_NesVecs );
  procedure get ( file : in file_type; v : in out Array_of_NesVecs );

  -- DESCRIPTION :
  --   Reads nested vectors from file or from standard input.
  --   The first line should always contain the dimensions n, a, and b.
  --   Reading occurs depth first, corresponding to the output format.

  procedure put ( v : in NesVec );
  procedure put ( file : in file_type; v : in NesVec );
  procedure put ( v : in Link_to_NesVec );
  procedure put ( file : in file_type; v : in Link_to_NesVec );
  procedure put ( v : in Array_of_NesVecs );
  procedure put ( file : in file_type; v : in Array_of_NesVecs );

  -- DESCRIPTION :
  --   The nested vectors are written on standard output or on file.
  --   First the dimensions, n, a, and b are written, followed on the
  --   next lines by either the actual numbers (if n = 1) on separate lines,
  --   or recursively, by the other nested vectors. 

end Generic_NesVecs_io;
