with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Abstract_Ring_io;
with Generic_Vectors;
with Generic_Vectors_io;
with Generic_VecVecs;
with Generic_Lists_of_Vectors;

generic

  with package Ring_io is new Abstract_Ring_io(<>);
  with package Vectors is new Generic_Vectors(Ring_io.Ring);
  with package Vectors_io is new Generic_Vectors_io(Ring_io,Vectors);
  with package VecVecs is new Generic_VecVecs(Ring_io.Ring,Vectors);
  with package Lists is
         new Generic_Lists_of_Vectors(Ring_io.Ring,Vectors,VecVecs);

package Generic_Lists_of_Vectors_io is

-- DESCRIPTION :
--   Input/Output of lists of links to vectors.

  use Lists;

  procedure get ( n,m : in natural32; L : out List );
  procedure get ( file : in file_type;
                  n,m : in natural32; L : out List );

  -- DESCRIPTION :
  --   Reads m vectors of length n from standard output or from file.

  procedure put ( L : in List );
  procedure put ( file : in file_type; L : in List );

  -- DESCRIPTION :
  --   Writes the vectors in l on standard output or on file.

end Generic_Lists_of_Vectors_io;
