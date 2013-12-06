with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Abstract_Ring_io;
with Generic_Vectors;
with Generic_Vectors_io;
with Generic_VecVecs;
with Generic_Lists_of_Vectors;
with Generic_Lists_of_Vectors_io;
with Generic_Arrays_of_Vector_Lists;

generic

  with package Ring_io is new Abstract_Ring_io(<>);
  with package Vectors is new Generic_Vectors(Ring_io.Ring);
  with package Vectors_io is new Generic_Vectors_io(Ring_io,Vectors);
  with package VecVecs is new Generic_VecVecs(Ring_io.Ring,Vectors);
  with package Lists is
         new Generic_Lists_of_Vectors(Ring_io.Ring,Vectors,VecVecs);
  with package Lists_io is
         new Generic_Lists_of_Vectors_io
               (Ring_io,Vectors,Vectors_io,VecVecs,Lists);
  with package Arrays is
         new Generic_Arrays_of_Vector_Lists(Ring_io.Ring,Vectors,VecVecs,Lists);

package Generic_Arrays_of_Vector_Lists_io is

-- DESCRIPTION :
--   Input/Output of arrays of lists of links to vectors.

  use Arrays;

  procedure get ( al : in out Link_to_Array_of_Lists );
  procedure get ( n : in natural32; m : in Standard_Natural_Vectors.Vector;
                  al : out Array_of_Lists );
  procedure get ( file : in file_type; al : in out Link_to_Array_of_Lists );
  procedure get ( file : in file_type;
                  n : in natural32; m : in Standard_Natural_Vectors.Vector;
                  al : out Array_of_Lists );

  -- DESCRIPTION :
  --   Reads a number of lists from standard input or from file;
  --   the ith list must contain m(i) vectors of length n.
  --   If n and m are not provided, then they are first read.

  procedure put ( al : in Array_of_Lists );
  procedure put ( file : in file_type; al : in Array_of_Lists );

  -- DESCRIPTION :
  --   Writes the lists in al on standard output or on file.

end Generic_Arrays_of_Vector_Lists_io;
