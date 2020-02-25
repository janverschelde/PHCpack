with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Abstract_Ring_io;
with Generic_Vectors;
with Generic_Matrices;
with Generic_Matrices_io;
with Generic_VecMats;

generic

  with package Ring_io is new Abstract_Ring_io(<>);
  with package Vectors is new Generic_Vectors(Ring_io.Ring);
  with package Matrices is new Generic_Matrices(Ring_io.Ring,Vectors);
  with package Matrices_io is new Generic_Matrices_io(Ring_io,Vectors,Matrices);
  with package VecMats is new Generic_VecMats(Ring_io.Ring,Vectors,Matrices);

package Generic_VecMats_io is

-- DESCRIPTION :
--   Provides input/output routines for vector of matrices with any entries.

  use VecMats;

  procedure get ( v : in out VecMat );
  procedure get ( file : in file_type; v : in out VecMat );

  -- DESCRIPTION :
  --   Numbers will be read from standard input or from file,
  --   until all entries of v are filled with n1-by-n2  matrices,
  --   where the dimensions n1 and n2 must be given on input. 
  --   The numbers must be separated by spaces or line breaks.

  procedure get ( n,n1,n2 : in natural32; v : in out Link_to_VecMat );
  procedure get ( file : in file_type; n,n1,n2 : in natural32;
                  v : in out Link_to_VecMat );

  -- DESCRIPTION :
  --   The vector on return will be of range 1..n and will be filled
  --   with matrices of ranges 1..n1,1..n2, with numbers read from standard
  --   input or from file.
  --   The numbers must be separated by spaces or line breaks.

  procedure put ( v : in VecMat );
  procedure put ( file : in file_type; v : in VecMat );
  procedure put ( v : in Link_to_VecMat );
  procedure put ( file : in file_type; v : in Link_to_VecMat );

  -- DESCRIPTION :
  --   The vector of matrices v is written on standard output or on file.
  --   The elements of v are written on separate lines, separated by 
  --   white lines.

  procedure put ( v : in VecMat; dp : in natural32 );
  procedure put ( file : in file_type; v : in VecMat; dp : in natural32 );
  procedure put ( v : in Link_to_VecMat; dp : in natural32 );
  procedure put ( file : in file_type;
                  v : in Link_to_VecMat; dp : in natural32 );

  -- DESCRIPTION :
  --   The vector of matrices v is written on standard output or on file.
  --   The elements of v are written on separate lines with dp decimal places.

end Generic_VecMats_io;
