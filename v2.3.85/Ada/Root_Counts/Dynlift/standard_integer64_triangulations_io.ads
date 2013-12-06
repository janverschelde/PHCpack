with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Lists_of_Integer64_Vectors;         use Lists_of_Integer64_Vectors;
with Standard_Integer64_Triangulations;  use Standard_Integer64_Triangulations;

package Standard_Integer64_Triangulations_io is

-- DESCRIPTION :
--   Input/output of triangulations of polytopes spanned by integer vertices.

  procedure get ( t : in out Triangulation );
  procedure get ( n,m : in natural32; t : in out Triangulation );
  procedure get ( file : in file_type; t : in out Triangulation );
  procedure get ( file : in file_type; n,m : in natural32;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Reads first the dimension n and the number of simplices m.
  --   if they are not specified as parameter.
  --   Either from standard input or from file, m times n integer vectors
  --   of length are read.

  procedure put ( n : in natural32; t : in Triangulation );
  procedure put ( n : in natural32; t : in Triangulation; v : out natural64 );
  procedure put ( file : in file_type;
                  n : in natural32; t : in Triangulation );
  procedure put ( file : in file_type;
                  n : in natural32; t : in Triangulation; v : out natural64 );

  -- DESCRIPTION :
  --   Writes the simplices in the triangulation on standard output
  --   or on file.  When the parameter `v' is supplied, the volume
  --   will be computed and returned.  Also, more text banners are provided.

  procedure put ( n : natural32; t : in Triangulation;
                  convecs : in out List; v : out natural64 );
  procedure put ( file : in file_type; n : natural32; t : in Triangulation;
                  convecs : in out List; v : out natural64 );

  -- DESCRIPTION :
  --   Also the connectivity vectors for each simplex will be written.
  --   A connectivity vector cv for a simplex s is defined as follows:
  --    cv(i) = 0 if Neighbor(s,i) = Null_Simplex
  --    cv(i) = k if Neighbor(s,i) /= Null_Simplex
  --         and Position(t,Neighbor(s,i)) = k.

end Standard_Integer64_Triangulations_io;
