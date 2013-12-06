with text_io;                            use text_io;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Triangulations;                     use Triangulations;

package Triangulations_io is

-- DESCRIPTION :
--   Input/output of triangulations of polytopes spanned by integer vertices.

  procedure get ( t : in out Triangulation );
  procedure get ( n,m : in natural; t : in out Triangulation );
  procedure get ( file : in file_type; t : in out Triangulation );
  procedure get ( file : in file_type; n,m : in natural;
                  t : in out Triangulation );

  -- DESCRIPTION :
  --   Reads first the dimension n and the number of simplices m.
  --   if they are not specified as parameter.
  --   Either from standard input or from file, m times n integer vectors
  --   of length are read.

  procedure put ( n : in natural; t : in Triangulation );
  procedure put ( n : in natural; t : in Triangulation; v : out natural );
  procedure put ( file : in file_type;
                  n : in natural; t : in Triangulation );
  procedure put ( file : in file_type;
                  n : in natural; t : in Triangulation; v : out natural );

  -- DESCRIPTION :
  --   Writes the simplices in the triangulation on standard output
  --   or on file.  When the parameter `v' is supplied, the volume
  --   will be computed and returned.  Also, more text banners are provided.

  procedure put ( n : natural; t : in Triangulation;
                  convecs : in out List; v : out natural );
  procedure put ( file : in file_type; n : natural; t : in Triangulation;
                  convecs : in out List; v : out natural );

  -- DESCRIPTION :
  --   Also the connectivity vectors for each simplex will be written.
  --   A connectivity vector cv for a simplex s is defined as follows:
  --    cv(i) = 0 if Neighbor(s,i) = Null_Simplex
  --    cv(i) = k if Neighbor(s,i) /= Null_Simplex
  --         and Position(t,Neighbor(s,i)) = k.

end Triangulations_io;
