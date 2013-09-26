with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer64_Simplices;       use Standard_Integer64_Simplices;

package Standard_Integer64_Simplices_io is

-- DESCRIPTION :
--   Input/output routines for simplices spanned by integer vertices.

  procedure get ( s : in out Simplex );
  procedure get ( n : in natural32; s : in out Simplex );
  procedure get ( file : in file_type; s : in out Simplex );
  procedure get ( file : in file_type; n : in natural32; s : in out Simplex );

   -- DESCRIPTION :
   --   Reads the dimension n if not specified as parameter,
   --   and then n integer vectors of length n, from standard input
   --   or from file.

  procedure put ( s : in Simplex );
  procedure put ( file : in file_type; s : in Simplex );

   -- DESCRIPTION :
   --   Writes the n vectors that span the simplex on standard output
   --   or on file.

end Standard_Integer64_Simplices_io;
