with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer32_Transformations; use Standard_Integer32_Transformations;

package Standard_Integer_Transformations_io is

-- DESCRIPTION :
--   Input/output routines for transformations.

  procedure get ( n : in natural32; t : out Transfo );
  procedure get ( file : in file_type; n : in natural32; t : out Transfo );

  -- DESCRIPTION :
  --   Reads n vectors from standard input or from file.
  --   These vectors are considered as the images of
  --   the basis vectors under the transformation t.

  procedure put ( t : in Transfo );
  procedure put ( file : in file_type; t : in Transfo );

  -- DESCRIPTION :
  --   Writes the images of the basis vectors under t.

end Standard_Integer_Transformations_io;
