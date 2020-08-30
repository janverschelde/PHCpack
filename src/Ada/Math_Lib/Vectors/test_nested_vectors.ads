with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_NesVecs;           use Standard_Complex_NesVecs;

package Test_Nested_Vectors is
 
-- DESCRIPTION :
--   Tests nested vectors, vectors of vectors of vectors, etc...

  procedure Interactive_Read ( nv : in out NesVec );

  -- DESCRIPTION :
  --   Interactive procedure to read a nested vector.

  procedure Interactive_Read ( n : in natural32; nv : in out Link_to_NesVec );

  -- DESCRIPTION :
  --   Interactive procedure to read a nested vector of dimension n.

  procedure Interactive_Read ( nv : in out Link_to_NesVec );

  -- DESCRIPTION :
  --   Interactive procedure to read a nested vector of any dimension.

  procedure Main;

  -- DESCRIPTION :
  --   Reads a multi-dimensional matrix and writes this matrix.

end Test_Nested_Vectors;
