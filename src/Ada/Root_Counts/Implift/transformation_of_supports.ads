with Standard_Integer_Matrices;
with Lists_of_Integer_Vectors;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;

package Transformation_of_Supports is

-- DESCRIPTION :
--   Supports with a span that has a positive corank
--   we can transform to eliminate as many variables as the corank.
--   The operations in this package are normally applied after the
--   calculation of the rank of the span of the supports.

  function Transform ( support : Lists_of_Integer_Vectors.List;
                       transfo : Standard_Integer_Matrices.Matrix )
                     return Lists_of_Integer_Vectors.List;

  -- DECRIPTION :
  --   Applies the transformation to the points in the list support,
  --   multiplying every point with the matrix transfo.
  --   On return is the transformed list.

  function Transform ( p : Standard_Complex_Laurentials.Poly;
                       transfo : Standard_Integer_Matrices.Matrix )
                     return Standard_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Applies the coordinate transformation defined by the matrix transfo
  --   to every monomial in the given polynomial p.
  --   On return is the transformed polynomial.

  function Transform ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                       transfo : Standard_Integer_Matrices.Matrix )
                     return Standard_Complex_Laur_Systems.Laur_Sys;

  -- DESCRIPTION :
  --   Applies the coordinate transformation defined by the matrix transfo
  --   to every monomial in the given Laurent system p.
  --   On return is the transformed Laurent system.

end Transformation_of_Supports;
