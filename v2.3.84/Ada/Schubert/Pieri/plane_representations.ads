with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Vectors;

package Plane_Representations is

-- DESCRIPTION :
--   This package allows to convert between various representations
--   of solution planes.

  function Localize ( locmap : Standard_Natural_Matrices.Matrix;
                      plamat : Standard_Complex_Matrices.Matrix )
                    return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the representation of the plane in local coordinates.

  function Vector_Rep ( plamat : Standard_Complex_Matrices.Matrix )
                      return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts the matrix representation of the plane into a vector
  --   representation, ordering the variables lexicographically.

  function Vector_Rep ( locmap : Standard_Natural_Matrices.Matrix;
                        plamat : Standard_Complex_Matrices.Matrix )
                      return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector representation of the plane, skipping the
  --   entries that correspond to zeros and ones in the localization.

  function Matrix_Rep ( locmap : Standard_Natural_Matrices.Matrix;
                        plavec : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Converts from vector to matrix representation of a plane.

  -- REQUIRED :
  --   The vector plavec contains exactly as many elements as there are free
  --   elements in the localization map.

  -- ON ENTRY :
  --   locmap           localization map, 0,1 for I, 2 for free variables;
  --   plavec           vector representation of a p-plane.

 -- ON RETURN : matrix representation as plane with the localization map.

end Plane_Representations;
