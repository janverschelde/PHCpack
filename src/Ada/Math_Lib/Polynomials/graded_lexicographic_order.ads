with Standard_Natural_Vectors;
with Standard_Integer_Vectors;

package Graded_Lexicographic_Order is

-- DESCRIPTION :
--   This package provides a graded lexicographical ordening on vectors
--   of standard natural and integer numbers;

  function "<" ( v1,v2 : Standard_Natural_Vectors.Vector ) return boolean;
  function "<" ( v1,v2 : Standard_Natural_Vectors.Link_to_Vector )
               return boolean;
  function ">" ( v1,v2 : Standard_Natural_Vectors.Vector ) return boolean;
  function ">" ( v1,v2 : Standard_Natural_Vectors.Link_to_Vector )
               return boolean;

  function ">" ( v1,v2 : Standard_Integer_Vectors.Vector ) return boolean;
  function ">" ( v1,v2 : Standard_Integer_Vectors.Link_to_Vector )
               return boolean;
  function "<" ( v1,v2 : Standard_Integer_Vectors.Vector ) return boolean;
  function "<" ( v1,v2 : Standard_Integer_Vectors.Link_to_Vector )
               return boolean;

end Graded_Lexicographic_Order;
