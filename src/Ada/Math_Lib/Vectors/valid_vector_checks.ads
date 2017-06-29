with Standard_Floating_Vectors;
with Standard_Complex_Vectors;

package Valid_Vector_Checks is

-- DESCRIPTION :
--   A vector is valid if all its components are valid numbers.

  function is_valid ( v : Standard_Floating_Vectors.Vector ) return boolean;
  function is_valid ( v : Standard_Complex_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if all components of v are valid numbers,
  --   returns false otherwise.

end Valid_Vector_Checks;
