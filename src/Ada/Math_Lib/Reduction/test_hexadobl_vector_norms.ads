package Test_HexaDobl_Vector_Norms is

-- DESCRIPTION :
--   Generates a random vector, computes a norm, divides every component 
--   of the vector by this norm, and then computes the norm of the 
--   normalized vector.  This second norm must equal one, 
--   within the accuracy of hexa double precision.

  procedure Test_Real_Two_Norm;

  -- DESCRIPTION :
  --   Prompts for a dimension, generates a real random vector, computes
  --   the 2-norm, divides every element in the vector by the 2-norm,
  --   and then computes the 2-norm of the normalized vector.

  procedure Test_Real_Sum_Norm;

  -- DESCRIPTION :
  --   Prompts for a dimension, generates a real random vector, computes
  --   the sum norm, divides every element in the vector by the sum norm,
  --   and then computes the sum norm of the normalized vector.

  procedure Test_Real_Max_Norm;

  -- DESCRIPTION :
  --   Prompts for a dimension, generates a real random vector, computes
  --   the max norm, divides every element in the vector by the max norm,
  --   and then computes the max norm of the normalized vector.

  procedure Test_Complex_Two_Norm;

  -- DESCRIPTION :
  --   Prompts for a dimension, generates a complex random vector, computes
  --   the 2-norm, divides every element in the vector by the 2-norm,
  --   and then computes the 2-norm of the normalized vector.

  procedure Test_Complex_Sum_Norm;

  -- DESCRIPTION :
  --   Prompts for a dimension, generates a complex random vector, computes
  --   the sum norm, divides every element in the vector by the sum norm,
  --   and then computes the sum norm of the normalized vector.

  procedure Test_Complex_Max_Norm;

  -- DESCRIPTION :
  --   Prompts for a dimension, generates a complex random vector, computes
  --   the max norm, divides every element in the vector by the max norm,
  --   and then computes the max norm of the normalized vector.

  procedure Main;

  -- DESCRIPTION :
  --   Displays a menu and prompts for a test.

end Test_HexaDobl_Vector_Norms;
