with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Test_Predictor_Convolutions is

  procedure Standard_Check_Solutions
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Given in chom is a system of convolution circuits for a homotopy
  --   and in sols a list of solutions, in double precision.
  --   Writes the result of the evaluation of each solution at the circuits.

  procedure DoblDobl_Check_Solutions
              ( chom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Given in chom is a system of convolution circuits for a homotopy
  --   and in sols a list of solutions, in double double precision.
  --   Writes the result of the evaluation of each solution at the circuits.

  procedure QuadDobl_Check_Solutions
              ( chom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Given in chom is a system of convolution circuits for a homotopy
  --   and in sols a list of solutions, in quad double precision.
  --   Writes the result of the evaluation of each solution at the circuits.

end Test_Predictor_Convolutions;
