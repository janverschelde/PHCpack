with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Homotopy_Convolution_Circuits is

-- DESCRIPTION :
--   A convolution circuit represents a homotopy if the coefficients
--   with positive powers in the series are nonzero.
--   The procedures in this package add a continuation parameter
--   to the coefficients of the convolution circuits.

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit );

  -- DESCRIPTION :
  --   To every coefficient adds 't' as the linear coefficient
  --   in all power series coefficients of the circuit.

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.
                       Convolution_Circuits );
  procedure Add_Continuation_Parameter
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Convolution_Circuits );
  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Convolution_Circuits );

  -- DESCRIPTION :
  --   To every coefficient adds 't' as the linear coefficient
  --   in all power series coefficients of the circuits.

  procedure Add_Parameter_to_Constant
              ( c : in Standard_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                deg : in integer32 );

  -- DESCRIPTION :
  --   Adds the continuation parameter t to the constant of c.

  -- REQUIRED : the coefficients in the power series of c
  --   have degree at least deg.

  procedure Add_Parameter_to_Constant
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System );
  procedure Add_Parameter_to_Constant
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Add_Parameter_to_Constant
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System );

  -- DESCRIPTION :
  --   Adds the continuation parameter to every circuit in s.

  -- REQUIRED : s /= null.

  procedure Set_Solution_Constant
              ( c : in Standard_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                z : in Standard_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                z : in DoblDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Convolution_Circuit;
                z : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the constant of c so that z becomes a solution.

  procedure Set_Solution_Constant
              ( c : in Standard_Speelpenning_Convolutions.Convolution_Circuits;
                z : in Standard_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.Convolution_Circuits;
                z : in DoblDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.Convolution_Circuits;
                z : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the constant of each circuit in c so that z becomes a solution.

end Homotopy_Convolution_Circuits;
