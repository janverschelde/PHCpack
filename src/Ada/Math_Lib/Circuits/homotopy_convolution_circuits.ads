with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with TripDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with PentDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors;
with DecaDobl_Complex_Vectors;
with HexaDobl_Complex_Vectors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with TripDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with PentDobl_Speelpenning_Convolutions;
with OctoDobl_Speelpenning_Convolutions;
with DecaDobl_Speelpenning_Convolutions;
with HexaDobl_Speelpenning_Convolutions;

package Homotopy_Convolution_Circuits is

-- DESCRIPTION :
--   A convolution circuit represents a homotopy if the coefficients
--   with positive powers in the series are nonzero.
--   The procedures in this package add a continuation parameter
--   to the coefficients of the convolution circuits.

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Add_Continuation_Parameter
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit );

  -- DESCRIPTION :
  --   To every coefficient adds 't' as the linear coefficient
  --   in all power series coefficients of the circuit.

  procedure Add_Continuation_Parameter
              ( c : in Standard_Speelpenning_Convolutions.Circuits );
  procedure Add_Continuation_Parameter
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits );
  procedure Add_Continuation_Parameter
              ( c : in TripDobl_Speelpenning_Convolutions.Circuits );
  procedure Add_Continuation_Parameter
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits );
  procedure Add_Continuation_Parameter
              ( c : in PentDobl_Speelpenning_Convolutions.Circuits );
  procedure Add_Continuation_Parameter
              ( c : in OctoDobl_Speelpenning_Convolutions.Circuits );
  procedure Add_Continuation_Parameter
              ( c : in DecaDobl_Speelpenning_Convolutions.Circuits );
  procedure Add_Continuation_Parameter
              ( c : in HexaDobl_Speelpenning_Convolutions.Circuits );

  -- DESCRIPTION :
  --   To every coefficient adds 't' as the linear coefficient
  --   in all power series coefficients of the circuits.

  procedure Add_Parameter_to_Constant
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                deg : in integer32 );
  procedure Add_Parameter_to_Constant
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit;
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
              ( s : in TripDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Add_Parameter_to_Constant
              ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Add_Parameter_to_Constant
              ( s : in PentDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Add_Parameter_to_Constant
              ( s : in OctoDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Add_Parameter_to_Constant
              ( s : in DecaDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Add_Parameter_to_Constant
              ( s : in HexaDobl_Speelpenning_Convolutions.Link_to_System );

  -- DESCRIPTION :
  --   Adds the continuation parameter to every circuit in s.

  -- REQUIRED : s /= null.

  procedure Set_Solution_Constant
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                z : in Standard_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DoblDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in TripDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in QuadDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in PentDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in OctoDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DecaDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in HexaDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the constant of c so that z becomes a solution.

  procedure Set_Solution_Constant
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                z : in Standard_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                z : in DoblDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in TripDobl_Speelpenning_Convolutions.Circuits;
                z : in TripDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                z : in QuadDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in PentDobl_Speelpenning_Convolutions.Circuits;
                z : in PentDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                z : in OctoDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                z : in DecaDobl_Complex_Vectors.Vector );
  procedure Set_Solution_Constant
              ( c : in HexaDobl_Speelpenning_Convolutions.Circuits;
                z : in HexaDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets the constant of each circuit in c so that z becomes a solution.

  procedure Newton_Homotopy
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                z : in Standard_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DoblDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in TripDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in TripDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in QuadDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in PentDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in PentDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in OctoDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in OctoDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in DecaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in DecaDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in HexaDobl_Speelpenning_Convolutions.Link_to_Circuit;
                z : in HexaDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Let y = Eval(c,z), then the constant of c becomes -y + t*y, so
  --   at t = 0, z is a solution of c, and at t = 1, z is the original c.

  -- REQUIRED : the degree of the power series in c is at least one.

  procedure Newton_Homotopy
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                z : in Standard_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                z : in DoblDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in TripDobl_Speelpenning_Convolutions.Circuits;
                z : in TripDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                z : in QuadDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in PentDobl_Speelpenning_Convolutions.Circuits;
                z : in PentDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                z : in OctoDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                z : in DecaDobl_Complex_Vectors.Vector );
  procedure Newton_Homotopy
              ( c : in HexaDobl_Speelpenning_Convolutions.Circuits;
                z : in HexaDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Constructs a Newton homotopy for all circuits in c.

  -- REQUIRED : the degree of all power series in c is at least one.

end Homotopy_Convolution_Circuits;
