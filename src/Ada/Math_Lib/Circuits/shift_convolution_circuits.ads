with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Shift_Convolution_Circuits is

-- DESCRIPTION :
--   Provides operations to shift the coefficients of convolution circuits.

  procedure Shift ( c,wrk : in out Standard_Complex_Vectors.Vector;
                    t : in double_float );
  procedure Shift ( c,wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float );
  procedure Shift ( c,wrk : in out Standard_Complex_Vectors.Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( c,wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( c,wrk : in out DoblDobl_Complex_Vectors.Vector;
                    t : in double_double );
  procedure Shift ( c,wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double );
  procedure Shift ( c,wrk : in out DoblDobl_Complex_Vectors.Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c,wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c,wrk : in out QuadDobl_Complex_Vectors.Vector;
                    t : in quad_double );
  procedure Shift ( c,wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double );
  procedure Shift ( c,wrk : in out QuadDobl_Complex_Vectors.Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c,wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   The coefficients in c of the series x(s) are shifted to 
  --   correspond to the coefficients of the series x(s-t),
  --   in double, double double, and quad double precision.

  -- ON ENTRY :
  --   c        coefficients of a power series;
  --   wrk      work space vector of same range as c;
  --   t        value of the shift.

  -- ON RETURN :
  --   c        shifted coefficients;
  --   wrk      work space equals c.

  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float );
  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float );
  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double );
  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double );
  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double );
  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double );
  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Shifts all the coefficients of the power series in c.

  -- ON ENTRY :
  --   c        a convolution circuit with coefficients as power series,
  --            all of the same degree;
  --   wrk      work space vector of same range as the coefficients in c;
  --   t        value of the shift.

  -- ON RETURN :
  --   c        circuit with shifted coefficients;
  --   wrk      modified work space.

  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float );
  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float );
  procedure Shift ( c : in out Standard_Speelpenning_Convolutions.Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double );
  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double );
  procedure Shift ( c : in out DoblDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double );
  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double );
  procedure Shift ( c : in out QuadDobl_Speelpenning_Convolutions.Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuits;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Shifts all the coefficients of the power series in c.

  -- ON ENTRY :
  --   c        convolution circuits with coefficients as power series,
  --            all of the same degree;
  --   wrk      work space vector of same range as the coefficients in c;
  --   t        value of the shift.

  -- ON RETURN :
  --   c        circuits with shifted coefficients;
  --   wrk      modified work space.

  procedure Shift ( s : in out Standard_Speelpenning_Convolutions.System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float );
  procedure Shift ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in double_float );
  procedure Shift ( s : in out Standard_Speelpenning_Convolutions.System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                    wrk : in Standard_Complex_Vectors.Link_to_Vector;
                    t : in Standard_Complex_Numbers.Complex_Number );
  procedure Shift ( s : in out DoblDobl_Speelpenning_Convolutions.System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double );
  procedure Shift ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in double_double );
  procedure Shift ( s : in out DoblDobl_Speelpenning_Convolutions.System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                    t : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( s : in out QuadDobl_Speelpenning_Convolutions.System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double );
  procedure Shift ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in quad_double );
  procedure Shift ( s : in out QuadDobl_Speelpenning_Convolutions.System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );
  procedure Shift ( s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                    wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                    t : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Shifts all the coefficients of the power series in s.

  -- ON ENTRY :
  --   s        system of circuits with coefficients as power series,
  --            all of the same degree;
  --   wrk      work space vector of same range as the coefficients in c;
  --   t        value of the shift.

  -- ON RETURN :
  --   s        circuit in s have shifted coefficients;
  --   wrk      modified work space.

end Shift_Convolution_Circuits;
