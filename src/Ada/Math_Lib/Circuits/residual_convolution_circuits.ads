with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;

package Residual_Convolution_Circuits is

-- DESCRIPTION :
--   To compute mixed residuals, mixing relative and absolute criteria,
--   the evaluation is required at a system where all coefficients are
--   replaced by their radius.  The functions in this package return copies
--   of convolution circuits with the radii of the coefficients of the input
--   circuits, in double, double double, and quad double precision.
  
  procedure AbsVal ( c : in out Standard_Complex_Numbers.Complex_Number );
  procedure AbsVal ( c : in out DoblDobl_Complex_Numbers.Complex_Number );
  procedure AbsVal ( c : in out QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Assigns to c its radius.

  procedure AbsVal ( v : in out Standard_Complex_Vectors.Vector );
  procedure AbsVal ( v : in Standard_Complex_Vectors.Link_to_Vector );
  procedure AbsVal ( v : in Standard_Complex_VecVecs.VecVec );
  procedure AbsVal ( v : in out DoblDobl_Complex_Vectors.Vector );
  procedure AbsVal ( v : in DoblDobl_Complex_Vectors.Link_to_Vector );
  procedure AbsVal ( v : in DoblDobl_Complex_VecVecs.VecVec );
  procedure AbsVal ( v : in out QuadDobl_Complex_Vectors.Vector );
  procedure AbsVal ( v : in QuadDobl_Complex_Vectors.Link_to_Vector );
  procedure AbsVal ( v : in QuadDobl_Complex_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Replaces all numbers in v by their radius.

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuit )
                  return Standard_Speelpenning_Convolutions.Circuit;
  function AbsVal ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit )
                  return Standard_Speelpenning_Convolutions.Link_to_Circuit;
  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuits )
                  return Standard_Speelpenning_Convolutions.Circuits;
  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
                  return DoblDobl_Speelpenning_Convolutions.Circuit;
  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit )
                  return DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuits )
                  return DoblDobl_Speelpenning_Convolutions.Circuits;
  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
                  return QuadDobl_Speelpenning_Convolutions.Circuit;
  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
                  return QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
                  return QuadDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns a copy of c. except for the coefficients,
  --   which are replaced by the radii of the coefficients of .

end Residual_Convolution_Circuits;
