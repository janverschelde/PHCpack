with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuit;
                    deg : integer32 := 0 )
                  return Standard_Speelpenning_Convolutions.Circuit;
  function AbsVal ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit;
                    deg : integer32 := 0 )
                  return Standard_Speelpenning_Convolutions.Link_to_Circuit;
  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuits;
                    deg : integer32 := 0 )
                  return Standard_Speelpenning_Convolutions.Circuits;
  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuit;
                    deg : integer32 := 0 )
                  return DoblDobl_Speelpenning_Convolutions.Circuit;
  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    deg : integer32 := 0 )
                  return DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuits;
                    deg : integer32 := 0 )
                  return DoblDobl_Speelpenning_Convolutions.Circuits;
  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuit;
                    deg : integer32 := 0 )
                  return QuadDobl_Speelpenning_Convolutions.Circuit;
  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    deg : integer32 := 0 )
                  return QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuits;
                    deg : integer32 := 0 )
                  return QuadDobl_Speelpenning_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns a copy of c, except for the coefficients,
  --   which are replaced by the radii of the coefficients of c,
  --   up to the given degree index deg.

  function Residual_Convolution_System
             ( s : Standard_Speelpenning_Convolutions.System;
               deg : integer32 := 0 )
             return Standard_Speelpenning_Convolutions.System;
  function Residual_Convolution_System
             ( s : Standard_Speelpenning_Convolutions.Link_to_System;
               deg : integer32 := 0 )
             return Standard_Speelpenning_Convolutions.Link_to_System;
  function Residual_Convolution_System
             ( s : DoblDobl_Speelpenning_Convolutions.System;
               deg : integer32 := 0 )
             return DoblDobl_Speelpenning_Convolutions.System;
  function Residual_Convolution_System
             ( s : DoblDobl_Speelpenning_Convolutions.Link_to_System;
               deg : integer32 := 0 )
             return DoblDobl_Speelpenning_Convolutions.Link_to_System;
  function Residual_Convolution_System
             ( s : QuadDobl_Speelpenning_Convolutions.System;
               deg : integer32 := 0 )
             return QuadDobl_Speelpenning_Convolutions.System;
  function Residual_Convolution_System
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System;
               deg : integer32 := 0 )
             return QuadDobl_Speelpenning_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Returns a copy of the system s where the circuits are the
  --   result of the above AbsVal function.
  --   The work space of s is not copied, but the system on return
  --   has allocated work space of the proper dimensions.

  procedure Update_Radii_of_Constants
              ( radc : in Standard_Speelpenning_Convolutions.Circuit; 
                c : in Standard_Speelpenning_Convolutions.Circuit );
  procedure Update_Radii_of_Constants
              ( radc : in DoblDobl_Speelpenning_Convolutions.Circuit; 
                c : in DoblDobl_Speelpenning_Convolutions.Circuit );
  procedure Update_Radii_of_Constants
              ( radc : in QuadDobl_Speelpenning_Convolutions.Circuit; 
                c : in QuadDobl_Speelpenning_Convolutions.Circuit );
  procedure Update_Radii_of_Constants
              ( radc : in Standard_Speelpenning_Convolutions.Link_to_Circuit; 
                c : in Standard_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Update_Radii_of_Constants
              ( radc : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit; 
                c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Update_Radii_of_Constants
              ( radc : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit; 
                c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit );
  procedure Update_Radii_of_Constants
              ( radc : in Standard_Speelpenning_Convolutions.Circuits; 
                c : in Standard_Speelpenning_Convolutions.Circuits );
  procedure Update_Radii_of_Constants
              ( radc : in DoblDobl_Speelpenning_Convolutions.Circuits; 
                c : in DoblDobl_Speelpenning_Convolutions.Circuits );
  procedure Update_Radii_of_Constants
              ( radc : in QuadDobl_Speelpenning_Convolutions.Circuits; 
                c : in QuadDobl_Speelpenning_Convolutions.Circuits );

  -- DESCRIPTION :
  --   Updates the constant coefficients of radc with the radii of c.

  procedure Update_Radii_of_Constants
              ( rads : in Standard_Speelpenning_Convolutions.Link_to_System; 
                s : in Standard_Speelpenning_Convolutions.Link_to_System );
  procedure Update_Radii_of_Constants
              ( rads : in DoblDobl_Speelpenning_Convolutions.Link_to_System; 
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System );
  procedure Update_Radii_of_Constants
              ( rads : in QuadDobl_Speelpenning_Convolutions.Link_to_System; 
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System );

  -- DESCRIPTION :
  --   Updates the constant coefficients of rads with the radii of s.

end Residual_Convolution_Circuits;
