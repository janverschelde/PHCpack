with Standard_Floating_VecVecVecs;
with QuadDobl_Speelpenning_Convolutions;
with QuadDobl_Coefficient_Convolutions;

package QuadDobl_Convolution_Splitters is

-- DESCRIPTION :
--   The functions in this package split the complex coefficients in
--   convolution circuits for more efficient evaluation and differentiation.

  function Split ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
                 return QuadDobl_Coefficient_Convolutions.Circuit;
  function Split ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
                 return QuadDobl_Coefficient_Convolutions.Link_to_Circuit;
  function Split ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
                 return QuadDobl_Coefficient_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns the circuit(s) with the same content as c, but with all
  --   complex coefficients separated into real high, imaginary high,
  --   real low, and imaginary low parts.

  procedure Split
              ( p : in QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                rp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ip : out Standard_Floating_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Splits the complex coefficients in the power table p in real and
  --   imaginary parts, returned respectively in the vectors rp and ip.

  function Split ( s : QuadDobl_Speelpenning_Convolutions.System )
                 return QuadDobl_Coefficient_Convolutions.System;
  function Split ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System )
                 return QuadDobl_Coefficient_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Returns the system with the same content as s, with workspace
  --   complex coefficients separated into real and imaginary parts.

end QuadDobl_Convolution_Splitters;
