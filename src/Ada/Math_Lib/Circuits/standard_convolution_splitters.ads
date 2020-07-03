with Standard_Floating_VecVecVecs;
with Standard_Speelpenning_Convolutions;
with Standard_Coefficient_Convolutions;

package Standard_Convolution_Splitters is

-- DESCRIPTION :
--   The functions in this package split the complex coefficients in
--   convolution circuits for more efficient evaluation and differentiation.

  function Split ( c : Standard_Speelpenning_Convolutions.Circuit )
                 return Standard_Coefficient_Convolutions.Circuit;
  function Split ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit )
                 return Standard_Coefficient_Convolutions.Link_to_Circuit;
  function Split ( c : Standard_Speelpenning_Convolutions.Circuits )
                 return Standard_Coefficient_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns the circuit(s) with the same content as c, but with all
  --   complex coefficients separated into real and imaginary parts.

  procedure Split
              ( p : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                rp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ip : out Standard_Floating_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Splits the complex coefficients in the power table p in real and
  --   imaginary parts, returned in rp and ip respectively.

  function Split ( s : Standard_Speelpenning_Convolutions.System )
                 return Standard_Coefficient_Convolutions.System;
  function Split ( s : Standard_Speelpenning_Convolutions.Link_to_System )
                 return Standard_Coefficient_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Returns the system with the same content as s, with workspace
  --   complex coefficients separated into real and imaginary parts.

end Standard_Convolution_Splitters;
