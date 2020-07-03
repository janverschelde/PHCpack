with Standard_Floating_VecVecVecs;
with DoblDobl_Speelpenning_Convolutions;
with DoblDobl_Coefficient_Convolutions;

package DoblDobl_Convolution_Splitters is

-- DESCRIPTION :
--   The functions in this package split the complex coefficients in
--   convolution circuits for more efficient evaluation and differentiation.

  function Split ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
                 return DoblDobl_Coefficient_Convolutions.Circuit;
  function Split ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit )
                 return DoblDobl_Coefficient_Convolutions.Link_to_Circuit;
  function Split ( c : DoblDobl_Speelpenning_Convolutions.Circuits )
                 return DoblDobl_Coefficient_Convolutions.Circuits;

  -- DESCRIPTION :
  --   Returns the circuit(s) with the same content as c, but with all
  --   complex coefficients separated into real high, imaginary high,
  --   real low, and imaginary low parts.

  procedure Split
              ( p : in DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                rhp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec );

  -- DESCRIPTION :
  --   Splits the complex coefficients in the power table p in real high,
  --   imaginary high, real low, imaginary low parts, returned respectively
  --   in the vectors rhp, ihp, rlp, ilp.

  function Split ( s : DoblDobl_Speelpenning_Convolutions.System )
                 return DoblDobl_Coefficient_Convolutions.System;
  function Split ( s : DoblDobl_Speelpenning_Convolutions.Link_to_System )
                 return DoblDobl_Coefficient_Convolutions.Link_to_System;

  -- DESCRIPTION :
  --   Returns the system with the same content as s, with workspace
  --   complex coefficients separated into real and imaginary parts.

end DoblDobl_Convolution_Splitters;
