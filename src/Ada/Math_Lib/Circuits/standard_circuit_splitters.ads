with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;

package Standard_Circuit_Splitters is

-- DESCRIPTION :
--   A circuit splitters takes a circuit, circuits, or a system
--   and splits the complex coefficients into real and imaginary
--   parts to form a coefficient circuit.

  function Split ( c : Standard_Complex_Circuits.Circuit )
                 return Standard_Coefficient_Circuits.Circuit;
  function Split ( c : Standard_Complex_Circuits.Link_to_Circuit )
                 return Standard_Coefficient_Circuits.Link_to_Circuit;

  -- DESCRIPTION :
  --   Splits the data in the circuit c 
  --   and returns the corresponding coefficient circuit.

  function Split ( c : Standard_Complex_Circuits.Circuits )
                 return Standard_Coefficient_Circuits.Circuits;

  -- DESCRIPTION :
  --   Splits the data in the circuit c
  --   and returns the corresponding coefficient circuit.

  function Split ( s : Standard_Complex_Circuits.System )
                 return Standard_Coefficient_Circuits.System;
  function Split ( s : Standard_Complex_Circuits.Link_to_System )
                 return Standard_Coefficient_Circuits.Link_to_System;

  -- DESCRIPTION :
  --   Splits the system into an equivalent coefficient system,
  --   allocates data for all work space.

end Standard_Circuit_Splitters;
