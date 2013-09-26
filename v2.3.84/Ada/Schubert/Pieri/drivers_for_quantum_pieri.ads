with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Drivers_for_Quantum_Pieri is

-- DESCRIPTION :
--   The quantum Pieri homotopies compute curves of p-planes which meet
--   given m-planes at prescribed interpolation points.
--   The two interactive drivers differ in the optional file argument.

   procedure Driver_for_Quantum_Pieri ( n,d,q : in natural32 );
   procedure Driver_for_Quantum_Pieri
              ( file : in file_type; n,d,q : in natural32 );

  -- DESCRIPTION :
  --   Call the quantum Pieri root counts and homotopies.

  -- ON ENTRY :
  --   file     output file (optional) must be opened for writing;
  --   n        ambient dimension;
  --   d        dimension of the solution planes;
  --   q        degree of the solution maps.

end Drivers_for_Quantum_Pieri;
