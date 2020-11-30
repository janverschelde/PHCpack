with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;

package Main_Quantum_Pieri is

-- DESCRIPTION :
--   The quantum Pieri homotopies compute curves of p-planes which meet
--   given m-planes at prescribed interpolation points.

   procedure Main ( n,d,q : in natural32 );
   procedure Main ( file : in file_type; n,d,q : in natural32 );

  -- DESCRIPTION :
  --   Calls the quantum Pieri root counts and homotopies.

  -- ON ENTRY :
  --   file     output file (optional) must be opened for writing;
  --   n        ambient dimension;
  --   d        dimension of the solution planes;
  --   q        degree of the solution maps.

end Main_Quantum_Pieri;
