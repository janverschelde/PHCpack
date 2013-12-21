with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Deflation_Trees;           use Standard_Deflation_Trees;

package Standard_Deflation_Trees_io is

-- DESCRIPTION :
--   Writes the deflation trees to files.

  procedure Add_Multiplier_Symbols ( k,n : in natural32 );

  -- Adds n symbols to the symbol table,
  -- of the form lm[k,i], for i in 1..n.

  procedure Write_Deflated_System
              ( file : in file_type; p,dp : in Poly_Sys );

  -- DESCRIPTION :
  --   Writes the deflated system dp to file, deflated from the original p.

  procedure Write ( file : in file_type; t : in Node );
  procedure Write ( file : in file_type; name : in string; t : in Node );

  -- DESCRIPTION :
  --   If no name is provided, then all systems in the tree are
  --   written to one single file.
  --   If a name is provided, then the deflated systems are written
  --   to separate files, with suffixed _dkRr, with k the deflation
  --   number and r the rank.

end Standard_Deflation_Trees_io;
