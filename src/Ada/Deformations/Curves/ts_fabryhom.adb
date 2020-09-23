with Communications_with_User;
with Standard_Fabry_on_Homotopy;
with DoblDobl_Fabry_on_Homotopy;
with QuadDobl_Fabry_on_Homotopy;

procedure ts_fabryhom is

-- DESCRIPTION :
--   Tests the Newton-Fabry convergence radius computation
--   for artificial or natural-parameter homotopies.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and then launches the computation in that precision.

    prc : character;

  begin
    prc := Communications_with_User.Prompt_for_Precision;
    case prc is
      when '0' => Standard_Fabry_on_Homotopy.Main;
      when '1' => DoblDobl_Fabry_on_Homotopy.Main;
      when '2' => QuadDobl_Fabry_on_Homotopy.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_fabryhom;
