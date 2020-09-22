with Communications_with_User;           use Communications_with_User;
with Test_Standard_Coeff_Convolutions;
with Test_DoblDobl_Coeff_Convolutions;
with Test_QuadDobl_Coeff_Convolutions;

procedure ts_speelcnv is

-- DESCRIPTION :
--   Tests the evaluation of the gradient of a polynomial in many variables,
--   in a power series of some fixed degree.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then calls the main test
  --   for that precision.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Test_Standard_Coeff_Convolutions.Main;
      when '1' => Test_DoblDobl_Coeff_Convolutions.Main;
      when '2' => Test_QuadDobl_Coeff_Convolutions.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_speelcnv;
