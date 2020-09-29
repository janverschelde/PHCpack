with Communications_with_User;          use Communications_with_User;
with Test_Standard_Pade_Approximants;
with Test_DoblDobl_Pade_Approximants;
with Test_QuadDobl_Pade_Approximants;

procedure ts_serpade is

-- DESCRIPTION :
--   Main test on rational approximants.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and 
  --   then calls the main test procedure for that precision.

    prc : constant character := Prompt_for_Precision;

  begin
    case prc is 
      when '0' => Test_Standard_Pade_Approximants.Main;
      when '1' => Test_DoblDobl_Pade_Approximants.Main;
      when '2' => Test_QuadDobl_Pade_Approximants.Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpade;
