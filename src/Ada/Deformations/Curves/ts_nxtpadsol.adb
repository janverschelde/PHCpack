with Communications_with_User;           use Communications_with_User;
with Interactive_Pade_Trackers;          use Interactive_Pade_Trackers;

procedure ts_nxtpadsol is

-- DESCRIPTION :
--   Interactive test on the get_next() method to track paths
--   with a series Pade predictor.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and launches the test.

    ans : constant character := Prompt_for_Precision;

  begin
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_nxtpadsol;
