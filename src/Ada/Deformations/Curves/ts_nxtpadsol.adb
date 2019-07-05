with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Interactive_Pade_Trackers;          use Interactive_Pade_Trackers;

procedure ts_nxtpadsol is

-- DESCRIPTION :
--   Interactive test on the get_next() method to track paths
--   with a series Pade predictor.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and then launches
  --   the corresponding test.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
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
