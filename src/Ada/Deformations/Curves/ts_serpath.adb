with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Series_Path_Trackers;

procedure ts_serpath is

-- DESCRIPTION :
--   Developing path trackers with power series.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the path tracker.

    ans : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Series_Path_Trackers.Standard_Main;
      when '1' => Series_Path_Trackers.DoblDobl_Main;
      when '2' => Series_Path_Trackers.QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpath;
