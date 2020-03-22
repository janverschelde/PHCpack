with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Series_Path_Trackers;

procedure ts_serpath is

-- DESCRIPTION :
--   Developing path trackers with power series.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision
  --   and then launches the path tracker.

    ans : constant character := Prompt_for_Precision;
    vrb : integer32 := 0;

  begin
    new_line;
    put("Verbose level ? (0 if no verbose level) : ");
    get(vrb);
    case ans is
      when '0' => Series_Path_Trackers.Standard_Main(vrb);
      when '1' => Series_Path_Trackers.DoblDobl_Main(vrb);
      when '2' => Series_Path_Trackers.QuadDobl_Main(vrb);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serpath;
