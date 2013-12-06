with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Drivers_to_Breakup_Components;      use Drivers_to_Breakup_Components;

procedure ts_breakup is

  ans : character;

begin
  new_line;
  put_line("Classifying generic points to break up components.");
  new_line;
  put_line("MENU for decomposing equidimensional solution set :");
  put_line("  1. incremental use interpolation filters to group points;");
  put_line("  2. predict grouping with monodromy group, verify with filters.");
  put("Type 1 or 2 to select : "); Ask_Alternative(ans,"12");
  case ans is
    when '1' => Breakup_with_Interpolation_Filters;
    when '2' => Breakup_with_Monodromy_Group_Actions;
    when others => null;
  end case;
end ts_breakup;
