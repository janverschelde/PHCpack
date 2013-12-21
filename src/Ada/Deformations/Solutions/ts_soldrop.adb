with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Symbol_Table,Symbol_Table_io;
with Solution_Drops;

procedure ts_soldrop is

  function Prompt_for_Drop_Variable return natural32 is

  -- DESCRIPTION :
  --   Returns the index of the coordinate to be dropped.
  --   If the symbol table is empty, then the user will be asked
  --   to enter a number, otherwise the user can enter the name
  --   of the variable for the drop.

    res : natural32 := 0;

  begin
    if Symbol_Table.Empty then
      put("give component to drop : "); get(res);
    else
      put("The variable names : "); Symbol_Table_io.Write; new_line;
      declare
        sb : Symbol_Table.symbol;
      begin
        loop
          put("give the symbol to be dropped : "); Symbol_Table_io.get(sb);
          res := Symbol_Table.Get(sb);
          exit when (res /= 0);
          put("symbol not found, please try again...");
        end loop;
      end;
    end if;
    return res;
  end Prompt_for_Drop_Variable;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Tests the dropping of a coordinate from a list of solutions
  --   in standard double precision.

    sols,dropped : Standard_Complex_Solutions.Solution_List;
    k : natural32;

  begin
    new_line;
    Read(sols);
    new_line;
    k := Prompt_for_Drop_Variable;
    dropped := Solution_Drops.Drop(sols,k);
    Symbol_Table.Remove(k);
    put("The solution list with coordinate "); put(k,1);
    put_line(" dropped :");
    put(standard_output,Standard_Complex_Solutions.Length_Of(dropped),
        natural32(Standard_Complex_Solutions.Head_Of(dropped).n),dropped);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Tests the dropping of a coordinate from a list of solutions
  --   in double double precision.

    sols,dropped : DoblDobl_Complex_Solutions.Solution_List;
    k : natural32;

  begin
    new_line;
    Read(sols);
    new_line;
    k := Prompt_for_Drop_Variable;
    dropped := Solution_Drops.Drop(sols,k);
    Symbol_Table.Remove(k);
    put("The solution list with coordinate "); put(k,1);
    put_line(" dropped :");
    put(standard_output,DoblDobl_Complex_Solutions.Length_Of(dropped),
        natural32(DoblDobl_Complex_Solutions.Head_Of(dropped).n),dropped);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Tests the dropping of a coordinate from a list of solutions
  --   in quad double precision.

    sols,dropped : QuadDobl_Complex_Solutions.Solution_List;
    k : natural32;

  begin
    new_line;
    Read(sols);
    new_line;
    k := Prompt_for_Drop_Variable;
    dropped := Solution_Drops.Drop(sols,k);
    Symbol_Table.Remove(k);
    put("The solution list with coordinate "); put(k,1);
    put_line(" dropped :");
    put(standard_output,QuadDobl_Complex_Solutions.Length_Of(dropped),
        natural32(QuadDobl_Complex_Solutions.Head_Of(dropped).n),dropped);
  end QuadDobl_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for dropping a coordinate of a solution...");
    put_line("  1. standard double precision;");
    put_line("  2. double double precision;");
    put_line("  3. quad double precision.");
    put("Type 1, 2, or 3 to select the precision level : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Standard_Test;
      when '2' => DoblDobl_Test;
      when '3' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_soldrop;
