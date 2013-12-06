with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Symbol_Table,Symbol_Table_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Polynomial_Drops;

procedure ts_poldrop is

-- DESCRIPTION :
--   Test on dropping a variable from a polynomial.

  function Prompt_for_Drop_Variable return integer32 is

  -- DESCRIPTION :
  --   Returns the index of the coordinate to be dropped.
  --   If the symbol table is empty, then the user will be asked
  --   to enter a number, otherwise the user can enter the name
  --   of the variable for the drop.

    res : integer32 := 0;

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
          res := integer32(Symbol_Table.Get(sb));
          exit when (res /= 0);
          put("symbol not found, please try again...");
        end loop;
      end;
    end if;
    return res;
  end Prompt_for_Drop_Variable;

  procedure Standard_Test is

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;
    k : integer32 := 0;

  begin
    new_line;
    get(lp);
    new_line;
    k := Prompt_for_Drop_Variable;
    declare
      dropped : Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,k);
    begin
      Symbol_Table.Remove(natural32(k));
      put("The polynomial system with coordinate "); put(k,1);
      put_line(" dropped :");
      put(dropped);
    end;
  end Standard_Test;

  procedure DoblDobl_Test is

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;
    k : integer32 := 0;

  begin
    new_line;
    get(lp);
    new_line;
    k := Prompt_for_Drop_Variable;
    declare
      dropped : Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,k);
    begin
      Symbol_Table.Remove(natural32(k));
      put("The polynomial system with coordinate "); put(k,1);
      put_line(" dropped :");
      put(dropped);
    end;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;
    k : integer32 := 0;

  begin
    new_line;
    get(lp);
    new_line;
    k := Prompt_for_Drop_Variable;
    declare
      dropped : Poly_Sys(lp'range) := Polynomial_Drops.Drop(lp.all,k);
    begin
      Symbol_Table.Remove(natural32(k));
      put("The polynomial system with coordinate "); put(k,1);
      put_line(" dropped :");
      put(dropped);
    end;
  end QuadDobl_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for dropping a variable from a polynomial :");
    put_line("  1. standard double precision;");
    put_line("  2. double double precision;");
    put_line("  3. quad double precision;");
    put("Type 1, 2, or 3 to select precision : ");
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
end ts_poldrop;
