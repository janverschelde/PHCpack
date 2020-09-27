with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Black_Box_Solvers;

procedure ts_bbsolve is

-- DESCRIPTION :
--   Tests the black box solver.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for a polynomial systems and calls the black box solver.

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    rc : natural32;
    vrb : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("-> your system : "); put(lp.all);
    new_line;
    put("Give the verbose level : "); get(vrb);
    new_line;
    put("Run only polyhedral homotopies ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Black_Box_Solvers.Polyhedral_Solve(lp.all,false,true,rc,sols,vrb);
     else Black_Box_Solvers.Solve(lp.all,false,true,rc,sols,vrb);
    end if;
    new_line;
    put("The root count : "); put(rc,1); new_line;
    if not Is_Null(sols) then
      new_line;
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Main;

begin
  Main;
end ts_bbsolve;
