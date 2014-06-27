with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Standard_Simpomial_Solvers;        use Standard_Simpomial_Solvers;
with DoblDobl_Simpomial_Solvers;        use DoblDobl_Simpomial_Solvers;
with QuadDobl_Simpomial_Solvers;        use QuadDobl_Simpomial_Solvers;

procedure ts_simposol is

-- DESCRIPTION :
--   Tests on the wrappers to parse and solve systems that can be
--   reduced to binomial systems.

  procedure Standard_Simpomial_Solver is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and
  --   tries to solve it, using simplex solvers
  --   in standard double complex arithmetic.

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    fail,zero_y : boolean;
    ans : character;
    file : file_type;
    tol : constant double_float := 1.0E-12;
    rcond : double_float;

  begin
    get(lp);
    new_line;
    if not Is_Simplex_System(lp.all) then
      put_line("The system is a not simplex system.");
    else
      put_line("The system is a simplex system.");
      Solve(lp.all,tol,rcond,sols,fail,zero_y);
      if fail then
        put_line("Simplex solver returned failure!");
      else
        put("Estimate for inverse condition number = ");
        put(rcond,3); new_line;
        if zero_y then
          put_line("No solutions with all components different from zero.");
        else
          put("Found "); put(Length_Of(sols),1); put_line(" solutions.");
          put("Do you want the solutions on file ? (y/n) ");          
          Ask_Yes_or_No(ans);
          if ans = 'y' then
            new_line;
            put_line("Reading the name of the output file.");
            Read_Name_and_Create_File(file);
            put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
            close(file);
          end if;
        end if;
      end if;
    end if;
  end Standard_Simpomial_Solver;

  procedure DoblDobl_Simpomial_Solver is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and
  --   tries to solve it, using simplex solvers
  --   in double double complex arithmetic.

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    fail,zero_y : boolean;
    ans : character;
    file : file_type;
    tol : constant double_double := create(1.0E-12);
    rcond : double_double;

  begin
    get(lp);
    new_line;
    if not Is_Simplex_System(lp.all) then
      put_line("The system is a not simplex system.");
    else
      put_line("The system is a simplex system.");
      Solve(lp.all,tol,rcond,sols,fail,zero_y);
      if fail then
        put_line("Simplex solver returned failure!");
      else
        put("Estimate for inverse condition number = ");
        put(rcond,3); new_line;
        if zero_y then
          put_line("No solutions with all components different from zero.");
        else
          put("Found "); put(Length_Of(sols),1); put_line(" solutions.");
          put("Do you want the solutions on file ? (y/n) ");          
          Ask_Yes_or_No(ans);
          if ans = 'y' then
            new_line;
            put_line("Reading the name of the output file.");
            Read_Name_and_Create_File(file);
            put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
            close(file);
          end if;
        end if;
      end if;
    end if;
  end DoblDobl_Simpomial_Solver;

  procedure QuadDobl_Simpomial_Solver is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and
  --   tries to solve it, using simplex solvers
  --   in double double complex arithmetic.

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    fail,zero_y : boolean;
    ans : character;
    file : file_type;
    tol : constant quad_double := create(1.0E-12);
    rcond : quad_double;

  begin
    get(lp);
    new_line;
    if not Is_Simplex_System(lp.all) then
      put_line("The system is a not simplex system.");
    else
      put_line("The system is a simplex system.");
      Solve(lp.all,tol,rcond,sols,fail,zero_y);
      if fail then
        put_line("Simplex solver returned failure!");
      else
        put("Estimate for inverse condition number = ");
        put(rcond,3); new_line;
        if zero_y then
          put_line("No solutions with all components different from zero.");
        else
          put("Found "); put(Length_Of(sols),1); put_line(" solutions.");
          put("Do you want the solutions on file ? (y/n) ");          
          Ask_Yes_or_No(ans);
          if ans = 'y' then
            new_line;
            put_line("Reading the name of the output file.");
            Read_Name_and_Create_File(file);
            put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
            close(file);
          end if;
        end if;
      end if;
    end if;
  end QuadDobl_Simpomial_Solver;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and
  --   tries to solve it, using simplex solvers.

    ans : character;

  begin
    new_line;
    put_line("MENU for the solving a simplex polynomial systems ...");
    put_line("  1. in standard double complex arithmetic;");
    put_line("  2. in double double complex arithmetic;");
    put_line("  3. in quad double complex arithmetic.");
    put("Type 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Standard_Simpomial_Solver;
      when '2' => DoblDobl_Simpomial_Solver;
      when '3' => QuadDobl_Simpomial_Solver;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_simposol;
