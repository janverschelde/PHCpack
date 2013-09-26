with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Simpomial_Solvers;        use Standard_Simpomial_Solvers;

procedure ts_simposol is

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system and
  --   tries to solve it, using simplex solvers.

    lp : Link_to_Poly_Sys;
    sols : Solution_List;
    fail,zero_y : boolean;
    ans : character;
    file : file_type;
    tol : constant double_float := 1.0E-12;
    rcond : double_float;

  begin
    new_line;
    put_line("Testing the solution of simplex polynomial systems ...");
    new_line;
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
  end Main;

begin
  Main;
end ts_simposol;
