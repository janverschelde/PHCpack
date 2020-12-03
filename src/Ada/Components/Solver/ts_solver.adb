with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with String_Splitters;                   use String_Splitters;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Main_Eqn_by_Eqn_Solvers;            use Main_Eqn_by_Eqn_Solvers;

procedure ts_solver is

-- DESCRIPTION :
--   Interactive stand-alone tester of the equation-by-equation solver.

  procedure Main is

    lp : Link_to_Poly_Sys;
    file : file_type;
    name : Link_to_String;
    ans : character;

  begin
    new_line;
    put_line("Testing the equation-by-equation solver ...");
    get(lp);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file,name);
    new_line;
    put("Do you want to test the version with generics ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Solver_using_Generics(file,name.all,lp.all);
     else Shuffle_Polynomials_and_Solve(file,name.all,lp.all);
    end if;
  end Main;

begin
  Main;
end ts_solver;
