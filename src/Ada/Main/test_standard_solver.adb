with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Black_Box_Solvers;
with Black_Box_Polyhedral_Solvers;

package body Test_Standard_Solver is

  procedure Solve_to_File
              ( nt : in natural32; mvonly : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    outfile : file_type;
    rc : natural32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    put(outfile,natural32(p'last),p);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if mvonly then
      if nt = 0
       then Black_Box_Polyhedral_Solvers.Solve(outfile,p,true,rc,sols,vrb);
       else Black_Box_Polyhedral_Solvers.Solve(outfile,nt,p,true,rc,sols,vrb);
      end if;
    else
      if nt = 0
       then Black_Box_Solvers.Solve(outfile,p,true,rc,sols,vrb);
       else Black_Box_Solvers.Solve(outfile,nt,p,true,rc,sols,vrb);
      end if;
    end if;
  end Solve_to_File;

  procedure Solve_without_File
              ( nt : in natural32; mvonly : in boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                vrb : in integer32 := 0 ) is

    rc : natural32 := 0;
    sols : Solution_List;

  begin
    if mvonly then
      if nt = 0
       then Black_Box_Polyhedral_Solvers.Solve(p,false,true,rc,sols,vrb);
       else Black_Box_Polyhedral_Solvers.Solve(nt,p,false,true,rc,sols,vrb);
      end if;
    else
      if nt = 0 
       then Black_Box_Solvers.Solve(p,false,true,rc,sols,vrb);
       else Black_Box_Solvers.Solve(nt,p,false,true,rc,sols,vrb);
      end if;
    end if;
    new_line;
    put("The root count : "); put(rc,1); new_line;
    if not Is_Null(sols) then
      new_line;
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Solve_without_File;

  procedure Main is

    lp : Link_to_Poly_Sys;
    nt : natural32 := 0;
    vrb : integer32 := 0;
    mvfocus : boolean;
    ans : character;

  begin
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    new_line;
    put_line("-> your system : "); put(lp.all);
    new_line;
    put("Give the verbose level : "); get(vrb);
    put("Give the number of tasks (0 for no multitasking) : "); get(nt);
    put("Focus on mixed volumes and polyhedral homotopies ? (y/n) ");
    Ask_Yes_or_No(ans); mvfocus := (ans = 'y');
    put("Write the output to file ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Solve_to_File(nt,mvfocus,lp.all,vrb);
     else Solve_without_File(nt,mvfocus,lp.all,vrb);
    end if;
  end Main;

end Test_Standard_Solver;
