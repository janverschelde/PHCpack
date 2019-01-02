with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;            use DoblDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;            use QuadDobl_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;   use Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;   use DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;   use QuadDobl_System_and_Solutions_io;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_BlackBox_Continuations;    use DoblDobl_BlackBox_Continuations;
with QuadDobl_BlackBox_Continuations;    use QuadDobl_BlackBox_Continuations;

procedure ts_bbpoco is

-- DESCRIPTION :
--   Calls the blackbox procedures for polynomial continuation.

  procedure Standard_Continuation is

  -- DESCRIPTION :
  --   Prompts the user for a target and start system
  --   and then calls the black box polynomal continuation
  --   in standard double precision.

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    file : file_type;
    lp,lq : Link_to_Poly_Sys;
    sols : Solution_List;
    gamma : constant Complex_Number := Random1;
    nt : integer32 := 0;
    elap : duration;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(lp);
    new_line;
    put_line("Reading the start system and start solutions ...");
    get(lq,sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Black_Box_Polynomial_Continuation
      (file,true,nt,lp.all,lq.all,gamma,sols,elap);
  end Standard_Continuation;

  procedure DoblDobl_Continuation is

  -- DESCRIPTION :
  --   Prompts the user for a target and start system
  --   and then calls the black box polynomal continuation
  --   in double double precision.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    file : file_type;
    lp,lq : Link_to_Poly_Sys;
    sols : Solution_List;
    gamma : constant Complex_Number := Random1;
    nt : integer32 := 0;
    elap : duration;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(lp);
    new_line;
    put_line("Reading the start system and start solutions ...");
    get(lq,sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Black_Box_Polynomial_Continuation(file,nt,lp.all,lq.all,gamma,sols,elap);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end DoblDobl_Continuation;

  procedure QuadDobl_Continuation is

  -- DESCRIPTION :
  --   Prompts the user for a target and start system
  --   and then calls the black box polynomal continuation
  --   in quad double precision.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    file : file_type;
    lp,lq : Link_to_Poly_Sys;
    sols : Solution_List;
    gamma : constant Complex_Number := Random1;
    nt : integer32 := 0;
    elap : duration;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(lp);
    new_line;
    put_line("Reading the start system and start solutions ...");
    get(lq,sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of tasks : "); get(nt);
    Black_Box_Polynomial_Continuation(file,nt,lp.all,lq.all,gamma,sols,elap);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end QuadDobl_Continuation;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision.

    ans : character;

  begin
    new_line;
    put_line("MENU to call the blackbox continuation procedures :");
    put_line("  1. use standard double float arithmetic;");
    put_line("  2. use double double float arithmetic;");
    put_line("  3. use quad double float arithmetic;");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Standard_Continuation;
      when '2' => DoblDobl_Continuation;
      when '3' => QuadDobl_Continuation;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_bbpoco;
