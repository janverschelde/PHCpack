with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io; 
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;

procedure ts_sols is

-- DESCRIPTION :
--   Test on input/output of solution lists.

  procedure Test_Plain_Read is

    use Standard_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    sols : Solution_List;

  begin
    new_line; 
    put_line("Reading the polynomial system.");
    get(lp);
    new_line;
    put_line("The polynomial system :");
    put(lp.all);
    new_line;
    Read(sols);
    new_line;
    put_line("The solutions : ");
    put(sols);
  end Test_Plain_Read;

  procedure Test_Interactive_Standard_Read_Write is

    use Standard_Complex_Solutions;

    solfile,outfile : file_type;
    sols,sols_last : Solution_List;
    n,len,dim,nbr,cnt : natural32;
    first : boolean := true;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the solutions file ...");
    Read_Name_and_Open_File(solfile);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Give the number of solutions you want to read : ");
    Read_Natural(n);
    Read_First(solfile,n,len,dim,nbr,sols,sols_last);
    put("The length of the list : "); put(len,1); new_line;
    put("The dimension of the vectors : "); put(dim,1); new_line;
    loop
      put("The number of solutions read : "); put(nbr,1); new_line;
      if Length_Of(sols) > 0 then
        new_line;
        put("Do you want to see the solutions ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then put(standard_output,Length_Of(sols),dim,sols);
        end if;
        new_line;
        put("Give the number of solutions to be written to file : ");
        Read_Natural(n);
        if n > 0 then
          if first then
            Write_First(outfile,n,len,dim,nbr,sols);
            first := false; cnt := nbr;
          else
            Write_Next(outfile,n,nbr,cnt,sols);
            put("Updated counter to "); put(cnt,1); put_line(".");
          end if;
          put("Written "); put(nbr,1); put_line(" solutions to file.");
        end if;
        new_line;
        put("Do you want to clear the solutions ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then Clear(sols); sols_last := sols;
        end if;
      end if;
      new_line;
      put("Give the number of solutions to be read (0 to exit) : ");
      Read_Natural(n);
      exit when (n = 0);
      Read_Next(solfile,n,dim,nbr,sols,sols_last);
    end loop;
  end Test_Interactive_Standard_Read_Write;

  procedure Test_Interactive_Multprec_Read is

    use Multprec_Complex_Solutions;

    solfile,outfile : file_type;
    sols,sols_last : Solution_List;
    n,len,dim,nbr : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the solutions file ...");
    Read_Name_and_Open_File(solfile);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Give the number of solutions you want to read : "); get(n);
    Read_First(solfile,n,len,dim,nbr,sols,sols_last);
    put("The length of the list : "); put(len,1); new_line;
    put("The dimension of the vectors : "); put(dim,1); new_line;
    loop
      put("The number of solutions read : "); put(nbr,1); new_line;
      if Length_Of(sols) > 0 then
        new_line;
        put("Do you wish to see the solutions ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then put(standard_output,Length_Of(sols),dim,sols);
        end if;
        new_line;
        put("Do you wish to clear the solutions ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y'
         then Clear(sols); sols_last := sols;
        end if;
      end if;
      new_line;
      put("Give the number of solutions to be read (0 to exit) : ");
      get(n);
      exit when (n = 0);
      Read_Next(solfile,n,dim,nbr,sols,sols_last);
    end loop;
  end Test_Interactive_Multprec_Read;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test the i/o of solution lists :");
    put_line("  1. plain reading of a solution list from file;");
    put_line("  2. incremental read/write of a standard solution list;");
    put_line("  3. incremental read/write of a multprec solution list;");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Test_Plain_Read;
      when '2' => Test_Interactive_Standard_Read_Write;
      when '3' => Test_Interactive_Multprec_Read;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sols;
