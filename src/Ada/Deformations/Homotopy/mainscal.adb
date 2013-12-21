with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Scaling;
with Drivers_for_Scaling;                use Drivers_for_Scaling;

procedure mainscal ( infilename,outfilename : in string ) is

  lp : Link_to_Poly_Sys;
  infile : file_type;
  outfile : file_type;
  basis : natural32;
  n : integer32;
  scalvec : Link_to_Vector;
  ans : character;
  sysonfile : boolean;

  procedure Read_System ( file : in out file_type; filename : in string ) is
  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      n := lp'length;
      sysonfile := true;
    else
      sysonfile := false;
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      sysonfile := false; lp := null; return;
  end Read_System;

  procedure Separate_File ( p : in Poly_Sys ) is

    scafile : file_type;
    nunk : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    new_line;
    put("Do you want the scaled system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(scafile);
      if nunk = natural32(p'length)
       then put(scafile,natural32(p'length),p);
       else put(scafile,natural32(p'length),nunk,p);
      end if;
      if basis /= 0 then
        new_line(scafile);
        put_line(scafile,"SCALING COEFFICIENTS :");
        new_line(scafile);
        put(scafile,basis,1); new_line(scafile);
        put_line(scafile,scalvec);
      end if;
      Close(scafile);
    end if;
  end Separate_File;

  procedure Rescale is
  
    sols : Solution_List;
    found : boolean;
    m : natural32;

  begin
    if sysonfile then                       -- scan for scaling coefficients
      Scan_and_Skip(infile,"SCALING COEFFICIENTS",found);
      if found then
        get(infile,basis);
        scalvec := new vector(1..2*n);
        get(infile,scalvec.all);
      end if;
    else 
      found := false;
    end if;
    if not found then
      put("Give the basis : "); get(basis);
      put("Give "); put(2*n,1); put_line(" complex scaling numbers : ");
      scalvec := new vector(1..2*n);
      for i in scalvec'range loop
        get(scalvec(i));
      end loop;
    end if;
    if sysonfile then                              -- scan for the solutions
      Reset(infile);
      Scan_and_Skip(infile,"SOLUTIONS",found);
      if found
       then get(infile,sols);
      end if;
      Close(infile);
    else
      found := false;
    end if;
    if not found then
      put_line("Reading the name of the file for the solutions.");
      Read_Name_and_Open_File(infile);
      get(infile,sols);
      Close(infile);
    end if;
    put_line(outfile,"THE SCALING COEFFICIENTS : ");
    new_line(outfile);
    put(outfile,basis,1); new_line(outfile);
    put_line(outfile,scalvec);
    new_line(outfile);
    Standard_Scaling.Scale(basis,scalvec.all,sols);
    m := Length_Of(sols);
    if m > 0 then
      put_line(outfile,"THE DE-SCALED SOLUTIONS : ");
      new_line(outfile);
      put(outfile,m,natural32(Head_Of(sols).n),sols);
    end if;
    Close(outfile);
  end Rescale;

  procedure Display_and_Dispatch_Menu
               ( file : in file_type; p : in out Poly_Sys ) is

  -- DESCRIPTION :
  --   Displays the menu and returns a choice, corresponding to one of the
  --   three available scaling procedures.

  begin
    loop
      new_line;
      put_line("MENU for Scaling Polynomial Systems :");
      put_line("  1 : Equation Scaling : divide by average coefficient      ");
      put_line("  2 : Variable Scaling : change of variables, as z = (2^c)*x");
      put_line("  3 : Solution Scaling : back to original coordinates       ");
      put("Type 1, 2, or 3 to select scaling, or i for info : ");
      Ask_Alternative(ans,"123i");
      if ans = 'i'
       then new_line; Drivers_for_Scaling.Display_Info; new_line;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(file,p); basis := 0;
      when '2' => Variable_Scaling(file,p,basis,scalvec);
      when '3' => Rescale;
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(file,p,basis,scalvec);
      when others    => null;
    end case;
    if ans /= '3'
     then Separate_File(p);
    end if;
  end Display_and_Dispatch_Menu;

begin
  Read_System(infile,infilename);
  if lp = null then
    loop
      new_line;
      put("Is the system on a file ? (y/n/i=info) ");
      Ask_Alternative(ans,"yni");
      if ans = 'i' then
        new_line;
        Standard_Complex_Poly_Systems_io.Display_Format;
        new_line;
      end if;
      exit when ans /= 'i';
    end loop;
    new_line;
    if ans = 'y' then
      put_line("Reading the name of the input file.");
      Read_Name_and_Open_File(infile);
      get(infile,lp);
      sysonfile := true;
      n := lp'length;
    else
      put("Give the dimension : "); get(n);
      lp := new Poly_Sys(1..n);
      put("Give "); put(n,1); put(" "); put(n,1);
      put_line("-variate polynomials :");
      get(natural32(n),lp.all);
      skip_line;  -- skip end_of_line symbol
      sysonfile := false;
    end if;
  end if;
  Create_Output_File(outfile,outfilename);
  put(outfile,natural32(lp'last),lp.all); new_line(outfile);
  Display_and_Dispatch_Menu(outfile,lp.all);
end mainscal;
