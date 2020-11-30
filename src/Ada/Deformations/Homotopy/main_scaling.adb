with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Scaling;
with DoblDobl_Scaling;
with QuadDobl_Scaling;
with Scaling_Methods;

package body Main_Scaling is

  procedure Standard_Read_System 
              ( file : in out file_type; filename : in string;
                dim : out integer32;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      dim := lp'length;
    else
      dim := 0;
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      dim := 0; lp := null; return;
  end Standard_Read_System;

  procedure DoblDobl_Read_System 
              ( file : in out file_type; filename : in string;
                dim : out integer32;
                lp : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      dim := lp'length;
    else
      dim := 0;
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      dim := 0; lp := null; return;
  end DoblDobl_Read_System;

  procedure QuadDobl_Read_System 
              ( file : in out file_type; filename : in string;
                dim : out integer32;
                lp : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      dim := lp'length;
    else
      dim := 0;
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      dim := 0; lp := null; return;
  end QuadDobl_Read_System;

  procedure Standard_Separate_File 
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scalvec : in Standard_Complex_Vectors.Link_to_Vector ) is

    ans : character;
    scafile : file_type;
    use Standard_Complex_Polynomials;
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
  end Standard_Separate_File;

  procedure DoblDobl_Separate_File 
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scalvec : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    ans : character;
    scafile : file_type;
   -- use DoblDobl_Complex_Polynomials;
   -- nunk : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    new_line;
    put("Do you want the scaled system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(scafile);
      put(scafile,p);
      if basis /= 0 then
        new_line(scafile);
        put_line(scafile,"SCALING COEFFICIENTS :");
        new_line(scafile);
        put(scafile,basis,1); new_line(scafile);
        put_line(scafile,scalvec);
      end if;
      Close(scafile);
    end if;
  end DoblDobl_Separate_File;

  procedure QuadDobl_Separate_File 
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                basis : in natural32;
                scalvec : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    ans : character;
    scafile : file_type;
   -- use QuadDobl_Complex_Polynomials;
   -- nunk : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    new_line;
    put("Do you want the scaled system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(scafile);
      put(scafile,p);
      if basis /= 0 then
        new_line(scafile);
        put_line(scafile,"SCALING COEFFICIENTS :");
        new_line(scafile);
        put(scafile,basis,1); new_line(scafile);
        put_line(scafile,scalvec);
      end if;
      Close(scafile);
    end if;
  end QuadDobl_Separate_File;

  procedure Standard_Rescale
                ( dim : in integer32;
                  infile,outfile : in out file_type;
                  sysonfile : in boolean ) is

    use Standard_Complex_Vectors;
    use Standard_Complex_Solutions;

    sols : Solution_List;
    found : boolean;
    basis : natural32 := 0;
    m : natural32;
    scalvec : Link_to_Vector;

  begin
    if sysonfile then                       -- scan for scaling coefficients
      Scan_and_Skip(infile,"SCALING COEFFICIENTS",found);
      if found then
        get(infile,basis);
        scalvec := new vector(1..2*dim);
        get(infile,scalvec.all);
      end if;
    else 
      found := false;
    end if;
    if not found then
      put("Give the basis : "); get(basis);
      put("Give "); put(2*dim,1); put_line(" complex scaling numbers : ");
      scalvec := new vector(1..2*dim);
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
      put_line(outfile,"THE RESCALED SOLUTIONS : ");
      new_line(outfile);
      put(outfile,m,natural32(Head_Of(sols).n),sols);
    end if;
    Close(outfile);
  end Standard_Rescale;

  procedure DoblDobl_Rescale
                ( dim : in integer32;
                  infile,outfile : in out file_type;
                  sysonfile : in boolean ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Solutions;

    sols : Solution_List;
    found : boolean;
    basis : natural32 := 0;
    m : natural32;
    scalvec : Link_to_Vector;

  begin
    if sysonfile then                       -- scan for scaling coefficients
      Scan_and_Skip(infile,"SCALING COEFFICIENTS",found);
      if found then
        get(infile,basis);
        scalvec := new vector(1..2*dim);
        get(infile,scalvec.all);
      end if;
    else 
      found := false;
    end if;
    if not found then
      put("Give the basis : "); get(basis);
      put("Give "); put(2*dim,1); put_line(" complex scaling numbers : ");
      scalvec := new vector(1..2*dim);
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
    DoblDobl_Scaling.Scale(basis,scalvec.all,sols);
    m := Length_Of(sols);
    if m > 0 then
      put_line(outfile,"THE RESCALED SOLUTIONS : ");
      new_line(outfile);
      put(outfile,m,natural32(Head_Of(sols).n),sols);
    end if;
    Close(outfile);
  end DoblDobl_Rescale;

  procedure QuadDobl_Rescale
                ( dim : in integer32;
                  infile,outfile : in out file_type;
                  sysonfile : in boolean ) is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Solutions;

    sols : Solution_List;
    found : boolean;
    basis : natural32 := 0;
    m : natural32;
    scalvec : Link_to_Vector;

  begin
    if sysonfile then                       -- scan for scaling coefficients
      Scan_and_Skip(infile,"SCALING COEFFICIENTS",found);
      if found then
        get(infile,basis);
        scalvec := new vector(1..2*dim);
        get(infile,scalvec.all);
      end if;
    else 
      found := false;
    end if;
    if not found then
      put("Give the basis : "); get(basis);
      put("Give "); put(2*dim,1); put_line(" complex scaling numbers : ");
      scalvec := new vector(1..2*dim);
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
    QuadDobl_Scaling.Scale(basis,scalvec.all,sols);
    m := Length_Of(sols);
    if m > 0 then
      put_line(outfile,"THE RESCALED SOLUTIONS : ");
      new_line(outfile);
      put(outfile,m,natural32(Head_Of(sols).n),sols);
    end if;
    Close(outfile);
  end QuadDobl_Rescale;

  procedure Standard_Display_and_Dispatch_Menu
               ( infile,outfile : in out file_type; dim : in integer32;
                 p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                 sysonfile : in boolean ) is

    ans : character;
    basis : natural32;
    scalvec : Standard_Complex_Vectors.Link_to_Vector;

    use Scaling_Methods;

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
       then new_line; Scaling_Methods.Display_Info; new_line;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(outfile,p); basis := 0;
      when '2' => Variable_Scaling(outfile,p,basis,scalvec);
      when '3' => Standard_Rescale(dim,infile,outfile,sysonfile);
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(outfile,p,basis,scalvec);
      when others    => null;
    end case;
    if ans /= '3'
     then Standard_Separate_File(p,basis,scalvec);
    end if;
  end Standard_Display_and_Dispatch_Menu;

  procedure DoblDobl_Display_and_Dispatch_Menu
               ( infile,outfile : in out file_type; dim : in integer32;
                 p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sysonfile : in boolean ) is

    ans : character;
    basis : natural32;
    scalvec : DoblDobl_Complex_Vectors.Link_to_Vector;

    use Scaling_Methods;

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
       then new_line; Scaling_Methods.Display_Info; new_line;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(outfile,p); basis := 0;
      when '2' => Variable_Scaling(outfile,p,basis,scalvec);
      when '3' => DoblDobl_Rescale(dim,infile,outfile,sysonfile);
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(outfile,p,basis,scalvec);
      when others    => null;
    end case;
    if ans /= '3'
     then DoblDobl_Separate_File(p,basis,scalvec);
    end if;
  end DoblDobl_Display_and_Dispatch_Menu;

  procedure QuadDobl_Display_and_Dispatch_Menu
               ( infile,outfile : in out file_type; dim : in integer32;
                 p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sysonfile : in boolean ) is

    ans : character;
    basis : natural32;
    scalvec : QuadDobl_Complex_Vectors.Link_to_Vector;

    use Scaling_Methods;

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
       then new_line; Scaling_Methods.Display_Info; new_line;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(outfile,p); basis := 0;
      when '2' => Variable_Scaling(outfile,p,basis,scalvec);
      when '3' => QuadDobl_Rescale(dim,infile,outfile,sysonfile);
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(outfile,p,basis,scalvec);
      when others    => null;
    end case;
    if ans /= '3'
     then QuadDobl_Separate_File(p,basis,scalvec);
    end if;
  end QuadDobl_Display_and_Dispatch_Menu;

  procedure Standard_Main ( infilename,outfilename : in string ) is

    ans : character;
    n : integer32;
    infile,outfile : file_type;
    sysonfile : boolean;
    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;
 
  begin
    Standard_Read_System(infile,infilename,n,lp);
    sysonfile := (n > 0);
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
    Standard_Display_and_Dispatch_Menu(infile,outfile,n,lp.all,sysonfile);
  end Standard_Main;

  procedure DoblDobl_Main ( infilename,outfilename : in string ) is

    ans : character;
    n : integer32;
    infile,outfile : file_type;
    sysonfile : boolean;
    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;
 
  begin
    DoblDobl_Read_System(infile,infilename,n,lp);
    sysonfile := (n > 0);
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
        get(lp);
        n := lp'last;
        sysonfile := false;
      end if;
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,lp.all); new_line(outfile);
    DoblDobl_Display_and_Dispatch_Menu(infile,outfile,n,lp.all,sysonfile);
  end DoblDobl_Main;

  procedure QuadDobl_Main ( infilename,outfilename : in string ) is

    ans : character;
    n : integer32;
    infile,outfile : file_type;
    sysonfile : boolean;
    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;
 
  begin
    QuadDobl_Read_System(infile,infilename,n,lp);
    sysonfile := (n > 0);
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
        get(lp);
        n := lp'last;
        sysonfile := false;
      end if;
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,lp.all); new_line(outfile);
    QuadDobl_Display_and_Dispatch_Menu(infile,outfile,n,lp.all,sysonfile);
  end QuadDobl_Main;

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    ans : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in mainscal.Main ...");
    end if;
    new_line;
    put_line("MENU for the precision of the scalers :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main(infilename,outfilename);
      when '1' => DoblDobl_Main(infilename,outfilename);
      when '2' => QuadDobl_Main(infilename,outfilename);
      when others => null;
    end case;
  end Main;

end Main_Scaling;
