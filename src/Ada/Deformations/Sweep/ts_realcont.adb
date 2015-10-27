with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Standard_Floating_Poly_Systems;
with Standard_Floating_Poly_Systems_io;  use Standard_Floating_Poly_Systems_io;
with Double_Double_Poly_Systems;
with Double_Double_Poly_Systems_io;      use Double_Double_Poly_Systems_io;
with Quad_Double_Poly_Systems;
with Quad_Double_Poly_Systems_io;        use Quad_Double_Poly_Systems_io;
with Standard_Complex_to_Real_Poly;      use Standard_Complex_to_Real_Poly;
with DoblDobl_Complex_to_Real_Poly;      use DoblDobl_Complex_to_Real_Poly;
with QuadDobl_Complex_to_Real_Poly;      use QuadDobl_Complex_to_Real_Poly;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Quad_Parameters;
with Standard_Quad_Sweepers;             use Standard_Quad_Sweepers;
with DoblDobl_Quad_Parameters;
with DoblDobl_Quad_Sweepers;             use DoblDobl_Quad_Sweepers;
with QuadDobl_Quad_Parameters;
with QuadDobl_Quad_Sweepers;             use QuadDobl_Quad_Sweepers;

procedure ts_realcont is

-- DESCRIPTION :
--   Real path following towards a quadratic bifurcation point.

  procedure Real_Sweep
              ( interactive : in boolean; nq,nv : in natural32; 
                p : in Standard_Floating_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Given a system with real coefficients, a real sweep starts,
  --   in standard double precision.

    use Standard_Complex_Solutions;

    eigval : boolean;
    solfile : file_type;
    repsols : Solution_List;
    ans : character;

  begin
    if interactive then
      new_line;
      put("Do you want to compute eigenvalues along the path ? (y/n) ");
      Ask_Yes_or_No(ans);
      eigval := (ans = 'y');
      new_line;
      put("Do you want to record solutions along the path ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the file name to write solutions along path...");
        Read_Name_and_Create_File(solfile);
        new_line;
        Standard_Quad_Parameters.Tune;
        Interactive_Real_Sweep(nq,nv,eigval,p,repsols);
        put(solfile,Length_Of(repsols),natural32(Head_Of(repsols).n),repsols);
        close(solfile);
      else
        new_line;
        Standard_Quad_Parameters.Tune;
        Interactive_Real_Sweep(nq,nv,eigval,p,sols);
      end if;
    else
      Run_Sweep(nq,nv,p,sols);
    end if;
  end Real_Sweep;

  procedure Real_Sweep
              ( interactive : in boolean; nq,nv : in natural32; 
                p : in Double_Double_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Given a system with real coefficients, a real sweep starts,
  --   in double double precision.

    use DoblDobl_Complex_Solutions;

    eigval : boolean;
    solfile : file_type;
    repsols : Solution_List;
    ans : character;

  begin
    if interactive then
      new_line;
      put("Do you want to compute eigenvalues along the path ? (y/n) ");
      Ask_Yes_or_No(ans);
      eigval := (ans = 'y');
      new_line;
      put("Do you want to record solutions along the path ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the file name to write solutions along path...");
        Read_Name_and_Create_File(solfile);
        new_line;
        Standard_Quad_Parameters.Tune;
        Interactive_Real_Sweep(nq,nv,eigval,p,repsols);
        put(solfile,Length_Of(repsols),natural32(Head_Of(repsols).n),repsols);
        close(solfile);
      else
        new_line;
        Standard_Quad_Parameters.Tune;
        Interactive_Real_Sweep(nq,nv,eigval,p,sols);
      end if;
    else
      Run_Sweep(nq,nv,p,sols);
    end if;
  end Real_Sweep;

  procedure Real_Sweep
              ( interactive : in boolean; nq,nv : in natural32; 
                p : in Quad_Double_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Given a system with real coefficients, a real sweep starts,
  --   in quad double precision.

    use QuadDobl_Complex_Solutions;

    eigval : boolean;
    solfile : file_type;
    repsols : Solution_List;
    ans : character;

  begin
    if interactive then
      new_line;
      put("Do you want to compute eigenvalues along the path ? (y/n) ");
      Ask_Yes_or_No(ans);
      eigval := (ans = 'y');
      new_line;
      put("Do you want to record solutions along the path ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the file name to write solutions along path...");
        Read_Name_and_Create_File(solfile);
        new_line;
        Standard_Quad_Parameters.Tune;
        Interactive_Real_Sweep(nq,nv,eigval,p,repsols);
        put(solfile,Length_Of(repsols),natural32(Head_Of(repsols).n),repsols);
        close(solfile);
      else
        new_line;
        Standard_Quad_Parameters.Tune;
        Interactive_Real_Sweep(nq,nv,eigval,p,sols);
      end if;
    else
      Run_Sweep(nq,nv,p,sols);
    end if;
  end Real_Sweep;

  procedure Complex_Sweep
              ( interactive : in boolean; nq,nv : in natural32; 
                q : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The system on input has some coefficients with nonzero
  --   imaginary parts for which a complex sweep will start,
  --   in standard double precision.

  begin
    if interactive then
      new_line;
      Standard_Quad_Parameters.Tune;
      Interactive_Complex_Sweep(nq,nv,q,sols);
    else
      Run_Sweep(nq,nv,q,sols);
    end if;
  end Complex_Sweep;

  procedure Complex_Sweep
              ( interactive : in boolean; nq,nv : in natural32; 
                q : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The system on input has some coefficients with nonzero
  --   imaginary parts for which a complex sweep will start,
  --   in double double precision.

  begin
    if interactive then
      new_line;
      DoblDobl_Quad_Parameters.Tune;
      Interactive_Complex_Sweep(nq,nv,q,sols);
    else
      Run_Sweep(nq,nv,q,sols);
    end if;
  end Complex_Sweep;

  procedure Complex_Sweep
              ( interactive : in boolean; nq,nv : in natural32; 
                q : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The system on input has some coefficients with nonzero
  --   imaginary parts for which a complex sweep will start,
  --   in quad double precision.

  begin
    if interactive then
      new_line;
      QuadDobl_Quad_Parameters.Tune;
      Interactive_Complex_Sweep(nq,nv,q,sols);
    else
      Run_Sweep(nq,nv,q,sols);
    end if;
  end Complex_Sweep;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Setups the data for a sweep in standard double precision.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    file : file_type;
    p : Standard_Floating_Poly_Systems.Link_to_Poly_Sys;
    q : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    nq,nv : natural32;
    ans : character;
    found,realcff,interactive : boolean;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Parameter continuation on polynomial systems.");
    new_line;
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(file);
    get(file,q);
    new_line;
    nq := natural32(q'last);
    nv := Number_of_Unknowns(q(q'first));
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    realcff := Is_Real(q.all);
    new_line;
    if realcff then
      p := Convert_Complex_to_Real(q);
      put_line("Your system : "); put(p.all);
      put_line("has all its coefficients real.");
    else
      put_line("Your system : "); put(q.all);
      put_line("has some coefficients with nonzero imaginary parts.");
    end if;
    new_line;
    put("Continue (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      File_Scanning.Scan_and_Skip(file,"SOLUTIONS",found);
      if found then
        get(file,sols);
        put_line("The solutions : ");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
      new_line;
      put("Interactive version for one solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      interactive := (ans = 'y');
      if realcff
       then Real_Sweep(interactive,nq,nv,p.all,sols);
       else Complex_Sweep(interactive,nq,nv,q.all,sols);
      end if;
    end if;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Setups the data for a sweep in double double precision.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;

    file : file_type;
    p : Double_Double_Poly_Systems.Link_to_Poly_Sys;
    q : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nq,nv : natural32;
    ans : character;
    found,realcff,interactive : boolean;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Parameter continuation on polynomial systems.");
    new_line;
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(file);
    get(file,q);
    new_line;
    nq := natural32(q'last);
    nv := Number_of_Unknowns(q(q'first));
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    realcff := Is_Real(q.all);
    new_line;
    if realcff then
      p := Convert_Complex_to_Real(q);
      put_line("Your system : "); put(p.all);
      put_line("has all its coefficients real.");
    else
      put_line("Your system : "); put(q.all);
      put_line("has some coefficients with nonzero imaginary parts.");
    end if;
    new_line;
    put("Continue (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      File_Scanning.Scan_and_Skip(file,"SOLUTIONS",found);
      if found then
        get(file,sols);
        put_line("The solutions : ");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
      new_line;
      put("Interactive version for one solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      interactive := (ans = 'y');
      if realcff
       then Real_Sweep(interactive,nq,nv,p.all,sols);
       else Complex_Sweep(interactive,nq,nv,q.all,sols);
      end if;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Setups the data for a sweep in double double precision.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;

    file : file_type;
    p : Quad_Double_Poly_Systems.Link_to_Poly_Sys;
    q : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nq,nv : natural32;
    ans : character;
    found,realcff,interactive : boolean;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("Parameter continuation on polynomial systems.");
    new_line;
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(file);
    get(file,q);
    new_line;
    nq := natural32(q'last);
    nv := Number_of_Unknowns(q(q'first));
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    realcff := Is_Real(q.all);
    new_line;
    if realcff then
      p := Convert_Complex_to_Real(q);
      put_line("Your system : "); put(p.all);
      put_line("has all its coefficients real.");
    else
      put_line("Your system : "); put(q.all);
      put_line("has some coefficients with nonzero imaginary parts.");
    end if;
    new_line;
    put("Continue (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      File_Scanning.Scan_and_Skip(file,"SOLUTIONS",found);
      if found then
        get(file,sols);
        put_line("The solutions : ");
        put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
      new_line;
      put("Interactive version for one solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      interactive := (ans = 'y');
      if realcff
       then Real_Sweep(interactive,nq,nv,p.all,sols);
       else Complex_Sweep(interactive,nq,nv,q.all,sols);
      end if;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then calls the corresponding test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_realcont;
