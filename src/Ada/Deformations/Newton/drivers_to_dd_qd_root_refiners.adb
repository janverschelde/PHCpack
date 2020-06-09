with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Laur_Systems;
with Multprec_Complex_Laurentials_io;
with Multprec_Complex_Laur_Systems;
with DoblDobl_Polynomial_Convertors;     use DoblDobl_Polynomial_Convertors;
with DoblDobl_Complex_Laurentials;       use DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Laurentials;       use QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Polynomial_Convertors;     use QuadDobl_Polynomial_Convertors;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Multprec_System_and_Solutions_io;
with Root_Refining_Parameters;           use Root_Refining_Parameters;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;             use QuadDobl_Root_Refiners;
with Multitasked_DD_QD_Refiners;         use Multitasked_DD_QD_Refiners;

package body Drivers_to_dd_qd_Root_Refiners is

  procedure Standard_to_DoblDobl_Complex
              ( p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out DoblDobl_Complex_Solutions.Solution_List ) is

    q : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Standard_System_and_Solutions_io.get(q,sols);
    declare
      sq : constant DoblDobl_Complex_Laur_Systems.Laur_Sys
         := Standard_Laur_Sys_to_DoblDobl_Complex(q.all);
    begin
      p := new DoblDobl_Complex_Laur_Systems.Laur_Sys'(sq);
    end;
    s := DoblDobl_Complex_Solutions.Create(sols);
  end Standard_to_DoblDobl_Complex;

  procedure Multprec_to_DoblDobl_Complex
              ( p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out DoblDobl_Complex_Solutions.Solution_List ) is

    q : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    Multprec_Complex_Laurentials_io.Set_Working_Precision(5);
    Multprec_System_and_Solutions_io.get(q,sols);
    declare
      sq : constant DoblDobl_Complex_Laur_Systems.Laur_Sys
         := Multprec_Laur_Sys_to_DoblDobl_Complex(q.all);
    begin
      p := new DoblDobl_Complex_Laur_Systems.Laur_Sys'(sq);
    end;
    s := DoblDobl_Complex_Solutions.Create(sols);
  end Multprec_to_DoblDobl_Complex;

  procedure Multprec_to_DoblDobl_Complex
              ( file : in file_type;
                p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out DoblDobl_Complex_Solutions.Solution_List ) is

    q : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    Multprec_Complex_Laurentials_io.Set_Working_Precision(5);
    Multprec_System_and_Solutions_io.get(file,q,sols);
    declare
      sq : constant DoblDobl_Complex_Laur_Systems.Laur_Sys
         := Multprec_Laur_Sys_to_DoblDobl_Complex(q.all);
    begin
      p := new DoblDobl_Complex_Laur_Systems.Laur_Sys'(sq);
    end;
    s := DoblDobl_Complex_Solutions.Create(sols);
  end Multprec_to_DoblDobl_Complex;

  procedure Standard_to_QuadDobl_Complex
              ( p : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out QuadDobl_Complex_Solutions.Solution_List ) is

    q : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Standard_System_and_Solutions_io.get(q,sols);
    declare
      sq : constant QuadDobl_Complex_Laur_Systems.Laur_Sys
         := Standard_Laur_Sys_to_QuadDobl_Complex(q.all);
    begin
      p := new QuadDobl_Complex_Laur_Systems.Laur_Sys'(sq);
    end;
    s := QuadDobl_Complex_Solutions.Create(sols);
  end Standard_to_QuadDobl_Complex;

  procedure Multprec_to_QuadDobl_Complex
              ( p : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out QuadDobl_Complex_Solutions.Solution_List ) is

    q : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    Multprec_Complex_Laurentials_io.Set_Working_Precision(10);
    Multprec_System_and_Solutions_io.get(q,sols);
    declare
      sq : constant QuadDobl_Complex_Laur_Systems.Laur_Sys
         := Multprec_Laur_Sys_to_QuadDobl_Complex(q.all);
    begin
      p := new QuadDobl_Complex_Laur_Systems.Laur_Sys'(sq);
    end;
    s := QuadDobl_Complex_Solutions.Create(sols);
  end Multprec_to_QuadDobl_Complex;

  procedure Multprec_to_QuadDobl_Complex
              ( file : in file_type;
                p : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                s : out QuadDobl_Complex_Solutions.Solution_List ) is

    q : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    Multprec_Complex_Laurentials_io.Set_Working_Precision(10);
    Multprec_System_and_Solutions_io.get(file,q,sols);
    declare
      sq : constant QuadDobl_Complex_Laur_Systems.Laur_Sys
         := Multprec_Laur_Sys_to_QuadDobl_Complex(q.all);
    begin
      p := new QuadDobl_Complex_Laur_Systems.Laur_Sys'(sq);
    end;
    s := QuadDobl_Complex_Solutions.Create(sols);
  end Multprec_to_QuadDobl_Complex;

  procedure Prompt_for_Tasks
             ( genuine : in boolean; nbeq,nvar : in natural32;
               nbtasks : out integer32; otp : out boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the number of tasks,
  --   but only if not genuine and nbeq = nvar.

  -- ON ENTRY :
  --   genuine   true if the system is a genuine Laurent system;
  --   nbeq      number of equations;
  --   nvar      number of variables.

  -- ON RETURN :
  --   nbtasks   the number of tasks;
  --   otp       flag to indicate if output is wanted.

    ans : character;

  begin
    nbtasks := 0;
    if (not genuine) and (nbeq = nvar) then
      new_line;
      put("Apply multitasking ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("-> give the number of tasks : "); get(nbtasks);
        put("Monitor progress of tasks ? (y/n) "); 
        Ask_Yes_or_No(ans); otp := (ans = 'y');
      end if;
    end if;
  end Prompt_for_Tasks;

  procedure DD_Root_Refinement ( infilename,outfilename : in string ) is

  -- DESCRIPTION :
  --   Driver dedicated to refining in double double precision.

  -- ON ENTRY :
  --   infilename    the name of the input file, could be empty;
  --   outfilename   the name of the output file, could be empty.

    timer : Timing_Widget;
    infile,outfile : file_type;
    dd_p : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    dd_s : DoblDobl_Complex_Solutions.Solution_List;
    nbeq,nvar : natural32;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout,genuine,otp : boolean;
    nbtasks : integer32 := 0;

  begin
    if infilename /= "" then
      Open_Input_File(infile,infilename);
    else
      new_line;
      put_line("Reading the name of the input file ...");
      Read_Name_and_Open_File(infile);
    end if;
    Multprec_to_DoblDobl_Complex(infile,dd_p,dd_s);
    genuine := DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(dd_p.all);
    Create_Output_File(outfile,outfilename);
    nvar := Number_of_Unknowns(dd_p(dd_p'first));
    nbeq := natural32(dd_p'last);
    if nbeq = nvar
     then put(outfile,nbeq,dd_p.all);
     else put(outfile,nbeq,nvar,dd_p.all);
    end if;
    DoblDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (outfile,epsxa,epsfa,tolsing,maxit,deflate,wout);
    Prompt_for_Tasks(genuine,nbeq,nvar,nbtasks,otp);
    new_line; put("Refining "); 
    put(DoblDobl_Complex_Solutions.Length_Of(dd_s),1);
    put(" solutions ...");
    put_line(" see the output file for results."); new_line;
    if genuine then
      tstart(timer);
     -- DoblDobl_Root_Refiner(dd_p.all,dd_s);
      Reporting_Root_Refiner -- no deflate for Laurent systems yet
        (outfile,dd_p.all,dd_s,epsxa,epsfa,tolsing,numit,maxit,wout);
      tstop(timer);
    else
      declare
        use DoblDobl_Laur_Poly_Convertors;         
        q : DoblDobl_Complex_Poly_Systems.Poly_Sys(dd_p'range)
          := Laurent_to_Polynomial_System(dd_p.all);
      begin
        tstart(timer);
       -- DoblDobl_Root_Refiner(q,dd_s);
        if nbtasks = 0 then
          Reporting_Root_Refiner
            (outfile,q,dd_s,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        else
          Multitasking_Refinement(q,dd_s,nbtasks,otp,epsxa,epsfa,maxit);
          new_line(outfile);
          put_line(outfile,"THE SOLUTIONS :");
          DoblDobl_Complex_Solutions_io.write(outfile,dd_s);
        end if;
        tstop(timer);
        DoblDobl_Complex_Poly_Systems.Clear(q);
      end;
    end if;
   -- DoblDobl_Complex_Solutions_io.write(outfile,dd_s);
    new_line(outfile);
    print_times(outfile,timer,"double double Newton refinement");
    close(outfile);
  end DD_Root_Refinement;

  procedure QD_Root_Refinement ( infilename,outfilename : in string ) is

  -- DESCRIPTION :
  --   Driver dedicated to root refinement in quad double precision.

  -- ON ENTRY :
  --   infilename    the name of the input file, could be empty;
  --   outfilename   the name of the output file, could be empty.

    timer : Timing_Widget;
    infile,outfile : file_type;
    qd_p : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    qd_s : QuadDobl_Complex_Solutions.Solution_List;
    nbeq,nvar : natural32;
    epsxa,epsfa,tolsing : double_float;
    numit,maxit : natural32 := 0;
    deflate,wout,genuine,otp : boolean;
    nbtasks : integer32 := 0;

  begin
    if infilename /= "" then
      Open_Input_File(infile,infilename);
    else
      new_line;
      put_line("Reading the name of the input file ...");
      Read_Name_and_Open_File(infile);
    end if;
    Multprec_to_QuadDobl_Complex(infile,qd_p,qd_s);
    genuine := QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(qd_p.all);
    Create_Output_File(outfile,outfilename);
    nvar := Number_of_Unknowns(qd_p(qd_p'first));
    nbeq := natural32(qd_p'last);
    if nbeq = nvar
     then put(outfile,nbeq,qd_p.all);
     else put(outfile,nbeq,nvar,qd_p.all);
    end if;
    QuadDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (outfile,epsxa,epsfa,tolsing,maxit,deflate,wout);
    Prompt_for_Tasks(genuine,nbeq,nvar,nbtasks,otp);
    new_line; put("Refining "); 
    put(QuadDobl_Complex_Solutions.Length_Of(qd_s),1);
    put(" solutions ...");
    put_line(" see the output file for results."); new_line;
    if genuine then
      tstart(timer);
     -- QuadDobl_Root_Refiner(qd_p.all,qd_s);
      Reporting_Root_Refiner  -- no deflation for Laurent systems yet
        (outfile,qd_p.all,qd_s,epsxa,epsfa,tolsing,numit,maxit,wout);
      tstop(timer);
    else
      declare
        use QuadDobl_Laur_Poly_Convertors;         
        q : QuadDobl_Complex_Poly_Systems.Poly_Sys(qd_p'range)
          := Laurent_to_Polynomial_System(qd_p.all);
      begin
        tstart(timer);
       -- QuadDobl_Root_Refiner(q,qd_s);
        if nbtasks = 0 then
          Reporting_Root_Refiner
            (outfile,q,qd_s,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        else
          Multitasking_Refinement(q,qd_s,nbtasks,otp,epsxa,epsfa,maxit);
          new_line(outfile);
          put_line(outfile,"THE SOLUTIONS :");
          QuadDobl_Complex_Solutions_io.write(outfile,qd_s);
        end if;
        tstop(timer);
        QuadDobl_Complex_Poly_Systems.Clear(q);
      end;
    end if;
   -- QuadDobl_Complex_Solutions_io.write(outfile,qd_s);
    new_line(outfile);
    print_times(outfile,timer,"quad double Newton refinement");
    close(outfile);
  end QD_Root_Refinement;

  procedure DD_QD_Root_Refinement ( infilename,outfilename : in string ) is

    ans : character;

  begin
    new_line;
    put_line("MENU for Newton in double double or quad double arithmetic :");
    put_line("  1. use complex double double arithmetic;");
    put_line("  2. use complex quad double arithmetic.");
    put("Type 1 or 2 to make a choice : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then DD_Root_Refinement(infilename,outfilename);
     else QD_Root_Refinement(infilename,outfilename);
    end if;
  end DD_QD_Root_Refinement;

end Drivers_to_dd_qd_Root_Refiners;
