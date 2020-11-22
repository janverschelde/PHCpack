with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Laurentials;
with DoblDobl_Complex_Laurentials;
with QuadDobl_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Laur_Poly_Convertors;
with Prompt_for_Systems;
with Prompt_for_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;             use QuadDobl_Root_Refiners;
with Root_Refining_Parameters;           use Root_Refining_Parameters;

package body Black_Box_Root_Refiners is

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;
    dim : constant integer32 := Head_Of(sols).n;

  begin
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    if p'last = dim then
      tstart(timer);
      Reporting_Root_Refiner
        (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
      tstop(timer);
    else
      tstart(timer);
      Reporting_Root_Sharpener
        (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
      tstop(timer);
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- NOTE : the deflate parameter is not used as deflation is not available
  --   for Laurent systems.  No Sharpener for overdetermined problems?

    use Standard_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;
    dim : constant integer32 := Head_Of(sols).n;

  begin
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false; -- not available
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    if p'last = dim then
      tstart(timer);
      Reporting_Root_Refiner
        (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,false);
      tstop(timer);
    else
      tstart(timer);
      Reporting_Root_Sharpener
        (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,false);
      tstop(timer);
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Refine_Roots
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;

  begin
    DoblDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    tstart(timer);
    Reporting_Root_Refiner
      (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Refine_Roots
               ( file : in file_type;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- NOTE : the deflate parameter is not used as deflation is not available
  --   for Laurent systems.  No Sharpener for overdetermined problems?

    use DoblDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;

  begin
    DoblDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false; -- not available
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    tstart(timer);
    Reporting_Root_Refiner
      (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,false);
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Refine_Roots
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;

  begin
    QuadDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    tstart(timer);
    Reporting_Root_Refiner
      (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Refine_Roots
               ( file : in file_type;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- NOTE : the deflate parameter is not used as deflation is not available
  --   for Laurent systems.  No Sharpener for overdetermined problems?

    use QuadDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;

  begin
    QuadDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false; -- not available
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    tstart(timer);
    Reporting_Root_Refiner
      (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,false);
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Standard_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    procedure Refine ( file : in out file_type; lp : in Link_to_Laur_Sys;
                       sysonfile : in boolean ) is

    -- DESCRIPTION :
    --   Calls the root refiner on the system lp,
    --   reading the solutions from file if sysonfile.

      outfile : file_type;
      sols : Solution_List;
      nbvar : constant natural32
            := Standard_Complex_Laurentials.Number_of_Unknowns(lp(lp'first));

    begin
      Create_Output_File(outfile,outfilename);
      if lp'last = integer32(nbvar)
       then put(outfile,natural32(lp'last),lp.all);
       else put(outfile,natural32(lp'last),nbvar,lp.all);
      end if;
      Prompt_for_Solutions.Read_Solutions(file,sysonfile,sols);
      if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
        Black_Box_Root_Refiners.Refine_Roots(outfile,lp.all,sols);
      else
        declare
          use Standard_Laur_Poly_Convertors;
          p : constant Poly_Sys(lp'range)
            := Positive_Laurent_Polynomial_System(lp.all);
        begin
          Black_Box_Root_Refiners.Refine_Roots(outfile,p,sols);
        end;
      end if;
    end Refine;

    procedure Main is

    -- DESCRIPTION :
    --   Reads the system, the solutions,
    --   and then calls the black box root refiner.

      infile : file_type;
      sysonfile : boolean;
      lp : Link_to_Laur_Sys := null;

    begin
      if verbose > 0 then
        put("At verbose level "); put(verbose,1);
        put_line(", in black_box_root_refiners.Standard_Main ...");
      end if;
      Prompt_for_Systems.Read_System(infile,infilename,lp,sysonfile);
      if lp /= null
       then Refine(infile,lp,sysonfile);
      end if;
    end Main;

  begin
    Main;
  end Standard_Main;

  procedure DoblDobl_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    procedure Refine ( file : in out file_type; lp : in Link_to_Laur_Sys;
                       sysonfile : in boolean ) is

    -- DESCRIPTION :
    --   Calls the root refiner on the system lp,
    --   reading the solutions from file if sysonfile.

      outfile : file_type;
      sols : Solution_List;
      nbvar : constant natural32
            := DoblDobl_Complex_Laurentials.Number_of_Unknowns(lp(lp'first));

    begin
      Create_Output_File(outfile,outfilename);
      if lp'last = integer32(nbvar)
       then put(outfile,natural32(lp'last),lp.all);
       else put(outfile,natural32(lp'last),nbvar,lp.all);
      end if;
      Prompt_for_Solutions.Read_Solutions(file,sysonfile,sols);
      if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
        Black_Box_Root_Refiners.Refine_Roots(outfile,lp.all,sols);
      else
        declare
          use DoblDobl_Laur_Poly_Convertors;
          p : constant Poly_Sys(lp'range)
            := Positive_Laurent_Polynomial_System(lp.all);
        begin
          Black_Box_Root_Refiners.Refine_Roots(outfile,p,sols);
        end;
      end if;
    end Refine;

    procedure Main is

    -- DESCRIPTION :
    --   Reads the system, the solutions,
    --   and then calls the black box root refiner.

      infile : file_type;
      lp : Link_to_Laur_Sys;
      sysonfile : boolean;

    begin
      if verbose > 0 then
        put("At verbose level "); put(verbose,1);
        put_line(", in black_box_root_refiners.DoblDobl_Main ...");
      end if;
      Prompt_for_Systems.Read_System(infile,infilename,lp,sysonfile);
      if lp /= null
       then Refine(infile,lp,sysonfile);
      end if;
    end Main;

  begin
    Main;
  end DoblDobl_Main;

  procedure QuadDobl_Main ( infilename,outfilename : in string;
                            verbose : in integer32 := 0 ) is


    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    procedure Refine ( file : in out file_type; lp : in Link_to_Laur_Sys;
                       sysonfile : in boolean ) is

    -- DESCRIPTION :
    --   Calls the root refiner on the system lp,
    --   reading the solutions from file if sysonfile.

      outfile : file_type;
      sols : Solution_List;
      nbvar : constant natural32
            := QuadDobl_Complex_Laurentials.Number_of_Unknowns(lp(lp'first));

    begin
      Create_Output_File(outfile,outfilename);
      if lp'last = integer32(nbvar)
       then put(outfile,natural32(lp'last),lp.all);
       else put(outfile,natural32(lp'last),nbvar,lp.all);
      end if;
      Prompt_for_Solutions.Read_Solutions(file,sysonfile,sols);
      if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
        Black_Box_Root_Refiners.Refine_Roots(outfile,lp.all,sols);
      else
        declare
          use QuadDobl_Laur_Poly_Convertors;
          p : constant Poly_Sys(lp'range)
            := Positive_Laurent_Polynomial_System(lp.all);
        begin
          Black_Box_Root_Refiners.Refine_Roots(outfile,p,sols);
        end;
      end if;
    end Refine;

    procedure Main is

    -- DESCRIPTION :
    --   Reads the system, the solutions,
    --   and then calls the black box root refiner.

      infile : file_type;
      lp : Link_to_Laur_Sys;
      sysonfile : boolean;

    begin
      if verbose > 0 then
        put("At verbose level "); put(verbose,1);
        put_line(", in black_box_root_refiners.QuadDobl_Main ...");
      end if;
      Prompt_for_Systems.Read_System(infile,infilename,lp,sysonfile);
      if lp /= null
       then Refine(infile,lp,sysonfile);
      end if;
    end Main;

  begin
    Main;
  end QuadDobl_Main;

end Black_Box_Root_Refiners;
