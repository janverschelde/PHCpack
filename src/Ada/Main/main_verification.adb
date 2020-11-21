with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;      use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;       use Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Poly_Strings;      use Multprec_Complex_Poly_Strings;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;   use Multprec_System_and_Solutions_io;
with Drivers_for_Condition_Tables;       use Drivers_for_Condition_Tables;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Multprec_Root_Refiners;             use Multprec_Root_Refiners;
with Multprec_Residual_Evaluations;      use Multprec_Residual_Evaluations;
with Symbolic_Symmetry_Group_io;
with Drivers_for_Symmetry_Group_io;      use Drivers_for_Symmetry_Group_io;
with Drivers_for_Orbits_of_Solutions;    use Drivers_for_Orbits_of_Solutions;
with Driver_for_Winding_Numbers;
with Drivers_to_Deflate_Singularities;   use Drivers_to_Deflate_Singularities;
with Standard_Multiplicity_Structure;
with DoblDobl_Multiplicity_Structure;
with QuadDobl_Multiplicity_Structure;
with Drivers_to_DD_QD_Root_Refiners;     use Drivers_to_DD_QD_Root_Refiners;
with Root_Refining_Parameters;           use Root_Refining_Parameters;
with Standard_Refiner_Circuits;
with valipoco;
with Bye_Bye_Message;
with Verification_of_Solutions;          use Verification_of_Solutions;
with Prompt_for_Systems;                 use Prompt_for_Systems;
with Prompt_for_Solutions;               use Prompt_for_Solutions;

package body Main_Verification is

  procedure Display_Verification_Info is

    i : array(1..13) of string(1..65);

  begin
    i(1):="A condition table for a list of solutions places the solutions in";
    i(2):="the list in various frequency tables, determining their  position";
    i(3):="based on logarithms of the corrector norms, condition number, and";
    i(4):="residual value.                                                  ";
    i(5):="Basic verification consists in the application of Newton's method";
    i(6):="on  the  list  of solutions.  There are facilities to extract the";
    i(7):="generating solutions when the symmetry group is submitted.       ";
    i(8):="  Winding  numbers  can  be  computed  by  homotopy  continuation";
    i(9):="methods.   The user must provide a start system with solutions at";
   i(10):="t < 1.                                                           ";
   i(11):="Polyhedral verification is based on  the  output  file  of  poco,";
   i(12):="where the polyhedral  end  game was turned on.  This verification";
   i(13):="puts up a frequency table of computed path directions.           ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Display_Verification_Info;

  procedure Refine_Roots
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                solsfile,invar,allperms,signsym : in boolean;
                v : in List_of_Permutations;
                epsxa,epsfa,tolsing : in double_float;
                maxit : in natural32; deflate : in out boolean;
                wout : in boolean;
                sols : in out Standard_Complex_Solutions.Solution_List;
                refsols : in out Standard_Complex_Solutions.Solution_List ) is

    numit : natural32 := 0;

  begin
    if solsfile or invar then
      if p'last /= Standard_Complex_Solutions.Head_Of(sols).n then
        Reporting_Root_Sharpener
          (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
      else 
        Reporting_Root_Refiner
          (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
        if invar
         then Driver_for_Orbits_of_Solutions
                (file,refsols,v,allperms,signsym,epsxa);
        end if;
      end if;
    else
      if p'last /= Standard_Complex_Solutions.Head_Of(sols).n then
        Reporting_Root_Sharpener
          (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
      else
        Reporting_Root_Refiner
          (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
      end if;
    end if;
  end Refine_Roots;

  procedure Refine_Roots
             ( file : in file_type;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               solsfile : in boolean;
               epsxa,epsfa,tolsing : in double_float; maxit : in natural32;
               deflate : in out boolean; wout : in boolean;
               sols : in out Standard_Complex_Solutions.Solution_List;
               refsols : in out Standard_Complex_Solutions.Solution_List ) is

    numit : natural32 := 0;

  begin
    if solsfile then
      if p'last /= Standard_Complex_Solutions.Head_Of(sols).n then
        Reporting_Root_Sharpener
          (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
      else
        Reporting_Root_Refiner
          (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
      end if;
    else
      if p'last /= Standard_Complex_Solutions.Head_Of(sols).n then
        Reporting_Root_Sharpener
          (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
      else
        Reporting_Root_Refiner
          (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
      end if;
    end if;
  end Refine_Roots;

  procedure Refine_Roots
             ( file : in file_type;
               p : in Standard_Complex_Laur_Systems.Laur_Sys;
               solsfile : in boolean;
               epsxa,epsfa,tolsing : in double_float; maxit : in natural32;
               wout : in boolean;
               sols : in out Standard_Complex_Solutions.Solution_List;
               refsols : in out Standard_Complex_Solutions.Solution_List ) is

    numit : natural32 := 0;

  begin
    if solsfile then
      if p'last /= Standard_Complex_Solutions.Head_Of(sols).n then
        Reporting_Root_Sharpener
          (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,wout);
      else
        Reporting_Root_Refiner
          (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,wout);
      end if;
    else
      if p'last /= Standard_Complex_Solutions.Head_Of(sols).n then
        Reporting_Root_Sharpener
          (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,wout);
      else
        Reporting_Root_Refiner
          (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,wout);
      end if;
    end if;
  end Refine_Roots;

  procedure Winding_Verification
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    infile,solsft,outfile : file_type;
    ans : character;
    sysonfile,solsfile,deflate,wout : boolean;
    sols,refsols: Standard_Complex_Solutions.Solution_List;
    epsxa,epsfa,tolsing : double_float;
    maxit : natural32;

  begin
    if vrblvl > 0
     then put_line("-> in main_verification.Winding_Verification ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all);
    Read_Solutions(infile,sysonfile,sols);
    new_line;
    put("Do you want the refined solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      solsfile := true;
      put_line("Reading the name of the file to write the solutions on.");
      Read_Name_and_Create_File(solsft);
    else
      solsfile := false;
    end if;
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (outfile,epsxa,epsfa,tolsing,maxit,deflate,wout);
    Driver_for_Winding_Numbers(outfile,lp.all,sols);
    tstart(timer);
    Refine_Roots(outfile,lp.all,solsfile,
                 epsxa,epsfa,tolsing,maxit,deflate,wout,sols,refsols);
    tstop(timer);
    if solsfile then
      put(solsft,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
      Close(solsft);
    end if;
    new_line(outfile);
    print_times(outfile,timer,"Root Refinement");
    Close(outfile);
  end Winding_Verification;

  procedure Standard_Weeding_Verification
              ( infile : in out file_type; outfilename : in string;
                lp : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sysonfile : in boolean ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    solsft,outfile : file_type;
    n,maxit : natural32;
    ans : character;
    solsfile,deflate,wout : boolean := false;
    invar,allperms,signsym,allsigns : boolean := false;
    g,v : List_of_Permutations;
    sols,refsols: Standard_Complex_Solutions.Solution_List;
    epsxa,epsfa,tolsing : double_float;
    nbequ : constant natural32 := natural32(lp'last);
    nbvar : constant natural32 := Number_of_Unknowns(lp(lp'first));

  begin
    Create_Output_File(outfile,outfilename);
    if nbequ = nbvar
     then put(outfile,nbequ,lp.all);
     else put(outfile,nbequ,nbvar,lp.all);
    end if;
    Read_Solutions(infile,sysonfile,sols);
    new_line;
    put("Is the system invariant under group actions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      invar := true; n := lp'length;
      Read_Symmetry_Group(n,g,v,allperms,signsym,allsigns);
      new_line(outfile);
      put_line(outfile,"THE SYMMETRY GROUP : ");
      new_line(outfile);
      Symbolic_Symmetry_Group_io.put(outfile,v);
      new_line(outfile);
    else
      invar := false;
    end if;
    new_line;
    put("Do you want the refined solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      solsfile := true;
      put_line("Reading the name of the file to write the solutions on.");
      Read_Name_and_Create_File(solsft);
    else
      solsfile := false;
    end if;
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    Standard_Menu_Root_Refining_Parameters
      (outfile,epsxa,epsfa,tolsing,maxit,deflate,wout);
    Communications_with_User.End_of_Input_Message;
   -- put("starting root refinement with #variables : ");
   -- put(Head_Of(sols).n,1); new_line;
    tstart(timer);
    Refine_Roots(outfile,lp.all,solsfile,invar,allperms,signsym,v,
                 epsxa,epsfa,tolsing,maxit,deflate,wout,sols,refsols);
    tstop(timer);
    if solsfile then
      put(solsft,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
      Close(solsft);
    end if;
    new_line(outfile);
    print_times(outfile,timer,"Root Refinement");
    Close(outfile);
  end Standard_Weeding_Verification;

  procedure Standard_Weeding_Verification
              ( infile : in out file_type; outfilename : in string;
                lp : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                sysonfile : in boolean ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Solutions;

    outfile,solsft : file_type;
    nbequ : constant natural32 := natural32(lp'last);
    nbvar : constant natural32 := Number_of_Unknowns(lp(lp'first));
    ans : character;
    sols,refsols : Solution_List;
    solsfile,deflate,wout : boolean;
    maxit : natural32;
    epsxa,epsfa,tolsing : double_float;
    timer : Timing_Widget;

  begin
    Create_Output_File(outfile,outfilename);
    if nbequ = nbvar
     then put(outfile,nbequ,lp.all);
     else put(outfile,nbequ,nbvar,lp.all);
    end if;
    Read_Solutions(infile,sysonfile,sols);
    new_line;
    put("Do you want the refined solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      solsfile := true;
      put_line("Reading the name of the file to write the solutions on.");
      Read_Name_and_Create_File(solsft);
    else
      solsfile := false;
    end if;
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    deflate := false; -- no deflation for Laurent systems
    Standard_Menu_Root_Refining_Parameters
      (outfile,epsxa,epsfa,tolsing,maxit,deflate,wout);
    Communications_with_User.End_of_Input_Message;
    tstart(timer);
    Refine_Roots(outfile,lp.all,solsfile,
                 epsxa,epsfa,tolsing,maxit,wout,sols,refsols);
    tstop(timer);
    if solsfile then
      put(solsft,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
      Close(solsft);
    end if;
    new_line(outfile);
    print_times(outfile,timer,"Root Refinement");
    Close(outfile);
  end Standard_Weeding_Verification;

  procedure Standard_Weeding_Verification
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    infile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sysonfile : boolean;

    use Standard_Laur_Poly_Convertors;

  begin
    if vrblvl > 0 then
      put_line("-> in main_verification.Standard_Weeding_Verification ...");
    end if;
    Read_System(infile,infilename,lq,sysonfile);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Standard_Weeding_Verification(infile,outfilename,lq,sysonfile);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Standard_Weeding_Verification(infile,outfilename,lp,sysonfile);
    end if;
  end Standard_Weeding_Verification;

  procedure Multprec_Residual_Evaluator
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    timer : Timing_Widget;
    infile,outfile : file_type;
    sysonfile : boolean;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put_line("-> in main_verification.Multprec_Residual_Evaluator ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all);
    Read_Solutions(infile,sysonfile,sols);
    new_line(outfile);
    put_line(outfile,"the solutions in standard precision : ");
    put(outfile,Standard_Complex_Solutions.Length_Of(sols),
        natural32(lp'last),sols);
    declare
      mpsols : Multprec_Complex_Solutions.Solution_List
             := Multprec_Complex_Solutions.Create(sols);
      mp : constant Multprec_Complex_Poly_Systems.Poly_Sys(lp'range)
         := Convert(lp.all);
      mp_eval : constant Multprec_Complex_Poly_SysFun.Eval_Poly_Sys(mp'range)
              := Create(mp);
      deci,size : natural32 := 0;
    begin
      new_line;
      put("Give the number of decimal places : "); get(deci);
      size := Decimal_to_Size(deci);
      Multprec_Complex_Solutions.Set_Size(mpsols,size);
      new_line(outfile);
      put(outfile,"THE RESIDUALS with "); put(outfile,deci,1);
      put_line(outfile," decimal places :");
      tstart(timer);
      Residuals(outfile,mp_eval,mpsols);            
      tstop(timer);
    end;
    new_line(outfile);
    print_times(outfile,timer,"Multi-Precision Residual Evaluation");
    Close(outfile);
  end Multprec_Residual_Evaluator;

  procedure Call_Multprec_Root_Refiner
               ( file : in file_type; n,m : in natural32;
                 ls : in Link_to_Array_of_Strings;
                 sols : in out Multprec_Complex_Solutions.Solution_List ) is

    timer : Timing_Widget;
    epsxa,epsfa,tolsing : Floating_Number;
    maxit,numit,deci,size : natural32 := 0;
    deflate,wout : boolean;
    p : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n));

  begin
    new_line;
    Multprec_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deci,deflate,wout);
    Multprec_Menu_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deci,deflate,wout);
    size := Decimal_to_Size(deci);
   -- put("Give the size of the numbers : "); get(size);
   -- Symbol_Table.Init(m);
    if Symbol_Table.Number < m
     then Symbol_Table.Init(m);
    end if;
   -- put("m = "); put(m,1); new_line;
   -- for i in ls'range loop
   --   put("p["); put(i,1); put("] = "); put_line(ls(i).all);
   -- end loop;
   -- put_line("parsing the system from strings ...");
    p := Parse(m,size,ls.all);
   -- put_line("done parsing the system");
    Multprec_Complex_Solutions.Set_Size(sols,size);
    Communications_with_User.End_of_Input_Message;
    numit := 0;
    tstart(timer);
    if p'last /= Multprec_Complex_Solutions.Head_Of(sols).n then
      Reporting_Root_Sharpener
        (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
    else
      Reporting_Root_Refiner
        (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,deflate,wout);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Multi-Precision Root Refinement");
  end Call_Multprec_Root_Refiner;

  procedure Call_Varbprec_Root_Refiner
               ( file : in file_type;
                 ls : in Link_to_Array_of_Strings;
                 sols : in out Multprec_Complex_Solutions.Solution_List ) is

    timer : Timing_Widget;
    wanted,maxitr,maxprc : natural32;
    verbose : boolean;

    use Multprec_Complex_Solutions;

  begin
    Menu_to_Set_Parameters(wanted,maxitr,maxprc,verbose);
    new_line(file);
    Write_Parameters(file,wanted,maxitr,maxprc,verbose);
    tstart(timer);
    if verbose then
      Verify_Solutions_of_Laurent_Polynomials
        (file,ls.all,sols,wanted,maxitr,maxprc);
    else
      Verify_Solutions_of_Laurent_Polynomials
        (ls.all,sols,wanted,maxitr,maxprc);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    new_line(file);
    print_times(file,timer,"Variable precision Newton steps");
  end Call_Varbprec_Root_Refiner;

  procedure Multprec_Weeding_Verification
              ( infilename,outfilename : in string;
                varbprec : in boolean ) is

    n,m : natural32;
    ls : Link_to_Array_of_Strings;
    infile,outfile : file_type;
    sols : Multprec_Complex_Solutions.Solution_List;

  begin
    if infilename = "" then
      get(n,m,ls,sols);
    else
      Open_Input_File(infile,infilename);
      get(infile,n,m,ls,sols);
      Close(infile);
    end if;
   -- put_line("reading the strings and solutions from file ...");
    Create_Output_File(outfile,outfilename);
    put(outfile,n,1);
    if m /= n
     then put(outfile,"  "); put(outfile,m,1);
    end if;
    new_line(outfile);
    for i in ls'range loop
      put_line(outfile,ls(i).all);
    end loop;
    if Multprec_Complex_Solutions.Length_Of(sols) > 0 then
      new_line(outfile);
      put_line(outfile,"the given solutions : ");
      put(outfile,Multprec_Complex_Solutions.Length_Of(sols),
          natural32(Multprec_Complex_Solutions.Head_Of(sols).n),sols);
      if varbprec
       then Call_Varbprec_Root_Refiner(outfile,ls,sols);
       else Call_Multprec_Root_Refiner(outfile,n,m,ls,sols);
      end if;
    end if;
    Close(outfile);
  end Multprec_Weeding_Verification;

  procedure Polyhedral_End_Game_Verification is

    pocofile,resultfile : file_type;

  begin
    new_line;
    put_line("Reading name of the output file of poco.");
    Read_Name_and_Open_File(pocofile);
    new_line;
    put_line("Reading name of output file.");
    Read_Name_and_Create_File(resultfile);
    Communications_with_User.End_of_Input_Message;
    valipoco(pocofile,resultfile);
    Close(pocofile);
    new_line(resultfile);
    put(resultfile,Bye_Bye_Message);
    Close(resultfile);
  end Polyhedral_End_Game_Verification;

  procedure Standard_Newton_with_Deflation
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    sysonfile : boolean;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    nbequ,nbvar : natural32;

  begin
    if vrblvl > 0 then
      put("-> in main_verification.");
      put_line("Standard_Newton_with_Deflation ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    nbequ := natural32(lp'last);
    nbvar := Number_of_Unknowns(lp(lp'first));
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file...");
      declare
        new_outfilename : constant string := Read_String;
      begin
        Create_Output_File(outfile,new_outfilename);
        if nbequ = nbvar
         then put(outfile,nbequ,lp.all);
         else put(outfile,nbequ,nbvar,lp.all);
        end if;
        Read_Solutions(infile,sysonfile,sols);
        new_line(outfile);
        put_line(outfile,"the solutions on input :");
        put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        Deflate_Singularities(outfile,new_outfilename,lp.all,sols,vrblvl-1);
      end;
    else
      Create_Output_File(outfile,outfilename);
      if nbequ = nbvar
       then put(outfile,nbequ,lp.all);
       else put(outfile,nbequ,nbvar,lp.all);
      end if;
      Read_Solutions(infile,sysonfile,sols);
      new_line(outfile);
      put_line(outfile,"the solutions on input :");
      put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Deflate_Singularities(outfile,outfilename,lp.all,sols,vrblvl-1);
    end if;
  end Standard_Newton_with_Deflation;

  procedure DoblDobl_Newton_with_Deflation
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;

    infile,outfile : file_type;
    sysonfile : boolean;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    nbequ,nbvar : natural32;

  begin
    if vrblvl > 0 then
      put("-> in main_verification.");
      put_line("DoblDobl_Newton_with_Deflation ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    nbequ := natural32(lp'last);
    nbvar := Number_of_Unknowns(lp(lp'first));
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file...");
      declare
        new_outfilename : constant string := Read_String;
      begin
        Create_Output_File(outfile,new_outfilename);
        if nbequ = nbvar
         then put(outfile,nbequ,lp.all);
         else put(outfile,nbequ,nbvar,lp.all);
        end if;
        Read_Solutions(infile,sysonfile,sols);
        new_line(outfile);
        put_line(outfile,"the solutions on input :");
        put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        Deflate_Singularities(outfile,new_outfilename,lp.all,sols,vrblvl-1);
      end;
    else
      Create_Output_File(outfile,outfilename);
      if nbequ = nbvar
       then put(outfile,nbequ,lp.all);
       else put(outfile,nbequ,nbvar,lp.all);
      end if;
      Read_Solutions(infile,sysonfile,sols);
      new_line(outfile);
      put_line(outfile,"the solutions on input :");
      put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Deflate_Singularities(outfile,outfilename,lp.all,sols,vrblvl-1);
    end if;
  end DoblDobl_Newton_with_Deflation;

  procedure QuadDobl_Newton_with_Deflation
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;

    infile,outfile : file_type;
    sysonfile : boolean;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    nbequ,nbvar : natural32;

  begin
    if vrblvl > 0 then
      put("-> in main_verification.");
      put_line("QuadDobl_Newton_with_Deflation ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    nbequ := natural32(lp'last);
    nbvar := Number_of_Unknowns(lp(lp'first));
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file...");
      declare
        new_outfilename : constant string := Read_String;
      begin
        Create_Output_File(outfile,new_outfilename);
        if nbequ = nbvar
         then put(outfile,nbequ,lp.all);
         else put(outfile,nbequ,nbvar,lp.all);
        end if;
        Read_Solutions(infile,sysonfile,sols);
        new_line(outfile);
        put_line(outfile,"the solutions on input :");
        put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        Deflate_Singularities(outfile,new_outfilename,lp.all,sols,vrblvl-1);
      end;
    else
      Create_Output_File(outfile,outfilename);
      if nbequ = nbvar
       then put(outfile,nbequ,lp.all);
       else put(outfile,nbequ,nbvar,lp.all);
      end if;
      Read_Solutions(infile,sysonfile,sols);
      new_line(outfile);
      put_line(outfile,"the solutions on input :");
      put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Deflate_Singularities(outfile,outfilename,lp.all,sols,vrblvl-1);
    end if;
  end QuadDobl_Newton_with_Deflation;

  procedure Newton_with_Deflation
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    prc : character;

  begin
    if vrblvl > 0
     then put_line("-> in mainvali.Newton_with_Deflation ...");
    end if;
    prc := Communications_with_User.Prompt_for_Precision;
    case prc is
      when '0' =>
        Standard_Newton_with_Deflation(infilename,outfilename,vrblvl-1);
      when '1' =>
        DoblDobl_Newton_with_Deflation(infilename,outfilename,vrblvl-1);
      when '2' =>
        QuadDobl_Newton_with_Deflation(infilename,outfilename,vrblvl-1);
      when others => null;
    end case;
  end Newton_with_Deflation;

  procedure Standard_Multiplicity
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    use Standard_Multiplicity_Structure;

    infile,outfile : file_type;
    sysonfile : boolean;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in main_verification.Standard_Multiplicity ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    put(outfile,natural32(lp'last),lp.all);
    Read_Solutions(infile,sysonfile,sols);
    Driver_to_Multiplicity_Structure(outfile,lp.all,sols);
  end Standard_Multiplicity;

  procedure DoblDobl_Multiplicity
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Multiplicity_Structure;

    infile,outfile : file_type;
    sysonfile : boolean;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in main_verification.DoblDobl_Multiplicity ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    put(outfile,natural32(lp'last),lp.all);
    Read_Solutions(infile,sysonfile,sols);
    Driver_to_Multiplicity_Structure(outfile,lp.all,sols);
  end DoblDobl_Multiplicity;

  procedure QuadDobl_Multiplicity
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Multiplicity_Structure;

    infile,outfile : file_type;
    sysonfile : boolean;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0
     then put_line("-> in main_verification.QuadDobl_Multiplicity ...");
    end if;
    Read_System(infile,infilename,lp,sysonfile);
    if outfilename = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
    else
      Create_Output_File(outfile,outfilename);
    end if;
    put(outfile,natural32(lp'last),lp.all);
    Read_Solutions(infile,sysonfile,sols);
    Driver_to_Multiplicity_Structure(outfile,lp.all,sols);
  end QuadDobl_Multiplicity;

  procedure Multiplicity_Structure
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    prc : character;

  begin
    if vrblvl > 0
     then put_line("-> in main_verification.Multiplicity_Structure ...");
    end if;
    prc := Communications_with_User.Prompt_for_Precision;
    case prc is
      when '0' => Standard_Multiplicity(infilename,outfilename,vrblvl-1);
      when '1' => DoblDobl_Multiplicity(infilename,outfilename,vrblvl-1);
      when '2' => QuadDobl_Multiplicity(infilename,outfilename,vrblvl-1);
      when others => null;
    end case;
  end Multiplicity_Structure;

  procedure Solution_Scanner
              ( infilename,outfilename : in string;
                vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0
     then put_line("-> in main_verification.Solution_Scanner ...");
    end if;
    new_line;
    put("Run Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Refiner_Circuits.Main(infilename,outfilename,vrblvl-1);
     else Main_Driver_to_Scan_Solution_Lists(infilename,outfilename,vrblvl-1);
    end if;
  end Solution_Scanner;

  procedure Main ( infilename,outfilename : in string;
                   verbose : in integer32 := 0 ) is

    ans : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in main_verification.Main ...");
    end if;
    loop
      new_line;
      put_line("MENU with Verification Methods : ");
      put_line
         ("  0. Scanning (huge) solution files and creating condition tables;");
      put_line
       ("  1. Basic Verification : refining and weeding out the solution set;");
      put_line
         ("  2. Evaluation of the residuals using multi-precision arithmetic;");
      put_line
         ("  3. Newton's method using multi-precision arithmetic;");
      put_line
         ("  4. Winding-Number Computation by homotopy continuation;");
      put_line
         ("  5. Polyhedral Verification : frequency table of path directions;");
      put_line
         ("  6. Newton's method with deflation for isolated singularities;");
      put_line
         ("  7. Multiplicity structure of isolated singular solutions;");
      put_line
         ("  8. Newton's method in double double or quad double arithmetic;");
      put_line
         ("  9. variable precision Newton steps to the desired accuracy.");
      put("Type 0, 1, 2, 3, 4, 5, 6, 7, 8, or 9 to select, or i for info : ");
      Ask_Alternative(ans,"0123456789i");
      case ans is
        when 'i' => new_line;
                    Display_Verification_Info;
        when '0' => Solution_Scanner(infilename,outfilename,verbose-1);
        when '1' =>
          Standard_Weeding_Verification(infilename,outfilename,verbose-1);
        when '2' =>
          Multprec_Residual_Evaluator(infilename,outfilename,verbose-1);
        when '3' =>
          Multprec_Weeding_Verification(infilename,outfilename,false);
        when '4' => Winding_Verification(infilename,outfilename,verbose-1);
        when '5' => Polyhedral_End_Game_Verification;
        when '6' => Newton_with_Deflation(infilename,outfilename,verbose-1);
        when '7' => Multiplicity_Structure(infilename,outfilename,verbose-1);
        when '8' => DD_QD_Root_Refinement(infilename,outfilename);
        when '9' =>
          Multprec_Weeding_Verification(infilename,outfilename,true);
        when others => null;
      end case;
      exit when ans /= 'i';
    end loop;
  end Main;

end Main_Verification;
