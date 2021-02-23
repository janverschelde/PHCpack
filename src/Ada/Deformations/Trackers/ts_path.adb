with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Double_Double_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Multprec_Floating_Vectors;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Multprec_Homotopy;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Multprec_System_and_Solutions_io;
with Continuation_Parameters;
with Standard_Continuation_Data;
with Standard_Path_Trackers;
with DoblDobl_Continuation_Data;
with DoblDobl_Path_Trackers;
with QuadDobl_Continuation_Data;
with QuadDobl_Path_Trackers;
with Multprec_Continuation_Data;
with Multprec_Path_Trackers;
with Main_Poly_Continuation;             use Main_Poly_Continuation;

procedure ts_path is

-- DESCRIPTION :
--   Test facility for the path trackers.

  tol : constant double_float := 10.0E-14;

  procedure Standard_Read_Homotopy
              ( lp,lq : in out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                qsols : in out Standard_Complex_Solutions.Solution_List;
                nbequ : out integer32 ) is

  -- DESCRIPTION :
  --   Reads the target system lp, start system lq and start solutions qsols,
  --   for standard complex numbers.
  --   The parameter nbequ equals lp'last if Gauss-Newton is needed.

    use Standard_Complex_Poly_Systems;

    targetfile,startfile : file_type;
    nbvar : integer32;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(targetfile);
    get(targetfile,lp);
    Close(targetfile);
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(startfile);
    Standard_System_and_Solutions_io.get(startfile,lq,qsols);
    nbvar := integer32(Number_of_Unknowns(lp(lp'first)));
    if lp'last > nbvar then
      nbequ := lp'last;
    else
      new_line;
      put("Apply Gauss-Newton correctors ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then nbequ := lp'last;
       else nbequ := 0;
      end if;
    end if;
    Close(startfile);
  end Standard_Read_Homotopy;

  procedure DoblDobl_Read_Homotopy
              ( lp,lq : in out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                qsols : in out DoblDobl_Complex_Solutions.Solution_List;
                nbequ : out integer32 ) is

  -- DESCRIPTION :
  --   Reads the target system lp, start system lq and start solutions qsols,
  --   for standard complex numbers.
  --   The parameter nbequ equals lp'last if Gauss-Newton is needed.

    use DoblDobl_Complex_Poly_Systems;

    targetfile,startfile : file_type;
    nbvar : integer32;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(targetfile);
    get(targetfile,lp);
    Close(targetfile);
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(startfile);
    DoblDobl_System_and_Solutions_io.get(startfile,lq,qsols);
    nbvar := integer32(Number_of_Unknowns(lp(lp'first)));
    if lp'last > nbvar then
      nbequ := lp'last;
    else
      new_line;
      put("Apply Gauss-Newton correctors ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then nbequ := lp'last;
       else nbequ := 0;
      end if;
    end if;
    Close(startfile);
  end DoblDobl_Read_Homotopy;

  procedure QuadDobl_Read_Homotopy
              ( lp,lq : in out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                qsols : in out QuadDobl_Complex_Solutions.Solution_List;
                nbequ : out integer32 ) is

  -- DESCRIPTION :
  --   Reads the target system lp, start system lq and start solutions qsols,
  --   for standard complex numbers.
  --   The parameter nbequ equals lp'last if Gauss-Newton is needed.

    use QuadDobl_Complex_Poly_Systems;

    targetfile,startfile : file_type;
    nbvar : integer32;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(targetfile);
    get(targetfile,lp);
    Close(targetfile);
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(startfile);
    QuadDobl_System_and_Solutions_io.get(startfile,lq,qsols);
    nbvar := integer32(Number_of_Unknowns(lp(lp'first)));
    if lp'last > nbvar then
      nbequ := lp'last;
    else
      new_line;
      put("Apply Gauss-Newton correctors ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then nbequ := lp'last;
       else nbequ := 0;
      end if;
    end if;
    Close(startfile);
  end QuadDobl_Read_Homotopy;

  procedure Multprec_Read_Homotopy
              ( lp,lq : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                qsols : in out Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Reads the target system lp, start system lq and start solutions qsols,
  --   for multi-precision numbers.

    targetfile,startfile : file_type;

  begin
    new_line;
    put_line("Reading the name of the file for the target system.");
    Read_Name_and_Open_File(targetfile);
    get(targetfile,lp);
    Close(targetfile);
    new_line;
    put_line("Reading the name of the file for the start system.");
    Read_Name_and_Open_File(startfile);
    Multprec_System_and_Solutions_io.get(startfile,lq,qsols);
    Close(startfile);
  end Multprec_Read_Homotopy;

  procedure Call_Standard_Path_Trackers
              ( file : in file_type;
                sols : in out Standard_Complex_Solutions.Solution_List;
                nbequ : in integer32 := 0 ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;
    use Standard_Continuation_Data;
    use Standard_Path_Trackers;

    tmp : Solution_List := sols;
    patpp : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_for_Path;
    endpp : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_End_Game;
    patcp : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_for_Path;
    endcp : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_End_Game;
    pathdir : Standard_Floating_Vectors.Link_to_Vector;
    wnd : integer32 := 1;
    errv : double_float := 0.0;
    order : constant integer32 
          := integer32(Continuation_Parameters.endext_order);
    timer : Timing_Widget;

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue
            (Max_Norm,Standard_Homotopy.Eval,
                      Standard_Homotopy.Diff,Standard_Homotopy.Diff);

    procedure Continue_in_End_Game is
      new Linear_Single_Conditioned_Reporting_Continue
            (Max_Norm,Standard_Homotopy.Eval,
                      Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Shallow_Create(ls);
      begin
        Continue_along_Path(file,s,Create(1.0),tol,false,patpp,patcp,nbequ);
        Continue_in_End_Game(file,s,Create(1.0),tol,false,
                             order,wnd,pathdir,errv,endpp,endcp,nbequ);
        ls.err := s.cora; ls.rco := s.rcond; ls.res := s.resa;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"path tracking with standard arithmetic");
  end Call_Standard_Path_Trackers;

  procedure Call_DoblDobl_Path_Trackers
              ( file : in file_type;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                nbequ : in integer32 := 0 ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Continuation_Data;
    use DoblDobl_Path_Trackers;

    one : constant double_double := create(1.0);
    target : constant Complex_Number := create(one);
    tmp : Solution_List := sols;
    patpp : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_for_Path;
    endpp : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_End_Game;
    patcp : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_for_Path;
    endcp : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_End_Game;
    w : integer32 := 1;
    pathdir : Double_Double_Vectors.Link_to_Vector;
    errv : double_double := create(0.0);
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);
    timer : Timing_Widget;

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue
            (Max_Norm,DoblDobl_Homotopy.Eval,
                      DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

    procedure Continue_in_End_Game is
      new Linear_Single_Conditioned_Reporting_Continue
            (Max_Norm,DoblDobl_Homotopy.Eval,
                      DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Shallow_Create(ls);
      begin
        Continue_along_Path(file,s,target,tol,false,patpp,patcp,nbequ);
        Continue_in_End_Game(file,s,target,tol,false,
                             order,w,pathdir,errv,endpp,endcp,nbequ);
        ls := Shallow_Create(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"path tracking with double double arithmetic");
  end Call_DoblDobl_Path_Trackers;

  procedure Call_QuadDobl_Path_Trackers
              ( file : in file_type;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                nbequ : in integer32 := 0 ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Continuation_Data;
    use QuadDobl_Path_Trackers;

    one : constant quad_double := create(1.0);
    target : constant Complex_Number := create(one);
    tmp : Solution_List := sols;
    patpp : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_for_Path;
    endpp : constant Continuation_Parameters.Pred_Pars
          := Continuation_Parameters.Create_End_Game;
    patcp : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_for_Path;
    endcp : constant Continuation_Parameters.Corr_Pars
          := Continuation_Parameters.Create_End_Game;
    w : integer32 := 1;
    pathdir : Quad_Double_Vectors.Link_to_Vector;
    errv : quad_double := create(0.0);
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);
    timer : Timing_Widget;

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue
            (Max_Norm,QuadDobl_Homotopy.Eval,
                      QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

    procedure Continue_in_End_Game is
      new Linear_Single_Conditioned_Reporting_Continue
            (Max_Norm,QuadDobl_Homotopy.Eval,
                      QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Shallow_Create(ls);
      begin
        Continue_along_Path(file,s,target,tol,false,patpp,patcp,nbequ);
        Continue_in_End_Game(file,s,target,tol,false,
                             order,w,pathdir,errv,endpp,endcp,nbequ);
        ls := Shallow_Create(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"path tracking with quad double arithmetic");
  end Call_QuadDobl_Path_Trackers;

  procedure Call_Multprec_Path_Trackers
              ( file : in file_type;
                sols : in out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;
    use Multprec_Continuation_Data;
    use Multprec_Path_Trackers;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    patpp : constant Pred_Pars
          := Convert(Continuation_Parameters.Create_for_Path);
    endpp : constant Pred_Pars
          := Convert(Continuation_Parameters.Create_End_Game);
    patcp : constant Corr_Pars
          := Convert(Continuation_Parameters.Create_for_Path);
    endcp : constant Corr_Pars
          := Convert(Continuation_Parameters.Create_End_Game);
    pathdir : Multprec_Floating_Vectors.Link_to_Vector;
    errv : Floating_Number;
    order : constant integer32
          := integer32(Continuation_Parameters.endext_order);
    timer : Timing_Widget;

    procedure Continue_along_Path is
      new Linear_Single_Normal_Reporting_Continue
            (Max_Norm,Multprec_Homotopy.Eval,
                      Multprec_Homotopy.Diff,Multprec_Homotopy.Diff);

    procedure Continue_in_End_Game is
      new Linear_Single_Conditioned_Reporting_Continue
            (Max_Norm,Multprec_Homotopy.Eval,
                      Multprec_Homotopy.Diff,Multprec_Homotopy.Diff);

  begin
    tstart(timer);
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
        s : Solu_Info := Deep_Create(ls.all);
      begin
        Continue_along_Path(file,s,Create(integer(1)),Create(tol),
                            false,patpp,patcp);
        Continue_in_End_Game(file,s,Create(integer(1)),Create(tol),false,
                             order,pathdir,errv,endpp,endcp);
        Append(res,res_last,Deep_Create(s));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,sols); Clear(res);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"path tracking in multiprecision arithmetic");
  end Call_Multprec_Path_Trackers;

  procedure Test_Standard_Path_Trackers is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    file : file_type;
    oc : natural32 := 0;
    lp,lq : Link_to_Poly_Sys;
    qsols : Solution_List;
    nbequ : integer32;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Standard_Read_Homotopy(lp,lq,qsols,nbequ);
    put_line(file,"TARGET SYSTEM : ");
    put(file,natural32(lp'last),lp.all);
    new_line(file);
    put_line(file,"START SYSTEM : ");
    put(file,natural32(lq'last),lq.all);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    if not Is_Null(qsols) then
      Check_Continuation_Parameter(qsols);
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      Standard_Homotopy.Create(lp.all,lq.all,1,Random1);
      new_line;
      Driver_for_Continuation_Parameters(file);
      new_line;
      Driver_for_Process_io(file,oc);
      new_line;
      put_line("See the output file for results ...");
      new_line;
      Call_Standard_Path_Trackers(file,qsols,nbequ);
      new_line(file);
      put_line(file,"THE COMPUTED SOLUTIONS : ");
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      Standard_Homotopy.Clear;
    end if;
  end Test_Standard_Path_Trackers;

  procedure Test_DoblDobl_Path_Trackers is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    file : file_type;
    oc : natural32 := 0;
    lp,lq : Link_to_Poly_Sys;
    qsols : Solution_List;
    ran : Complex_Number;
    nbequ : integer32;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    DoblDobl_Read_Homotopy(lp,lq,qsols,nbequ);
    put_line(file,"TARGET SYSTEM : ");
    put(file,lp.all);
    new_line(file);
    put_line(file,"START SYSTEM : ");
    put(file,lq.all);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    if not Is_Null(qsols) then
      Check_Continuation_Parameter(qsols);
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      ran := DoblDobl_Random_Numbers.Random1;
      DoblDobl_Homotopy.Create(lp.all,lq.all,1,ran);
      new_line;
      Continuation_Parameters.Tune(0,32);
      Driver_for_Continuation_Parameters(file);
      new_line;
      Driver_for_Process_io(file,oc);
      new_line;
      put_line("See the output file for results ...");
      new_line;
      Call_DoblDobl_Path_Trackers(file,qsols,nbequ);
      new_line(file);
      put_line(file,"THE COMPUTED SOLUTIONS : ");
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      DoblDobl_Homotopy.Clear;
    end if;
  end Test_DoblDobl_Path_Trackers;

  procedure Test_QuadDobl_Path_Trackers is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    file : file_type;
    oc : natural32 := 0;
    lp,lq : Link_to_Poly_Sys;
    qsols : Solution_List;
    one : constant quad_double := create(1.0);
    ran : Complex_Number := create(one);
    nbequ : integer32;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    QuadDobl_Read_Homotopy(lp,lq,qsols,nbequ);
    put_line(file,"TARGET SYSTEM : ");
    put(file,lp.all);
    new_line(file);
    put_line(file,"START SYSTEM : ");
    put(file,lq.all);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    if not Is_Null(qsols) then
      Check_Continuation_Parameter(qsols);
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      new_line(file);
      ran := QuadDobl_Random_Numbers.Random1;
      QuadDobl_Homotopy.Create(lp.all,lq.all,1,ran);
      new_line;
      Continuation_Parameters.Tune(0,64);
      Driver_for_Continuation_Parameters(file);
      new_line;
      Driver_for_Process_io(file,oc);
      new_line;
      put_line("See the output file for results ...");
      new_line;
      Call_QuadDobl_Path_Trackers(file,qsols,nbequ);
      new_line(file);
      put_line(file,"THE COMPUTED SOLUTIONS : ");
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      QuadDobl_Homotopy.Clear;
    end if;
  end Test_QuadDobl_Path_Trackers;

  procedure Test_Multprec_Path_Trackers is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Solutions;

    file : file_type;
    oc,deci,size : natural32 := 0;
    lp,lq : Link_to_Poly_Sys;
    qsols : Solution_List;
    ran : Complex_Number;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Multprec_Read_Homotopy(lp,lq,qsols);
    put_line(file,"TARGET SYSTEM : ");
    put(file,lp.all);
    new_line(file);
    put_line(file,"START SYSTEM : ");
    put(file,lq.all);
    new_line(file);
    put_line(file,"START SOLUTIONS : ");
    new_line;
    put("Give the number of decimal places : "); get(deci); skip_line;
    size := Decimal_to_Size(deci); 
    if not Is_Null(qsols) then
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      Set_Size(lp.all,size);
      Set_Size(lq.all,size);
      Set_Size(qsols,size);
      new_line(file);
      ran := Random(size);
      Multprec_Homotopy.Create(lp.all,lq.all,1,ran);
      new_line;
      Continuation_Parameters.Tune(0,deci);
      Driver_for_Continuation_Parameters(file);
      new_line;
      Driver_for_Process_io(file,oc);
      Call_Multprec_Path_Trackers(file,qsols);
      new_line(file);
      put_line(file,"THE COMPUTED SOLUTIONS : ");
      put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
      Multprec_Homotopy.Clear;
    end if;
  end Test_Multprec_Path_Trackers;

  procedure Main is
  
    ans : character;

  begin
    new_line;
    put_line("Testing the path trackers");
   -- loop
      new_line;
      put_line("Choose one of the following : ");
     -- put_line("  0. Exit this program.");
      put_line("  1. Test standard predictor-corrector method;");
      put_line("  2. run predictor-corrector with double double numbers;");
      put_line("  3. run predictor-corrector with quad double numbers;");
      put_line("  4. Test multi-precision predictor-corrector method.");
     -- put("Type 0, 1, 2, 3, or 4 to select : ");
      put("Type 1, 2, 3, or 4 to select : ");
      Ask_Alternative(ans,"01234");
     -- exit when ans = '0';
      case ans is
        when '1' => Test_Standard_Path_Trackers;
        when '2' => Test_DoblDobl_Path_Trackers;
        when '3' => Test_QuadDobl_Path_Trackers;
        when '4' => Test_Multprec_Path_Trackers;
        when others => null;
      end case;
   -- end loop;
  end Main;

begin
  Main;
end ts_path;
