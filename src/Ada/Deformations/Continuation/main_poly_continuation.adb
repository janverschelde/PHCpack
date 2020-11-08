with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;        use QuadDobl_Complex_Numbers_cv;
with Multprec_Floating_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Numbers_io;                         use Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with Quad_Double_Vectors;
with Quad_Double_VecVecs;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions_io;      use Multprec_Complex_Solutions_io;
with Standard_Homotopy;
with Standard_Laurent_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Laurent_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Laurent_Homotopy;
with Multprec_Homotopy;
with Drivers_for_Homotopy_Creation;      use Drivers_for_Homotopy_Creation;
with Continuation_Parameters;
with Continuation_Parameters_io;
with Standard_IncFix_Continuation;
with DoblDobl_IncFix_Continuation;
with QuadDobl_IncFix_Continuation;
with Multprec_IncFix_Continuation;
with Process_io;                         use Process_io;
with Drivers_for_Path_Directions;        use Drivers_for_Path_Directions;

package body Main_Poly_Continuation is

-- AUXILIARIES :

  procedure LaurCont
              ( file : in file_type;
                sols : in out Standard_Complex_Solutions.Solution_List;
                proj,report : in boolean; nbq : in integer32 := 0;
                target : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path trackers for standard complex solutions,
  --   for Laurent homotopies in standard double precision.

    use Standard_IncFix_Continuation;
   
    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Standard_Laurent_Homotopy.Eval,
                          Standard_Laurent_Homotopy.Diff,
                          Standard_Laurent_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Standard_Laurent_Homotopy.Eval,
                             Standard_Laurent_Homotopy.Diff,
                             Standard_Laurent_Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,nbq,target=>target);
     else Sil_Cont(sols,proj,nbq,target=>target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end LaurCont;

  procedure LaurCont
              ( file : in file_type;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean; nbq : in integer32 := 0;
                target : in DoblDobl_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path trackers for double double complex solutions,
  --   for Laurent homotopies in double double precision.

    use DoblDobl_IncFix_Continuation;
   
    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,DoblDobl_Laurent_Homotopy.Eval,
                          DoblDobl_Laurent_Homotopy.Diff,
                          DoblDobl_Laurent_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,DoblDobl_Laurent_Homotopy.Eval,
                             DoblDobl_Laurent_Homotopy.Diff,
                             DoblDobl_Laurent_Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,nbq,target=>target);
     else Sil_Cont(sols,nbq,target=>target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end LaurCont;

  procedure LaurCont
              ( file : in file_type;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                report : in boolean; nbq : in integer32 := 0;
                target : in QuadDobl_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path trackers for quad double complex solutions,
  --   for Laurent homotopies in quad double precision.

    use QuadDobl_IncFix_Continuation;
   
    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,QuadDobl_Laurent_Homotopy.Eval,
                          QuadDobl_Laurent_Homotopy.Diff,
                          QuadDobl_Laurent_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,QuadDobl_Laurent_Homotopy.Eval,
                             QuadDobl_Laurent_Homotopy.Diff,
                             QuadDobl_Laurent_Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,nbq,target=>target);
     else Sil_Cont(sols,nbq,target=>target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end LaurCont;

  procedure Continue
              ( file : in file_type;
                sols : in out Standard_Complex_Solutions.Solution_List;
                proj,report : in boolean; nbq : in integer32 := 0;
                target : in Standard_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path trackers for standard complex solutions,
  --   for polynomial homotopies in standard double precision.

    use Standard_IncFix_Continuation;
   
    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);
     -- new Silent_Continue(Max_Norm,Standard_Coefficient_Homotopy.Eval,
     --                     Standard_Homotopy.Diff,
     --                     Standard_Coefficient_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);
     -- new Reporting_Continue(Max_Norm,Standard_Coefficient_Homotopy.Eval,
     --                        Standard_Homotopy.Diff,
     --                        Standard_Coefficient_Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,nbq,target=>target);
     else Sil_Cont(sols,proj,nbq,target=>target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end Continue;

  procedure Continue
              ( file : in file_type;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                report : in boolean; nbq : in integer32 := 0;
                target : in DoblDobl_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path trackers for double double complex solutions,
  --   for polynomial homotopies in double double precision.

    use DoblDobl_IncFix_Continuation;
   
    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                          DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,DoblDobl_Homotopy.Eval,
                             DoblDobl_Homotopy.Diff,DoblDobl_Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,nbq,target=>target);
     else Sil_Cont(sols,nbq,target=>target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end Continue;

  procedure Continue
              ( file : in file_type;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                report : in boolean; nbq : in integer32 := 0;
                target : in QuadDobl_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path trackers for quad double complex solutions,
  --   for polynomial homotopies in quad double precision.

    use QuadDobl_IncFix_Continuation;
   
    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                          QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,QuadDobl_Homotopy.Eval,
                             QuadDobl_Homotopy.Diff,QuadDobl_Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,nbq,target=>target);
     else Sil_Cont(sols,nbq,target=>target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end Continue;

  procedure Continue
              ( file : in file_type;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                proj,report : in boolean;
                target : in Multprec_Complex_Numbers.Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path trackers for multiprecision complex solutions,
  --   for polynomial homotopies in arbitrary multiprecision.

    use Multprec_IncFix_Continuation;

    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Multprec_Homotopy.Eval,
                          Multprec_Homotopy.Diff,Multprec_Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Multprec_Homotopy.Eval,
                             Multprec_Homotopy.Diff,Multprec_Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,target);
     else Sil_Cont(sols,proj,target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end Continue;

  procedure Ask_Symbol is

    sb : Symbol;

  begin
    put("Give symbol to display additional unknown : ");
    sb := (sb'range => ' ');
    Symbol_Table.Enlarge(1);
    Symbol_Table_io.Get(sb);
    Symbol_Table.Add(sb);
  end Ask_Symbol;

-- TARGET ROUTINES :

  procedure Driver_for_Process_io ( file : in file_type ) is

    ans : character;
    m : array(0..8) of string(1..65);

  begin
    put_line("MENU for Output Information during Continuation : ");
    m(0):="  0 : no intermediate output information during continuation     ";
    m(1):="  1 : only the final solutions at the end of the paths           ";
    m(2):="  2 : intermediate solutions at each step along the paths        ";
    m(3):="  3 : information of the predictor: t and step length            ";
    m(4):="  4 : information of the corrector: corrections and residuals    ";
    m(5):="  5 : intermediate solutions and information of the predictor    ";
    m(6):="  6 : intermediate solutions and information of the corrector    ";
    m(7):="  7 : information of predictor and corrector                     ";
    m(8):="  8 : intermediate solutions, info of predictor and corrector    ";
    for i in m'range loop
      put_line(m(i));
    end loop;
    put("Type a number between 0 and 8 to select output information : ");
    Ask_Alternative(ans,"012345678");
    new_line(file);
    put_line(file,"OUTPUT INFORMATION DURING CONTINUATION :");
    case ans is
      when '0' => Set_output_code(nil); put_line(file,m(0));
      when '1' => Set_output_code(nil); put_line(file,m(1));
      when '2' => Set_output_code(s);   put_line(file,m(2));
      when '3' => Set_output_code(p);   put_line(file,m(3));
      when '4' => Set_output_code(c);   put_line(file,m(4));
      when '5' => Set_output_code(sp);  put_line(file,m(5));
      when '6' => Set_output_code(sc);  put_line(file,m(6));
      when '7' => Set_output_code(pc);  put_line(file,m(7));
      when '8' => Set_output_code(spc); put_line(file,m(8));
      when others => null;
    end case;
  end Driver_for_Process_io;

  procedure Driver_for_Process_io ( file : in file_type; oc : out natural32 ) is

    ans : character;
    m : array(0..8) of string(1..65);

  begin
    put_line("MENU for Output Information during Continuation : ");
    m(0):="  0 : no intermediate output information during continuation     ";
    m(1):="  1 : only the final solutions at the end of the paths           ";
    m(2):="  2 : intermediate solutions at each step along the paths        ";
    m(3):="  3 : information of the predictor: t and step length            ";
    m(4):="  4 : information of the corrector: corrections and residuals    ";
    m(5):="  5 : intermediate solutions and information of the predictor    ";
    m(6):="  6 : intermediate solutions and information of the corrector    ";
    m(7):="  7 : information of predictor and corrector                     ";
    m(8):="  8 : intermediate solutions, info of predictor and corrector    ";
    for i in m'range loop
      put_line(m(i));
    end loop;
    put("Type a number between 0 and 8 to select output information : ");
    Ask_Alternative(ans,"012345678");
    new_line(file);
    put_line(file,"OUTPUT INFORMATION DURING CONTINUATION :");
    case ans is
      when '0' => Set_output_code(nil); oc := 0; put_line(file,m(0));
      when '1' => Set_output_code(nil); oc := 1; put_line(file,m(1));
      when '2' => Set_output_code(s);   oc := 2; put_line(file,m(2));
      when '3' => Set_output_code(p);   oc := 3; put_line(file,m(3));
      when '4' => Set_output_code(c);   oc := 4; put_line(file,m(4));
      when '5' => Set_output_code(sp);  oc := 5; put_line(file,m(5));
      when '6' => Set_output_code(sc);  oc := 6; put_line(file,m(6));
      when '7' => Set_output_code(pc);  oc := 7; put_line(file,m(7));
      when '8' => Set_output_code(spc); oc := 8; put_line(file,m(8));
      when others => null;
    end case;
  end Driver_for_Process_io;

  procedure Begin_Banner ( file : in file_type ) is
  begin
    put_line(file,"****************** CURRENT CONTINUATION PARAMETERS "
                   & "*****************");
  end Begin_Banner;

  procedure End_Banner ( file : in file_type ) is
  begin
    put_line(file,"***************************************************"
                   & "*****************");
  end End_Banner;

  procedure Driver_for_Continuation_Parameters is

    nb : natural32 := 0;

  begin
    loop
      Begin_Banner(Standard_Output);
      Continuation_Parameters_io.put;
      End_Banner(Standard_Output);
      Continuation_Parameters_io.get(nb);
      exit when (nb = 0);
    end loop;
  end Driver_for_Continuation_Parameters;

  procedure Driver_for_Continuation_Parameters ( file : in file_type ) is
  begin
    Driver_for_Continuation_Parameters;
    new_line(file);
    Begin_Banner(file);
    Continuation_Parameters_io.put(file);
    End_Banner(file);
  end Driver_for_Continuation_Parameters;

  procedure Driver_for_Continuation_Parameters
              ( precision : in natural32 ) is

    nb : natural32 := 0;

  begin
    Continuation_Parameters.Tune(0); -- ,precision); -- leave default
    loop
      Begin_Banner(Standard_Output);
      Continuation_Parameters_io.put;
      End_Banner(Standard_Output);
      Continuation_Parameters_io.get(nb);
      exit when (nb = 0);
    end loop;
  end Driver_for_Continuation_Parameters;

  procedure Driver_for_Continuation_Parameters
              ( file : in file_type; precision : in natural32 ) is
  begin
    Driver_for_Continuation_Parameters(precision);
    new_line(file);
    Begin_Banner(file);
    Continuation_Parameters_io.put(file);
    End_Banner(file);
  end Driver_for_Continuation_Parameters;

  procedure Check_Continuation_Parameter
                ( sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    ans : character;
    tre,tim : double_float;

  begin
    if not Is_Null(sols) then
      if Head_Of(sols).t = Create(1.0) then
        put_line("The first solution has continuation parameter t = 1.0.");
        put_line("By default, the continuation goes from t = 0.0 to 1.0.");
        put("Do you want to change t ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          put("Give real part of t (by default, type 0) : ");
          Read_Double_Float(tre);
          put("Give imaginary part of t (by default, type 0) : ");
          Read_Double_Float(tim);
          Set_Continuation_Parameter(sols,Create(tre,tim));
        end if;
      end if;
    end if;
  end Check_Continuation_Parameter;

  procedure Check_Continuation_Parameter
                ( sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    ans : character;
    tre,tim : double_float;
    dd_one : constant double_double := create(1.0);
    one : constant DoblDobl_Complex_Numbers.Complex_Number := Create(dd_one);

  begin
    if not Is_Null(sols) then
      if Head_Of(sols).t = one then
        put_line("The first solution has continuation parameter t = 1.0.");
        put_line("By default, the continuation goes from t = 0.0 to 1.0.");
        put("Do you want to change t ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          put("Give real part of t (by default, type 0) : ");
          Read_Double_Float(tre);
          put("Give imaginary part of t (by default, type 0) : ");
          Read_Double_Float(tim);
          declare
            dd_tre : constant double_double := create(tre);
            dd_tim : constant double_double := create(tim);
          begin
            Set_Continuation_Parameter(sols,Create(dd_tre,dd_tim));
          end;
        end if;
      end if;
    end if;
  end Check_Continuation_Parameter;

  procedure Check_Continuation_Parameter
                ( sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    ans : character;
    tre,tim : double_float;
    qd_zero : constant quad_double := create(0.0);
    zero : constant QuadDobl_Complex_Numbers.Complex_Number := Create(qd_zero);

  begin
    if not Is_Null(sols) then
      if not (Head_Of(sols).t = zero) then
        put_line("The first solution has continuation parameter t != 0.0.");
        put_line("By default, the continuation goes from t = 0.0 to 1.0.");
        put("Do you want to change t ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          put("Give real part of t (by default, type 0) : ");
          Read_Double_Float(tre);
          put("Give imaginary part of t (by default, type 0) : ");
          Read_Double_Float(tim);
          declare
            dd_tre : constant quad_double := create(tre);
            dd_tim : constant quad_double := create(tim);
          begin
            Set_Continuation_Parameter(sols,Create(dd_tre,dd_tim));
          end;
        end if;
      end if;
    end if;
  end Check_Continuation_Parameter;

  procedure Check_Continuation_Parameter
                ( sols : in out Multprec_Complex_Solutions.Solution_List ) is

    use Multprec_Floating_Numbers,Multprec_Complex_Numbers;
    use Multprec_Complex_Solutions;

    ans : character;
    tre,tim : double_float;
    zero : constant Multprec_Complex_Numbers.Complex_Number
         := create(integer(0));

  begin
    if not Is_Null(sols) then
      if not Equal(Head_Of(sols).t,zero) then
        put_line("The first solution has continuation parameter t != 0.0.");
        put_line("By default, the continuation goes from t = 0.0 to 1.0.");
        put("Do you want to change t ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y' then
          put("Give real part of t (by default, type 0) : ");
          Read_Double_Float(tre);
          put("Give imaginary part of t (by default, type 0) : ");
          Read_Double_Float(tim);
          declare
            mp_tre : constant Floating_Number := create(tre);
            mp_tim : constant Floating_Number := create(tim);
          begin
            Set_Continuation_Parameter(sols,Create(mp_tre,mp_tim));
          end;
        end if;
      end if;
    end if;
  end Check_Continuation_Parameter;

  procedure Read_Start_Solutions
                ( infile,outfile : in file_type;
                  qsols : out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Scans the input file infile for the solutions banner
  --   and reads the solutions into the list qsols.
  --   If there is no banner, then the user has the opportunity
  --   to provide another file name with start solutions.
  --   The continuation parameter is checked.
  --   The solution list is written to the output file outfile.

    use Standard_Complex_Solutions;

    found : boolean;

  begin
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found 
     then get(infile,qsols);
     else new_line; Read(qsols);
    end if;
    Check_Continuation_Parameter(qsols);
    put_line(outfile,"THE START SOLUTIONS : ");
    put(outfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    new_line(outfile);
  end Read_Start_Solutions;

  procedure Read_Start_Solutions
                ( infile,outfile : in file_type;
                  qsols : out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Scans the input file infile for the solutions banner
  --   and reads the solutions into the list qsols.
  --   If there is no banner, then the user has the opportunity
  --   to provide another file name with start solutions.
  --   The continuation parameter is checked.
  --   The solution list is written to the output file outfile.

    use DoblDobl_Complex_Solutions;

    found : boolean;

  begin
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found 
     then get(infile,qsols);
     else new_line; Read(qsols);
    end if;
    Check_Continuation_Parameter(qsols);
    put_line(outfile,"THE START SOLUTIONS : ");
    put(outfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    new_line(outfile);
  end Read_Start_Solutions;

  procedure Read_Start_Solutions
                ( infile,outfile : in file_type;
                  qsols : out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Scans the input file infile for the solutions banner
  --   and reads the solutions into the list qsols.
  --   If there is no banner, then the user has the opportunity
  --   to provide another file name with start solutions.
  --   The continuation parameter is checked.
  --   The solution list is written to the output file outfile.

    use QuadDobl_Complex_Solutions;

    found : boolean;

  begin
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found 
     then get(infile,qsols);
     else new_line; Read(qsols);
    end if;
    Check_Continuation_Parameter(qsols);
    put_line(outfile,"THE START SOLUTIONS : ");
    put(outfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    new_line(outfile);
  end Read_Start_Solutions;

  procedure Read_Start_Solutions
                ( infile,outfile : in file_type;
                  qsols : out Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Scans the input file infile for the solutions banner
  --   and reads the solutions into the list qsols.
  --   If there is no banner, then the user has the opportunity
  --   to provide another file name with start solutions.
  --   The continuation parameter is checked.
  --   The solution list is written to the output file outfile.

    use Multprec_Complex_Solutions;

    found : boolean;

  begin
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found 
     then get(infile,qsols);
     else new_line; Read(qsols);
    end if;
    Check_Continuation_Parameter(qsols);
    put_line(outfile,"THE START SOLUTIONS : ");
    put(outfile,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    new_line(outfile);
  end Read_Start_Solutions;

  procedure Read_Start_System
                ( file : in file_type;
                  q : out Standard_Complex_Poly_Systems.Poly_Sys;
                  qsols : out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Interactive reading of start system with start solutions.
  --   The system will be written to file.

    infile : file_type;
    readfail : boolean := true;
    nbvar,nbequ : natural32;

    procedure Try_to_Read is

    -- DESCRIPTION :
    --   Prompts the user for a start system and allows to retry
    --   in case of failures.

      m : natural32 := 0;

    begin
      put_line("Reading the name of the file for start system.");
      Read_Name_and_Open_File(infile);
      m := natural32(Scan_Line_for_Number(infile));
      get(infile,q);
      readfail := false;
    exception
      when others => put("The system on the file is not correct,");
                     put(" m = "); put(m,1); put(".");
                     readfail := true;
                     put_line("  Try again..."); close(infile);
    end Try_to_Read;

  begin
    new_line;
    loop
      Try_to_Read;
      exit when not readfail;
    end loop;
    nbvar := Number_of_Unknowns(q(q'first));
    nbequ := natural32(q'last);
    put_line(file,"THE START SYSTEM : ");
    if nbequ = nbvar 
     then put(file,nbequ,q);
     else put(file,nbequ,nbvar,q);
    end if;
    new_line(file);
    Read_Start_Solutions(infile,file,qsols);
    Close(infile);
  end Read_Start_System;

  procedure Read_Start_System
                ( file : in file_type;
                  q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                  qsols : out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Interactive reading of start system with start solutions.
  --   The system will be written to file.

    infile : file_type;
    readfail : boolean := true;
   -- nbequ,nbvar : natural32;  

    procedure Try_to_Read is

    -- DESCRIPTION :
    --   Prompts the user for a start system and allows to retry
    --   in case of failures.

      m : natural32 := 0;

    begin
      put_line("Reading the name of the file for start system.");
      Read_Name_and_Open_File(infile);
      m := natural32(Scan_Line_for_Number(infile));
      get(infile,q);
      readfail := false;
    exception
      when others => put("The system on the file is not correct,");
                     put(" m = "); put(m,1); put(".");
                     readfail := true;
                     put_line("  Try again..."); close(infile);
    end Try_to_Read;

  begin
    new_line;
    loop
      Try_to_Read;
      exit when not readfail;
    end loop;
   -- nbvar := Number_of_Unknowns(q(q'first));
   -- nbequ := natural32(q'last);
    put_line(file,"THE START SYSTEM : ");
   -- if nbequ = nbvar
   --  then
    put(file,q);
   --  else put(file,nbequ,nbvar,q);
   -- end if;
    new_line(file);
    Read_Start_Solutions(infile,file,qsols);
    Close(infile);
  end Read_Start_System;

  procedure Read_Start_System
                ( file : in file_type;
                  q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                  qsols : out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Interactive reading of start system with start solutions.
  --   The system will be written to file.

    infile : file_type;
    readfail : boolean := true;

    procedure Try_to_Read is

    -- DESCRIPTION :
    --   Prompts the user for a start system and allows to retry
    --   in case of failures.

      m : natural32 := 0;

    begin
      put_line("Reading the name of the file for start system.");
      Read_Name_and_Open_File(infile);
      m := natural32(Scan_Line_for_Number(infile));
      get(infile,q);
      readfail := false;
    exception
      when others => put("The system on the file is not correct,");
                     put(" m = "); put(m,1); put(".");
                     readfail := true;
                     put_line("  Try again..."); close(infile);
    end Try_to_Read;

  begin
    new_line;
    loop
      Try_to_Read;
      exit when not readfail;
    end loop;
    put_line(file,"THE START SYSTEM : ");
    put(file,q); new_line(file);
    Read_Start_Solutions(infile,file,qsols);
    Close(infile);
  end Read_Start_System;

  procedure Read_Start_System
                ( file : in file_type;
                  q : out Multprec_Complex_Poly_Systems.Poly_Sys;
                  qsols : out Multprec_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Interactive reading of start system with start solutions.
  --   The system will be written to file.

    infile : file_type;
    readfail : boolean := true;

    procedure Try_to_Read is

    -- DESCRIPTION :
    --   Prompts the user for a start system and allows to retry
    --   in case of failures.

      m : natural32 := 0;

    begin
      put_line("Reading the name of the file for start system.");
      Read_Name_and_Open_File(infile);
      m := natural32(Scan_Line_for_Number(infile));
      get(infile,q);
      readfail := false;
    exception
      when others => put("The system on the file is not correct,");
                     put(" m = "); put(m,1); put(".");
                     readfail := true;
                     put_line("  Try again..."); close(infile);
    end Try_to_Read;

  begin
    new_line;
    loop
      Try_to_Read;
      exit when not readfail;
    end loop;
    put_line(file,"THE START SYSTEM : ");
    put(file,q); new_line(file);
    Read_Start_Solutions(infile,file,qsols);
    Close(infile);
  end Read_Start_System;

  procedure Read_Start_System
                ( file : in file_type;
                  q : out Standard_Complex_Laur_Systems.Laur_Sys;
                  qsols : out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Interactive reading of start system with start solutions.
  --   The system will be written to file.

    infile : file_type;
    readfail : boolean := true;
    nbequ,nbvar : natural32;

    procedure Try_to_Read is

    -- DESCRIPTION :
    --   Prompts the user for a start system and allows to retry
    --   in case of failures.

      m : natural32 := 0;

    begin
      put_line("Reading the name of the file for start system.");
      Read_Name_and_Open_File(infile);
      m := natural32(Scan_Line_for_Number(infile));
      get(infile,q);
      readfail := false;
    exception
      when others => put("The system on the file is not correct,");
                     put(" m = "); put(m,1); put(".");
                     readfail := true;
                     put_line("  Try again..."); close(infile);
    end Try_to_Read;

  begin
    new_line;
    loop
      Try_to_Read;
      exit when not readfail;
    end loop;
    put_line(file,"THE START SYSTEM : ");
    nbvar := Number_of_Unknowns(q(q'first));
    nbequ := natural32(q'last);
    if nbequ = nbvar
     then put(file,nbequ,q);
     else put(file,nbequ,nbvar,q);
    end if;
    new_line(file);
    Read_Start_Solutions(infile,file,qsols);
    Close(infile);
  end Read_Start_System;

  procedure Driver_for_Polynomial_Continuation 
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys; 
                  prclvl : in natural32;
                  ls : in String_Splitters.Link_to_Array_of_Strings;
                  sols : out Standard_Complex_Solutions.Solution_List;
                  ddsols : out DoblDobl_Complex_Solutions.Solution_List;
                  qdsols : out QuadDobl_Complex_Solutions.Solution_List;
                  mpsols : out Multprec_Complex_Solutions.Solution_List;
                  target : out Complex_Number; verbose : in integer32 := 0 ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    pp,q : Poly_Sys(p'range);
    t : Complex_Number;
    qsols : Solution_List;
    mqsols : Multprec_Complex_Solutions.Solution_List;
    proj : boolean;
    deci,size : natural32 := 0;
    nbequ,nbvar : integer32;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Polynomial_Continuation 1 ...");
    end if;
    Read_Start_System(file,q,qsols);
    nbequ := q'last;
    nbvar := Head_Of(qsols).n;
    Copy(p,pp);
    if prclvl <= 1 then
      deci := 16;
    elsif prclvl = 2 then
      deci := 32;
    else
      deci := 64;
    end if;
    Driver_for_Homotopy_Construction(file,ls,pp,q,qsols,t,deci);
    proj := (Number_of_Unknowns(q(q'first)) > natural32(q'last));
    if proj
     then Ask_Symbol;
    end if;
    new_line;
    if deci <= 16 then
      if nbequ = nbvar
       then Driver_for_Standard_Continuation(file,qsols,proj,target=>t);
       else Driver_for_Standard_Continuation(file,qsols,proj,nbequ,target=>t);
      end if;
      sols := qsols;
    elsif deci <= 32 then
      ddsols := DoblDobl_Complex_Solutions.Create(qsols);
      if nbequ = nbvar
       then Driver_for_DoblDobl_Continuation(file,ddsols,target=>t);
       else Driver_for_DoblDobl_Continuation(file,ddsols,nbequ,target=>t);
      end if;
    elsif deci <= 64 then
      qdsols := QuadDobl_Complex_Solutions.Create(qsols);
      if nbequ = nbvar
       then Driver_for_QuadDobl_Continuation(file,qdsols,target=>t);
       else Driver_for_QuadDobl_Continuation(file,qdsols,nbequ,target=>t);
      end if;
    else
      mqsols := Multprec_Complex_Solutions.Create(qsols);
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Multprec_Complex_Solutions.Set_Size(mqsols,size);
      Driver_for_Multprec_Continuation(file,mqsols,proj,deci,t);
      mpsols := mqsols;
    end if;
   -- Homotopy.Clear;  --> clearing here creates difficulties for root refiner
    target := t;
  end Driver_for_Polynomial_Continuation;

  procedure Driver_for_Laurent_Continuation
                ( file : in file_type;
                  p : in Standard_Complex_Laur_Systems.Laur_Sys;
                  prclvl : in natural32;
                  ls : in String_Splitters.Link_to_Array_of_Strings;
                  sols : out Standard_Complex_Solutions.Solution_list;
                  ddsols : out DoblDobl_Complex_Solutions.Solution_list;
                  qdsols : out QuadDobl_Complex_Solutions.Solution_list;
                  mpsols : out Multprec_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    pp,q : Laur_Sys(p'range);
    t : Complex_Number;
    qsols : Solution_List;
    proj : boolean;
    deci : natural32 := 0;
    nbequ,nbvar : integer32;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Laurent_Continuation ...");
    end if;
    Read_Start_System(file,q,qsols);
    nbequ := q'last;
    nbvar := integer32(Number_of_Unknowns(q(q'first)));
    Copy(p,pp);
    if prclvl <= 1 then
      deci := 16;
    elsif prclvl = 2 then
      deci := 32;
    else
      deci := 64;
    end if;
    Driver_for_Homotopy_Construction(file,ls,pp,q,qsols,t,deci);
    proj := (Number_of_Unknowns(q(q'first)) > natural32(q'last));
    if proj
     then Ask_Symbol;
    end if;
    new_line;
    if deci <= 16 then
      if nbequ = nbvar
       then Driver_for_Standard_Laurent_Continuation(file,qsols,proj,target=>t);
       else Driver_for_Standard_Laurent_Continuation(file,qsols,proj,nbequ,t);
      end if;
    elsif deci <= 32 then
      ddsols := DoblDobl_Complex_Solutions.Create(qsols);
      if nbequ = nbvar
       then Driver_for_DoblDobl_Laurent_Continuation(file,ddsols,target=>t);
       else Driver_for_DoblDobl_Laurent_Continuation(file,ddsols,nbequ,t);
      end if;
    elsif deci <= 64 then
      qdsols := QuadDobl_Complex_Solutions.Create(qsols);
      if nbequ = nbvar
       then Driver_for_QuadDobl_Laurent_Continuation(file,qdsols,target=>t);
       else Driver_for_QuadDobl_Laurent_Continuation(file,qdsols,nbequ,t);
      end if;
    else
      mpsols := Multprec_Complex_Solutions.Create(qsols);
   --   size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
   --   Multprec_Complex_Solutions.Set_Size(mqsols,size);
   --   Driver_for_Multprec_Continuation(file,mqsols,proj,deci,t);
   --   mpsols := mqsols;
    end if;
    sols := qsols;
   -- Homotopy.Clear;  --> clearing here creates difficulties for root refiner
    target := t;
  end Driver_for_Laurent_Continuation;

  procedure Driver_for_Parameter_Continuation
                ( file : in file_type;
                  p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  k : in natural32; target : in Complex_Number;
                  sols : out Standard_Complex_Solutions.Solution_list;
                  verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
    qsols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Parameter_Continuation ...");
    end if;
    new_line; Read(qsols);
    put_line(file,"THE START SOLUTIONS :");
    put(file,Length_Of(qsols),natural32(Head_Of(qsols).n),qsols);
    new_line(file);
    Standard_Homotopy.Create(p,integer32(k));
    put_line(file,"HOMOTOPY PARAMETERS :");
    put(file,"  k : "); put(file,k,2); new_line(file);
    put(file,"  a : "); put(file,target); new_line(file);
    Driver_for_Standard_Continuation(file,qsols,false,target=>target);
   -- Homotopy.Clear; --> clearing here creates difficulties for root refiner
    sols := qsols;
  end Driver_for_Parameter_Continuation;

-- DRIVERS FOR DOUBLE DOUBLE, QUAD DOUBLE & MULTIPRECISION CONTINUATION :

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type;
                  p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                  sols : out DoblDobl_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;
    q : Poly_Sys(p'range);
    qsols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Polynomial_Continuation 2 ...");
    end if;
    Read_Start_System(file,q,qsols);
    Driver_for_Homotopy_Construction(file,p,q,target);
    Driver_for_DoblDobl_Continuation(file,qsols,target=>target);
    sols := qsols;
  end Driver_for_Polynomial_Continuation;

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type;
                  p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                  sols : out QuadDobl_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;
    q : Poly_Sys(p'range);
    qsols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Polynomial_Continuation 3 ...");
    end if;
    Read_Start_System(file,q,qsols);
    Driver_for_Homotopy_Construction(file,p,q,target);
    Driver_for_QuadDobl_Continuation(file,qsols,target=>target);
    sols := qsols;
  end Driver_for_Polynomial_Continuation;

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type; dp : in natural32;
                  p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                  sols : out Multprec_Complex_Solutions.Solution_list;
                  target : out Complex_Number; verbose : in integer32 := 0 ) is

    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Solutions;
    q : Poly_Sys(p'range);
    qsols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Polynomial_Continuation 4 ...");
    end if;
    Read_Start_System(file,q,qsols);
    Driver_for_Homotopy_Construction(file,dp,p,q,target);
    Driver_for_Multprec_Continuation(file,qsols,false,dp,target);
    sols := qsols;
  end Driver_for_Polynomial_Continuation;

-- REDEFINING ARTIFIICIAL-PARAMETER HOMOTOPIES

  procedure Standard_Redefine_Homotopy is

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Homotopy.Accessibility_Constant;
    dim : constant integer32 := Standard_Homotopy.Dimension;
    p : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Standard_Homotopy.Target_System;
    q : constant Standard_Complex_Poly_Systems.Poly_Sys(1..dim)
      := Standard_Homotopy.Start_System;
    start,target : Standard_Complex_Poly_Systems.Poly_Sys(1..dim);

  begin
    Standard_Complex_Poly_Systems.Copy(p,target);
    Standard_Complex_Poly_Systems.Copy(q,start);
    Standard_Homotopy.Clear;
    Standard_Homotopy.Create(target,start,1,gamma);
  end Standard_Redefine_Homotopy;

  procedure DoblDobl_Redefine_Homotopy is

    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Homotopy.Accessibility_Constant;
    dim : constant integer32 := DoblDobl_Homotopy.Dimension;
    p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..dim)
      := DoblDobl_Homotopy.Target_System;
    q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..dim)
      := DoblDobl_Homotopy.Start_System;
    start,target : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..dim);

  begin
    DoblDobl_Complex_Poly_Systems.Copy(p,target);
    DoblDobl_Complex_Poly_Systems.Copy(q,start);
    DoblDobl_Homotopy.Clear;
    DoblDobl_Homotopy.Create(target,start,1,gamma);
  end DoblDobl_Redefine_Homotopy;

  procedure QuadDobl_Redefine_Homotopy is

    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Homotopy.Accessibility_Constant;
    dim : constant integer32 := QuadDobl_Homotopy.Dimension;
    p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..dim)
      := QuadDobl_Homotopy.Target_System;
    q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..dim)
      := QuadDobl_Homotopy.Start_System;
    start,target : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..dim);

  begin
    QuadDobl_Complex_Poly_Systems.Copy(p,target);
    QuadDobl_Complex_Poly_Systems.Copy(q,start);
    QuadDobl_Homotopy.Clear;
    QuadDobl_Homotopy.Create(target,start,1,gamma);
  end QuadDobl_Redefine_Homotopy;

-- CALLING THE PATH TRACKERS :

  procedure Driver_for_Standard_Continuation
                ( file : in file_type;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  proj : in boolean; nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    oc : natural32;
    report : boolean;
    n : constant natural32 := natural32(Head_Of(sols).n);
    nv : constant natural32 := Length_Of(sols);
    w : Standard_Integer_Vectors.Vector(1..integer32(nv));
    v : Standard_Floating_VecVecs.Link_to_VecVec;
    errv : Standard_Floating_Vectors.Link_to_Vector;
    k : natural32;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Standard_Continuation ...");
    end if;
    new_line;
    Driver_for_Continuation_Parameters(file);
    if Continuation_Parameters.endext_order > 0
     then Init_Path_Directions(n,nv,v,errv);
    end if;
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    if Continuation_Parameters.endext_order > 0 then
      k := Standard_Homotopy.Relaxation_Power;
      --put("the original relaxation power : "); put(k,1); new_line;
      if k /= 1
       then Standard_Redefine_Homotopy;
      end if;
     -- k := Standard_Homotopy.Relaxation_Power;
     -- put("the redefined relaxation power : "); put(k,1); new_line;
      w := (w'range => 1);
      Toric_Continue(file,sols,proj,report,w,v.all,errv.all,target);
      Write_Directions(file,w,v.all,errv.all);
    else
      Continue(file,sols,proj,report,nbq,target=>target);
    end if;
  end Driver_for_Standard_Continuation;

  procedure Driver_for_Standard_Laurent_Continuation
                ( file : in file_type;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  proj : in boolean; nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 ) is

    oc : natural32;
    report : boolean;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Standard_Laurent_Continuation ...");
    end if;
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    LaurCont(file,sols,proj,report,nbq,target);
  end Driver_for_Standard_Laurent_Continuation;

  procedure Driver_for_DoblDobl_Laurent_Continuation
              ( file : in file_type;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                nbq : in integer32 := 0;
                target : in Complex_Number := Create(1.0);
                verbose : in integer32 := 0 ) is

    oc : natural32;
    report : boolean;
    dd_target : constant DoblDobl_Complex_Numbers.Complex_Number
              := Standard_to_DoblDobl_Complex(target);

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_DoblDobl_Laurent_Continuation ...");
    end if;
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    LaurCont(file,sols,report,nbq,dd_target);
  end Driver_for_DoblDobl_Laurent_Continuation;

  procedure Driver_for_QuadDobl_Laurent_Continuation
              ( file : in file_type;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                nbq : in integer32 := 0;
                target : in Complex_Number := Create(1.0);
                verbose : in integer32 := 0 ) is

    oc : natural32;
    report : boolean;
    qd_target : constant QuadDobl_Complex_Numbers.Complex_Number
              := Standard_to_QuadDobl_Complex(target);

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_QuadDobl_Laurent_Continuation ...");
    end if;
    new_line;
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    LaurCont(file,sols,report,nbq,qd_target);
  end Driver_for_QuadDobl_Laurent_Continuation;

  procedure Driver_for_Multprec_Continuation
              ( file : in file_type;
                sols : in out Multprec_Complex_Solutions.Solution_List;
                proj : in boolean; deci : in natural32;
                target : in Complex_Number := Create(1.0);
                verbose : in integer32 := 0 ) is

    oc : natural32;
    report : boolean;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_Multprec_Continuation ...");
    end if;
    new_line;
    Continuation_Parameters.Tune(0); -- ,deci); -- just leave default ...
    Driver_for_Continuation_Parameters(file);
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    Continue(file,sols,proj,report,Create(target));
  end Driver_for_Multprec_Continuation;

  procedure Driver_for_DoblDobl_Continuation
                ( file : in file_type;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List;
                  nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    oc : natural32;
    report : boolean;
    dd_targ : constant DoblDobl_Complex_Numbers.Complex_Number
            := Standard_to_DoblDobl_Complex(target);
    n : constant natural32 := natural32(Head_Of(sols).n);
    nv : constant natural32 := Length_Of(sols);
    w : Standard_Integer_Vectors.Vector(1..integer32(nv));
    v : Double_Double_VecVecs.Link_to_VecVec;
    errv : Double_Double_Vectors.Link_to_Vector;
    k : natural32;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_DoblDobl_Continuation ...");
    end if;
    new_line;
    Continuation_Parameters.Tune(0); -- ,32); -- too severe !!!
    Driver_for_Continuation_Parameters(file);
    if Continuation_Parameters.endext_order > 0
     then Init_Path_Directions(n,nv,v,errv);
    end if;
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    if Continuation_Parameters.endext_order > 0 then
      k := DoblDobl_Homotopy.Relaxation_Power;
      --put("the original relaxation power : "); put(k,1); new_line;
      if k /= 1
       then DoblDobl_Redefine_Homotopy;
      end if;
     -- k := DoblDobl_Homotopy.Relaxation_Power;
     -- put("the redefined relaxation power : "); put(k,1); new_line;
      w := (w'range => 1);
      Toric_Continue(file,sols,false,report,w,v.all,errv.all,dd_targ);
      Write_Directions(file,w,v.all,errv.all);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    else
      Continue(file,sols,report,nbq,dd_targ);
    end if;
  end Driver_for_DoblDobl_Continuation;

  procedure Driver_for_QuadDobl_Continuation
                ( file : in file_type;
                  sols : in out QuadDobl_Complex_Solutions.Solution_List;
                  nbq : in integer32 := 0;
                  target : in Complex_Number := Create(1.0);
                  verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    oc : natural32;
    report : boolean;
    qd_targ : constant QuadDobl_Complex_Numbers.Complex_Number
             := Standard_to_QuadDobl_Complex(target);
    n : constant natural32 := natural32(Head_Of(sols).n);
    nv : constant natural32 := Length_Of(sols);
    w : Standard_Integer_Vectors.Vector(1..integer32(nv));
    v : Quad_Double_VecVecs.Link_to_VecVec;
    errv : Quad_Double_Vectors.Link_to_Vector;
    k : natural32;

  begin
    if verbose > 0 then
      put("-> in drivers_for_poly_continuation.");
      put_line("Driver_for_QuadDobl_Continuation ...");
    end if;
    new_line;
    Continuation_Parameters.Tune(0); -- ,64); -- too severe !!!
    Driver_for_Continuation_Parameters(file);
    if Continuation_Parameters.endext_order > 0
     then Init_Path_Directions(n,nv,v,errv);
    end if;
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    if Continuation_Parameters.endext_order > 0 then
      w := (w'range => 1);
      k := QuadDobl_Homotopy.Relaxation_Power;
      --put("the original relaxation power : "); put(k,1); new_line;
      if k /= 1
       then QuadDobl_Redefine_Homotopy;
      end if;
     -- k := QuadDobl_Homotopy.Relaxation_Power;
     -- put("the redefined relaxation power : "); put(k,1); new_line;
      Toric_Continue(file,sols,false,report,w,v.all,errv.all,qd_targ);
      Write_Directions(file,w,v.all,errv.all);
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    else
      Continue(file,sols,report,nbq,qd_targ);
    end if;
  end Driver_for_QuadDobl_Continuation;

end Main_Poly_Continuation;
