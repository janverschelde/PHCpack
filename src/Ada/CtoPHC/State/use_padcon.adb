with Interfaces.C;
with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Homotopy;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with DoblDobl_Homotopy;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with QuadDobl_Homotopy;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Standard_Pade_Approximants;
with DoblDobl_Pade_Approximants;
with QuadDobl_Pade_Approximants;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Drivers_to_Series_Trackers;
with Standard_SeriesPade_Tracker;
with DoblDobl_SeriesPade_Tracker;
with QuadDobl_SeriesPade_Tracker;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;

function use_padcon ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is -- set default values

    pars : constant Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;

  begin
   -- put_line("setting default values for parameters ...");
    Homotopy_Continuation_Parameters.Construct(pars);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 0 of use_padcon.");
      return 735;
  end Job0;

  function Job1 return integer32 is -- clear parameter values
  begin
    Homotopy_Continuation_Parameters.Destruct;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 1 of use_padcon.");
      return 736;
  end Job1;

  function Job2 return integer32 is -- get a value

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    idx : constant natural32 := natural32(v_a(v_a'first));
    fail : integer32 := 0;
    homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    use Homotopy_Continuation_Parameters;

  begin
    if homconpars = null then
      fail := Job0; -- initialize with default values
      homconpars := Homotopy_Continuation_Parameters.Retrieve;
    end if;
    case idx is
      when  1 => -- put(" 1. gamma : "); put(homconpars.gamma); new_line; 
                 Assign(homconpars.gamma,c);
      when  2 => Assign(integer32(homconpars.numdeg),b);
      when  3 => Assign(integer32(homconpars.dendeg),b);
      when  4 => Assign(homconpars.maxsize,c);
      when  5 => Assign(homconpars.minsize,c);
      when  6 => Assign(homconpars.sbeta,c);
      when  7 => Assign(homconpars.pbeta,c);
      when  8 => Assign(homconpars.cbeta,c);
      when  9 => Assign(homconpars.alpha,c);
      when 10 => Assign(homconpars.tolres,c);
      when 11 => Assign(homconpars.epsilon,c);
      when 12 => Assign(integer32(homconpars.corsteps),b);
      when 13 => Assign(integer32(homconpars.maxsteps),b);
      when others => 
        put_line("Index value for the parameter is out of range.");
        fail := 737;
    end case;
    return fail;
  exception
    when others =>
      put_line("Exception raised in job 2 of use_padcon.");
      return 737;
  end Job2;

  function Job3 return integer32 is -- set a value

    fail : integer32 := 0;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    idx : constant natural32 := natural32(v_a(v_a'first));
    v_b : C_Integer_Array(0..0);
    v_c : C_Double_Array(0..0);
    v_gamma : C_Double_Array(0..1);
    regamma,imgamma : double_float;

    homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    use Interfaces.C;
    use Homotopy_Continuation_Parameters;

  begin
    if homconpars = null then
      fail := Job0; -- initialize with default values
      homconpars := Homotopy_Continuation_Parameters.Retrieve;
    end if;
    case idx is
      when  1 =>
        v_gamma := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
        regamma := double_float(v_gamma(v_gamma'first));
        imgamma := double_float(v_gamma(v_gamma'first+1));
        homconpars.gamma := Standard_Complex_Numbers.Create(regamma,imgamma);
      when  2 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.numdeg := natural32(v_b(v_b'first));
      when  3 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.dendeg := natural32(v_b(v_b'first));
      when  4 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.maxsize := double_float(v_c(v_c'first));
      when  5 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.minsize := double_float(v_c(v_c'first));
      when  6 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.sbeta := double_float(v_c(v_c'first));
      when  7 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.pbeta := double_float(v_c(v_c'first));
      when  8 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.cbeta := double_float(v_c(v_c'first));
      when  9 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.alpha := double_float(v_c(v_c'first));
      when 10 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.tolres := double_float(v_c(v_c'first));
      when 11 => v_c := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
                 homconpars.epsilon := double_float(v_c(v_c'first));
      when 12 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.corsteps := natural32(v_b(v_b'first));
      when 13 => v_b := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
                 homconpars.maxsteps := natural32(v_b(v_b'first));
      when others =>
        put_line("Index value for the parameter is out of range.");
        fail := 738;
    end case;
    return fail;
  exception
    when others =>
      put_line("Exception raised in job 3 of use_padcon.");
      return 738;
  end Job3;

  procedure Standard_Track ( name : in string; verbose : in boolean ) is

  -- DESCRIPTION :
  --   Tracks the solution paths in standard precision,
  --   writing no output if name is empty,
  --   otherwise, creates an output file with the given name.

    file : file_type;
    start,target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
   -- put_line("Retrieving the data from PHCpack_Operations ...");
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- put_line("Entering the creation of the homotopy ...");
   -- declare
   --   cp,cq : Standard_Complex_Poly_Systems.Poly_Sys(target'range);
   -- begin
   --   Standard_Complex_Poly_Systems.Copy(target.all,cp);
   --   Standard_Complex_Poly_Systems.Copy(start.all,cq);
     -- Standard_Homotopy.Clear;
   --   Standard_Homotopy.Create(cp,cq,tpow,homconpars.gamma);
    Standard_Homotopy.Create(target.all,start.all,tpow,homconpars.gamma);
   -- end;
   -- put_line("Done with the creation of the homotopy, call trackers ...");
    if name = "" then
      if not verbose then
        Drivers_to_Series_Trackers.Standard_Track
          (target'last,sols,homconpars.all); --,verbose);
      else
        Homotopy_Continuation_Parameters_io.put(homconpars.all);
        Drivers_to_Series_Trackers.Standard_Track
          (standard_output,target'last,sols,homconpars.all,verbose);
      end if;
    else
      Create(file,out_file,name);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.Standard_Track
        (file,target'last,sols,homconpars.all,verbose);
      close(file);
    end if;
   -- put_line("Clearing the solutions container ...");
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(sols);
  end Standard_Track;

  procedure DoblDobl_Track ( name : in string; verbose : in boolean ) is

  -- DESCRIPTION :
  --   Tracks the solution paths in standard precision,
  --   writing no output if name is empty,
  --   otherwise, creates an output file with the given name.

    use DoblDobl_Complex_Numbers_cv;

    file : file_type;
    start,target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := homconpars.gamma;
    dd_gamma : constant DoblDobl_Complex_Numbers.Complex_Number
             := Standard_to_DoblDobl_Complex(gamma);

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    PHCpack_Operations.Retrieve_Target_System(target);
    DoblDobl_Homotopy.Create(target.all,start.all,tpow,dd_gamma);
    if name = "" then
      if not verbose then
        Drivers_to_Series_Trackers.DoblDobl_Track
          (target'last,sols,homconpars.all); -- ,verbose);
      else
        Homotopy_Continuation_Parameters_io.put(homconpars.all);
        Drivers_to_Series_Trackers.DoblDobl_Track
          (standard_output,target'last,sols,homconpars.all,verbose);
      end if;
    else
      Create(file,out_file,name);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.DoblDobl_Track
        (file,target'last,sols,homconpars.all,verbose);
      close(file);
    end if;
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(sols);
  end DoblDobl_Track;

  procedure QuadDobl_Track ( name : in string; verbose : in boolean ) is

  -- DESCRIPTION :
  --   Tracks the solution paths in standard precision,
  --   writing no output if name is empty,
  --   otherwise, creates an output file with the given name.

    use QuadDobl_Complex_Numbers_cv;

    file : file_type;
    start,target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    tpow : constant natural32 := 2;

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

    gamma : constant Standard_Complex_Numbers.Complex_Number
          := homconpars.gamma;
    qd_gamma : constant QuadDobl_Complex_Numbers.Complex_Number
             := Standard_to_QuadDobl_Complex(gamma);

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    PHCpack_Operations.Retrieve_Target_System(target);
    QuadDobl_Homotopy.Create(target.all,start.all,tpow,qd_gamma);
    if name = "" then
      if not verbose then
        Drivers_to_Series_Trackers.QuadDobl_Track
          (target'last,sols,homconpars.all); -- ,verbose);
      else
        Homotopy_Continuation_Parameters_io.put(homconpars.all);
        Drivers_to_Series_Trackers.QuadDobl_Track
          (standard_output,target'last,sols,homconpars.all,verbose);
      end if;
    else
      Create(file,out_file,name);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.QuadDobl_Track
        (file,target'last,sols,homconpars.all,verbose);
      close(file);
    end if;
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(sols);
  end QuadDobl_Track;

  function Job4 return integer32 is -- track paths

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    prc : constant natural32 := natural32(v_a(v_a'first));
    nbc : constant natural32 := natural32(v_a(v_a'first+1));
    vrb : constant natural32 := natural32(v_a(v_a'first+2));
    verbose : constant boolean := (vrb > 0);

  begin
    if nbc = 0 then
      case prc is
        when 0 => Standard_Track("",verbose);
        when 1 => DoblDobl_Track("",verbose);
        when 2 => QuadDobl_Track("",verbose);
        when others => null;
      end case;
    else
      declare
        v_b : constant C_Integer_Array(0..Interfaces.C.size_T(nbc-1))
            := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
        name : constant string := C_Integer_Array_to_String(nbc,v_b);
      begin
        case prc is
          when 0 => Standard_Track(name,verbose);
          when 1 => DoblDobl_Track(name,verbose);
          when 2 => QuadDobl_Track(name,verbose);
          when others => null;
        end case;
      end;
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 4 of use_padcon.");
      return 739;
  end Job4;

  procedure Standard_Initialize_Tracker is

  -- DESCRIPTION :
  --   Retrieves target and start system and initializes the
  --   Series-Pade tracker in standard double precision.

    start,target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    Standard_SeriesPade_Tracker.Init(target,start);
  end Standard_Initialize_Tracker;

  procedure DoblDobl_Initialize_Tracker is

  -- DESCRIPTION :
  --   Retrieves target and start system and initializes the
  --   Series-Pade tracker in double double precision.

    start,target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    DoblDobl_SeriesPade_Tracker.Init(target,start);
  end DoblDobl_Initialize_Tracker;

  procedure QuadDobl_Initialize_Tracker is

  -- DESCRIPTION :
  --   Retrieves target and start system and initializes the
  --   Series-Pade tracker in quad double precision.

    start,target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    QuadDobl_SeriesPade_Tracker.Init(target,start);
  end QuadDobl_Initialize_Tracker;

  procedure Standard_Initialize_Tracker ( idx : in integer32 ) is

  -- DESCRIPTION :
  --   Retrieves the target system and initializes the
  --   Series-Pade tracker in standard double precision,
  --   with the continuation parameter defined in idx.

    target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(target);
    Standard_SeriesPade_Tracker.Init(target,idx);
  end Standard_Initialize_Tracker;

  procedure DoblDobl_Initialize_Tracker ( idx : in integer32 ) is

  -- DESCRIPTION :
  --   Retrieves the target system and initializes the
  --   Series-Pade tracker in double double precision,
  --   with the continuation parameter defined in idx.

    target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(target);
    DoblDobl_SeriesPade_Tracker.Init(target,idx);
  end DoblDobl_Initialize_Tracker;

  procedure QuadDobl_Initialize_Tracker ( idx : in integer32 ) is

  -- DESCRIPTION :
  --   Retrieves the target system and initializes the
  --   Series-Pade tracker in quad double precision,
  --   with the continuation parameter defined in idx.

    target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(target);
    QuadDobl_SeriesPade_Tracker.Init(target,idx);
  end QuadDobl_Initialize_Tracker;

  function Job5 return integer32 is -- initialize seriespade tracker

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
    case prc is
      when 0 =>
        if verbose then
          put("Initializing homotopy in Series-Pade tracker ");
          put("in double precision ...");
        end if;
        Standard_SeriesPade_Tracker.Init(homconpars.all);
        Standard_Initialize_Tracker;
      when 1 =>
        if verbose then
          put("Initializing homotopy in Series-Pade tracker ");
          put_line("in double double precision ...");
        end if;
        DoblDobl_SeriesPade_Tracker.Init(homconpars.all);
        DoblDobl_Initialize_Tracker;
      when 2 =>
        if verbose then
          put("Initializing homotopy in Series-Pade tracker ");
          put_line("in quad double precision ...");
        end if;
        QuadDobl_SeriesPade_Tracker.Init(homconpars.all);
        QuadDobl_Initialize_Tracker;
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    if verbose
     then put_line(" done!");
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 5 of use_padcon.");
      return 860;
  end Job5;

  procedure Standard_Initialize_Solution
              ( idx : in natural32; verbose : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Retrieves the solution at index idx from the solutions container
  --   in standard double precision and writes output if verbose.

    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    if verbose
     then put("initializing solution "); put(idx,1); put_line(" :");
    end if;
    Standard_Solutions_Container.Retrieve(idx,ls,fail);
    Standard_SeriesPade_Tracker.Init(ls);
    if verbose
     then put(ls.all); new_line;
    end if;
  end Standard_Initialize_Solution;

  procedure DoblDobl_Initialize_Solution
              ( idx : in natural32; verbose : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Retrieves the solution at index idx from the solutions container
  --   in double double precision and writes output if verbose.

    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    if verbose
     then put("initializing solution "); put(idx,1); put_line(" :");
    end if;
    DoblDobl_Solutions_Container.Retrieve(idx,ls,fail);
    DoblDobl_SeriesPade_Tracker.Init(ls);
    if verbose
     then put(ls.all); new_line;
    end if;
  end DoblDobl_Initialize_Solution;

  procedure QuadDobl_Initialize_Solution
              ( idx : in natural32; verbose : in boolean;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Retrieves the solution at index idx from the solutions container
  --   in quad double precision and writes output if verbose.
  --   The failure code of the retrieval is returned in fail.

    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    if verbose
     then put("initializing solution "); put(idx,1); put_line(" :");
    end if;
    QuadDobl_Solutions_Container.Retrieve(idx,ls,fail);
    QuadDobl_SeriesPade_Tracker.Init(ls);
    if verbose
     then put(ls.all); new_line;
    end if;
  end QuadDobl_Initialize_Solution;

  function Job6 return integer32 is -- initialize next start solution

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    prc : constant natural32 := natural32(v_a(v_a'first));
    idx : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    fail : boolean;

  begin
    case prc is
      when 0 => Standard_Initialize_Solution(idx,verbose,fail);
      when 1 => DoblDobl_Initialize_Solution(idx,verbose,fail);
      when 2 => QuadDobl_Initialize_Solution(idx,verbose,fail);
      when others => put_line("Wrong value for the precision.");
    end case;
    if fail then
      if verbose
       then put_line("The initialization of the start solution failed.");
      end if;
      assign(1,a);
    else
      if verbose
       then put_line("The initialization of the start solution succeeded.");
      end if;
      assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 6 of use_padcon.");
      return 861;
  end Job6;

  function Job7 return integer32 is -- run next predict-correct step

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    fail : boolean;

  begin
    case prc is
      when 0 => Standard_SeriesPade_Tracker.Predict_and_Correct(fail,verbose);
      when 1 => DoblDobl_SeriesPade_Tracker.Predict_and_Correct(fail,verbose);
      when 2 => QuadDobl_SeriesPade_Tracker.Predict_and_Correct(fail,verbose);
      when others => put_line("Wrong value for the precision.");
    end case;
    if fail then
      if verbose
       then put_line("The predict-correct step failed.");
      end if;
      assign(1,a);
    else
      if verbose
       then put_line("The predict-correct step succeeded.");
      end if;
      assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 7 of use_padcon.");
      return 862;
  end Job7;

  function Job8 return integer32 is -- get the current solution

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    prc : constant natural32 := natural32(v_a(v_a'first));
    idx : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_ls : Standard_Complex_Solutions.Link_to_Solution;
    dd_ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    qd_ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if verbose then
      put("Retrieving the current solution, at index "); put(idx,1);
      put_line(" ...");
    end if;
    case prc is
      when 0 =>
        st_ls := Standard_SeriesPade_Tracker.Get_Current_Solution;
        Standard_Solutions_Container.Replace(idx,st_ls,fail);
      when 1 =>
        dd_ls := DoblDobl_SeriesPade_Tracker.Get_Current_Solution;
        DoblDobl_Solutions_Container.Replace(idx,dd_ls,fail);
      when 2 =>
        qd_ls := QuadDobl_SeriesPade_Tracker.Get_Current_Solution;
        QuadDobl_Solutions_Container.Replace(idx,qd_ls,fail);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    if fail then
      if verbose
       then put_line("Placement of the current solution failed.");
      end if;
      assign(1,a);
    else
      if verbose
       then put_line("Placement of the current solution succeeded.");
      end if;
      assign(0,a);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 8 of use_padcon.");
      return 863;
  end Job8;

  function Job9 return integer32 is -- clear data

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));

  begin
    case prc is
      when 0 => Standard_SeriesPade_Tracker.Clear;
      when 1 => DoblDobl_SeriesPade_Tracker.Clear;
      when 2 => QuadDobl_SeriesPade_Tracker.Clear;
      when others => put_line("Wrong value for the precision.");
    end case;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 9 of use_padcon.");
      return 864;
  end Job9;

  function Job10 return integer32 is -- get pole radius

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    frp : double_float;

  begin
    case prc is
      when 0 =>
        frp := Standard_SeriesPade_Tracker.Get_Current_Pole_Radius;
      when 1 =>
        frp := hi_part(DoblDobl_SeriesPade_Tracker.Get_Current_Pole_Radius);
      when 2 =>
        frp := hihi_part(QuadDobl_SeriesPade_Tracker.Get_Current_Pole_Radius);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    Assign(frp,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 10 of use_padcon.");
      return 865;
  end Job10;

  function Job11 return integer32 is -- get closest pole

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    st_cfp : Standard_Complex_Numbers.Complex_Number;
    dd_cfp : DoblDobl_Complex_Numbers.Complex_Number;
    qd_cfp : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    case prc is
      when 0 =>
        st_cfp := Standard_SeriesPade_Tracker.Get_Current_Closest_Pole;
      when 1 =>
        dd_cfp := DoblDobl_SeriesPade_Tracker.Get_Current_Closest_Pole;
        st_cfp := DoblDobl_Complex_to_Standard(dd_cfp);
      when 2 =>
        qd_cfp := QuadDobl_SeriesPade_Tracker.Get_Current_Closest_Pole;
        st_cfp := QuadDobl_Complex_to_Standard(qd_cfp);
      when others =>
        put_line("Wrong value for the precision.");
    end case;
    Assign(st_cfp,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 11 of use_padcon.");
      return 866;
  end Job11;

  function Job12 return integer32 is -- get current t value

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    tval : double_float;

  begin
    case prc is
      when 0 => tval := Standard_SeriesPade_Tracker.Get_Current_t_Value;
      when 1 => tval := DoblDobl_SeriesPade_Tracker.Get_Current_t_Value;
      when 2 => tval := QuadDobl_SeriesPade_Tracker.Get_Current_t_Value;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(tval,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 12 of use_padcon.");
      return 867;
  end Job12;

  function Job13 return integer32 is -- get current step size

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Step_Size;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Step_Size;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Step_Size;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 13 of use_padcon.");
      return 868;
  end Job13;

  function Job19 return integer32 is -- get current series step size

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Series_Step;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Series_Step;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Series_Step;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 19 of use_padcon.");
      return 885;
  end Job19;

  function Job20 return integer32 is -- get current pole step size

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Pole_Step;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Pole_Step;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Pole_Step;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 20 of use_padcon.");
      return 886;
  end Job20;

  function Job21 return integer32 is -- get current estimated distance

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    dist : double_float;
    dd_dist : double_double;
    qd_dist : quad_double;

  begin
    case prc is
      when 0 =>
        dist := Standard_SeriesPade_Tracker.Get_Current_Estimated_Distance;
      when 1 =>
        dd_dist := DoblDobl_SeriesPade_Tracker.Get_Current_Estimated_Distance;
        dist := hi_part(dd_dist);
      when 2 => 
        qd_dist := QuadDobl_SeriesPade_Tracker.Get_Current_Estimated_Distance;
        dist := hihi_part(qd_dist);
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(dist,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 21 of use_padcon.");
      return 887;
  end Job21;

  function Job22 return integer32 is -- get current Hessian step size

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    step : double_float;

  begin
    case prc is
      when 0 => step := Standard_SeriesPade_Tracker.Get_Current_Hessian_Step;
      when 1 => step := DoblDobl_SeriesPade_Tracker.Get_Current_Hessian_Step;
      when 2 => step := QuadDobl_SeriesPade_Tracker.Get_Current_Hessian_Step;
      when others => put_line("Wrong value for the precision.");
    end case;
    Assign(step,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 22 of use_padcon.");
      return 888;
  end Job22;

  function Standard_Series_Coefficient
             ( leadidx,cffidx : integer32; verbose : boolean )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current series vectors in standard double precision.
  --   If verbose, then extra output is written to screen.

    res : Standard_Complex_Numbers.Complex_Number;
    srv : constant Standard_Complex_Series_Vectors.Link_to_Vector
        := Standard_SeriesPade_Tracker.Get_Current_Series_Vector;

  begin
    res := srv(leadidx).cff(cffidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end Standard_Series_Coefficient;

  function DoblDobl_Series_Coefficient
             ( leadidx,cffidx : integer32; verbose : boolean )
             return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current series vectors in double double precision.
  --   If verbose, then extra output is written to screen.

    res : DoblDobl_Complex_Numbers.Complex_Number;
    srv : constant DoblDobl_Complex_Series_Vectors.Link_to_Vector
        := DoblDobl_SeriesPade_Tracker.Get_Current_Series_Vector;

  begin
    res := srv(leadidx).cff(cffidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end DoblDobl_Series_Coefficient;

  function QuadDobl_Series_Coefficient
             ( leadidx,cffidx : integer32; verbose : boolean )
             return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current series vectors in quad double precision.
  --   If verbose, then extra output is written to screen.

    res : QuadDobl_Complex_Numbers.Complex_Number;
    srv : constant QuadDobl_Complex_Series_Vectors.Link_to_Vector
        := QuadDobl_SeriesPade_Tracker.Get_Current_Series_Vector;

  begin
    res := srv(leadidx).cff(cffidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end QuadDobl_Series_Coefficient;

  function Job14 return integer32 is -- get series coefficient

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    prc : constant natural32 := natural32(v_a(v_a'first));
    leadidx : constant integer32 := integer32(v_a(v_a'first+1));
    cffidx : constant integer32 := integer32(v_a(v_a'first+2));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_cff : Standard_Complex_Numbers.Complex_Number;
    dd_cff : DoblDobl_Complex_Numbers.Complex_Number;
    qd_cff : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    if verbose then
      put("Precision : ");
      case prc is
        when 0 => put("double");
        when 1 => put("double double");
        when 2 => put("quad double");
        when others => put("invalid!");
      end case;
      put("  lead idx : "); put(leadidx,1);
      put("  cff idx : "); put(cffidx,1); new_line;
    end if;
    case prc is
      when 0 => st_cff := Standard_Series_Coefficient(leadidx,cffidx,verbose);
      when 1 => dd_cff := DoblDobl_Series_Coefficient(leadidx,cffidx,verbose);
                st_cff := DoblDobl_Complex_to_Standard(dd_cff);
      when 2 => qd_cff := QuadDobl_Series_Coefficient(leadidx,cffidx,verbose);
                st_cff := QuadDobl_Complex_to_Standard(qd_cff);
      when others => put_line("Invalid precision.");
    end case;
    Assign(st_cff,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 14 of use_padcon.");
      return 869;
  end Job14;

  function Standard_Pade_Coefficient
             ( leadidx,cffidx : integer32; num : natural32;
               verbose : boolean )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current Pade vectors in standard double precision.
  --   If num = 0, then the denominator coefficient will be returned,
  --   otherwise, the numerator coefficient will be returned.
  --   If verbose, then extra output is written to screen.

    res : Standard_Complex_Numbers.Complex_Number;
    pv : constant Standard_Pade_Approximants.Link_to_Pade_Vector
       := Standard_SeriesPade_Tracker.Get_Current_Pade_Vector;

    use Standard_Pade_Approximants;

  begin
    if num = 0 then
      declare
        c : constant Standard_Complex_Vectors.Vector
          := Denominator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    else
      declare
        c : constant Standard_Complex_Vectors.Vector
          := Numerator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    end if;
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end Standard_Pade_Coefficient;

  function DoblDobl_Pade_Coefficient
             ( leadidx,cffidx : integer32; num : natural32;
               verbose : boolean )
             return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current Pade vectors in double double precision.
  --   If num = 0, then the denominator coefficient will be returned,
  --   otherwise, the numerator coefficient will be returned.
  --   If verbose, then extra output is written to screen.

    res : DoblDobl_Complex_Numbers.Complex_Number;
    pv : constant DoblDobl_Pade_Approximants.Link_to_Pade_Vector
       := DoblDobl_SeriesPade_Tracker.Get_Current_Pade_Vector;

    use DoblDobl_Pade_Approximants;

  begin
    if num = 0 then
      declare
        c : constant DoblDobl_Complex_Vectors.Vector
          := Denominator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    else
      declare
        c : constant DoblDobl_Complex_Vectors.Vector
          := Numerator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    end if;
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end DoblDobl_Pade_Coefficient;

  function QuadDobl_Pade_Coefficient
             ( leadidx,cffidx : integer32; num : natural32;
               verbose : boolean )
             return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the coefficient at position cffidx in component leadidx
  --   of the current Pade vectors in quad double precision.
  --   If num = 0, then the denominator coefficient will be returned,
  --   otherwise, the numerator coefficient will be returned.
  --   If verbose, then extra output is written to screen.

    res : QuadDobl_Complex_Numbers.Complex_Number;
    pv : constant QuadDobl_Pade_Approximants.Link_to_Pade_Vector
       := QuadDobl_SeriesPade_Tracker.Get_Current_Pade_Vector;

    use QuadDobl_Pade_Approximants;

  begin
    if num = 0 then
      declare
        c : constant QuadDobl_Complex_Vectors.Vector
          := Denominator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    else
      declare
        c : constant QuadDobl_Complex_Vectors.Vector
          := Numerator_Coefficients(pv(leadidx));
      begin
        res := c(cffidx);
      end;
    end if;
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end QuadDobl_Pade_Coefficient;

  function Job15 return integer32 is -- get Pade coefficient

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    prc : constant natural32 := natural32(v_a(v_a'first));
    num : constant natural32 := natural32(v_a(v_a'first+1));
    leadidx : constant integer32 := integer32(v_a(v_a'first+2));
    cffidx : constant integer32 := integer32(v_a(v_a'first+3));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_cff : Standard_Complex_Numbers.Complex_Number;
    dd_cff : DoblDobl_Complex_Numbers.Complex_Number;
    qd_cff : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    if verbose then
      put("Precision : ");
      case prc is
        when 0 => put("double");
        when 1 => put("double double");
        when 2 => put("quad double");
        when others => put("invalid!");
      end case;
      put("  lead idx : "); put(leadidx,1);
      put("  cff idx : "); put(cffidx,1);
      if num = 0
       then put_line("  denominator");
       else put_line("  numernator");
      end if;
    end if;
    case prc is
      when 0 =>
        st_cff := Standard_Pade_Coefficient(leadidx,cffidx,num,verbose);
      when 1 =>
        dd_cff := DoblDobl_Pade_Coefficient(leadidx,cffidx,num,verbose);
        st_cff := DoblDobl_Complex_to_Standard(dd_cff);
      when 2 =>
        qd_cff := QuadDobl_Pade_Coefficient(leadidx,cffidx,num,verbose);
        st_cff := QuadDobl_Complex_to_Standard(qd_cff);
      when others =>
        put_line("Invalid value for the precision.");
    end case;
    Assign(st_cff,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 15 of use_padcon.");
      return 870;
  end Job15;

  function Standard_Pole
             ( leadidx,poleidx : integer32; verbose : boolean )
             return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the pole at position cffidx in component leadidx
  --   of the current poles in standard double precision.
  --   If verbose, then extra output is written to screen.

    res : Standard_Complex_Numbers.Complex_Number;
    poles: constant Standard_Complex_VecVecs.Link_to_VecVec
        := Standard_SeriesPade_Tracker.Get_Current_Poles;

  begin
    res := poles(leadidx)(poleidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end Standard_Pole;

  function DoblDobl_Pole
             ( leadidx,poleidx : integer32; verbose : boolean )
             return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the pole at position cffidx in component leadidx
  --   of the current poles in double double precision.
  --   If verbose, then extra output is written to screen.

    res : DoblDobl_Complex_Numbers.Complex_Number;
    poles: constant DoblDobl_Complex_VecVecs.Link_to_VecVec
        := DoblDobl_SeriesPade_Tracker.Get_Current_Poles;

  begin
    res := poles(leadidx)(poleidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end DoblDobl_Pole;

  function QuadDobl_Pole
             ( leadidx,poleidx : integer32; verbose : boolean )
             return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the pole at position cffidx in component leadidx
  --   of the current poles in quad double precision.
  --   If verbose, then extra output is written to screen.

    res : QuadDobl_Complex_Numbers.Complex_Number;
    poles: constant QuadDobl_Complex_VecVecs.Link_to_VecVec
        := QuadDobl_SeriesPade_Tracker.Get_Current_Poles;

  begin
    res := poles(leadidx)(poleidx);
    if verbose
     then put("Returning "); put(res); new_line;
    end if;
    return res;
  end QuadDobl_Pole;

  function Job16 return integer32 is  -- get pole

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    prc : constant natural32 := natural32(v_a(v_a'first));
    leadidx : constant integer32 := integer32(v_a(v_a'first+1));
    poleidx : constant integer32 := integer32(v_a(v_a'first+2));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    vrb : constant natural32 := natural32(v_b(v_b'first));
    verbose : constant boolean := (vrb = 1);
    st_pole : Standard_Complex_Numbers.Complex_Number;
    dd_pole : DoblDobl_Complex_Numbers.Complex_Number;
    qd_pole : QuadDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;
    use QuadDobl_Complex_Numbers_cv;

  begin
    if verbose then
      put("Precision : ");
      case prc is
        when 0 => put("double");
        when 1 => put("double double");
        when 2 => put("quad double");
        when others => put("invalid!");
      end case;
      put("  lead idx : "); put(leadidx,1);
      put("  pole idx : "); put(poleidx,1); new_line;
    end if;
    case prc is
      when 0 => st_pole := Standard_Pole(leadidx,poleidx,verbose);
      when 1 => dd_pole := DoblDobl_Pole(leadidx,poleidx,verbose);
                st_pole := DoblDobl_Complex_to_Standard(dd_pole);
      when 2 => qd_pole := QuadDobl_Pole(leadidx,poleidx,verbose);
                st_pole := QuadDobl_Complex_to_Standard(qd_pole);
      when others => put_line("Invalid value for the precision.");
    end case;
    Assign(st_pole,c);
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 16 of use_padcon.");
      return 871;
  end Job16;

  function Job17 return integer32 is -- write parameters to defined file

    pars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
         := Homotopy_Continuation_Parameters.Retrieve;

  begin
    if PHCpack_Operations.Is_File_Defined then
      new_line(PHCpack_Operations.output_file);
      Homotopy_Continuation_Parameters_io.put
        (PHCpack_Operations.output_file,pars.all);
      text_io.flush(PHCpack_Operations.output_file);
    else
      Homotopy_Continuation_Parameters_io.put(pars.all);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 17 of use_padcon.");
      return 874;
  end Job17;

  function Job18 return integer32 is -- initializes natural parameter homotopy

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    prc : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    idx : constant integer32 := integer32(v_b(v_b'first));
    vrb : constant natural32 := natural32(v_b(v_b'first+1));
    verbose : constant boolean := (vrb = 1);

    homconpars : constant Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
    case prc is
      when 0 =>
        if verbose then
          put("Initializing parameter homotopy in Series-Pade tracker ");
          put("in double precision ...");
        end if;
        Standard_SeriesPade_Tracker.Init(homconpars.all);
        Standard_Initialize_Tracker(idx);
      when 1 =>
        if verbose then
          put("Initializing parameter homotopy in Series-Pade tracker ");
          put("in double double precision ...");
        end if;
        DoblDobl_SeriesPade_Tracker.Init(homconpars.all);
        DoblDobl_Initialize_Tracker(idx);
      when 2 =>
        if verbose then
          put("Initializing parameter homotopy in Series-Pade tracker ");
          put("in quad double precision ...");
        end if;
        QuadDobl_SeriesPade_Tracker.Init(homconpars.all);
        QuadDobl_Initialize_Tracker(idx);
      when others => null;
    end case;
    return 0;
  exception
    when others =>
      put_line("Exception raised in job 18 of use_padcon.");
      return 878;
  end Job18;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- default values
      when 1 => return Job1; -- clear parameter values
      when 2 => return Job2; -- get value
      when 3 => return Job3; -- set value
      when 4 => return Job4; -- track paths
      when 5 => return Job5; -- initialize seriespade homotopy tracker
      when 6 => return Job6; -- initialize next start solution
      when 7 => return Job7; -- run next predict-correct step
      when 8 => return Job8; -- get the current solution
      when 9 => return Job9; -- deallocates seriespade homotopy tracker data
      when 10 => return Job10; -- get the smallest forward pole radius
      when 11 => return Job11; -- get the closest pole
      when 12 => return Job12; -- get the current t value
      when 13 => return Job13; -- get the current step size
      when 14 => return Job14; -- get series coefficient
      when 15 => return Job15; -- get Pade coefficient
      when 16 => return Job16; -- get pole
      when 17 => return Job17; -- write parameters to defined output file
      when 18 => return Job18; -- initializes natural parameter homotopy
      when 19 => return Job19; -- get current series step
      when 20 => return Job20; -- get current pole step
      when 21 => return Job21; -- get current estimated distance
      when 22 => return Job22; -- get current Hessian step
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  exception
    when others => put("Exception raised in use_padcon handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_padcon;
