with Interfaces.C;
with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Homotopy;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Solutions;
with DoblDobl_Homotopy;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Solutions;
with QuadDobl_Homotopy;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Drivers_to_Series_Trackers;
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
  end Job0;

  function Job1 return integer32 is -- clear parameter values
  begin
    Homotopy_Continuation_Parameters.Destruct;
    return 0;
  end Job1;

  function Job2 return integer32 is -- get a value

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    idx : constant natural32 := natural32(v_a(v_a'first));
    fail : integer32;
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
      when  8 => Assign(homconpars.alpha,c);
      when  9 => Assign(homconpars.tolres,c);
      when 10 => Assign(homconpars.epsilon,c);
      when 11 => Assign(integer32(homconpars.corsteps),b);
      when 12 => Assign(integer32(homconpars.maxsteps),b);
      when others => 
        put_line("Index value for the parameter is out of range.");
    end case;
    return 0;
  end Job2;

  function Job3 return integer32 is -- set a value

    fail : integer32;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
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
      when  2 => v_b := C_intarrs.Value(b);
                 homconpars.numdeg := natural32(v_b(v_b'first));
      when  3 => v_b := C_intarrs.Value(b);
                 homconpars.dendeg := natural32(v_b(v_b'first));
      when  4 => v_c := C_dblarrs.Value(c);
                 homconpars.maxsize := double_float(v_c(v_c'first));
      when  5 => v_c := C_dblarrs.Value(c);
                 homconpars.minsize := double_float(v_c(v_c'first));
      when  6 => v_c := C_dblarrs.Value(c);
                 homconpars.sbeta := double_float(v_c(v_c'first));
      when  7 => v_c := C_dblarrs.Value(c);
                 homconpars.pbeta := double_float(v_c(v_c'first));
      when  8 => v_c := C_dblarrs.Value(c);
                 homconpars.alpha := double_float(v_c(v_c'first));
      when  9 => v_c := C_dblarrs.Value(c);
                 homconpars.tolres := double_float(v_c(v_c'first));
      when 10 => v_c := C_dblarrs.Value(c);
                 homconpars.epsilon := double_float(v_c(v_c'first));
      when 11 => v_b := C_intarrs.Value(b);
                 homconpars.corsteps := natural32(v_b(v_b'first));
      when 12 => v_b := C_intarrs.Value(b);
                 homconpars.maxsteps := natural32(v_b(v_b'first));
      when others =>
        put_line("Index value for the parameter is out of range.");
    end case;
    return 0;
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

    homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters
               := Homotopy_Continuation_Parameters.Retrieve;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    PHCpack_Operations.Retrieve_Target_System(target);
    Standard_Homotopy.Create(target.all,start.all,tpow,homconpars.gamma);
    if name = "" then
      Drivers_to_Series_Trackers.Standard_Track
        (target'last,sols,homconpars.all); --,verbose);
    else
      Create(file,out_file,name);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.Standard_Track
        (file,target'last,sols,homconpars.all,verbose);
    end if;
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

    homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters
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
      Drivers_to_Series_Trackers.DoblDobl_Track
        (target'last,sols,homconpars.all); -- ,verbose);
    else
      Create(file,out_file,name);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.DoblDobl_Track
        (file,target'last,sols,homconpars.all,verbose);
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

    homconpars : Homotopy_Continuation_Parameters.Link_to_Parameters
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
      Drivers_to_Series_Trackers.QuadDobl_Track
        (target'last,sols,homconpars.all); -- ,verbose);
    else
      Create(file,out_file,name);
      Homotopy_Continuation_Parameters_io.put(file,homconpars.all);
      Drivers_to_Series_Trackers.QuadDobl_Track
        (file,target'last,sols,homconpars.all,verbose);
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
        v_b : C_Integer_Array(0..Interfaces.C.size_T(nbc-1))
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
  end Job4;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- default values
      when 1 => return Job1; -- clear parameter values
      when 2 => return Job2; -- get value
      when 3 => return Job3; -- set value
      when 4 => return Job4; -- track paths
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
