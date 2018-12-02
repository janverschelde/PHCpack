with Interfaces.C;
with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
--with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Homotopy_Continuation_Parameters;

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
      when  1 => v_gamma := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
                 regamma := double_float(v_gamma(v_gamma'first));
                 imgamma := double_float(v_gamma(v_gamma'first+1));
                 homconpars.gamma := Create(regamma,imgamma);
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

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0; -- default values
      when 1 => return Job1; -- clear parameter values
      when 2 => return Job2; -- get value
      when 3 => return Job3; -- set value
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  exception
    when others => put("Exception raised in use_numbtrop handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_padcon;
