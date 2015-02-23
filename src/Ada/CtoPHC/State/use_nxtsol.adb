with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;
--with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
--with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
--with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
--with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Solutions_Container;
with DoblDobl_Complex_Solutions;
--with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with DoblDobl_Solutions_Container;
with QuadDobl_Complex_Solutions;
--with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with QuadDobl_Solutions_Container;
with Multprec_Complex_Solutions;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Multprec_Solutions_Container;
with PHCpack_Operations;
with Standard_Path_Tracker;
with DoblDobl_Path_Tracker;
with QuadDobl_Path_Tracker;
with Multprec_Path_Tracker;
with Varbprec_Path_Tracker;

function use_nxtsol ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32 is
 
  function Job0 return integer32 is -- initialize standard homotopy

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    start,target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    if fixed = 1
     then Standard_Path_Tracker.Init(target,start,true);
     else Standard_Path_Tracker.Init(target,start,false);
    end if;
    return 0;
  exception
    when others => return 500;
  end Job0;
 
  function Job1 return integer32 is -- initialize dobldobl homotopy

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    start,target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    if fixed = 1
     then DoblDobl_Path_Tracker.Init(target,start,true);
     else DoblDobl_Path_Tracker.Init(target,start,false);
    end if;
    return 0;
  exception
    when others => return 501;
  end Job1;
 
  function Job2 return integer32 is -- initialize quaddobl homotopy

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    start,target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    if fixed = 1
     then QuadDobl_Path_Tracker.Init(target,start,true);
     else QuadDobl_Path_Tracker.Init(target,start,false);
    end if;
    return 0;
  exception
    when others => return 502;
  end Job2;

  function Job3 return integer32 is -- initialize standard solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("initializing solution "); put(k,1); put_line(" :");
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    Standard_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => return 503;
  end Job3;

  function Job4 return integer32 is -- initialize dobldobl solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("initializing solution "); put(k,1); put_line(" :");
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    DoblDobl_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => return 504;
  end Job4;

  function Job5 return integer32 is -- initialize quaddobl solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("initializing solution "); put(k,1); put_line(" :");
    QuadDobl_Solutions_Container.Retrieve(k,ls,fail);
    QuadDobl_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => return 505;
  end Job5;

  function Job6 return integer32 is -- next standard solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := Standard_Path_Tracker.get_next;
    Standard_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => return 506;
  end Job6;

  function Job7 return integer32 is -- next dobldobl solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := DoblDobl_Path_Tracker.get_next;
    DoblDobl_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => return 507;
  end Job7;

  function Job8 return integer32 is -- next quaddobl solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := QuadDobl_Path_Tracker.get_next;
    QuadDobl_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => return 508;
  end Job8;

  function Job12 return integer32 is -- initialize multiprecision tracker


    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    deci : constant natural32 := natural32(v_b(v_b'first));
    start,target : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    if fixed = 1
     then Multprec_Path_Tracker.Init(target,start,true,deci);
     else Multprec_Path_Tracker.Init(target,start,false,deci);
    end if;
    return 0;
  exception
    when others => return 512;
  end Job12;

  function Job13 return integer32 is -- initialize multiprecision solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Multprec_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("initializing solution "); put(k,1); put_line(" :");
    Multprec_Solutions_Container.Retrieve(k,ls,fail);
    Multprec_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => return 513;
  end Job13;

  function Job14 return integer32 is -- next multiprecision solution

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Multprec_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := Multprec_Path_Tracker.get_next;
    Multprec_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => return 514;
  end Job14;

  function Job16 return integer32 is -- init variable precision homotopy

  -- DESCRIPTION :
  --   The function expects three values in a:
  --   a[0] : whether a fixed gamma is to be used or not (1 or 0);
  --   a[1] : the total number of characters in the string b;
  --   a[2] : the number of characters in the first string that
  --   contains the target system.
  --   The parameter b holds the string representations of two systems,
  --   respectively the target and start system in the homotopy.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    fix : constant natural32 := natural32(v_a(v_a'first));
    fixed : constant boolean := (fix = 1);
    len : constant integer := integer(v_a(v_a'first+1));
    sec : constant integer := integer(v_a(v_a'first+2));
    v_b : constant C_Integer_Array(0..Interfaces.C.size_t(len))
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(len+1));
    hom : constant string(1..len)
        := C_Integer_Array_to_String(natural32(len),v_b);
    str_target : constant string(1..sec) := hom(1..sec);
    str_starts : constant string(sec+1..len) := hom(sec+1..len);
    tnq : constant natural := Count_Delimiters(str_target,';');
    snq : constant natural := Count_Delimiters(str_starts,';');
    target : constant Array_of_Strings := Split(tnq,str_target,';');
    starts : constant Array_of_Strings := Split(snq,str_starts,';');
    lp : constant Link_to_Array_of_Strings := new Array_of_Strings'(target);
    lq : constant Link_to_Array_of_Strings := new Array_of_Strings'(starts);

  begin
   -- put_line("The target system : "); put_line(str_target);
   -- put_line("The start system : "); put_line(str_starts);
   -- put_line("The list of target polynomials : ");
   -- for i in lp'range loop
   --   put_line(lp(i).all);
   -- end loop;
   -- put_line("The list of start polynomials : ");
   -- for i in lq'range loop
   --   put_line(lq(i).all);
   -- end loop;
    Varbprec_Path_Tracker.Init(lp,lq,fixed);
    return 0;
  exception
    when others => return 516;
  end Job16;

  function Job17 return integer32 is -- init variable precision solution

  -- DESCRIPTION :
  --   The function expects two values in a:
  --   a[0] : the length of the vector in b;
  --   a[1] : the number of variables in the solution represented by b.
  --   In b is the string representation of a solution in PHCpack format.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    len : constant integer := integer(v_a(v_a'first));
    dim : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array(0..Interfaces.C.size_t(len))
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(len+1));
    sol : constant string(1..len)
        := C_Integer_Array_to_String(natural32(len),v_b);
    lsl : constant Link_to_String := new string'(sol);

  begin
   -- put_line("initializing the variable precision path tracker ...");
   -- put_line("with the solution : "); put_line(lsl.all);
    Varbprec_Path_Tracker.Init(lsl,dim);
    return 0;
  exception
    when others => return 517;
  end Job17;

  function Job18 return integer32 is -- next variable precision solution

  -- DESCRIPTION :
  --   The function expects four values in a:
  --   a[0] : the wanted number of accurate decimal places;
  --   a[1] : the maximum precision to be used;
  --   a[2] : the maximum number of corrector steps;
  --   a[3] : whether intermediate output is wanted or not,
  --   where 1 is true and 0 is false.
  --   Because the size of the solution string can fluctuate as the
  --   demands on the precision vary, job 20 must be called immediately
  --   after this job to retrieve the solution.
  --   The size of the solution string is returned in b.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    want : constant natural32 := natural32(v_a(v_a'first));
    maxprc : constant natural32 := natural32(v_a(v_a'first+1));
    maxitr : constant natural32 := natural32(v_a(v_a'first+2));
    output : constant natural32 := natural32(v_a(v_a'first+3));
    otp : constant boolean := (output = 1);
    sol : Link_to_String;

  begin
   -- put_line("inside use_nxtsol, calling get_current ...");
    sol := Varbprec_Path_Tracker.get_current;
   -- put_line("inside use_nxtsol, calling get_next ...");
    sol := Varbprec_Path_Tracker.get_next(want,maxprc,maxitr,otp);
   -- put_line("The solution returned : ");
   -- put_line(sol.all);
    Assign(integer32(sol'last),b);
    return 0;
  exception
    when others => put_line("some exception occurred"); raise; return 518;
  end Job18;

  function Job20 return integer32 is -- current variable precision solution

  -- DESCRIPTION :
  --   Returns the current variable precision solution as a string,
  --   as computed previously by Job 18.
  --   Assign in a[0] the number of characters in the string b.

    sol : Link_to_String;

  begin
    sol := Varbprec_Path_Tracker.get_current;
    declare
      sv : constant Standard_Integer_Vectors.Vector
         := String_to_integer_Vector(sol.all);
    begin
      Assign(sv'last,a);
      Assign(sv,b);
    end;
    return 0;
  exception
    when others => return 520;
  end Job20;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Job0;   -- initialize standard homotopy
      when 1 => return Job1;   -- initialize dobldobl homotopy
      when 2 => return Job2;   -- initialize quaddobl homotopy
      when 3 => return Job3;   -- initialize standard solution
      when 4 => return Job4;   -- initialize dobldobl solution
      when 5 => return Job5;   -- initialize quaddobl solution
      when 6 => return Job6;   -- next standard solution
      when 7 => return Job7;   -- next dobldobl solution
      when 8 => return Job8;   -- next quaddobl solution
      when 9 => Standard_Path_Tracker.Clear; return 0;
      when 10 => DoblDobl_Path_Tracker.Clear; return 0;
      when 11 => QuadDobl_Path_Tracker.Clear; return 0;
      when 12 => return Job12; -- initialize multiprecision tracker
      when 13 => return Job13; -- initialize multiprecision solution
      when 14 => return Job14; -- next multiprecision solution
      when 15 => Multprec_Path_Tracker.Clear; return 0;
      when 16 => return Job16; -- initialize variable precision homotopy 
      when 17 => return Job17; -- initialize variable precision solution
      when 18 => return Job18; -- next variable precision solution
      when 19 => Varbprec_Path_Tracker.Clear; return 0;
      when 20 => return Job20; -- current variable precision solution
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_nxtsol;
