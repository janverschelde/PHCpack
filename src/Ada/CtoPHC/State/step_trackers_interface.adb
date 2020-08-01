with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with String_Splitters;                  use String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
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

package body Step_Trackers_Interface is
 
  function Step_Trackers_Standard_Homotopy
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    start,target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Standard_Homotopy ...");
    end if;
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    if fixed = 1 then
      Standard_Path_Tracker.Init(target,start,true);
    else
      declare
        v_c : constant C_Double_Array
            := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
        regamma : constant double_float := double_float(v_c(v_c'first));
        imgamma : constant double_float := double_float(v_c(v_c'first+1));
        gamma : Standard_Complex_Numbers.Complex_Number;
      begin
        if regamma = 0.0 and imgamma = 0.0 then
          Standard_Path_Tracker.Init(target,start,false);
        else
          gamma := Standard_Complex_Numbers.Create(regamma,imgamma);
          Standard_Path_Tracker.Init(target,start,gamma,2,0);
        end if;
      end;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Standard_Homotopy.");
      end if;
      return 500;
  end Step_Trackers_Standard_Homotopy;
 
  function Step_Trackers_DoblDobl_Homotopy
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    start,target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_DoblDobl_Homotopy ...");
    end if;
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    if fixed = 1 then
      DoblDobl_Path_Tracker.Init(target,start,true);
    else
      declare
        v_c : constant C_Double_Array
            := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
        regamma : constant double_float := double_float(v_c(v_c'first));
        imgamma : constant double_float := double_float(v_c(v_c'first+1));
        gamma : DoblDobl_Complex_Numbers.Complex_Number;
        dd_regamma,dd_imgamma : double_double;
      begin
        if regamma = 0.0 and imgamma = 0.0 then
          DoblDobl_Path_Tracker.Init(target,start,false);
        else
          dd_regamma := create(regamma);
          dd_imgamma := create(imgamma);
          gamma := DoblDobl_Complex_Numbers.Create(dd_regamma,dd_imgamma);
          DoblDobl_Path_Tracker.Init(target,start,gamma,2,0);
        end if;
      end;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_DoblDobl_Homotopy.");
      end if;
      return 501;
  end Step_Trackers_DoblDobl_Homotopy;
 
  function Step_Trackers_QuadDobl_Homotopy
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    start,target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_QuadDobl_Homotopy ...");
    end if;
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    if fixed = 1 then
      QuadDobl_Path_Tracker.Init(target,start,true);
    else
      declare
        v_c : constant C_Double_Array
            := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2));
        regamma : constant double_float := double_float(v_c(v_c'first));
        imgamma : constant double_float := double_float(v_c(v_c'first+1));
        gamma : QuadDobl_Complex_Numbers.Complex_Number;
        qd_regamma,qd_imgamma : quad_double;
      begin
        if regamma = 0.0 and imgamma = 0.0 then
          QuadDobl_Path_Tracker.Init(target,start,false);
        else
          qd_regamma := create(regamma);
          qd_imgamma := create(imgamma);
          gamma := QuadDobl_Complex_Numbers.Create(qd_regamma,qd_imgamma);
          QuadDobl_Path_Tracker.Init(target,start,gamma,2,0);
        end if;
      end;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_DoblDobl_Homotopy.");
      end if;
      return 502;
  end Step_Trackers_QuadDobl_Homotopy;

  function Step_Trackers_Set_Standard_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Set_Standard_Solution ...");
    end if;
   -- put("initializing solution "); put(k,1); put_line(" :");
    Standard_Solutions_Container.Retrieve(k,ls,fail);
    Standard_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Set_Standard_Solution.");
      end if;
      return 503;
  end Step_Trackers_Set_Standard_Solution;

  function Step_Trackers_Set_DoblDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Set_DoblDobl_Solution ...");
    end if;
   -- put("initializing solution "); put(k,1); put_line(" :");
    DoblDobl_Solutions_Container.Retrieve(k,ls,fail);
    DoblDobl_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Set_DoblDobl_Solution.");
      end if;
      return 504;
  end Step_Trackers_Set_DoblDobl_Solution;

  function Step_Trackers_Set_QuadDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Set_QuadDobl_Solution ...");
    end if;
   -- put("initializing solution "); put(k,1); put_line(" :");
    QuadDobl_Solutions_Container.Retrieve(k,ls,fail);
    QuadDobl_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Set_QuadDobl_Solution.");
      end if;
      return 505;
  end Step_Trackers_Set_QuadDobl_Solution;

  function Step_Trackers_Next_Standard_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Standard_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Next_Standard_Solution ...");
    end if;
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := Standard_Path_Tracker.get_next;
    Standard_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Next_Standard_Solution.");
      end if;
      return 506;
  end Step_Trackers_Next_Standard_Solution;

  function Step_Trackers_Next_DoblDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Next_DoblDobl_Solution ...");
    end if;
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := DoblDobl_Path_Tracker.get_next;
    DoblDobl_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Next_DoblDobl_Solution.");
      end if;
      return 507;
  end Step_Trackers_Next_DoblDobl_Solution;

  function Step_Trackers_Next_QuadDobl_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Next_QuadDobl_Solution ...");
    end if;
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := QuadDobl_Path_Tracker.get_next;
    QuadDobl_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Next_DoblDobl_Solution.");
      end if;
      return 508;
  end Step_Trackers_Next_QuadDobl_Solution;

  function Step_Trackers_Multprec_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    deci : constant natural32 := natural32(v_b(v_b'first));
    start,target : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Multprec_Homotopy ...");
    end if;
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    if fixed = 1
     then Multprec_Path_Tracker.Init(target,start,true,deci);
     else Multprec_Path_Tracker.Init(target,start,false,deci);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Multprec_Homotopy.");
      end if;
      return 512;
  end Step_Trackers_Multprec_Homotopy;

  function Step_Trackers_Set_Multprec_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Multprec_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Set_Multprec_Solution ...");
    end if;
   -- put("initializing solution "); put(k,1); put_line(" :");
    Multprec_Solutions_Container.Retrieve(k,ls,fail);
    Multprec_Path_Tracker.Init(ls);
   -- put(ls.all);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Set_Multprec_Solution.");
      end if;
      return 513;
  end Step_Trackers_Set_Multprec_Solution;

  function Step_Trackers_Next_Multprec_Solution
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    ls : Multprec_Complex_Solutions.Link_to_Solution;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Next_Multprec_Solution ...");
    end if;
   -- put("predictor-corrector step on solution "); put(k,1); new_line;
    ls := Multprec_Path_Tracker.get_next;
    Multprec_Solutions_Container.Replace(k,ls,fail);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Next_Multprec_Solution.");
      end if;
      return 514;
  end Step_Trackers_Next_Multprec_Solution;

  function Step_Trackers_Varbprec_Homotopy
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Varbprec_Homotopy ...");
    end if;
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
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Varbprec_Homotopy.");
      end if;
      return 516;
  end Step_Trackers_Varbprec_Homotopy;

  function Step_Trackers_Set_Varbprec_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

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
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Set_Varbprec_Solution ...");
    end if;
   -- put_line("initializing the variable precision path tracker ...");
   -- put_line("with the solution : "); put_line(lsl.all);
    Varbprec_Path_Tracker.Init(lsl,dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Set_Varbprec_Solution.");
      end if;
      return 517;
  end Step_Trackers_Set_Varbprec_Solution;

  function Step_Trackers_Next_Varbprec_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));
    want : constant natural32 := natural32(v_a(v_a'first));
    maxprc : constant natural32 := natural32(v_a(v_a'first+1));
    maxitr : constant natural32 := natural32(v_a(v_a'first+2));
    output : constant natural32 := natural32(v_a(v_a'first+3));
    otp : constant boolean := (output = 1);
    sol : Link_to_String;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Next_Varbprec_Solution ...");
    end if;
   -- put_line("inside use_nxtsol, calling get_current ...");
   -- sol := Varbprec_Path_Tracker.get_current;
   -- put_line("inside use_nxtsol, calling get_next ...");
    sol := Varbprec_Path_Tracker.get_next(want,maxprc,maxitr,otp);
   -- put_line("The solution returned : ");
   -- put_line(sol.all);
    Assign(integer32(sol'last),b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Next_Varbprec_Solution.");
      end if;
      return 518;
  end Step_Trackers_Next_Varbprec_Solution;

  function Step_Trackers_Get_Varbprec_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    sol : Link_to_String;

  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Get_Varbprec_Solution ...");
    end if;
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
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Step_Trackers_Get_Varbprec_Solution.");
      end if;
     return 520;
  end Step_Trackers_Get_Varbprec_Solution;

  function Step_Trackers_Standard_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Standard_Clear ...");
    end if;
    Standard_Path_Tracker.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Standard_Clear.");
      end if;
      return 509;
  end Step_Trackers_Standard_Clear;

  function Step_Trackers_DoblDobl_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_DoblDobl_Clear ...");
    end if;
    DoblDobl_Path_Tracker.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("DoblDobl_Clear.");
      end if;
      return 510;
  end Step_Trackers_DoblDobl_Clear;

  function Step_Trackers_QuadDobl_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_QuadDobl_Clear ...");
    end if;
    QuadDobl_Path_Tracker.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("QuadDobl_Clear.");
      end if;
      return 511;
  end Step_Trackers_QuadDobl_Clear;

  function Step_Trackers_Multprec_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Multprec_Clear ...");
    end if;
    Multprec_Path_Tracker.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Multprec_Clear.");
      end if;
      return 512;
  end Step_Trackers_Multprec_Clear;

  function Step_Trackers_Varbprec_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in step_trackers_interface.");
      put_line("Step_Trackers_Varbprec_Clear ...");
    end if;
    Varbprec_Path_Tracker.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in step_trackers_interface.");
        put_line("Varbprec_Clear.");
      end if;
      return 515;
  end Step_Trackers_Varbprec_Clear;

end Step_Trackers_Interface;
