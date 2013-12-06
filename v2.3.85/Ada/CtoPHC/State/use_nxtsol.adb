with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
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
with Multprec_Solutions_Container;
with PHCpack_Operations;
with Standard_Path_Tracker;
with DoblDobl_Path_Tracker;
with QuadDobl_Path_Tracker;
with Multprec_Path_Tracker;

function use_nxtsol ( job : integer32;
                      a : C_intarrs.Pointer;
                      b : C_intarrs.Pointer;
                      c : C_dblarrs.Pointer ) return integer32 is
 
  function Job0 return integer32 is -- initialize standard homotopy

    start,target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    Standard_Path_Tracker.Init(target,start);
    return 0;
  end Job0;
 
  function Job1 return integer32 is -- initialize dobldobl homotopy

    start,target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    DoblDobl_Path_Tracker.Init(target,start);
    return 0;
  end Job1;
 
  function Job2 return integer32 is -- initialize quaddobl homotopy

    start,target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
   -- new_line;
   -- put_line("The target system : "); put(target.all);
   -- put_line("The start system : "); put(start.all);
    QuadDobl_Path_Tracker.Init(target,start);
    return 0;
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
  end Job8;

  function Job9 return integer32 is -- clear standard path tracker
  begin
    Standard_Path_Tracker.Clear;
    return 0;
  end Job9;

  function Job10 return integer32 is -- clear standard path tracker
  begin
    DoblDobl_Path_Tracker.Clear;
    return 0;
  end Job10;

  function Job11 return integer32 is -- clear standard path tracker
  begin
    QuadDobl_Path_Tracker.Clear;
    return 0;
  end Job11;

  function Job12 return integer32 is -- initialize multiprecision tracker

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    deci : constant natural32 := natural32(v_a(v_a'first));
    start,target : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(start);
    PHCpack_Operations.Retrieve_Target_System(target);
    Multprec_Path_Tracker.Init(target,start,deci);
    return 0;
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
  end Job14;

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
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_nxtsol;
