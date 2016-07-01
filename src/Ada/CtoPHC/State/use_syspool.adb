with Characters_and_Numbers; use Characters_and_Numbers;

-- with GNAT.Float_Control;

with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with PHCpack_Operations;                use PHCpack_Operations;
with Standard_PolySys_Container;
with Standard_Systems_Pool;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_PolySys_Container;
with DoblDobl_Systems_Pool;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_PolySys_Container;
with QuadDobl_Systems_Pool;
with Solutions_Pool;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

function use_syspool ( job : integer32;
                       a : C_intarrs.Pointer;
		       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  function Job0 return integer32 is -- initialize pool with n = a[0]

    v : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant integer32 := integer32(v(v'first));

  begin
    Standard_Systems_Pool.Initialize(n);
    return 0;
  end Job0;

  function Job1 return integer32 is -- size of the standard pool

    n : constant natural32 := Standard_Systems_Pool.Size;

  begin
    Assign(integer32(n),a);
    return 0;
  end Job1;

  function Job9 return integer32 is -- size of the dobldobl pool

    n : constant natural32 := DoblDobl_Systems_Pool.Size;

  begin
    Assign(integer32(n),a);
    return 0;
  end Job9;

  function Job10 return integer32 is -- size of the quaddobl pool

    n : constant natural32 := QuadDobl_Systems_Pool.Size;

  begin
    Assign(integer32(n),a);
    return 0;
  end Job10;

  function Job2 return integer32 is -- read and create k-th system, k = a[0]

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : Link_to_Poly_Sys;
    
  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(p);
    Standard_Systems_Pool.Create(k,p.all);
    return 0;
  end Job2;

  function Job3 return integer32 is -- write k-th system, k = a[0]

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : constant Link_to_Poly_Sys := Standard_Systems_Pool.Retrieve(k);

  begin
    if p /= null then
      if PHCpack_Operations.Is_File_Defined
       then put(PHCpack_Operations.output_file,natural32(p'last),p.all);
       else put(standard_output,natural32(p'last),p.all);
      end if;
    end if;
    return 0;
  end Job3;

  function Job4 return integer32 is -- creates k-th system from container

    v : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant integer32 := integer32(v(v'first));
    p : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if p /= null
     then Standard_Systems_Pool.Create(k,p.all);
    end if;
    return 0;
  end Job4;

  function Job5 return integer32 is -- refines a solution using k-th system

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant integer32 := integer32(v_a(v_a'first));
   -- n : constant integer32 := integer32(v_a(v_a'first+1));
    f : constant Link_to_Eval_Poly_Sys := Standard_Systems_Pool.Evaluator(k);
    jf : constant Link_to_Eval_Jaco_Mat
       := Standard_Systems_Pool.Jacobian_Evaluator(k);
    sols : constant Solution_List := Solutions_Pool.Retrieve(k);
    len : constant natural32 := Solutions_Pool.Length(k);
    tmp : Solution_List;
   -- sol : Solution(n) := Convert_to_Solution(b,c);
    ls : Link_to_Solution; -- := Convert_to_Solution(b,c);
    epsxa : constant double_float := 1.0E-14;
    epsfa : constant double_float := 1.0E-14;
    max : constant natural32 := 3;
    cnt : natural32 := 0;
    numit : natural32 := 0;
    fail : boolean;
   -- x : Standard_Complex_Vectors.Vector(1..n) := (1..n => Create(1.0));
   -- y : Standard_Complex_Vectors.Vector(f'range);

  begin
   -- GNAT.Float_Control.Reset;
    put_line("Thread " & convert(k) & " starts refining " 
                       & convert(integer32(len)) & " solutions.");
    tmp := sols;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp); cnt := cnt + 1;
      put_line("Thread " & convert(k) & " refines solution "
                         & convert(integer32(cnt)));
      Silent_Newton(f.all,jf.all,ls.all,epsxa,epsfa,numit,max,fail);
     -- y := Eval(f.all,ls.v);
      tmp := Tail_Of(tmp);
    end loop;
   -- Assign_Solution(sol,b,c);
   -- Assign_Solution(ls,b,c);
    put_line(" done");
    return 0;
  exception
    when others => put_line("exception occurred in root refiner...");
                   return 305;
  end Job5;

  function Job6 return integer32 is -- copy to standard container

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant integer32 := integer32(v_a(v_a'first));
    sys : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
        := Standard_Systems_Pool.Retrieve(ind);

  begin
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(sys.all);
    return 0;
  exception
    when others => put_line("exception when copying from standard pool");
                   return 313;
  end Job6;

  function Job7 return integer32 is -- copy to dobldobl container

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant integer32 := integer32(v_a(v_a'first));
    sys : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
        := DoblDobl_Systems_Pool.Retrieve(ind);

  begin
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(sys.all);
    return 0;
  exception
    when others => put_line("exception when copying from dobldobl pool");
                   return 314;
  end Job7;

  function Job8 return integer32 is -- copy to quaddobl container

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    ind : constant integer32 := integer32(v_a(v_a'first));
    sys : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
        := QuadDobl_Systems_Pool.Retrieve(ind);

  begin
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(sys.all);
    return 0;
  exception
    when others => put_line("exception when copying from quaddobl pool");
                   return 315;
  end Job8;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when  0 => return Job0;  -- initialize pool with a[0]
      when  1 => return Job1;  -- returns size of standard pool in a[0]
      when  2 => return Job2;  -- read and create k-th system, k = a[0]
      when  3 => return Job3;  -- write k-th system, k = a[0]
      when  4 => return Job4;  -- creates k-th system from container
      when  5 => return Job5;  -- refine a root with k-th system
      when  6 => return Job6;  -- copy to standard systems container
      when  7 => return Job7;  -- copy to dobldobl systems container
      when  8 => return Job8;  -- copy to quaddobl systems container
      when  9 => return Job9;  -- returns size of dobldobl pool in a[0]
      when 10 => return Job10; -- returns size of quaddobl pool in a[0]
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_syspool;
