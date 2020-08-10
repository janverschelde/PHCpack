with Characters_and_Numbers; use Characters_and_Numbers;

-- with GNAT.Float_Control;

with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Standard_Systems_Pool;
with Solutions_Pool;

with Standard_SysPool_Interface;
with DoblDobl_SysPool_Interface;
with QuadDobl_SysPool_Interface;

function use_syspool ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

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

  function Handle_Jobs return integer32 is

    use Standard_SysPool_Interface;
    use DoblDobl_SysPool_Interface;
    use QuadDobl_SysPool_Interface;

  begin
    case job is
      when  0 => return Standard_SysPool_Initialize(a,vrblvl);
      when  1 => return Standard_SysPool_Size(a,vrblvl);
      when  2 => return Standard_SysPool_Read(a,vrblvl);
      when  3 => return Standard_SysPool_Write(a,vrblvl);
      when  4 => return Standard_SysPool_from_Container(a,vrblvl);
      when  5 => return Job5;  -- refine a root with k-th system
      when  6 => return Standard_SysPool_into_Container(a,vrblvl);
      when  7 => return DoblDobl_SysPool_into_Container(a,vrblvl);
      when  8 => return QuadDobl_SysPool_into_Container(a,vrblvl);
      when  9 => return DoblDobl_SysPool_Size(a,vrblvl);
      when 10 => return QuadDobl_SysPool_Size(a,vrblvl);
      when 11 => return DoblDobl_SysPool_Initialize(a,vrblvl);
      when 12 => return QuadDobl_SysPool_Initialize(a,vrblvl);
      when 13 => return Standard_SysPool_Clear(vrblvl);
      when 14 => return DoblDobl_SysPool_Clear(vrblvl);
      when 15 => return QuadDobl_SysPool_Clear(vrblvl);
      when 16 => return DoblDobl_SysPool_from_Container(a,vrblvl);
      when 17 => return QuadDobl_SysPool_from_Container(a,vrblvl);
      when others => put_line("invalid operation"); return 1;
    end case;
  end Handle_Jobs;

begin
  return Handle_Jobs;
end use_syspool;
