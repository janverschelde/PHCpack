with text_io,integer_io;                use text_io,integer_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Homotopy;
with Increment_and_Fix_Continuation;    use Increment_and_Fix_Continuation;
with Standard_Root_Refiners;            use Standard_Root_Refiners;
with Interfaces.C;
with C_Integer_Arrays_io;               use C_Integer_Arrays_io;
with C_Double_Arrays_io;                use C_Double_Arrays_io;
with C_to_Ada_Arrays;                   use C_to_Ada_Arrays;
with Coefficient_Support_Poly_Systems;  use Coefficient_Support_Poly_Systems;
with Coefficient_Solution_Vectors;      use Coefficient_Solution_Vectors;

procedure Track_Paths ( n,m : in integer;
                        p_mc : in C_intarrs.Pointer;
                        p_ns : in integer; p_s : in C_intarrs.Pointer;
                        p_nc : in integer; p_c : in C_dblarrs.Pointer;
                        q_mc : in C_intarrs.Pointer;
                        q_ns : in integer; q_s : in C_intarrs.Pointer;
                        q_nc : in integer; q_c : in C_dblarrs.Pointer;
                        nbmul : in integer; mul : in C_intarrs.Pointer;
                        nbcfs : in integer; cfs : in C_dblarrs.Pointer ) is 

  p_mva : C_Integer_Array(0..Interfaces.C.size_T(m-1))
      := C_intarrs.Value(p_mc,Interfaces.C.ptrdiff_T(m));
  p_sva : C_Integer_Array(0..Interfaces.C.size_T(p_ns-1))
      := C_intarrs.Value(p_s,Interfaces.C.ptrdiff_T(p_ns));
  p_cva : C_Double_Array(0..Interfaces.C.size_T(p_nc-1)) 
      := C_dblarrs.Value(p_c,Interfaces.C.ptrdiff_T(p_nc));
  p : Poly_Sys(1..m) := Create(n,p_mva,p_cva,p_sva);
  q_mva : C_Integer_Array(0..Interfaces.C.size_T(m-1))
      := C_intarrs.Value(q_mc,Interfaces.C.ptrdiff_T(m));
  q_sva : C_Integer_Array(0..Interfaces.C.size_T(q_ns-1))
      := C_intarrs.Value(q_s,Interfaces.C.ptrdiff_T(q_ns));
  q_cva : C_Double_Array(0..Interfaces.C.size_T(q_nc-1)) 
      := C_dblarrs.Value(q_c,Interfaces.C.ptrdiff_T(q_nc));
  q : Poly_Sys(1..m) := Create(n,q_mva,q_cva,q_sva);
  mul_va : C_Integer_Array(0..Interfaces.C.size_T(nbmul-1))
         := C_intarrs.Value(mul,Interfaces.C.ptrdiff_T(nbmul));
  cfs_va : C_Double_Array(0..Interfaces.C.size_T(nbcfs-1))
         := C_dblarrs.Value(cfs,Interfaces.C.ptrdiff_T(nbcfs));
 -- s_mva : constant Standard_Integer_Vectors.Vector := Convert(mul_va);
 -- s_cva : constant Standard_Floating_Vectors.Vector := Convert(cfs_va);
 -- sols : Solution_List := Create(n,s_mva,s_cva);
  sols : Solution_List := Create(n,mul_va,cfs_va);

  function path_sols ( n,nbsols : integer; mul : C_Integer_Array;
                       nbcff : integer; cff : C_Double_Array )
                     return integer;
  pragma Import(C,path_sols,"path_sols");

  procedure Export_Results ( sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Calls the routine path_sols to export the calculated solutions
  --   to the C routines.

    n : constant natural := Head_Of(sols).n;
    m : constant C_Integer_Array := Multiplicities(sols);
    c : constant C_Double_Array := Coefficients(sols);
   -- m : constant Standard_Integer_Vectors.Vector := Multiplicities(sols);
   -- c : constant Standard_Floating_Vectors.Vector := Coefficients(sols);
   -- c_m : constant C_Integer_Array := Convert(m);
   -- c_c : constant C_Double_Array := Convert(c);
    result_of_call : integer;

  begin
   -- result_of_call := path_sols(n,m'length,c_m,c'length,c_c);
    result_of_call := path_sols(n,m'length,m,c'length,c);
    if result_of_call = 0
     then put_line("Call to path_sols ended normally.");
     else put_line("Call to path_sols ended abnormally.");
    end if;
  end Export_Results;

  procedure Refine_Roots ( file : in file_type;
                           p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the root refiner, for square polynomial systems only.

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural := 3;
    numit : natural := 0;

  begin
    Reporting_Root_Refiner
      (file,p,sols,epsxa,epsfa,tolsing,numit,max,false,false);
  end Refine_Roots;

  procedure Refine_Roots ( p : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Calls the root refiner, for square polynomial systems only.

    epsxa : constant double_float := 10.0**(-8);
    epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural := 3;
    numit : natural := 0;

  begin
    Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max,false);
  end Refine_Roots;

  procedure Continuation is

    a : Complex_Number := Random1;
    k : natural := 2;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                          Standard_Homotopy.Diff,Standard_Homotopy.Diff);

    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Standard_Homotopy.Eval,
                             Standard_Homotopy.Diff,Standard_Homotopy.Diff);

  begin
    Standard_Homotopy.Create(p,q,k,a);
    Sil_Cont(sols,false,Create(1.0));
    put_line("The solutions of the target system :");
    Refine_Roots(Standard_Output,p,sols);
    Export_Results(sols);
  end Continuation;

begin
 -- put_line("in the Ada routine...");
 -- put("n = "); put(n,1); new_line;
 -- put("m = "); put(m,1); new_line;
 -- put("mc = "); put(mva'length,mva); new_line;
 -- put("ns = "); put(ns,1); new_line;
 -- put("s = "); put(sva'length,sva); new_line;
 -- put("nc = "); put(nc,1); new_line;
 -- put_line("coefficients : "); put(cva'length,cva);
  put_line("The target system : "); put(p);
  put_line("The start system : "); put_line(q);
  put_line("The start solutions : ");
  put(Standard_Output,Length_Of(sols),Head_Of(sols).n,sols);
  Continuation;
 -- put_line("leaving the Ada routine...");
end Track_Paths;
