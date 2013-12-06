with text_io,integer_io;                use text_io,integer_io;
with Interfaces.C;
with Standard_Complex_Vectors;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Coefficient_Support_Poly_Systems;  use Coefficient_Support_Poly_Systems;
with PHCpack;
with C_Double_Arrays;                   use C_Double_Arrays;
with C_to_Ada_Arrays;                   use C_to_Ada_Arrays;

procedure phc_solver ( n,m : in integer; mc : in C_intarrs.Pointer;
                       ns : in integer; s : in C_intarrs.Pointer;
                       nc : in integer; c : in C_dblarrs.Pointer ) is

  mva : C_Integer_Array(0..Interfaces.C.size_T(m-1))
      := C_intarrs.Value(mc,Interfaces.C.ptrdiff_T(m));
  sva : C_Integer_Array(0..Interfaces.C.size_T(ns-1))
      := C_intarrs.Value(s,Interfaces.C.ptrdiff_T(ns));
  cva : C_Double_Array(0..Interfaces.C.size_T(nc-1)) 
      := C_dblarrs.Value(c,Interfaces.C.ptrdiff_T(nc));
  p : Poly_Sys(1..m) := Create(n,mva,cva,sva);
  rc : natural;
  q : Poly_Sys(p'range);
  sols : Solution_List;

  function phc_sols ( n,r,ns : integer; s : C_Double_Array ) return integer;
  pragma Import(C,phc_sols,"phc_sols");

  -- DESCRIPTION :
  --   The function phc_sols takes the solutions from Ada to C.

  function Concat_Solution_Vectors
             ( sols : Solution_List )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Concatenates all solutions vectors into one long vector.

    res : Standard_Complex_Vectors.Vector(1..n*Length_Of(sols));
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : natural := 0;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      for i in ls.v'range loop
        cnt := cnt + 1;
        res(cnt) := ls.v(i);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Concat_Solution_Vectors;

  procedure Output_Solutions ( sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Converts the solution list into an array of C doubles and
  --   passes the solutions to the C function phc_sols.

    sv : constant Standard_Complex_Vectors.Vector
       := Concat_Solution_Vectors(sols);
    dv : constant C_Double_Array := Convert(sv);
    result_of_call : integer;

  begin
    result_of_call := phc_sols(n,rc,dv'length,dv);
  end Output_Solutions;

begin
  put_line("The polynomial system defined by coefficients and support : ");
  put_line(p);
  PHCpack.Count_Roots(p,rc,q,sols);
  PHCpack.Artificial_Parameter_Continuation(p,q,sols);
  PHCpack.Refine_Roots(Standard_Output,p,sols);
  Output_Solutions(sols);
end phc_solver;
