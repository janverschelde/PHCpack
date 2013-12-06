with Interfaces.C;                      use Interfaces.C;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

package body Assignments_of_Solutions is

-- AUXILIARIES :

  function DoblDobl_Create
              ( x : C_Double_Array; i : integer32 ) return double_double is

  -- DESCRIPTION :
  --   Returns a double double from x(i) and x(i+1).
 
    res : constant double_double
        := create(double_float(x(Interfaces.C.size_t(i))),
                  double_float(x(Interfaces.C.size_t(i+1))));

  begin
    return res;
  end DoblDobl_Create;

  function QuadDobl_Create
              ( x : C_Double_Array; i : integer32 ) return quad_double is

  -- DESCRIPTION :
  --   Returns a quad double from x(i), x(i+1), x(i+2), x(i+3).
 
    res_hi : constant double_double := DoblDobl_Create(x,i);
    res_lo : constant double_double := DoblDobl_Create(x,i+2);
    res : constant quad_double := create(res_hi,res_lo);

  begin
    return res;
  end QuadDobl_Create;

  function DoblDobl_Complex_Create
              ( x : C_Double_Array; i : integer32 )
              return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns a complex double double with real part in x(i) and x(i+1)
  --   and imaginary part in x(i+2) and x(i+3).
 
    res_re : constant double_double := DoblDobl_Create(x,i);
    res_im : constant double_double := DoblDobl_Create(x,i+2);
    res : constant DoblDobl_Complex_Numbers.Complex_Number
        := DoblDobl_Complex_Numbers.Create(res_re,res_im);

  begin
    return res;
  end DoblDobl_Complex_Create;

  function QuadDobl_Complex_Create
              ( x : C_Double_Array; i : integer32 )
              return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns a complex quad double 
  --   with real part in x(i), x(i+1), x(i+2), x(i+3);
  --   and imaginary part in x(i+4), x(i+5), x(i+6), x(i+7).
 
    res_re : constant quad_double := QuadDobl_Create(x,i);
    res_im : constant quad_double := QuadDobl_Create(x,i+4);
    res : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Complex_Numbers.Create(res_re,res_im);

  begin
    return res;
  end QuadDobl_Complex_Create;

  procedure Assign_Double_Double
               ( x : in double_double; i : in integer32;
                 y : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns the high and low part of x to y(i) and y(i+1).

  begin
    y(i) := hi_part(x);
    y(i+1) := lo_part(x);
  end Assign_Double_Double;

  procedure Assign_DoblDobl_Complex
               ( x : in DoblDobl_Complex_Numbers.Complex_Number;
                 i : in integer32;
                 y : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns the high and low part 
  --   of the real part of x to y(i), y(i+1); 
  --   and the high and low part
  --   of the imaginary part of x to y(i+2), y(i+3).

    x_re : constant double_double
         := DoblDobl_Complex_Numbers.REAL_PART(x);
    x_im : constant double_double
         := DoblDobl_Complex_Numbers.IMAG_PART(x);

  begin
    Assign_Double_Double(x_re,i,y);
    Assign_Double_Double(x_im,i+2,y);
  end Assign_DoblDobl_Complex;

  procedure Assign_Quad_Double
               ( x : in quad_double; i : in integer32;
                 y : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns x to y(i), y(i+1), y(i+2), y(i+3).

    x_hi : constant double_double := hi_part(x);
    x_lo : constant double_double := lo_part(x);

  begin
    Assign_Double_Double(x_hi,i,y);
    Assign_Double_Double(x_lo,i+2,y);
  end Assign_Quad_Double;

  procedure Assign_QuadDobl_Complex
               ( x : in QuadDobl_Complex_Numbers.Complex_Number;
                 i : in integer32;
                 y : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Assigns the real part of x to y(i), y(i+1), y(i+2), y(i+3); 
  --   and the imaginary part of x to y(i+4), y(i+5), y(i+6), y(i+7).

    x_re : constant quad_double
         := QuadDobl_Complex_Numbers.REAL_PART(x);
    x_im : constant quad_double
         := QuadDobl_Complex_Numbers.IMAG_PART(x);

  begin
    Assign_Quad_Double(x_re,i,y);
    Assign_Quad_Double(x_im,i+4,y);
  end Assign_QuadDobl_Complex;

-- TARGET FUNCTIONS :

  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Numbers;

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    nv : constant integer32 := integer32(v(0));
    dim : constant Interfaces.C.size_t := Interfaces.C.size_t(2*nv+5);
    sol : C_Double_Array(0..dim-1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(dim));
    res : Standard_Complex_Solutions.Solution(nv);
    ind : Interfaces.C.size_t := 2;

  begin
    res.t := Create(double_float(sol(0)),double_float(sol(1)));
    res.m := integer32(v(1));
    for i in res.v'range loop
      res.v(i) := Create(double_float(sol(ind)),double_float(sol(ind+1)));
      ind := ind + 2;
    end loop;
    res.err := double_float(sol(Interfaces.C.size_t(2*nv+2)));
    res.rco := double_float(sol(Interfaces.C.size_t(2*nv+3)));
    res.res := double_float(sol(Interfaces.C.size_t(2*nv+4)));
    return res;
  end Convert_to_Solution;

  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Numbers;

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    nv : constant integer32 := integer32(v(0));
    dim : constant Interfaces.C.size_t := Interfaces.C.size_t(4*nv+10);
    sol : C_Double_Array(0..dim-1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(dim));
    res : DoblDobl_Complex_Solutions.Solution(nv);
    ind : integer32 := 4;

  begin
    res.t := DoblDobl_Complex_Create(sol,0);
    res.m := integer32(v(1));
    for i in res.v'range loop
      res.v(i) := DoblDobl_Complex_Create(sol,ind);
      ind := ind + 4;
    end loop;
    res.err := DoblDobl_Create(sol,ind);
    res.rco := DoblDobl_Create(sol,ind+2);
    res.res := DoblDobl_Create(sol,ind+4);
    return res;
  end Convert_to_Solution;

  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Numbers;

    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    nv : constant integer32 := integer32(v(0));
    dim : constant Interfaces.C.size_t := Interfaces.C.size_t(8*nv+20);
    sol : C_Double_Array(0..dim-1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(dim));
    res : QuadDobl_Complex_Solutions.Solution(nv);
    ind : integer32 := 8;

  begin
    res.t := QuadDobl_Complex_Create(sol,0);
    res.m := integer32(v(1));
    for i in res.v'range loop
      res.v(i) := QuadDobl_Complex_Create(sol,ind);
      ind := ind + 8;
    end loop;
    res.err := QuadDobl_Create(sol,ind);
    res.rco := QuadDobl_Create(sol,ind+4);
    res.res := QuadDobl_Create(sol,ind+8);
    return res;
  end Convert_to_Solution;

  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return Standard_Complex_Solutions.Link_to_Solution is

    use Standard_Complex_Numbers;

    res : Standard_Complex_Solutions.Link_to_Solution;
    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    nv : constant integer32 := integer32(v(0));
    dim : constant Interfaces.C.size_t := Interfaces.C.size_t(2*nv+5);
    sol : C_Double_Array(0..dim-1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(dim));
    s : Standard_Complex_Solutions.Solution(nv);
    ind : Interfaces.C.size_t := 2;

  begin
    s.t := Create(double_float(sol(0)),double_float(sol(1)));
    s.m := integer32(v(1));
    for i in s.v'range loop
      s.v(i) := Create(double_float(sol(ind)),double_float(sol(ind+1)));
      ind := ind + 2;
    end loop;
    s.err := double_float(sol(Interfaces.C.size_t(2*nv+2)));
    s.rco := double_float(sol(Interfaces.C.size_t(2*nv+3)));
    s.res := double_float(sol(Interfaces.C.size_t(2*nv+4)));
    res := new Standard_Complex_Solutions.Solution'(s);
    return res;
  end Convert_to_Solution;

  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return DoblDobl_Complex_Solutions.Link_to_Solution is

    use DoblDobl_Complex_Numbers;

    res : DoblDobl_Complex_Solutions.Link_to_Solution;
    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    nv : constant integer32 := integer32(v(0));
    dim : constant Interfaces.C.size_t := Interfaces.C.size_t(4*nv+10);
    sol : C_Double_Array(0..dim-1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(dim));
    s : DoblDobl_Complex_Solutions.Solution(nv);
    ind : integer32 := 4;

  begin
    s.t := DoblDobl_Complex_Create(sol,0);
    s.m := integer32(v(1));
    for i in s.v'range loop
      s.v(i) := DoblDobl_Complex_Create(sol,ind);
      ind := ind + 4;
    end loop;
    s.err := DoblDobl_Create(sol,ind);
    s.rco := DoblDobl_Create(sol,ind+2);
    s.res := DoblDobl_Create(sol,ind+4);
    res := new DoblDobl_Complex_Solutions.Solution'(s);
    return res;
  end Convert_to_Solution;

  function Convert_to_Solution
             ( b : C_intarrs.Pointer; c : C_dblarrs.Pointer )
             return QuadDobl_Complex_Solutions.Link_to_Solution is

    use QuadDobl_Complex_Numbers;

    res : QuadDobl_Complex_Solutions.Link_to_Solution;
    v : constant C_Integer_Array(0..1)
      := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    nv : constant integer32 := integer32(v(0));
    dim : constant Interfaces.C.size_t := Interfaces.C.size_t(8*nv+20);
    sol : C_Double_Array(0..dim-1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(dim));
    s : QuadDobl_Complex_Solutions.Solution(nv);
    ind : integer32 := 8;

  begin
    s.t := QuadDobl_Complex_Create(sol,0);
    s.m := integer32(v(1));
    for i in s.v'range loop
      s.v(i) := QuadDobl_Complex_Create(sol,ind);
      ind := ind + 8;
    end loop;
    s.err := QuadDobl_Create(sol,ind);
    s.rco := QuadDobl_Create(sol,ind+4);
    s.res := QuadDobl_Create(sol,ind+8);
    res := new QuadDobl_Complex_Solutions.Solution'(s);
    return res;
  end Convert_to_Solution;

  procedure Assign_Solution
             ( s : in Standard_Complex_Solutions.Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := 2*s.n + 5;
    sol : Standard_Floating_Vectors.Vector(1..dim);
    ind : integer32 := 3;

  begin
    sol(1) := REAL_PART(s.t);
    sol(2) := IMAG_PART(s.t);
    Assign(s.m,b);
    for i in s.v'range loop
      sol(ind) := REAL_PART(s.v(i));
      sol(ind+1) := IMAG_PART(s.v(i));
      ind := ind + 2;
    end loop;
    sol(dim-2) := s.err;
    sol(dim-1) := s.rco;
    sol(dim)   := s.res;
    Assign(sol,c);
  end Assign_Solution;

  procedure Assign_Solution
             ( s : in DoblDobl_Complex_Solutions.Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := 4*s.n + 10;
    sol : Standard_Floating_Vectors.Vector(1..dim);
    ind : integer32 := 5;

  begin
    Assign(s.m,b);
    Assign_DoblDobl_Complex(s.t,1,sol);
    for i in s.v'range loop
      Assign_DoblDobl_Complex(s.v(i),ind,sol);
      ind := ind + 4;
    end loop;
    Assign_Double_Double(s.err,ind,sol);
    Assign_Double_Double(s.rco,ind+2,sol);
    Assign_Double_Double(s.res,ind+4,sol);
    Assign(sol,c);
  end Assign_Solution;

  procedure Assign_Solution
             ( s : in QuadDobl_Complex_Solutions.Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := 8*s.n + 20;
    sol : Standard_Floating_Vectors.Vector(1..dim);
    ind : integer32 := 9;

  begin
    Assign(s.m,b);
    Assign_QuadDobl_Complex(s.t,1,sol);
    for i in s.v'range loop
      Assign_QuadDobl_Complex(s.v(i),ind,sol);
      ind := ind + 8;
    end loop;
    Assign_Quad_Double(s.err,ind,sol);
    Assign_Quad_Double(s.rco,ind+4,sol);
    Assign_Quad_Double(s.res,ind+8,sol);
    Assign(sol,c);
  end Assign_Solution;

  procedure Assign_Solution
             ( ls : in Standard_Complex_Solutions.Link_to_Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer ) is

    use Standard_Complex_Numbers;

    dim : constant integer32 := 2*ls.n + 5;
    sol : Standard_Floating_Vectors.Vector(1..dim);
    ind : integer32 := 3;

  begin
    sol(1) := REAL_PART(ls.t);
    sol(2) := IMAG_PART(ls.t);
    Assign(ls.m,b);
    for i in ls.v'range loop
      sol(ind) := REAL_PART(ls.v(i));
      sol(ind+1) := IMAG_PART(ls.v(i));
      ind := ind + 2;
    end loop;
    sol(dim-2) := ls.err;
    sol(dim-1) := ls.rco;
    sol(dim)   := ls.res;
    Assign(sol,c);
  end Assign_Solution;

  procedure Assign_Solution
             ( ls : in DoblDobl_Complex_Solutions.Link_to_Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer ) is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := 4*ls.n + 10;
    sol : Standard_Floating_Vectors.Vector(1..dim);
    ind : integer32 := 5;

  begin
    Assign(ls.m,b);
    Assign_DoblDobl_Complex(ls.t,1,sol);
    for i in ls.v'range loop
      Assign_DoblDobl_Complex(ls.v(i),ind,sol);
      ind := ind + 4;
    end loop;
    Assign_Double_Double(ls.err,ind,sol);
    Assign_Double_Double(ls.rco,ind+2,sol);
    Assign_Double_Double(ls.res,ind+4,sol);
    Assign(sol,c);
  end Assign_Solution;

  procedure Assign_Solution
             ( ls : in QuadDobl_Complex_Solutions.Link_to_Solution;
               b : in C_intarrs.Pointer; c : in C_dblarrs.Pointer ) is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := 8*ls.n + 20;
    sol : Standard_Floating_Vectors.Vector(1..dim);
    ind : integer32 := 9;

  begin
    Assign(ls.m,b);
    Assign_QuadDobl_Complex(ls.t,1,sol);
    for i in ls.v'range loop
      Assign_QuadDobl_Complex(ls.v(i),ind,sol);
      ind := ind + 8;
    end loop;
    Assign_Quad_Double(ls.err,ind,sol);
    Assign_Quad_Double(ls.rco,ind+4,sol);
    Assign_Quad_Double(ls.res,ind+8,sol);
    Assign(sol,c);
  end Assign_Solution;

end Assignments_of_Solutions;
