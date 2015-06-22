with text_io;                           use text_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Random_Vectors;
with Multprec_Complex_Poly_SysFun;      use Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Laur_SysFun;      use Multprec_Complex_Laur_SysFun;
with Lexicographical_Supports;

package body Multprec_Polynomial_Flatteners is

-- CONSTRUCTORS for SUPPORTS :

  function Number_of_Terms ( p : Poly_Sys ) return natural32 is

    res : natural32 := 0;

  begin
    for i in p'range loop
      res := res + Multprec_Complex_Polynomials.Number_of_Terms(p(i));
    end loop;
    return res;
  end Number_of_Terms;

  function Number_of_Terms ( p : Laur_Sys ) return natural32 is

    res : natural32 := 0;

  begin
    for i in p'range loop
      res := res + Multprec_Complex_Laurentials.Number_of_Terms(p(i));
    end loop;
    return res;
  end Number_of_Terms;

  procedure Update_Supports ( s,s_last : in out List;
                              p : in Multprec_Complex_Polynomials.Poly ) is

    use Multprec_Complex_Polynomials;

    procedure Scan ( t : in Term; continue : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      if not Is_In(s,v)
       then Append(s,s_last,v);
      end if;
      continue := true;
    end Scan;
    procedure Scan_Support is new Visiting_Iterator(Scan);

  begin
    Scan_Support(p);
  end Update_Supports;

  procedure Update_Supports ( s,s_last : in out List;
                              p : in Multprec_Complex_Laurentials.Poly ) is

    use Multprec_Complex_Laurentials;

    procedure Scan ( t : in Term; continue : out boolean ) is

      v : constant Standard_Integer_Vectors.Vector := t.dg.all;

    begin
      if not Is_In(s,v)
       then Append(s,s_last,v);
      end if;
      continue := true;
    end Scan;
    procedure Scan_Support is new Visiting_Iterator(Scan);

  begin
    Scan_Support(p);
  end Update_Supports;

  function Distinct_Supports ( p : Poly_Sys ) return List is

    res,res_last : List;

  begin
    for i in p'range loop
      Update_Supports(res,res_last,p(i));
    end loop;
    return res;
  end Distinct_Supports;

  function Distinct_Supports ( p : Laur_Sys ) return List is

    res,res_last : List;

  begin
    for i in p'range loop
      Update_Supports(res,res_last,p(i));
    end loop;
    return res;
  end Distinct_Supports;

-- CONSTRUCTORS for the DENSE case :

  procedure Update_Coefficient_Matrix
              ( A : in out Multprec_Complex_Matrices.Matrix;
                i : in integer32;
                s : in Standard_Integer_VecVecs.VecVec;
                p : in Multprec_Complex_Polynomials.Poly ) is

    use Multprec_Complex_Polynomials;

    procedure Scan ( t : in Term; continue : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);
      k : integer32;

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      k := Lexicographical_Supports.Index(s,v);
      if k > 0
       then Copy(t.cf,A(i,k));
      end if;
      continue := true;
    end Scan;
    procedure Scan_Terms is new Visiting_Iterator(Scan);

  begin
    for j in A'range(2) loop
      A(i,j) := Create(integer(0));
    end loop;
    Scan_Terms(p);
  end Update_Coefficient_Matrix;

  procedure Update_Coefficient_Matrix
              ( A : in out Multprec_Complex_Matrices.Matrix;
                i : in integer32;
                s : in Standard_Integer_VecVecs.VecVec;
                p : in Multprec_Complex_Laurentials.Poly ) is

    use Multprec_Complex_Laurentials;

    procedure Scan ( t : in Term; continue : out boolean ) is

      v : constant Standard_Integer_Vectors.Vector := t.dg.all;
      k : integer32;

    begin
      k := Lexicographical_Supports.Index(s,v);
      if k > 0
       then Copy(t.cf,A(i,k));
      end if;
      continue := true;
    end Scan;
    procedure Scan_Terms is new Visiting_Iterator(Scan);

  begin
    for j in A'range(2) loop
      A(i,j) := Create(integer(0));
    end loop;
    Scan_Terms(p);
  end Update_Coefficient_Matrix;

  function Coefficient_Matrix 
             ( p : Poly_Sys; s : Standard_Integer_VecVecs.VecVec ) 
             return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(p'range,s'range);

  begin
    for i in p'range loop
      Update_Coefficient_Matrix(res,i,s,p(i));
    end loop;
    return res;
  end Coefficient_Matrix;

  function Coefficient_Matrix 
             ( p : Laur_Sys; s : Standard_Integer_VecVecs.VecVec ) 
             return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(p'range,s'range);

  begin
    for i in p'range loop
      Update_Coefficient_Matrix(res,i,s,p(i));
    end loop;
    return res;
  end Coefficient_Matrix;

  procedure Flatten ( p : in Poly_Sys;
                      s : out Standard_Integer_VecVecs.Link_to_VecVec;
                      c : out Multprec_Complex_Matrices.Link_to_Matrix ) is

    ds : List := Distinct_Supports(p);
    sp : List := Lexicographical_Supports.Sort(ds);
    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(sp);
    A : constant Multprec_Complex_Matrices.Matrix(p'range,v'range)
      := Coefficient_Matrix(p,v);

  begin
    s := new Standard_Integer_VecVecs.VecVec'(v);
    c := new Multprec_Complex_Matrices.Matrix'(A);
    Lists_of_Integer_Vectors.Shallow_Clear(ds);
    Lists_of_Integer_Vectors.Shallow_Clear(sp);
  end Flatten;

  procedure Flatten ( p : in Laur_Sys;
                      s : out Standard_Integer_VecVecs.Link_to_VecVec;
                      c : out Multprec_Complex_Matrices.Link_to_Matrix ) is

    ds : List := Distinct_Supports(p);
    sp : List := Lexicographical_Supports.Sort(ds);
    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(sp);
    A : constant Multprec_Complex_Matrices.Matrix(p'range,v'range)
      := Coefficient_Matrix(p,v);

  begin
    s := new Standard_Integer_VecVecs.VecVec'(v);
    c := new Multprec_Complex_Matrices.Matrix'(A);
    Lists_of_Integer_Vectors.Shallow_Clear(ds);
    Lists_of_Integer_Vectors.Shallow_Clear(sp);
  end Flatten;

-- CONSTRUCTORS for the SPARSE case :

  procedure Coefficients_of_Support
              ( p : in Multprec_Complex_Polynomials.Poly;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out Multprec_Complex_Vectors.Link_to_Vector;
                k : out Standard_Natural_Vectors.Link_to_Vector ) is

    use Multprec_Complex_Polynomials;
    ntp : constant integer32 := integer32(Number_of_Terms(p));
    ind : integer32 := 0;

    procedure Scan ( t : in Term; continue : out boolean ) is

      v : Standard_Integer_Vectors.Vector(t.dg'range);
      vk : integer32;

    begin
      for i in v'range loop
        v(i) := integer32(t.dg(i));
      end loop;
      ind := ind + 1;
      c(ind) := t.cf;
      vk := Lexicographical_Supports.Index(s,v);
      if vk > 0
       then k(ind) := natural32(vk);
      end if;
      continue := true;
    end Scan;
    procedure Scan_Terms is new Visiting_Iterator(Scan);

  begin
    c := new Multprec_Complex_Vectors.Vector(1..ntp);
    k := new Standard_Natural_Vectors.Vector(1..ntp);
    Scan_Terms(p);
  end Coefficients_of_Support;

  procedure Coefficients_of_Support
              ( p : in Multprec_Complex_Laurentials.Poly;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out Multprec_Complex_Vectors.Link_to_Vector;
                k : out Standard_Natural_Vectors.Link_to_Vector ) is

    use Multprec_Complex_Laurentials;
    ntp : constant integer32 := integer32(Number_of_Terms(p));
    ind : integer32 := 0;

    procedure Scan ( t : in Term; continue : out boolean ) is

      v : constant Standard_Integer_Vectors.Vector := t.dg.all;
      vk : integer32;

    begin
      vk := Lexicographical_Supports.Index(s,v);
      ind := ind + 1;
      c(ind) := t.cf;
      if vk > 0
       then k(ind) := natural32(vk); 
      end if;
      continue := true;
    end Scan;
    procedure Scan_Terms is new Visiting_Iterator(Scan);

  begin
    c := new Multprec_Complex_Vectors.Vector(1..ntp);
    k := new Standard_Natural_Vectors.Vector(1..ntp);
    Scan_Terms(p);
  end Coefficients_of_Support;

  procedure Coefficients_of_Supports
              ( p : in Poly_Sys;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out Multprec_Complex_VecVecs.VecVec;
                k : out Standard_Natural_VecVecs.VecVec ) is
  begin
    for i in p'range loop
      Coefficients_of_Support(p(i),s,c(i),k(i));
    end loop;
  end Coefficients_of_Supports;

  procedure Coefficients_of_Supports
              ( p : in Laur_Sys;
                s : in Standard_Integer_VecVecs.VecVec;
                c : out Multprec_Complex_VecVecs.VecVec;
                k : out Standard_Natural_VecVecs.VecVec ) is
  begin
    for i in p'range loop
      Coefficients_of_Support(p(i),s,c(i),k(i));
    end loop;
  end Coefficients_of_Supports;

-- EVALUATORS :

  function Eval ( v : Standard_Integer_Vectors.Vector;
                  x : Multprec_Complex_Vectors.Vector )
                return Complex_Number is

    res : Complex_Number := Create(integer(1));

  begin
    for i in v'range loop
      if v(i) /= 0
       then res := res*(x(i)**integer(v(i)));
      end if;
    end loop;
    return res;
  end Eval;

  function Compressed_Eval
             ( v : Standard_Integer_Vectors.Vector;
               x : Multprec_Complex_Vectors.Vector )
             return Complex_Number is

    res : Complex_Number := Create(integer(1));
    ind : integer32 := 1;
    pow : integer32;
    first : boolean := true;

  begin
    while ind < v'last loop
      pow := v(ind+1);
      if first then
        if pow = 1
         then res := x(v(ind));
         else res := x(v(ind))**integer(pow);
        end if;
        first := false;
      else
        if pow = 1
         then res := res*x(v(ind));
         else res := res*(x(v(ind))**integer(pow));
        end if;
      end if;
      ind := ind + 2;
    end loop;
    return res;
  end Compressed_Eval;

  function Eval ( v : Standard_Integer_VecVecs.VecVec;
                  x : Multprec_Complex_Vectors.Vector )
                return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Eval(v(i).all,x);
    end loop;
    return res;
  end Eval;

  function Factored_Eval
              ( v : Standard_Integer_VecVecs.VecVec;
                x : Multprec_Complex_Vectors.Vector )
              return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);
    ptr : Standard_Integer_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      ptr := v(i);
      res(i) := Eval(ptr(x'range),x);
      if ptr(0) > 0
       then res(i) := res(i)*res(ptr(0));
      end if;
    end loop;
    return res;
  end Factored_Eval;

  function Factored_Compressed_Eval
              ( v : Standard_Integer_VecVecs.VecVec;
                x : Multprec_Complex_Vectors.Vector )
              return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);
    ptr : Standard_Integer_Vectors.Link_to_Vector;

  begin
    for i in v'range loop
      ptr := v(i);
      res(i) := Compressed_Eval(ptr.all,x);
      if ptr(0) > 0
       then res(i) := res(i)*res(ptr(0));
      end if;
    end loop;
    return res;
  end Factored_Compressed_Eval;

  function Eval ( A : in Multprec_Complex_Matrices.Matrix;
                  v : in Standard_Integer_VecVecs.VecVec;
                  x : in Multprec_Complex_Vectors.Vector ) 
                return Multprec_Complex_Vectors.Vector is

    y : constant Multprec_Complex_Vectors.Vector := Eval(v,x);
    use Multprec_Complex_Matrices;
    res : constant Multprec_Complex_Vectors.Vector := A*y;

  begin
    return res;
  end Eval;

  function Eval ( c : Multprec_Complex_Vectors.Vector;
                  v : Standard_Integer_VecVecs.VecVec;
                  k : Standard_Natural_Vectors.Vector;
                  x : in Multprec_Complex_Vectors.Vector )
                return Complex_Number is

    res : Complex_Number := Create(integer(0));
    one : constant Complex_Number := Create(integer(1));

  begin
    for i in c'range loop
      if c(i) = one
       then res := res + Eval(v(integer32(k(i))).all,x);
       else res := res + c(i)*Eval(v(integer32(k(i))).all,x);
      end if;
    end loop;
    return res;
  end Eval;

  function Eval ( c,vx : Multprec_Complex_Vectors.Vector;
                  k : Standard_Natural_Vectors.Vector )
                return Complex_Number is

    res : Complex_Number := Create(integer(0));
    one : constant Complex_Number := Create(integer(1));

  begin
    for i in c'range loop
      if c(i) = one 
       then res := res + vx(integer32(k(i)));
       else res := res + c(i)*vx(integer32(k(i)));
      end if;
    end loop;
    return res;
  end Eval;

  function Eval ( c : Multprec_Complex_VecVecs.VecVec;
                  v : Standard_Integer_VecVecs.VecVec;
                  k : Standard_Natural_VecVecs.VecVec;
                  x : Multprec_Complex_Vectors.Vector )
                return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(c'range);
    vx : constant Multprec_Complex_Vectors.Vector := Eval(v,x);

  begin
    for i in res'range loop
      res(i) := Eval(c(i).all,vx,k(i).all);
    end loop;
    return res;
  end Eval;

  function Factored_Eval
             ( c : Multprec_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(c'range);
    vx : constant Multprec_Complex_Vectors.Vector := Factored_Eval(v,x);

  begin
    for i in res'range loop
      res(i) := Eval(c(i).all,vx,k(i).all);
    end loop;
    return res;
  end Factored_Eval;

  function Factored_Compressed_Eval
             ( c : Multprec_Complex_VecVecs.VecVec;
               v : Standard_Integer_VecVecs.VecVec;
               k : Standard_Natural_VecVecs.VecVec;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(c'range);
    vx : constant Multprec_Complex_Vectors.Vector
       := Factored_Compressed_Eval(v,x);

  begin
    for i in res'range loop
      res(i) := Eval(c(i).all,vx,k(i).all);
    end loop;
    return res;
  end Factored_Compressed_Eval;

-- TEST PROCEDURES :

  procedure Spy ( A : in Multprec_Complex_Matrices.Matrix; 
                  v : Standard_Integer_VecVecs.VecVec ) is

    zero : constant Complex_Number := Create(integer(0));

  begin
    for k in v'range loop
      put(v(k).all);
      for i in A'range(1) loop
        put(" ");
        if A(i,k) = zero
         then put("0");
         else put("*");
        end if;
      end loop;
      new_line;
    end loop;
  end Spy;

  procedure Test_Eval ( p : in Poly_Sys; sz : in natural32;
                        A : in Multprec_Complex_Matrices.Matrix;
                        v : in Standard_Integer_VecVecs.VecVec ) is

    f : constant Standard_Integer_Vectors.Link_to_Vector := v(v'first);
    x : Multprec_Complex_Vectors.Vector(f'range);
    y,z : Multprec_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    loop
      x := Multprec_Random_Vectors.Random_Vector(1,f'last,sz);
      y := Eval(p,x);
      z := Eval(A,v,x);
      put_line("value at the polynomial system :"); put_line(y);
      put_line("value at the flattened representation :"); put_line(z);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Eval;

  procedure Test_Eval ( p : in Laur_Sys; sz : in natural32;
                        A : in Multprec_Complex_Matrices.Matrix;
                        v : in Standard_Integer_VecVecs.VecVec ) is

    f : constant Standard_Integer_Vectors.Link_to_Vector := v(v'first);
    x : Multprec_Complex_Vectors.Vector(f'range);
    y,z : Multprec_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    loop
      x := Multprec_Random_Vectors.Random_Vector(1,f'last,sz);
      y := Eval(p,x);
      z := Eval(A,v,x);
      put_line("value at the polynomial system :"); put_line(y);
      put_line("value at the flattened representation :"); put_line(z);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Eval;

  procedure Test_Eval ( p : in Poly_Sys; sz : in natural32;
                        c : in Multprec_Complex_VecVecs.VecVec;
                        v,fv,cfv : in Standard_Integer_VecVecs.VecVec;
                        k : in Standard_Natural_VecVecs.VecVec ) is

    f : constant Standard_Integer_Vectors.Link_to_Vector := v(v'first);
    x : Multprec_Complex_Vectors.Vector(f'range);
    y,z,w1,w2 : Multprec_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    loop
      x := Multprec_Random_Vectors.Random_Vector(1,f'last,sz);
      y := Eval(p,x);
      z := Eval(c,v,k,x);
      w1 := Factored_Eval(c,fv,k,x);
      w2 := Factored_Compressed_Eval(c,cfv,k,x);
      put_line("value at the polynomial system :"); put_line(y);
      put_line("value at the flattened representation :"); put_line(z);
      put_line("value with factored form of monomials :"); put_line(w1);
      put_line("value with compressed factored form :"); put_line(w2);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Eval;

  procedure Test_Eval ( p : in Laur_Sys; sz : in natural32;
                        c : in Multprec_Complex_VecVecs.VecVec;
                        v : in Standard_Integer_VecVecs.VecVec;
                        k : in Standard_Natural_VecVecs.VecVec ) is

    f : constant Standard_Integer_Vectors.Link_to_Vector := v(v'first);
    x : Multprec_Complex_Vectors.Vector(f'range);
    y,z : Multprec_Complex_Vectors.Vector(p'range);
    ans : character;

  begin
    loop
      x := Multprec_Random_Vectors.Random_Vector(1,f'last,sz);
      y := Eval(p,x);
      z := Eval(c,v,k,x);
      put_line("value at the polynomial system :"); put_line(y);
      put_line("value at the flattened representation :"); put_line(z);
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Eval;

end Multprec_Polynomial_Flatteners;
