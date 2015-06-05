with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Natural_Vectors;

package body Cyclic_Roots_System is

  function Number_of_Monomials ( n,i : integer32 ) return integer32 is
  begin
    if i = n
     then return 2;
     else return n;
    end if;
  end Number_of_Monomials;

  function Support_of_Cyclic
             ( n,i : integer32 )
             return Standard_Natural_VecVecs.VecVec is

    dim : constant integer32 := Number_of_Monomials(n,i);
    res : Standard_Natural_VecVecs.VecVec(1..dim);
    deg : Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);
    tmp : natural32;

  begin
    for k in 1..i loop
      deg(k) := 1;
    end loop;
    res(1) := new Standard_Natural_Vectors.Vector'(deg);
    if dim = 2 then
      deg := (1..n => 0);
      res(2) := new Standard_Natural_Vectors.Vector'(deg);
    else
      for k in 2..res'last loop
        tmp := deg(n);
        for j in reverse 1..n-1 loop
          deg(j+1) := deg(j);
        end loop;
        deg(1) := tmp;
        res(k) := new Standard_Natural_Vectors.Vector'(deg);
      end loop;
    end if;
    return res;
  end Support_of_Cyclic;

  function Supports_of_Cyclic
             ( n : integer32 )
             return Standard_Natural_VecVecs.Array_of_VecVecs is

    res : Standard_Natural_VecVecs.Array_of_VecVecs(1..n);

  begin
    for i in 1..n loop
      declare
        s : constant Standard_Natural_VecVecs.VecVec
          := Support_of_Cyclic(n,i);
      begin
        res(i) := new Standard_Natural_VecVecs.VecVec'(s);
      end;
    end loop;
    return res;
  end Supports_of_Cyclic;

  function Standard_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    t : Term;

  begin
    t.cf := Standard_Complex_Numbers.Create(1.0);
    t.dg := new Standard_Natural_Vectors.Vector'(s(s'first).all);
    Add(res,t);
    if s'last = s'first+1 then
      t.dg.all := s(s'last).all;
      Sub(res,t);
    else
      for i in s'first+1..s'last loop
        t.dg.all := s(i).all;
        Add(res,t);
      end loop;
    end if;
    Clear(t);
    return res;
  end Standard_Cyclic_Polynomial;

  function DoblDobl_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    one : constant double_double := create(1.0);
    t : Term;

  begin
    t.cf := DoblDobl_Complex_Numbers.Create(one);
    t.dg := new Standard_Natural_Vectors.Vector'(s(s'first).all);
    Add(res,t);
    if s'last = s'first+1 then
      t.dg.all := s(s'last).all;
      Sub(res,t);
    else
      for i in s'first+1..s'last loop
        t.dg.all := s(i).all;
        Add(res,t);
      end loop;
    end if;
    Clear(t);
    return res;
  end DoblDobl_Cyclic_Polynomial;

  function QuadDobl_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;
    one : constant quad_double := create(1.0);
    t : Term;

  begin
    t.cf := QuadDobl_Complex_Numbers.Create(one);
    t.dg := new Standard_Natural_Vectors.Vector'(s(s'first).all);
    Add(res,t);
    if s'last = s'first+1 then
      t.dg.all := s(s'last).all;
      Sub(res,t);
    else
      for i in s'first+1..s'last loop
        t.dg.all := s(i).all;
        Add(res,t);
      end loop;
    end if;
    Clear(t);
    return res;
  end QuadDobl_Cyclic_Polynomial;

  function Multprec_Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec )
             return Multprec_Complex_Polynomials.Poly is

    use Multprec_Complex_Polynomials;

    res : Poly := Null_Poly;
    one : Floating_Number := create(1.0);
    t : Term;

  begin
    t.cf := Multprec_Complex_Numbers.Create(one);
    t.dg := new Standard_Natural_Vectors.Vector'(s(s'first).all);
    Add(res,t);
    if s'last = s'first+1 then
      t.dg.all := s(s'last).all;
      Sub(res,t);
    else
      for i in s'first+1..s'last loop
        t.dg.all := s(i).all;
        Add(res,t);
      end loop;
    end if;
    Clear(t);
    Clear(one);
    return res;
  end Multprec_Cyclic_Polynomial;

  function Standard_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Standard_Cyclic_Polynomial(s(i).all);
    end loop;
    return res;
  end Standard_Cyclic_Polynomial_System;

  function DoblDobl_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := DoblDobl_Cyclic_Polynomial(s(i).all);
    end loop;
    return res;
  end DoblDobl_Cyclic_Polynomial_System;

  function QuadDobl_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := QuadDobl_Cyclic_Polynomial(s(i).all);
    end loop;
    return res;
  end QuadDobl_Cyclic_Polynomial_System;

  function Multprec_Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Multprec_Cyclic_Polynomial(s(i).all);
    end loop;
    return res;
  end Multprec_Cyclic_Polynomial_System;

end Cyclic_Roots_System;
