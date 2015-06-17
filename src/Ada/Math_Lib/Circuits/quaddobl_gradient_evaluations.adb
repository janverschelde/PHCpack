with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with QuadDobl_Speelpenning_Products;
with QuadDobl_Monomial_Evaluations;

package body QuadDobl_Gradient_Evaluations is

  procedure Split_Common_Factor
              ( e : in Standard_Natural_Vectors.Vector;
                f,b : out Standard_Natural_Vectors.Vector ) is
  begin
    for i in e'range loop
      if e(i) = 0 then
        b(i) := 0; f(i) := 0;
      elsif e(i) = 1 then
        b(i) := 1; f(i) := 0;
      else -- e(i) > 1
        b(i) := 1; f(i) := e(i) - 1;
      end if;
    end loop;
  end Split_Common_Factor;

  procedure Split_Common_Factors
              ( e : in Standard_Natural_VecVecs.VecVec;
                f,b : out Standard_Natural_VecVecs.VecVec ) is
  begin
    for i in e'range loop
      declare
        ee : constant Standard_Natural_Vectors.Vector := e(i).all;
        ff,bb : Standard_Natural_Vectors.Vector(ee'range);
      begin
        Split_Common_Factor(ee,ff,bb);
        f(i) := new Standard_Natural_Vectors.Vector'(ff);
        b(i) := new Standard_Natural_Vectors.Vector'(bb);
      end;
    end loop;
  end Split_Common_Factors;

  function Reverse_Speel
             ( b : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(b'range);
    z : QuadDobl_Complex_Vectors.Vector(0..x'last);

  begin
    for i in b'range loop
      z := QuadDobl_Speelpenning_Products.Reverse_Speel(b(i).all,x);
      res(i) := new QuadDobl_Complex_Vectors.Vector'(z);
    end loop;
    return res;
  end Reverse_Speel;

  procedure Reverse_Speel
             ( b : in Standard_Natural_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               s : in out QuadDobl_Complex_VecVecs.VecVec ) is

    z : QuadDobl_Complex_Vectors.Vector(0..x'last);
    sind : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    for i in b'range loop
      z := QuadDobl_Speelpenning_Products.Reverse_Speel(b(i).all,x);
      sind := s(i);
      for j in z'range loop
        sind(j) := z(j);
      end loop;
    end loop;
  end Reverse_Speel;

  function Gradient_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec is

    y : QuadDobl_Complex_Vectors.Vector(f'range);
    s : QuadDobl_Complex_VecVecs.VecVec(b'range);
    z : Complex_Number;
    zero : constant Complex_Number := Create(integer(0));
    m : natural32;
    m_qd : Quad_Double_Numbers.quad_double;

  begin
    y := QuadDobl_Monomial_Evaluations.Eval_with_Power_Table(f,x);
    s := Reverse_Speel(b,x);
    for i in s'range loop
      s(i)(0) := y(i)*s(i)(0);
      for j in 1..x'last loop
        z := s(i)(j);
        if z /= zero then
          m := f(i)(j) + 1;
          if m > 1 
           then m_qd := Quad_Double_Numbers.Create(m);
                s(i)(j) := m_qd*y(i)*z;
           else s(i)(j) := y(i)*z;
          end if;
        end if;
      end loop;
    end loop;
    return s;
  end Gradient_Monomials;

  procedure Gradient_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               s : in out QuadDobl_Complex_VecVecs.VecVec ) is

    y : QuadDobl_Complex_Vectors.Vector(f'range);
    z : Complex_Number;
    zero : constant Complex_Number := Create(integer(0));
    m : natural32;
    find : Standard_Natural_Vectors.Link_to_Vector;
    sind : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    y := QuadDobl_Monomial_Evaluations.Eval_with_Power_Table(f,x);
    Reverse_Speel(b,x,s);
    for i in s'range loop
      sind := s(i);
      sind(0) := y(i)*sind(0);
      find := f(i);
      for j in 1..x'last loop
        z := sind(j);
        if z /= zero then
          m := find(j) + 1;
          if m > 1 
           then sind(j) := Quad_Double_Numbers.create(m)*y(i)*z;
           else sind(j) := y(i)*z;
          end if;
        end if;
      end loop;
    end loop;
  end Gradient_Monomials;

  function Gradient_Sum_of_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : QuadDobl_Complex_Vectors.Vector(0..n)
        := (0..n => Create(integer(0)));
    y : QuadDobl_Complex_VecVecs.VecVec(b'range);

  begin
    y := Gradient_Monomials(f,b,x);
    for i in y'range loop
      for j in res'range loop
        res(j) := res(j) + y(i)(j);
      end loop;
    end loop;
    return res;
  end Gradient_Sum_of_Monomials;

  procedure Gradient_Sum_of_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               y : in out QuadDobl_Complex_VecVecs.VecVec;
               r : out QuadDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(f,b,x,y);
    r := (0..n => Create(integer(0)));
    for i in y'range loop
      yind := y(i);
      for j in r'range loop
        r(j) := r(j) + yind(j);
      end loop;
    end loop;
  end Gradient_Sum_of_Monomials;

  function Gradient_of_Polynomial
             ( f,b : Standard_Natural_VecVecs.VecVec;
               c,x : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : QuadDobl_Complex_Vectors.Vector(0..n)
        := (0..n => Create(integer(0)));
    y : QuadDobl_Complex_VecVecs.VecVec(b'range);

  begin
    y := Gradient_Monomials(f,b,x);
    for i in y'range loop
      for j in res'range loop
        res(j) := res(j) + c(i)*y(i)(j);
      end loop;
    end loop;
    return res;
  end Gradient_of_Polynomial;

  procedure Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.VecVec;
               ydx : out QuadDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : QuadDobl_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(f,b,x,wrk);
    ydx := (0..n => Create(integer(0)));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ydx(j) := ydx(j) + c(i)*yind(j);
      end loop;
    end loop;
  end Gradient_of_Polynomial;

  procedure Conditioned_Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.VecVec;
               ydx : out QuadDobl_Complex_Vectors.Vector;
               numcnd : out Quad_Double_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : QuadDobl_Complex_Vectors.Link_to_Vector;
    ciyj : Complex_Number;

    use Quad_Double_Numbers;

  begin
    Gradient_Monomials(f,b,x,wrk);
    ydx := (0..n => Create(integer(0)));
    numcnd := (0..n => create(integer(0)));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ciyj := c(i)*yind(j);
        ydx(j) := ydx(j) + ciyj;
        numcnd(j) := numcnd(j) + Absval(ciyj); 
      end loop;
    end loop;
  end Conditioned_Gradient_of_Polynomial;

end QuadDobl_Gradient_Evaluations;
