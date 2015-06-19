with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;            use Multprec_Complex_Numbers;
with Multprec_Speelpenning_Products;
with Multprec_Monomial_Evaluations;

package body Multprec_Gradient_Evaluations is

  function Reverse_Speel
             ( b : Standard_Natural_VecVecs.VecVec;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(b'range);
    z : Multprec_Complex_Vectors.Vector(0..x'last);

  begin
    for i in b'range loop
      z := Multprec_Speelpenning_Products.Reverse_Speel(b(i).all,x);
      res(i) := new Multprec_Complex_Vectors.Vector'(z);
    end loop;
    return res;
  end Reverse_Speel;

  procedure Reverse_Speel
             ( b : in Standard_Natural_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               s : in out Multprec_Complex_VecVecs.VecVec ) is

    z : Multprec_Complex_Vectors.Vector(0..x'last);
    sind : Multprec_Complex_Vectors.Link_to_Vector;

  begin
    for i in b'range loop
      z := Multprec_Speelpenning_Products.Reverse_Speel(b(i).all,x);
      sind := s(i);
      for j in z'range loop
        Copy(z(j),sind(j));
      end loop;
    end loop;
  end Reverse_Speel;

  function Gradient_Monomials
             ( b : Standard_Natural_VecVecs.VecVec;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_VecVecs.VecVec is

    y : Multprec_Complex_Vectors.Vector(b'range);
    s : Multprec_Complex_VecVecs.VecVec(b'range);
    zero : constant Complex_Number := Create(integer(0));

  begin
    for i in y'range loop
      y(i) := Create(integer(1));
    end loop;
    s := Reverse_Speel(b,x);
    for i in s'range loop
      Mul(s(i)(0),y(i));
      for j in 1..x'last loop
        if not Equal(s(i)(j),zero)
         then Multprec_Complex_Numbers.Mul(s(i)(j),y(i));
        end if;
      end loop;
    end loop;
    return s;
  end Gradient_Monomials;

  procedure Gradient_Monomials
             ( b : in Standard_Natural_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               s : in out Multprec_Complex_VecVecs.VecVec ) is

    y : Multprec_Complex_Vectors.Vector(b'range);
    zero : constant Complex_Number := Create(integer(0));
    sind : Multprec_Complex_Vectors.Link_to_Vector;

  begin
    for i in y'range loop
      y(i) := Create(integer(1));
    end loop;
    Reverse_Speel(b,x,s);
    for i in s'range loop
      sind := s(i);
      Mul(sind(0),y(i));
      for j in 1..x'last loop
        if not Equal(sind(j),zero)
         then Multprec_Complex_Numbers.Mul(sind(j),y(i));
        end if;
      end loop;
    end loop;
  end Gradient_Monomials;

  function Gradient_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_VecVecs.VecVec is

    y : Multprec_Complex_Vectors.Vector(f'range);
    s : Multprec_Complex_VecVecs.VecVec(b'range);
    zero : constant Complex_Number := Create(integer(0));
    m : natural32;
    mpm : Multprec_Floating_Numbers.Floating_Number;

  begin
    y := Multprec_Monomial_Evaluations.Eval_with_Power_Table(f,x);
    s := Reverse_Speel(b,x);
    for i in s'range loop
      Mul(s(i)(0),y(i));
      for j in 1..x'last loop
        if not Equal(s(i)(j),zero) then
          m := f(i)(j) + 1;
          if m > 1 then
            mpm := Multprec_Floating_Numbers.Create(m);
            Multprec_Complex_Numbers.Mul(s(i)(j),mpm);
            Multprec_Floating_Numbers.Clear(mpm);
          end if;
          Multprec_Complex_Numbers.Mul(s(i)(j),y(i));
        end if;
      end loop;
    end loop;
    return s;
  end Gradient_Monomials;

  procedure Gradient_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               s : in out Multprec_Complex_VecVecs.VecVec ) is

    y : Multprec_Complex_Vectors.Vector(f'range);
    zero : constant Complex_Number := Create(integer(0));
    m : natural32;
    mpm : Multprec_Floating_Numbers.Floating_Number;
    find : Standard_Natural_Vectors.Link_to_Vector;
    sind : Multprec_Complex_Vectors.Link_to_Vector;

  begin
    y := Multprec_Monomial_Evaluations.Eval_with_Power_Table(f,x);
    Reverse_Speel(b,x,s);
    for i in s'range loop
      sind := s(i);
      Mul(sind(0),y(i));
      find := f(i);
      for j in 1..x'last loop
        if not Equal(sind(j),zero) then
          m := find(j) + 1;
          if m > 1 then
            mpm := Multprec_Floating_Numbers.create(m);
            Multprec_Complex_Numbers.Mul(sind(j),mpm);
            Multprec_Floating_Numbers.Clear(mpm);
          end if;
          Multprec_Complex_Numbers.Mul(sind(j),y(i));
        end if;
      end loop;
    end loop;
  end Gradient_Monomials;

  function Gradient_Sum_of_Monomials
             ( b : Standard_Natural_VecVecs.VecVec;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n)
        := (0..n => Create(integer(0)));
    y : Multprec_Complex_VecVecs.VecVec(b'range);

  begin
    y := Gradient_Monomials(b,x);
    for i in y'range loop
      for j in res'range loop
        Multprec_Complex_Numbers.Add(res(j),y(i)(j));
      end loop;
    end loop;
    return res;
  end Gradient_Sum_of_Monomials;

  procedure Gradient_Sum_of_Monomials
             ( b : in Standard_Natural_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               y : in out Multprec_Complex_VecVecs.VecVec;
               r : out Multprec_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Multprec_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(b,x,y);
    r := (0..n => Create(integer(0)));
    for i in y'range loop
      yind := y(i);
      for j in r'range loop
        Multprec_Complex_Numbers.Add(r(j),yind(j));
      end loop;
    end loop;
  end Gradient_Sum_of_Monomials;

  function Gradient_Sum_of_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n)
        := (0..n => Create(integer(0)));
    y : Multprec_Complex_VecVecs.VecVec(b'range);

  begin
    y := Gradient_Monomials(f,b,x);
    for i in y'range loop
      for j in res'range loop
        Multprec_Complex_Numbers.Add(res(j),y(i)(j));
      end loop;
    end loop;
    return res;
  end Gradient_Sum_of_Monomials;

  procedure Gradient_Sum_of_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               y : in out Multprec_Complex_VecVecs.VecVec;
               r : out Multprec_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Multprec_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(f,b,x,y);
    r := (0..n => Create(integer(0)));
    for i in y'range loop
      yind := y(i);
      for j in r'range loop
        Multprec_Complex_Numbers.Add(r(j),yind(j));
      end loop;
    end loop;
  end Gradient_Sum_of_Monomials;

  function Gradient_of_Polynomial
             ( b : Standard_Natural_VecVecs.VecVec;
               c,x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n)
        := (0..n => Create(integer(0)));
    y : Multprec_Complex_VecVecs.VecVec(b'range);
    prd : Complex_Number;

  begin
    y := Gradient_Monomials(b,x);
    for i in y'range loop
      for j in res'range loop
        prd := c(i)*y(i)(j);
        Multprec_Complex_Numbers.Add(res(j),prd);
        Multprec_Complex_Numbers.Clear(prd);
      end loop;
    end loop;
    return res;
  end Gradient_of_Polynomial;

  procedure Gradient_of_Polynomial
             ( b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Multprec_Complex_Vectors.Link_to_Vector;
    prd : Complex_Number;

  begin
    Gradient_Monomials(b,x,wrk);
    ydx := (0..n => Create(integer(0)));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        prd := c(i)*yind(j);
        Multprec_Complex_Numbers.Add(ydx(j),prd);
        Multprec_Complex_Numbers.Clear(prd);
      end loop;
    end loop;
  end Gradient_of_Polynomial;

  function Gradient_of_Polynomial
             ( f,b : Standard_Natural_VecVecs.VecVec;
               c,x : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Multprec_Complex_Vectors.Vector(0..n)
        := (0..n => Create(integer(0)));
    y : Multprec_Complex_VecVecs.VecVec(b'range);
    prd : Complex_Number;

  begin
    y := Gradient_Monomials(f,b,x);
    for i in y'range loop
      for j in res'range loop
        prd := c(i)*y(i)(j);
        Multprec_Complex_Numbers.Add(res(j),prd);
        Multprec_Complex_Numbers.Clear(prd);
      end loop;
    end loop;
    return res;
  end Gradient_of_Polynomial;

  procedure Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Multprec_Complex_Vectors.Link_to_Vector;
    prd : Complex_Number;

  begin
    Gradient_Monomials(f,b,x,wrk);
    ydx := (0..n => Create(integer(0)));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        prd := c(i)*yind(j);
        Multprec_Complex_Numbers.Add(ydx(j),prd);
        Multprec_Complex_Numbers.Clear(prd);
      end loop;
    end loop;
  end Gradient_of_Polynomial;

  procedure Conditioned_Gradient_of_Polynomial
             ( b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector;
               numcnd : out Multprec_Floating_Vectors.Vector ) is

    use Multprec_Floating_Numbers;

    n : constant integer32 := x'last;
    yind : Multprec_Complex_Vectors.Link_to_Vector;
    ciyj : Complex_Number;
    absc : Floating_Number;

  begin
    Gradient_Monomials(b,x,wrk);
    ydx := (0..n => Create(integer(0)));
    numcnd := (0..n => create(integer(0)));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ciyj := c(i)*yind(j);
        Multprec_Complex_Numbers.Add(ydx(j),ciyj);
        absc := Absval(ciyj); 
        Multprec_Floating_Numbers.Add(numcnd(j),absc);
        Multprec_Floating_Numbers.Clear(absc);
        Multprec_Complex_Numbers.Clear(ciyj);
      end loop;
    end loop;
  end Conditioned_Gradient_of_Polynomial;

  procedure Conditioned_Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector;
               numcnd : out Multprec_Floating_Vectors.Vector ) is

    use Multprec_Floating_Numbers;

    n : constant integer32 := x'last;
    yind : Multprec_Complex_Vectors.Link_to_Vector;
    ciyj : Complex_Number;
    absc : Floating_Number;

  begin
    Gradient_Monomials(f,b,x,wrk);
    ydx := (0..n => Create(integer(0)));
    numcnd := (0..n => create(integer(0)));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ciyj := c(i)*yind(j);
        Multprec_Complex_Numbers.Add(ydx(j),ciyj);
        absc := Absval(ciyj); 
        Multprec_Floating_Numbers.Add(numcnd(j),absc);
        Multprec_Floating_Numbers.Clear(absc);
        Multprec_Complex_Numbers.Clear(ciyj);
      end loop;
    end loop;
  end Conditioned_Gradient_of_Polynomial;

end Multprec_Gradient_Evaluations;
