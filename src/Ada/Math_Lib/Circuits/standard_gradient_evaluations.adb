with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Speelpenning_Products;
with Standard_Monomial_Evaluations;

package body Standard_Gradient_Evaluations is

  function Reverse_Speel
             ( b : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(b'range);
    z : Standard_Complex_Vectors.Vector(0..x'last);

  begin
    for i in b'range loop
      z := Standard_Speelpenning_Products.Reverse_Speel(b(i).all,x);
      res(i) := new Standard_Complex_Vectors.Vector'(z);
    end loop;
    return res;
  end Reverse_Speel;

  procedure Reverse_Speel
             ( b : in Standard_Natural_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               s : in out Standard_Complex_VecVecs.VecVec ) is

    z : Standard_Complex_Vectors.Vector(0..x'last);
    sind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in b'range loop
      z := Standard_Speelpenning_Products.Reverse_Speel(b(i).all,x);
      sind := s(i);
      for j in z'range loop
        sind(j) := z(j);
      end loop;
    end loop;
  end Reverse_Speel;

  function Gradient_Monomials
             ( b : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    y : Standard_Complex_Vectors.Vector(b'range);
    s : Standard_Complex_VecVecs.VecVec(b'range);
    z : Complex_Number;
    zero : constant Complex_Number := Create(0.0);

  begin
    for i in y'range loop
      y(i) := Create(1.0);
    end loop;
    s := Reverse_Speel(b,x);
    for i in s'range loop
      s(i)(0) := y(i)*s(i)(0);
      for j in 1..x'last loop
        z := s(i)(j);
        if z /= zero
         then s(i)(j) := y(i)*z;
        end if;
      end loop;
    end loop;
    return s;
  end Gradient_Monomials;

  procedure Gradient_Monomials
             ( b : in Standard_Natural_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               s : in out Standard_Complex_VecVecs.VecVec ) is

    y : Standard_Complex_Vectors.Vector(b'range);
    z : Complex_Number;
    zero : constant Complex_Number := Create(0.0);
    sind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in y'range loop
      y(i) := Create(1.0);
    end loop;
    Reverse_Speel(b,x,s);
    for i in s'range loop
      sind := s(i);
      sind(0) := y(i)*sind(0);
      for j in 1..x'last loop
        z := sind(j);
        if z /= zero
         then sind(j) := y(i)*z;
        end if;
      end loop;
    end loop;
  end Gradient_Monomials;

  function Gradient_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    y : Standard_Complex_Vectors.Vector(f'range);
    s : Standard_Complex_VecVecs.VecVec(b'range);
    z : Complex_Number;
    zero : constant Complex_Number := Create(0.0);
    m : natural32;

  begin
    y := Standard_Monomial_Evaluations.Eval_with_Power_Table(f,x);
    s := Reverse_Speel(b,x);
    for i in s'range loop
      s(i)(0) := y(i)*s(i)(0);
      for j in 1..x'last loop
        z := s(i)(j);
        if z /= zero then
          m := f(i)(j) + 1;
          if m > 1 
           then s(i)(j) := double_float(m)*y(i)*z;
           else s(i)(j) := y(i)*z;
          end if;
        end if;
      end loop;
    end loop;
    return s;
  end Gradient_Monomials;

  procedure Gradient_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               s : in out Standard_Complex_VecVecs.VecVec ) is

    y : Standard_Complex_Vectors.Vector(f'range);
    z : Complex_Number;
    zero : constant Complex_Number := Create(0.0);
    m : natural32;
    find : Standard_Natural_Vectors.Link_to_Vector;
    sind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    y := Standard_Monomial_Evaluations.Eval_with_Power_Table(f,x);
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
           then sind(j) := double_float(m)*y(i)*z;
           else sind(j) := y(i)*z;
          end if;
        end if;
      end loop;
    end loop;
  end Gradient_Monomials;

  function Gradient_Sum_of_Monomials
             ( b : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Standard_Complex_Vectors.Vector(0..n) := (0..n => Create(0.0));
    y : Standard_Complex_VecVecs.VecVec(b'range);
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    y := Gradient_Monomials(b,x);
    for i in y'range loop
      yind := y(i);
      for j in res'range loop
        res(j) := res(j) + yind(j);
      end loop;
    end loop;
    return res;
  end Gradient_Sum_of_Monomials;

  procedure Gradient_Sum_of_Monomials
             ( b : in Standard_Natural_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               y : in out Standard_Complex_VecVecs.VecVec;
               r : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(b,x,y);
    r := (0..n => Create(0.0));
    for i in y'range loop
      yind := y(i);
      for j in r'range loop
        r(j) := r(j) + yind(j);
      end loop;
    end loop;
  end Gradient_Sum_of_Monomials;

  function Gradient_Sum_of_Monomials
             ( f,b : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Standard_Complex_Vectors.Vector(0..n) := (0..n => Create(0.0));
    y : Standard_Complex_VecVecs.VecVec(b'range);
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    y := Gradient_Monomials(f,b,x);
    for i in y'range loop
      yind := y(i);
      for j in res'range loop
        res(j) := res(j) + yind(j);
      end loop;
    end loop;
    return res;
  end Gradient_Sum_of_Monomials;

  procedure Gradient_Sum_of_Monomials
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               y : in out Standard_Complex_VecVecs.VecVec;
               r : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(f,b,x,y);
    r := (0..n => Create(0.0));
    for i in y'range loop
      yind := y(i);
      for j in r'range loop
        r(j) := r(j) + yind(j);
      end loop;
    end loop;
  end Gradient_Sum_of_Monomials;

  function Gradient_of_Polynomial
             ( b : Standard_Natural_VecVecs.VecVec;
               c,x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Standard_Complex_Vectors.Vector(0..n) := (0..n => Create(0.0));
    y : Standard_Complex_VecVecs.VecVec(b'range);
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    y := Gradient_Monomials(b,x);
    for i in y'range loop
      yind := y(i);
      for j in res'range loop
        res(j) := res(j) + c(i)*yind(j);
      end loop;
    end loop;
    return res;
  end Gradient_of_Polynomial;

  procedure Gradient_of_Polynomial
             ( b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.VecVec;
               ydx : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(b,x,wrk);
    ydx := (0..n => Create(0.0));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ydx(j) := ydx(j) + c(i)*yind(j);
      end loop;
    end loop;
  end Gradient_of_Polynomial;

  function Gradient_of_Polynomial
             ( f,b : Standard_Natural_VecVecs.VecVec;
               c,x : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    n : constant integer32 := x'last;
    res : Standard_Complex_Vectors.Vector(0..n) := (0..n => Create(0.0));
    y : Standard_Complex_VecVecs.VecVec(b'range);
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    y := Gradient_Monomials(f,b,x);
    for i in y'range loop
      yind := y(i);
      for j in res'range loop
        res(j) := res(j) + c(i)*yind(j);
      end loop;
    end loop;
    return res;
  end Gradient_of_Polynomial;

  procedure Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.VecVec;
               ydx : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Standard_Complex_Vectors.Link_to_Vector;

  begin
    Gradient_Monomials(f,b,x,wrk);
    ydx := (0..n => Create(0.0));
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ydx(j) := ydx(j) + c(i)*yind(j);
      end loop;
    end loop;
  end Gradient_of_Polynomial;

  procedure Conditioned_Gradient_of_Polynomial
             ( b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.VecVec;
               ydx : out Standard_Complex_Vectors.Vector;
               numcnd : out Standard_Floating_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Standard_Complex_Vectors.Link_to_Vector;
    ciyj : Complex_Number;

  begin
    Gradient_Monomials(b,x,wrk);
    ydx := (0..n => Create(0.0));
    numcnd := (0..n => 0.0);
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ciyj := c(i)*yind(j);
        ydx(j) := ydx(j) + ciyj;
        numcnd(j) := numcnd(j) + Absval(ciyj); 
      end loop;
    end loop;
  end Conditioned_Gradient_of_Polynomial;

  procedure Conditioned_Gradient_of_Polynomial
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.VecVec;
               ydx : out Standard_Complex_Vectors.Vector;
               numcnd : out Standard_Floating_Vectors.Vector ) is

    n : constant integer32 := x'last;
    yind : Standard_Complex_Vectors.Link_to_Vector;
    ciyj : Complex_Number;

  begin
    Gradient_Monomials(f,b,x,wrk);
    ydx := (0..n => Create(0.0));
    numcnd := (0..n => 0.0);
    for i in wrk'range loop
      yind := wrk(i);
      for j in ydx'range loop
        ciyj := c(i)*yind(j);
        ydx(j) := ydx(j) + ciyj;
        numcnd(j) := numcnd(j) + Absval(ciyj); 
      end loop;
    end loop;
  end Conditioned_Gradient_of_Polynomial;

end Standard_Gradient_Evaluations;
