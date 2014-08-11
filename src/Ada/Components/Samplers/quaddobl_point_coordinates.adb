package body QuadDobl_Point_Coordinates is

  function Affine_Coordinates ( x : Vector ) return Vector is

    res : Vector(1..x'last);

  begin
    for i in res'range loop
      res(i) := x(i)/x(0);
    end loop;
    return res;
  end Affine_Coordinates;

  function Projective_Coordinates ( x : Vector ) return Vector is

    res : Vector(0..x'last);
    one : constant quad_double := create(1.0);

  begin
    res(0) := Create(one);
    res(x'range) := x;
    return res;
  end Projective_Coordinates;

  procedure Max_Norm ( x : in Vector;
                       ind : out integer32; nrm : out quad_double ) is

    tmp : quad_double;

  begin
    ind := x'first;
    nrm := AbsVal(x(ind));
    for i in x'first+1..x'last loop
      tmp := AbsVal(x(i));
      if tmp > nrm then
        nrm := tmp;
        ind := i;
      end if;
    end loop;
  end Max_Norm;

  procedure Scale ( x : in out Vector; ind : in integer32 ) is

    one : constant quad_double := create(1.0);

  begin
    for k in x'first..ind-1 loop
      x(k) := x(k)/x(ind);
    end loop;
    for k in ind+1..x'last loop
      x(k) := x(k)/x(ind);
    end loop;
    x(ind) := Create(one);
  end Scale;

  function Affine_Expand ( c : Complex_Number; b,v : Vector ) return Vector is

    res : constant Vector(b'range) := b + c*v;

  begin
    return res;
  end Affine_Expand;

  function Affine_Expand ( c,b : Vector; v : VecVec ) return Vector is

    res : Vector(b'range) := b;

  begin
    for i in v'range loop
      res := res + c(i)*v(i).all;
    end loop;
    return res;
  end Affine_Expand;

  function Affine_Expand ( c : Vector; p : Matrix ) return Vector is

    res : Vector(p'range(1));

  begin
    for i in res'range loop
      res(i) := p(i,0);
      for j in 1..p'last(2) loop
        res(i) := res(i) + c(j)*p(i,j);
      end loop;
    end loop;
    return res;
  end Affine_Expand;

  function Projective_Expand ( c : Vector; p : Matrix ) return Vector is

    res : Vector(p'range(1));

  begin
    for i in res'range loop
      res(i) := c(0)*p(i,0);
      for j in 1..p'last(2) loop
        res(i) := res(i) + c(j)*p(i,j);
      end loop;
    end loop;
    return res;
  end Projective_Expand;

  function Inner_Product ( u,v : Vector ) return Complex_Number is

    zero : constant quad_double := create(0.0);
    res : Complex_Number := Create(zero);

  begin
    for i in u'range loop
      res := res + u(i)*Conjugate(v(i));
    end loop;
    return res;
  end Inner_Product;

  function Project ( x,b : Vector; v : VecVec ) return Vector is

    y : constant Vector(b'range) := x-b;
    res : Vector(v'range);

  begin
    for i in res'range loop
      res(i) := Inner_Product(y,v(i).all);
    end loop;
    return res;
  end Project;

  function Project ( x : Vector; p : Matrix ) return Vector is

    res : Vector(1..p'last(2));
    y : Vector(p'range(1));
    zero : constant quad_double := create(0.0);

  begin
    for i in p'range(1) loop
      y(i) := x(i) - p(i,0);
    end loop;
    for i in res'range loop
      res(i) := Create(zero);
      for j in y'range loop
        res(i) := res(i) + y(j)*Conjugate(p(j,i));
      end loop;
    end loop;
    return res;
  end Project;

end QuadDobl_Point_Coordinates;
