with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;

package body Standard_Recentering_Coordinates is

  function Linear_Offset_Shift
             ( a,b : Vector; t : Complex_Number ) return Vector is

    res : Vector(b'range);
    one : constant Complex_Number := Create(1.0);
    s : constant Complex_Number := one - t;

  begin
    for i in b'range loop
      res(i) := s*a(i) + t*b(i);
    end loop;
    return res;
  end Linear_Offset_Shift;

  function Complement_of_Projection
             ( p : Matrix; v : Vector ) return Vector is

    res : Vector(v'range) := v;
    w : Vector(v'range);
    c : Complex_Number;

  begin
    for j in 1..p'last(2) loop
      for i in w'range loop
        w(i) := p(i,j);
      end loop;
      c := Conjugated_Inner_Product(w,res);
      for i in v'range loop
        res(i) := res(i) - c*w(i);
      end loop;
    end loop;
    return res;
  end Complement_of_Projection;

  function Distance ( p : Matrix; x : Vector ) return double_float is

    bx,c : Vector(p'range(1));

  begin
    for i in bx'range loop
      bx(i) := p(i,0) - x(i);
    end loop;
    c := Complement_of_Projection(p,bx);
    return Norm2(c);
  end Distance;

end Standard_Recentering_Coordinates;
