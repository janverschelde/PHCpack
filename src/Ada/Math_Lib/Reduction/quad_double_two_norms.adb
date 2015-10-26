with QuadDobl_Mathematical_Functions;    use QuadDobl_Mathematical_Functions;

package body Quad_Double_Two_Norms is

  function Norm2 ( v : Vector ) return quad_double is

    res : quad_double := create(0.0);

  begin
    for i in v'range loop
      res := res + v(i)*v(i);
    end loop;
    res := SQRT(res);
    return res;
  end Norm2;

  procedure Normalize ( v : in out Vector ) is

    nrm : constant quad_double := Norm2(v);

  begin
    for i in v'range loop
      v(i) := v(i)/nrm;
    end loop;
  end Normalize;

  procedure Normalize( v : in out VecVec ) is
  begin
    for i in v'range loop
      Normalize(v(i).all);
    end loop;
  end Normalize;

end Quad_Double_Two_Norms;
