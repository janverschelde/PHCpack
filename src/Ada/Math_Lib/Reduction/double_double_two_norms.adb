with DoblDobl_Mathematical_Functions;    use DoblDobl_Mathematical_Functions;

package body Double_Double_Two_Norms is

  function Norm2 ( v : Vector ) return double_double is

    res : double_double := create(0.0);

  begin
    for i in v'range loop
      res := res + v(i)*v(i);
    end loop;
    res := SQRT(res);
    return res;
  end Norm2;

  procedure Normalize ( v : in out Vector ) is

    nrm : constant double_double := Norm2(v);

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

end Double_Double_Two_Norms;
