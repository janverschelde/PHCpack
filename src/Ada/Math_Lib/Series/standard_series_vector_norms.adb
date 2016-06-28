with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Algebraic_Series;
with Standard_Dense_Series_Norms;

package body Standard_Series_Vector_Norms is

  function Conjugate ( v : Vector ) return Vector is

    res : Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_Dense_Series.Conjugate(v(i));
    end loop;
    return res;
  end Conjugate;

  function Inner_Product ( u,v : Vector ) return Series is

    w : constant Vector(u'range) := Conjugate(u);
    res : Series := Create(0);

  begin
    for i in w'range loop
      res := res + w(i)*v(i);
    end loop;
    return res;
  end Inner_Product;

  function Square_of_Norm ( v : Vector ) return Series is

    c : constant Vector(v'range) := Conjugate(v);
    res : constant Series := c*v;

  begin
    return res;
  end Square_of_Norm;

  function Norm ( v : Vector ) return Series is

    sn : constant Series := Square_of_Norm(v);
    res : constant Series := Standard_Algebraic_Series.sqrt(sn,0);

  begin
    return res;
  end Norm;

  procedure Normalize ( v : in out Vector ) is

    sn : constant Series := Norm(v);
    invsn : constant Series := Inverse(sn);

  begin
    Mul(v,invsn);
  end Normalize;

  function Normalize ( v : Vector ) return Vector is

    sn : constant Series := Norm(v);
    invsn : constant Series := Inverse(sn);
    res : constant Vector(v'range) := invsn*v;

  begin
    return res;
  end Normalize;

  function Max_Norm ( v : Vector ) return double_float is

    res : double_float := Standard_Dense_Series_Norms.Max_Norm(v(v'first));
    nrm : double_float;

  begin
    for i in v'first+1..v'last loop
      nrm := Standard_Dense_Series_Norms.Max_Norm(v(i));
      if nrm > res
       then res := nrm;
      end if;
    end loop;
    return res;
  end Max_Norm;
  
end Standard_Series_Vector_Norms;
