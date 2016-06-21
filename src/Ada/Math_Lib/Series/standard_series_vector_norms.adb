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
  
end Standard_Series_Vector_Norms;
