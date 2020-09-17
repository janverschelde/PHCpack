with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with PentDobl_Complex_Algebraic_Series;
with PentDobl_Complex_Series_Norms;

package body PentDobl_CSeries_Vector_Norms is

  function Conjugate ( v : Vector ) return Vector is

    res : Vector(v'range);

  begin
    for i in v'range loop
      res(i) := PentDobl_Complex_Series_Norms.Conjugate(v(i));
    end loop;
    return res;
  end Conjugate;

  function Inner_Product ( u,v : Vector ) return Series is

    w : constant Vector(u'range) := Conjugate(u);
    deg : constant integer32 := u(u'first).cff'last;
    res : Series := Create(0,deg);

  begin
    for i in w'range loop
      res := res + w(i).all*v(i).all;
    end loop;
    return res;
  end Inner_Product;

  function Square_of_Norm ( v : Vector ) return Series is

    c : constant Vector(v'range) := Conjugate(v);
    res : constant Link_to_Series := c*v;

  begin
    return res.all;
  end Square_of_Norm;

  function Norm ( v : Vector ) return Series is

    sn : constant Series := Square_of_Norm(v);
    res : constant Series := PentDobl_Complex_Algebraic_Series.sqrt(sn,0);

  begin
    return res;
  end Norm;

  procedure Normalize ( v : in out Vector ) is

    sn : constant Series := Norm(v);
    invsn : constant Series := Inverse(sn);

  begin
    for i in v'range loop
      Mul(v(i).all,invsn);
    end loop;
  end Normalize;

  function Normalize ( v : Vector ) return Vector is

    sn : constant Series := Norm(v);
    invsn : constant Series := Inverse(sn);
    res : Vector(v'range);

  begin
    for i in v'range loop
      declare
        prd : constant Series := invsn*v(i).all;
      begin
        res(i) := new Series'(prd);
      end;
    end loop;
    return res;
  end Normalize;

  function Max_Norm ( v : Vector ) return penta_double is

    res : penta_double
        := PentDobl_Complex_Series_Norms.Max_Norm(v(v'first).all);
    nrm : penta_double;

  begin
    for i in v'first+1..v'last loop
      nrm := PentDobl_Complex_Series_Norms.Max_Norm(v(i).all);
      if nrm > res
       then res := nrm;
      end if;
    end loop;
    return res;
  end Max_Norm;
  
end PentDobl_CSeries_Vector_Norms;
