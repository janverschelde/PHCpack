with QuadDobl_Complex_Numbers_cv;        use QuadDobl_Complex_Numbers_cv;

package body QuadDobl_Complex_Vectors_cv is

  function Standard_to_QuadDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_to_QuadDobl_Complex(v(i));
    end loop;
    return res;
  end Standard_to_QuadDobl_Complex;

  function Multprec_to_QuadDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Multprec_to_QuadDobl_Complex(v(i));
    end loop;
    return res;
  end Multprec_to_QuadDobl_Complex;

  function QuadDobl_Complex_to_Standard
             ( v : QuadDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := QuadDobl_Complex_to_Standard(v(i));
    end loop;
    return res;
  end QuadDobl_Complex_to_Standard;

  function QuadDobl_Complex_to_DoblDobl
             ( v : QuadDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := QuadDobl_Complex_to_DoblDobl(v(i));
    end loop;
    return res;
  end QuadDobl_Complex_to_DoblDobl;

  function QuadDobl_Complex_to_TripDobl
             ( v : QuadDobl_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector is

    res : TripDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := QuadDobl_Complex_to_TripDobl(v(i));
    end loop;
    return res;
  end QuadDobl_Complex_to_TripDobl;

  function QuadDobl_Complex_to_Multprec
             ( v : QuadDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := QuadDobl_Complex_to_Multprec(v(i));
    end loop;
    return res;
  end QuadDobl_Complex_to_Multprec;

  function to_triple_double
             ( v : QuadDobl_Complex_VecVecs.VecVec )
             return TripDobl_Complex_VecVecs.VecVec is

    res : TripDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant TripDobl_Complex_Vectors.Vector
            := QuadDobl_Complex_to_TripDobl(v(i).all);
      begin
        res(i) := new TripDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_triple_double;

  function to_triple_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return TripDobl_Complex_VecVecs.Link_to_VecVec is

    res : TripDobl_Complex_VecVecs.Link_to_VecVec;
    tdv : TripDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in tdv'range loop
      declare
        vec : constant TripDobl_Complex_Vectors.Vector
            := QuadDobl_Complex_to_TripDobl(v(i).all);
      begin
        tdv(i) := new TripDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new TripDobl_Complex_VecVecs.VecVec'(tdv);
    return res;
  end to_triple_double;

  function to_double_double
             ( v : QuadDobl_Complex_VecVecs.VecVec )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant DoblDobl_Complex_Vectors.Vector
            := QuadDobl_Complex_to_DoblDobl(v(i).all);
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_double_double;

  function to_double_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return DoblDobl_Complex_VecVecs.Link_to_VecVec is

    res : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    ddv : DoblDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in ddv'range loop
      declare
        vec : constant DoblDobl_Complex_Vectors.Vector
            := QuadDobl_Complex_to_DoblDobl(v(i).all);
      begin
        ddv(i) := new DoblDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new DoblDobl_Complex_VecVecs.VecVec'(ddv);
    return res;
  end to_double_double;

  function to_double
             ( v : QuadDobl_Complex_VecVecs.VecVec )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant Standard_Complex_Vectors.Vector
            := QuadDobl_Complex_to_Standard(v(i).all);
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_double;

  function to_double
             ( v : QuadDobl_Complex_VecVecs.Link_to_VecVec )
             return Standard_Complex_VecVecs.Link_to_VecVec is

    res : Standard_Complex_VecVecs.Link_to_VecVec;
    dv : Standard_Complex_VecVecs.VecVec(v'range);

  begin
    for i in dv'range loop
      declare
        vec : constant Standard_Complex_Vectors.Vector
            := QuadDobl_Complex_to_Standard(v(i).all);
      begin
        dv(i) := new Standard_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new Standard_Complex_VecVecs.VecVec'(dv);
    return res;
  end to_double;

end QuadDobl_Complex_Vectors_cv;
