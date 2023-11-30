with HexaDobl_Complex_Numbers_cv;        use HexaDobl_Complex_Numbers_cv;

package body HexaDobl_Complex_Vectors_cv is

  function Standard_to_HexaDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return HexaDobl_Complex_Vectors.Vector is

    res : HexaDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_to_HexaDobl_Complex(v(i));
    end loop;
    return res;
  end Standard_to_HexaDobl_Complex;

  function Multprec_to_HexaDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return HexaDobl_Complex_Vectors.Vector is

    res : HexaDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Multprec_to_HexaDobl_Complex(v(i));
    end loop;
    return res;
  end Multprec_to_HexaDobl_Complex;

  function HexaDobl_Complex_to_Standard
             ( v : HexaDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_Standard(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_Standard;

  function HexaDobl_Complex_to_Multprec
             ( v : HexaDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_Multprec(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_Multprec;

  function HexaDobl_Complex_to_DoblDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_DoblDobl(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_DoblDobl;

  function HexaDobl_Complex_to_TripDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector is

    res : TripDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_TripDobl(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_TripDobl;

  function HexaDobl_Complex_to_QuadDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_QuadDobl(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_QuadDobl;

  function HexaDobl_Complex_to_PentDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_PentDobl(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_PentDobl;

  function HexaDobl_Complex_to_OctoDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_OctoDobl(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_OctoDobl;

  function HexaDobl_Complex_to_DecaDobl
             ( v : HexaDobl_Complex_Vectors.Vector )
             return DecaDobl_Complex_Vectors.Vector is

    res : DecaDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_to_DecaDobl(v(i));
    end loop;
    return res;
  end HexaDobl_Complex_to_DecaDobl;

  function to_deca_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return DecaDobl_Complex_VecVecs.VecVec is

    res : DecaDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant DecaDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_DecaDobl(v(i).all);
      begin
        res(i) := new DecaDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_deca_double;

  function to_deca_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return DecaDobl_Complex_VecVecs.Link_to_VecVec is

    res : DecaDobl_Complex_VecVecs.Link_to_VecVec;
    tdv : DecaDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in tdv'range loop
      declare
        vec : constant DecaDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_DecaDobl(v(i).all);
      begin
        tdv(i) := new DecaDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new DecaDobl_Complex_VecVecs.VecVec'(tdv);
    return res;
  end to_deca_double;

  function to_octo_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return OctoDobl_Complex_VecVecs.VecVec is

    res : OctoDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant OctoDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_OctoDobl(v(i).all);
      begin
        res(i) := new OctoDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_octo_double;

  function to_octo_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return OctoDobl_Complex_VecVecs.Link_to_VecVec is

    res : OctoDobl_Complex_VecVecs.Link_to_VecVec;
    tdv : OctoDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in tdv'range loop
      declare
        vec : constant OctoDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_OctoDobl(v(i).all);
      begin
        tdv(i) := new OctoDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new OctoDobl_Complex_VecVecs.VecVec'(tdv);
    return res;
  end to_octo_double;

  function to_penta_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return PentDobl_Complex_VecVecs.VecVec is

    res : PentDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant PentDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_PentDobl(v(i).all);
      begin
        res(i) := new PentDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_penta_double;

  function to_penta_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return PentDobl_Complex_VecVecs.Link_to_VecVec is

    res : PentDobl_Complex_VecVecs.Link_to_VecVec;
    tdv : PentDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in tdv'range loop
      declare
        vec : constant PentDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_PentDobl(v(i).all);
      begin
        tdv(i) := new PentDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new PentDobl_Complex_VecVecs.VecVec'(tdv);
    return res;
  end to_penta_double;

  function to_quad_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant QuadDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_QuadDobl(v(i).all);
      begin
        res(i) := new QuadDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_quad_double;

  function to_quad_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return QuadDobl_Complex_VecVecs.Link_to_VecVec is

    res : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    tdv : QuadDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in tdv'range loop
      declare
        vec : constant QuadDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_QuadDobl(v(i).all);
      begin
        tdv(i) := new QuadDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new QuadDobl_Complex_VecVecs.VecVec'(tdv);
    return res;
  end to_quad_double;

  function to_triple_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return TripDobl_Complex_VecVecs.VecVec is

    res : TripDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant TripDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_TripDobl(v(i).all);
      begin
        res(i) := new TripDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_triple_double;

  function to_triple_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return TripDobl_Complex_VecVecs.Link_to_VecVec is

    res : TripDobl_Complex_VecVecs.Link_to_VecVec;
    tdv : TripDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in tdv'range loop
      declare
        vec : constant TripDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_TripDobl(v(i).all);
      begin
        tdv(i) := new TripDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new TripDobl_Complex_VecVecs.VecVec'(tdv);
    return res;
  end to_triple_double;

  function to_double_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant DoblDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_DoblDobl(v(i).all);
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_double_double;

  function to_double_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return DoblDobl_Complex_VecVecs.Link_to_VecVec is

    res : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    ddv : DoblDobl_Complex_VecVecs.VecVec(v'range);

  begin
    for i in ddv'range loop
      declare
        vec : constant DoblDobl_Complex_Vectors.Vector
            := HexaDobl_Complex_to_DoblDobl(v(i).all);
      begin
        ddv(i) := new DoblDobl_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new DoblDobl_Complex_VecVecs.VecVec'(ddv);
    return res;
  end to_double_double;

  function to_double
             ( v : HexaDobl_Complex_VecVecs.VecVec )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(v'range);

  begin
    for i in res'range loop
      declare
        vec : constant Standard_Complex_Vectors.Vector
            := HexaDobl_Complex_to_Standard(v(i).all);
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    return res;
  end to_double;

  function to_double
             ( v : HexaDobl_Complex_VecVecs.Link_to_VecVec )
             return Standard_Complex_VecVecs.Link_to_VecVec is

    res : Standard_Complex_VecVecs.Link_to_VecVec;
    dv : Standard_Complex_VecVecs.VecVec(v'range);

  begin
    for i in dv'range loop
      declare
        vec : constant Standard_Complex_Vectors.Vector
            := HexaDobl_Complex_to_Standard(v(i).all);
      begin
        dv(i) := new Standard_Complex_Vectors.Vector'(vec);
      end;
    end loop;
    res := new Standard_Complex_VecVecs.VecVec'(dv);
    return res;
  end to_double;

end HexaDobl_Complex_Vectors_cv;
