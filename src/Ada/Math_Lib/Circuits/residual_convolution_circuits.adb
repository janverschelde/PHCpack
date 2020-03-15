with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers_Polar;

package body Residual_Convolution_Circuits is

  procedure AbsVal ( c : in out Standard_Complex_Numbers.Complex_Number ) is

    rad : constant double_float := Standard_Complex_Numbers_Polar.Radius(c);

  begin
    c := Standard_Complex_Numbers.Create(rad);
  end AbsVal;

  procedure AbsVal ( c : in out DoblDobl_Complex_Numbers.Complex_Number ) is

    rad : constant double_double := DoblDobl_Complex_Numbers_Polar.Radius(c);

  begin
    c := DoblDobl_Complex_Numbers.Create(rad);
  end AbsVal;

  procedure AbsVal ( c : in out QuadDobl_Complex_Numbers.Complex_Number ) is

    rad : constant quad_double := QuadDobl_Complex_Numbers_Polar.Radius(c);

  begin
    c := QuadDobl_Complex_Numbers.Create(rad);
  end AbsVal;

  procedure AbsVal ( v : in out Standard_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      AbsVal(v(i));
    end loop;
  end AbsVal;

  procedure AbsVal ( v : in out DoblDobl_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      AbsVal(v(i));
    end loop;
  end AbsVal;

  procedure AbsVal ( v : in out QuadDobl_Complex_Vectors.Vector ) is
  begin
    for i in v'range loop
      AbsVal(v(i));
    end loop;
  end AbsVal;

  procedure AbsVal ( v : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Vectors;

  begin
    if v /= null
     then AbsVal(v.all);
    end if;
  end AbsVal;

  procedure AbsVal ( v : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    use DoblDobl_Complex_Vectors;

  begin
    if v /= null
     then AbsVal(v.all);
    end if;
  end AbsVal;

  procedure AbsVal ( v : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    use QuadDobl_Complex_Vectors;

  begin
    if v /= null
     then AbsVal(v.all);
    end if;
  end AbsVal;

  procedure AbsVal ( v : in Standard_Complex_VecVecs.VecVec ) is
  begin
    for i in v'range loop
      AbsVal(v(i));
    end loop;
  end AbsVal;

  procedure AbsVal ( v : in DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in v'range loop
      AbsVal(v(i));
    end loop;
  end AbsVal;

  procedure AbsVal ( v : in QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in v'range loop
      AbsVal(v(i));
    end loop;
  end AbsVal;

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuit )
                  return Standard_Speelpenning_Convolutions.Circuit is

    use Standard_Speelpenning_Convolutions;

    res : Circuit(c.nbr,c.dim,c.dim1,c.dim2);

  begin
    Copy(c,res);
    AbsVal(res.cff);
    AbsVal(res.cst);
    return res;
  end AbsVal;

  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
                  return DoblDobl_Speelpenning_Convolutions.Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    res : Circuit(c.nbr,c.dim,c.dim1,c.dim2);

  begin
    Copy(c,res);
    AbsVal(res.cff);
    AbsVal(res.cst);
    return res;
  end AbsVal;

  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
                  return QuadDobl_Speelpenning_Convolutions.Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    res : Circuit(c.nbr,c.dim,c.dim1,c.dim2);

  begin
    Copy(c,res);
    AbsVal(res.cff);
    AbsVal(res.cst);
    return res;
  end AbsVal;

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit )
                  return Standard_Speelpenning_Convolutions.Link_to_Circuit is

    use Standard_Speelpenning_Convolutions;

    res : Link_to_Circuit;

  begin
    if c /= null then
      declare
        crc : constant Circuit(c.nbr,c.dim,c.dim1,c.dim2) := AbsVal(c.all);
      begin
        res := new Circuit'(crc);
      end;
    end if;
    return res;
  end AbsVal;

  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit )
                  return DoblDobl_Speelpenning_Convolutions.Link_to_Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    res : Link_to_Circuit;

  begin
    if c /= null then
      declare
        crc : constant Circuit(c.nbr,c.dim,c.dim1,c.dim2) := AbsVal(c.all);
      begin
        res := new Circuit'(crc);
      end;
    end if;
    return res;
  end AbsVal;

  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit )
                  return QuadDobl_Speelpenning_Convolutions.Link_to_Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    res : Link_to_Circuit;

  begin
    if c /= null then
      declare
        crc : constant Circuit(c.nbr,c.dim,c.dim1,c.dim2) := AbsVal(c.all);
      begin
        res := new Circuit'(crc);
      end;
    end if;
    return res;
  end AbsVal;

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuits )
                  return Standard_Speelpenning_Convolutions.Circuits is

    res : Standard_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for i in c'range loop
      res(i) := AbsVal(c(i));
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuits )
                  return DoblDobl_Speelpenning_Convolutions.Circuits is

    res : DoblDobl_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for i in c'range loop
      res(i) := AbsVal(c(i));
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
                  return QuadDobl_Speelpenning_Convolutions.Circuits is

    res : QuadDobl_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for i in c'range loop
      res(i) := AbsVal(c(i));
    end loop;
    return res;
  end AbsVal;

end Residual_Convolution_Circuits;
