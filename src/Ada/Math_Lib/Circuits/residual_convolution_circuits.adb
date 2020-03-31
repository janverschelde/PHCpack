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

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuit;
                    deg : integer32 := 0 )
                  return Standard_Speelpenning_Convolutions.Circuit is

    use Standard_Speelpenning_Convolutions;

    res : Circuit(c.nbr,c.dim,c.dim1,c.dim2);

  begin
    res := Copy(c,deg); -- Copy(c,res);
    AbsVal(res.cff);
    AbsVal(res.cst);
    return res;
  end AbsVal;

  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuit;
                    deg : integer32 := 0 )
                  return DoblDobl_Speelpenning_Convolutions.Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    res : Circuit(c.nbr,c.dim,c.dim1,c.dim2);

  begin
    res := Copy(c,deg); -- Copy(c,res);
    AbsVal(res.cff);
    AbsVal(res.cst);
    return res;
  end AbsVal;

  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuit;
                    deg : integer32 := 0 )
                  return QuadDobl_Speelpenning_Convolutions.Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    res : Circuit(c.nbr,c.dim,c.dim1,c.dim2);

  begin
    res := Copy(c,deg); -- Copy(c,res);
    AbsVal(res.cff);
    AbsVal(res.cst);
    return res;
  end AbsVal;

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit;
                    deg : integer32 := 0 )
                  return Standard_Speelpenning_Convolutions.Link_to_Circuit is

    use Standard_Speelpenning_Convolutions;

    res : Link_to_Circuit;

  begin
    if c /= null then
      declare
        crc : constant Circuit(c.nbr,c.dim,c.dim1,c.dim2)
            := AbsVal(c.all,deg);
      begin
        res := new Circuit'(crc);
      end;
    end if;
    return res;
  end AbsVal;

  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    deg : integer32 := 0 )
                  return DoblDobl_Speelpenning_Convolutions.Link_to_Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    res : Link_to_Circuit;

  begin
    if c /= null then
      declare
        crc : constant Circuit(c.nbr,c.dim,c.dim1,c.dim2)
            := AbsVal(c.all,deg);
      begin
        res := new Circuit'(crc);
      end;
    end if;
    return res;
  end AbsVal;

  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                    deg : integer32 := 0 )
                  return QuadDobl_Speelpenning_Convolutions.Link_to_Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    res : Link_to_Circuit;

  begin
    if c /= null then
      declare
        crc : constant Circuit(c.nbr,c.dim,c.dim1,c.dim2)
            := AbsVal(c.all,deg);
      begin
        res := new Circuit'(crc);
      end;
    end if;
    return res;
  end AbsVal;

  function AbsVal ( c : Standard_Speelpenning_Convolutions.Circuits;
                    deg : integer32 := 0 )
                  return Standard_Speelpenning_Convolutions.Circuits is

    res : Standard_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for i in c'range loop
      res(i) := AbsVal(c(i),deg);
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( c : DoblDobl_Speelpenning_Convolutions.Circuits;
                    deg : integer32 := 0 )
                  return DoblDobl_Speelpenning_Convolutions.Circuits is

    res : DoblDobl_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for i in c'range loop
      res(i) := AbsVal(c(i),deg);
    end loop;
    return res;
  end AbsVal;

  function AbsVal ( c : QuadDobl_Speelpenning_Convolutions.Circuits;
                    deg : integer32 := 0 )
                  return QuadDobl_Speelpenning_Convolutions.Circuits is

    res : QuadDobl_Speelpenning_Convolutions.Circuits(c'range);

  begin
    for i in c'range loop
      res(i) := AbsVal(c(i),deg);
    end loop;
    return res;
  end AbsVal;

  function Residual_Convolution_System
             ( s : Standard_Speelpenning_Convolutions.System;
               deg : integer32 := 0 )
             return Standard_Speelpenning_Convolutions.System is

    use  Standard_Speelpenning_Convolutions;

    crc : constant Circuits(s.crc'range) := AbsVal(s.crc,deg);
    res : constant System(s.neq,s.neq1,s.dim,s.dim1,deg)
        := Create(crc,s.dim,deg);

  begin
    return res;
  end Residual_Convolution_System;

  function Residual_Convolution_System
             ( s : DoblDobl_Speelpenning_Convolutions.System;
               deg : integer32 := 0 )
             return DoblDobl_Speelpenning_Convolutions.System is

    use  DoblDobl_Speelpenning_Convolutions;

    crc : constant Circuits(s.crc'range) := AbsVal(s.crc,deg);
    res : constant System(s.neq,s.neq1,s.dim,s.dim1,deg)
        := Create(crc,s.dim,deg);

  begin
    return res;
  end Residual_Convolution_System;

  function Residual_Convolution_System
             ( s : QuadDobl_Speelpenning_Convolutions.System;
               deg : integer32 := 0 )
             return QuadDobl_Speelpenning_Convolutions.System is

    use  QuadDobl_Speelpenning_Convolutions;

    crc : constant Circuits(s.crc'range) := AbsVal(s.crc,deg);
    res : constant System(s.neq,s.neq1,s.dim,s.dim1,deg)
        := Create(crc,s.dim,deg);

  begin
    return res;
  end Residual_Convolution_System;

  function Residual_Convolution_System
             ( s : Standard_Speelpenning_Convolutions.Link_to_System;
               deg : integer32 := 0 )
             return Standard_Speelpenning_Convolutions.Link_to_System is

    use  Standard_Speelpenning_Convolutions;

    res : Link_to_System;

  begin
    if s /= null then
      declare
        sys : constant System(s.neq,s.neq1,s.dim,s.dim1,deg)
            := Residual_Convolution_System(s.all,deg);
      begin
        res := new System'(sys);
      end;
    end if;
    return res;
  end Residual_Convolution_System;

  function Residual_Convolution_System
             ( s : DoblDobl_Speelpenning_Convolutions.Link_to_System;
               deg : integer32 := 0 )
             return DoblDobl_Speelpenning_Convolutions.Link_to_System is

    use  DoblDobl_Speelpenning_Convolutions;

    res : Link_to_System;

  begin
    if s /= null then
      declare
        sys : constant System(s.neq,s.neq1,s.dim,s.dim1,deg)
            := Residual_Convolution_System(s.all,deg);
      begin
        res := new System'(sys);
      end;
    end if;
    return res;
  end Residual_Convolution_System;

  function Residual_Convolution_System
             ( s : QuadDobl_Speelpenning_Convolutions.Link_to_System;
               deg : integer32 := 0 )
             return QuadDobl_Speelpenning_Convolutions.Link_to_System is

    use  QuadDobl_Speelpenning_Convolutions;

    res : Link_to_System;

  begin
    if s /= null then
      declare
        sys : constant System(s.neq,s.neq1,s.dim,s.dim1,deg)
            := Residual_Convolution_System(s.all,deg);
      begin
        res := new System'(sys);
      end;
    end if;
    return res;
  end Residual_Convolution_System;

  procedure Update_Radii_of_Constants
              ( radc : in Standard_Speelpenning_Convolutions.Circuit; 
                c : in Standard_Speelpenning_Convolutions.Circuit ) is

    radlnk,lnk : Standard_Complex_Vectors.Link_to_Vector;
    rad : double_float;

    use Standard_Complex_Vectors;

  begin
    for k in radc.cff'range loop
      radlnk := radc.cff(k); lnk := c.cff(k);
      rad := Standard_Complex_Numbers_Polar.Radius(lnk(0));
      radlnk(0) := Standard_Complex_Numbers.Create(rad);
    end loop;
    if radc.cst /= null and c.cst /= null then
      rad := Standard_Complex_Numbers_Polar.Radius(c.cst(0));
      radc.cst(0) := Standard_Complex_Numbers.Create(rad);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in DoblDobl_Speelpenning_Convolutions.Circuit; 
                c : in DoblDobl_Speelpenning_Convolutions.Circuit ) is

    radlnk,lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    rad : double_double;

    use DoblDobl_Complex_Vectors;

  begin
    for k in radc.cff'range loop
      radlnk := radc.cff(k); lnk := c.cff(k);
      rad := DoblDobl_Complex_Numbers_Polar.Radius(lnk(0));
      radlnk(0) := DoblDobl_Complex_Numbers.Create(rad);
    end loop;
    if radc.cst /= null and c.cst /= null then
      rad := DoblDobl_Complex_Numbers_Polar.Radius(c.cst(0));
      radc.cst(0) := DoblDobl_Complex_Numbers.Create(rad);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in QuadDobl_Speelpenning_Convolutions.Circuit; 
                c : in QuadDobl_Speelpenning_Convolutions.Circuit ) is

    radlnk,lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    rad : quad_double;

    use QuadDobl_Complex_Vectors;

  begin
    for k in radc.cff'range loop
      radlnk := radc.cff(k); lnk := c.cff(k);
      rad := QuadDobl_Complex_Numbers_Polar.Radius(lnk(0));
      radlnk(0) := QuadDobl_Complex_Numbers.Create(rad);
    end loop;
    if radc.cst /= null and c.cst /= null then
      rad := QuadDobl_Complex_Numbers_Polar.Radius(c.cst(0));
      radc.cst(0) := QuadDobl_Complex_Numbers.Create(rad);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in Standard_Speelpenning_Convolutions.Link_to_Circuit; 
                c : in Standard_Speelpenning_Convolutions.Link_to_Circuit ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if radc /= null and c /= null
     then Update_Radii_of_Constants(radc.all,c.all);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit; 
                c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if radc /= null and c /= null
     then Update_Radii_of_Constants(radc.all,c.all);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit; 
                c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if radc /= null and c /= null
     then Update_Radii_of_Constants(radc.all,c.all);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in Standard_Speelpenning_Convolutions.Circuits; 
                c : in Standard_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in radc'range loop
      Update_Radii_of_Constants(radc(k),c(k));
    end loop;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in DoblDobl_Speelpenning_Convolutions.Circuits; 
                c : in DoblDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in radc'range loop
      Update_Radii_of_Constants(radc(k),c(k));
    end loop;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( radc : in QuadDobl_Speelpenning_Convolutions.Circuits; 
                c : in QuadDobl_Speelpenning_Convolutions.Circuits ) is
  begin
    for k in radc'range loop
      Update_Radii_of_Constants(radc(k),c(k));
    end loop;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( rads : in Standard_Speelpenning_Convolutions.Link_to_System; 
                s : in Standard_Speelpenning_Convolutions.Link_to_System ) is

    use Standard_Speelpenning_Convolutions;

  begin
    if rads /= null and s /= null
     then Update_Radii_of_Constants(rads.crc,s.crc);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( rads : in DoblDobl_Speelpenning_Convolutions.Link_to_System; 
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if rads /= null and s /= null
     then Update_Radii_of_Constants(rads.crc,s.crc);
    end if;
  end Update_Radii_of_Constants;

  procedure Update_Radii_of_Constants
              ( rads : in QuadDobl_Speelpenning_Convolutions.Link_to_System; 
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if rads /= null and s /= null
     then Update_Radii_of_Constants(rads.crc,s.crc);
    end if;
  end Update_Radii_of_Constants;

end Residual_Convolution_Circuits;
