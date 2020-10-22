with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with TripDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with PentDobl_Complex_Vectors;
with OctoDobl_Complex_Vectors;
with DecaDobl_Complex_Vectors;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with TripDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with PentDobl_Random_Vectors;
with OctoDobl_Random_Vectors;
with DecaDobl_Random_Vectors;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Random_Series;
with Standard_Random_Series_Vectors;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Random_Series;
with DoblDobl_Random_Series_Vectors;
with TripDobl_Complex_Series;
with TripDobl_Complex_Series_Vectors;
with TripDobl_Complex_Random_Series;
with TripDobl_Random_Series_Vectors;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Random_Series;
with QuadDobl_Random_Series_Vectors;
with PentDobl_Complex_Series;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Random_Series;
with PentDobl_Random_Series_Vectors;
with OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Random_Series;
with OctoDobl_Random_Series_Vectors;
with DecaDobl_Complex_Series;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Random_Series;
with DecaDobl_Random_Series_Vectors;
with Exponent_Indices;
with Series_Coefficient_Vectors;
with Homotopy_Convolution_Circuits;
with Standard_Newton_Convolutions;
with DoblDobl_Newton_Convolutions;
with TripDobl_Newton_Convolutions;
with QuadDobl_Newton_Convolutions;
with PentDobl_Newton_Convolutions;
with OctoDobl_Newton_Convolutions;
with DecaDobl_Newton_Convolutions;

package body Random_Convolution_Circuits is

  function Random_Exponents
             ( dim,nbr : integer32; pwr : integer32 := 1 ) 
             return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..nbr);
    nz : integer32;

  begin
    for i in 1..nbr loop
      loop
        declare
          xp : constant Standard_Integer_Vectors.Vector(1..dim)
             := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
        begin
          nz := Standard_Integer_Vectors.Sum(xp);
          if nz > 0
           then res(i) := new Standard_Integer_Vectors.Vector'(xp);
          end if;
        end;
        exit when (nz > 0);
      end loop;
    end loop;
    return res;
  end Random_Exponents;

  function Standard_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return Standard_Speelpenning_Convolutions.Circuit is

    use Standard_Speelpenning_Convolutions;

    res : Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant Standard_Complex_Series_Vectors.Vector(1..nbr)
           := Standard_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    rancst : constant Standard_Complex_Series.Series
           := Standard_Complex_Random_Series.Random_Series(deg);
    cstcff : constant Standard_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Series_Coefficient_Vectors.Standard_Series_Coefficients(polcff);
    res.cst := new Standard_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end Standard_Random_Convolution_Circuit;

  function DoblDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return DoblDobl_Speelpenning_Convolutions.Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    res : Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant DoblDobl_Complex_Series_Vectors.Vector(1..nbr)
           := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    rancst : constant DoblDobl_Complex_Series.Series
           := DoblDobl_Complex_Random_Series.Random_Series(deg);
    cstcff : constant DoblDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Series_Coefficient_Vectors.DoblDobl_Series_Coefficients(polcff);
    res.cst := new DoblDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end DoblDobl_Random_Convolution_Circuit;

  function TripDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return TripDobl_Speelpenning_Convolutions.Circuit is

    use TripDobl_Speelpenning_Convolutions;

    res : Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant TripDobl_Complex_Series_Vectors.Vector(1..nbr)
           := TripDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    rancst : constant TripDobl_Complex_Series.Series
           := TripDobl_Complex_Random_Series.Random_Series(deg);
    cstcff : constant TripDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Series_Coefficient_Vectors.TripDobl_Series_Coefficients(polcff);
    res.cst := new TripDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end TripDobl_Random_Convolution_Circuit;

  function QuadDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    res : Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant QuadDobl_Complex_Series_Vectors.Vector(1..nbr)
           := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    rancst : constant QuadDobl_Complex_Series.Series
           := QuadDobl_Complex_Random_Series.Random_Series(deg);
    cstcff : constant QuadDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Series_Coefficient_Vectors.QuadDobl_Series_Coefficients(polcff);
    res.cst := new QuadDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end QuadDobl_Random_Convolution_Circuit;

  function PentDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return PentDobl_Speelpenning_Convolutions.Circuit is

    use PentDobl_Speelpenning_Convolutions;

    res : Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant PentDobl_Complex_Series_Vectors.Vector(1..nbr)
           := PentDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    rancst : constant PentDobl_Complex_Series.Series
           := PentDobl_Complex_Random_Series.Random_Series(deg);
    cstcff : constant PentDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Series_Coefficient_Vectors.PentDobl_Series_Coefficients(polcff);
    res.cst := new PentDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end PentDobl_Random_Convolution_Circuit;

  function OctoDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return OctoDobl_Speelpenning_Convolutions.Circuit is

    use OctoDobl_Speelpenning_Convolutions;

    res : Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant OctoDobl_Complex_Series_Vectors.Vector(1..nbr)
           := OctoDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    rancst : constant OctoDobl_Complex_Series.Series
           := OctoDobl_Complex_Random_Series.Random_Series(deg);
    cstcff : constant OctoDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Series_Coefficient_Vectors.OctoDobl_Series_Coefficients(polcff);
    res.cst := new OctoDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end OctoDobl_Random_Convolution_Circuit;

  function DecaDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return DecaDobl_Speelpenning_Convolutions.Circuit is

    use DecaDobl_Speelpenning_Convolutions;

    res : Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant DecaDobl_Complex_Series_Vectors.Vector(1..nbr)
           := DecaDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    rancst : constant DecaDobl_Complex_Series.Series
           := DecaDobl_Complex_Random_Series.Random_Series(deg);
    cstcff : constant DecaDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Series_Coefficient_Vectors.DecaDobl_Series_Coefficients(polcff);
    res.cst := new DecaDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end DecaDobl_Random_Convolution_Circuit;

  function Standard_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return Standard_Speelpenning_Convolutions.Circuits is

    use Standard_Speelpenning_Convolutions;

    res : Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Circuit(nbr,dim,dim-1,dim-2)
          := Standard_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Circuit'(c);
      end;
    end loop;
    return res;
  end Standard_Random_Convolution_Circuits;

  function DoblDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return DoblDobl_Speelpenning_Convolutions.Circuits is

    use DoblDobl_Speelpenning_Convolutions;

    res : Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Circuit(nbr,dim,dim-1,dim-2)
          := DoblDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Circuit'(c);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Convolution_Circuits;

  function TripDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return TripDobl_Speelpenning_Convolutions.Circuits is

    use TripDobl_Speelpenning_Convolutions;

    res : Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Circuit(nbr,dim,dim-1,dim-2)
          := TripDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Circuit'(c);
      end;
    end loop;
    return res;
  end TripDobl_Random_Convolution_Circuits;

  function QuadDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Circuits is

    use QuadDobl_Speelpenning_Convolutions;

    res : Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Circuit(nbr,dim,dim-1,dim-2)
          := QuadDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Circuit'(c);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Convolution_Circuits;

  function PentDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return PentDobl_Speelpenning_Convolutions.Circuits is

    use PentDobl_Speelpenning_Convolutions;

    res : Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Circuit(nbr,dim,dim-1,dim-2)
          := PentDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Circuit'(c);
      end;
    end loop;
    return res;
  end PentDobl_Random_Convolution_Circuits;

  function OctoDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return OctoDobl_Speelpenning_Convolutions.Circuits is

    use OctoDobl_Speelpenning_Convolutions;

    res : Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Circuit(nbr,dim,dim-1,dim-2)
          := OctoDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Circuit'(c);
      end;
    end loop;
    return res;
  end OctoDobl_Random_Convolution_Circuits;

  function DecaDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return DecaDobl_Speelpenning_Convolutions.Circuits is

    use DecaDobl_Speelpenning_Convolutions;

    res : Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Circuit(nbr,dim,dim-1,dim-2)
          := DecaDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Circuit'(c);
      end;
    end loop;
    return res;
  end DecaDobl_Random_Convolution_Circuits;

  function Standard_Random_System
             ( dim,deg,nbr,pwr : integer32 )
             return Standard_Speelpenning_Convolutions.Link_to_System is

    c : constant Standard_Speelpenning_Convolutions.Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    return Standard_Speelpenning_Convolutions.Create(c,dim,deg);
  end Standard_Random_System;

  function DoblDobl_Random_System
             ( dim,deg,nbr,pwr : integer32 )
             return DoblDobl_Speelpenning_Convolutions.Link_to_System is

    c : constant DoblDobl_Speelpenning_Convolutions.Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    return DoblDobl_Speelpenning_Convolutions.Create(c,dim,deg);
  end DoblDobl_Random_System;

  function TripDobl_Random_System
             ( dim,deg,nbr,pwr : integer32 )
             return TripDobl_Speelpenning_Convolutions.Link_to_System is

    c : constant TripDobl_Speelpenning_Convolutions.Circuits
      := TripDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    return TripDobl_Speelpenning_Convolutions.Create(c,dim,deg);
  end TripDobl_Random_System;

  function QuadDobl_Random_System
             ( dim,deg,nbr,pwr : integer32 )
             return QuadDobl_Speelpenning_Convolutions.Link_to_System is

    c : constant QuadDobl_Speelpenning_Convolutions.Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    return QuadDobl_Speelpenning_Convolutions.Create(c,dim,deg);
  end QuadDobl_Random_System;

  function PentDobl_Random_System
             ( dim,deg,nbr,pwr : integer32 )
             return PentDobl_Speelpenning_Convolutions.Link_to_System is

    c : constant PentDobl_Speelpenning_Convolutions.Circuits
      := PentDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    return PentDobl_Speelpenning_Convolutions.Create(c,dim,deg);
  end PentDobl_Random_System;

  function OctoDobl_Random_System
             ( dim,deg,nbr,pwr : integer32 )
             return OctoDobl_Speelpenning_Convolutions.Link_to_System is

    c : constant OctoDobl_Speelpenning_Convolutions.Circuits
      := OctoDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    return OctoDobl_Speelpenning_Convolutions.Create(c,dim,deg);
  end OctoDobl_Random_System;

  function DecaDobl_Random_System
             ( dim,deg,nbr,pwr : integer32 )
             return DecaDobl_Speelpenning_Convolutions.Link_to_System is

    c : constant DecaDobl_Speelpenning_Convolutions.Circuits
      := DecaDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);

  begin
    return DecaDobl_Speelpenning_Convolutions.Create(c,dim,deg);
  end DecaDobl_Random_System;

  procedure Standard_Random_Newton_Homotopy
             ( dim,deg,nbr,pwr : in integer32;
               s : out Standard_Speelpenning_Convolutions.Link_to_System;
               x : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    z : constant Standard_Complex_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    sz : constant Standard_Complex_VecVecs.VecVec(1..dim)
       := Standard_Newton_Convolutions.Series_Coefficients(z,deg);

  begin
    x := new Standard_Complex_VecVecs.VecVec'(sz);
    s := Standard_Random_System(dim,deg,nbr,pwr);
    Homotopy_Convolution_Circuits.Newton_Homotopy(s.crc,z);
  end Standard_Random_Newton_Homotopy;

  procedure DoblDobl_Random_Newton_Homotopy
             ( dim,deg,nbr,pwr : in integer32;
               s : out DoblDobl_Speelpenning_Convolutions.Link_to_System;
               x : out DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    z : constant DoblDobl_Complex_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    sz : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
       := DoblDobl_Newton_Convolutions.Series_Coefficients(z,deg);

  begin
    x := new DoblDobl_Complex_VecVecs.VecVec'(sz);
    s := DoblDobl_Random_System(dim,deg,nbr,pwr);
    Homotopy_Convolution_Circuits.Newton_Homotopy(s.crc,z);
  end DoblDobl_Random_Newton_Homotopy;

  procedure TripDobl_Random_Newton_Homotopy
             ( dim,deg,nbr,pwr : in integer32;
               s : out TripDobl_Speelpenning_Convolutions.Link_to_System;
               x : out TripDobl_Complex_VecVecs.Link_to_VecVec ) is

    z : constant TripDobl_Complex_Vectors.Vector(1..dim)
      := TripDobl_Random_Vectors.Random_Vector(1,dim);
    sz : constant TripDobl_Complex_VecVecs.VecVec(1..dim)
       := TripDobl_Newton_Convolutions.Series_Coefficients(z,deg);

  begin
    x := new TripDobl_Complex_VecVecs.VecVec'(sz);
    s := TripDobl_Random_System(dim,deg,nbr,pwr);
    Homotopy_Convolution_Circuits.Newton_Homotopy(s.crc,z);
  end TripDobl_Random_Newton_Homotopy;

  procedure QuadDobl_Random_Newton_Homotopy
             ( dim,deg,nbr,pwr : in integer32;
               s : out QuadDobl_Speelpenning_Convolutions.Link_to_System;
               x : out QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    z : constant QuadDobl_Complex_Vectors.Vector(1..dim)
      := QuadDobl_Random_Vectors.Random_Vector(1,dim);
    sz : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
       := QuadDobl_Newton_Convolutions.Series_Coefficients(z,deg);

  begin
    x := new QuadDobl_Complex_VecVecs.VecVec'(sz);
    s := QuadDobl_Random_System(dim,deg,nbr,pwr);
    Homotopy_Convolution_Circuits.Newton_Homotopy(s.crc,z);
  end QuadDobl_Random_Newton_Homotopy;

  procedure PentDobl_Random_Newton_Homotopy
             ( dim,deg,nbr,pwr : in integer32;
               s : out PentDobl_Speelpenning_Convolutions.Link_to_System;
               x : out PentDobl_Complex_VecVecs.Link_to_VecVec ) is

    z : constant PentDobl_Complex_Vectors.Vector(1..dim)
      := PentDobl_Random_Vectors.Random_Vector(1,dim);
    sz : constant PentDobl_Complex_VecVecs.VecVec(1..dim)
       := PentDobl_Newton_Convolutions.Series_Coefficients(z,deg);

  begin
    x := new PentDobl_Complex_VecVecs.VecVec'(sz);
    s := PentDobl_Random_System(dim,deg,nbr,pwr);
    Homotopy_Convolution_Circuits.Newton_Homotopy(s.crc,z);
  end PentDobl_Random_Newton_Homotopy;

  procedure OctoDobl_Random_Newton_Homotopy
             ( dim,deg,nbr,pwr : in integer32;
               s : out OctoDobl_Speelpenning_Convolutions.Link_to_System;
               x : out OctoDobl_Complex_VecVecs.Link_to_VecVec ) is

    z : constant OctoDobl_Complex_Vectors.Vector(1..dim)
      := OctoDobl_Random_Vectors.Random_Vector(1,dim);
    sz : constant OctoDobl_Complex_VecVecs.VecVec(1..dim)
       := OctoDobl_Newton_Convolutions.Series_Coefficients(z,deg);

  begin
    x := new OctoDobl_Complex_VecVecs.VecVec'(sz);
    s := OctoDobl_Random_System(dim,deg,nbr,pwr);
    Homotopy_Convolution_Circuits.Newton_Homotopy(s.crc,z);
  end OctoDobl_Random_Newton_Homotopy;

  procedure DecaDobl_Random_Newton_Homotopy
             ( dim,deg,nbr,pwr : in integer32;
               s : out DecaDobl_Speelpenning_Convolutions.Link_to_System;
               x : out DecaDobl_Complex_VecVecs.Link_to_VecVec ) is

    z : constant DecaDobl_Complex_Vectors.Vector(1..dim)
      := DecaDobl_Random_Vectors.Random_Vector(1,dim);
    sz : constant DecaDobl_Complex_VecVecs.VecVec(1..dim)
       := DecaDobl_Newton_Convolutions.Series_Coefficients(z,deg);

  begin
    x := new DecaDobl_Complex_VecVecs.VecVec'(sz);
    s := DecaDobl_Random_System(dim,deg,nbr,pwr);
    Homotopy_Convolution_Circuits.Newton_Homotopy(s.crc,z);
  end DecaDobl_Random_Newton_Homotopy;

end Random_Convolution_Circuits;
