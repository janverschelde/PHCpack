with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Series;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Random_Series;
with Standard_Random_Series_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Random_Series;
with DoblDobl_Random_Series_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Random_Series;
with QuadDobl_Random_Series_Vectors;
with Exponent_Indices;
with Series_Polynomial_Gradients;

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
             return Standard_Speelpenning_Convolutions.Convolution_Circuit is

    use Standard_Speelpenning_Convolutions;

    res : Convolution_Circuit(nbr,dim,dim-1,dim-2);
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
    res.cff := Series_Polynomial_Gradients.Standard_Series_Coefficients(polcff);
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
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuit is

    use DoblDobl_Speelpenning_Convolutions;

    res : Convolution_Circuit(nbr,dim,dim-1,dim-2);
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
    res.cff := Series_Polynomial_Gradients.DoblDobl_Series_Coefficients(polcff);
    res.cst := new DoblDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end DoblDobl_Random_Convolution_Circuit;

  function QuadDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuit is

    use QuadDobl_Speelpenning_Convolutions;

    res : Convolution_Circuit(nbr,dim,dim-1,dim-2);
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
    res.cff := Series_Polynomial_Gradients.QuadDobl_Series_Coefficients(polcff);
    res.cst := new QuadDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end QuadDobl_Random_Convolution_Circuit;

  function Standard_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return Standard_Speelpenning_Convolutions.Convolution_Circuits is

    use Standard_Speelpenning_Convolutions;

    res : Convolution_Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Convolution_Circuit(nbr,dim,dim-1,dim-2)
          := Standard_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Convolution_Circuit'(c);
      end;
    end loop;
    return res;
  end Standard_Random_Convolution_Circuits;

  function DoblDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuits is

    use DoblDobl_Speelpenning_Convolutions;

    res : Convolution_Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Convolution_Circuit(nbr,dim,dim-1,dim-2)
          := DoblDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Convolution_Circuit'(c);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Convolution_Circuits;

  function QuadDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuits is

    use QuadDobl_Speelpenning_Convolutions;

    res : Convolution_Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Convolution_Circuit(nbr,dim,dim-1,dim-2)
          := QuadDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Convolution_Circuit'(c);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Convolution_Circuits;

end Random_Convolution_Circuits;
