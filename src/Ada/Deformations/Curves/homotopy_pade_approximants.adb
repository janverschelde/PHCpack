with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;
with Homotopy_Series_Readers;           use Homotopy_Series_Readers;

package body Homotopy_Pade_Approximants is

  procedure Standard_Pade_Approximant
              ( sol : in Standard_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out Standard_Dense_Series_Vectors.Vector;
                pv : out Standard_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false ) is

    nbterms : constant natural32 := natural32(numdeg+dendeg+1);

  begin
    Standard_Series_Newton(sol,idx,nbequ,nbterms,nbiters,srv,eva);
    pv := Standard_Pade_Approximants.Create(numdeg,dendeg,srv,verbose);
  end Standard_Pade_Approximant;

  procedure DoblDobl_Pade_Approximant
              ( sol : in DoblDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out DoblDobl_Dense_Series_Vectors.Vector;
                pv : out DoblDobl_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false ) is

    nbterms : constant natural32 := natural32(numdeg+dendeg+1);

  begin
    DoblDobl_Series_Newton(sol,idx,nbequ,nbterms,nbiters,srv,eva);
    pv := DoblDobl_Pade_Approximants.Create(numdeg,dendeg,srv);
  end DoblDobl_Pade_Approximant;

  procedure QuadDobl_Pade_Approximant
              ( sol : in QuadDobl_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out QuadDobl_Dense_Series_Vectors.Vector;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false ) is

    nbterms : constant natural32 := natural32(numdeg+dendeg+1);

  begin
    QuadDobl_Series_Newton(sol,idx,nbequ,nbterms,nbiters,srv,eva);
    pv := QuadDobl_Pade_Approximants.Create(numdeg,dendeg,srv);
  end QuadDobl_Pade_Approximant;

  function Standard_Poles
              ( p : Standard_Pade_Approximants.Pade )
              return Standard_Complex_Vectors.Vector is

    deg : constant integer32
        := Standard_Pade_Approximants.Denominator_Degree(p);
    res : Standard_Complex_Vectors.Vector(1..deg);
    cff : constant Standard_Complex_Vectors.Vector
        := Standard_Pade_Approximants.Denominator_Coefficients(p);
    dsc,sqrtdsc,den2cff2 : Standard_Complex_Numbers.Complex_Number;
    two : constant double_float := 2.0;
    four : constant double_float := 4.0;

  begin
    if deg = 1 then
      res(1) := -cff(0)/cff(1);
    elsif deg = 2 then
      dsc := cff(1)**2 - four*cff(0)*cff(2);
      sqrtdsc := Standard_Complex_Numbers_Polar.Root(dsc,2,1);
      den2cff2 := two*cff(2);
      res(1) := (-cff(1) + sqrtdsc)/den2cff2;
      res(2) := (-cff(1) - sqrtdsc)/den2cff2;
    end if;
    return res;
  end Standard_Poles;

  function DoblDobl_Poles
              ( p : DoblDobl_Pade_Approximants.Pade )
              return DoblDobl_Complex_Vectors.Vector is

    deg : constant integer32
        := DoblDobl_Pade_Approximants.Denominator_Degree(p);
    res : DoblDobl_Complex_Vectors.Vector(1..deg);
    cff : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Pade_Approximants.Denominator_Coefficients(p);
    dsc,sqrtdsc,den2cff2 : DoblDobl_Complex_Numbers.Complex_Number;
    two : constant double_double := create(2.0);
    four : constant double_double := create(4.0);

  begin
    if deg = 1 then
      res(1) := -cff(0)/cff(1);
    elsif deg = 2 then
      dsc := cff(1)*cff(1) - four*cff(0)*cff(2);
      sqrtdsc := DoblDobl_Complex_Numbers_Polar.Root(dsc,2,1);
      den2cff2 := two*cff(2);
      res(1) := (-cff(1) + sqrtdsc)/den2cff2;
      res(2) := (-cff(1) - sqrtdsc)/den2cff2;
    end if;
    return res;
  end DoblDobl_Poles;

  function QuadDobl_Poles
              ( p : QuadDobl_Pade_Approximants.Pade )
              return QuadDobl_Complex_Vectors.Vector is

    deg : constant integer32
        := QuadDobl_Pade_Approximants.Denominator_Degree(p);
    res : QuadDobl_Complex_Vectors.Vector(1..deg);
    cff : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Pade_Approximants.Denominator_Coefficients(p);
    dsc,sqrtdsc,den2cff2 : QuadDobl_Complex_Numbers.Complex_Number;
    two : constant quad_double := create(2.0);
    four : constant quad_double := create(4.0);

  begin
    if deg = 1 then
      res(1) := -cff(0)/cff(1);
    elsif deg = 2 then
      dsc := cff(1)**2 - four*cff(0)*cff(2);
      sqrtdsc := QuadDobl_Complex_Numbers_Polar.Root(dsc,2,1);
      den2cff2 := two*cff(2);
      res(1) := (-cff(1) + sqrtdsc)/den2cff2;
      res(2) := (-cff(1) - sqrtdsc)/den2cff2;
    end if;
    return res;
  end QuadDobl_Poles;

  function Standard_Poles
              ( pv : Standard_Pade_Approximants.Pade_Vector )
              return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(pv'range);

  begin
    for i in pv'range loop
      declare
        poles : constant Standard_Complex_Vectors.Vector
              := Standard_Poles(pv(i));
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(poles);
      end;
    end loop;
    return res;
  end Standard_Poles;

  function DoblDobl_Poles
              ( pv : DoblDobl_Pade_Approximants.Pade_Vector )
              return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(pv'range);

  begin
    for i in pv'range loop
      declare
        poles : constant DoblDobl_Complex_Vectors.Vector
              := DoblDobl_Poles(pv(i));
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector'(poles);
      end;
    end loop;
    return res;
  end DoblDobl_Poles;

  function QuadDobl_Poles
              ( pv : QuadDobl_Pade_Approximants.Pade_Vector )
              return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(pv'range);

  begin
    for i in pv'range loop
      declare
        poles : constant QuadDobl_Complex_Vectors.Vector
              := QuadDobl_Poles(pv(i));
      begin
        res(i) := new QuadDobl_Complex_Vectors.Vector'(poles);
      end;
    end loop;
    return res;
  end QuadDobl_Poles;

end Homotopy_Pade_Approximants;
