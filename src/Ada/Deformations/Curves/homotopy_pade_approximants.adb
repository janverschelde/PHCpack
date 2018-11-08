with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with Double_Double_Numbers;             use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_Polar;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_Polar;
with Standard_Floating_Vectors;
with Black_Box_Univariate_Solvers;
with Homotopy_Series_Readers;           use Homotopy_Series_Readers;

package body Homotopy_Pade_Approximants is

  procedure Standard_Pade_Approximant
              ( sol : in Standard_Complex_Vectors.Vector;
                idx,nbequ,numdeg,dendeg : in integer32;
                nbiters : in natural32;
                srv,eva : out Standard_Complex_Series_Vectors.Vector;
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
                srv,eva : out DoblDobl_Complex_Series_Vectors.Vector;
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
                srv,eva : out QuadDobl_Complex_Series_Vectors.Vector;
                pv : out QuadDobl_Pade_Approximants.Pade_Vector;
                verbose : in boolean := false ) is

    nbterms : constant natural32 := natural32(numdeg+dendeg+1);

  begin
    QuadDobl_Series_Newton(sol,idx,nbequ,nbterms,nbiters,srv,eva);
    pv := QuadDobl_Pade_Approximants.Create(numdeg,dendeg,srv);
  end QuadDobl_Pade_Approximant;

  function Numerical_Degree
              ( p : Standard_Complex_Vectors.Vector;
                tol : double_float ) return integer32 is

    val : double_float;

  begin
    for i in reverse p'range loop
      val := Standard_Complex_Numbers.AbsVal(p(i));
      if val > tol
       then return i;
      end if;
    end loop;
    return -1;
  end Numerical_Degree;

  function Numerical_Degree
              ( p : DoblDobl_Complex_Vectors.Vector;
                tol : double_float ) return integer32 is

    val : double_double;

  begin
    for i in reverse p'range loop
      val := DoblDobl_Complex_Numbers.AbsVal(p(i));
      if val > tol
       then return i;
      end if;
    end loop;
    return -1;
  end Numerical_Degree;

  function Numerical_Degree
              ( p : QuadDobl_Complex_Vectors.Vector;
                tol : double_float ) return integer32 is

    val : quad_double;

  begin
    for i in reverse p'range loop
      val := QuadDobl_Complex_Numbers.AbsVal(p(i));
      if val > tol
       then return i;
      end if;
    end loop;
    return -1;
  end Numerical_Degree;

  function Standard_Poles
              ( p : Standard_Pade_Approximants.Pade )
              return Standard_Complex_Vectors.Vector is

    deg : constant integer32
        := Standard_Pade_Approximants.Denominator_Degree(p);
    minone : constant Standard_Complex_Numbers.Complex_Number
           := Standard_Complex_Numbers.Create(-1.0);
    res : Standard_Complex_Vectors.Vector(1..deg) := (1..deg => minone);
    cff : constant Standard_Complex_Vectors.Vector
        := Standard_Pade_Approximants.Denominator_Coefficients(p);
    numdeg : constant integer32 := Numerical_Degree(cff,1.0e-14);

  begin
    if numdeg > 0 then
      if numdeg = 1 then
        res(1) := -cff(0)/cff(1);
      elsif numdeg = 2 then
        declare
          dsc,sqrtdsc,den2cff2 : Standard_Complex_Numbers.Complex_Number;
          two : constant double_float := 2.0;
          four : constant double_float := 4.0;
        begin
          dsc := cff(1)**2 - four*cff(0)*cff(2);
          sqrtdsc := Standard_Complex_Numbers_Polar.Root(dsc,2,1);
          den2cff2 := two*cff(2);
          res(1) := (-cff(1) + sqrtdsc)/den2cff2;
          res(2) := (-cff(1) - sqrtdsc)/den2cff2;
        end;
      else
        declare
          err,rco,rsi : Standard_Floating_Vectors.Vector(1..numdeg);
          fail : boolean;
          use Black_Box_Univariate_Solvers;
        begin
          Standard_Compute_Roots(cff(0..numdeg),res(1..numdeg),
                                 err,rco,rsi,fail);
        end;
      end if;
    end if;
    return res;
  end Standard_Poles;

  function DoblDobl_Poles
              ( p : DoblDobl_Pade_Approximants.Pade )
              return DoblDobl_Complex_Vectors.Vector is

    deg : constant integer32
        := DoblDobl_Pade_Approximants.Denominator_Degree(p);
    dd_minone : constant double_double := create(-1.0);
    minone : constant DoblDobl_Complex_Numbers.Complex_Number
           := DoblDobl_Complex_Numbers.Create(dd_minone);
    res : DoblDobl_Complex_Vectors.Vector(1..deg) := (1..deg => minone);
    cff : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Pade_Approximants.Denominator_Coefficients(p);
    numdeg : constant integer32 := Numerical_Degree(cff,1.0e-30);

  begin
    if numdeg > 0 then
      if numdeg = 1 then
        res(1) := -cff(0)/cff(1);
      elsif numdeg = 2 then
        declare
          dsc,sqrtdsc,den2cff2 : DoblDobl_Complex_Numbers.Complex_Number;
          two : constant double_double := create(2.0);
          four : constant double_double := create(4.0);
        begin
          dsc := cff(1)*cff(1) - four*cff(0)*cff(2);
          sqrtdsc := DoblDobl_Complex_Numbers_Polar.Root(dsc,2,1);
          den2cff2 := two*cff(2);
          res(1) := (-cff(1) + sqrtdsc)/den2cff2;
          res(2) := (-cff(1) - sqrtdsc)/den2cff2;
        end;
      else
        declare
          err,rco,rsi : Standard_Floating_Vectors.Vector(1..numdeg);
          fail : boolean;
          use Black_Box_Univariate_Solvers;
        begin
          DoblDobl_Compute_Roots(cff(0..numdeg),res(1..numdeg),
                                 err,rco,rsi,fail);
        end;
      end if;
    end if;
    return res;
  end DoblDobl_Poles;

  function QuadDobl_Poles
              ( p : QuadDobl_Pade_Approximants.Pade )
              return QuadDobl_Complex_Vectors.Vector is

    deg : constant integer32
        := QuadDobl_Pade_Approximants.Denominator_Degree(p);
    qd_minone : constant quad_double := create(-1.0);
    minone : constant QuadDobl_Complex_Numbers.Complex_Number
           := QuadDobl_Complex_Numbers.Create(qd_minone);
    res : QuadDobl_Complex_Vectors.Vector(1..deg) := (1..deg => minone);
    cff : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Pade_Approximants.Denominator_Coefficients(p);
    numdeg : constant integer32 := Numerical_Degree(cff,1.0e-60);

  begin
    if numdeg > 0 then
      if numdeg = 1 then
        res(1) := -cff(0)/cff(1);
      elsif numdeg = 2 then
        declare
          dsc,sqrtdsc,den2cff2 : QuadDobl_Complex_Numbers.Complex_Number;
          two : constant quad_double := create(2.0);
          four : constant quad_double := create(4.0);
        begin
          dsc := cff(1)**2 - four*cff(0)*cff(2);
          sqrtdsc := QuadDobl_Complex_Numbers_Polar.Root(dsc,2,1);
          den2cff2 := two*cff(2);
          res(1) := (-cff(1) + sqrtdsc)/den2cff2;
          res(2) := (-cff(1) - sqrtdsc)/den2cff2;
        end;
      else
        declare
          err,rco,rsi : Standard_Floating_Vectors.Vector(1..numdeg);
          fail : boolean;
          use Black_Box_Univariate_Solvers;
        begin
          QuadDobl_Compute_Roots(cff(0..numdeg),res(1..numdeg),
                                 err,rco,rsi,fail);
        end;
      end if;
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

  procedure Smallest_Forward_Pole
              ( v : in Standard_Complex_Vectors.Vector;
                idx : out integer32; minval : out double_float ) is

    radval : double_float;

  begin
    idx := v'first;
    minval := Standard_Complex_Numbers_Polar.Radius(v(idx));
    if Standard_Complex_Numbers.REAL_Part(v(idx)) <= 0.0
     then minval := -minval;
    end if;
    for k in v'first+1..v'last loop
      radval := Standard_Complex_Numbers_Polar.Radius(v(k));
      if Standard_Complex_Numbers.REAL_PART(v(idx)) <= 0.0 then
        if minval < 0.0 then
          if radval < -minval
           then minval := -radval; idx := k; -- smallest backward pole
          end if;
        end if;
      elsif minval < 0.0 then
        minval := radval; idx := k; -- first forward pole
      elsif radval < minval then
        minval := radval; idx := k; -- found smaller forward pole
      end if;
    end loop;
  end Smallest_Forward_Pole;

  procedure Smallest_Forward_Pole
              ( v : in DoblDobl_Complex_Vectors.Vector;
                idx : out integer32; minval : out double_double ) is

    radval : double_double;

  begin
    idx := v'first;
    minval := DoblDobl_Complex_Numbers_Polar.Radius(v(idx));
    if DoblDobl_Complex_Numbers.REAL_Part(v(idx)) <= 0.0
     then minval := -minval;
    end if;
    for k in v'first+1..v'last loop
      radval := DoblDobl_Complex_Numbers_Polar.Radius(v(k));
      if DoblDobl_Complex_Numbers.REAL_PART(v(idx)) <= 0.0 then
        if minval < 0.0 then
          if radval < -minval
           then minval := -radval; idx := k; -- smallest backward pole
          end if;
        end if;
      elsif minval < 0.0 then
        minval := radval; idx := k; -- first forward pole
      elsif radval < minval then
        minval := radval; idx := k; -- found smaller forward pole
      end if;
    end loop;
  end Smallest_Forward_Pole;

  procedure Smallest_Forward_Pole
              ( v : in QuadDobl_Complex_Vectors.Vector;
                idx : out integer32; minval : out quad_double ) is

    radval : quad_double;

  begin
    idx := v'first;
    minval := QuadDobl_Complex_Numbers_Polar.Radius(v(idx));
    if QuadDobl_Complex_Numbers.REAL_Part(v(idx)) <= 0.0
     then minval := -minval;
    end if;
    for k in v'first+1..v'last loop
      radval := QuadDobl_Complex_Numbers_Polar.Radius(v(k));
      if QuadDobl_Complex_Numbers.REAL_PART(v(idx)) <= 0.0 then
        if minval < 0.0 then
          if radval < -minval
           then minval := -radval; idx := k; -- smallest backward pole
          end if;
        end if;
      elsif minval < 0.0 then
        minval := radval; idx := k; -- first forward pole
      elsif radval < minval then
        minval := radval; idx := k; -- found smaller forward pole
      end if;
    end loop;
  end Smallest_Forward_Pole;

  procedure Smallest_Forward_Pole
              ( v : in Standard_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out double_float ) is

    vkidx : integer32;
    radval : double_float;

  begin
    leadidx := v'first;
    Smallest_Forward_Pole(v(leadidx).all,idx,minval);
    for k in v'first+1..v'last loop
      Smallest_Forward_Pole(v(k).all,vkidx,radval);
      if minval < 0.0 then
        if radval < 0.0 then
          if -radval < -minval -- found smaller backward pole
           then minval := radval; leadidx := k; idx := vkidx;
          end if;
        else
          minval := radval; leadidx := k; idx := vkidx; -- first forward pole
        end if;
      end if;
      if radval > 0.0 then -- check if smaller forward pole found
        if radval < minval
         then minval := radval; leadidx := k; idx := vkidx;
        end if;
      end if;
    end loop;
  end Smallest_Forward_Pole;

  procedure Smallest_Forward_Pole
              ( v : in DoblDobl_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out double_double ) is

    vkidx : integer32;
    radval : double_double;

  begin
    leadidx := v'first;
    Smallest_Forward_Pole(v(leadidx).all,idx,minval);
    for k in v'first+1..v'last loop
      Smallest_Forward_Pole(v(k).all,vkidx,radval);
      if minval < 0.0 then
        if radval < 0.0 then
          if -radval < -minval -- found smaller backward pole
           then minval := radval; leadidx := k; idx := vkidx;
          end if;
        else
          minval := radval; leadidx := k; idx := vkidx; -- first forward pole
        end if;
      end if;
      if radval > 0.0 then -- check if smaller forward pole found
        if radval < minval
         then minval := radval; leadidx := k; idx := vkidx;
        end if;
      end if;
    end loop;
  end Smallest_Forward_Pole;

  procedure Smallest_Forward_Pole
              ( v : in QuadDobl_Complex_VecVecs.VecVec;
                leadidx,idx : out integer32; minval : out quad_double ) is

    vkidx : integer32;
    radval : quad_double;

  begin
    leadidx := v'first;
    Smallest_Forward_Pole(v(leadidx).all,idx,minval);
    for k in v'first+1..v'last loop
      Smallest_Forward_Pole(v(k).all,vkidx,radval);
      if minval < 0.0 then
        if radval < 0.0 then
          if -radval < -minval -- found smaller backward pole
           then minval := radval; leadidx := k; idx := vkidx;
          end if;
        else
          minval := radval; leadidx := k; idx := vkidx; -- first forward pole
        end if;
      end if;
      if radval > 0.0 then -- check if smaller forward pole found
        if radval < minval
         then minval := radval; leadidx := k; idx := vkidx;
        end if;
      end if;
    end loop;
  end Smallest_Forward_Pole;

end Homotopy_Pade_Approximants;
