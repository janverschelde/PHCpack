with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Dense_Series;
with Standard_Series_Poly_Functions;
with DoblDobl_Dense_Series;
with DoblDobl_Series_Poly_Functions;
with QuadDobl_Dense_Series;
with QuadDobl_Series_Poly_Functions;

package body Series_Polynomial_Gradients is

  function Standard_Polynomial
             ( c : Standard_Speelpenning_Convolutions.Convolution_Circuit )
             return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    for k in c.xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
      for i in 1..c.dim loop
        trm.dg(i) := natural32(c.xps(k)(i));
      end loop;
      trm.cf := Standard_Dense_Series.Create(c.cff(k).all);
      Standard_Series_Polynomials.Add(res,trm);
      Standard_Series_Polynomials.Clear(trm);
    end loop;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
    trm.cf := Standard_Dense_Series.Create(c.cst.all);
    Standard_Series_Polynomials.Add(res,trm);
    Standard_Series_Polynomials.Clear(trm);
    return res;
  end Standard_Polynomial;

  function DoblDobl_Polynomial
             ( c : DoblDobl_Speelpenning_Convolutions.Convolution_Circuit )
             return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Series_Polynomials.Null_Poly;
    trm : DoblDobl_Series_Polynomials.Term;

  begin
    for k in c.xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
      for i in 1..c.dim loop
        trm.dg(i) := natural32(c.xps(k)(i));
      end loop;
      trm.cf := DoblDobl_Dense_Series.Create(c.cff(k).all);
      DoblDobl_Series_Polynomials.Add(res,trm);
      DoblDobl_Series_Polynomials.Clear(trm);
    end loop;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
    trm.cf := DoblDobl_Dense_Series.Create(c.cst.all);
    DoblDobl_Series_Polynomials.Add(res,trm);
    DoblDobl_Series_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Polynomial;

  function QuadDobl_Polynomial
             ( c : QuadDobl_Speelpenning_Convolutions.Convolution_Circuit )
             return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Series_Polynomials.Null_Poly;
    trm : QuadDobl_Series_Polynomials.Term;

  begin
    for k in c.xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
      for i in 1..c.dim loop
        trm.dg(i) := natural32(c.xps(k)(i));
      end loop;
      trm.cf := QuadDobl_Dense_Series.Create(c.cff(k).all);
      QuadDobl_Series_Polynomials.Add(res,trm);
      QuadDobl_Series_Polynomials.Clear(trm);
    end loop;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
    trm.cf := QuadDobl_Dense_Series.Create(c.cst.all);
    QuadDobl_Series_Polynomials.Add(res,trm);
    QuadDobl_Series_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Polynomial;

  function Standard_System
             ( c : Standard_Speelpenning_Convolutions.Convolution_Circuits )
             return Standard_Series_Poly_Systems.Poly_Sys is

    res : Standard_Series_Poly_Systems.Poly_Sys(c'range);

  begin
    for i in c'range loop
      res(i) := Standard_Polynomial(c(i).all);
    end loop;
    return res;
  end Standard_System;

  function DoblDobl_System
             ( c : DoblDobl_Speelpenning_Convolutions.Convolution_Circuits )
             return DoblDobl_Series_Poly_Systems.Poly_Sys is

    res : DoblDobl_Series_Poly_Systems.Poly_Sys(c'range);

  begin
    for i in c'range loop
      res(i) := DoblDobl_Polynomial(c(i).all);
    end loop;
    return res;
  end DoblDobl_System;

  function QuadDobl_System
             ( c : QuadDobl_Speelpenning_Convolutions.Convolution_Circuits )
             return QuadDobl_Series_Poly_Systems.Poly_Sys is

    res : QuadDobl_Series_Poly_Systems.Poly_Sys(c'range);

  begin
    for i in c'range loop
      res(i) := QuadDobl_Polynomial(c(i).all);
    end loop;
    return res;
  end QuadDobl_System;

  function Standard_Product
             ( dim,deg : in integer32 )
             return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := Standard_Dense_Series.Create(1.0,deg);
    res := Standard_Series_Polynomials.Create(trm);
    Standard_Series_Polynomials.Clear(trm);
    return res;
  end Standard_Product;

  function DoblDobl_Product
             ( dim,deg : in integer32 )
             return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly;
    trm : DoblDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := DoblDobl_Dense_Series.Create(1.0,deg);
    res := DoblDobl_Series_Polynomials.Create(trm);
    DoblDobl_Series_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Product;

  function QuadDobl_Product
             ( dim,deg : in integer32 )
             return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly;
    trm : QuadDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := QuadDobl_Dense_Series.Create(1.0,deg);
    res := QuadDobl_Series_Polynomials.Create(trm);
    QuadDobl_Series_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Product;

  function Standard_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := Standard_Dense_Series.Create(1.0,deg);
    res := Standard_Series_Polynomials.Create(trm);
    Standard_Series_Polynomials.Clear(trm);
    return res;
  end Standard_Product;

  function DoblDobl_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly;
    trm : DoblDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := DoblDobl_Dense_Series.Create(1.0,deg);
    res := DoblDobl_Series_Polynomials.Create(trm);
    DoblDobl_Series_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Product;

  function QuadDobl_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly;
    trm : QuadDobl_Series_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := QuadDobl_Dense_Series.Create(1.0,deg);
    res := QuadDobl_Series_Polynomials.Create(trm);
    QuadDobl_Series_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Product;

  function Standard_Polynomial
             ( dim,deg : in integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      for i in xps(k)'range loop
        trm.dg(xps(k)(i)) := 1;
      end loop;
      trm.cf := Standard_Dense_Series.Create(1.0,deg);
      Standard_Series_Polynomials.Add(res,trm);
      Standard_Series_Polynomials.Clear(trm);
    end loop;
    return res;
  end Standard_Polynomial;

  function Standard_Polynomial
             ( dim : in integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : Standard_Dense_Series_Vectors.Vector;
               isxidx : boolean := true )
             return Standard_Series_Polynomials.Poly is

    res : Standard_Series_Polynomials.Poly
        := Standard_Series_Polynomials.Null_Poly;
    trm : Standard_Series_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      if isxidx then
        for i in xps(k)'range loop
          trm.dg(xps(k)(i)) := 1;
        end loop;
      else
        for i in 1..dim loop
          trm.dg(i) := natural32(xps(k)(i));
        end loop;
      end if;
      Standard_Dense_Series.Copy(cff(k),trm.cf);
      Standard_Series_Polynomials.Add(res,trm);
      Standard_Series_Polynomials.Clear(trm);
    end loop;
    return res;
  end Standard_Polynomial;

  function DoblDobl_Polynomial
             ( dim,deg : in integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Series_Polynomials.Null_Poly;
    trm : DoblDobl_Series_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      for i in xps(k)'range loop
        trm.dg(xps(k)(i)) := 1;
      end loop;
      trm.cf := DoblDobl_Dense_Series.Create(1.0,deg);
      DoblDobl_Series_Polynomials.Add(res,trm);
      DoblDobl_Series_Polynomials.Clear(trm);
    end loop;
    return res;
  end DoblDobl_Polynomial;

  function DoblDobl_Polynomial
             ( dim : in integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : DoblDobl_Dense_Series_Vectors.Vector;
               isxidx : boolean := true )
             return DoblDobl_Series_Polynomials.Poly is

    res : DoblDobl_Series_Polynomials.Poly
        := DoblDobl_Series_Polynomials.Null_Poly;
    trm : DoblDobl_Series_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      if isxidx then
        for i in xps(k)'range loop
          trm.dg(xps(k)(i)) := 1;
        end loop;
      else
        for i in 1..dim loop
          trm.dg(i) := natural32(xps(k)(i));
        end loop;
      end if;
      DoblDobl_Dense_Series.Copy(cff(k),trm.cf);
      DoblDobl_Series_Polynomials.Add(res,trm);
      DoblDobl_Series_Polynomials.Clear(trm);
    end loop;
    return res;
  end DoblDobl_Polynomial;

  function QuadDobl_Polynomial
             ( dim,deg : in integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Series_Polynomials.Null_Poly;
    trm : QuadDobl_Series_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      for i in xps(k)'range loop
        trm.dg(xps(k)(i)) := 1;
      end loop;
      trm.cf := QuadDobl_Dense_Series.Create(1.0,deg);
      QuadDobl_Series_Polynomials.Add(res,trm);
      QuadDobl_Series_Polynomials.Clear(trm);
    end loop;
    return res;
  end QuadDobl_Polynomial;

  function QuadDobl_Polynomial
             ( dim : in integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : QuadDobl_Dense_Series_Vectors.Vector;
               isxidx : boolean := true )
             return QuadDobl_Series_Polynomials.Poly is

    res : QuadDobl_Series_Polynomials.Poly
        := QuadDobl_Series_Polynomials.Null_Poly;
    trm : QuadDobl_Series_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      if isxidx then
        for i in xps(k)'range loop
          trm.dg(xps(k)(i)) := 1;
        end loop;
      else
        for i in 1..dim loop
          trm.dg(i) := natural32(xps(k)(i));
        end loop;
      end if;
      QuadDobl_Dense_Series.Copy(cff(k),trm.cf);
      QuadDobl_Series_Polynomials.Add(res,trm);
      QuadDobl_Series_Polynomials.Clear(trm);
    end loop;
    return res;
  end QuadDobl_Polynomial;

  function Standard_Gradient
             ( p : Standard_Series_Polynomials.Poly;
               x : Standard_Dense_Series_Vectors.Vector )
             return Standard_Dense_Series_Vectors.Vector is

    res : Standard_Dense_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : Standard_Series_Polynomials.Poly
            := Standard_Series_Polynomials.Diff(p,k);
        val : constant Standard_Dense_Series.Series
            := Standard_Series_Poly_Functions.Eval(dpk,x);
      begin
        Standard_Series_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end Standard_Gradient;

  function DoblDobl_Gradient
             ( p : DoblDobl_Series_Polynomials.Poly;
               x : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Dense_Series_Vectors.Vector is

    res : DoblDobl_Dense_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : DoblDobl_Series_Polynomials.Poly
            := DoblDobl_Series_Polynomials.Diff(p,k);
        val : constant DoblDobl_Dense_Series.Series
            := DoblDobl_Series_Poly_Functions.Eval(dpk,x);
      begin
        DoblDobl_Series_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end DoblDobl_Gradient;

  function QuadDobl_Gradient
             ( p : QuadDobl_Series_Polynomials.Poly;
               x : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Dense_Series_Vectors.Vector is

    res : QuadDobl_Dense_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : QuadDobl_Series_Polynomials.Poly
            := QuadDobl_Series_Polynomials.Diff(p,k);
        val : constant QuadDobl_Dense_Series.Series
            := QuadDobl_Series_Poly_Functions.Eval(dpk,x);
      begin
        QuadDobl_Series_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end QuadDobl_Gradient;

  function Standard_Series_Coefficients
             ( s : Standard_Dense_Series_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant Standard_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Standard_Series_Coefficients;

  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Dense_Series_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant DoblDobl_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end DoblDobl_Series_Coefficients;

  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Dense_Series_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end QuadDobl_Series_Coefficients;

end Series_Polynomial_Gradients;
