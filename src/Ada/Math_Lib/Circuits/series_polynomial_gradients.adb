with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
-- with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Series;
with Standard_CSeries_Poly_Functions;
with DoblDobl_Complex_Series;
with DoblDobl_CSeries_Poly_Functions;
with QuadDobl_Complex_Series;
with QuadDobl_CSeries_Poly_Functions;

package body Series_Polynomial_Gradients is

  function Standard_Polynomial
             ( c : Standard_Speelpenning_Convolutions.Circuit )
             return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly
        := Standard_CSeries_Polynomials.Null_Poly;
    trm : Standard_CSeries_Polynomials.Term;

  begin
    for k in c.xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
      for i in 1..c.dim loop
        trm.dg(i) := natural32(c.xps(k)(i));
      end loop;
      trm.cf := Standard_Complex_Series.Create(c.cff(k).all);
      Standard_CSeries_Polynomials.Add(res,trm);
      Standard_CSeries_Polynomials.Clear(trm);
    end loop;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
    trm.cf := Standard_Complex_Series.Create(c.cst.all);
    Standard_CSeries_Polynomials.Add(res,trm);
    Standard_CSeries_Polynomials.Clear(trm);
    return res;
  end Standard_Polynomial;

  function DoblDobl_Polynomial
             ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
             return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly
        := DoblDobl_CSeries_Polynomials.Null_Poly;
    trm : DoblDobl_CSeries_Polynomials.Term;

  begin
    for k in c.xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
      for i in 1..c.dim loop
        trm.dg(i) := natural32(c.xps(k)(i));
      end loop;
      trm.cf := DoblDobl_Complex_Series.Create(c.cff(k).all);
      DoblDobl_CSeries_Polynomials.Add(res,trm);
      DoblDobl_CSeries_Polynomials.Clear(trm);
    end loop;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
    trm.cf := DoblDobl_Complex_Series.Create(c.cst.all);
    DoblDobl_CSeries_Polynomials.Add(res,trm);
    DoblDobl_CSeries_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Polynomial;

  function QuadDobl_Polynomial
             ( c : QuadDobl_Speelpenning_Convolutions.Circuit )
             return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly
        := QuadDobl_CSeries_Polynomials.Null_Poly;
    trm : QuadDobl_CSeries_Polynomials.Term;

  begin
    for k in c.xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
      for i in 1..c.dim loop
        trm.dg(i) := natural32(c.xps(k)(i));
      end loop;
      trm.cf := QuadDobl_Complex_Series.Create(c.cff(k).all);
      QuadDobl_CSeries_Polynomials.Add(res,trm);
      QuadDobl_CSeries_Polynomials.Clear(trm);
    end loop;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..c.dim => 0);
    trm.cf := QuadDobl_Complex_Series.Create(c.cst.all);
    QuadDobl_CSeries_Polynomials.Add(res,trm);
    QuadDobl_CSeries_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Polynomial;

  function Standard_System
             ( c : Standard_Speelpenning_Convolutions.Circuits )
             return Standard_CSeries_Poly_Systems.Poly_Sys is

    res : Standard_CSeries_Poly_Systems.Poly_Sys(c'range);

  begin
    for i in c'range loop
      res(i) := Standard_Polynomial(c(i).all);
    end loop;
    return res;
  end Standard_System;

  function DoblDobl_System
             ( c : DoblDobl_Speelpenning_Convolutions.Circuits )
             return DoblDobl_CSeries_Poly_Systems.Poly_Sys is

    res : DoblDobl_CSeries_Poly_Systems.Poly_Sys(c'range);

  begin
    for i in c'range loop
      res(i) := DoblDobl_Polynomial(c(i).all);
    end loop;
    return res;
  end DoblDobl_System;

  function QuadDobl_System
             ( c : QuadDobl_Speelpenning_Convolutions.Circuits )
             return QuadDobl_CSeries_Poly_Systems.Poly_Sys is

    res : QuadDobl_CSeries_Poly_Systems.Poly_Sys(c'range);

  begin
    for i in c'range loop
      res(i) := QuadDobl_Polynomial(c(i).all);
    end loop;
    return res;
  end QuadDobl_System;

  function Standard_Product
             ( dim,deg : in integer32 )
             return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly;
    trm : Standard_CSeries_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := Standard_Complex_Series.Create(1,deg);
    res := Standard_CSeries_Polynomials.Create(trm);
    Standard_CSeries_Polynomials.Clear(trm);
    return res;
  end Standard_Product;

  function DoblDobl_Product
             ( dim,deg : in integer32 )
             return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly;
    trm : DoblDobl_CSeries_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := DoblDobl_Complex_Series.Create(1,deg);
    res := DoblDobl_CSeries_Polynomials.Create(trm);
    DoblDobl_CSeries_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Product;

  function QuadDobl_Product
             ( dim,deg : in integer32 )
             return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly;
    trm : QuadDobl_CSeries_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 1);
    trm.cf := QuadDobl_Complex_Series.Create(1,deg);
    res := QuadDobl_CSeries_Polynomials.Create(trm);
    QuadDobl_CSeries_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Product;

  function Standard_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly;
    trm : Standard_CSeries_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := Standard_Complex_Series.Create(1,deg);
    res := Standard_CSeries_Polynomials.Create(trm);
    Standard_CSeries_Polynomials.Clear(trm);
    return res;
  end Standard_Product;

  function DoblDobl_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly;
    trm : DoblDobl_CSeries_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := DoblDobl_Complex_Series.Create(1,deg);
    res := DoblDobl_CSeries_Polynomials.Create(trm);
    DoblDobl_CSeries_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Product;

  function QuadDobl_Product
             ( deg : in integer32;
               xp : Standard_Integer_Vectors.Vector )
             return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly;
    trm : QuadDobl_CSeries_Polynomials.Term;

  begin
    trm.dg := new Standard_Natural_Vectors.Vector(1..xp'last);
    for i in xp'range loop
      trm.dg(i) := natural32(xp(i));
    end loop;
    trm.cf := QuadDobl_Complex_Series.Create(1,deg);
    res := QuadDobl_CSeries_Polynomials.Create(trm);
    QuadDobl_CSeries_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Product;

  -- function Is_Zero ( s : Standard_Complex_Series.Link_to_Series )
  --                  return boolean is

  -- DESCRIPTION :
  --   Returns true if the series is zero.

  --  use Standard_Complex_Numbers;
  --  use Standard_Complex_Series;

  --  zero : constant Complex_Number := Create(0.0);

  -- begin
  --   if s = null then
  --     return true;
  --   else
  --     for i in s.cff'range loop
  --       if s.cff(i) /= zero
  --        then return false;
  --       end if;
  --     end loop;
  --     return true;
  --   end if;
  -- end Is_Zero;

  function Standard_Polynomial
             ( dim,deg : in integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               isidx : boolean := true )
             return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly
        := Standard_CSeries_Polynomials.Null_Poly;
    trm : Standard_CSeries_Polynomials.Term;
   -- cff : Standard_Complex_Series.Link_to_Series;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      if isidx then
        for i in xps(k)'range loop
          trm.dg(xps(k)(i)) := 1;
        end loop;
      else
        for i in 1..dim loop
          trm.dg(i) := natural32(xps(k)(i));
        end loop;
      end if;
     -- cff := Standard_CSeries_Polynomials.Coeff(res,trm.dg);
     -- if Is_Zero(cff) then
      trm.cf := Standard_Complex_Series.Create(1,deg);
      Standard_CSeries_Polynomials.Add(res,trm);
     -- end if;
      Standard_CSeries_Polynomials.Clear(trm);
    end loop;
    return res;
  end Standard_Polynomial;

  function Standard_Polynomial
             ( dim : in integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : Standard_Complex_Series_Vectors.Vector;
               isxidx : boolean := true )
             return Standard_CSeries_Polynomials.Poly is

    res : Standard_CSeries_Polynomials.Poly
        := Standard_CSeries_Polynomials.Null_Poly;
    trm : Standard_CSeries_Polynomials.Term;

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
      Standard_Complex_Series.Copy(cff(k),trm.cf);
      Standard_CSeries_Polynomials.Add(res,trm);
      Standard_CSeries_Polynomials.Clear(trm);
    end loop;
    return res;
  end Standard_Polynomial;

  function DoblDobl_Polynomial
             ( dim,deg : in integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly
        := DoblDobl_CSeries_Polynomials.Null_Poly;
    trm : DoblDobl_CSeries_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      for i in xps(k)'range loop
        trm.dg(xps(k)(i)) := 1;
      end loop;
      trm.cf := DoblDobl_Complex_Series.Create(1,deg);
      DoblDobl_CSeries_Polynomials.Add(res,trm);
      DoblDobl_CSeries_Polynomials.Clear(trm);
    end loop;
    return res;
  end DoblDobl_Polynomial;

  function DoblDobl_Polynomial
             ( dim : in integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : DoblDobl_Complex_Series_Vectors.Vector;
               isxidx : boolean := true )
             return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly
        := DoblDobl_CSeries_Polynomials.Null_Poly;
    trm : DoblDobl_CSeries_Polynomials.Term;

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
      DoblDobl_Complex_Series.Copy(cff(k),trm.cf);
      DoblDobl_CSeries_Polynomials.Add(res,trm);
      DoblDobl_CSeries_Polynomials.Clear(trm);
    end loop;
    return res;
  end DoblDobl_Polynomial;

  function QuadDobl_Polynomial
             ( dim,deg : in integer32;
               xps : Standard_Integer_VecVecs.VecVec )
             return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly
        := QuadDobl_CSeries_Polynomials.Null_Poly;
    trm : QuadDobl_CSeries_Polynomials.Term;

  begin
    for k in xps'range loop
      trm.dg := new Standard_Natural_Vectors.Vector'(1..dim => 0);
      for i in xps(k)'range loop
        trm.dg(xps(k)(i)) := 1;
      end loop;
      trm.cf := QuadDobl_Complex_Series.Create(1,deg);
      QuadDobl_CSeries_Polynomials.Add(res,trm);
      QuadDobl_CSeries_Polynomials.Clear(trm);
    end loop;
    return res;
  end QuadDobl_Polynomial;

  function QuadDobl_Polynomial
             ( dim : in integer32;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : QuadDobl_Complex_Series_Vectors.Vector;
               isxidx : boolean := true )
             return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly
        := QuadDobl_CSeries_Polynomials.Null_Poly;
    trm : QuadDobl_CSeries_Polynomials.Term;

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
      QuadDobl_Complex_Series.Copy(cff(k),trm.cf);
      QuadDobl_CSeries_Polynomials.Add(res,trm);
      QuadDobl_CSeries_Polynomials.Clear(trm);
    end loop;
    return res;
  end QuadDobl_Polynomial;

  function Standard_Gradient
             ( p : Standard_CSeries_Polynomials.Poly;
               x : Standard_Complex_Series_Vectors.Vector )
             return Standard_Complex_Series_Vectors.Vector is

    res : Standard_Complex_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : Standard_CSeries_Polynomials.Poly
            := Standard_CSeries_Polynomials.Diff(p,k);
        val : constant Standard_Complex_Series.Link_to_Series
            := Standard_CSeries_Poly_Functions.Eval(dpk,x);
      begin
        Standard_CSeries_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end Standard_Gradient;

  function DoblDobl_Gradient
             ( p : DoblDobl_CSeries_Polynomials.Poly;
               x : DoblDobl_Complex_Series_Vectors.Vector )
             return DoblDobl_Complex_Series_Vectors.Vector is

    res : DoblDobl_Complex_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : DoblDobl_CSeries_Polynomials.Poly
            := DoblDobl_CSeries_Polynomials.Diff(p,k);
        val : constant DoblDobl_Complex_Series.Link_to_Series
            := DoblDobl_CSeries_Poly_Functions.Eval(dpk,x);
      begin
        DoblDobl_CSeries_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end DoblDobl_Gradient;

  function QuadDobl_Gradient
             ( p : QuadDobl_CSeries_Polynomials.Poly;
               x : QuadDobl_Complex_Series_Vectors.Vector )
             return QuadDobl_Complex_Series_Vectors.Vector is

    res : QuadDobl_Complex_Series_Vectors.Vector(x'range);

  begin
    for k in x'range loop
      declare
        dpk : QuadDobl_CSeries_Polynomials.Poly
            := QuadDobl_CSeries_Polynomials.Diff(p,k);
        val : constant QuadDobl_Complex_Series.Link_to_Series
            := QuadDobl_CSeries_Poly_Functions.Eval(dpk,x);
      begin
        QuadDobl_CSeries_Polynomials.Clear(dpk);
        res(k) := val;
      end;
    end loop;
    return res;
  end QuadDobl_Gradient;

end Series_Polynomial_Gradients;
