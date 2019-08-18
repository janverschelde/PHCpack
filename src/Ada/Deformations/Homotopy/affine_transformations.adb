with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;

package body Affine_Transformations is

  function Make_Affine
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

   res : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last-1);
   one : constant Standard_Complex_Numbers.Complex_Number
       := Standard_Complex_Numbers.Create(1.0);

  begin
    for i in res'range loop
      res(i) := Standard_Complex_Poly_Functions.Eval(p(i),one,p'last);
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

   res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-1);
   dd_one : constant double_double := create(1.0);
   one : constant DoblDobl_Complex_Numbers.Complex_Number
       := DoblDobl_Complex_Numbers.Create(dd_one);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Poly_Functions.Eval(p(i),one,p'last);
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

   res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-1);
   qd_one : constant quad_double := create(1.0);
   one : constant QuadDobl_Complex_Numbers.Complex_Number
       := QuadDobl_Complex_Numbers.Create(qd_one);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Poly_Functions.Eval(p(i),one,p'last);
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

   im : constant integer32 := integer32(m);
   res : Standard_Complex_Poly_Systems.Poly_Sys(p'first..p'last-im);
   one : constant Standard_Complex_Numbers.Complex_Number
       := Standard_Complex_Numbers.Create(1.0);
   eva : Standard_Complex_Polynomials.Poly;
   idx : integer32;

  begin
    for i in res'range loop
      idx := p'last;
      for k in 1..m loop
        if k = 1
         then eva := Standard_Complex_Poly_Functions.Eval(p(i),one,idx);
         else eva := Standard_Complex_Poly_Functions.Eval(res(i),one,idx);
        end if;
        Standard_Complex_Polynomials.Copy(eva,res(i));
        Standard_Complex_Polynomials.Clear(eva);
        idx := idx - 1;
      end loop;
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

   im : constant integer32 := integer32(m);
   res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-im);
   dd_one : constant double_double := create(1.0);
   one : constant DoblDobl_Complex_Numbers.Complex_Number
       := DoblDobl_Complex_Numbers.Create(dd_one);
   eva : DoblDobl_Complex_Polynomials.Poly;
   idx : integer32;

  begin
    for i in res'range loop
      idx := p'last;
      for k in 1..m loop
        if k = 1
         then eva := DoblDobl_Complex_Poly_Functions.Eval(p(i),one,idx);
         else eva := DoblDobl_Complex_Poly_Functions.Eval(res(i),one,idx);
        end if;
        DoblDobl_Complex_Polynomials.Copy(eva,res(i));
        DoblDobl_Complex_Polynomials.Clear(eva);
        idx := idx - 1;
      end loop;
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

   im : constant integer32 := integer32(m);
   res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'first..p'last-im);
   qd_one : constant quad_double := create(1.0);
   one : constant QuadDobl_Complex_Numbers.Complex_Number
       := QuadDobl_Complex_Numbers.Create(qd_one);
   eva : QuadDobl_Complex_Polynomials.Poly;
   idx : integer32;

  begin
    for i in res'range loop
      idx := p'last;
      for k in 1..m loop
        if k = 1
         then eva := QuadDobl_Complex_Poly_Functions.Eval(p(i),one,idx);
         else eva := QuadDobl_Complex_Poly_Functions.Eval(res(i),one,idx);
        end if;
        QuadDobl_Complex_Polynomials.Copy(eva,res(i));
        QuadDobl_Complex_Polynomials.Clear(eva);
        idx := idx - 1;
      end loop;
    end loop;
    return res;
  end Make_Affine;

end Affine_Transformations;
