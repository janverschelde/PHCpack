with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Functions;
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

end Affine_Transformations;
