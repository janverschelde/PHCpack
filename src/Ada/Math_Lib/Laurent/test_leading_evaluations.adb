with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Numbers_Polar;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Double_Leading_Evaluations;

package body Test_Leading_Evaluations is

  function Random_Monomial
             ( dim,low,upp : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..dim);

  begin
    for i in 1..dim loop
      res(i) := Standard_Random_Numbers.Random(low,upp);
    end loop;
    return res;
  end Random_Monomial;

  function Random_Polynomial
             ( nbr,dim,low,upp : integer32 )
             return Standard_Integer_VecVecs.VecVec is

    res : Standard_Integer_VecVecs.VecVec(1..nbr);

  begin
    for i in 1..nbr loop
      declare
        mon : constant Standard_Integer_Vectors.Vector(1..dim)
            := Random_Monomial(dim,low,upp);
      begin
        res(i) := new Standard_Integer_Vectors.Vector'(mon);
      end;
    end loop;
    return res;
  end Random_Polynomial;

  function Random_Leading_Powers
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..dim);
    rnd : double_float;

  begin
    for i in 1..dim loop
      rnd := Standard_Random_Numbers.Random;
      res(i) := abs(rnd);
    end loop;
    return res;
  end Random_Leading_Powers;

  procedure Test_Monomial_Derivative
              ( deg : in Standard_Integer_Vectors.Vector;
                pwr : in Standard_Floating_Vectors.Vector;
                cff : in Standard_Complex_Vectors.Vector;
                idx : in integer32; err : out double_float ) is

    lpr : constant double_float
        := Double_Leading_Evaluations.Leading_Power(deg,pwr,idx);
    lcf : constant Complex_Number
        := Double_Leading_Evaluations.Leading_Coefficient(deg,cff,idx);
   -- variables for the random point test
    t_rnd : constant double_float := abs(Standard_Random_Numbers.Random);
    t_pwr : constant double_float := t_rnd**lpr;
    t_val : constant Complex_Number := t_pwr*lcf;
    t_cff : Standard_Complex_Vectors.Vector(cff'range);
    t_tst,error : Complex_Number;

  begin
    put("-> testing derivative with respect to variable "); put(idx,1);
    put_line(" ...");
    put("leading power : "); put(lpr); new_line;
    put("leading coefficient : "); put(lcf); new_line;
    put("value at random point t = "); put(t_rnd); put_line(" :");
    put(t_val); new_line;
    for i in t_cff'range loop
      t_cff(i) := (t_rnd**pwr(i))*cff(i);
    end loop;
    t_tst := Double_Leading_Evaluations.Leading_Coefficient(deg,t_cff,idx);
    put_line("comparing to alternative evaluation :");
    put(t_tst); new_line;
    error := t_tst - t_val;
    err := AbsVal(error);
    put("error :"); put(err,3);
    if err < 1.0E-10
     then put_line(", okay");
     else put_line(" >= 1.0E-10, bug!");
    end if;
  end Test_Monomial_Derivative;

  procedure Test_Monomial ( dim : in integer32 ) is
               
    deg : constant Standard_Integer_Vectors.Vector(1..dim)
        := Random_Monomial(dim,-9,9);
    pwr : constant Standard_Floating_Vectors.Vector(1..dim)
        := Random_Leading_Powers(dim);
    cff : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    lpr : constant double_float
        := Double_Leading_Evaluations.Leading_Power(deg,pwr);
    lcf : constant Complex_Number
        := Double_Leading_Evaluations.Leading_Coefficient(deg,cff);
   -- variables for the random point test
    t_rnd : constant double_float := abs(Standard_Random_Numbers.Random);
    t_pwr : constant double_float := t_rnd**lpr;
    t_val : constant Complex_Number := t_pwr*lcf;
    t_cff : Standard_Complex_Vectors.Vector(1..dim);
    t_tst,error : Complex_Number;
    err,diferr : double_float;

  begin
    put("the degrees :"); put(deg); new_line;
    put_line("leading powers of the series :"); put_line(pwr);
    put_line("leading coefficients of the series :"); put_line(cff);
    put("power value : "); put(lpr); new_line;
    put("coeff value : "); put(lcf); new_line;
    put(t_rnd); put("**"); put(lpr); put(" = ");
    put(t_pwr); new_line;
    put("evaluating at t ="); put(t_rnd); put_line(" :");
    put(t_val); new_line;
    put("|value| :");
    put(Standard_Complex_Numbers_Polar.radius(t_val)); new_line;
    put("  angle : ");
    put(Standard_Complex_Numbers_Polar.angle(t_val)); new_line;
    for i in t_cff'range loop
      t_cff(i) := (t_rnd**pwr(i))*cff(i);
    end loop;
    put_line("leading coefficients multiplied by the point :");
    put_line(t_cff);
    t_tst := Double_Leading_Evaluations.Leading_Coefficient(deg,t_cff);
    put_line("value at the random point :");
    put(t_tst); new_line;
    put("|value| :");
    put(Standard_Complex_Numbers_Polar.radius(t_tst)); new_line;
    put("  angle : ");
    put(Standard_Complex_Numbers_Polar.angle(t_tst)); new_line;
    error := t_tst - t_val;
    err := AbsVal(error);
    put("error :"); put(err,3);
    if err < 1.0E-10
     then put_line(", okay");
     else put_line(" >= 1.0E-10, bug!");
    end if;
    for i in 1..dim loop
      if deg(i) /= 0 then
        Test_Monomial_Derivative(deg,pwr,cff,i,diferr);
        err := err + diferr;
      end if;
    end loop;
    put("sum of all errors :"); put(err,3);
    if err < 1.0E-10
     then put_line(", okay");
     else put_line(" >= 1.0E-10, bug!");
    end if;
  end Test_Monomial;

  procedure Test_Polynomial ( nbr,dim : in integer32 ) is

    deg : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Polynomial(nbr,dim,-9,9);
    pwr : constant Standard_Floating_Vectors.Vector(1..dim)
        := Random_Leading_Powers(dim);
    lpr : double_float;
    idx : integer32;

  begin
    put_line("The degrees : ");
    for i in deg'range loop
      put(deg(i).all); new_line;
    end loop;
    put_line("leading powers of the series :"); put_line(pwr);
    Double_Leading_Evaluations.Leading_Power(deg,pwr,lpr,idx,1);
    put("power value : "); put(lpr);
    put(" at index "); put(idx,1); new_line;
  end Test_Polynomial;

  procedure Main is

    nbr,dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables : "); get(dim);
    put("Give the number of monomials : "); get(nbr);
    if nbr = 1
     then Test_Monomial(dim);
     else Test_Polynomial(nbr,dim);
    end if;
  end Main;

end Test_Leading_Evaluations;
