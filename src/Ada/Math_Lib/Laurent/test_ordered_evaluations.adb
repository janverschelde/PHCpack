with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_IO;       use Standard_Integer_Vectors_IO;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;
with Random_Laurent_Homotopy;
with Double_Real_Powered_Series;
with Double_Ordered_Evaluations;

package body Test_Ordered_Evaluations is

  procedure Test_First_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Tests the first derivative first order evaluation of a polynomial 
  --   represented by (pcf, pct, pdg), evaluated at a series
  --   with coefficients in cff and powers in pwr.

    use Double_Ordered_Evaluations;

    dim : constant integer32 := cff'last;
    nbr : constant integer32 := pdg'last;
    size : constant integer32 := Size_Evaluation(dim,1,1,nbr);
    ycf1,ycf2 : Standard_Complex_Vectors.Vector(1..size);
    ydg1,ydg2 : Standard_Floating_Vectors.Vector(1..size);
    err,sumerr : double_float;
    cf0,cf1 : Standard_Complex_Vectors.Vector(1..dim);
    pw1 : Standard_Floating_Vectors.Vector(1..dim);

  begin
    First_Derivative_First_Order(pcf,pct,pdg,cff,pwr,ycf1,ydg1,1);
    for i in cf0'range loop
      cf0(i) := cff(i)(cff(i)'first);
      cf1(i) := cff(i)(cff(i)'first+1);
      pw1(i) := pwr(i)(pwr(i)'first);
    end loop;
    First_Order_Evaluation(pcf,pct,pdg,cf0,cf1,pw1,1,ycf2,ydg2,1);
    put_line("the first derivative first order evaluation :");
    sumerr := 0.0;
    for i in ycf1'range loop
      put(ycf1(i)); put("  t^"); put(ydg1(i)); new_line;
      put(ycf2(i)); put("  t^"); put(ydg2(i)); new_line;
      err := AbsVal(ycf1(i) - ycf2(i)) + abs(ydg1(i) - ydg2(i));
      put("error :"); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("sum of errors :"); put(sumerr,3); new_line;
  end Test_First_Derivative_First_Order;

  procedure Test_Second_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Tests the second derivative first order evaluation of a polynomial
  --   represented by (pcf, pct, pdg), evaluated at a series
  --   with coefficients in cff and powers in pwr.

    use Double_Ordered_Evaluations;

    dim : constant integer32 := cff'last;
    nbr : constant integer32 := pdg'last;
    size : constant integer32 := Size_Evaluation(dim,2,1,nbr);
    ycf1,ycf2 : Standard_Complex_Vectors.Vector(1..size);
    ydg1,ydg2 : Standard_Floating_Vectors.Vector(1..size);
    err,sumerr : double_float;
    cf0,cf1 : Standard_Complex_Vectors.Vector(1..dim);
    pw1 : Standard_Floating_Vectors.Vector(1..dim);

  begin
    Second_Derivative_First_Order(pcf,pct,pdg,cff,pwr,ycf1,ydg1,1);
    for i in cf0'range loop
      cf0(i) := cff(i)(cff(i)'first);
      cf1(i) := cff(i)(cff(i)'first+1);
      pw1(i) := pwr(i)(pwr(i)'first);
    end loop;
    First_Order_Evaluation(pcf,pct,pdg,cf0,cf1,pw1,2,ycf2,ydg2,1);
    put_line("the second derivative first order evaluation :");
    sumerr := 0.0;
    for i in ycf1'range loop
      put(ycf1(i)); put("  t^"); put(ydg1(i)); new_line;
      put(ycf2(i)); put("  t^"); put(ydg2(i)); new_line;
      err := AbsVal(ycf1(i) - ycf2(i)) + abs(ydg1(i) - ydg2(i));
      put("error :"); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("sum of errors :"); put(sumerr,3); new_line;
  end Test_Second_Derivative_First_Order;

  procedure Test_Third_Derivative_First_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Tests the third derivative first order evaluation of a polynomial
  --   represented by (pdg, pcf, pct), evaluated at a series
  --   with coefficients in cff and powers in pwr.

    use Double_Ordered_Evaluations;

    dim : constant integer32 := cff'last;
    nbr : constant integer32 := pdg'last;
    size : constant integer32 := Size_Evaluation(dim,3,1,nbr);
    ycf1,ycf2 : Standard_Complex_Vectors.Vector(1..size);
    ydg1,ydg2 : Standard_Floating_Vectors.Vector(1..size);
    err,sumerr : double_float;
    cf0,cf1 : Standard_Complex_Vectors.Vector(1..dim);
    pw1 : Standard_Floating_Vectors.Vector(1..dim);

  begin
    Third_Derivative_First_Order(pcf,pct,pdg,cff,pwr,ycf1,ydg1,1);
    for i in cf0'range loop
      cf0(i) := cff(i)(cff(i)'first);
      cf1(i) := cff(i)(cff(i)'first+1);
      pw1(i) := pwr(i)(pwr(i)'first);
    end loop;
    First_Order_Evaluation(pcf,pct,pdg,cf0,cf1,pw1,3,ycf2,ydg2,1);
    put_line("the third derivative first order evaluation :");
    sumerr := 0.0;
    for i in ycf1'range loop
      put(ycf1(i)); put("  t^"); put(ydg1(i)); new_line;
      put(ycf2(i)); put("  t^"); put(ydg2(i)); new_line;
      err := AbsVal(ycf1(i) - ycf2(i)) + abs(ydg1(i) - ydg2(i));
      put("error :"); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("sum of errors :"); put(sumerr,3); new_line;
  end Test_Third_Derivative_First_Order;

  procedure Test_First_Order_Evaluation
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec;
                difmax : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the first order evaluation of a polynomial
  --   represented by (pdg, pcf, pct), evaluated at a series
  --   with coefficients in cff and powers in pwr,
  --   using all derivatives up to difmax.

    use Double_Ordered_Evaluations;

    dim : constant integer32 := cff'last;
    nbr : constant integer32 := pdg'last;
    size : constant integer32 := Size_Evaluation(dim,difmax,1,nbr);
    ycf : Standard_Complex_Vectors.Vector(1..size);
    ydg : Standard_Floating_Vectors.Vector(1..size);
    cf0,cf1 : Standard_Complex_Vectors.Vector(1..dim);
    pw1 : Standard_Floating_Vectors.Vector(1..dim);

  begin
    for i in cf0'range loop
      cf0(i) := cff(i)(cff(i)'first);
      cf1(i) := cff(i)(cff(i)'first+1);
      pw1(i) := pwr(i)(pwr(i)'first);
    end loop;
    First_Order_Evaluation(pcf,pct,pdg,cf0,cf1,pw1,difmax,ycf,ydg,1);
    put("derivative "); put(difmax,1);
    put_line(" first order evaluation :");
    for i in ycf'range loop
      put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
    end loop;
  end Test_First_Order_Evaluation;

  procedure Test_Second_Derivative_Second_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Tests the 2nd derivative 2nd order evaluation of a polynomial
  --   represented by (pcf, pct, pdg), evaluated at a series
  --   with coefficients in cff and powers in pwr.

    use Double_Ordered_Evaluations;

    dim : constant integer32 := cff'last;
    nbr : constant integer32 := pdg'last;
    size : constant integer32 := Size_Evaluation(dim,2,2,nbr);
    ycf : Standard_Complex_Vectors.Vector(1..size);
    ydg : Standard_Floating_Vectors.Vector(1..size);
    cf0,cf1,cf2 : Standard_Complex_Vectors.Vector(1..dim);
    pw1,pw2 : Standard_Floating_Vectors.Vector(1..dim);

  begin
    for i in cf0'range loop
      cf0(i) := cff(i)(cff(i)'first);
      cf1(i) := cff(i)(cff(i)'first+1);
      cf2(i) := cff(i)(cff(i)'first+2);
      pw1(i) := pwr(i)(pwr(i)'first);
      pw2(i) := pwr(i)(pwr(i)'first+1);
    end loop;
    Second_Derivative_Second_Order(pcf,pct,pdg,cf0,cf1,cf2,pw1,pw2,ycf,ydg,1);
    put_line("second derivative second order evaluation :");
    for i in ycf'range loop
      put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
    end loop;
  end Test_Second_Derivative_Second_Order;

  procedure Test_Third_Derivative_Second_Order
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                cff : in Standard_Complex_VecVecs.VecVec;
                pwr : in Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Tests the 3rd derivative 2nd order evaluation of a polynomial
  --   represented by (pcf, pct, pdg), evaluated at a series
  --   with coefficients in cff and powers in pwr.

    use Double_Ordered_Evaluations;

    dim : constant integer32 := cff'last;
    nbr : constant integer32 := pdg'last;
    size : constant integer32 := Size_Evaluation(dim,3,2,nbr);
    ycf : Standard_Complex_Vectors.Vector(1..size);
    ydg : Standard_Floating_Vectors.Vector(1..size);
    cf0,cf1,cf2 : Standard_Complex_Vectors.Vector(1..dim);
    pw1,pw2 : Standard_Floating_Vectors.Vector(1..dim);

  begin
    for i in cf0'range loop
      cf0(i) := cff(i)(cff(i)'first);
      cf1(i) := cff(i)(cff(i)'first+1);
      cf2(i) := cff(i)(cff(i)'first+2);
      pw1(i) := pwr(i)(pwr(i)'first);
      pw2(i) := pwr(i)(pwr(i)'first+1);
    end loop;
    Third_Derivative_Second_Order(pcf,pct,pdg,cf0,cf1,cf2,pw1,pw2,ycf,ydg,1);
    put_line("third derivative second order evaluation :");
    for i in ycf'range loop
      put(ycf(i)); put("  t^"); put(ydg(i)); new_line;
    end loop;
  end Test_Third_Derivative_Second_Order;

  procedure Generate_Input
              ( dim,nbr,ord : in integer32;
                hcf : out Standard_Complex_Vectors.Vector;
                hct : out Standard_Floating_Vectors.Vector;
                hdg : out Standard_Integer_VecVecs.VecVec;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Given the dimension dim, number nbr of monomials in a polynomial,
  --   and the order ord of the series, generates a Laurent homotopy
  --   and a power series solution.

    nbt : constant Standard_Integer_Vectors.Vector(1..dim)
        := (1..dim => ord);
    pdg : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Laurent_Homotopy.Random_Polynomial_Support(nbr,dim,-9,9);
    pcf : constant Standard_Complex_Vectors.Vector(1..nbr)
        := Standard_Random_Vectors.Random_Vector(1,nbr);
    pct : constant Standard_Floating_Vectors.Vector(1..nbr)
        := Double_Real_Powered_Series.Random_Leading_Powers(nbr);

  begin
    Double_Real_Powered_Series.Random_Power_Series(dim,nbt,cff,pwr);
    for i in cff'range loop
      put("-> a random power series "); put(i,1); put_line(" : ");
      put(cff(i)(0)); new_line;
      for j in 1..nbt(i) loop
        put(cff(i)(j)); put("  t^"); put(pwr(i)(j)); new_line;
      end loop;
    end loop;
    put_line("A random polynomial :");
    for i in 1..nbr loop
      put(pcf(i)); 
      put(" t^"); put(pct(i));
      put("  "); put(pdg(i)); new_line;
    end loop;
    Random_Laurent_Homotopy.Random_Homotopy_Polynomial
      (pdg,pcf,pct,cff,pwr,1,hdg,hcf,hct);
    Random_Laurent_Homotopy.Scale_Homotopy_Powers(hct);
    put_line("the homotopy polynomial :");
    for i in hcf'range loop
      put(hcf(i)); 
      put(" t^"); put(hct(i));
      put("  "); put(hdg(i)); new_line;
    end loop;
  end Generate_Input;

  procedure Test_First_Order ( dim,nbr,ord : in integer32 ) is

    cff : Standard_Complex_VecVecs.VecVec(1..dim);
    pwr : Standard_Floating_VecVecs.VecVec(1..dim);
    size : constant integer32 := (ord+2)*nbr;
    hcf : Standard_Complex_Vectors.Vector(1..size);
    hct : Standard_Floating_Vectors.Vector(1..size);
    hdg : Standard_Integer_VecVecs.VecVec(1..size);

  begin
    Generate_Input(dim,nbr,ord,hcf,hct,hdg,cff,pwr);
   -- for derivatives 1, 2, 3, we have two different ways to compute
    Test_First_Derivative_First_Order(hcf,hct,hdg,cff,pwr);
    Test_Second_Derivative_First_Order(hcf,hct,hdg,cff,pwr);
    Test_Third_Derivative_First_Order(hcf,hct,hdg,cff,pwr);
   -- for derivatives 4, 5, 6, 7, 8, there is only one way to compute
    for k in 4..8 loop
      Test_First_Order_Evaluation(hcf,hct,hdg,cff,pwr,integer32(k));
    end loop;
  end Test_First_Order;

  procedure Test_Second_Order ( dim,nbr,ord : in integer32 ) is

    cff : Standard_Complex_VecVecs.VecVec(1..dim);
    pwr : Standard_Floating_VecVecs.VecVec(1..dim);
    size : constant integer32 := (ord+2)*nbr;
    hcf : Standard_Complex_Vectors.Vector(1..size);
    hct : Standard_Floating_Vectors.Vector(1..size);
    hdg : Standard_Integer_VecVecs.VecVec(1..size);

  begin
    Generate_Input(dim,nbr,ord,hcf,hct,hdg,cff,pwr);
    Test_Second_Derivative_Second_Order(hcf,hct,hdg,cff,pwr);
    Test_Third_Derivative_Second_Order(hcf,hct,hdg,cff,pwr);
  end Test_Second_Order;

  procedure Main is

    dim,nbr,ord : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the number of monomials : "); get(nbr);
    put("Give the order of the series : "); get(ord);
   -- if ord > 0
   --  then Test_First_Order(dim,nbr,ord);
   -- end if;
    if ord > 1
     then Test_Second_Order(dim,nbr,ord);
    end if;
  end Main;

end Test_Ordered_Evaluations;
