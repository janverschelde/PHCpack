with Ada.Calendar;                       use Ada.Calendar;
with Ada.unchecked_conversion;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;      use Standard_Complex_Vector_Norms;
with Double_Exponential_Arithmetic;      use Double_Exponential_Arithmetic;

package body Test_Double_Exponentials is

  procedure Make_Random_Exponentials
              ( deg : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                sxp : out Standard_Floating_Vectors.Vector ) is
  begin
    sxp(0) := 0.0;
    cff(0) := Standard_Random_Numbers.Random1;
    sxp(1) := 1.0 + abs(Standard_Random_Numbers.Random);
    cff(1) := Standard_Random_Numbers.Random1;
    for i in 2..deg loop
      cff(i) := Standard_Random_Numbers.Random1;
      sxp(i) := sxp(i-1) + abs(Standard_Random_Numbers.Random); 
    end loop;
  end Make_Random_Exponentials;

  function Is_Sorted ( xp : Standard_Floating_Vectors.Vector )
                     return boolean is
  begin
    for i in xp'first+1..xp'last loop
      if xp(i) < xp(i-1)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Sorted;

  procedure Normalize
              ( cff : in out Standard_Complex_Vectors.Vector;
                sxp : in out Standard_Floating_Vectors.Vector ) is

    cfftmp : Complex_Number;
    sxptmp : double_float;
    swapped : boolean;

  begin
    loop
      swapped := false;
      for i in sxp'first+1..sxp'last loop
        if sxp(i) < sxp(i-1) then
          sxptmp := sxp(i); sxp(i) := sxp(i-1); sxp(i-1) := sxptmp;
          cfftmp := cff(i); cff(i) := cff(i-1); cff(i-1) := cfftmp;
          swapped := true;
        end if;
      end loop;
      exit when not swapped;
    end loop;
  end Normalize;

  procedure Write_Exponential_Series
              ( file : in file_type;
                cff : in Standard_Complex_Vectors.Vector;
                sxp : in Standard_Floating_Vectors.Vector ) is
  begin
    if sxp(0) /= 0.0 then
      put("  (");
      if REAL_PART(cff(0)) >= 0.0
       then put(file," ");
      end if;
      put(file,REAL_PART(cff(0)),1,16,3);
      if IMAG_PART(cff(0)) >= 0.0
       then put(file," + ");
       else put(file," - ");
      end if;
      put(file,abs(IMAG_PART(cff(0))),1,16,3);
      put(file,"*I)*t^");
      put(file,sxp(0),1,16,3); new_line(file);
    else
      if REAL_PART(cff(0)) >= 0.0
       then put(file,"    ");
       else put(file,"   ");
      end if;
      put(file,REAL_PART(cff(0)),1,16,3);
      if IMAG_PART(cff(0)) >= 0.0
       then put(file," + ");
       else put(file," - ");
      end if;
      put(file,abs(IMAG_PART(cff(0))),1,16,3);
      put_line(file,"*I");
    end if;
    for i in 1..cff'last loop
      put(file,"+ (");
      if REAL_PART(cff(i)) >= 0.0
       then put(file," ");
      end if;
      put(file,REAL_PART(cff(i)),1,16,3);
      if IMAG_PART(cff(i)) >= 0.0
       then put(file," + ");
       else put(file," - ");
      end if;
      put(file,abs(IMAG_PART(cff(i))),1,16,3);
      put(file,"*I)*t^");
      put(file,sxp(i),1,16,3); new_line(file);
    end loop;
  end Write_Exponential_Series;

  function Extension_Degree ( alpha,beta : double_float ) return integer32 is

    res : integer32 := 0;

  begin
    if beta = 0.0 then
      return 0;
    else
      res := integer32(alpha/beta);
      while double_float(res)*beta < alpha loop
        res := res + 1;
      end loop;
      return res;
    end if;
  end Extension_Degree;

  function Extension_Degree
	     ( alpha : double_float;
               beta : Standard_Floating_Vectors.Vector) return integer32 is

    res,deg : integer32 := 0;

  begin
    for i in beta'range loop
      deg := Extension_Degree(alpha,beta(i));
      if deg > res
       then res := deg;
      end if;
    end loop;
    return res;
  end Extension_Degree;

  procedure Test_Inverse
              ( cff : in Standard_Complex_Vectors.Vector;
                sxp : in Standard_Floating_Vectors.Vector ) is                 

    deg : constant integer32 := cff'last;
    invcff : Standard_Complex_Vectors.Vector(0..deg);
    prd : Standard_Complex_Vectors.Vector(0..deg);
    nrm : double_float;

  begin
    put_line("An exponential series :");
    Write_Exponential_Series(standard_output,cff,sxp);
    invcff := Inverse(cff);
    put_line("The inverse of the exponential series :");
    Write_Exponential_Series(standard_output,invcff,sxp);
    prd := Convolute(cff,invcff);
    put_line("The product with the inverse of the exponential series :");
    Write_Exponential_Series(standard_output,prd,sxp);
    nrm := Max_Norm(prd);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
  end Test_Inverse;

  procedure Extend ( deg,extdeg : in integer32;
                     cff : in Standard_Complex_Vectors.Vector;
                     sxp : in Standard_Floating_Vectors.Vector;
                     extcff : out Standard_Complex_Vectors.Vector;
                     extsxp : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Extends the series with coefficients in cff and corresponding
  --   exponents of truncation degree to the new degree extdeg,
  --   returning in extcff and extsxp the extended series.

    newdeg : constant integer32 := deg + extdeg;
    idx,degidx : integer32;

  begin
    extcff(cff'range) := cff;
    extsxp(sxp'range) := sxp;
    extcff(cff'last+1..newdeg) := (cff'last+1..newdeg => create(0.0));
    idx := deg + 1;
    while idx <= newdeg loop
      if idx <= 2*deg then
        extsxp(idx) := 2.0*extsxp(idx-deg);
      else
        degidx := (idx mod deg) + 1;
        extsxp(idx) := extsxp(idx-deg) + extsxp(degidx);
      end if;
      idx := idx + 1;
    end loop;
    Normalize(extcff,extsxp);
  end Extend;

  procedure Test_Inverse ( deg : in integer32 ) is

    cff : Standard_Complex_Vectors.Vector(0..deg);
    sxp : Standard_Floating_Vectors.Vector(0..deg);
    extdeg : integer32 := 0;

  begin
    Make_Random_Exponentials(deg,cff,sxp);
    if not Is_Sorted(sxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("Random coefficients :"); put_line(cff);
    put_line("Random exponentials :"); put_line(sxp);
    Test_Inverse(cff,sxp);
    put("Give the extension degree : "); get(extdeg);
    if extdeg > 0 then
      declare
        newdeg : constant integer32 := deg + extdeg;
        extcff : Standard_Complex_Vectors.Vector(0..newdeg);
        extsxp : Standard_Floating_Vectors.Vector(0..newdeg);
      begin
        Extend(deg,extdeg,cff,sxp,extcff,extsxp);
        Test_Inverse(extcff,extsxp);
      end;
    end if;
  end Test_Inverse;

  procedure Test_Sum ( adeg,bdeg : in integer32 ) is

    sumdeg : constant integer32 := adeg + bdeg;
    acf : Standard_Complex_Vectors.Vector(0..adeg);
    axp : Standard_Floating_Vectors.Vector(0..adeg);
    bcf : Standard_Complex_Vectors.Vector(0..bdeg);
    bxp : Standard_Floating_Vectors.Vector(0..bdeg);
    sumcf,difcf : Standard_Complex_Vectors.Vector(0..sumdeg);
    sumxp,difxp : Standard_Floating_Vectors.Vector(0..sumdeg);
    nrm : double_float;

  begin
    Make_Random_Exponentials(adeg,acf,axp);
    if not Is_Sorted(axp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("The first series :");
    Write_Exponential_Series(standard_output,acf,axp);
    Make_Random_Exponentials(bdeg,bcf,bxp);
    if not Is_Sorted(bxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("The second series :");
    Write_Exponential_Series(standard_output,bcf,bxp);
    Add(adeg,bdeg,sumdeg,acf,bcf,axp,bxp,sumcf,sumxp);
    if not Is_Sorted(sumxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("The sum of the two series :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    Sub(sumdeg,bdeg,sumdeg,sumcf,bcf,sumxp,bxp,difcf,difxp);
    if not Is_Sorted(difxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("After subtracting second series from the sum :");
    Write_Exponential_Series(standard_output,difcf,difxp);
    Sub(sumdeg,adeg,sumdeg,difcf,acf,difxp,axp,sumcf,sumxp);
    if not Is_Sorted(sumxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("After subtracting first series from the difference :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    nrm := Max_Norm(sumcf);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
    Add(adeg,bdeg,sumdeg,acf,bcf,axp,bxp,sumcf,sumxp);
    if not Is_Sorted(sumxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("The sum of the two series :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    Sub(sumdeg,adeg,sumdeg,sumcf,acf,sumxp,axp,difcf,difxp);
    if not Is_Sorted(difxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("After subtracting first series from the sum :");
    Write_Exponential_Series(standard_output,difcf,difxp);
    Sub(sumdeg,bdeg,sumdeg,difcf,bcf,difxp,bxp,sumcf,sumxp);
    if not Is_Sorted(sumxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("After subtracting second series from the difference :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    nrm := Max_Norm(sumcf);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
  end Test_Sum;

  procedure Test_Multiplicative_Commutativity
              ( adeg,bdeg : in integer32;
                acf,bcf : in Standard_Complex_Vectors.Vector;
                axp,bxp : in Standard_Floating_Vectors.Vector ) is

    proddeg : constant integer32 := (adeg+1)*(bdeg+1) - 1;
    abprodcf,baprodcf : Standard_Complex_Vectors.Vector(0..proddeg);
    abprodxp,baprodxp : Standard_Floating_Vectors.Vector(0..proddeg);
    prdcf,wrkcf,difcf : Standard_Complex_Vectors.Vector(0..proddeg);
    prdxp,wrkxp,difxp : Standard_Floating_Vectors.Vector(0..proddeg);
    nrm : double_float;

  begin
    Mul(adeg,bdeg,proddeg,acf,bcf,axp,bxp,
        abprodcf,abprodxp,prdcf,wrkcf,prdxp,wrkxp);
    if not Is_Sorted(abprodxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("The product of the first with the second series :");
    Write_Exponential_Series
      (standard_output,abprodcf(0..proddeg-1),abprodxp(0..proddeg-1));
    Mul(bdeg,adeg,proddeg,bcf,acf,bxp,axp,
        baprodcf,baprodxp,prdcf,wrkcf,prdxp,wrkxp);
    if not Is_Sorted(baprodxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("The product of the second with the first series :");
    Write_Exponential_Series(standard_output,baprodcf,baprodxp);
    Sub(proddeg,proddeg,proddeg,abprodcf,baprodcf,abprodxp,baprodxp,
        difcf,difxp);
    if not Is_Sorted(difxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("After subtracting the products :");
    Write_Exponential_Series(standard_output,difcf,difxp);
    nrm := Max_Norm(difcf);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
  end Test_Multiplicative_Commutativity;

  procedure Test_Product
              ( adeg,bdeg : in integer32;
                acf,bcf : in Standard_Complex_Vectors.Vector;
                axp,bxp : in Standard_Floating_Vectors.Vector ) is

    extdeg : constant integer32 
           := Extension_Degree(axp(axp'last),bxp(0..bdeg));
    proddeg : constant integer32 := (adeg+1)*(bdeg+extdeg+1) - 1;
    ebcf : Standard_Complex_Vectors.Vector(0..bdeg + extdeg);
    ebxp : Standard_Floating_Vectors.Vector(0..bdeg + extdeg);
    prodcf : Standard_Complex_Vectors.Vector(0..proddeg);
    prodxp : Standard_Floating_Vectors.Vector(0..proddeg);
    quotdeg : constant integer32 := (proddeg+1)*(bdeg+extdeg+1) - 1;
    quotcf,difcf : Standard_Complex_Vectors.Vector(0..quotdeg);
    quotxp,difxp : Standard_Floating_Vectors.Vector(0..quotdeg);
    invbcf : Standard_Complex_Vectors.Vector(0..bdeg+extdeg);
    prdcf,wrkcf : Standard_Complex_Vectors.Vector(0..quotdeg);
    prdxp,wrkxp : Standard_Floating_Vectors.Vector(0..quotdeg);
    nrm : double_float;
    ans : character;

  begin
    put("-> the computed extension degree : "); put(extdeg,1); new_line;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    Extend(bdeg,extdeg,bcf,bxp,ebcf,ebxp);
    put_line("The second series, extended :");
    Write_Exponential_Series(standard_output,ebcf,ebxp);
    if not Is_Sorted(ebxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    Mul(adeg,bdeg+extdeg,proddeg,acf,ebcf,axp,ebxp,
        prodcf,prodxp,prdcf,wrkcf,prdxp,wrkxp);
    if not Is_Sorted(prodxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("The product of the two series :");
    Write_Exponential_Series(standard_output,prodcf,prodxp);
    Div(proddeg-1,bdeg+extdeg,quotdeg,prodcf,ebcf,prodxp,ebxp,
        quotcf,quotxp,invbcf,prdcf,wrkcf,prdxp,wrkxp);
    if not Is_Sorted(quotxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("the inverse of the second series :");
    Write_Exponential_Series(standard_output,invbcf,ebxp);
    put_line("After dividing second series from the product :");
    Write_Exponential_Series(standard_output,quotcf,quotxp);
    Sub(adeg,adeg,adeg,quotcf(0..adeg),acf,quotxp(0..adeg),axp,difcf,difxp);
    if not Is_Sorted(difxp(0..adeg))
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("After subtracting first series from the difference :");
    Write_Exponential_Series(standard_output,difcf(0..adeg),difxp(0..adeg));
    nrm := Max_Norm(difcf(0..adeg));
    put("-> max norm of the coefficients :"); put(nrm); new_line;
   -- put_line("the first series :");
   -- Write_Exponential_Series(standard_output,acf,axp);
  end Test_Product;

  procedure Test_Product ( adeg,bdeg : in integer32 ) is

    acf : Standard_Complex_Vectors.Vector(0..adeg);
    axp : Standard_Floating_Vectors.Vector(0..adeg);
    bcf : Standard_Complex_Vectors.Vector(0..bdeg);
    bxp : Standard_Floating_Vectors.Vector(0..bdeg);
    ans : character;

  begin
    Make_Random_Exponentials(adeg,acf,axp);
    put_line("The first series :");
    Write_Exponential_Series(standard_output,acf,axp);
    if not Is_Sorted(axp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    Make_Random_Exponentials(bdeg,bcf,bxp);
   -- make sure exponents of second series are large enough
   -- for i in 1..bxp'last loop
   --   bxp(i) := bxp(i) + axp(axp'last);
   -- end loop;
    put_line("The second series :");
    Write_Exponential_Series(standard_output,bcf,bxp);
    if not Is_Sorted(bxp)
     then put_line("Exponents are NOT in increasing order!");
    end if;
    put_line("Testing commutativity ...");
    Test_Multiplicative_Commutativity(adeg,bdeg,acf,bcf,axp,bxp);
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    Test_Product(adeg,bdeg,acf,bcf,axp,bxp);
  end Test_Product;

  procedure Main is

    ans : character;
    deg,adeg,bdeg : integer32 := 0;
    now : constant Time := Clock;
    year,month,day : integer;
    seconds : Duration;
    natsecs,fliprem,flipquo : natural32;
    modulus1 : constant natural32 := 32768;

  begin
    Split(now,year,month,day,seconds);
    natsecs := natural32(seconds);
    fliprem := natsecs rem modulus1;
    flipquo := natsecs/modulus1;
    natsecs := fliprem*modulus1 + flipquo;
    Standard_Random_Numbers.Set_Seed(natsecs);
    new_line;
    put("Give the truncation degree of the exponential series : ");
    get(deg);
    put("-> generating series of degree "); put(deg,1); put_line(" ...");
    put_line("****** testing the inverse ******");
    Test_Inverse(deg);
    put_line("****** testing the sum ******");
    put("Give the degree of the first series : "); get(adeg);
    put("Give the degree of the second series : "); get(bdeg);
    Test_Sum(adeg,bdeg);
    put("Continue with the test on the product ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("****** testing the product ******");
      Test_Product(adeg,bdeg);
    end if;
  end Main;

end Test_Double_Exponentials;
