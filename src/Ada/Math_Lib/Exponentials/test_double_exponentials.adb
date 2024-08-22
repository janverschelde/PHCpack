with Ada.Calendar;                       use Ada.Calendar;
with Ada.unchecked_conversion;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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
                exp : out Standard_Floating_Vectors.Vector ) is
  begin
    exp(0) := 0.0;
    cff(0) := Standard_Random_Numbers.Random1;
    for i in 1..deg loop
      cff(i) := Standard_Random_Numbers.Random1;
      exp(i) := exp(i-1) + abs(Standard_Random_Numbers.Random); 
    end loop;
  end Make_Random_Exponentials;

  procedure Write_Exponential_Series
              ( file : in file_type;
                cff : in Standard_Complex_Vectors.Vector;
                exp : in Standard_Floating_Vectors.Vector ) is
  begin
    if exp(0) /= 0.0 then
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
      put(file,exp(0),1,16,3); new_line(file);
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
      put(file,exp(i),1,16,3); new_line(file);
    end loop;
  end Write_Exponential_Series;

  procedure Test_Inverse ( deg : in integer32 ) is

    cff : Standard_Complex_Vectors.Vector(0..deg);
    invcff : Standard_Complex_Vectors.Vector(0..deg);
    prd : Standard_Complex_Vectors.Vector(0..deg);
    exp : Standard_Floating_Vectors.Vector(0..deg);
    nrm : double_float;

  begin
    Make_Random_Exponentials(deg,cff,exp);
    put_line("Random coefficients :"); put_line(cff);
    put_line("Random exponentials :"); put_line(exp);
    put_line("A random exponential series :");
    Write_Exponential_Series(standard_output,cff,exp);
    invcff := Inverse(cff);
    put_line("The inverse of the exponential series :");
    Write_Exponential_Series(standard_output,invcff,exp);
    prd := Convolute(cff,invcff);
    put_line("The product with the inverse of the exponential series :");
    Write_Exponential_Series(standard_output,prd,exp);
    nrm := Max_Norm(prd);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
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
    put_line("The first series :");
    Write_Exponential_Series(standard_output,acf,axp);
    Make_Random_Exponentials(bdeg,bcf,bxp);
    put_line("The second series :");
    Write_Exponential_Series(standard_output,bcf,bxp);
    Add(acf,bcf,axp,bxp,sumcf,sumxp);
    put_line("The sum of the two series :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    Sub(sumcf,bcf,sumxp,bxp,difcf,difxp);
    put_line("After subtracting second series from the sum :");
    Write_Exponential_Series(standard_output,difcf,difxp);
    Sub(difcf,acf,difxp,axp,sumcf,sumxp);
    put_line("After subtracting first series from the difference :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    nrm := Max_Norm(sumcf);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
    Add(acf,bcf,axp,bxp,sumcf,sumxp);
    put_line("The sum of the two series :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    Sub(sumcf,acf,sumxp,axp,difcf,difxp);
    put_line("After subtracting first series from the sum :");
    Write_Exponential_Series(standard_output,difcf,difxp);
    Sub(difcf,bcf,difxp,bxp,sumcf,sumxp);
    put_line("After subtracting second series from the difference :");
    Write_Exponential_Series(standard_output,sumcf,sumxp);
    nrm := Max_Norm(sumcf);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
  end Test_Sum;

  procedure Test_Product ( adeg,bdeg : in integer32 ) is

    maxdeg : constant integer32 := (adeg+1)*bdeg;
    acf : Standard_Complex_Vectors.Vector(0..adeg);
    axp : Standard_Floating_Vectors.Vector(0..adeg);
    bcf : Standard_Complex_Vectors.Vector(0..bdeg);
    bxp : Standard_Floating_Vectors.Vector(0..bdeg);
    prodcf,quotcf,difcf : Standard_Complex_Vectors.Vector(0..maxdeg);
    prodxp,quotxp,difxp : Standard_Floating_Vectors.Vector(0..maxdeg);
    invbcf : Standard_Complex_Vectors.Vector(0..bdeg);
    prdcf,wrkcf : Standard_Complex_Vectors.Vector(0..maxdeg);
    prdxp,wrkxp : Standard_Floating_Vectors.Vector(0..maxdeg);
    nrm : double_float;

  begin
    Make_Random_Exponentials(adeg,acf,axp);
    put_line("The first series :");
    Write_Exponential_Series(standard_output,acf,axp);
    Make_Random_Exponentials(bdeg,bcf,bxp);
   -- make sure exponents of second series are large enough
    for i in 1..bxp'last loop
      bxp(i) := bxp(i) + axp(axp'last);
    end loop;
    put_line("The second series :");
    Write_Exponential_Series(standard_output,bcf,bxp);
    Mul(acf,bcf,axp,bxp,prodcf,prodxp,prdcf,wrkcf,prdxp,wrkxp);
    put_line("The product of the two series :");
    Write_Exponential_Series(standard_output,prodcf,prodxp);
    Div(prodcf(0..adeg),bcf,prodxp(0..adeg),bxp,
        quotcf,quotxp,invbcf,prdcf,wrkcf,prdxp,wrkxp);
    put_line("After dividing second series from the product :");
    Write_Exponential_Series(standard_output,quotcf,quotxp);
    Sub(quotcf(0..adeg),acf,quotxp(0..adeg),axp,difcf,difxp);
    put_line("After subtracting first series from the difference :");
    Write_Exponential_Series(standard_output,difcf(0..adeg),difxp(0..adeg));
    nrm := Max_Norm(difcf);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
  end Test_Product;

  procedure Main is

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
    put("Give the degree of the first series : "); get(adeg);
    put("Give the degree of the second series : "); get(bdeg);
    put_line("****** testing the sum ******");
    Test_Sum(adeg,bdeg);
    put_line("****** testing the product ******");
    Test_Product(adeg,bdeg);
  end Main;

end Test_Double_Exponentials;
