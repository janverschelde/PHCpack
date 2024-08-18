with Ada.Calendar;                       use Ada.Calendar;
with Ada.unchecked_conversion;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;

procedure ts_expops is

-- DESCRIPTION :
--   Tests arithmetic operations on complex exponential series.

  procedure Make_Random_Exponentials
              ( deg : in integer32;
                cff : out Standard_Complex_Vectors.Vector;
                exp : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Generates deg+1 random complex coefficients and
  --   corresponding exponents for the terms in the series.

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

  -- DESCRIPTION :
  --   Writes the exponential series to file.

  begin
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

  function Inverse ( cff : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the coefficients of the inverse series.

  -- REQUIRED : cff(0) is nonzero.

    res : Standard_Complex_Vectors.Vector(cff'range);

  begin
    res(0) := 1.0/cff(0);
    for i in 1..res'last loop
      res(i) := -cff(1)*res(i-1);
      for j in 2..i loop
        res(i) := res(i) - cff(j)*res(i-j);
      end loop;
      res(i) := res(i)/cff(0);
    end loop;
    return res;
  end Inverse;

  function Convolute ( a,b : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the product of the series with coefficients
  --   with the same exponents and same length.

  -- REQUIRED : a'last = b'last.

    res : Standard_Complex_Vectors.Vector(0..a'last);

  begin
    for i in 0..res'last loop
      res(i) := a(0)*b(i);
      for j in 1..i loop
        res(i) := res(i) + a(j)*b(i-j);
      end loop;
    end loop;
    return res;
  end Convolute;

  procedure Add ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the sum of two exponential series,
  --   truncated at the same degree.

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series.

  -- ON RETURN :
  --   ccf        coefficients of the sum;
  --   cxp        exponents of the sum.

    aix : integer32 := acf'first;
    bix : integer32 := bcf'first;

  begin
    for i in ccf'range loop
      if axp(aix) < bxp(bix) then
        ccf(i) := acf(aix);
        cxp(i) := axp(aix);
        aix := aix + 1;
      elsif axp(aix) > bxp(bix) then
        ccf(i) := bcf(bix);
        cxp(i) := bxp(bix);
        bix := bix + 1;
      else -- axp(aix) = bxp(bix) 
        ccf(i) := acf(aix) + bcf(bix);
        cxp(i) := axp(aix);
        aix := aix + 1;
        bix := bix + 1;
      end if;
    end loop;
  end Add;

  procedure Sub ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the difference of two exponential series,
  --   truncated at the same degree.

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series.

  -- ON RETURN :
  --   ccf        coefficients of the difference;
  --   cxp        exponents of the difference.

    aix : integer32 := acf'first;
    bix : integer32 := bcf'first;

  begin
    for i in ccf'range loop
      if axp(aix) < bxp(bix) then
        ccf(i) := acf(aix);
        cxp(i) := axp(aix);
        aix := aix + 1;
      elsif axp(aix) > bxp(bix) then
        ccf(i) := -bcf(bix);
        cxp(i) := bxp(bix);
        bix := bix + 1;
      else -- axp(aix) = bxp(bix) 
        ccf(i) := acf(aix) - bcf(bix);
        cxp(i) := axp(aix);
        aix := aix + 1;
        bix := bix + 1;
      end if;
    end loop;
  end Sub;

  procedure Inc ( acf : in out Standard_Complex_Vectors.Vector;
                  bcf : in Standard_Complex_Vectors.Vector;
                  axp : in out Standard_Floating_Vectors.Vector;
                  bxp : in Standard_Floating_Vectors.Vector;
                  wrkcf : in out Standard_Complex_Vectors.Vector;
                  wrkxp : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Increments the first series (acf, axp) with (bcf, bxp),
  --   truncated at the same degree.

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series.

  -- ON RETURN :
  --   acf        coefficients of the sum;
  --   axp        exponents of the sum.
  --   wrkcf      work space coefficients;
  --   wrkxp      work space exponents.

  begin
    Add(acf,bcf,axp,bxp,wrkcf,wrkxp);
    for i in acf'range loop
      acf(i) := wrkcf(i);
      axp(i) := wrkxp(i);
    end loop;
  end Inc;

  procedure Mul ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  prdcf,wrkcf : in out Standard_Complex_Vectors.Vector;
                  prdxp,wrkxp : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the product of two exponential series,
  --   truncated at the same degree.

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series.

  -- ON RETURN :
  --   ccf        coefficients of the product;
  --   cxp        exponents of the product.
  --   prdcf      work space coefficients for product term;
  --   prdxp      work space exponents for product term;
  --   wrkcf      work space coefficients for increment;
  --   wrkxp      work space exponents for increment.

  begin
    for i in ccf'range loop
      ccf(i) := acf(0)*bcf(i);
      cxp(i) := axp(0)+bxp(i);
    end loop;
    for i in 1..ccf'last loop
      for j in bcf'range loop
        prdcf(j) := acf(i)*bcf(j);
        prdxp(j) := axp(i)+bxp(j);
      end loop;
      Inc(ccf,prdcf,cxp,prdxp,wrkcf,wrkxp);
    end loop;
  end Mul;

  procedure Div ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  axp,bxp : in Standard_Floating_Vectors.Vector;
                  ccf : out Standard_Complex_Vectors.Vector;
                  cxp : out Standard_Floating_Vectors.Vector;
                  invbcf,prdcf,wrkcf : in out Standard_Complex_Vectors.Vector;
                  prdxp,wrkxp : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the quotient of two exponential series,
  --   truncated at the same degree.

  -- REQUIRED : bcf(0) is nonzero.

  -- ON ENTRY :
  --   acf        coefficients of the first series;
  --   bcf        coefficients of the second series;
  --   axp        exponents of the first series;
  --   bxp        exponents of the second series.

  -- ON RETURN :
  --   ccf        coefficients of the quotient;
  --   cxp        exponents of the quotient;
  --   invbcf     coefficients of inverse of second series;
  --   prdcf      work space coefficients for product term;
  --   prdxp      work space exponents for product term;
  --   wrkcf      work space coefficients for increment;
  --   wrkxp      work space exponents for increment.

  begin
    invbcf := Inverse(bcf);
    Mul(acf,invbcf,axp,bxp,ccf,cxp,prdcf,wrkcf,prdxp,wrkxp);
  end Div;

  procedure Test_Inverse ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Makes a random series truncated at degree deg,
  --   computes its inverse and the product to check.

    cff : Standard_Complex_Vectors.Vector(0..deg);
    invcff : Standard_Complex_Vectors.Vector(0..deg);
    prd : Standard_Complex_Vectors.Vector(0..deg);
    exp : Standard_Floating_Vectors.Vector(0..deg);

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
  end Test_Inverse;

  procedure Test_Sum ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Makes two random series truncated at degree deg,
  --   computes their sum and the difference to check.

    acf,bcf,sumcf,difcf : Standard_Complex_Vectors.Vector(0..deg);
    axp,bxp,sumxp,difxp : Standard_Floating_Vectors.Vector(0..deg);

  begin
    Make_Random_Exponentials(deg,acf,axp);
    put_line("The first series :");
    Write_Exponential_Series(standard_output,acf,axp);
    Make_Random_Exponentials(deg,bcf,bxp);
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
  end Test_Sum;

  procedure Test_Product ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Makes two random series truncated at degree deg,
  --   computes their product and the quotient to check.

    acf,bcf,prodcf,quotcf,difcf : Standard_Complex_Vectors.Vector(0..deg);
    axp,bxp,prodxp,quotxp,difxp : Standard_Floating_Vectors.Vector(0..deg);
    prdcf,wrkcf,invbcf : Standard_Complex_Vectors.Vector(0..deg);
    prdxp,wrkxp : Standard_Floating_Vectors.Vector(0..deg);

  begin
    Make_Random_Exponentials(deg,acf,axp);
    put_line("The first series :");
    Write_Exponential_Series(standard_output,acf,axp);
    Make_Random_Exponentials(deg,bcf,bxp);
    put_line("The second series :");
    Write_Exponential_Series(standard_output,bcf,bxp);
    Mul(acf,bcf,axp,bxp,prodcf,prodxp,prdcf,wrkcf,prdxp,wrkxp);
    put_line("The product of the two series :");
    Write_Exponential_Series(standard_output,prodcf,prodxp);
    Div(prodcf,bcf,prodxp,bxp,quotcf,quotxp,invbcf,prdcf,wrkcf,prdxp,wrkxp);
    put_line("After dividing second series from the product :");
    Write_Exponential_Series(standard_output,quotcf,quotxp);
    Sub(quotcf,acf,quotxp,axp,difcf,difxp);
    put_line("After subtracting first series from the difference :");
    Write_Exponential_Series(standard_output,difcf,difxp);
  end Test_Product;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of terms,
  --   generates random exponential series,
  --   and tests the arithmetical operations.

    deg : integer32 := 0;
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
    Test_Inverse(deg);
    Test_Sum(deg);
    Test_Product(deg);
  end Main;

begin
  Main;
end ts_expops;
