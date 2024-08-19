with Ada.Calendar;                       use Ada.Calendar;
with Ada.unchecked_conversion;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;      use Standard_Complex_Vector_Norms;

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
    cix : integer32 := ccf'first;

  begin
    while cix <= ccf'last loop
      if axp(aix) < bxp(bix) then
        ccf(cix) := acf(aix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        cix := cix + 1;
      elsif axp(aix) > bxp(bix) then
        ccf(cix) := bcf(bix);
        cxp(cix) := bxp(bix);
        bix := bix + 1;
        cix := cix + 1;
      else -- axp(aix) = bxp(bix) 
        ccf(cix) := acf(aix) + bcf(bix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        bix := bix + 1;
        if AbsVal(ccf(cix)) > 1.0e-14
         then cix := cix + 1;
        end if;
      end if;
      exit when (aix > acf'last) or (bix > bcf'last);
    end loop;
    if cix <= ccf'last then
      if aix <= acf'last then
        while cix <= ccf'last loop
          ccf(cix) := acf(aix);
          cxp(cix) := axp(aix);
          cix := cix + 1;
          aix := aix + 1;
          exit when (aix > acf'last);
        end loop;
      elsif bix <= bcf'last then
        while cix <= ccf'last loop
          ccf(cix) := bcf(bix);
          cxp(cix) := bxp(bix);
          cix := cix + 1;
          bix := bix + 1;
          exit when (bix > bcf'last);
        end loop;
      end if;
      while cix <= ccf'last loop
        ccf(cix) := create(0.0);
        cxp(cix) := cxp(cix-1) + 1.0;
        cix := cix + 1;
      end loop;
    end if;
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
    cix : integer32 := ccf'first;

  begin
    while cix <= ccf'last loop
      if axp(aix) < bxp(bix) then
        ccf(cix) := acf(aix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        cix := cix + 1;
      elsif axp(aix) > bxp(bix) then
        ccf(cix) := -bcf(bix);
        cxp(cix) := bxp(bix);
        bix := bix + 1;
        cix := cix + 1;
      else -- axp(aix) = bxp(bix) 
        ccf(cix) := acf(aix) - bcf(bix);
        cxp(cix) := axp(aix);
        aix := aix + 1;
        bix := bix + 1;
        if AbsVal(ccf(cix)) > 1.0e-14 then
          cix := cix + 1;
        else
          put("aix : "); put(aix,1);
          put("  bix : "); put(bix,1);
          put("  cix : "); put(cix,1); new_line;
          put(ccf(cix)); new_line;
        end if;
      end if;
      exit when (aix > acf'last) or (bix > bcf'last);
    end loop;
    if cix <= ccf'last then
     -- put("at end, aix : "); put(aix,1);
     -- put("  bix : "); put(bix,1);
     -- put("  cix : "); put(cix,1); new_line;
      if aix <= acf'last then
       -- put_line("copying from a ...");
        while cix <= ccf'last loop
          ccf(cix) := acf(aix);
          cxp(cix) := axp(aix);
          cix := cix + 1;
          aix := aix + 1;
          exit when (aix > acf'last);
        end loop;
      elsif bix <= bcf'last then
       -- put_line("copying from b ...");
        while cix <= ccf'last loop
          ccf(cix) := -bcf(bix);
          cxp(cix) := bxp(bix);
          cix := cix + 1;
          bix := bix + 1;
          exit when (bix > bcf'last);
        end loop;
      end if;
      while cix <= ccf'last loop
        ccf(cix) := create(0.0);
        cxp(cix) := cxp(cix-1) + 1.0;
        cix := cix + 1;
      end loop;
    end if;
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

    deg : constant integer32 := acf'last;
    cix,pix,wix : integer32;

  begin
    for i in ccf'range loop
      exit when (i > bcf'last);
      ccf(i) := acf(0)*bcf(i);
      cxp(i) := axp(0)+bxp(i);
    end loop;
    for i in 1..deg loop
      for j in bcf'range loop
        exit when (i > acf'last);
        prdcf(j) := acf(i)*bcf(j);
        prdxp(j) := axp(i)+bxp(j);
      end loop;
      cix := ccf'first;
      pix := prdcf'first;
      wix := wrkcf'first;
      while wix <= (i+1)*deg loop
        if cxp(cix) < prdxp(pix) then
          wrkcf(wix) := ccf(cix);
          wrkxp(wix) := cxp(cix);
          cix := cix + 1;
          wix := wix + 1;
        elsif cxp(cix) > prdxp(pix) then
          wrkcf(wix) := prdcf(pix);
          wrkxp(wix) := prdxp(pix);
          pix := pix + 1;
          wix := wix + 1;
        else -- cxp(cix) = prdxp(pix)
          wrkcf(wix) := ccf(cix) + prdcf(pix);
          wrkxp(wix) := cxp(cix);
          cix := cix + 1;
          pix := pix + 1;
          wix := wix + 1;
        end if;
        exit when (pix > deg) or (cix > i*deg);
      end loop;
      if wix <= (i+1)*deg then
        if pix <= deg then
          while wix <= (i+1)*deg loop
            wrkcf(wix) := prdcf(pix);
            wrkxp(wix) := prdxp(pix);
            wix := wix + 1;
            pix := pix + 1;
            exit when (pix > deg);
          end loop;
        elsif cix <= i*deg then
          while wix <= (i+1)*deg loop
            wrkcf(wix) := ccf(cix);
            wrkxp(wix) := cxp(cix);
            wix := wix + 1;
            cix := cix + 1;
            exit when (cix > i*deg);
          end loop;
        end if;
      end if;
      put("i = "); put(i,1); 
      put("  deg = "); put(deg,1); 
      put("  (i+1)*deg = "); put((i+1)*deg,1); new_line;
      for j in 0..(i+1)*deg loop
        ccf(j) := wrkcf(j);
        cxp(j) := wrkxp(j);
      end loop;
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

  procedure Test_Sum ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Makes two random series truncated at degree deg,
  --   computes their sum and the difference to check.

    acf,bcf : Standard_Complex_Vectors.Vector(0..deg);
    sumcf,difcf : Standard_Complex_Vectors.Vector(0..2*deg);
    axp,bxp : Standard_Floating_Vectors.Vector(0..deg);
    sumxp,difxp : Standard_Floating_Vectors.Vector(0..2*deg);
    nrm : double_float;

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

  procedure Test_Product ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Makes two random series truncated at degree deg,
  --   computes their product and the quotient to check.

    maxdeg : constant integer32 := (deg+1)*deg;
    acf,bcf : Standard_Complex_Vectors.Vector(0..deg);
    prodcf,quotcf,difcf : Standard_Complex_Vectors.Vector(0..maxdeg);
    axp,bxp : Standard_Floating_Vectors.Vector(0..deg);
    prodxp,quotxp,difxp : Standard_Floating_Vectors.Vector(0..maxdeg);
    invbcf : Standard_Complex_Vectors.Vector(0..deg);
    prdcf,wrkcf : Standard_Complex_Vectors.Vector(0..maxdeg);
    prdxp,wrkxp : Standard_Floating_Vectors.Vector(0..maxdeg);
    nrm : double_float;

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
    Div(prodcf(0..deg),bcf,prodxp(0..deg),bxp,
        quotcf,quotxp,invbcf,prdcf,wrkcf,prdxp,wrkxp);
    put_line("After dividing second series from the product :");
    Write_Exponential_Series(standard_output,quotcf,quotxp);
    Sub(quotcf(0..deg),acf,quotxp(0..deg),axp,difcf,difxp);
    put_line("After subtracting first series from the difference :");
    Write_Exponential_Series(standard_output,difcf(0..deg),difxp(0..deg));
    nrm := Max_Norm(difcf);
    put("-> max norm of the coefficients :"); put(nrm); new_line;
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
    put_line("****** testing the inverse ******");
    Test_Inverse(deg);
    put_line("****** testing the sum ******");
    Test_Sum(deg);
    put_line("****** testing the product ******");
    Test_Product(deg);
  end Main;

begin
  Main;
end ts_expops;
