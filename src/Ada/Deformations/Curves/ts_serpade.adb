with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;

procedure ts_serpade is

-- DESCRIPTION :
--   Interactive development of the rational approximation of a function,
--   given the power series at the origin.

  function log1plusx ( dim : integer32 )
                     return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the first dim+1 coefficients of the series of log(1+x)
  --   as a vector of range 0..dim.

    res : Standard_Complex_Vectors.Vector(0..dim);

  begin
    res(0) := Create(0.0);
    for k in 1..dim loop
      res(k) := Create(1.0/double_float(k));
      if k mod 2 = 0
       then Min(res(k));
      end if;
    end loop;
    return res;
  end log1plusx;

  procedure Standard_Denominator_System
              ( numdeg,dendeg : in integer32; 
                cff : in Standard_Complex_Vectors.Vector;
                mat : out Standard_Complex_Matrices.Matrix;
                rhs : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Defines the coefficient matrix and the right hand side vector
  --   for the denominator of a rational approximation.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   cff      coefficients of the series, of range 0..numdeg+dendeg.

  -- ON RETURN :
  --   mat      square coefficient matrix of range 1..dendeg.
  --   rhs      right hand side vector of range 1..dendeg.

    dim : constant integer32 := numdeg+dendeg;
    idx : integer32 := 0;

  begin
    for i in 1..dendeg loop
      idx := numdeg - dendeg + i;
      for j in 1..dendeg loop
        if idx < 1
         then mat(i,j) := Create(0.0);
         else mat(i,j) := cff(idx);
        end if;
        idx := idx + 1;
      end loop;
    end loop;
    idx := 0;
    for k in numdeg+1..dim loop
      idx := idx + 1;
      rhs(idx) := -cff(k);
    end loop;
  end Standard_Denominator_System;

  function Standard_Numerator_Coefficients
              ( numdeg,dendeg : integer32;
                dencff,sercff : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the coefficients of the numerator in the rational
  --   approximation as a vector of range 0..numdeg.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   dencff   coefficients of the denominator;
  --   sercff   coefficients of the power series.

    res : Standard_Complex_Vectors.Vector(0..numdeg);
    mindeg : integer32;

  begin
    res(0) := sercff(0);
    if dendeg <= numdeg
     then mindeg := dendeg;
     else mindeg := numdeg;
    end if;
    for i in 1..numdeg loop
      res(i) := sercff(i);
      for j in 1..i loop
        exit when (j > mindeg);
        res(i) := res(i) + dencff(j)*sercff(i-j);
      end loop; 
    end loop;
    return res;
  end Standard_Numerator_Coefficients;

  function Standard_Evaluate
              ( p : Standard_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Evaluates the polynomial with coefficients in p at x.

    res : Complex_Number := p(p'last);

  begin
    for k in reverse 0..p'last-1 loop
      res := res*x + p(k);
    end loop;
    return res;
  end Standard_Evaluate;

  function Standard_Evaluate
              ( num,den : Standard_Complex_Vectors.Vector;
                x : Complex_Number ) return Complex_Number is

  -- DESCRIPTION :
  --   Given the coefficients of the numerator in num
  --   and the coefficients of the denominator in den,
  --   return the value of the rational approximation at x.

    numval : constant Complex_Number := Standard_Evaluate(num,x);
    denval : constant Complex_Number := Standard_Evaluate(den,x);
    res : constant Complex_Number := numval/denval;

  begin
    return res;
  end Standard_Evaluate;

  procedure Standard_Test
              ( numdeg,dendeg : in integer32;
                cff : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   for a series with coefficients given in cff.
  --   The simplest textbook definition is applied in the computation.

  -- ON ENTRY :
  --   numdeg   degree of the numerator;
  --   dendeg   degree of the denominator;
  --   cff      coefficients of the power series.

  -- REQUIRED : cff'range = 0..numdeg+dendeg.

    dim : constant integer32 := numdeg + dendeg;
    mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
    rhs : Standard_Complex_Vectors.Vector(1..dendeg);
    ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
    info : integer32;
    dencff : Standard_Complex_Vectors.Vector(0..dendeg);
    numcff : Standard_Complex_Vectors.Vector(0..numdeg);
    pnt,eva : Complex_Number;
    chkpnt : constant double_float
           := Standard_Mathematical_Functions.LN(1.1);

  begin
    put_line("The coefficient vector for log(1+x) :"); put_line(cff);
    Standard_Denominator_System(numdeg,dendeg,cff,mat,rhs);
    lufac(mat,dendeg,ipvt,info);
    lusolve(mat,dendeg,ipvt,rhs);
    dencff(0) := Create(1.0);
    for i in 1..dendeg loop
      dencff(i) := rhs(dendeg-i+1);
    end loop;
    put_line("The coefficients of the denominator :"); put_line(dencff);
    numcff := Standard_Numerator_Coefficients(numdeg,dendeg,dencff,cff);
    put_line("The coefficients of the numerator :"); put_line(numcff);
    pnt := Create(0.1);
    eva := Standard_Evaluate(numcff,dencff,pnt);
    put("The value at 1.1      :"); put(eva); new_line;
    put("The value of log(1.1) :"); put(chkpnt); new_line;
  end Standard_Test;

  procedure Standard_Logarithm_Test ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the construction in standard floating point arithmetic
  --   on the natural logarithm of 1 + x.

    dim : constant integer32 := numdeg + dendeg;
    cff : constant Standard_Complex_Vectors.Vector(0..dim)
        := log1plusx(dim);

  begin
    Standard_Test(numdeg,dendeg,cff);
  end Standard_Logarithm_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension
  --   and the coefficient vector of the series.

    degnum : integer32 := 0; -- degree numerator
    degden : integer32 := 0; -- degree denominator
    dim : integer32 := 0;    -- degnum + degden

  begin
    new_line;
    put_line("Rational approximation for a function ...");
    put("Give the degree of the numerator : "); get(degnum);
    put("Give the degree of the denominator : "); get(degden);
    dim := degnum + degden;
    put("The dimension : "); put(dim,1); new_line;
    Standard_Logarithm_Test(degnum,degden);
  end Main;

begin
  Main;
end ts_serpade;
