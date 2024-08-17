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
    declare
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
    end;
  end Main;

begin
  Main;
end ts_expops;
