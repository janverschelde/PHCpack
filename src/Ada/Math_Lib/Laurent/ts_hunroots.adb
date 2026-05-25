with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Double_Real_Powered_Series;
with Test_Real_Powered_Series;

procedure ts_hunroots is

-- DESCRIPTION :
--   Tests the square root computation of real power series.

  procedure Add ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the sum of (acf, apw) and (bcf, bpw).
  --   The addition of two random series is defined by the
  --   sum of their constants and the merge sort of the terms.
  --   For the result to be correct, the size of (cff, pwt)
  --   must be the sum of the sizes of (acf, apw) and (bcf, bpw).

    idx : integer32 := pwt'first;
    adx : integer32 := apw'first;
    bdx : integer32 := bpw'first;

  begin
    cff(0) := acf(0) + bcf(0);
    while (adx <= apw'last) and (bdx <= bpw'last) loop
      exit when (idx > pwt'last);
      if apw(adx) < bpw(bdx) then
        pwt(idx) := apw(adx);
        cff(idx) := acf(adx);
        adx := adx + 1;
      else
        pwt(idx) := bpw(bdx);
        cff(idx) := bcf(bdx);
        bdx := bdx + 1;
      end if;
      idx := idx + 1;
    end loop;
    for i in adx..apw'last loop
      exit when (idx > pwt'last);
      pwt(idx) := apw(i);
      cff(idx) := acf(i);
      idx := idx + 1;
    end loop;
    for i in bdx..bpw'last loop
      exit when (idx > pwt'last);
      pwt(idx) := bpw(i);
      cff(idx) := bcf(i);
      idx := idx + 1;
    end loop;
    Double_Real_Powered_Series.normalize(cff,pwt);
  end Add;

  procedure Sub ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the result of (acf, apw) - (bcf, bpw).
  --   For the result to be correct, the size of (cff, pwt)
  --   must be the sum of the sizes of (acf, apw) and (bcf, bpw).

    minbcf : Standard_Complex_Vectors.Vector(bcf'range);

  begin
    for i in bcf'range loop
      minbcf(i) := -bcf(i);
    end loop;
    Add(acf,minbcf,apw,bpw,cff,pwt);
  end Sub;

  procedure Mul ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the product of (acf, apw) and (bcf, bpw).
  --   For the result to be correct, the size of (cff, pwt)
  --   must be the product of the one plus the sizes of (acf, apw)
  --   and (bcf, bpw) minus one (for the constant).

    idx : integer32 := 0;

  begin
    cff(0) := acf(0)*bcf(0);  -- multiply with constant of first series
    for j in bpw'range loop
      exit when (j > cff'last);
      cff(j) := acf(0)*bcf(j);
      pwt(j) := bpw(j);
    end loop;
    idx := bpw'last;
    for i in apw'range loop  -- multiply with i-th term of first series
      idx := idx + 1;
      exit when (idx > cff'last);
      cff(idx) := acf(i)*bcf(0);
      pwt(idx) := apw(i);
      for j in bpw'range loop
        idx := idx + 1;
        exit when (idx > cff'last);
        cff(idx) := acf(i)*bcf(j);
        pwt(idx) := apw(i)+bpw(j);
      end loop;
    end loop;
    Double_Real_Powered_Series.sort(pwt,cff);
    Double_Real_Powered_Series.normalize(cff,pwt);
  end Mul;

  function Equal ( acf,bcf : Standard_Complex_Vectors.Vector;
                   apw,bpw : Standard_Floating_Vectors.Vector;
                   tol : double_float := 1.0E-12 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the sum of the componentwise errors of
  --   the coefficients and the powers is less than the tolerance.
  --   Both series are expected to have the same size.

    sumerr : double_float := AbsVal(acf(0) - bcf(0));

  begin
    if sumerr > tol then
      return false;
    else
      for i in apw'range loop
        sumerr := sumerr + AbsVal(acf(i) - bcf(i));
        if sumerr > tol
         then return false;
        end if;
        sumerr := sumerr + abs(apw(i) - bpw(i));
        if sumerr > tol
         then return false;
        end if;
      end loop;
    end if;
    return true;
  end Equal;

  procedure Test_Random_Series ( size : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random real power series of the given size.

    cff : Standard_Complex_Vectors.Vector(0..size);
    pwt : Standard_Floating_Vectors.Vector(1..size);

  begin
    put("-> generating a series of size "); put(size,1); put_line(" ...");
    Test_Real_Powered_Series.random_series(size,cff,pwt);
    Test_Real_Powered_Series.write(cff,pwt);
  end Test_Random_Series;

  procedure Test_Addition ( size : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random series of the given size,
  --   tests their addition and the subtraction.

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    sumcf : Standard_Complex_Vectors.Vector(0..2*size);
    sumpw : Standard_Floating_Vectors.Vector(1..2*size);
    difcf : Standard_Complex_Vectors.Vector(0..3*size);
    difpw : Standard_Floating_Vectors.Vector(1..3*size);
    equ : boolean;

  begin
    put("-> adding two series of size "); put(size,1); put_line(" ...");
    Test_Real_Powered_Series.random_series(size,acf,apw);
    Test_Real_Powered_Series.random_series(size,bcf,bpw);
    put_line("The sum of ");
    Test_Real_Powered_Series.write(acf,apw);
    put_line("and");
    Test_Real_Powered_Series.write(bcf,bpw);
    put_line("is");
    Add(acf,bcf,apw,bpw,sumcf,sumpw);
    Test_Real_Powered_Series.write(sumcf,sumpw);
    Sub(sumcf,bcf,sumpw,bpw,difcf,difpw);
    put_line("After subtracting the second series from the sum :");
    Test_Real_Powered_Series.write(difcf(0..size),difpw(1..size));
    put_line("the first series :");
    Test_Real_Powered_Series.write(acf,apw);
    equ := Equal(difcf(0..size),acf,difpw(1..size),apw);
    put("Component wise equal ? ");
    if equ
     then put_line("Yes.");
     else put_line("No, bug!?");
    end if;
  end Test_Addition;

  procedure Test_Multiplication ( size : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random series of the given size,
  --   tests their multiplication and the division.

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    prdsize : constant integer32 := (size+1)*(size+1) - 1;
    prdcf : Standard_Complex_Vectors.Vector(0..prdsize);
    prdpw : Standard_Floating_Vectors.Vector(1..prdsize);

  begin
    put("-> multiplying two series of size "); put(size,1); put_line(" ...");
    Test_Real_Powered_Series.random_series(size,acf,apw);
    Test_Real_Powered_Series.random_series(size,bcf,bpw);
    put_line("The product of ");
    Test_Real_Powered_Series.write(acf,apw);
    put_line("and");
    Test_Real_Powered_Series.write(bcf,bpw);
    put_line("is");
    Mul(acf,bcf,apw,bpw,prdcf,prdpw);
    Test_Real_Powered_Series.write(prdcf,prdpw);
  end Test_Multiplication;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts to select a test, the size of a real power series,
  --   and then runs the test.

    ans : character;
    size : integer32 := 0;

  begin
    new_line;
    put_line("MENU to test real powered series arithmetic ...");
    put_line("  0. generate a random series");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication and division");
    put("Type 0, 1, or 2 to select a test : ");
    Communications_with_User.Ask_Alternative(ans,"012");
    new_line;
    put("Give the size of the series : "); get(size);
    case ans is
      when '0' => Test_Random_Series(size);
      when '1' => Test_Addition(size);
      when '2' => Test_Multiplication(size);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hunroots;
