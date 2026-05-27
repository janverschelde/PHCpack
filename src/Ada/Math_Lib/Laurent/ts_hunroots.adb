with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Double_Real_Powered_Series;
with Test_Real_Powered_Series;

procedure ts_hunroots is

-- DESCRIPTION :
--   Tests the square root computation of real power series.

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

  procedure Inv ( acf : in Standard_Complex_Vectors.Vector;
                  apw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the multiplicative inverse of the
  --   series given in (acf, apw).

  -- REQUIRED : acf(0) /= 0.

    sqrdiv : Complex_Number;
    minidx : Standard_Integer_Vectors.Vector(apw'range) := (apw'range => 1);
    apwidx : integer32 := 1;
    sumpwr : double_float;
    cnvsum : Complex_Number;

  begin
    cff(0) := 1.0/acf(0);               -- constant term
    sqrdiv := cff(0)/acf(0);
    pwt(1) := apw(1);                   -- first order term
    cff(1) := -acf(1)*sqrdiv;
    apwidx := 2;                        -- apw(1) term canceled
    if pwt'last > 1 then
      if apw(2) < 2.0*apw(1) then       -- second order term
        pwt(2) := apw(2);
        cff(2) := -acf(2)*sqrdiv;
        apwidx := 3;                    -- apw(2) term canceled
      else
        pwt(2) := 2.0*apw(1);
        minidx(1) := 2;                 -- apw(1) + pwt(1) canceled
        if apw(2) > 2.0*apw(1) then
          cff(2) := -acf(1)*cff(1)/acf(0);
        else -- apw(2) = 2.0*apw(1)
          cff(2) := -(acf(1)*cff(1) + acf(2)*cff(0))/acf(0);
          apwidx := 3;                  -- apw(2) term canceled
        end if;
      end if;
      for k in 3..apw'last loop
        exit when ((k > pwt'last) or (k > cff'last));
        pwt(k) := apw(k);
        for i in 1..k-1 loop
          sumpwr := apw(i) + pwt(minidx(i));
          if sumpwr < pwt(k)
           then pwt(k) := sumpwr;
          end if;
        end loop;
        if pwt(k) = apw(k) then
          apwidx := apwidx + 1;
          cnvsum := cff(0)*acf(k);
        else
          cnvsum := create(0.0);
        end if;
        for i in 1..k-1 loop
          sumpwr := apw(i) + pwt(minidx(i));
          if sumpwr = pwt(k) then
            cnvsum := cnvsum + acf(i)*cff(minidx(i));
            minidx(i) := minidx(i) + 1;
          end if;
        end loop;
        cff(k) := -cnvsum/acf(0);
      end loop;
    end if;
  end Inv;

  procedure Div ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Returns in (cff, pwt) the series (acf, apw) divided by (bcf, bpw),
  --   via multiplication of the inverse of (bcf, bpw).

  -- REQUIRED : bcf(0) /= 0.

    size : constant integer32 := bcf'last;
    invbcf : Standard_Complex_Vectors.Vector(0..size);
    invbpw : Standard_Floating_Vectors.Vector(1..size);
    prdsize : constant integer32 := (size+1)*(size+1) - 1;
    prdcf : Standard_Complex_Vectors.Vector(0..prdsize);
    prdpw : Standard_Floating_Vectors.Vector(1..prdsize);

  begin
    Inv(bcf,bpw,invbcf,invbpw);
    Mul(acf,invbcf,apw,invbpw,prdcf,prdpw);
    cff(0) := prdcf(0);
    for i in pwt'range loop
      cff(i) := prdcf(i);
      pwt(i) := prdpw(i);
    end loop;
  end Div;

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
  --   tests their multiplication.

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

  procedure Test_Inverse ( size : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random series of the given size,
  --   computes the inverse and multiplies the random series
  --   with the inverse to verify the correctness.

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    prdsize : constant integer32 := (size+1)*(size+1) - 1;
    prdcf : Standard_Complex_Vectors.Vector(0..prdsize);
    prdpw : Standard_Floating_Vectors.Vector(1..prdsize);

  begin
    Test_Real_Powered_Series.random_series(size,acf,apw);
    for i in apw'range loop
      apw(i) := apw(i) + 1.0;
    end loop;
    put("-> the inverse of a random series of size "); put(size,1);
    put_line(" :");
    Test_Real_Powered_Series.write(acf,apw);
    put_line("is");
    Inv(acf,apw,bcf,bpw);
    Test_Real_Powered_Series.write(bcf,bpw);
    put_line("the product of a random series with its inverse :");
    Mul(acf,bcf,apw,bpw,prdcf,prdpw);
    Test_Real_Powered_Series.write(prdcf(0..size),prdpw(1..size));
  end Test_Inverse;

  procedure Test_Division ( size : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two random series of the given size,
  --   tests their division.

    acf,bcf : Standard_Complex_Vectors.Vector(0..size);
    apw,bpw : Standard_Floating_Vectors.Vector(1..size);
    qsize : constant integer32 := (size+1)*(size+1) - 1;
    qcf : Standard_Complex_Vectors.Vector(0..qsize);
    qpw : Standard_Floating_Vectors.Vector(1..qsize);
    psize : constant integer32 := (size+1)*(qsize+1) - 1;
    pcf : Standard_Complex_Vectors.Vector(0..psize);
    ppw : Standard_Floating_Vectors.Vector(1..psize);
    equ : boolean;

  begin
    put("-> dividing two series of size "); put(size,1); put_line(" ...");
    Test_Real_Powered_Series.random_series(size,acf,apw);
    Test_Real_Powered_Series.random_series(size,bcf,bpw);
    for i in apw'range loop
      apw(i) := apw(i) + 1.0;
      bpw(i) := bpw(i) + 1.0;
    end loop;
    put_line("The division of ");
    Test_Real_Powered_Series.write(acf,apw);
    put_line("by");
    Test_Real_Powered_Series.write(bcf,bpw);
    put_line("is");
    Div(acf,bcf,apw,bpw,qcf,qpw);
    Test_Real_Powered_Series.write(qcf(0..size),qpw(1..size));
    put_line("multiplying the quotient with the second series ...");
    Mul(bcf,qcf,bpw,qpw,pcf,ppw);
    Test_Real_Powered_Series.write(pcf(0..size),ppw(1..size));
    put_line("The first series :");
    Test_Real_Powered_Series.write(acf,apw);
    equ := Equal(pcf(0..size),acf,ppw(1..size),apw);
    put("Component wise equal ? ");
    if equ
     then put_line("Yes.");
     else put_line("No, bug!?");
    end if;
  end Test_Division;

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
    put_line("  2. test multiplication");
    put_line("  3. test multiplicative inverse");
    put_line("  4. test division");
    put("Type 0, 1, 2, 3, or 4 to select a test : ");
    Communications_with_User.Ask_Alternative(ans,"01234");
    new_line;
    put("Give the size of the series : "); get(size);
    case ans is
      when '0' => Test_Random_Series(size);
      when '1' => Test_Addition(size);
      when '2' => Test_Multiplication(size);
      when '3' => Test_Inverse(size);
      when '4' => Test_Division(size);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_hunroots;
