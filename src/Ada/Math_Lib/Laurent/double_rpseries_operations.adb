with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Numbers;
with Standard_Random_Vectors;

package body Double_rpSeries_Operations is

  procedure Sort ( x : in out Standard_Floating_Vectors.Vector ) is

    val : double_float;

  begin
    for i in x'range loop
      for j in i+1..x'last loop
        val := x(j);
        if val < x(i) then -- x(j) is the new minimum
          x(j) := x(i);    -- swap x(i) and x(j)
          x(i) := val;     -- x(i) is the minimum
        end if;
      end loop;
    end loop;
  end Sort;

  procedure Sort ( x : in out Standard_Floating_Vectors.Vector;
                   y : in out Standard_Complex_Vectors.Vector ) is

    val : double_float;
    tmp : Complex_Number;

  begin
    for i in x'range loop
      for j in i+1..x'last loop
        val := x(j);
        if val < x(i) then -- x(j) is the new minimum
          x(j) := x(i);    -- swap x(i) and x(j)
          x(i) := val;     -- x(i) is the minimum
          tmp := y(i);     -- swap y(i) and y(j)
          y(i) := y(j);
          y(j) := tmp;
        end if;
      end loop;
    end loop;
  end Sort;

  procedure Shift_Zeros ( cff : in out Standard_Complex_Vectors.Vector;
                          pwt : in out Standard_Floating_Vectors.Vector;
                          tol : in double_float := 1.0E-12 ) is 

    idx : integer32 := pwt'first;
    nzidx : integer32;

  begin
    while idx <= pwt'last loop
      if AbsVal(cff(idx)) > tol then
         nzidx := idx + 1;
      else -- if AbsVal(cff(idx)) < tol then
        if idx < pwt'last then
          nzidx := idx + 1; -- index of next nonzero coefficient
          while AbsVal(cff(nzidx)) < tol loop 
            nzidx := nzidx + 1;
            exit when (nzidx > pwt'last);
          end loop;
          if nzidx <= pwt'last then
            cff(idx) := cff(nzidx);
            pwt(idx) := pwt(nzidx);
            cff(nzidx) := create(0.0);
          end if;
        end if;
      end if;
      exit when (nzidx > pwt'last); -- no more nonzero coefficients
      idx := idx + 1;
    end loop;
  end Shift_Zeros;

  procedure Normalize ( cff : in out Standard_Complex_Vectors.Vector;
                        pwt : in out Standard_Floating_Vectors.Vector;
                        tol : in double_float := 1.0E-12 ) is 

    dif : double_float;

  begin
    for i in pwt'first..pwt'last-1 loop
      for j in i+1..pwt'last loop
        dif := abs(pwt(i) - pwt(j));
        exit when (dif > tol);
        cff(i) := cff(i) + cff(j);
        cff(j) := create(0.0);
      end loop;
    end loop;
    Shift_Zeros(cff,pwt,tol);
  end Normalize;

  procedure Normalize ( cf : in Standard_Complex_VecVecs.VecVec;
                        dg : in Standard_Floating_VecVecs.VecVec;
                        tol : in double_float := 1.0E-12 ) is
  begin
    for i in cf'range loop
      Normalize(cf(i).all,dg(i).all,tol);
    end loop;
  end Normalize;

  function Equal ( acf,bcf : Standard_Complex_Vectors.Vector;
                   apw,bpw : Standard_Floating_Vectors.Vector;
                   tol : double_float := 1.0E-12 ) return boolean is

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

-- BASIC ARITHMETIC :

  procedure Add ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

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
    Normalize(cff,pwt);
  end Add;

  procedure Sub ( acf,bcf : in Standard_Complex_Vectors.Vector;
                  apw,bpw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

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
    Sort(pwt,cff);
    Normalize(cff,pwt);
  end Mul;

  procedure Inv ( acf : in Standard_Complex_Vectors.Vector;
                  apw : in Standard_Floating_Vectors.Vector;
                  cff : out Standard_Complex_Vectors.Vector;
                  pwt : out Standard_Floating_Vectors.Vector ) is

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

-- USEFUL FUNCTIONS :

  function Positive_Minimum_Index
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector;
               tol : double_float := 1.0E-12 ) return integer32 is

    res,idx : integer32;
    psm : double_float;

  begin
    for i in v'range loop -- find first positive number
      if AbsVal(c(i)) > tol then
        if v(i) > tol then
          idx := i;
          psm := v(i);
          exit;
        end if;
      end if;
    end loop;
    res := idx;
    for i in idx+1..v'last loop
      if AbsVal(c(i)) > tol then
        if v(i) > tol and then v(i) < psm
         then psm := v(i); res := i;
        end if;
      end if;
    end loop;
    return res;
  end Positive_Minimum_Index;

  function Positive_Minimum
             ( v : Standard_Floating_Vectors.Vector;
               tol : double_float := 1.0E-12 ) return double_float is

    res : double_float;
    idx : integer32;

  begin
    for i in v'range loop -- find first positive number
      if v(i) > tol then
        idx := i;
        res := v(i);
        exit;
      end if;
    end loop;
    for i in idx+1..v'last loop
      if v(i) > tol and then v(i) < res
       then res := v(i);
      end if;
    end loop;
    return res;
  end Positive_Minimum;

  function Positive_Minimum
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector;
               tol : double_float := 1.0E-12 ) return double_float is

    res : double_float;
    idx : integer32;

  begin
    for i in v'range loop -- find first positive number
      if AbsVal(c(i)) > tol then
        if v(i) > tol then
          idx := i;
          res := v(i);
          exit;
        end if;
      end if;
    end loop;
    for i in idx+1..v'last loop
      if AbsVal(c(i)) > tol then
        if v(i) > tol and then v(i) < res
         then res := v(i);
        end if;
      end if;
    end loop;
    return res;
  end Positive_Minimum;

  function Positive_Minima
             ( c : Standard_Complex_VecVecs.VecVec;
               v : Standard_Floating_VecVecs.VecVec;
               tol : double_float := 1.0E-12 )
             return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := Positive_Minimum(c(i).all,v(i).all,tol);
    end loop;
    return res;
  end Positive_Minima;

  function Coefficient ( c : Standard_Complex_Vectors.Vector;
                         e : Standard_Floating_Vectors.Vector;
                         p : double_float; tol : double_float := 1.0E-12 )
                       return Complex_Number is
  begin
    for i in e'range loop
      if abs(e(i) - p) < tol
       then return c(i);
      end if;
    end loop;
    return create(0.0);
  end Coefficient;

  function Coefficients ( c : Standard_Complex_VecVecs.VecVec;
                          e : Standard_Floating_VecVecs.VecVec;
                          p : Standard_Floating_Vectors.Vector;
                          tol : double_float := 1.0E-12 )
                       return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(p'range);

  begin
    for i in res'range loop
      res(i) := Coefficient(c(i).all,e(i).all,p(i),tol);
    end loop;
    return res;
  end Coefficients;

  function Random_Leading_Powers
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..dim);
    rnd : double_float;

  begin
    for i in 1..dim loop
      rnd := Standard_Random_Numbers.Random;
      res(i) := abs(rnd);
      while res(i) < 1.0 loop -- make larger than 1.0
        res(i) := res(i) + 0.1;
      end loop;
    end loop;
    return res;
  end Random_Leading_Powers;

  procedure Random_Power_Series
              ( dim : in integer32;
                nbt : in Standard_Integer_Vectors.Vector;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec ) is

    rnd : double_float;

  begin
    for i in 1..dim loop
      declare
        cfi : constant Standard_Complex_Vectors.Vector(0..nbt(i))
            := Standard_Random_Vectors.Random_Vector(0,nbt(i));
        pwi : Standard_Floating_Vectors.Vector(1..nbt(i));
      begin
        cff(i) := new Standard_Complex_Vectors.Vector'(cfi);
        rnd := Standard_Random_Numbers.Random;
        pwi(1) := abs(rnd);
        while pwi(1) < 1.0 loop  -- make larger than 1.0
          pwi(1) := pwi(1) + 0.1;
        end loop;
        for j in 2..nbt(i) loop
          rnd := abs(Standard_Random_Numbers.Random); -- rnd in [0,1]
          pwi(j) := pwi(j-1) + rnd;
          while pwi(j) > 2.0*pwi(j-1) loop -- make pwi(j) < 2*pwi(j-1)
            rnd := rnd/2.0;
            pwi(j) := pwi(j-1) + rnd;
          end loop;
        end loop;
        pwr(i) := new Standard_Floating_Vectors.Vector'(pwi);
      end;
    end loop;
  end Random_Power_Series;

  function Evaluate_Series
             ( cff : Standard_Complex_VecVecs.VecVec;
               pwr : Standard_Floating_VecVecs.VecVec; tpt : double_float )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(cff'range);

  begin
    for i in res'range loop
      res(i) := cff(i)(0);
      for j in 1..cff(i)'last loop
        res(i) := res(i) + cff(i)(j)*(tpt**pwr(i)(j));
      end loop;
    end loop;
    return res;
  end Evaluate_Series;

end Double_rpSeries_Operations;
