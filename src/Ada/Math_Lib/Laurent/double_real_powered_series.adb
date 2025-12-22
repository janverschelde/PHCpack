with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Random_Numbers;
with Standard_Random_Vectors;

package body Double_Real_Powered_Series is

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

  procedure Normalize ( cf : in out Standard_Complex_Vectors.Vector;
                        dg : in Standard_Floating_Vectors.Vector ) is

    tol : constant double_float := 1.0E-12;
    dif : double_float;

  begin
    for i in cf'first..cf'last-1 loop
      for j in i+1..cf'last loop
        dif := abs(dg(i) - dg(j));
        exit when (dif > tol);
        cf(i) := cf(i) + cf(j);
        cf(j) := create(0.0);
      end loop;
    end loop;
  end Normalize;

  function Positive_Minimum_Index
             ( c : Standard_Complex_Vectors.Vector;
               v : Standard_Floating_Vectors.Vector ) return integer32 is

    res,idx : integer32;
    psm : double_float;
    tol : constant double_float := 1.0E-12;

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

  function Random_Leading_Powers
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..dim);
    rnd : double_float;

  begin
    for i in 1..dim loop
      rnd := Standard_Random_Numbers.Random;
      res(i) := abs(rnd);
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

end Double_Real_Powered_Series;
