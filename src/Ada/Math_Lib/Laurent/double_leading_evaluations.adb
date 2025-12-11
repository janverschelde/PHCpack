with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Floating_Vectors_IO;      use Standard_Floating_Vectors_IO;

package body Double_Leading_Evaluations is

  function Leading_Power 
             ( deg : Standard_Integer_Vectors.Vector;
               pwr : Standard_Floating_Vectors.Vector;
               difidx : integer32 := 0; vrblvl : integer32 := 0 )
             return double_float is

    res : double_float := 0.0;

  begin
    if vrblvl > 0 then
      put("-> in Double_Leading_Evaluations.leading_power 1, difidx : ");
      put(difidx,1); put_line(" ...");
    end if;
    for i in deg'range loop
      if deg(i) /= 0 then
        if difidx = i
         then res := res + double_float(deg(i)-1)*pwr(i);
         else res := res + double_float(deg(i))*pwr(i);
        end if;
      end if;
    end loop;
    return res;
  end Leading_Power;

  procedure Leading_Power 
              ( deg : in Standard_Integer_VecVecs.VecVec;
                pwr : in Standard_Floating_Vectors.Vector;
                val : out double_float; idx : out integer32;
                vrblvl : in integer32 := 0 ) is

    lpr : double_float;

  begin
    if vrblvl > 0 then
      put_line("-> in Double_Leading_Evaluations.leading_power 2 ...");
    end if;
    idx := deg'first;
    val := Leading_Power(deg(idx).all,pwr,0,vrblvl-1);
    if vrblvl > 0 then
      put("leading power at i = "); put(idx,1); put(" : ");
      put(val); new_line;
    end if;
    for i in deg'first+1..deg'last loop
      lpr := Leading_Power(deg(i).all,pwr,0,vrblvl-1);
      if vrblvl > 0 then
        put("leading power at i = "); put(i,1); put(" : ");
        put(lpr); new_line;
      end if;
      if lpr < val
       then val := lpr; idx := i;
      end if;
    end loop;
  end Leading_Power;

  procedure Sort ( x : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Sorts the numbers in x in increasing order.

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

  -- DESCRIPTION :
  --   Sorts the numbers in x in increasing order
  --   and swaps the corresponding numbers in y accordingly.

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

  procedure Evaluate_Powers
              ( deg : in Standard_Integer_VecVecs.VecVec;
                pwr : in Standard_Floating_Vectors.Vector;
                val : out Standard_Floating_Vectors.Vector;
                idx : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in Double_Leading_Evaluations.evaluate_powers ...");
    end if;
    idx := deg'first;
    val(idx) := Leading_Power(deg(idx).all,pwr,0,vrblvl-1);
    if vrblvl > 0 then
      put("leading power at i = "); put(idx,1); put(" : ");
      put(val(idx)); new_line;
    end if;
    for i in deg'first+1..deg'last loop
      val(i) := Leading_Power(deg(i).all,pwr,0,vrblvl-1);
      if vrblvl > 0 then
        put("leading power at i = "); put(i,1); put(" : ");
        put(val(i)); new_line;
      end if;
      if val(i) < val(idx)
       then idx := i;
      end if;
    end loop;
    Sort(val);
  end Evaluate_Powers;

  function Leading_Coefficient
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               difidx : integer32 := 0; vrblvl : integer32 := 0 )
             return Complex_Number is

    res : Complex_Number := create(1.0);

  begin
    if vrblvl > 0 then
      put("-> in Double_Leaving_Evaluations.leading_coefficient, difidx : ");
      put(difidx,1); put_line(" ...");
    end if;
    if difidx > 0 then
      if deg(difidx) = 0
       then return create(0.0);
      end if;
    end if;
    for i in deg'range loop
      if deg(i) > 0 then
        if difidx /= i then
          for j in 1..deg(i) loop
            res := res*cff(i);
          end loop;
        elsif deg(i) > 1 then
          for j in 1..(deg(i)-1) loop
            res := res*cff(i);
          end loop;
          res := double_float(deg(i))*res;
        end if;
      elsif deg(i) < 0 then
        for j in 1..(-deg(i)) loop
          res := res/cff(i);
        end loop;
        if difidx = i then
          res := res/cff(i);
          res := double_float(deg(i))*res;
        end if;
      end if;
    end loop;
    return res;
  end Leading_Coefficient;

  procedure Evaluate_Polynomial
              ( pcf : in Standard_Complex_Vectors.Vector;
                pct : in Standard_Floating_Vectors.Vector;
                pdg : in Standard_Integer_VecVecs.VecVec;
                xcf : in Standard_Complex_Vectors.Vector;
                xdg : in Standard_Floating_Vectors.Vector;
                ycf : out Standard_Complex_Vectors.Vector;
                ydg : out Standard_Floating_Vectors.Vector;
                idx : out integer32; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in Double_Leading_Evaluations.evaluate_polynomial ...");
    end if;
    idx := pdg'first;
    ydg(idx) := pct(idx) + Leading_Power(pdg(idx).all,xdg,0,vrblvl-1);
    ycf(idx) := pcf(idx)*Leading_Coefficient(pdg(idx).all,xcf,0,vrblvl-1);
    if vrblvl > 0 then
      put("leading power at i = "); put(idx,1); put(" : ");
      put(ydg(idx)); new_line;
    end if;
    for i in pdg'first+1..pdg'last loop
      ydg(i) := Leading_Power(pdg(i).all,xdg,0,vrblvl-1);
      ycf(i) := pcf(i)*Leading_Coefficient(pdg(i).all,xcf,0,vrblvl-1);
      if vrblvl > 0 then
        put("leading power at i = "); put(i,1); put(" : ");
        put(ydg(i)); new_line;
      end if;
      if ydg(i) < ydg(idx)
       then idx := i;
      end if;
    end loop;
    sort(ydg,ycf);
  end Evaluate_Polynomial;

  procedure Evaluate_System
              ( pcf : in Standard_Complex_VecVecs.VecVec;
                pct : in Standard_Floating_VecVecs.VecVec;
                pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                xcf : in Standard_Complex_Vectors.Vector;
                xdg : in Standard_Floating_Vectors.Vector;
                ycf : in Standard_Complex_VecVecs.VecVec;
                ydg : in Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    idx : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in Double_Leading_Evaluations.evaluate_system ...");
    end if;
    for i in pcf'range loop
      if vrblvl > 0 then
        put("Evaluating polynomial "); put(i,1); put_line(" ...");
      end if;
      Double_Leading_Evaluations.Evaluate_Polynomial
        (pcf(i).all,pct(i).all,pdg(i).all,xcf,xdg,ycf(i).all,ydg(i).all,idx,
         vrblvl-1);
      if vrblvl > 0 then
        put_line("evaluated powers :"); put_line(ydg(i).all);
        put("power value : "); put(ydg(i)(ydg(i)'first));
        put(" at index "); put(idx,1); new_line;
      end if;
    end loop;
  end Evaluate_System;

end Double_Leading_Evaluations;
