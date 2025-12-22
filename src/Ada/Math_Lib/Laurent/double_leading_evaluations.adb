with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Floating_Vectors_IO;      use Standard_Floating_Vectors_IO;
with Double_Real_Powered_Series;

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
    Double_Real_Powered_Series.Sort(val);
  end Evaluate_Powers;

  function Leading_Coefficient
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               difidx : integer32 := 0; vrblvl : integer32 := 0 )
             return Complex_Number is

    res : Complex_Number := create(1.0);

  begin
    if vrblvl > 0 then
      put("-> in Double_Leading_Evaluations.leading_coefficient, difidx : ");
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

  function Second_Derivative
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               i : integer32; vrblvl : integer32 := 0 )
             return Complex_Number is

    res : Complex_Number := create(1.0);

  begin
    if vrblvl > 0 then
      put("-> in Double_Leading_Evaluations.second_derivative, i : ");
      put(i,1); put_line(" ...");
    end if;
    if deg(i) = 0 then
      return create(0.0);
    elsif deg(i) = 1 then
      return create(0.0);
    end if;
    for k in deg'range loop
      if deg(k) > 0 then
        if k /= i then
          for j in 1..deg(k) loop
            res := res*cff(k);
          end loop;
        elsif deg(k) > 1 then
          for j in 1..(deg(k)-2) loop
            res := res*cff(k);
          end loop;
          res := double_float(deg(k)*(deg(k)-1))*res;
        end if;
      elsif deg(k) < 0 then
        for j in 1..(-deg(k)) loop
          res := res/cff(k);
        end loop;
        if k = i then
          res := res/cff(k);
          res := res/cff(k);
          res := double_float(deg(k)*(deg(k)-1))*res;
        end if;
      end if;
    end loop;
    return res;
  end Second_Derivative;

  function Second_Mixed_Derivative
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               i,j : integer32; vrblvl : integer32 := 0 )
             return Complex_Number is

    res : Complex_Number := create(1.0);

  begin
    if vrblvl > 0 then
      put("-> in Double_Leading_Evaluations.second_mixed_derivative, i : ");
      put(i,1); put(", j : "); put(j,1); put_line(" ...");
    end if;
    if deg(i) = 0 then
      return create(0.0);
    elsif deg(j) = 0 then
      return create(0.0);
    end if;
    for k in deg'range loop
      if k = i or k = j then
        if deg(k) > 1 then
          for kk in 1..(deg(k)-1) loop   
            res := res*cff(k);
          end loop;
        elsif deg(k) < 0 then
          for kk in 1..(-deg(k)+1) loop
            res := res/cff(k);
          end loop;
        end if;
        res := double_float(deg(k))*res;
      else
        if deg(k) > 0 then
          for kk in 1..deg(k) loop
            res := res*cff(k);
          end loop;
        elsif deg(k) < 0 then
          for kk in 1..(-deg(k)) loop
            res := res/cff(k);
          end loop;
        end if;
      end if;
    end loop;
    return res;
  end Second_Mixed_Derivative;

  function Third_Derivative
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               i : integer32; vrblvl : integer32 := 0 )
             return Complex_Number is

    res : Complex_Number := create(1.0);

  begin
    if vrblvl > 0 then
      put("-> in Double_Leading_Evaluations.third_derivative, i : ");
      put(i,1); put_line(" ...");
    end if;
    if deg(i) = 0 then
      return create(0.0);
    elsif deg(i) = 1 then
      return create(0.0);
    elsif deg(i) = 2 then
      return create(0.0);
    end if;
    for k in deg'range loop
      if deg(k) > 0 then
        if k /= i then
          for j in 1..deg(k) loop
            res := res*cff(k);
          end loop;
        elsif deg(k) > 2 then
          for j in 1..(deg(k)-3) loop
            res := res*cff(k);
          end loop;
          res := double_float(deg(k)*(deg(k)-1)*(deg(k)-2))*res;
        end if;
      elsif deg(k) < 0 then
        for j in 1..(-deg(k)) loop
          res := res/cff(k);
        end loop;
        if k = i then
          res := res/cff(k);
          res := res/cff(k);
          res := res/cff(k);
          res := double_float(deg(k)*(deg(k)-1)*(deg(k)-2))*res;
        end if;
      end if;
    end loop;
    return res;
  end Third_Derivative;

  function Third_Semi_Mixed_Derivative
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               i,j : integer32; vrblvl : integer32 := 0 )
             return Complex_Number is

    res : Complex_Number := create(1.0);

  begin
    if vrblvl > 0 then
      put("-> in Double_Leading_Evaluations.");
      put_line("third_semi_mixed_derivative ...");
      put("  i : "); put(i,1);
      put(", j : "); put(j,1); new_line;
    end if;
    if deg(i) = 0 then
      return create(0.0);
    elsif deg(i) = 1 then
      return create(0.0);
    elsif deg(j) = 0 then
      return create(0.0);
    end if;
    for k in deg'range loop
      if k = i then -- take second derivative w.r.t. i
        if deg(k) > 1 then
          for kk in 1..(deg(k)-2) loop
            res := res*cff(k);
          end loop;
        elsif deg(k) < 0 then
          for kk in 1..(-deg(k)+2) loop
            res := res/cff(k);
          end loop;
        end if;
        res := double_float(deg(k)*(deg(k)-1))*res;
      elsif k = j then -- take first derivative w.r.t. j
        if deg(k) > 1 then
          for kk in 1..(deg(k)-1) loop
            res := res*cff(k);
          end loop;
        else
          for kk in 1..(-deg(k)+1) loop
            res := res/cff(k);
          end loop;
        end if;
        res := double_float(deg(k))*res;
      else -- just evaluate 
        if deg(k) > 0 then
          for kk in 1..deg(k) loop
            res := res*cff(k);
          end loop;
        elsif deg(k) < 0 then
          for kk in 1..(-deg(k)) loop
            res := res/cff(k);
          end loop;
        end if;
      end if;
    end loop;
    return res;
  end Third_Semi_Mixed_Derivative;

  function Third_Fully_Mixed_Derivative
             ( deg : Standard_Integer_Vectors.Vector;
               cff : Standard_Complex_Vectors.Vector;
               i,j,k : integer32; vrblvl : integer32 := 0 )
             return Complex_Number is

    res : Complex_Number := create(1.0);

  begin
    if vrblvl > 0 then
      put("-> in Double_Leading_Evaluations.");
      put_line("third_fully_mixed_derivative ...");
      put("  i : "); put(i,1);
      put(", j : "); put(j,1);
      put(", k : "); put(k,1); new_line;
    end if;
    if deg(i) = 0 then
      return create(0.0);
    elsif deg(j) = 0 then
      return create(0.0);
    elsif deg(k) = 0 then
      return create(0.0);
    end if;
    for L in deg'range loop
      if L = i or L = j or L = k then
        if deg(L) > 1 then
          for kk in 1..(deg(L)-1) loop   
            res := res*cff(L);
          end loop;
        elsif deg(L) < 0 then
          for kk in 1..(-deg(L)+1) loop
            res := res/cff(L);
          end loop;
        end if;
        res := double_float(deg(L))*res;
      else
        if deg(L) > 0 then
          for kk in 1..deg(L) loop
            res := res*cff(L);
          end loop;
        elsif deg(L) < 0 then
          for kk in 1..(-deg(L)) loop
            res := res/cff(L);
          end loop;
        end if;
      end if;
    end loop;
    return res;
  end Third_Fully_Mixed_Derivative;

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
    Double_Real_Powered_Series.Sort(ydg,ycf);
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
