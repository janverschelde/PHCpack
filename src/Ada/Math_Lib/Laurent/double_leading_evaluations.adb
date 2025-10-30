with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;

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
    val := Leading_Power(deg(deg'first).all,pwr,0,vrblvl-1);
    if vrblvl > 0 then
      put("leading power at i = "); put(deg'first,1); put(" : ");
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

end Double_Leading_Evaluations;
