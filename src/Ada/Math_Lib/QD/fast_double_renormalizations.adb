with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Basics;

package body Fast_Double_Renormalizations is

  procedure fast_renorm
              ( x0,x1,x2,x3 : in double_float;
                r0,r1,r2 : out double_float ) is

    f0,f1,f2,f3,pr : double_float;
    ptr : integer32;

  begin
    Double_Double_Basics.quick_two_sum(x2,x3,pr,f3);  
    Double_Double_Basics.quick_two_sum(x1,pr,pr,f2);  
    Double_Double_Basics.quick_two_sum(x0,pr,f0,f1);  
    if f1 = 0.0 then
      pr := f0;
      ptr := 0;
      Double_Double_Basics.quick_two_sum(pr,f2,r0,pr);
    else
      r0 := f0;
      pr := f1;
      ptr := 1;
      Double_Double_Basics.quick_two_sum(pr,f2,r1,pr);
    end if;
    if pr = 0.0 then
      if ptr = 0
       then pr := r0;
       else pr := r1;
      end if;
    else
      ptr := ptr + 1;
    end if;
    if ptr = 0 then
      Double_Double_Basics.quick_two_sum(pr,f3,r0,pr);
    elsif ptr = 1 then
      Double_Double_Basics.quick_two_sum(pr,f3,r1,pr);
    else
      Double_Double_Basics.quick_two_sum(pr,f3,r2,pr);
    end if;
    if pr = 0.0 then
      if ptr = 0 then
        pr := r0;
      elsif ptr = 1 then
        pr := r1;
      else
        pr := r2;
      end if;
    else
      ptr := ptr + 1;
    end if;
    if ptr < 3 and pr /= 0.0 then
      if ptr = 0 then
        r0 := pr;
      elsif ptr = 1 then
        r1 := pr;
      else
        r2 := pr;
      end if;
      ptr := ptr + 1;
    end if;
    if ptr < 1 then
      r0 := 0.0; r1 := 0.0; r2 := 0.0;
    elsif ptr < 2 then
      r1 := 0.0; r2 := 0.0;
    elsif ptr < 3 then
      r2 := 0.0;
    end if;
  end fast_renorm;

end Fast_Double_Renormalizations;
