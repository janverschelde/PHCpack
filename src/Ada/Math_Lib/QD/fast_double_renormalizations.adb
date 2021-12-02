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

  procedure renorm_add1
              ( x0,x1,x2 : in double_float; y : in double_float;
                r0,r1,r2 : out double_float ) is

    f0,f1,f2,f3,pr : double_float;
    ptr : integer32;

  begin
    Double_Double_Basics.two_sum(x2,y,pr,f3);
    Double_Double_Basics.two_sum(x1,pr,pr,f2);
    Double_Double_Basics.two_sum(x0,pr,f0,f1);
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
  end renorm_add1;

  procedure renorm5
              ( f0,f1,f2,f3,f4,f5 : in double_float;
                pr : in out double_float;
                r0,r1,r2,r3,r4 : out double_float ) is

    ptr : integer;

  begin
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
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f3,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f3,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f3,r2,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f4,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f4,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f4,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f4,r3,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f5,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f5,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f5,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f5,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f5,r4,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    if ptr < 5 and pr /= 0.0 then
      case ptr is
        when 0 => r0 := pr;
        when 1 => r1 := pr;
        when 2 => r2 := pr;
        when 3 => r3 := pr;
        when 4 => r4 := pr;
        when others => null;
      end case;
      ptr := ptr + 1;
    end if;
    if ptr < 1 then
      r4 := 0.0; r3 := 0.0; r2 := 0.0; r1 := 0.0; r0 := 0.0;
    elsif ptr < 2 then
      r4 := 0.0; r3 := 0.0; r2 := 0.0; r1 := 0.0;
    elsif ptr < 3 then
      r4 := 0.0; r3 := 0.0; r2 := 0.0;
    elsif ptr < 4 then
      r4 := 0.0; r3 := 0.0;
    elsif ptr < 5 then
      r4 := 0.0;
    end if;
  end renorm5;

  procedure fast_renorm
              ( x0,x1,x2,x3,x4,x5 : in double_float;
                r0,r1,r2,r3,r4 : out double_float ) is

    f0,f1,f2,f3,f4,f5,pr : double_float;

  begin
    Double_Double_Basics.quick_two_sum(x4,x5,pr,f5);
    Double_Double_Basics.quick_two_sum(x3,pr,pr,f4);
    Double_Double_Basics.quick_two_sum(x2,pr,pr,f3);
    Double_Double_Basics.quick_two_sum(x1,pr,pr,f2);
    Double_Double_Basics.quick_two_sum(x0,pr,f0,f1);
    renorm5(f0,f1,f2,f3,f4,f5,pr,r0,r1,r2,r3,r4);
  end fast_renorm;

  procedure renorm_add1
              ( x0,x1,x2,x3,x4 : in double_float;
                y : in double_float;
                r0,r1,r2,r3,r4 : out double_float ) is

    f0,f1,f2,f3,f4,f5,pr : double_float;

  begin
    Double_Double_Basics.two_sum(x4,y,pr,f5);
    Double_Double_Basics.two_sum(x3,pr,pr,f4);
    Double_Double_Basics.two_sum(x2,pr,pr,f3);
    Double_Double_Basics.two_sum(x1,pr,pr,f2);
    Double_Double_Basics.two_sum(x0,pr,f0,f1);
    renorm5(f0,f1,f2,f3,f4,f5,pr,r0,r1,r2,r3,r4);
  end renorm_add1;

  procedure renorm8
              ( f0,f1,f2,f3,f4,f5,f6,f7,f8 : in double_float;
                pr : in out double_float;
                r0,r1,r2,r3,r4,r5,r6,r7 : out double_float ) is

    ptr : integer;

  begin
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
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f3,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f3,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f3,r2,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f4,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f4,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f4,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f4,r3,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f5,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f5,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f5,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f5,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f5,r4,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f6,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f6,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f6,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f6,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f6,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f6,r5,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f7,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f7,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f7,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f7,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f7,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f7,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f7,r6,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f8,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f8,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f8,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f8,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f8,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f8,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f8,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f8,r7,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    if ptr < 8 and pr /= 0.0 then
      case ptr is
        when 0 => r0 := pr;
        when 1 => r1 := pr;
        when 2 => r2 := pr;
        when 3 => r3 := pr;
        when 4 => r4 := pr;
        when 5 => r5 := pr;
        when 6 => r6 := pr;
        when 7 => r7 := pr;
        when others => null;
      end case;
      ptr := ptr + 1;
    end if;
    if ptr < 1 then
      r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0; r3 := 0.0; r2 := 0.0;
      r1 := 0.0; r0 := 0.0;
    elsif ptr < 2 then
      r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0; r3 := 0.0; r2 := 0.0;
      r1 := 0.0;
    elsif ptr < 3 then
      r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0; r3 := 0.0; r2 := 0.0;
    elsif ptr < 4 then
      r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0; r3 := 0.0;
    elsif ptr < 5 then
      r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
    elsif ptr < 6 then
      r7 := 0.0; r6 := 0.0; r5 := 0.0;
    elsif ptr < 7 then
      r7 := 0.0; r6 := 0.0;
    elsif ptr < 8 then
      r7 := 0.0;
    end if;
  end renorm8;

  procedure fast_renorm
              ( x0,x1,x2,x3,x4,x5,x6,x7,x8 : in double_float;
                r0,r1,r2,r3,r4,r5,r6,r7 : out double_float ) is

    f0,f1,f2,f3,f4,f5,f6,f7,f8,pr : double_float;

  begin
    Double_Double_Basics.quick_two_sum(x7,x8,pr,f8);
    Double_Double_Basics.quick_two_sum(x6,pr,pr,f7);
    Double_Double_Basics.quick_two_sum(x5,pr,pr,f6);
    Double_Double_Basics.quick_two_sum(x4,pr,pr,f5);
    Double_Double_Basics.quick_two_sum(x3,pr,pr,f4);
    Double_Double_Basics.quick_two_sum(x2,pr,pr,f3);
    Double_Double_Basics.quick_two_sum(x1,pr,pr,f2);
    Double_Double_Basics.quick_two_sum(x0,pr,f0,f1);
    renorm8(f0,f1,f2,f3,f4,f5,f6,f7,f8,pr,r0,r1,r2,r3,r4,r5,r6,r7);
  end fast_renorm;

  procedure renorm_add1
              ( x0,x1,x2,x3,x4,x5,x6,x7 : in double_float;
                y : in double_float;
                r0,r1,r2,r3,r4,r5,r6,r7 : out double_float ) is

    f0,f1,f2,f3,f4,f5,f6,f7,f8,pr : double_float;

  begin
    Double_Double_Basics.two_sum(x7,y,pr,f8);
    Double_Double_Basics.two_sum(x6,pr,pr,f7);
    Double_Double_Basics.two_sum(x5,pr,pr,f6);
    Double_Double_Basics.two_sum(x4,pr,pr,f5);
    Double_Double_Basics.two_sum(x3,pr,pr,f4);
    Double_Double_Basics.two_sum(x2,pr,pr,f3);
    Double_Double_Basics.two_sum(x1,pr,pr,f2);
    Double_Double_Basics.two_sum(x0,pr,f0,f1);
    renorm8(f0,f1,f2,f3,f4,f5,f6,f7,f8,pr,r0,r1,r2,r3,r4,r5,r6,r7);
  end renorm_add1;

  procedure renorm10
              ( f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10 : in double_float;
                pr : in out double_float;
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 : out double_float ) is

    ptr : integer;

  begin
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
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f3,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f3,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f3,r2,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f4,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f4,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f4,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f4,r3,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f5,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f5,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f5,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f5,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f5,r4,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f6,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f6,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f6,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f6,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f6,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f6,r5,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f7,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f7,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f7,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f7,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f7,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f7,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f7,r6,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f8,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f8,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f8,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f8,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f8,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f8,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f8,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f8,r7,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f9,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f9,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f9,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f9,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f9,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f9,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f9,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f9,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f9,r8,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is
      when 0 => Double_Double_Basics.quick_two_sum(pr,f10,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f10,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f10,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f10,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f10,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f10,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f10,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f10,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f10,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f10,r9,pr);
      when others => null;
    end case;
    if pr = 0.0 then
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    if ptr < 10 and pr /= 0.0 then
      case ptr is
        when 0 => r0 := pr;
        when 1 => r1 := pr;
        when 2 => r2 := pr;
        when 3 => r3 := pr;
        when 4 => r4 := pr;
        when 5 => r5 := pr;
        when 6 => r6 := pr;
        when 7 => r7 := pr;
        when 8 => r8 := pr;
        when 9 => r9 := pr;
        when others => null;
      end case;
      ptr := ptr + 1;
    end if;
    if ptr < 1 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0; r2 := 0.0; r1 := 0.0; r0 := 0.0;
    elsif ptr < 2 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0; r2 := 0.0; r1 := 0.0;
    elsif ptr < 3 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0; r2 := 0.0;
    elsif ptr < 4 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0;
    elsif ptr < 5 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
    elsif ptr < 6 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0;
    elsif ptr < 7 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0;
    elsif ptr < 8 then
      r9 := 0.0; r8 := 0.0; r7 := 0.0;
    elsif ptr < 9 then
      r9 := 0.0; r8 := 0.0;
    elsif ptr < 10 then
      r9 := 0.0;
    end if;
  end renorm10;

  procedure fast_renorm
              ( x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 : in double_float;
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 : out double_float ) is

    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr : double_float;

  begin
    Double_Double_Basics.quick_two_sum(x9,x10,pr,f10);
    Double_Double_Basics.quick_two_sum(x8,pr,pr,f9);
    Double_Double_Basics.quick_two_sum(x7,pr,pr,f8);
    Double_Double_Basics.quick_two_sum(x6,pr,pr,f7);
    Double_Double_Basics.quick_two_sum(x5,pr,pr,f6);
    Double_Double_Basics.quick_two_sum(x4,pr,pr,f5);
    Double_Double_Basics.quick_two_sum(x3,pr,pr,f4);
    Double_Double_Basics.quick_two_sum(x2,pr,pr,f3);
    Double_Double_Basics.quick_two_sum(x1,pr,pr,f2);
    Double_Double_Basics.quick_two_sum(x0,pr,f0,f1);
    renorm10(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr,
             r0,r1,r2,r3,r4,r5,r6,r7,r8,r9);
  end fast_renorm;

  procedure renorm_add1
              ( x0,x1,x2,x3,x4,x5,x6,x7,x8,x9 : in double_float;
                y : in double_float;
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 : out double_float ) is

    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr : double_float;

  begin
    Double_Double_Basics.two_sum(x9,y,pr,f10);
    Double_Double_Basics.two_sum(x8,pr,pr,f9);
    Double_Double_Basics.two_sum(x7,pr,pr,f8);
    Double_Double_Basics.two_sum(x6,pr,pr,f7);
    Double_Double_Basics.two_sum(x5,pr,pr,f6);
    Double_Double_Basics.two_sum(x4,pr,pr,f5);
    Double_Double_Basics.two_sum(x3,pr,pr,f4);
    Double_Double_Basics.two_sum(x2,pr,pr,f3);
    Double_Double_Basics.two_sum(x1,pr,pr,f2);
    Double_Double_Basics.two_sum(x0,pr,f0,f1);
    renorm10(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,pr,
             r0,r1,r2,r3,r4,r5,r6,r7,r8,r9);
  end renorm_add1;

  procedure renorm16
              ( f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10 : in double_float;
                f11,f12,f13,f14,f15,f16 : in double_float;
                pr : in out double_float;
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 : out double_float;
                r10,r11,r12,r13,r14,r15 : out double_float ) is

    ptr : integer;

  begin
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
    if pr = 0.0 then -- first test on pr
      if ptr = 0
       then pr := r0;
       else pr := r1;
      end if;
    else
      ptr := ptr + 1;
    end if;
    case ptr is      -- first test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f3,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f3,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f3,r2,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- second test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is      -- second test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f4,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f4,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f4,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f4,r3,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- third test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is      -- third test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f5,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f5,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f5,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f5,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f5,r4,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- fourth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is      -- fourth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f6,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f6,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f6,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f6,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f6,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f6,r5,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- fifth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is      -- fifth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f7,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f7,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f7,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f7,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f7,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f7,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f7,r6,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- sixth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is      -- sixth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f8,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f8,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f8,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f8,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f8,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f8,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f8,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f8,r7,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- seventh test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- seventh test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f9,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f9,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f9,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f9,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f9,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f9,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f9,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f9,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f9,r8,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- eighth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- eighth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f10,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f10,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f10,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f10,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f10,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f10,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f10,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f10,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f10,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f10,r9,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- nineth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- nineth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f11,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f11,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f11,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f11,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f11,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f11,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f11,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f11,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f11,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f11,r9,pr);
      when 10 => Double_Double_Basics.quick_two_sum(pr,f11,r10,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- tenth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when 10 => pr := r10;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- tenth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f12,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f12,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f12,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f12,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f12,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f12,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f12,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f12,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f12,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f12,r9,pr);
      when 10 => Double_Double_Basics.quick_two_sum(pr,f12,r10,pr);
      when 11 => Double_Double_Basics.quick_two_sum(pr,f12,r11,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- eleventh test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when 10 => pr := r10;
        when 11 => pr := r11;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- eleventh test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f13,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f13,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f13,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f13,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f13,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f13,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f13,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f13,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f13,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f13,r9,pr);
      when 10 => Double_Double_Basics.quick_two_sum(pr,f13,r10,pr);
      when 11 => Double_Double_Basics.quick_two_sum(pr,f13,r11,pr);
      when 12 => Double_Double_Basics.quick_two_sum(pr,f13,r12,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- twelveth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when 10 => pr := r10;
        when 11 => pr := r11;
        when 12 => pr := r12;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- twelveth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f14,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f14,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f14,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f14,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f14,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f14,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f14,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f14,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f14,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f14,r9,pr);
      when 10 => Double_Double_Basics.quick_two_sum(pr,f14,r10,pr);
      when 11 => Double_Double_Basics.quick_two_sum(pr,f14,r11,pr);
      when 12 => Double_Double_Basics.quick_two_sum(pr,f14,r12,pr);
      when 13 => Double_Double_Basics.quick_two_sum(pr,f14,r13,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- thirteenth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when 10 => pr := r10;
        when 11 => pr := r11;
        when 12 => pr := r12;
        when 13 => pr := r13;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- thirteenth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f15,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f15,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f15,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f15,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f15,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f15,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f15,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f15,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f15,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f15,r9,pr);
      when 10 => Double_Double_Basics.quick_two_sum(pr,f15,r10,pr);
      when 11 => Double_Double_Basics.quick_two_sum(pr,f15,r11,pr);
      when 12 => Double_Double_Basics.quick_two_sum(pr,f15,r12,pr);
      when 13 => Double_Double_Basics.quick_two_sum(pr,f15,r13,pr);
      when 14 => Double_Double_Basics.quick_two_sum(pr,f15,r14,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- forteenth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when 10 => pr := r10;
        when 11 => pr := r11;
        when 12 => pr := r12;
        when 13 => pr := r13;
        when 14 => pr := r14;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
    case ptr is     -- forteenth test on ptr
      when 0 => Double_Double_Basics.quick_two_sum(pr,f16,r0,pr);
      when 1 => Double_Double_Basics.quick_two_sum(pr,f16,r1,pr);
      when 2 => Double_Double_Basics.quick_two_sum(pr,f16,r2,pr);
      when 3 => Double_Double_Basics.quick_two_sum(pr,f16,r3,pr);
      when 4 => Double_Double_Basics.quick_two_sum(pr,f16,r4,pr);
      when 5 => Double_Double_Basics.quick_two_sum(pr,f16,r5,pr);
      when 6 => Double_Double_Basics.quick_two_sum(pr,f16,r6,pr);
      when 7 => Double_Double_Basics.quick_two_sum(pr,f16,r7,pr);
      when 8 => Double_Double_Basics.quick_two_sum(pr,f16,r8,pr);
      when 9 => Double_Double_Basics.quick_two_sum(pr,f16,r9,pr);
      when 10 => Double_Double_Basics.quick_two_sum(pr,f16,r10,pr);
      when 11 => Double_Double_Basics.quick_two_sum(pr,f16,r11,pr);
      when 12 => Double_Double_Basics.quick_two_sum(pr,f16,r12,pr);
      when 13 => Double_Double_Basics.quick_two_sum(pr,f16,r13,pr);
      when 14 => Double_Double_Basics.quick_two_sum(pr,f16,r14,pr);
      when 15 => Double_Double_Basics.quick_two_sum(pr,f16,r15,pr);
      when others => null;
    end case;
    if pr = 0.0 then -- fifteenth test on pr
      case ptr is
        when 0 => pr := r0;
        when 1 => pr := r1;
        when 2 => pr := r2;
        when 3 => pr := r3;
        when 4 => pr := r4;
        when 5 => pr := r5;
        when 6 => pr := r6;
        when 7 => pr := r7;
        when 8 => pr := r8;
        when 9 => pr := r9;
        when 10 => pr := r10;
        when 11 => pr := r11;
        when 12 => pr := r12;
        when 13 => pr := r13;
        when 14 => pr := r14;
        when 15 => pr := r15;
        when others => null;
      end case;
    else
      ptr := ptr + 1;
    end if;
  -- beginning of the end ...
    if ptr < 16 and pr /= 0.0 then
      case ptr is
        when 0 => r0 := pr;
        when 1 => r1 := pr;
        when 2 => r2 := pr;
        when 3 => r3 := pr;
        when 4 => r4 := pr;
        when 5 => r5 := pr;
        when 6 => r6 := pr;
        when 7 => r7 := pr;
        when 8 => r8 := pr;
        when 9 => r9 := pr;
        when 10 => r10 := pr;
        when 11 => r11 := pr;
        when 12 => r12 := pr;
        when 13 => r13 := pr;
        when 14 => r14 := pr;
        when 15 => r15 := pr;
        when others => null;
      end case;
      ptr := ptr + 1;
    end if;
    if ptr < 1 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0; r2 := 0.0; r1 := 0.0; r0 := 0.0;
    elsif ptr < 2 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0; r2 := 0.0; r1 := 0.0;
    elsif ptr < 3 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0; r2 := 0.0;
    elsif ptr < 4 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
      r3 := 0.0;
    elsif ptr < 5 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0; r4 := 0.0;
    elsif ptr < 6 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0; r5 := 0.0;
    elsif ptr < 7 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0; r6 := 0.0;
    elsif ptr < 8 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0; r7 := 0.0;
    elsif ptr < 9 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0; r8 := 0.0;
    elsif ptr < 10 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
      r9 := 0.0;
    elsif ptr < 11 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0; r10 := 0.0;
    elsif ptr < 12 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0; r11 := 0.0;
    elsif ptr < 13 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0; r12 := 0.0;
    elsif ptr < 14 then
      r15 := 0.0; r14 := 0.0; r13 := 0.0;
    elsif ptr < 15 then
      r15 := 0.0; r14 := 0.0;
    elsif ptr < 16 then
      r15 := 0.0;
    end if;
  end renorm16;

  procedure fast_renorm
              ( x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 : in double_float;
                x11,x12,x13,x14,x15,x16 : in double_float;
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 : out double_float;
                r10,r11,r12,r13,r14,r15 : out double_float ) is

    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10 : double_float;
    f11,f12,f13,f14,f15,f16,pr : double_float;

  begin
    Double_Double_Basics.quick_two_sum(x15,x16,pr,f16);
    Double_Double_Basics.quick_two_sum(x14,pr,pr,f15);
    Double_Double_Basics.quick_two_sum(x13,pr,pr,f14);
    Double_Double_Basics.quick_two_sum(x12,pr,pr,f13);
    Double_Double_Basics.quick_two_sum(x11,pr,pr,f12);
    Double_Double_Basics.quick_two_sum(x10,pr,pr,f11);
    Double_Double_Basics.quick_two_sum(x9,pr,pr,f10);
    Double_Double_Basics.quick_two_sum(x8,pr,pr,f9);
    Double_Double_Basics.quick_two_sum(x7,pr,pr,f8);
    Double_Double_Basics.quick_two_sum(x6,pr,pr,f7);
    Double_Double_Basics.quick_two_sum(x5,pr,pr,f6);
    Double_Double_Basics.quick_two_sum(x4,pr,pr,f5);
    Double_Double_Basics.quick_two_sum(x3,pr,pr,f4);
    Double_Double_Basics.quick_two_sum(x2,pr,pr,f3);
    Double_Double_Basics.quick_two_sum(x1,pr,pr,f2);
    Double_Double_Basics.quick_two_sum(x0,pr,f0,f1);
    renorm16(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,pr,
             r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
  end fast_renorm;

  procedure renorm_add1
              ( x0,x1,x2,x3,x4,x5,x6,x7,x8,x9 : in double_float;
                x10,x11,x12,x13,x14,x15 : in double_float;
                y : in double_float;
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 : out double_float;
                r10,r11,r12,r13,r14,r15 : out double_float ) is

    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10 : double_float;
    f11,f12,f13,f14,f15,f16,pr : double_float;

  begin
    Double_Double_Basics.two_sum(x15,y,pr,f16);
    Double_Double_Basics.two_sum(x14,pr,pr,f15);
    Double_Double_Basics.two_sum(x13,pr,pr,f14);
    Double_Double_Basics.two_sum(x12,pr,pr,f13);
    Double_Double_Basics.two_sum(x11,pr,pr,f12);
    Double_Double_Basics.two_sum(x10,pr,pr,f11);
    Double_Double_Basics.two_sum(x9,pr,pr,f10);
    Double_Double_Basics.two_sum(x8,pr,pr,f9);
    Double_Double_Basics.two_sum(x7,pr,pr,f8);
    Double_Double_Basics.two_sum(x6,pr,pr,f7);
    Double_Double_Basics.two_sum(x5,pr,pr,f6);
    Double_Double_Basics.two_sum(x4,pr,pr,f5);
    Double_Double_Basics.two_sum(x3,pr,pr,f4);
    Double_Double_Basics.two_sum(x2,pr,pr,f3);
    Double_Double_Basics.two_sum(x1,pr,pr,f2);
    Double_Double_Basics.two_sum(x0,pr,f0,f1);
    renorm16(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,pr,
             r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15);
  end renorm_add1;

end Fast_Double_Renormalizations;
