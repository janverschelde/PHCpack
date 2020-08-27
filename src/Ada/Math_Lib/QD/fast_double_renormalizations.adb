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

end Fast_Double_Renormalizations;
