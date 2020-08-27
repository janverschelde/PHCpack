with text_io;                            use text_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Double_Double_Basics;

procedure ts_decadbl is

-- DESCRIPTION :
--   Test development procedure for deca double arithmetic,
--   extending the range of the quad double libary,
--   with code generated from the software CAMPARY.

-- DATA STRUCTURE :
--   A deca double is stored as a record of ten doubles,
--   ordered from high to low significance.
--   The naming of the ten doubles follows the typical count
--   of a right handed person, starting with the fingers on the
--   right hand and then moving on to the fingers of the left hand.

  type deca_double is record
    right_thumb : double_float;  -- most significant double
    right_index : double_float;  -- second most significant double
    right_middle : double_float; -- third most significant double
    right_ring : double_float;   -- fourth most significant double
    right_pink : double_float;   -- fifth most significant double
    left_thumb : double_float;   -- fifth least significant double
    left_index : double_float;   -- fourth least significant double
    left_middle : double_float;  -- third least significant double
    left_ring : double_float;    -- second least significant double
    left_pink : double_float;    -- least significant double
  end record;

  function create ( x : double_float ) return deca_double is

  -- DESCRIPTION :
  --   The most significant part of the deca double on return is x.

    res : deca_double;

  begin
    res.right_thumb := x;
    res.right_index := 0.0;
    res.right_middle := 0.0;
    res.right_ring := 0.0;
    res.right_pink := 0.0;
    res.left_thumb := 0.0;
    res.left_index := 0.0;
    res.left_middle := 0.0;
    res.left_ring := 0.0;
    res.left_pink := 0.0;
    return res;
  end create;

  function thumb_right ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the most significant double of x.

  begin
    return x.right_thumb;
  end thumb_right;

  function index_right ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the second most significant double of x.

  begin
    return x.right_index;
  end index_right;

  function middle_right ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the third significant double of x.

  begin
    return x.right_middle;
  end middle_right;

  function ring_right ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the fourth most significant double of x.

  begin
    return x.right_ring;
  end ring_right;

  function pink_right ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the fifth most significant double of x.

  begin
    return x.right_pink;
  end pink_right;

  function thumb_left ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the fifth least significant double of x.

  begin
    return x.left_thumb;
  end thumb_left;

  function index_left ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the fourth least significant double of x.

  begin
    return x.left_index;
  end index_left;

  function middle_left ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the third least significant double of x.

  begin
    return x.left_middle;
  end middle_left;

  function ring_left ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the second least significant double of x.

  begin
    return x.left_ring;
  end ring_left;

  function pink_left ( x : deca_double ) return double_float is

  -- DESCRIPTION :
  --   Returns the least significant double of x.

  begin
    return x.left_pink;
  end pink_left;

  procedure write ( x : deca_double ) is

  -- DESCRIPTION :
  --   Writes all parts of x, one part per line.

  begin
    put("  thumb right  : "); put(thumb_right(x),2,17,3); new_line;
    put("  index right  : "); put(index_right(x),2,17,3); new_line;
    put("  middle right : "); put(middle_right(x),2,17,3); new_line;
    put("  ring right   : "); put(ring_right(x),2,17,3); new_line;
    put("  pink right   : "); put(pink_right(x),2,17,3); new_line;
    put("  thumb left   : "); put(thumb_left(x),2,17,3); new_line;
    put("  index left   : "); put(index_left(x),2,17,3); new_line;
    put("  middle left  : "); put(middle_left(x),2,17,3); new_line;
    put("  ring left    : "); put(ring_left(x),2,17,3); new_line;
    put("  pink left    : "); put(pink_left(x),2,17,3); new_line;
  end Write;

  procedure renorm10
              ( f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10 : in double_float;
                pr : in out double_float;
                r0,r1,r2,r3,r4,r5,r6,r7,r8,r9 : out double_float ) is

  -- DESCRIPTION :
  --   Definitions common to fast_renorm2L<11,10> and renorm2L_4Add1<10,10>
  --   from code of the specRenorm.h, generated by the CAMPARY library,
  --   in an unrolled form (Valentina Popescu, Mioara Joldes), with
  --   double double basics of QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  -- ON ENTRY :
  --   x0       most significant word;
  --   x1       second most significant word;
  --   x2       third most significant word;
  --   x3       fourth most significant word;
  --   x4       fifth most significant word;
  --   x5       sixth most significant word;
  --   x6       seventh most significant word;
  --   x7       eighth most significant word;
  --   x8       nineth most significant word;
  --   x9       tenth most significant word;
  --   x10      least significant word;
  --   pr       computed by the start of the renormalization.

  -- ON RETURN :
  --   pr       updated value by renormalization;
  --   r0       highest part of a deca double number;
  --   r1       second highest part of a deca double number;
  --   r2       third highest part of a deca double number.
  --   r3       fourth highest part of a deca double number;
  --   r4       fifth highest part of a deca double number;
  --   r5       fifth lowest part of a deca double number;
  --   r6       fourth lowest part of a deca double number;
  --   r7       third lowest part of a deca double number;
  --   r8       second lowest part of a deca double number;
  --   r9       lowest part of a deca double number.

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

  -- DESCRIPTION :
  --   The definition is based on the fast_renorm2L<11,10>,
  --   from code of the specRenorm.h, generated by the CAMPARY library,
  --   in an unrolled form (Valentina Popescu, Mioara Joldes), with
  --   double double basics of QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  -- ON ENTRY :
  --   x0       most significant word;
  --   x1       second most significant word;
  --   x2       third most significant word;
  --   x3       fourth most significant word;
  --   x4       fifth most significant word;
  --   x5       sixth most significant word;
  --   x6       seventh most significant word;
  --   x7       eighth most significant word;
  --   x8       nineth most significant word;
  --   x9       tenth most significant word;
  --   x10      least significant word.

  -- ON RETURN :
  --   r0       highest part of a deca double number;
  --   r1       second highest part of a deca double number;
  --   r2       third highest part of a deca double number.
  --   r3       fourth highest part of a deca double number;
  --   r4       fifth highest part of a deca double number;
  --   r5       fifth lowest part of a deca double number;
  --   r6       fourth lowest part of a deca double number;
  --   r7       third lowest part of a deca double number;
  --   r8       second lowest part of a deca double number;
  --   r9       lowest part of a deca double number.

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

  -- DESCRIPTION :
  --   The definition is based on renorm2L_4Add1<10,10>
  --   from code of the specRenorm.h, generated by the CAMPARY library,
  --   in an unrolled form (Valentina Popescu, Mioara Joldes), with
  --   double double basics of QD-2.3.9 (Y. Hida, X.S. Li, and D.H. Bailey).

  -- ON ENTRY :
  --   x0       most significant word of a deca double x;
  --   x1       second most significant word of a deca double x;
  --   x2       third most significant word of a deca double x;
  --   x3       fourth most significant word of a deca double x;
  --   x4       fifth most significant word of a deca double x;
  --   x5       fifth least significant word of a deca double x;
  --   x6       fourth least significant word of a deca double x;
  --   x7       third least significant word of a deca double x;
  --   x8       second least significant word of a deca double x;
  --   x9       least significant word of a deca double x;
  --   y        double to be added to x.

  -- ON RETURN :
  --   r0       most significant word of x + y;
  --   r1       second most significant word of x + y;
  --   r2       third most significant word of x + y;
  --   r3       fourth most significant word of x + y;
  --   r4       fifth most significant word of x + y;
  --   r5       fifth least significant word of x + y;
  --   r6       fourth least significant word of x + y;
  --   r7       third least significant word of x + y;
  --   r8       second least significant word of x + y;
  --   r9       least significant word of x + y.

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

  function "+" ( x,y : deca_double ) return deca_double is

  -- ALGORITHM : baileyAdd_fast<10,10,10> generated by CAMPARY.
   
    res : deca_double;
    f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,e : double_float;

  begin
    f10 := 0.0;
    Double_Double_Basics.two_sum(x.left_pink,y.left_pink,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_ring,y.left_ring,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_middle,y.left_middle,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_index,y.left_index,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.left_thumb,y.left_thumb,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_pink,y.right_pink,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_ring,y.right_ring,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_middle,y.right_middle,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_index,y.right_index,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    Double_Double_Basics.two_sum(x.right_thumb,y.right_thumb,f0,e);
    Double_Double_Basics.two_sum(f1,e,f1,e);
    Double_Double_Basics.two_sum(f2,e,f2,e);
    Double_Double_Basics.two_sum(f3,e,f3,e);
    Double_Double_Basics.two_sum(f4,e,f4,e);
    Double_Double_Basics.two_sum(f5,e,f5,e);
    Double_Double_Basics.two_sum(f6,e,f6,e);
    Double_Double_Basics.two_sum(f7,e,f7,e);
    Double_Double_Basics.two_sum(f8,e,f8,e);
    Double_Double_Basics.two_sum(f9,e,f9,e);
    f10 := f10 + e;
    fast_renorm(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end "+";

  function "+" ( x : deca_double; y : double_float ) return deca_double is

    res : deca_double;

  begin
    renorm_add1(x.right_thumb,x.right_index,x.right_middle,
                x.right_ring,x.right_pink,
                x.left_thumb,x.left_index,x.left_middle,
                x.left_ring,x.left_pink,y,
                res.right_thumb,res.right_index,res.right_middle,
                res.right_ring,res.right_pink,
                res.left_thumb,res.left_index,res.left_middle,
                res.left_ring,res.left_pink);
    return res;
  end "+";

  function "+" ( x : deca_double ) return deca_double is

  -- DESCRIPTION :
  --   Returns a copy of x.

    res : deca_double;

  begin
    res.right_thumb := x.right_thumb;
    res.right_index := x.right_index;
    res.right_middle := x.right_middle;
    res.right_ring := x.right_ring;
    res.right_pink := x.right_pink;
    res.left_thumb := x.left_thumb;
    res.left_index := x.left_index;
    res.left_middle := x.left_middle;
    res.left_ring := x.left_ring;
    res.left_pink := x.left_pink;
    return res;
  end "+";

  function "-" ( x : deca_double ) return deca_double is

  -- DESCRIPTION :
  --   Returns a copy of x, with all signs flipped.

    res : deca_double;

  begin
    res.right_thumb := -x.right_thumb;
    res.right_index := -x.right_index;
    res.right_middle := -x.right_middle;
    res.right_ring := -x.right_ring;
    res.right_pink := -x.right_pink;
    res.left_thumb := -x.left_thumb;
    res.left_index := -x.left_index;
    res.left_middle := -x.left_middle;
    res.left_ring := -x.left_ring;
    res.left_pink := -x.left_pink;
    return res;
  end "-";

  function "-" ( x,y : deca_double ) return deca_double is

    mny : constant deca_double := -y;
    res : constant deca_double := x + mny;

  begin
    return res;
  end "-";

  function random return deca_double is

  -- DESCRIPTION :
  --   Returns a random deca double number from adding
  --   random double numbers in [-1,+1].

    res : deca_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..10 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end random;

  procedure Test_Addition_and_Subtraction is

  -- DESCRIPTION :
  --   Generates two random numbers, adds and subtracts.

    x : constant deca_double := random;
    y : constant deca_double := random;
    z : constant deca_double := x + y;
    v : constant deca_double := z - y;

  begin
   put_line("All parts of a random deca double x :"); Write(x);
   put_line("All parts of a random deca double y :"); Write(y);
   put_line("All parts of x + y :"); Write(z);
   put_line("All parts of (x + y) - y :"); Write(v);
  end Test_Addition_and_Subtraction;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for a test.

  begin
    new_line;
    put_line("Testing deca double arithmetic ...");
    Test_Addition_and_Subtraction;
  end Main;

begin
  Main;
end ts_decadbl;
