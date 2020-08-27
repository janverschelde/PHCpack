with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;

procedure ts_decadbl is

-- DESCRIPTION :
--   Test development procedure for deca double arithmetic,
--   extending the range of the quad double libary,
--   with code generated from the software CAMPARY.

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

  procedure Test_Multiplication_and_Division is

  -- DESCRIPTION :
  --   Generates two random numbers, multiplies and divides.

    x : constant deca_double := random;
    y : constant deca_double := random;
    z : constant deca_double := x*y;
    v : constant deca_double := z/y;

  begin
   put_line("All parts of a random deca double x :"); Write(x);
   put_line("All parts of a random deca double y :"); Write(y);
   put_line("All parts of x * y :"); Write(z);
   put_line("All parts of (x * y) / y :"); Write(v);
  end Test_Multiplication_and_Division;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for a test.

    ans : character;

  begin
    new_line;
    put_line("Testing deca double arithmetic ...");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication and division");
    put("Type 1 or 2 to select a test : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Test_Addition_and_Subtraction;
      when '2' => Test_Multiplication_and_Division;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_decadbl;
