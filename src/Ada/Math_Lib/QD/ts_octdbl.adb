with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;

procedure ts_octdbl is

-- DESCRIPTION :
--   Development of arithmetic with 8-fold doubles "octo doubles"
--   with the help of generated code of the CAMPARY library.

  function random return octo_double is

  -- DESCRIPTION :
  --   Returns a random octo double number from adding
  --   random double numbers in [-1,+1].

    res : octo_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..8 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end random;

  procedure Write ( x : in octo_double ) is

  -- DESCRIPTION :
  --   Writes all parts of x, one part per line.

  begin
    put("  hihihi : "); put(hihihi_part(x)); new_line;
    put("  lohihi : "); put(lohihi_part(x)); new_line;
    put("  hilohi : "); put(hilohi_part(x)); new_line;
    put("  lolohi : "); put(lolohi_part(x)); new_line;
    put("  hihilo : "); put(hihilo_part(x)); new_line;
    put("  lohilo : "); put(lohilo_part(x)); new_line;
    put("  hilolo : "); put(hilolo_part(x)); new_line;
    put("  lololo : "); put(lololo_part(x)); new_line;
  end Write;

  procedure Test_Add_and_Subtract is

  -- DESCRIPTION :
  --   Generates two random numbers, adds and subtracts.

    x : constant octo_double := random;
    y : constant octo_double := random;
    z : constant octo_double := x + y;
    v : constant octo_double := z - y;

  begin
   put_line("All parts of a random octo double x :"); Write(x);
   put_line("All parts of a random octo double y :"); Write(y);
   put_line("All parts of x + y :"); Write(z);
   put_line("All parts of (x + y) - y :"); Write(v);
  end Test_Add_and_Subtract;

  procedure Test_Multiplication_and_Division is

  -- DESCRIPTION :
  --   Generates two random numbers, multiplies and divides.

    x : constant octo_double := random;
    y : constant octo_double := random;
    z : constant octo_double := x*y;
    v : constant octo_double := z/y;

  begin
   put_line("All parts of a random octo double x :"); Write(x);
   put_line("All parts of a random octo double y :"); Write(y);
   put_line("All parts of x * y :"); Write(z);
   put_line("All parts of (x * y) / y :"); Write(v);
  end Test_Multiplication_and_Division;

  procedure Main is

  -- DESCRIPTION :
  --   Launches tests.

    ans : character;

  begin
    new_line;
    put_line("Testing octo double arithmetic ...");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication and division");
    put("Type 1 or 2 to select a test : "); Ask_Alternative(ans,"12");
    case ans is
      when '1' => Test_Add_and_Subtract;
      when '2' => Test_Multiplication_and_Division;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_octdbl;
