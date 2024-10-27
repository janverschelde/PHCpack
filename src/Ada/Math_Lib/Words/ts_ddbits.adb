with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
-- with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Double_Double_Constants;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Bits_of_Doubles;                    use Bits_of_Doubles;
with Mask_Bits_of_Doubles;

procedure ts_ddbits is

-- DESCRIPTION :
--   Test on the bits of double doubles.

  procedure Mod_Split ( x : in double_float;
                        xhi,xlo : out double_float ) is

  -- DESCRIPTION :
  --   Splits the fraction of the double x in two equal halves,
  --   using modular arithmetic, without vectors of natural numbers.

    e : constant integer32 := integer32(double_float'exponent(x));
    f : constant double_float := double_float'fraction(x);
    s : constant double_float := double_float'compose(f, 52);
    m : constant integer64 := integer64(double_float'truncation(s));
    mlast,mchop : integer64;
    cnt : integer32 := 0;

  begin
    mlast := Mask_Bits_of_Doubles.last_bits(m,26);
    mchop := m - mlast;
   -- put("expo : "); put(e,1); new_line;
   -- put("frac : "); put(m,1,b=>2); new_line;
   -- put("chop : "); put(mchop,1,b=>2); new_line;
   -- put("last : "); put(mlast,1,b=>2); new_line;
    xhi := double_float'compose(double_float(mchop),e);
    while mlast < 2**natural(25-cnt) loop
      cnt := cnt + 1;
    end loop;
   -- put(" cnt : "); put(cnt,1); new_line;
    xlo := double_float'compose(double_float(mlast),e-26-cnt);
  end Mod_Split;

  procedure Vec_Split ( x : in double_float;
                        xhi,xlo : out double_float ) is

  -- DESCRIPTION :
  --   Splits the fraction of the double x in two equal halves,
  --   using vectors of natural numbers.

    head,tail : Standard_Natural_Vectors.Vector(0..51);
    val : integer64;
    expo : constant integer32 := integer32(double_float'exponent(x));
    cnt : integer32 := 0;

  begin
    xhi := x;
    chop_last_bits(xhi,26,head,tail);
    val := value_52bits(tail);
    while tail(cnt+26) = 0 loop
      cnt := cnt + 1;
    end loop;
    xlo := double_float'compose(double_float(val),expo-26-cnt);
  end Vec_Split;

  procedure Test_Mod_Split is

  -- DESCRIPTION :
  --   Chops the doubles of a double double in half.

    ddpi : constant double_double := Double_Double_Constants.pi;
    hipi : constant double_float := hi_part(ddpi);
    lopi : constant double_float := lo_part(ddpi);
    pi0,pi1,pi2,pi3 : double_float;
    ddnbr : double_double; 
    qdnbr : quad_double; 

  begin
    put("dd pi : "); put(ddpi); new_line;
    Mod_Split(hipi,pi0,pi1);
    put("hi pi :"); put(hipi); new_line;
    put("pi0+1 :"); standard_floating_numbers_io.put(pi0+pi1); new_line;
    Mod_Split(lopi,pi2,pi3);
    put("lo pi :"); put(lopi); new_line;
    put("pi2+3 :"); standard_floating_numbers_io.put(pi2+pi3); new_line;
    put("pi 0 :"); put(pi0); new_line;
    put("pi 1 :"); put(pi1); new_line;
    put("pi 2 :"); put(pi2); new_line;
    put("pi 3 :"); put(pi3); new_line;
    ddnbr := create(pi0+pi1,pi2+pi3);
    qdnbr := create(pi0,pi1,pi2,pi3);
    put(ddnbr); new_line;
    put(ddpi); new_line;
    put(qdnbr,32); new_line;
    put("error : "); double_double_numbers_io.put(abs(ddpi - ddnbr),2);
    new_line;
    put("error : "); quad_double_numbers_io.put(abs(ddpi - qdnbr),2);
    new_line;
  end Test_Mod_Split;

  procedure Test_Vec_Split is

  -- DESCRIPTION :
  --   Chops the doubles of a double double in half,
  --   using a vector of natural numbers.

    ddpi : constant double_double := Double_Double_Constants.pi;
    hipi : constant double_float := hi_part(ddpi);
    lopi : constant double_float := lo_part(ddpi);
    pi0,pi1,pi2,pi3 : double_float;
    ddnbr : double_double; 
    qdnbr : quad_double; 

  begin
    put("dd pi : "); put(ddpi); new_line;
    Vec_Split(hipi,pi0,pi1);
    put("hi pi :"); put(hipi); new_line;
    put("pi0+1 :"); standard_floating_numbers_io.put(pi0+pi1); new_line;
    Vec_Split(lopi,pi2,pi3);
    put("lo pi :"); put(lopi); new_line;
    put("pi2+3 :"); standard_floating_numbers_io.put(pi2+pi3); new_line;
    put("pi 0 :"); put(pi0); new_line;
    put("pi 1 :"); put(pi1); new_line;
    put("pi 2 :"); put(pi2); new_line;
    put("pi 3 :"); put(pi3); new_line;
    ddnbr := create(pi0+pi1,pi2+pi3);
    qdnbr := create(pi0,pi1,pi2,pi3);
    put(ddnbr); new_line;
    put(ddpi); new_line;
    put(qdnbr,32); new_line;
    put("error : "); double_double_numbers_io.put(abs(ddpi - ddnbr),2);
    new_line;
    put("error : "); quad_double_numbers_io.put(abs(ddpi - qdnbr),2);
    new_line;
  end Test_Vec_Split;

begin
  Test_Mod_Split;
  Test_Vec_Split;
end ts_ddbits;
