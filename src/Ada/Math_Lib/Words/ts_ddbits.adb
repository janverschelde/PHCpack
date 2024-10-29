with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Constants;
with Double_Double_Basics;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with DoblDobl_Random_Numbers;
with Double_Double_Vectors;
with DoblDobl_Random_Vectors;
with Bits_of_Doubles;                    use Bits_of_Doubles;

procedure ts_ddbits is

-- DESCRIPTION :
--   Test on the bits of double doubles.

  procedure write_fraction_bits ( nbr : in double_float ) is

  -- DESCRIPTION :
  --   Writes the bits of the fraction of nbr in binary.

    frc : constant double_float := double_float'fraction(nbr);
    sfr : constant double_float := double_float'compose(frc, 52);
    mfr : constant integer64 := integer64(double_float'truncation(sfr));

  begin
    put(mfr,1,b=>2); new_line;
  end write_fraction_bits;

  procedure Test_Mod_Split ( nbr : in double_double ) is

  -- DESCRIPTION :
  --   Chops the doubles of a double double nbr in half.

    hinb : constant double_float := hi_part(nbr);
    lonb : constant double_float := lo_part(nbr);
    nb0,nb1,nb2,nb3 : double_float;
    ddnbr : double_double; 
    qdnbr : quad_double; 

  begin
    put("dd nb : "); put(nbr); new_line;
   -- Mod_Split(hinb,nb0,nb1);
    Double_Double_Basics.split(hinb,nb0,nb1);
    put("hi nb :"); put(hinb); new_line;
    put("hi nb :"); write_fraction_bits(hinb);
    put("nb0+1 :"); standard_floating_numbers_io.put(nb0+nb1); new_line;
   -- Mod_Split(lonb,nb2,nb3);
    Double_Double_Basics.split(lonb,nb2,nb3);
    put("lo nb :"); put(lonb); new_line;
    put("nb2+3 :"); standard_floating_numbers_io.put(nb2+nb3); new_line;
    put("hinb0 :"); put(nb0); new_line;
    put("hinbi0 :"); write_fraction_bits(nb0);
    put("hinb1 :"); put(nb1); new_line;
    put("hinbi1 :"); write_fraction_bits(nb1);
    put("nb 2 :"); put(nb2); new_line;
    put("nb 3 :"); put(nb3); new_line;
    ddnbr := create(nb0+nb1,nb2+nb3);
    qdnbr := create(nb0,nb1,nb2,nb3);
    put(ddnbr); new_line;
    put(nbr); new_line;
    put(qdnbr,32); new_line;
    put("error : "); double_double_numbers_io.put(abs(nbr - ddnbr),2);
    new_line;
    put("error : "); quad_double_numbers_io.put(abs(nbr - qdnbr),2);
    new_line;
  end Test_Mod_Split;

  procedure Test_Vec_Split ( nbr : in double_double ) is

  -- DESCRIPTION :
  --   Chops the doubles of a double double nbr in half,
  --   using a vector of natural numbers.

    hinb : constant double_float := hi_part(nbr);
    lonb : constant double_float := lo_part(nbr);
    nb0,nb1,nb2,nb3 : double_float;
    ddnbr : double_double; 
    qdnbr : quad_double; 

  begin
    put("dd nb : "); put(nbr); new_line;
    Vec_Split(hinb,nb0,nb1);
    put("hi nb :"); put(hinb); new_line;
    put("nb0+1 :"); standard_floating_numbers_io.put(nb0+nb1); new_line;
    Vec_Split(lonb,nb2,nb3);
    put("lo nb :"); put(lonb); new_line;
    put("nb2+3 :"); standard_floating_numbers_io.put(nb2+nb3); new_line;
    put("nb 0 :"); put(nb0); new_line;
    put("nb 1 :"); put(nb1); new_line;
    put("nb 2 :"); put(nb2); new_line;
    put("nb 3 :"); put(nb3); new_line;
    ddnbr := create(nb0+nb1,nb2+nb3);
    qdnbr := create(nb0,nb1,nb2,nb3);
    put(ddnbr); new_line;
    put(nbr); new_line;
    put(qdnbr,32); new_line;
    put("error : "); double_double_numbers_io.put(abs(nbr - ddnbr),2);
    new_line;
    put("error : "); quad_double_numbers_io.put(abs(nbr - qdnbr),2);
    new_line;
  end Test_Vec_Split;

  procedure Test_Split is

  -- DESCRIPTION :
  --   Tests the splitting of the fraction of a double
  --   in two equal halves.

    ddpi : constant double_double := Double_Double_Constants.pi;
    ddnb : constant double_double := abs(DoblDobl_Random_Numbers.Random);

  begin
   -- put_line("testing pi ...");
   -- Test_Mod_Split(ddpi);
   -- Test_Vec_Split(ddpi);
    put_line("testing a random number ...");
    Test_Mod_Split(ddnb);
   -- Test_Vec_Split(ddnb);
  end Test_Split;

  procedure Test_Sum ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the sum of the elements in a random vector
  --   of doubles of dimension dim.

    x : Double_Double_Vectors.Vector(1..dim)
      := DoblDobl_Random_Vectors.Random_Vector(1,dim);
    ddsum0 : double_double := create(0.0);
    x0,x1,x2,x3,sum,err : double_float;
    s0,s1,s2,s3 : double_float := 0.0;
    ddsum1 : double_double;

  begin
    for i in x'range loop
      if x(i) < 0.0
       then Min(x(i));
      end if;
      ddsum0 := ddsum0 + x(i);
     -- Mod_Split(hi_part(x(i)),x0,x1);
      Double_Double_Basics.split(hi_part(x(i)),x0,x1);
      put("x0 : "); write_fraction_bits(x0);
      put("x1 : "); write_fraction_bits(x1);
     -- Mod_Split(lo_part(x(i)),x2,x3);
      Double_Double_Basics.split(lo_part(x(i)),x2,x3);
     -- put("x2 : "); write_fraction_bits(x2);
     -- put("x3 : "); write_fraction_bits(x3);
      s0 := s0 + x0;
      s1 := s1 + x1;
      put("s0 : "); write_fraction_bits(s0);
      put("s1 : "); write_fraction_bits(s1);
      s2 := s2 + x2;
      s3 := s3 + x3;
    end loop;
    put("s0 :"); put(s0); new_line;
    put("s1 :"); put(s1); new_line;
    put("s2 :"); put(s2); new_line;
    put("s3 :"); put(s3); new_line;
    put("dd sum 0 : "); put(ddsum0); new_line;
    put("dd sum 0 hi :"); put(hi_part(ddsum0)); new_line;
    put("dd sum 0 lo :"); put(lo_part(ddsum0)); new_line;
    Double_Double_Basics.quick_two_sum(s0,s1,sum,err);
    ddsum1 := create(sum,err);
    put("dd sum 1 : "); put(ddsum1); new_line;
    put("dd sum 1 hi :"); put(hi_part(ddsum1)); new_line;
    put("dd sum 1 lo :"); put(lo_part(ddsum1)); new_line;
    put("error : ");
    double_double_numbers_io.put(abs(ddsum0 - ddsum1),2);
    new_line;
  end Test_Sum;

  procedure Main is

  -- DESCRIPTION :
  --   Runs the tests.

  begin
    Test_Split;
    Test_Sum(1000);
  end Main;

begin
  Main;
end ts_ddbits;
