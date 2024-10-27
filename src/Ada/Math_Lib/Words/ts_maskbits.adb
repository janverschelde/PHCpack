with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Mask_Bits_of_Doubles;               use Mask_Bits_of_Doubles;

procedure ts_maskbits is

-- DESCRIPTION :
--   Test on bit masking.

  procedure test_last_bits ( lastbits : in natural32 ) is

  -- DESCRIPTION :
  --   Tests the extraction of the last bits,
  --   equal to the value of lastbits.

    r : constant unsigned_integer64 := 2**natural(lastbits);
    x : unsigned_integer64 := r;
    m : constant unsigned_integer64 := last_mask(lastbits);
    y : unsigned_integer64;

  begin
    put("   x : "); put(integer64(x),1,b=>2); new_line;
    put("mask :  "); put(integer64(m),1,b=>2); new_line;
    for i in 1..r loop
      y := x and m;
      put("last bits of "); put(integer64(x)); put(" : ");
      put(integer64(y),1,b=>2); put(" = ");
      put(last_bits(integer64(x),lastbits),1,b=>2);
      if y = last_bits(x,lastbits)
       then put_line(" true");
       else put_line(" false");
      end if;
      x := x + 1;
    end loop;
  end test_last_bits;

  procedure test_last_bits is

  -- DESCRIPTION :
  --   Tests the extraction of the last bits.

  begin
    for i in 1..4 loop
      put("testing the extraction of the last ");
      put(integer32(i),1); put_line(" bits ...");
      test_last_bits(natural32(i));
    end loop;
  end test_last_bits;

  procedure test_first_bits ( firstbits : in natural32 ) is

  -- DESCRITPION :
  --   Tests the extraction of the first bits,
  --   equal to the value of first bits.

    r : constant unsigned_integer64 := 2**natural(firstbits);
    xexp : constant integer32 := 51 - integer32(firstbits);
    xfloor : constant unsigned_integer64 := 2**natural(xexp);
    modulus : constant unsigned_integer64 := 2*xfloor;
    x : unsigned_integer64 := xfloor;
    y,z : unsigned_integer64;

  begin
    put("2^"); put(xexp); put(" : ");
    put(integer64(xfloor),1,b=>2); new_line;
    put("   x : "); put(integer64(x),1,b=>2); new_line;
    for i in 0..(r-1) loop
     -- y := (x - last_bits(x,52 - firstbits));
      y := first_bits(x,firstbits);
      z := y/modulus;
      put(integer64(x),59,b=>2); new_line;
      put(integer64(y),59,b=>2); new_line;
      put(integer64(z),10,b=>2);
      put(" = "); put(integer32(i));
      if z = i
       then put_line(" true");
       else put_line(" false");
      end if;
      x := x + modulus;
    end loop;
  end test_first_bits;

  procedure test_first_bits is

  -- DESCRIPTION :
  --   Tests the extraction of the first bits.

  begin
    for i in 1..4 loop
      put("testing the extraction of the first ");
      put(integer32(i),1); put_line(" bits ...");
      test_first_bits(natural32(i));
    end loop;
  end test_first_bits;

begin
  test_last_bits;
  test_first_bits;
end ts_maskbits;
