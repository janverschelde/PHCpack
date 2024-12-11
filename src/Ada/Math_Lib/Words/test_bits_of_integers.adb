with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Characters_and_Numbers;
with Bits_of_Integers;

package body Test_Bits_of_Integers is

  procedure Test_Signed_Split is

    rnd,rhi,rlo,chk : integer64;
    maxint : constant integer64 := 2**62;
    threshold : constant integer64 := 2**32;

  begin
    put_line("Splitting a random signed 64-bit integer ...");
    rnd := Standard_Random_Numbers.Random(0,maxint);
    put("-> x : "); put(rnd); new_line;
    declare
      s : constant string := Characters_and_Numbers.Convert(rnd);
    begin
      put("-> s : "); put(s); new_line;
      put("-> S : "); put(integer32(s'length),1); new_line;
    end;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    Bits_of_Integers.Split(rnd,rhi,rlo);
    put("-> H : "); put(rhi); new_line;
    put("->   : "); put(rhi,20,b=>16); new_line;
    put("->   : "); put(rhi,65,b=>2); new_line;
    put("-> L : "); put(rlo); new_line;
    put("->   : "); put(rlo,20,b=>16); new_line;
    put("->   : "); put(rlo,65,b=>2); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    put("-> H : "); put(rhi,65,b=>2); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    put("-> L : "); put(rlo,65,b=>2); new_line;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> H : "); put(rhi,20,b=>16); new_line;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> L : "); put(rlo,20,b=>16); new_line;
    chk := rhi + rlo;
    put("-> c : "); put(chk); new_line;
    put("-> x : "); put(rnd);
    if chk = rnd
     then put_line(" Okay, H + L = x");
     else put_line(" H + L /= x !?");
    end if;
    put("-> rlo < 2**32 ? ");
    if rlo < threshold
     then put_line("True");
     else put_line("False");
    end if;
    if rhi > 0 then
      put("-> rhi >= 2**32 ? ");
      if rhi >= threshold
       then put_line("True");
       else put_line("False");
      end if;
    end if;
  end Test_Signed_Split;

  procedure Test_Unsigned_Split is

    use Bits_of_Integers;

    rnd,rhi,rlo,chk : unsigned_integer64;
    srnd : integer64;
    maxint : constant integer64 := 2**62;
    threshold : constant unsigned_integer64 := 2**32;

  begin
    put_line("Splitting a random unsigned 64-bit integer ...");
    srnd := Standard_Random_Numbers.Random(0,maxint);
    put("-> x : "); put(srnd); new_line;
    declare
      s : constant string := Characters_and_Numbers.Convert(srnd);
    begin
      put("-> s : "); put(s); new_line;
      put("-> S : "); put(integer32(s'length),1); new_line;
    end;
    rnd := unsigned_integer64(srnd);
    put("-> h : "); put(integer64(rnd),20,b=>16); new_line;
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    Bits_of_Integers.Split(rnd,rhi,rlo);
    put("-> H : "); put(integer64(rhi)); new_line;
    put("->   : "); put(integer64(rhi),20,b=>16); new_line;
    put("->   : "); put(integer64(rhi),65,b=>2); new_line;
    put("-> L : "); put(integer64(rlo)); new_line;
    put("->   : "); put(integer64(rlo),20,b=>16); new_line;
    put("->   : "); put(integer64(rlo),65,b=>2); new_line;
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    put("-> H : "); put(integer64(rhi),65,b=>2); new_line;
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    put("-> L : "); put(integer64(rlo),65,b=>2); new_line;
    put("-> h : "); put(integer64(rnd),20,b=>16); new_line;
    put("-> H : "); put(integer64(rhi),20,b=>16); new_line;
    put("-> h : "); put(integer64(rnd),20,b=>16); new_line;
    put("-> L : "); put(integer64(rlo),20,b=>16); new_line;
    chk := rhi + rlo;
    put("-> c : "); put(integer64(chk)); new_line;
    put("-> x : "); put(srnd);
    if chk = rnd
     then put_line(" Okay, H + L = x");
     else put_line(" H + L /= x !?");
    end if;
    put("-> rlo < 2**32 ? ");
    if rlo < threshold
     then put_line("True");
     else put_line("False");
    end if;
    if rhi > 0 then
      put("-> rhi >= 2**32 ? ");
      if rhi >= threshold
       then put_line("True");
       else put_line("False");
      end if;
    end if;
  end Test_Unsigned_Split;

  procedure Test_30_Bit_Split is

    use Bits_of_Integers;

    rnd,rhi,rlo,carry,chk : unsigned_integer64;
    rhi32,rlo32,carry32 : integer32;
    srnd : integer64;
    maxint : constant integer64 := 2**62;
    threshold30 : constant unsigned_integer64 := 2**30;
    threshold60 : constant unsigned_integer64 := 2**60;

  begin
    put_line("Splitting a random 64-bit integer in 30 bits ...");
    srnd := Standard_Random_Numbers.Random(0,maxint);
    while unsigned_integer64(srnd) < threshold60 loop
      srnd := 2*srnd;
    end loop;
    put("-> x : "); put(srnd); new_line;
    declare
      s : constant string := Characters_and_Numbers.Convert(srnd);
    begin
      put("-> s : "); put(s); new_line;
      put("-> S : "); put(integer32(s'length),1); new_line;
    end;
    rnd := unsigned_integer64(srnd);
    put("-> h : "); put(integer64(rnd),20,b=>16); new_line;
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    Bits_of_Integers.Split_30_Bits(rnd,rhi,rlo,carry);
    put("-> C : "); put(integer64(carry)); new_line;
    put("->   : "); put(integer64(carry),20,b=>16); new_line;
    put("->   : "); put(integer64(carry),65,b=>2); new_line;
    put("-> H : "); put(integer64(rhi)); new_line;
    put("->   : "); put(integer64(rhi),20,b=>16); new_line;
    put("->   : "); put(integer64(rhi),65,b=>2); new_line;
    put("-> L : "); put(integer64(rlo)); new_line;
    put("->   : "); put(integer64(rlo),20,b=>16); new_line;
    put("->   : "); put(integer64(rlo),65,b=>2); new_line;
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    put("-> H : "); put(integer64(rhi),65,b=>2); new_line;
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    put("-> C : "); put(integer64(carry),65,b=>2); new_line;
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    put("-> L : "); put(integer64(rlo),65,b=>2); new_line;
    put("-> h : "); put(integer64(rnd),20,b=>16); new_line;
    put("-> C : "); put(integer64(carry),20,b=>16); new_line;
    put("-> h : "); put(integer64(rnd),20,b=>16); new_line;
    put("-> H : "); put(integer64(rhi),20,b=>16); new_line;
    put("-> h : "); put(integer64(rnd),20,b=>16); new_line;
    put("-> L : "); put(integer64(rlo),20,b=>16); new_line;
    chk := rhi + rlo + carry;
    put("-> c : "); put(integer64(chk)); new_line;
    put("-> x : "); put(srnd);
    if chk = rnd
     then put_line(" Okay, C + H + L = x");
     else put_line(" C + H + L /= x !?");
    end if;
    put("-> rlo < 2**30 ? ");
    if rlo < threshold30
     then put_line("True");
     else put_line("False");
    end if;
    if rhi > 0 then
      put("-> rhi >= 2**30 ? ");
      if rhi >= threshold30
       then put_line("True");
       else put_line("False");
      end if;
      put("-> rhi < 2**60 ? ");
      if rhi < threshold60
       then put_line("True");
       else put_line("False");
      end if;
    end if;
    if carry > 0 then
      put("-> carry >= 2**60 ? ");
      if carry >= threshold60
       then put_line("True");
       else put_line("False");
      end if;
    end if;
    Bits_of_Integers.Split_30_Bit_Words(rnd,rhi32,rlo32,carry32);
    put("-> b : "); put(integer64(rnd),65,b=>2); new_line;
    put("-> C : "); put(integer64(carry32),5,b=>2); new_line;
    put("-> H : "); put(integer64(rhi32),35,b=>2); new_line;
    put("-> L : "); put(integer64(rlo32),65,b=>2); new_line;
  end Test_30_Bit_Split;

  procedure Test_Signed_Quarter is

    rnd,rhihi,rlohi,rhilo,rlolo,chk : integer64;
    maxint : constant integer64 := 2**62;

  begin
    put_line("Quartering a random signed 64-bit integer ...");
    rnd := 2*Standard_Random_Numbers.Random(0,maxint);
    put("-> x : "); put(rnd,1); new_line;
    declare
      s : constant string := Characters_and_Numbers.Convert(rnd);
    begin
      put("-> s : "); put(s); new_line;
      put("-> S : "); put(integer32(s'length),1); new_line;
    end;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    Bits_of_Integers.Quarter(rnd,rhihi,rlohi,rhilo,rlolo);
    put("-> A : "); put(rhihi); new_line;
    put("->   : "); put(rhihi,20,b=>16); new_line;
    put("->   : "); put(rhihi,65,b=>2); new_line;
    put("-> B : "); put(rlohi); new_line;
    put("->   : "); put(rlohi,20,b=>16); new_line;
    put("->   : "); put(rlohi,65,b=>2); new_line;
    put("-> C : "); put(rhilo); new_line;
    put("->   : "); put(rhilo,20,b=>16); new_line;
    put("->   : "); put(rhilo,65,b=>2); new_line;
    put("-> D : "); put(rlolo); new_line;
    put("->   : "); put(rlolo,20,b=>16); new_line;
    put("->   : "); put(rlolo,65,b=>2); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    put("-> A : "); put(rhihi,65,b=>2); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    put("-> B : "); put(rlohi,65,b=>2); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    put("-> C : "); put(rhilo,65,b=>2); new_line;
    put("-> b : "); put(rnd,65,b=>2); new_line;
    put("-> D : "); put(rlolo,65,b=>2); new_line;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> A : "); put(rhihi,20,b=>16); new_line;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> B : "); put(rlohi,20,b=>16); new_line;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> C : "); put(rhilo,20,b=>16); new_line;
    put("-> h : "); put(rnd,20,b=>16); new_line;
    put("-> D : "); put(rlolo,20,b=>16); new_line;
    chk := rhihi + rlohi + rhilo + rlolo;
    put("-> c : "); put(chk); new_line;
    put("-> x : "); put(rnd);
    if chk = rnd
     then put_line(" Okay, H + L = x");
     else put_line(" H + L /= x !?");
    end if;
  end Test_Signed_Quarter;

  procedure Test_Bit_Split is

    maxint : constant integer64 := 2**62;
    srnd,high,low : integer64;
    expo : integer32 := 0;
    bsz : integer32;

  begin
    put_line("Splitting a random 64-bit integer in bits ...");
    srnd := 2*Standard_Random_Numbers.Random(0,maxint);
    put("-> x : "); put(srnd,1); new_line;
    bsz := Bits_of_Integers.Bit_Size(srnd)+3;
    put("-> b : "); put(srnd,natural32(bsz),b=>2); new_line;
    put("Give the number of last bits : "); get(expo);
    Bits_of_Integers.Split_Bits(srnd,expo,high,low);
    put("Separating the last ");put(expo,1); put_line(" bits :");
    put("-> H : "); put(high,1); new_line;
    put("-> L : "); put(low,1); new_line;
    put("-> b : "); put(srnd,natural32(bsz),b=>2); new_line;
    if expo > 0 then
      put("-> H : "); put(high,natural32(bsz-expo),b=>2); new_line;
    else
      put("-> H : "); put(high,natural32(bsz+expo),b=>2); new_line;
    end if;
    put("-> L : "); put(low,natural32(bsz),b=>2); new_line;
  end Test_Bit_Split;

  procedure Main is
  begin
    Test_Signed_Split;
    Test_Unsigned_Split;
    Test_30_Bit_Split;
    Test_Signed_Quarter;
    Test_Bit_Split;
  end Main;

end Test_Bits_of_Integers;
