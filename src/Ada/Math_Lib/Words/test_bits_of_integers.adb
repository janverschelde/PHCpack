with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Characters_and_Numbers;
with Bits_of_Integers;

package body Test_Bits_of_Integers is

  procedure Test_Split is

    rnd,rhi,rlo,chk : integer64;
    maxint : constant integer64 := 2**62;
    threshold : constant integer64 := 2**32;

  begin
    put_line("Splitting a random 64-bit integer ...");
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
  end Test_Split;

  procedure Test_Quarter is

    rnd,rhihi,rlohi,rhilo,rlolo,chk : integer64;
    maxint : constant integer64 := 2**62;

  begin
    put_line("Quartering a random 64-bit integer ...");
    rnd := Standard_Random_Numbers.Random(0,maxint);
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
  end Test_Quarter;

  procedure Main is
  begin
    Test_Split;
    Test_Quarter;
  end Main;

end Test_Bits_of_Integers;
