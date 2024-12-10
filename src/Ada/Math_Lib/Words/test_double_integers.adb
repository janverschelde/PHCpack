with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Integer64_Numbers_io;      use Multprec_Integer64_Numbers_io;
with Standard_Random_Numbers;
with Bits_of_Integers;
with Double_Integer_Arithmetic;          use Double_Integer_Arithmetic;

package body Test_Double_Integers is

  maxint52 : constant integer64 := 2**52;
  maxint62 : constant integer64 := 2**62;

  procedure Random_Double52_Integer ( hi,lo : out integer64 ) is
  begin
    hi := Standard_Random_Numbers.Random(0,maxint52);
    lo := Standard_Random_Numbers.Random(0,maxint52);
  end Random_Double52_Integer;

  procedure Random_Double64_Integer ( hi,lo : out integer64 ) is
  begin
    hi := Standard_Random_Numbers.Random(0,maxint62);
    lo := Standard_Random_Numbers.Random(0,maxint62);
  end Random_Double64_Integer;

  function Value52 ( hi,lo : integer64 ) return Integer_Number is

    res : Integer_Number := create64(hi);
    base : constant integer64 := 2**52;
    mpbase : Integer_Number := create64(base);

  begin
    Mul(res,mpbase);
    Add(res,lo);
    Clear(mpbase);
    return res;
  end Value52;

  function Value60 ( hi,lo : integer64; verbose : boolean := true )
                   return Integer_Number is

    res : Integer_Number;
    hihi,lohi,hilo,lolo,hicarry,locarry : integer32;
    mphihi,mplohi,mphilo,mplolo : Integer_Number;
    base : constant integer32 := 2**30;
    mpbase : Integer_Number := create32(base);

  begin
    Bits_of_Integers.Split_30_Bit_Words(hi,hihi,lohi,hicarry);
    if verbose then
      put("hihi : "); put(hihi); new_line;
      put("lohi : "); put(lohi); new_line;
    end if;
    Bits_of_Integers.Split_30_Bit_Words(lo,hilo,lolo,locarry);
    if locarry > 0 then
      lohi := lohi + locarry;
      if verbose
       then put("lohi : "); put(lohi); new_line;
      end if;
    end if;
    if verbose then
      put("hilo : "); put(hilo); new_line;
      put("lolo : "); put(lolo); new_line;
    end if;
    mphihi := Multprec_Integer64_Numbers.create32(hihi);
    mplohi := Multprec_Integer64_Numbers.create32(lohi);
    mphilo := Multprec_Integer64_Numbers.create32(hilo);
    mplolo := Multprec_Integer64_Numbers.create32(lolo);
    if hicarry = 0 then
      res := integer64(hihi)*mpbase;
    else
      res := integer64(hicarry)*mpbase;
      Add(res,mphihi);
      Mul(res,mpbase);
    end if;
    Add(res,mplohi);
    Mul(res,mpbase);
    Add(res,mphilo);
    Mul(res,mpbase);
    Add(res,mplolo);
    Clear(mphihi); Clear(mplohi); Clear(mphilo); Clear(mplolo);
    Clear(mpbase);
    return res;
  end Value60;

  function Value60 ( hihi,lohi,hilo,lolo : integer64;
                     verbose : boolean := true ) return Integer_Number is

    res : Integer_Number := Value60(hihi,lohi,verbose);
    inc : Integer_Number := Value60(hilo,lolo,verbose);
    base : constant integer64 := 2**60;
    mpbase60 : Integer_Number := create64(base);
    mpbase120 : Integer_Number := mpbase60*mpbase60;

  begin
    Mul(res,mpbase120);
    Add(res,inc);
    Clear(inc); Clear(mpbase60); Clear(mpbase120);
    return res;
  end Value60;

  procedure Test_Double52_Sum is

    xhi,xlo,yhi,ylo,zhi,zlo,carry : integer64;
    mpx,mpy,mpsum,mpz,err,mpzhi,mpzlo : Integer_Number;
    mpcarry,mpbase30,mpbase60,mpbase120 : Integer_Number;

  begin
    new_line;
    put_line("-> testing sum of double integers, base 2^52 ...");
    Random_Double52_Integer(xhi,xlo);
    Random_Double52_Integer(yhi,ylo);
    mpx := Value52(xhi,xlo);
    mpy := Value52(yhi,ylo);
    put("->   x : "); put(mpx); new_line;
    put("->   y : "); put(mpy); new_line;
    mpsum := mpx + mpy;
    put("-> x+y : "); put(mpsum); new_line;
    Add60(xhi,xlo,yhi,ylo,zhi,zlo,carry);
    mpz := Value52(zhi,zlo);
    mpbase30 := Multprec_Integer64_Numbers.Create32(integer32(2**30));
    mpbase60 := mpbase30*mpbase30;
    mpbase120 := mpbase60*mpbase60;
    mpzlo := Rmd(mpsum,mpbase60);
    mpzhi := mpsum - mpzlo;
    Div(mpzhi,mpbase60);
    put("->   zlo : "); put(zlo); new_line;
    put("-> mpzlo : "); put(mpzlo); new_line;
    put("->   zhi : "); put(zhi); new_line;
    put("-> mpzhi : "); put(mpzhi); new_line;
    if carry > 0 then
      mpcarry := Multprec_Integer64_Numbers.create32(integer32(carry));
      Multprec_Integer64_Numbers.Mul(mpcarry,mpbase120);
      Multprec_Integer64_Numbers.Add(mpz,mpcarry);
      Clear(mpcarry); Clear(mpbase30); Clear(mpbase60); Clear(mpbase120);
    end if;
    put("->   z : "); put(mpz); new_line;
    put("->   c : "); put(carry); new_line;
    err := mpsum - mpz;
    put("-> err : "); put(err); new_line;
    Clear(mpx); Clear(mpy); Clear(mpsum); Clear(mpz); Clear(err);
  end Test_Double52_Sum;

  procedure Test_Double60_Sum is

    xhi,xlo,yhi,ylo,zhi,zlo,carry : integer64;
    mpx,mpy,mpsum,mpz,err,mpzhi,mpzlo : Integer_Number;
    mpcarry,mpbase30,mpbase60,mpbase120 : Integer_Number;

  begin
    new_line;
    put_line("-> testing sum of double integers, base 2^60 ...");
    Random_Double64_Integer(xhi,xlo);
    Random_Double64_Integer(yhi,ylo);
    mpx := Value60(xhi,xlo,false);
    mpy := Value60(yhi,ylo,false);
    put("->   x : "); put(mpx); new_line;
    put("->   y : "); put(mpy); new_line;
    mpsum := mpx + mpy;
    put("-> x+y : "); put(mpsum); new_line;
    Add60(xhi,xlo,yhi,ylo,zhi,zlo,carry);
    mpz := Value60(zhi,zlo,false);
    mpbase30 := Multprec_Integer64_Numbers.Create32(integer32(2**30));
    mpbase60 := mpbase30*mpbase30;
    mpbase120 := mpbase60*mpbase60;
    mpzlo := Rmd(mpsum,mpbase60);
    mpzhi := mpsum - mpzlo;
    Div(mpzhi,mpbase60);
    put("->   zlo : "); put(zlo); new_line;
    put("-> mpzlo : "); put(mpzlo); new_line;
    put("->   zhi : "); put(zhi); new_line;
    put("-> mpzhi : "); put(mpzhi); new_line;
    if carry > 0 then
      mpcarry := Multprec_Integer64_Numbers.create32(integer32(carry));
      Multprec_Integer64_Numbers.Mul(mpcarry,mpbase120);
      Multprec_Integer64_Numbers.Add(mpz,mpcarry);
      Clear(mpcarry);
    end if;
    Clear(mpbase30); Clear(mpbase60); Clear(mpbase120);
    put("->   z : "); put(mpz); new_line;
    put("->   c : "); put(carry); new_line;
    err := mpsum - mpz;
    put("-> err : "); put(err); new_line;
    Clear(mpx); Clear(mpy); Clear(mpsum); Clear(mpz); Clear(err);
  end Test_Double60_Sum;

  procedure Test_Product is

    x,y,zhi,zlo,carry : integer64;
    mpx,mpy,mpprd,mpz,err : Integer_Number;
    sx,sy : integer32;

  begin
    new_line;
    put_line("-> testing product of two 60-bit integers ...");
    Random_Double64_Integer(x,y);
   -- x := 1152921504606846975; -- largest number of 60 bits
   -- y := 1152921504606846975;
    x := x/4; -- ensure the number of bits is 60 or less
    sx := Bits_of_Integers.Bit_Size(x);
    y := y/4;
    sy := Bits_of_Integers.Bit_Size(y);
    mpx := Value60(0,x,false);
    mpy := Value60(0,y,false);
    put("->   x : "); put(mpx); put(" size : "); put(sx); new_line;
    put("->   y : "); put(mpy); put(" size : "); put(sy); new_line;
    mpprd := mpx * mpy;
    put("-> x*y : "); put(mpprd); new_line;
    Mul60(x,y,zhi,zlo,carry);
    mpz := Value60(zhi,zlo,false);
    put("->   z : "); put(mpz); new_line;
    err := mpprd - mpz;
    put("-> err : "); put(err); new_line;
    Clear(mpx); Clear(mpy); Clear(mpz);
  end Test_Product;

  procedure Test_Double_Product is

    xhi,xlo,yhi,ylo,zhihi,zlohi,zhilo,zlolo,carry : integer64;
    mpx,mpy,mpprd,mpz,err : Integer_Number;

  begin
    new_line;
    put_line("-> testing product of two double integers ...");
    Random_Double64_Integer(xhi,xlo);
    xhi := xhi/4; xlo := xlo/4;
    Random_Double64_Integer(yhi,ylo);
    yhi := yhi/4; ylo := ylo/4;
   -- xhi := 1152921504606846975; -- largest number of 60 bits
   -- xlo := 1152921504606846975; 
   -- yhi := 1152921504606846975;
   -- ylo := 1152921504606846975;
    mpx := Value60(xhi,xlo,false);
    mpy := Value60(yhi,ylo,false);
    put("->   x : "); put(mpx); new_line;
    put("->   y : "); put(mpy); new_line;
    mpprd := mpx * mpy;
    put("-> x*y : "); put(mpprd); new_line;
    Dbl_Mul60(xhi,xlo,yhi,ylo,zhihi,zlohi,zhilo,zlolo,carry);
    mpz := Value60(zhihi,zlohi,zhilo,zlolo,false);
    put("-> x*y : "); put(mpprd); new_line;
    put("->   z : "); put(mpz); new_line;
    err := mpprd - mpz;
    put("-> err : "); put(err); new_line;
    Clear(mpx); Clear(mpy); Clear(mpz);
  end Test_Double_Product;

  procedure Main is
  begin
    Test_Double52_Sum;
    Test_Double60_Sum;
    Test_Product;
    Test_Double_Product;
  end Main;

end Test_Double_Integers;
