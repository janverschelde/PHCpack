with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Integer_Numbers_io;        use Multprec_Integer_Numbers_io;
with Standard_Random_Numbers;
with Bits_of_Integers;
with Double_Integer_Arithmetic;          use Double_Integer_Arithmetic;

package body Test_Double_Integers is

  maxint : constant integer64 := 2**62;

  procedure Random_Double_Integer ( hi,lo : out integer64 ) is
  begin
    hi := Standard_Random_Numbers.Random(0,maxint);
    lo := Standard_Random_Numbers.Random(0,maxint);
  end Random_Double_Integer;

  function Value ( hi,lo : integer64; verbose : boolean := true )
                 return Integer_Number is

    result : Integer_Number;
    hihi,lohi,hilo,lolo,hicarry,locarry : integer32;
    mphihi,mplohi,mphilo,mplolo : Integer_Number;
    base : constant integer32 := 2**30;
    mpbase : Integer_Number := create(base);

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
    mphihi := Create(hihi);
    mplohi := Create(lohi);
    mphilo := Create(hilo);
    mplolo := Create(lolo);
    if hicarry = 0 then
      result := hihi*mpbase;
    else
      result := hicarry*mpbase;
      Add(result,mphihi);
      Mul(result,mpbase);
    end if;
    Add(result,mplohi);
    Mul(result,mpbase);
    Add(result,mphilo);
    Mul(result,mpbase);
    Add(result,mplolo);
    Clear(mphihi); Clear(mplohi); Clear(mphilo); Clear(mplolo);
    Clear(mpbase);
    return result;
  end Value;

  procedure Test_Double_Sum is

    xhi,xlo,yhi,ylo,zhi,zlo,carry : integer64;
    mpx,mpy,mpsum,mpz,err,mpzhi,mpzlo : Integer_Number;
    mpcarry,mpbase30,mpbase60,mpbase120 : Integer_Number;

  begin
    put_line("-> testing sum of double integers ...");
    Random_Double_Integer(xhi,xlo);
    Random_Double_Integer(yhi,ylo);
    mpx := Value(xhi,xlo,false);
    mpy := Value(yhi,ylo,false);
    put("->   x : "); put(mpx); new_line;
    put("->   y : "); put(mpy); new_line;
    mpsum := mpx + mpy;
    put("-> x+y : "); put(mpsum); new_line;
    Add(xhi,xlo,yhi,ylo,zhi,zlo,carry);
    mpz := Value(zhi,zlo,false);
    mpbase30 := Multprec_Integer_Numbers.Create(integer32(2**30));
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
      mpcarry := Multprec_Integer_Numbers.Create(integer32(carry));
      Multprec_Integer_Numbers.Mul(mpcarry,mpbase120);
      Multprec_Integer_Numbers.Add(mpz,mpcarry);
      Clear(mpcarry); Clear(mpbase30); Clear(mpbase60); Clear(mpbase120);
    end if;
    put("->   z : "); put(mpz); new_line;
    put("->   c : "); put(carry); new_line;
    err := mpsum - mpz;
    put("-> err : "); put(err); new_line;
    Clear(mpx); Clear(mpy); Clear(mpsum); Clear(mpz); Clear(err);
  end Test_Double_Sum;

  procedure Test_Double_Product is

    x,y,zhi,zlo,carry : integer64;
    mpx,mpy,mpprd,mpz,err : Integer_Number;
    sx,sy : integer32;

  begin
    put_line("-> testing product of double integers ...");
    Random_Double_Integer(x,y);
   -- x := 1152921504606846975; -- largest number of 60 bits
   -- y := 1152921504606846975;
    x := x/4; -- ensure the number of bits is 60 or less
    sx := Bits_of_Integers.Bit_Size(x);
    y := y/4;
    sy := Bits_of_Integers.Bit_Size(y);
    mpx := Value(0,x,false);
    mpy := Value(0,y,false);
    put("->   x : "); put(mpx); put(" size : "); put(sx); new_line;
    put("->   y : "); put(mpy); put(" size : "); put(sy); new_line;
    mpprd := mpx * mpy;
    put("-> x*y : "); put(mpprd); new_line;
    Mul(x,y,zhi,zlo,carry);
    mpz := Value(zhi,zlo,false);
    put("->   z : "); put(mpz); new_line;
    err := mpprd - mpz;
    put("-> err : "); put(err); new_line;
    Clear(mpx); Clear(mpy); Clear(mpz);
  end Test_Double_Product;

  procedure Main is
  begin
    Test_Double_Sum;
    Test_Double_Product;
  end Main;

end Test_Double_Integers;
