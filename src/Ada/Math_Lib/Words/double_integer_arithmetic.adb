with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Bits_of_Integers;                   use Bits_of_Integers;

package body Double_Integer_Arithmetic is

  maxint60 : constant unsigned_integer64 := 2**60;
  b30 : constant integer64 := 2**30;

  procedure Add ( xhi,xlo,yhi,ylo : in integer64; 
                  zhi,zlo,carry : out integer64;
                  verbose : in boolean := true ) is

    mask60 : constant unsigned_integer64 := 1152921504606846975;
    uxhi : constant unsigned_integer64 := unsigned_integer64(xhi);
    uxlo : constant unsigned_integer64 := unsigned_integer64(xlo);
    uyhi : constant unsigned_integer64 := unsigned_integer64(yhi);
    uylo : constant unsigned_integer64 := unsigned_integer64(ylo);
    uzhi,uzlo,ucarry : unsigned_integer64;

  begin
    uzlo := uxlo + uylo;
    uzhi := uxhi + uyhi;
    if verbose then
      put("-> uzlo : "); put(integer64(uzlo)); new_line;
      put("-> uzhi : "); put(integer64(uzhi)); new_line;
    end if;
    if uzlo <= maxint60 then
      zlo := integer64(uzlo);
    else
      zlo := integer64(uzlo and mask60);
      ucarry := uzlo - unsigned_integer64(zlo);
      uzhi := uzhi + ucarry/maxint60;
      if verbose then
        put("-> locarry : "); put(integer64(ucarry/maxint60)); new_line;
      end if;
    end if;
    if uzhi <= maxint60 then
      zhi := integer64(uzhi);
      carry := 0;
    else
      zhi := integer64(uzhi and mask60);
      ucarry := uzhi - unsigned_integer64(zhi);
      carry := integer64(ucarry/maxint60);
      if verbose then
        put("-> hicarry : "); put(carry); new_line;
      end if;
    end if;
  end Add;

  procedure Mul ( x,y : in integer64; 
                  zhi,zlo,carry : out integer64;
                  verbose : in boolean := true ) is

    xhi32,xlo32,xcarry32,yhi32,ylo32,ycarry32 : integer32;
    zhihi,zlohi,zhilo,zlolo : unsigned_integer64;
    phi32,plo32,pcarry32,qhi32,qlo32,qcarry32 : integer32;

  begin
    Split_30_Bit_Words(x,xhi32,xlo32,xcarry32);
    Split_30_Bit_Words(y,yhi32,ylo32,ycarry32);
    if verbose then
      put("->  xhi32 : "); put(xhi32); new_line;
      put("->  xlo32 : "); put(xlo32); new_line;
      put("-> xcarry : "); put(xcarry32); new_line;
      put("->  yhi32 : "); put(yhi32); new_line;
      put("->  ylo32 : "); put(ylo32); new_line;
      put("-> ycarry : "); put(ycarry32); new_line;
    end if;
    zlolo := unsigned_integer64(xlo32)*unsigned_integer64(ylo32);
    zlohi := unsigned_integer64(xlo32)*unsigned_integer64(yhi32);
    zhilo := unsigned_integer64(xhi32)*unsigned_integer64(ylo32);
    zhihi := unsigned_integer64(xhi32)*unsigned_integer64(yhi32);
    if verbose then
      put("-> zlolo : "); put(integer64(zlolo)); new_line;
      put("-> zlohi : "); put(integer64(zlohi)); new_line;
      put("-> zhilo : "); put(integer64(zhilo)); new_line;
      put("-> zhihi : "); put(integer64(zhihi)); new_line;
    end if;
    Split_30_Bit_Words(zlohi,phi32,plo32,pcarry32);
    Split_30_Bit_Words(zhilo,qhi32,qlo32,qcarry32);
    if verbose then
      put("->  phi32 : "); put(phi32); new_line;
      put("->  plo32 : "); put(plo32); new_line;
      put("-> pcarry : "); put(pcarry32); new_line;
      put("->  qhi32 : "); put(qhi32); new_line;
      put("->  qlo32 : "); put(qlo32); new_line;
      put("-> qcarry : "); put(qcarry32); new_line;
    end if;
    zlo := integer64(zlolo) + integer64(plo32)*b30 + integer64(qlo32)*b30;
    zhi := integer64(zhihi) + integer64(phi32) + integer64(qhi32);
    carry := 0;
  end Mul;

  procedure Dbl_Mul ( xhi,xlo,yhi,ylo : in integer64; 
                      zhihi,zlohi,zhilo,zlolo,carry : out integer64;
                      verbose : in boolean := true ) is

    ahi,alo,acr,bhi,blo,bcr,chi,clo,ccr,dhi,dlo,dcr : integer64;
    ehi,elo,ecr : integer64;

  begin
    Mul(xhi,yhi,ahi,alo,acr,verbose);
    Mul(xhi,ylo,bhi,blo,bcr,verbose);
    Mul(xlo,yhi,chi,clo,ccr,verbose);
    Mul(xlo,ylo,dhi,dlo,dcr,verbose);
    if verbose then
      put("-> ahi : "); put(ahi); new_line;
      put("-> alo : "); put(alo); new_line;
      put("-> acr : "); put(acr); new_line;
      put("-> bhi : "); put(bhi); new_line;
      put("-> blo : "); put(blo); new_line;
      put("-> bcr : "); put(bcr); new_line;
      put("-> chi : "); put(chi); new_line;
      put("-> clo : "); put(clo); new_line;
      put("-> ccr : "); put(ccr); new_line;
      put("-> dhi : "); put(dhi); new_line;
      put("-> dlo : "); put(dlo); new_line;
      put("-> dcr : "); put(dcr); new_line;
    end if;
    Add(dhi,dlo,blo,0,ehi,elo,ecr,verbose);
    Add(ehi,elo,clo,0,zhilo,zlolo,carry,verbose);
    if verbose then
      put("-> ehi : "); put(ehi); new_line;
      put("-> elo : "); put(elo); new_line;
      put("-> ecr : "); put(ecr); new_line;
      put("-> zhilo : "); put(zhilo); new_line;
      put("-> zlolo : "); put(zlolo); new_line;
      put("-> carry : "); put(carry); new_line;
    end if;
    carry := carry + ecr;
    Add(ahi,alo,0,bhi+carry,ehi,elo,ecr,verbose);
    Add(ehi,elo,0,chi,zhihi,zlohi,carry,verbose);
    if verbose then
      put("-> ehi : "); put(ehi); new_line;
      put("-> elo : "); put(elo); new_line;
      put("-> ecr : "); put(ecr); new_line;
      put("-> zhihi : "); put(zhihi); new_line;
      put("-> zlohi : "); put(zlohi); new_line;
      put("-> carry : "); put(carry); new_line;
    end if;
    carry := carry + ecr;
  end Dbl_Mul;

end Double_Integer_Arithmetic;
