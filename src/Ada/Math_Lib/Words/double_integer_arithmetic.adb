with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Bits_of_Integers;                   use Bits_of_Integers;

package body Double_Integer_Arithmetic is

  maxint60 : constant unsigned_integer64 := 2**60;

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
        put("-> carry : "); put(integer64(ucarry/maxint60)); new_line;
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
        put("-> carry : "); put(carry); new_line;
      end if;
    end if;
  end Add;

  procedure Mul ( xhi,xlo,yhi,ylo : in integer64; 
                  zhi,zlo,carry : out integer64 ) is
  begin
    zhi := 0;
    zlo := 0;
    carry := 0;
  end Mul;

end Double_Integer_Arithmetic;
