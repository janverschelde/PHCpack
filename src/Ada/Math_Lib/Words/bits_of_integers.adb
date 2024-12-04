package body Bits_of_Integers is

  type unsigned_integer64 is mod 2**integer64'size;

  procedure Split ( nbr : in integer64; high,low : out integer64 ) is

    mask : constant unsigned_integer64 := 4294967295;
    unbr : constant unsigned_integer64 := unsigned_integer64(nbr);

  begin
    low := integer64(unbr and mask);
    high := nbr - low;
  end Split;

  procedure Quarter ( nbr : in integer64;
                      hihi,lohi,hilo,lolo : out integer64 ) is

    high,low : integer64;
    lowmask : constant unsigned_integer64 := 65535;
    highmask : constant unsigned_integer64 := 281474976710655;
    ulow,uhigh : unsigned_integer64;

  begin
    Split(nbr,high,low);
    ulow := unsigned_integer64(low);
    lolo := integer64(ulow and lowmask);
    hilo := low - lolo;
    uhigh := unsigned_integer64(high);
    lohi := integer64(uhigh and highmask);
    hihi := high - lohi;
  end Quarter;

end Bits_of_Integers;
