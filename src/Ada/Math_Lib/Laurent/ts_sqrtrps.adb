with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Double_Real_Power_Series;
with Test_Real_Powered_Series;

procedure ts_sqrtrps is

-- DESCRIPTION :
--   Tests the encapsulated arithmetic on real power series
--   via a square root computation.

  function Random_Series
             ( size : integer32; vrblvl : integer32 := 0 )
             return Double_Real_Power_Series.Link_to_Series is

   -- DESCRIPTION :
   --   Returns a random series of the given size.

    res : Double_Real_Power_Series.Link_to_Series;
    cff : Standard_Complex_Vectors.Vector(0..size);
    pwt : Standard_Floating_Vectors.Vector(1..size);

  begin
    if vrblvl > 0 then
      put("-> generating a series of size "); put(size,1); put_line(" ...");
    end if;
    Test_Real_Powered_Series.Random_Series(size,cff,pwt);
    if vrblvl > 0
     then Test_Real_Powered_Series.Write(cff,pwt);
    end if;
    res := Double_Real_Power_Series.make(cff,pwt);
    return res;
  end Random_Series;

  procedure Write ( s : in Double_Real_Power_Series.Link_to_Series ) is

  -- DESCRIPTION :
  --   Writes coefficients and powers of the series s.

  begin
    Test_Real_Powered_Series.Write(s.cff,s.pwt);
  end Write;

  procedure Test_SQRT ( size : in integer32 ) is

    x : Double_Real_Power_Series.Link_to_Series := Random_Series(size);

  begin
    Write(x);
  end Test_SQRT;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the size of a real power series,
  --   generates a random series, and then applies Newton's method
  --   to compute its square root.

    size : integer32 := 0;

  begin
    new_line;
    put("Give the size of the series : "); get(size);
    Test_SQRT(size);
  end Main;

begin
  Main;
end ts_sqrtrps;
