with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Binomial_Coefficients;              use Binomial_Coefficients;

procedure ts_bincff is

-- DESCRIPTION :
--   Test on the computation of binomial coefficients.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for n and k
  --   and prints n choose k.

    n,k,b,s : integer32 := 0;
    pw : natural32;

  begin
    put("Give n : "); get(n);
    put("Give k : "); get(k);
    b := Binomial(n,k);
    put(n,1); put(" choose "); put(k,1);
    put(" : "); put(b); new_line;
    new_line;
    for i in 0..n loop
      b := binomial(n,i);
      put(n,1); put(" choose "); put(i,1); 
      put(" : "); put(b,1); new_line;
      s := s + b;
    end loop;
    put("sum : "); put(s,1); new_line;
    put("2**"); put(n,1); put(" : ");
    pw := 2**natural(n);
    put(pw); new_line;
  end Main;

begin
  Main;
end ts_bincff;
