with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Random_Numbers;
with Test_Double_rpSeries_Polynomials;

package body Test_Double_rpSeries_PolySystems is

  procedure Write ( p : in Laur_Sys ) is
  begin
    for i in p'range loop
      put("polynomial "); put(i,1); put_line(" :");
      Test_Double_rpSeries_Polynomials.Write(p(i));
    end loop;
  end Write;

  procedure Test_Random_Make ( dim,tdx,nbr : in integer32 ) is

    p : Laur_Sys(1..dim);

    use Test_Double_rpSeries_Polynomials;

  begin
    for i in 1..dim loop
      p(i) := Random_Polynomial(dim,tdx,nbr,0,9);
    end loop;
    Write(p);
  end Test_Random_Make;

  procedure Main is

    seed,dim,tdx,nbr : integer32 := 0;

  begin
    new_line;
    put("Give the seed (0 for any seed) : "); get(seed);
    if seed = 0
     then seed := Standard_Random_Numbers.get_Seed;
     else Standard_Random_Numbers.Set_Seed(natural32(seed));
    end if;
    put("-> seed used : "); put(seed,1); new_line;
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the truncation index : "); get(tdx);
    put("Give the number of monomials : "); get(nbr);
    Test_Random_Make(dim,tdx,nbr);
    new_line;
    put("seed used : "); put(seed,1); new_line;
  end Main;

end Test_Double_rpSeries_PolySystems;
