with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
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

    dim,tdx,nbr : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    put("Give the truncation index : "); get(tdx);
    put("Give the number of monomials : "); get(nbr);
    Test_Random_Make(dim,tdx,nbr);
  end Main;

end Test_Double_rpSeries_PolySystems;
