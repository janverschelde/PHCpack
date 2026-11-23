with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Communications_with_User;
with Standard_Integer_Vectors_IO;
with Standard_Random_Numbers;
with Standard_Random_Vectors;
with Double_Real_Power_Series;
with Double_rpSeries_Vectors;
with Test_Double_rpSeries_Arithmetic;

package body Test_Double_rpSeries_Monomials is

  function Random_Monomial
             ( dim,tdx,low,upp : integer32;
               incpwr : double_float := 0.0 ) return Monomial is

    res : Monomial(dim,tdx);

  begin
    res.cff := Test_Double_rpSeries_Arithmetic.Random_Series(tdx,incpwr);
    res.deg := Standard_Random_Vectors.Random_Vector(1,dim,low,upp);
    return res;
  end Random_Monomial;

  function Random_Monomial
             ( dim,tdx,low,upp : integer32;
               incpwr : double_float := 0.0 ) return Link_to_Monomial is

    res_rep : constant Monomial(dim,tdx)
            := Random_Monomial(dim,tdx,low,upp,incpwr);
    res : constant Link_to_Monomial := new Monomial'(res_rep);

  begin
    return res;
  end Random_Monomial;

  function Random_Monomials
             ( nbr,dim,tdx,low,upp : integer32;
               incpwr : double_float := 0.0 ) return Monomial_Array is

    res : Monomial_Array(1..nbr);

  begin
    for i in 1..nbr loop
      res(i) := Random_Monomial(dim,tdx,low,upp,incpwr);
    end loop;
    return res;
  end Random_Monomials;

  procedure Write ( mon : in Monomial ) is
  begin
    Test_Double_rpSeries_Arithmetic.Write(mon.cff);
    Standard_Integer_Vectors_IO.put(mon.deg); new_line;
  end Write;

  procedure Write ( mon : in Link_to_Monomial ) is
  begin
    if mon /= null
     then Write(mon.all);
    end if;
  end Write;

  procedure Write ( mons : in Monomial_Array ) is
  begin
    for i in mons'range loop
      Write(mons(i));
    end loop;
  end Write;

  procedure Test_Gradient ( mon : in Link_to_Monomial ) is
  begin
    for i in mon.deg'range loop
      put("-> component "); put(i,1); put_line(" of the gradient :");
      declare
        difmon : Link_to_Monomial := Diff(mon,i);
      begin
        Write(difmon); Clear(difmon);
      end;
    end loop;
  end Test_Gradient;

  procedure Test_Random_Make ( dim,tdx : in integer32 ) is

    mon : constant Link_to_Monomial := Random_Monomial(dim,tdx,0,5);

  begin
    put_line("-> a random coefficient :");
    Test_Double_rpSeries_Arithmetic.Write(mon.cff);
    put("-> random exponents :");
    Standard_Integer_Vectors_IO.put(mon.deg); new_line;
    Test_Gradient(mon);
  end Test_Random_Make;

  procedure Test_Random_Eval ( dim,tdx : in integer32 ) is

    use Test_Double_rpSeries_Arithmetic;

    mon : Link_to_Monomial := Random_Monomial(dim,tdx,0,5,1.0);
    x : Double_rpSeries_Vectors.Vector(1..dim)
      := Random_Series_Vector(1,dim,tdx,1.0);
    y : Double_Real_Power_Series.Link_to_Series := Eval(mon,x);

  begin
    put_line("-> a random coefficient :");
    Test_Double_rpSeries_Arithmetic.Write(mon.cff);
    put("-> random exponents :");
    Standard_Integer_Vectors_IO.put(mon.deg); new_line;
    put_line("a random vector :"); Write(x);
    put_line("evaluated at the monomial :"); Write(y);
    clear(mon);
    Double_rpSeries_Vectors.clear(x);
    Double_Real_Power_Series.clear(y);
  end Test_Random_Eval;

  procedure Main is

    dim,tdx,seed : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the seed (0 for any seed) : "); get(seed);
    if seed = 0
     then seed := Standard_Random_Numbers.get_Seed;
     else Standard_Random_Numbers.Set_Seed(natural32(seed));
    end if;
    put("Give the dimension : "); get(dim);
    put("Give the truncation index : "); get(tdx);
    new_line;
    put_line("MENU to test monomials with series coefficients :");
    put_line("  1. make a random monomial");
    put_line("  2. evaluate a random monomial");
    put("Type 1 or 2 to select a test : ");
    Communications_with_User.Ask_Alternative(ans,"12");
    new_line;
    case ans is
      when '1' => Test_Random_Make(dim,tdx);
      when '2' => Test_Random_Eval(dim,tdx);
      when others => null;
    end case;
    new_line;
    put("seed used : "); put(seed,1); new_line;
  end Main;

end Test_Double_rpSeries_Monomials;
