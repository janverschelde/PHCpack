with Ada.Text_IO;                       use Ada.Text_IO;
with Communications_with_User;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Integer_Vectors;
with Double_Real_Power_Series;
with Double_rpSeries_Vectors;
with Double_rpSeries_PolyFunctions;
with Test_Double_rpSeries_Arithmetic;
with Test_Double_rpSeries_Monomials;

package body Test_Double_rpSeries_Polynomials is

  procedure Add ( p : in out Poly; mon : in Monomial ) is
  
    trm : Term;

  begin
    trm.dg := new Standard_Integer_Vectors.Vector'(mon.deg);
    trm.cf := new Double_Real_Power_Series.Series'(mon.cff);
    Add(p,trm);
    Clear(trm);
  end Add;

  procedure Add ( p : in out Poly; mon : in Link_to_Monomial ) is
  begin
    if mon /= null
     then Add(p,mon.all);
    end if;
  end Add;

  function Make ( mons : Monomial_Array ) return Poly is

    res : Poly := Null_Poly;

  begin
    for i in mons'range loop
      Add(res,mons(i));
    end loop;
    return res;
  end Make;

  function Make ( p : Poly ) return Monomial_Array is

    res : Monomial_Array(1..integer32(Number_of_Terms(p)));
    idx : integer32 := 0;

    procedure Visit_Term ( trm : in Term; continue : out boolean ) is
    begin
      idx := idx + 1;
      declare
        cff : constant Double_Real_Power_Series.Series := trm.cf.all;
        deg : constant Standard_Integer_Vectors.Vector := trm.dg.all;
      begin
        res(idx) := Make(cff,deg);
      end;
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Make;

  procedure Write ( p : in Poly ) is

    procedure Write_Term ( trm : in Term; continue : out boolean ) is
    begin
      declare
        cff : constant Double_Real_Power_Series.Series := trm.cf.all;
        deg : constant Standard_Integer_Vectors.Vector := trm.dg.all;
        mon : constant Monomial := Make(cff,deg);
      begin
        Test_Double_rpSeries_Monomials.Write(mon);
      end;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
  end Write;

  function Random_Polynomial ( dim,tdx,nbr,low,upp : integer32 ) return Poly is

    use Test_Double_rpSeries_Monomials;

    mons : Monomial_Array(1..nbr) := Random_Monomials(nbr,dim,tdx,low,upp);
    res : constant Poly := Make(mons);
  
  begin
    Clear(mons);
    return res;
  end Random_Polynomial;

  procedure Test_Random_Make ( dim,tdx,nbr : in integer32 ) is

    p : Poly := Random_Polynomial(nbr,dim,tdx,0,9);
  
  begin
    put_line("a random polynomial :"); Write(p);
    Clear(p);
  end Test_Random_Make;

  procedure Test_Random_Eval ( dim,tdx,nbr : in integer32 ) is

    use Test_Double_rpSeries_Arithmetic;
    use Test_Double_rpSeries_Monomials;

    mons : Monomial_Array(1..nbr) := Random_Monomials(nbr,dim,tdx,1,9);
    p : Poly := Make(mons);
    x : Double_rpSeries_Vectors.Vector(1..nbr)
      := Random_Series_Vector(1,nbr,tdx);
    y : Double_Real_Power_Series.Link_to_Series
      := Double_rpSeries_PolyFunctions.Eval(p,x);
    z : Double_Real_Power_Series.Link_to_Series := Eval(mons,x);

  begin
    put_line("a random polynomial :"); Write(p);
    put_line("a random vector :"); Write(x);
    put_line("evaluated at the polynomial :"); Write(y);
    put_line("evaluated at the monomial array :"); Write(z);
    clear(mons); clear(p);
    Double_rpSeries_Vectors.clear(x);
    Double_Real_Power_Series.clear(y);
    Double_Real_Power_Series.clear(z);
  end Test_Random_Eval;

  procedure Main is

    dim,tdx,nbr : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the truncation index : "); get(tdx);
    put("Give the number of monomials : "); get(nbr);
    new_line;
    put_line("MENU to test polynomials with series coefficients :");
    put_line("  1. make a random polynomial");
    put_line("  2. evaluate a random polynomial");
    put("Type 1 or 2 to select a test : ");
    Communications_with_User.Ask_Alternative(ans,"12");
    new_line;
    case ans is
      when '1' => Test_Random_Make(dim,tdx,nbr);
      when '2' => Test_Random_Eval(dim,tdx,nbr);
      when others => null;
    end case;
  end Main;

end Test_Double_rpSeries_Polynomials;
