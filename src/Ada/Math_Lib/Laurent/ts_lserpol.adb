with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Laurent_Series;
with Random_Laurent_Series;             use Random_Laurent_Series;
with Test_Standard_Lseries_Matrices;
with Standard_Lseries_Polynomials;      use Standard_Lseries_Polynomials;

procedure ts_lserpol is

-- DESCRIPTION :
--   Development of the evaluation and differentiation of a polynomial
--   at a sequence of Laurent series.

  procedure Test ( dim,nbr,deg,pwr,low,upp : in integer32 ) is

  -- DESCRIPTION :
  --   Generates random data and runs tests.

  -- ON ENTRY :
  --   dim      the dimension is the number of variables;
  --   nbr      number of monomials in the polynomial;
  --   deg      degree of the series;
  --   pwr      largest power for every variable;
  --   low      lower bound on leading exponents of the series;
  --   upp      upper bound on leading exponents of the series.

    plead : Standard_Integer_Vectors.Vector(1..nbr);
    pcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    pmons : Standard_Integer_VecVecs.VecVec(1..nbr);
    xlead : Standard_Integer_Vectors.Vector(1..dim);
    xcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    ye : integer32;
    yc : Standard_Complex_Vectors.Vector(0..deg);

  begin
    Make_Random_Polynomial(dim,nbr,deg,pwr,low,upp,plead,pcffs,pmons);
    put_line("A random polynomial with Laurent series coefficients :");
    Write(plead,pcffs,pmons);
    Random_Vector(dim,deg,low,upp,xlead,xcffs);
    put_line("A random vector of Laurent series :");
    Test_Standard_Lseries_Matrices.Write(xlead,xcffs,"x");
    Eval(deg,plead,pcffs,pmons,xlead,xcffs,ye,yc);
    put_line("The result of the evaluation :");
    Standard_Laurent_Series.Write(ye,yc);
  end Test;

  procedure Test_Input is

  -- DESCRIPTION :
  --   Prompts for a Laurent polynomial system and then constructs the data
  --   to evaluate the system at a vector of Laurent series.

    p : Link_to_Laur_Sys;
    neq,dim,tdx,deg : integer32 := 0;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ..."); get(p);
    new_line;
    put_line("-> your system :"); put(p.all);
    new_line;
    neq := p'last;
    dim := integer32(Number_of_Unknowns(p(p'first)));
    put("Read "); put(neq,1); put(" polynomials in "); put(dim,1);
    put(" variables :"); Symbol_Table_io.Write; new_line;
    if neq /= dim then
      for k in 1..natural32(dim) loop
        declare
          sb : constant Symbol_Table.Symbol := Symbol_Table.get(k);
        begin
          if sb(1) = 't'
           then tdx := integer32(k); exit;
          end if;
        end;
      end loop;
      put("-> index of t : "); put(tdx,1); new_line;
    end if;
    new_line;
    put("Give the degree of the series : "); get(deg);
    new_line;
    if tdx = 0 then
      for k in p'range loop
        Make_Series_Polynomial(p(k),dim,dim,0,deg);
      end loop;
    else
      for k in p'range loop
        Make_Series_Polynomial(p(k),dim,dim-1,tdx,deg);
      end loop;
    end if;
  end Test_Input;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the parameters of the tests and then runs tests.

    dim,nbr,deg,pwr,low,upp : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Generate a random polynomial ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then
      Test_Input;
    else
      new_line;
      put("Give the number of variables : "); get(dim);
      put("Give the number of monomials : "); get(nbr);
      put("Give the degree of the series : "); get(deg);
      put("Give the largest power of the variables : "); get(pwr);
      put("Give the lower bound on the leading exponents : "); get(low);
      put("Give the upper bound on the leading exponents : "); get(upp);
      Test(dim,nbr,deg,pwr,low,upp);
    end if;
  end Main;

begin
  Main;
end ts_lserpol;
