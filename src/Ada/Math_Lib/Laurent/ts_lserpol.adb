with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;
with Standard_Laurent_Series;
with Random_Laurent_Series;             use Random_Laurent_Series;
with Test_Standard_Lseries_Matrices;

procedure ts_lserpol is

-- DESCRIPTION :
--   Development of the evaluation and differentiation of a polynomial
--   at a sequence of Laurent series.

  procedure Write ( plead : in Standard_Integer_Vectors.Vector;
                    pcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                    pmons : in Standard_Integer_VecVecs.VecVec;
                    s : in string := "p" ) is
  begin
    for k in plead'range loop
      put(s & "("); put(k,1); put(") :"); put(pmons(k)); new_line;
      Standard_Laurent_Series.Write(plead(k),pcffs(k).all);
    end loop;
  end Write;

  procedure Make_Random_Polynomial
              ( dim,nbr,deg,pwr,low,upp : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec;
                mons : out Standard_Integer_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Makes a random polynomial with Laurent series as coefficients.

  -- ON ENTRY :
  --   dim      the dimension is the number of variables;
  --   nbr      number of monomials in the polynomial;
  --   deg      degree of the series;
  --   pwr      largest power for every variable;
  --   low      lower bound on leading exponents of the series;
  --   upp      upper bound on leading exponents of the series.

  -- ON RETURN :
  --   lead     an array of range 1..nbr with the leading exponents
  --            of the power series coefficients;
  --   cffs     coefficient vectors of the power series coefficients;
  --   mons     exponents of the monomials in the polynomial.

  begin
    for k in 1..nbr loop
      declare
        mon : constant Standard_Integer_Vectors.Vector(1..dim)
            := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
      begin
        mons(k) := new Standard_Integer_Vectors.Vector'(mon);
      end;
    end loop;
    Random_Vector(nbr,deg,low,upp,lead,cffs);
  end Make_Random_Polynomial;

  procedure Eval ( deg,mlead : in integer32;
                   cff : in Standard_Complex_Vectors.Link_to_Vector;
                   mon : in Standard_Integer_Vectors.Link_to_Vector;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Evaluates a monomial at a Laurent series.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   mlead    leading exponent of the Laurent series coefficient
  --            of the monomial;
  --   cff      Laurent series coefficient of the monomial;
  --   mon      exponents of the monomial;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation.

  -- ON RETURN :
  --   ye       leading exponent of the result of the evaluation;
  --   yc       coefficient vector of the value of the monomial.

    ze : integer32;
    zc : Standard_Complex_Vectors.Vector(0..deg);

  begin
    ye := mlead;
    for k in 0..deg loop -- initialize result with monomial coefficient
      yc(k) := cff(k);
    end loop;
    for i in mon'range loop -- mon(i) is the power of the i-th variable
      if mon(i) > 0 then
        ye := ye + xlead(i)*mon(i);
        for j in 1..mon(i) loop
          Standard_Laurent_Series.Multiply
            (deg,ye,xlead(i),yc,xcffs(i).all,ze,zc);
          ye := ze;
          for k in 0..deg loop
            yc(k) := zc(k);
          end loop;
        end loop;
      end if;
    end loop;
  end Eval;

  procedure Eval ( deg : in integer32;
                   plead : in Standard_Integer_Vectors.Vector;
                   pcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   pmons : in Standard_Integer_VecVecs.VecVec;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Evaluates a polynomial at a Laurent series.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   plead    leading exponents of the Laurent series coefficients;
  --   pcffs    coefficient vectors of the Laurent series coefficients;
  --   pmons    exponents of the monomials in the polynomial;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation.

  -- ON RETURN :
  --   ye       leading exponent of the result of the evaluation;
  --   yc       coefficient vector of the value of the polynomial.

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..deg);

  begin
    Eval(deg,plead(1),pcffs(1),pmons(1),xlead,xcffs,ye,yc);
    for i in 2..plead'last loop
      Eval(deg,plead(i),pcffs(i),pmons(i),xlead,xcffs,ze,zc);
      Standard_Laurent_Series.Add(deg,ye,ze,yc,zc,ewrk,cwrk);
      ye := ewrk;
      for k in 0..deg loop
        yc(k) := cwrk(k);
      end loop;
    end loop;
  end Eval;

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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the parameters of the tests and then runs tests.

    dim,nbr,deg,pwr,low,upp : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables : "); get(dim);
    put("Give the number of monomials : "); get(nbr);
    put("Give the degree of the series : "); get(deg);
    put("Give the largest power of the variables : "); get(pwr);
    put("Give the lower bound on the leading exponents : "); get(low);
    put("Give the upper bound on the leading exponents : "); get(upp);
    Test(dim,nbr,deg,pwr,low,upp);
  end Main;

begin
  Main;
end ts_lserpol;
