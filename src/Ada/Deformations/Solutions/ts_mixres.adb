with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Random_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Standard_Mixed_Residuals;

procedure ts_mixres is

-- DESCRIPTION :
--   Test the computation of mixed residuals.

  procedure Standard_Random_Test
              ( dim : in integer32; deg,nbr : in natural32 ) is

  -- DESCRIPTION :
  --   Runs a test on random polynomials and vectors
  --   in standard double precision.

  -- ON ENTRY :
  --   dim     number of variables in the random polynomial;
  --   deg     degree of the random polynomial;
  --   nbr     number of monomials in the random polynomial.

    rpt : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    apt : constant Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Mixed_Residuals.AbsVal(rpt);
    nvr : constant natural32 := natural32(dim);
    pol : constant Standard_Complex_Polynomials.Poly
        := Standard_Random_Polynomials.Random_Sparse_Poly(nvr,deg,nbr,0);
    abp : constant Standard_Complex_Polynomials.Poly
        := Standard_Mixed_Residuals.AbsVal(pol);
    res : double_float;

  begin
    put_line("a random point :"); put_line(rpt);
    put_line("its absolute values :"); put_line(apt);
    Symbol_Table.Init(nvr);
    put_line("a random polynomial :"); put(pol); new_line;
    put_line("its absolute version :"); put(abp); new_line;
    res := Standard_Mixed_Residuals.Residual(pol,abp,rpt);
    put("The mixed residual : "); put(res); new_line;
  end Standard_Random_Test;

  procedure Standard_Test
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Computes the mixed residuals of all solution points at p.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    abp : constant Poly_Sys(p'range)
        := Standard_Mixed_Residuals.AbsVal(p.all);
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    res : double_float;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      res := Standard_Mixed_Residuals.Residual(p.all,abp,ls.v);
      put("residual of solution "); put(i,1);
      put(" :"); put(res,3); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a dimension
  --   and then launches the tests.

    ans : character;
    dim : integer32 := 0;
    deg,nbr : natural32 := 0;
    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    put("Run a random test ? ");
    Ask_yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put("Give the dimension : "); get(dim);
      put("Give the largest degree : "); get(deg);
      put("Give the number of monomials : "); get(nbr);
      Standard_Random_Test(dim,deg,nbr);
    else
      new_line;
      put_line("Reading a system and its solutions ...");
      Standard_System_and_Solutions_io.get(p,sols);
      Standard_Test(p,sols);
    end if;
  end Main;

begin
  Main;
end ts_mixres;
