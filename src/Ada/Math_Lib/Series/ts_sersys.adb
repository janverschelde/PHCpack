with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;
with Standard_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Series_and_Polynomials;
with Series_and_Polynomials_io;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with Standard_Series_Jaco_Matrices;

procedure ts_sersys is

-- DESCRIPTION :
--   Tests the methods on systems of series polynomials.

  procedure Test_Evaluation
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                idx : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a vector of series
  --   and evaluates the system p.
  --   The idx is the index of the variable used as series variable.

    use Standard_Series_Polynomials;

    n : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    x : Standard_Dense_Series_Vectors.Link_to_Vector;
    v : integer32 := 1;
    y : Standard_Dense_Series_Vectors.Vector(p'range);
    order : integer32 := 0;

  begin
    if idx = 0 then
      Symbol_Table.Enlarge(1);
      v := n+1;
    end if;
    new_line;
    put("Reading a vector of "); put(n,1); put_line(" series ...");
    Series_and_Polynomials_io.get(x,v);
    put_line("The vector x of series :");
    Series_and_Polynomials_io.put(x.all);
    new_line;
    put("Give the order of the evaluation : "); get(order);
    Series_and_Polynomials.Set_Order(x.all,order);
    put_line("Evaluating the series ...");
    y := Standard_Series_Poly_SysFun.Eval(p,x.all);
    put_line("The value of the system at x :");
    Series_and_Polynomials_io.put(y);
  end Test_Evaluation;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a system of series polynomials.

    ls : Standard_Series_Poly_Systems.Link_to_Poly_Sys;
    ix : integer32 := 0;

  begin
    new_line;
    put("Give the index of the series variable : "); get(ix);
    new_line;
    Series_and_Polynomials_io.get(ls,ix);
    new_line;
    put_line("The polynomial system : ");
    Series_and_Polynomials_io.put(ls.all,ix);
    Test_Evaluation(ls.all,ix);
  end Main;

begin
  Main;
end ts_sersys;
