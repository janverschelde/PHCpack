with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Supports_of_Polynomial_Systems;    use Supports_of_Polynomial_Systems;
with Integer_Mixed_Subdivisions;        use Integer_Mixed_Subdivisions;
with Drivers_for_Static_Lifting;        use Drivers_for_Static_Lifting;

procedure ts_puiseux is

-- DESCRIPTION :
--   Development of the Newton-Puiseux algorithm.

  procedure Tropisms ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Given a system of n Laurent polynomials in n+1 variables,
  --   computes the tropisms, where the last variable is the parameter.

    sup : Array_of_Lists(p'range) := Create(p);
    dim : constant integer32 := p'last;
    mix : Standard_Integer_Vectors.Vector(1..dim);
    mcc : Mixed_Subdivision;
    mv : natural32;

  begin
    put_line("The supports : "); put(sup);
    mix := (mix'range => 1);
    Integer_Create_Mixed_Cells(standard_output,dim,mix,false,sup,mcc);
    Integer_Volume_Computation(standard_output,dim,mix,true,sup,mcc,mv);
  end Tropisms;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a Laurent polynomial system
  --   and checks whether the n polynomials have n+1 variables.

    lp : Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    new_line;
    put_line("Reading a Laurent polynomial system ...");
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Number of polynomials : "); put(nq,1); new_line;
    put("Number of variables : "); put(nv,1); new_line;
    if nv /= nq+1 then
      put(nv,1); put(" /= "); put(nq,1); put(" + 1");
    else
      Tropisms(lp.all);
    end if;
  end Main;

begin
  Main;
end ts_puiseux;
