with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Cayley_Trick;                       use Cayley_Trick;
with Driver_for_Minkowski_Polynomials;

procedure ts_drivmink is

-- DESCRIPTION :
--   This procedure tests the computation of the Minkowski polynomial.

  lp : Link_to_Poly_Sys;

begin
  new_line;
  put_line("Interactive testing of power lists.");
  new_line;
  get(lp);
  declare
    supports : Array_of_Lists(lp'range) := Create(lp.all);
    lifted : Array_of_Lists(supports'range);
    n : constant integer32 := lp'last;
    mix,perms : Link_to_Vector;
    t : Triangulation;
    mixsub : Mixed_Subdivision;
  begin
    put_line("The supports of the system : "); put(supports);
    Compute_Mixture(supports,mix,perms);
    Dynamic_Cayley(n,mix.all,supports,false,true,0,lifted,t);
    Driver_for_Minkowski_Polynomials(Standard_Output,n,mix.all,t,true,mixsub);
  end;
end ts_drivmink;
