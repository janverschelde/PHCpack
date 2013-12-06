with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Linear_Product_System;
with Black_Box_Solvers;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Tropisms_of_Supports;

procedure ts_degcurve is

-- DESCRIPTION :
--   Attempt to compute the degree of a curve via polyhedral methods.

  function Add_Hyperplane ( p : Poly_Sys ) return Poly_Sys is

  -- DESCRIPTION :
  --   Adds a random hyperplane to the system p.

    res : Poly_Sys(p'first..p'last+1);
    n : constant natural32 := Number_of_Unknowns(p(p'first));
    cff : constant Standard_Complex_Vectors.Vector(0..integer32(n))
        := Standard_Random_Vectors.Random_Vector(0,integer32(n));
    hyp : constant Poly := Standard_Linear_Product_System.Polynomial(cff);

  begin
    res(p'range) := p;
    res(res'last) := hyp;
    return res;
  end Add_Hyperplane;

  function Degree ( p : Poly_Sys) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the curve by solving the system
  --   augmented with a random hyperplane.

  -- REQUIRED : p'last + 1 = number of variables.

    res : natural32;
    ph : Poly_Sys(p'first..p'last+1) := Add_Hyperplane(p);
    rc : natural32;
    sols : Solution_List;

  begin
    Black_Box_Solvers.Solve(ph,true,rc,sols);
    res := Length_Of(sols);
    Clear(sols);
    Clear(ph(ph'last));
    return res;
  end Degree;

  procedure Compute_Tropisms ( p : in Poly_Sys ) is

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    s : Array_of_Lists(p'range) := Create(p);
    e : Array_of_Faces(s'range); 

  begin
    put_line("The supports :"); put(s);
    e := Tropisms_of_Supports.Edges(n,s);
    Tropisms_of_Supports.Show_Parallel_Edges(n,s,e);
    Tropisms_of_Supports.Sort_Parallel_Edges(n,s,e);
  end Compute_Tropisms;

  procedure Main is

    lp : Link_to_Poly_Sys;
   -- d : natural;

  begin
    new_line;
   -- put_line("Reading a system of n equations in n+1 unknowns ...");
    put_line("Reading a polynomial system ...");
    new_line;
    get(lp);
    put_line("Your system : "); put(lp.all);
   -- d := Degree(lp.all);
   -- put("Degree of the curve : "); put(d,1); new_line;
    Compute_Tropisms(lp.all);
  end Main;

begin
  Main;
end ts_degcurve;
