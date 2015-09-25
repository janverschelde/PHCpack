with Standard_Poly_Laur_Convertors;     use Standard_Poly_Laur_Convertors;
with Standard_Permanent_Factors;        use Standard_Permanent_Factors;
with Standard_Monomial_Map_Filters;     use Standard_Monomial_Map_Filters;

package body Black_Box_Binomial_Solvers is

  procedure Black_Box_Binomial_Solver
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pure : in boolean;
         sols : out Standard_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is

    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
      := Polynomial_to_Laurent_System(p);

  begin
    Black_Box_Binomial_Solver(q,pure,sols,fail); 
    Standard_Complex_Laur_Systems.Clear(q);
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pure : in boolean;
         sols : out DoblDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is

  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pure : in boolean;
         sols : out QuadDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is

  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                pure : in boolean;
         sols : out Standard_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is

    maps : Standard_Monomial_Maps.Monomial_Map_List;

  begin
    Silent_Affine_Solutions_with_Iterator(p,pure,maps,fail);
    if not fail
     then Silent_Filter(p,maps,sols);
    end if;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pure : in boolean;
         sols : out DoblDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is
  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pure : in boolean;
         sols : out QuadDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is
  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
         sols : out Standard_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is

    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range)
      := Polynomial_to_Laurent_System(p);

  begin
    Black_Box_Binomial_Solver(file,q,sols,fail);
    Standard_Complex_Laur_Systems.Clear(q);
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
         sols : out DoblDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is
  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
         sols : out QuadDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is
  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
         sols : out Standard_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is

    maps : Standard_Monomial_Maps.Monomial_Map_List;

  begin
    Interactive_Affine_Solutions_with_Iterator(p,maps,fail);
    if not fail
     then Reporting_Filter(file,p,maps,sols);
    end if;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
         sols : out DoblDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is
  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

  procedure Black_Box_Binomial_Solver
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
         sols : out QuadDobl_Monomial_Maps.Link_to_Array_of_Monomial_Map_Lists;
                fail : out boolean ) is
  begin
    sols := null;
    fail := true;
  end Black_Box_Binomial_Solver;

end Black_Box_Binomial_Solvers;
