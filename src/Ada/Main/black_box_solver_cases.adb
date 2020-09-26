with text_io;                            use text_io;
with Standard_Linear_Poly_Solvers;
with DoblDobl_Linear_Poly_Solvers;
with QuadDobl_Linear_Poly_Solvers;
with Black_Box_Univariate_Solvers;       use Black_Box_Univariate_Solvers;
with Black_Box_Simplex_Solvers;          use Black_Box_Simplex_Solvers;
with Black_Box_Helpers;

package body Black_Box_Solver_Cases is

  procedure Solve_for_Special_Cases
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                rc : out natural32;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solver_cases.Solve_for_Special_Cases 1 ...");
    end if;
    if p'first = p'last then
      n := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(Standard_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols); fail := false;
         else fail := true;
        end if;
      end if;
    else
      Standard_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Black_Box_Helpers.Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
      end if;
    end if;
    if fail
     then rc := 0;
     else rc := Length_Of(sols);
    end if;
  end Solve_for_Special_Cases;

  procedure Solve_for_Special_Cases
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                rc : out natural32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solver_cases.Solve_for_Special_Cases 2 ...");
    end if;
    if p'first = p'last then
      n := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(DoblDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols); fail := false;
         else fail := true;
        end if;
      end if;
    else
      DoblDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Black_Box_Helpers.Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
      end if;
    end if;
    if fail
     then rc := 0;
     else rc := Length_Of(sols);
    end if;
  end Solve_for_Special_Cases;

  procedure Solve_for_Special_Cases
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                rc : out natural32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solver_cases.Solve_for_Special_Cases 3 ...");
    end if;
    if p'first = p'last then
      n := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(QuadDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols); fail := false;
         else fail := true;
        end if;
      end if;
    else
      QuadDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Black_Box_Helpers.Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
      end if;
    end if;
    if fail
     then rc := 0;
     else rc := Length_Of(sols);
    end if;
  end Solve_for_Special_Cases;

end Black_Box_Solver_Cases;
