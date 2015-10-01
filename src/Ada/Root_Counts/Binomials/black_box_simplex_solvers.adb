with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;             use QuadDobl_Root_Refiners;
with Standard_Simpomial_Solvers;
with DoblDobl_Simpomial_Solvers;
with QuadDobl_Simpomial_Solvers;

package body Black_Box_Simplex_Solvers is

  procedure Black_Box_Simplex_Solver
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use Standard_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    numit,max : natural32;
    tol_zero : constant double_float := 1.0E-12;
    zero_y : boolean;
    deflate : boolean := false;

  begin
    Standard_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if (not fail and then (Length_Of(sols) > 0)) then
      epsxa := 1.0E-12;
      epsfa := 1.0E-12;
      tolsing := 1.0E-08;
      max := 5;
      numit := 0;
      Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max,deflate);
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use DoblDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    numit,max : natural32;
    tol_zero : constant double_double := create(1.0E-24);
    zero_y : boolean;

  begin
    DoblDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if (not fail and then (Length_Of(sols) > 0)) then
      epsxa := 1.0E-24;
      epsfa := 1.0E-24;
      tolsing := 1.0E-16;
      max := 5;
      numit := 0;
      Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max);
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use QuadDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    numit,max : natural32;
    tol_zero : constant quad_double := create(1.0E-48);
    zero_y : boolean;

  begin
    QuadDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if (not fail and then (Length_Of(sols) > 0)) then
      epsxa := 1.0E-48;
      epsfa := 1.0E-48;
      tolsing := 1.0E-24;
      max := 5;
      numit := 0;
      Silent_Root_Refiner(p,sols,epsxa,epsfa,tolsing,numit,max);
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use Standard_Complex_Solutions;

    tol_zero : constant double_float := 1.0E-12;
    zero_y : boolean;

  begin
    Standard_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use DoblDobl_Complex_Solutions;

    tol_zero : constant double_double := create(1.0E-24);
    zero_y : boolean;

  begin
    DoblDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use QuadDobl_Complex_Solutions;

    tol_zero : constant quad_double := create(1.0E-48);
    zero_y : boolean;

  begin
    QuadDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use Standard_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    numit,max : natural32;
    tol_zero : constant double_float := 1.0E-12;
    zero_y : boolean;
    deflate : boolean := false;

  begin
    Standard_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if not fail then
      if zero_y then
        new_line(file);
        put_line(file,
    "The system has no solutions with all components different from zero.");
        put_line(file,
    "Try perturbing the constant terms, solve the perturbed system, and");
        put_line(file,
    "remove the perturbations then again via homotopy continuation.");
      elsif Length_Of(sols) > 0 then
        epsxa := 1.0E-12;
        epsfa := 1.0E-12;
        tolsing := 1.0E-08;
        max := 5;
        numit := 0;
        Reporting_Root_Refiner
          (file,p,sols,epsxa,epsfa,tolsing,numit,max,deflate,false);
      end if;
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use DoblDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    numit,max : natural32;
    tol_zero : constant double_double := create(1.0E-24);
    zero_y : boolean;
    deflate : boolean := false;

  begin
    DoblDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if not fail then
      if zero_y then
        new_line(file);
        put_line(file,
    "The system has no solutions with all components different from zero.");
        put_line(file,
    "Try perturbing the constant terms, solve the perturbed system, and");
        put_line(file,
    "remove the perturbations then again via homotopy continuation.");
      elsif Length_Of(sols) > 0 then
        epsxa := 1.0E-24;
        epsfa := 1.0E-24;
        tolsing := 1.0E-16;
        max := 5;
        numit := 0;
        Reporting_Root_Refiner(file,p,sols,epsxa,epsfa,tolsing,numit,max,false);
      end if;
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use QuadDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    numit,max : natural32;
    tol_zero : constant quad_double := create(1.0E-48);
    zero_y : boolean;
    deflate : boolean := false;

  begin
    QuadDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if not fail then
      if zero_y then
        new_line(file);
        put_line(file,
    "The system has no solutions with all components different from zero.");
        put_line(file,
    "Try perturbing the constant terms, solve the perturbed system, and");
        put_line(file,
    "remove the perturbations then again via homotopy continuation.");
      elsif Length_Of(sols) > 0 then
        epsxa := 1.0E-48;
        epsfa := 1.0E-48;
        tolsing := 1.0E-24;
        max := 5;
        numit := 0;
        Reporting_Root_Refiner(file,p,sols,epsxa,epsfa,tolsing,numit,max,false);
      end if;
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use Standard_Complex_Solutions;

    tol_zero : constant double_float := 1.0E-12;
    zero_y : boolean;
    epsxa,epsfa,tolsing : double_float;
    numit,max : natural32;
    ref_sols : Solution_List;

  begin
    Standard_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if not fail then
      if Length_Of(sols) > 0 then
        epsxa := 1.0E-12;
        epsfa := 1.0E-12;
        tolsing := 1.0E-8;
        max := 5; numit := 0;
        Reporting_Root_Refiner(file,p,sols,ref_sols,
          epsxa,epsfa,tolsing,numit,max,false);
      end if;
    end if;
    if fail
     then put_line(file,"Blackbox simpomial solver reports failure.");
    end if;
    if zero_y
     then put_line(file,"Blackbox simpomial solver reports zero component.");
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use DoblDobl_Complex_Solutions;

    tol_zero : constant double_double := create(1.0E-24);
    zero_y : boolean;

  begin
    DoblDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if fail
     then put_line(file,"Blackbox simpomial solver reports failure.");
    end if;
    if zero_y
     then put_line(file,"Blackbox simpomial solver reports zero component.");
    end if;
  end Black_Box_Simplex_Solver;

  procedure Black_Box_Simplex_Solver
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                fail : out boolean ) is

    use QuadDobl_Complex_Solutions;

    tol_zero : constant quad_double := create(1.0E-48);
    zero_y : boolean;

  begin
    QuadDobl_Simpomial_Solvers.Solve(p,tol_zero,sols,fail,zero_y);
    if fail
     then put_line(file,"Blackbox simpomial solver reports failure.");
    end if;
    if zero_y
     then put_line(file,"Blackbox simpomial solver reports zero component.");
    end if;
  end Black_Box_Simplex_Solver;

end Black_Box_Simplex_Solvers;
