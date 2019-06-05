with Standard_Linear_Poly_Solvers;
with DoblDobl_Linear_Poly_Solvers;
with QuadDobl_Linear_Poly_Solvers;
with Standard_Scaling;
with DoblDobl_Scaling;
with QuadDobl_Scaling;
with Black_Box_Univariate_Solvers;       use Black_Box_Univariate_Solvers;
with Black_Box_Simplex_Solvers;          use Black_Box_Simplex_Solvers;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_Blackbox_Continuations;    use DoblDobl_Blackbox_Continuations;
with QuadDobl_Blackbox_Continuations;    use QuadDobl_Blackbox_Continuations;
with Black_Box_Solver_Cases;             use Black_Box_Solver_Cases;

package body Black_Box_Solvers is

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 1,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,false,rc,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        Standard_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation(deflate,pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 2,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,false,rc,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        DoblDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 3,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,false,rc,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        QuadDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 4,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        Standard_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation(deflate,pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 5,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        DoblDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 6,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        QuadDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation(pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 7,");
      put_line("for polynomial systems in double precision ...");
    end if;
    if p'first = p'last then
      n := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(Standard_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      Standard_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          Standard_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(file,0,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            Standard_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (file,deflate,pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 8,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    if p'first = p'last then
      n := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(DoblDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      DoblDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          DoblDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(file,0,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            DoblDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(file,pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 9,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    if p'first = p'last then
      n := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(QuadDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      QuadDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          QuadDobl_Complex_Poly_Systems.Copy(p,pp);
          Black_Box_Root_Counting(file,0,pp,false,rc,q,sols,sols0,roco,hoco);
          if rc /= 0 then
            QuadDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation(file,pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 10,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 11,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 12,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 13,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,rc,rocos,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 14,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,rc,rocos,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 15,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,rc,rocos,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 16,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(file,pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 17,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(file,pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 18,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco);
      if rc /= 0
       then Black_Box_Polynomial_Continuation(file,pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 19,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then
      Standard_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco);
      else
        Black_Box_Root_Counting
          (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco);
      end if;
      if rc /= 0 then
        Standard_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (deflate,integer32(nt),pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 20,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        DoblDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 21,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco);
      if rc /= 0 then
        QuadDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 22,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then
      Standard_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      end if;
      if rc /= 0 then
        Standard_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (deflate,integer32(nt),pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 23,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      end if;
      if rc /= 0 then
        DoblDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 24,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco);
      end if;
      if rc /= 0 then
        QuadDobl_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,sols,sols0,poco);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    use Standard_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 25,");
      put_line("for polynomial systems in double precision ...");
    end if;
    if p'first = p'last then
      n := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(Standard_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      Standard_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          Standard_Complex_Poly_Systems.Copy(p,pp);
          if nt >= 2 then
            Pipelined_Root_Counting
              (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          end if;
          if rc /= 0 then
            Standard_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (file,deflate,integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    use DoblDobl_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 26,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    if p'first = p'last then
      n := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(DoblDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      DoblDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          DoblDobl_Complex_Poly_Systems.Copy(p,pp);
          if nt >= 2 then
            Pipelined_Root_Counting
              (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          end if;
          if rc /= 0 then
            DoblDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (file,integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    use QuadDobl_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 27,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    if p'first = p'last then
      n := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      if n > 0 then
        rc := natural32(QuadDobl_Complex_Polynomials.Degree(p(p'first)));
        if n = 1
         then Black_Box_Durand_Kerner(p(p'first),sols);
        end if;
      end if;
    else
      QuadDobl_Linear_Poly_Solvers.Solve(p,s,fail);
      if not fail then
        rc := 1;
        ls := new Solution'(s);
        Construct(ls,sols);
      else
        if Are_Constants_In(p) then
          Black_Box_Simplex_Solver(p,sols,fail);
          fail := (fail or (Length_Of(sols) = 0));
        else
          fail := true;
        end if;
        if not fail then
          rc := Length_Of(sols);
        else
          QuadDobl_Complex_Poly_Systems.Copy(p,pp);
          if nt >= 2 then
            Pipelined_Root_Counting
              (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          else
            Black_Box_Root_Counting
              (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco);
          end if;
          if rc /= 0 then
            QuadDobl_Scaling.Scale(pp);
            Black_Box_Polynomial_Continuation
              (file,integer32(nt),pp,q,sols,sols0,poco);
            Push(sols0,sols);
          end if;
        end if;
      end if;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 28,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco);
      else
        Black_Box_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 29,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco);
      else
        Black_Box_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 30,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco);
      else
        Black_Box_Root_Counting(integer32(nt),silent,pp,rc,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 31,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(integer32(nt),pp,rc,rocos,q,sols,roco);
      else
        Black_Box_Root_Counting(integer32(nt),pp,rc,rocos,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 32,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(integer32(nt),pp,rc,rocos,q,sols,roco);
      else
        Black_Box_Root_Counting(integer32(nt),pp,rc,rocos,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 33,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(integer32(nt),pp,rc,rocos,q,sols,roco);
      else
        Black_Box_Root_Counting(integer32(nt),pp,rc,rocos,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 34,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco);
      else
        Black_Box_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(file,integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 35,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco);
      else
        Black_Box_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(file,integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 36,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco);
      else
        Black_Box_Root_Counting(file,integer32(nt),pp,rc,q,sols,roco,hoco);
      end if;
      if rc /= 0 then
        Black_Box_Polynomial_Continuation(file,integer32(nt),pp,q,sols,poco);
      end if;
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

end Black_Box_Solvers;
