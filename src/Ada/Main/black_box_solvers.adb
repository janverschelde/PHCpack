with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Linear_Poly_Solvers;
with DoblDobl_Linear_Poly_Solvers;
with QuadDobl_Linear_Poly_Solvers;
with Standard_Scaling;
with DoblDobl_Scaling;
with QuadDobl_Scaling;
with Black_Box_Univariate_Solvers;       use Black_Box_Univariate_Solvers;
with Black_Box_Simplex_Solvers;          use Black_Box_Simplex_Solvers;
with Standard_Monomial_Maps;
with Standard_Monomial_Maps_io;          use Standard_Monomial_Maps_io;
with DoblDobl_Monomial_Maps;
with QuadDobl_Monomial_Maps;
with Black_Box_Binomial_Solvers;         use Black_Box_Binomial_Solvers;
with Black_Box_Factorization;            use Black_Box_Factorization;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_Blackbox_Continuations;    use DoblDobl_Blackbox_Continuations;
with QuadDobl_Blackbox_Continuations;    use QuadDobl_Blackbox_Continuations;
with Greeting_Banners;
with Black_Box_Solver_Cases;             use Black_Box_Solver_Cases;

package body Black_Box_Solvers is

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;
    deflate : boolean := true;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;
    deflate : boolean := true;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;
    deflate : boolean := true;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;
    deflate : boolean := true;

  begin
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
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;
    deflate : boolean := true;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    use Standard_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;
    deflate : boolean := true;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    use DoblDobl_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    use QuadDobl_Complex_Solutions;

    fail : boolean;
    ls : Link_to_Solution;
    s : Solution(p'last);
    n : natural32;
    pp,q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
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
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    Black_Box_Simplex_Solver(p,sols,fail);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 0 then
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out Standard_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out DoblDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
                    sols : out QuadDobl_Complex_Solutions.Solution_List ) is
 
    fail : boolean;
    pp,q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
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
