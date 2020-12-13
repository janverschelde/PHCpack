with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Scaling;
with DoblDobl_Scaling;
with QuadDobl_Scaling;
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
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 1,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (0,silent,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        Standard_Scaling.Scale(pp);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (deflate,pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      Standard_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 2,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve(p,silent,deflate,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 3,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (0,silent,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        DoblDobl_Scaling.Scale(pp);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      DoblDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 4,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve(p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 5,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (0,silent,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        QuadDobl_Scaling.Scale(pp);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
        QuadDobl_Complex_Poly_Systems.Clear(pp);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 6,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve(p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 7,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (0,pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        Standard_Scaling.Scale(pp);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (deflate,pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
        Standard_Complex_Poly_Systems.Clear(pp);
      end if;
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 8,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve(p,deflate,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 9,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then -- not a special case
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (0,pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        DoblDobl_Scaling.Scale(pp);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      DoblDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 10,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve(p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 11,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (0,pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        QuadDobl_Scaling.Scale(pp);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      QuadDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 12,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve(p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 13,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (file,0,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;     
        Standard_Scaling.Scale(pp);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,deflate,pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      Standard_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 14,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve(file,p,deflate,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 15,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (file,0,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        DoblDobl_Scaling.Scale(pp);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,pp,q,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      DoblDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 16,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve(file,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 17,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (file,0,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        QuadDobl_Scaling.Scale(pp);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      QuadDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 18,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve(file,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 19,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        Standard_Complex_Solutions.Copy(sols,qsols);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation(pp,q,gamma,sols,poco,verbose-1);
      end if;
      Standard_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 20,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Solve(p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 21,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        DoblDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation(pp,q,gamma,sols,poco,verbose-1);
      end if;
      DoblDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 22,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Solve(p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 23,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,silent,pp,rc,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        QuadDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation(pp,q,gamma,sols,poco,verbose-1);
      end if;
      QuadDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 24,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Solve(p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 25,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,rc,rocos,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        Standard_Complex_Solutions.Copy(sols,qsols);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation(pp,q,gamma,sols,poco,verbose-1);
      end if;
      Standard_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 26,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Solve(p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 27,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,rc,rocos,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        DoblDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation(pp,q,gamma,sols,poco,verbose-1);
      end if;
      DoblDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 28,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Solve(p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 29,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(0,pp,rc,rocos,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        QuadDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation(pp,q,gamma,sols,poco,verbose-1);
      end if;
      QuadDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 30,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Solve(p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 31,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        Standard_Complex_Solutions.Copy(sols,qsols);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,pp,q,gamma,sols,poco,verbose-1);
      end if;
      Standard_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 32,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Solve(file,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 33,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        DoblDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,pp,q,gamma,sols,poco,verbose-1);
      end if;
      DoblDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 34,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Solve(file,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 35,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      Black_Box_Root_Counting(file,0,pp,rc,q,sols,roco,hoco,verbose-1);
      if rc /= 0 then
        QuadDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,pp,q,gamma,sols,poco,verbose-1);
      end if;
      QuadDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 36,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Solve(file,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 37,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail);
    if fail then
      Standard_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        Standard_Scaling.Scale(pp);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (deflate,integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      Standard_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 38,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve(nt,p,silent,deflate,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 39,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        DoblDobl_Scaling.Scale(pp);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      DoblDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 40,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve(nt,p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 41,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Root_Counting
        (integer32(nt),silent,pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        QuadDobl_Scaling.Scale(pp);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      QuadDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 42,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve(nt,p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    fail : boolean;
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 43,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then
      Standard_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        Standard_Scaling.Scale(pp);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (deflate,integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      Standard_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 44,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve(nt,p,deflate,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;
 
    fail : boolean;
    pp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 45,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        DoblDobl_Scaling.Scale(pp);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      DoblDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 46,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve(nt,p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;
 
    fail : boolean;
    pp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 47,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,false,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        QuadDobl_Scaling.Scale(pp);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,sols0,poco);
        Push(sols0,sols);
      end if;
      QuadDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 48,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve(nt,p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    use Standard_Complex_Solutions;

    fail : boolean;
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 49,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      Standard_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      else
        Black_Box_Root_Counting
          (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        Standard_Scaling.Scale(pp);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,deflate,integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      Standard_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 50,");
      put_line("for polynomial systems in double precision ...");
    end if;
    Solve(file,nt,p,deflate,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    use DoblDobl_Complex_Solutions;

    fail : boolean;
    pp : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 51,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      DoblDobl_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      else
        Black_Box_Root_Counting
          (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(sols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        DoblDobl_Scaling.Scale(pp);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      DoblDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 52,");
      put_line("for polynomial systems in double double precision ...");
    end if;
    Solve(file,nt,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    use QuadDobl_Complex_Solutions;

    fail : boolean;
    pp : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 53,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve_for_Special_Cases(p,rc,sols,fail,verbose-1);
    if fail then -- not a special case
      QuadDobl_Complex_Poly_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      else
        Black_Box_Root_Counting
          (file,integer32(nt),pp,false,rc,q,sols,sols0,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Copy(sols,qsols);
        if not Is_Null(qsols0) then
          Copy(sols0,qsols0);
          Push(qsols0,qsols);
        end if;
        QuadDobl_Scaling.Scale(pp);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,integer32(nt),pp,q,gamma,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
      QuadDobl_Complex_Poly_Systems.Clear(pp);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 54,");
      put_line("for polynomial systems in quad double precision ...");
    end if;
    Solve(file,nt,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 55,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),silent,pp,rc,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),silent,pp,rc,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Standard_Complex_Solutions.Copy(sols,qsols);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      Standard_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 56,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Solve(nt,p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 57,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),silent,pp,rc,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),silent,pp,rc,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        DoblDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      DoblDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 58,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Solve(nt,p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 59,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),silent,pp,rc,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),silent,pp,rc,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        QuadDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      QuadDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    silent : in boolean; rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 60,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Solve(nt,p,silent,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 61,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,rc,rocos,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,rc,rocos,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Standard_Complex_Solutions.Copy(sols,qsols);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      Standard_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 62,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Solve(nt,p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 63,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,rc,rocos,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,rc,rocos,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        DoblDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      DoblDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 64,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Solve(nt,p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 65,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (integer32(nt),pp,rc,rocos,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (integer32(nt),pp,rc,rocos,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        QuadDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      QuadDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 66,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Solve(nt,p,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Laur_Systems.Laur_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 67,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := Standard_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      Standard_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (file,integer32(nt),pp,rc,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (file,integer32(nt),pp,rc,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        Standard_Complex_Solutions.Copy(sols,qsols);
        gamma := Standard_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      Standard_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 68,");
      put_line("for Laurent polynomial systems in double precision ...");
    end if;
    Solve(file,nt,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out DoblDobl_Complex_Numbers.Complex_Number;
                    q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out DoblDobl_Complex_Solutions.Solution_List;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 69,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := DoblDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      DoblDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (file,integer32(nt),pp,rc,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (file,integer32(nt),pp,rc,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        DoblDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := DoblDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      DoblDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out DoblDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 70,");
      put_line("for Laurent polynomial systems in double double precision ...");
    end if;
    Solve(file,nt,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      DoblDobl_Complex_Solutions.Deep_Clear(qsols);
      DoblDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    gamma : out QuadDobl_Complex_Numbers.Complex_Number;
                    q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    qsols : out QuadDobl_Complex_Solutions.Solution_List;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    fail : boolean;
    pp : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    roco,hoco,poco : duration;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 71,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Black_Box_Simplex_Solver(p,sols,fail,verbose-1);
    rc := QuadDobl_Complex_Solutions.Length_Of(sols);
    if fail or (rc = 0) then
      QuadDobl_Complex_Laur_Systems.Copy(p,pp);
      if nt >= 2 then
        Pipelined_Root_Counting
          (file,integer32(nt),pp,rc,q,sols,roco,verbose-1);
      else
        Black_Box_Root_Counting
          (file,integer32(nt),pp,rc,q,sols,roco,hoco,verbose-1);
      end if;
      if rc /= 0 then
        QuadDobl_Complex_Solutions.Copy(sols,qsols);
        gamma := QuadDobl_Random_Numbers.Random1;
        Black_Box_Polynomial_Continuation
          (file,integer32(nt),pp,q,gamma,sols,poco,verbose-1);
      end if;
      QuadDobl_Complex_Laur_Systems.Clear(pp);
    else
      roco := 0.0; hoco := 0.0; poco := 0.0;
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    rc : out natural32;
                    sols : out QuadDobl_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    qsols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_solvers.Solve 72,");
      put_line("for Laurent polynomial systems in quad double precision ...");
    end if;
    Solve(file,nt,p,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      QuadDobl_Complex_Solutions.Deep_Clear(qsols);
      QuadDobl_Complex_Laur_Systems.Clear(q);
    end if;
  end Solve;

end Black_Box_Solvers;
