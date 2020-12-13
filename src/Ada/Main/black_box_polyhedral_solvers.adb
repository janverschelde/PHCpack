with Standard_Random_Numbers;
with Standard_Scaling;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with Black_Box_Mixed_Volumes;

package body Black_Box_Polyhedral_Solvers is

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 1 ...");
    end if;
    Standard_Complex_Poly_Systems.Copy(p,pp);
    Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
      (silent,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
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
      put_line("-> in black_box_polyhedral_solvers.Solve 2 ...");
    end if;
    Solve(p,silent,deflate,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
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
 
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 3 ...");
    end if;
    Standard_Complex_Poly_Systems.Copy(p,pp);
    Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
      (pp,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
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
      put_line("-> in black_box_polyhedral_solvers.Solve 4 ...");
    end if;
    Solve(p,deflate,rc,rocos,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    gamma : out Standard_Complex_Numbers.Complex_Number;
                    q : out Standard_Complex_Poly_Systems.Poly_Sys;
                    qsols : out Standard_Complex_Solutions.Solution_List;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 5 ...");
    end if;
    if nt < 2 then
      Solve(p,silent,deflate,rc,gamma,q,qsols,sols,verbose);
    else
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
        (nt,silent,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
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
                    silent,deflate : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is
 
    q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    qsols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 6 ...");
    end if;
    if nt < 2
     then Solve(p,silent,deflate,rc,gamma,q,qsols,sols,verbose);
     else Solve(nt,p,silent,deflate,rc,gamma,q,qsols,sols,verbose);
    end if;
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
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
 
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 7 ...");
    end if;
    if nt < 2 then
      Solve(p,deflate,rc,rocos,gamma,q,qsols,sols,verbose);
    else
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
        (nt,pp,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
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
      put_line("-> in black_box_polyhedral_solvers.Solve 8 ...");
    end if;
    if nt < 2
     then Solve(p,deflate,rc,rocos,gamma,q,qsols,sols,verbose);
     else Solve(nt,p,deflate,rc,rocos,gamma,q,qsols,sols,verbose);
    end if;
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
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
 
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 9 ...");
    end if;
    Standard_Complex_Poly_Systems.Copy(p,pp);
    Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
      (file,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
    if rc /= 0 then
      Copy(sols,qsols);
      if not Is_Null(sols0) then
        Copy(sols0,qsols0);
        Push(qsols0,qsols);
      end if;
      Standard_Scaling.Scale(pp);
      gamma := Standard_Random_Numbers.Random1;
      Black_Box_Polynomial_Continuation
        (file,deflate,pp,q,sols,sols0,poco,verbose-1);
      Push(sols0,sols);
    end if;
    Standard_Complex_Poly_Systems.Clear(pp);
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
      put_line("-> in black_box_polyhedral_solvers.Solve 10 ...");
    end if;
    Solve(file,p,deflate,rc,gamma,q,qsols,sols,verbose);
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
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
 
    pp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0,qsols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 11 ...");
    end if;
    if nt < 2 then
      Solve(file,p,deflate,rc,gamma,q,qsols,sols,verbose);
    else
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
        (file,nt,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
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
      put_line("-> in black_box_polyhedral_solvers.Solve 12 ...");
    end if;
    if nt < 2
     then Solve(file,p,deflate,rc,gamma,q,qsols,sols,verbose);
     else Solve(file,nt,p,deflate,rc,gamma,q,qsols,sols,verbose);
    end if;
    if rc /= 0 then
      Standard_Complex_Solutions.Deep_Clear(qsols);
      Standard_Complex_Poly_Systems.Clear(q);
    end if;
  end Solve;

end Black_Box_Polyhedral_Solvers;
