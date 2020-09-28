with Standard_Scaling;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with Black_Box_Mixed_Volumes;

package body Black_Box_Polyhedral_Solvers is

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 1 ...");
    end if;
    Standard_Complex_Poly_Systems.Copy(p,pp);
    Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
      (silent,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
    if rc /= 0 then
      Standard_Scaling.Scale(pp);
      Black_Box_Polynomial_Continuation
        (deflate,pp,q,sols,sols0,poco,verbose-1);
      Push(sols0,sols);
    end if;
  end Solve;

  procedure Solve ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean;
                    rc : out natural32; rocos : out Link_to_String;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 2 ...");
    end if;
    Standard_Complex_Poly_Systems.Copy(p,pp);
    Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
      (pp,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
    if rc /= 0 then
      Standard_Scaling.Scale(pp);
      Black_Box_Polynomial_Continuation
        (deflate,pp,q,sols,sols0,poco,verbose-1);
      Push(sols0,sols);
    end if;
  end Solve;

  procedure Solve ( nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    silent,deflate : in boolean;
                    rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 3 ...");
    end if;
    if nt < 2 then
      Solve(p,silent,deflate,rc,sols,verbose);
    else
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
        (nt,silent,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Standard_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (deflate,integer32(nt),pp,q,sols,sols0,poco,verbose-1);
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
 
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 4 ...");
    end if;
    if nt < 2 then
      Solve(p,deflate,rc,rocos,sols,verbose);
    else
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
        (nt,pp,rc,rocos,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Standard_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (deflate,integer32(nt),pp,q,sols,sols0,poco,verbose-1);
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
 
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 5 ...");
    end if;
    Standard_Complex_Poly_Systems.Copy(p,pp);
    Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
      (file,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
    if rc /= 0 then
      Standard_Scaling.Scale(pp);
      Black_Box_Polynomial_Continuation
        (file,deflate,pp,q,sols,sols0,poco,verbose-1);
      Push(sols0,sols);
    end if;
  end Solve;

  procedure Solve ( file : in file_type; nt : in natural32;
                    p : in Standard_Complex_Poly_Systems.Poly_Sys;
                    deflate : in boolean; rc : out natural32;
                    sols : out Standard_Complex_Solutions.Solution_List;
                    verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;
 
    pp,q : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    roco,hoco,poco : duration;
    sols0 : Solution_List;

  begin
    if verbose > 0 then
      put_line("-> in black_box_polyhedral_solvers.Solve 6 ...");
    end if;
    if nt < 2 then
      Solve(file,p,deflate,rc,sols,verbose);
    else
      Standard_Complex_Poly_Systems.Copy(p,pp);
      Black_Box_Mixed_Volumes.Black_Box_Polyhedral_Homotopies
        (file,nt,pp,rc,q,sols,sols0,roco,hoco,verbose-1);
      if rc /= 0 then
        Standard_Scaling.Scale(pp);
        Black_Box_Polynomial_Continuation
          (file,deflate,integer32(nt),pp,q,sols,sols0,poco,verbose-1);
        Push(sols0,sols);
      end if;
    end if;
  end Solve;

end Black_Box_Polyhedral_Solvers;
