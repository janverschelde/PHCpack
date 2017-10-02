with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with DoblDobl_Solution_Splitters;        use DoblDobl_Solution_Splitters;
with QuadDobl_Solution_Splitters;        use QuadDobl_Solution_Splitters;
with Continuation_Parameters;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_BlackBox_Continuations;    use DoblDobl_BlackBox_Continuations;
with QuadDobl_BlackBox_Continuations;    use QuadDobl_BlackBox_Continuations;
with Witness_Sets;

package body Cascade_Homotopy_Steps is

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation(target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation(target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation(target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation(target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation(target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation(target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
    end if;
    Filter_and_Split_Solutions
      (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
  end Down_Continuation;

end Cascade_Homotopy_Steps;
