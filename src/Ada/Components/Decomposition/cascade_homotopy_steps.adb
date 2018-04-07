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
with Standard_Solution_Filters;          use Standard_Solution_Filters;
with DoblDobl_Solution_Filters;          use DoblDobl_Solution_Filters;
with QuadDobl_Solution_Filters;          use QuadDobl_Solution_Filters;
with Standard_Solution_Manipulators;     use Standard_Solution_Manipulators;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with DoblDobl_Solution_Manipulators;     use DoblDobl_Solution_Manipulators;
with DoblDobl_Solution_Splitters;        use DoblDobl_Solution_Splitters;
with QuadDobl_Solution_Manipulators;     use QuadDobl_Solution_Manipulators;
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
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;
    deflate : constant boolean := false; -- := (level = 1);
    -- deflate only at the last stage, when removing the last slack variable
    -- deflation must wait after the homotopy membership test!

  begin
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation(deflate,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (deflate,integer32(nt),target,embsys,sols,time);
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then
      Filter_and_Split_Solutions
        (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (sols,integer32(n),integer32(level)-1,zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;
    deflate : constant boolean := false; -- (level = 1);
    -- deflate only at the last stage, when removing the last slack variable
    -- deflate must wait after the homotopy membership test

  begin
    new_line(file);
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(sols),1);
    put(file," paths at level "); put(file,level,1);
    put_line(file," in the cascade ...");
    flush(file);
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,deflate,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,deflate,integer32(nt),target,embsys,sols,time);
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (file,sols,integer32(n),integer32(level)-1,
     --    zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
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
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (sols,integer32(n),integer32(level)-1,zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List;
                sols0,sols1 : out Standard_Complex_Solutions.Solution_List;
                time : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    new_line(file);
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(sols),1);
    put(file," paths at level "); put(file,level,1);
    put_line(file," in the cascade ...");
    flush(file);
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (file,sols,integer32(n),integer32(level)-1,
     --    zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
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
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (sols,integer32(n),integer32(level)-1,zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    new_line(file);
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(sols),1);
    put(file," paths at level "); put(file,level,1);
    put_line(file," in the cascade ...");
    flush(file);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (file,sols,integer32(n),integer32(level)-1,
     --    zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
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
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (sols,integer32(n),integer32(level)-1,zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out DoblDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    new_line(file);
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(sols),1);
    put(file," paths at level "); put(file,level,1);
    put_line(file," in the cascade ...");
    flush(file);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (file,sols,integer32(n),integer32(level)-1,
     --    zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
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
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (sols,integer32(n),integer32(level)-1,zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    new_line(file);
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(sols),1);
    put(file," paths at level "); put(file,level,1);
    put_line(file," in the cascade ...");
    flush(file);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (file,sols,integer32(n),integer32(level)-1,
     --    zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
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
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (sols,integer32(n),integer32(level)-1,zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32; zerotol,tolsing : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                sols0,sols1 : out QuadDobl_Complex_Solutions.Solution_List;
                time : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Witness_Sets.Remove_Slice(embsys);
    n : constant natural32 := natural32(embsys'last)-level;

  begin
    new_line(file);
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    new_line(file);
    put(file,"Tracking "); put(file,Length_Of(sols),1);
    put(file," paths at level "); put(file,level,1);
    put_line(file," in the cascade ...");
    flush(file);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,time);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,time);
      Remove_Imaginary_Target(sols);
    end if;
    if level > 1 then -- at least one slack variable left
      Filter_and_Split_Solutions
        (file,sols,integer32(n),integer32(level)-1,zerotol,sols0,sols1);
     -- Zero_Singular_Split_Filter
     --   (file,sols,integer32(n),integer32(level)-1,
     --    zerotol,tolsing,sols0,sols1);
    else
      sols0 := Vanishing_Filter(sols,zerotol);
    end if;
  end Down_Continuation;

end Cascade_Homotopy_Steps;
