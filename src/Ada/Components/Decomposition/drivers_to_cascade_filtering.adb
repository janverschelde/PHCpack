with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Laur_Poly_Convertors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with DoblDobl_Solution_Splitters;        use DoblDobl_Solution_Splitters;
with QuadDobl_Solution_Splitters;        use QuadDobl_Solution_Splitters;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Prompt_for_Systems;
with Standard_BlackBox_Continuations;    use Standard_BlackBox_Continuations;
with DoblDobl_BlackBox_Continuations;    use DoblDobl_BlackBox_Continuations;
with QuadDobl_BlackBox_Continuations;    use QuadDobl_BlackBox_Continuations;
with Black_Box_Solvers;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Square_and_Embed_Systems;           use Square_and_Embed_Systems;
with Continuation_Parameters;
with Drivers_for_Poly_Continuation;      use Drivers_for_Poly_Continuation;
with Path_Counts_Table;                  use Path_Counts_Table;
with Add_and_Remove_Embedding;           use Add_and_Remove_Embedding;
with Greeting_Banners;
with Write_Seed_Number;

package body Drivers_to_Cascade_Filtering is

  procedure Write_Witness_Points
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Witness_Points;

  procedure Write_Witness_Points
              ( file : in file_type;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Witness_Points;

  procedure Write_Witness_Points
              ( file : in file_type;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Witness_Points;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,pocotime);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,pocotime);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                level : in natural32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    target : constant Standard_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(0.0));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,pocotime);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,pocotime);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,pocotime);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,pocotime);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    target : constant DoblDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,pocotime);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,pocotime);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                level : in natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,pocotime);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,pocotime);
    end if;
  end Down_Continuation;

  procedure Down_Continuation
              ( file : in file_type; nt : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                level : in natural32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pocotime : out duration ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    target : constant QuadDobl_Complex_Laur_Systems.Laur_Sys(embsys'range)
           := Remove_Slice(embsys);

  begin
    put(file,"START SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,embsys);
    new_line(file);
    put(file,"TARGET SYSTEM at level "); put(file,level,1);
    put_line(file," :"); put_line(file,target);
    Set_Continuation_Parameter(sols,Create(integer(0)));
    if nt = 0 then
      Black_Box_Polynomial_Continuation
        (file,target,embsys,sols,pocotime);
    else
      Black_Box_Polynomial_Continuation
        (file,integer32(nt),target,embsys,sols,pocotime);
    end if;
  end Down_Continuation;

  procedure Witness_Generate
              ( outfile,resfile : in file_type; nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type; nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(resfile);
        put_line(resfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Points(resfile,rsols0);
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  function Append_ck ( name : string; k : natural32 ) return string is

  -- DESCRIPTION :
  --   Appends "_swk" to the name.

    nbk : constant string := Convert(integer32(k));

  begin
     return name & "_sw" & nbk;
  end Append_ck;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
   -- put_line("Writing to file" & filename);
    create(file,out_file,filename);
    put_line(file,p);
    Write_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim;
    pocotime : duration;
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse 1..integer32(topdim) loop
        Down_Continuation(outfile,nt,embsys(i).all,natural32(i),wsols,pocotime);
        Clear(sols0); Clear(sols1);
        Filter_and_Split_Solutions
          (outfile,wsols,integer32(n),i-1,zerotol,sols0,sols1);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
      Clear(sols0);
    end if;
    tstop(timer);
    Write_Path_Counts(outfile,pathcnts);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Standard_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Standard_Laur_Poly_Convertors;
    use Standard_Complex_Solutions;

    infile,outfile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      Standard_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         Standard_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,1.0E-8);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,1.0E-8);
    end if;
  end Standard_Witness_Generate;

  procedure DoblDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Laur_Poly_Convertors;
    use DoblDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      DoblDobl_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         DoblDobl_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,1.0E-8);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,1.0E-8);
    end if;
  end DoblDobl_Witness_Generate;

  procedure QuadDobl_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Laur_Poly_Convertors;
    use QuadDobl_Complex_Solutions;

    infile,outfile : file_type;
    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    filename : Link_to_String;

  begin
    if inpname /= "" then
      Open_Input_File(infile,inpname);
      QuadDobl_Read_Embedding(infile,lq,sols,dim);
      filename := new string'(inpname);
    end if;
    if lq = null then
       new_line;
       declare
         name : constant string := Read_String;
       begin
         Open_Input_File(infile,name);
         QuadDobl_Read_Embedding(infile,lq,sols,dim);
         filename := new string'(name);
       end;
    end if;
    Create_Output_File(outfile,outname);
    new_line;
    Continuation_Parameters.Tune(0);
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    put_line("See the input and output file for results ...");
    new_line;
    close(infile);
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Witness_Generate(filename.all,outfile,nt,lq.all,sols,dim,1.0E-8);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Witness_Generate(filename.all,outfile,nt,lp.all,sols,dim,1.0E-8);
    end if;
  end QuadDobl_Witness_Generate;

  procedure Driver_to_Witness_Generate
              ( nt : in natural32; inpname,outname : in string ) is

    p : constant character := Prompt_for_Precision;

  begin
    case p is
      when '0' => Standard_Witness_Generate(nt,inpname,outname);
      when '1' => DoblDobl_Witness_Generate(nt,inpname,outname);
      when '2' => QuadDobl_Witness_Generate(nt,inpname,outname);
      when others => null;
    end case;
  end Driver_to_Witness_Generate;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    ep : Link_to_Poly_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    ep : Link_to_Poly_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    ep : Link_to_Poly_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    ep : Link_to_Laur_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    ep : Link_to_Laur_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( file : in file_type; name : in string; nt : in natural32; 
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    ep : Link_to_Laur_Sys;
    timer : timing_widget;
    sols : Solution_List;
    tol : constant double_float := 1.0E-8;
    k,rc : natural32;
                
  begin
    new_line;
    Interactive_Square_and_Embed(file,p,ep,k);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    if nt = 0 then
      tstart(timer);
      Black_Box_Solvers.Solve(file,ep.all,rc,sols);
      tstop(timer);
    else
      tstart(timer);
      Black_Box_Solvers.Solve(file,nt,ep.all,rc,sols);
      tstop(timer);
    end if;
    new_line(file);
    print_times(file,timer,"calling the blackbox solver for the top");
    if not Is_Null(sols) then
      Witness_Generate(name,file,nt,ep.all,sols,k,tol);
    end if;
  end QuadDobl_Embed_and_Cascade;

  procedure Standard_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    use Standard_Laur_Poly_Convertors;

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      Standard_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all);
    else
      lp := new Standard_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      Standard_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Standard_Embed_and_Cascade;

  procedure DoblDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    use DoblDobl_Laur_Poly_Convertors;

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      DoblDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all);
    else
      lp := new DoblDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      DoblDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end DoblDobl_Embed_and_Cascade;

  procedure QuadDobl_Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    use QuadDobl_Laur_Poly_Convertors;

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    infile,outfile : file_type;
    sysonfile : boolean;
    outfilename : Link_to_String;

  begin
    Prompt_for_Systems.Read_System(infile,inpname,lq,sysonfile);
    Create_Output_File(outfile,outname,outfilename);
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lq.all) then
      QuadDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lq.all);
    else
      lp := new QuadDobl_Complex_Poly_Systems.Poly_Sys'
                  (Positive_Laurent_Polynomial_System(lq.all));
      QuadDobl_Embed_and_Cascade(outfile,outfilename.all,nt,lp.all);
    end if;
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end QuadDobl_Embed_and_Cascade;

  procedure Embed_and_Cascade
              ( nt : in natural32; inpname,outname : in string ) is

    prc : constant character := Prompt_for_Precision;


  begin
    case prc is
      when '0' => Standard_Embed_and_Cascade(nt,inpname,outname);
      when '1' => DoblDobl_Embed_and_Cascade(nt,inpname,outname);
      when '2' => QuadDobl_Embed_and_Cascade(nt,inpname,outname);
      when others => null;
    end case;
  end Embed_and_Cascade;

end Drivers_to_Cascade_Filtering;
