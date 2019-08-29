with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Solution_Splitters;        use Standard_Solution_Splitters;
with DoblDobl_Solution_Splitters;        use DoblDobl_Solution_Splitters;
with QuadDobl_Solution_Splitters;        use QuadDobl_Solution_Splitters;
with Witness_Sets;                       use Witness_Sets;
with Path_Counts_Table;                  use Path_Counts_Table;
with Greeting_Banners;
with Write_Seed_Number;
with Cascade_Homotopy_Steps;             use Cascade_Homotopy_Steps;
with Cascade_Homotopies_io;              use Cascade_Homotopies_io;

package body Cascade_Homotopies is

  procedure Witness_Generate
              ( outfile : in file_type; nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
   -- nbq : constant natural32 := natural32(ep'last); -- #equations
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(outfile,ep);
    Write_Super_Witness_Points(outfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(outfile);
        put_line(outfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Super_Witness_Points(outfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(outfile,rsols0);
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
   -- new_line(outfile);
   -- Write_Seed_Number(outfile);
   -- put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile : in file_type; nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(outfile,ep);
    Write_Super_Witness_Points(outfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(outfile);
        put_line(outfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Super_Witness_Points(outfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(outfile,rsols0);
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
   -- new_line(outfile);
   -- Write_Seed_Number(outfile);
   -- put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(outfile,ep);
    Write_Super_Witness_Points(outfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(outfile);
        put_line(outfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Super_Witness_Points(outfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(outfile,rsols0);
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
   -- new_line(outfile);
   -- Write_Seed_Number(outfile);
   -- put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(outfile,ep);
    Write_Super_Witness_Points(outfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(outfile);
        put_line(outfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Super_Witness_Points(outfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(outfile,rsols0);
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
   -- new_line(outfile);
   -- Write_Seed_Number(outfile);
   -- put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(outfile,ep);
    Write_Super_Witness_Points(outfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(outfile);
        put_line(outfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Super_Witness_Points(outfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(outfile,rsols0);
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
   -- new_line(outfile);
   -- Write_Seed_Number(outfile);
   -- put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(outfile,ep);
    Write_Super_Witness_Points(outfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(sols0),Length_Of(sols1));
        new_line(outfile);
        put_line(outfile,embsys(i-1).all);
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Super_Witness_Points(outfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(outfile,rsols0);
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
   -- new_line(outfile);
   -- Write_Seed_Number(outfile);
   -- put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( outfile,resfile : in file_type; nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
   -- nbq : constant natural32 := natural32(ep'last); -- #equations
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Super_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
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
              Write_Super_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(resfile,rsols0);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Super_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
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
              Write_Super_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(resfile,rsols0);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Super_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
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
              Write_Super_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(resfile,rsols0);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Super_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
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
              Write_Super_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(resfile,rsols0);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Super_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
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
              Write_Super_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(resfile,rsols0);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
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
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    put_line(resfile,ep);
    Write_Super_Witness_Points(resfile,sols0);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
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
              Write_Super_Witness_Points(resfile,rsols1);
            end;
          end if;
        elsif not Is_Null(sols0) then
          declare
            rsols0 : constant Solution_List := Remove_Component(sols0);
          begin
            Write_Super_Witness_Points(resfile,rsols0);
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
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    Write_Witness_Superset(name,ep,esols0(integer32(topdim)),topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,esols0(i-1),sols1,pocotime);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        if i = 1 then
          if not Is_Null(sols1) then
            declare
              rsols1 : constant Solution_List := Remove_Component(sols1);
            begin
              Write_Witness_Superset(name,embsys(i-1).all,rsols1,0);
            end;
          end if;
        elsif not Is_Null(esols0(i-1)) then
          declare
            rsols0 : constant Solution_List := Remove_Component(esols0(i-1));
          begin
            Write_Witness_Superset(name,embsys(i-1).all,rsols0,natural32(i-1));
          end;
        end if;
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
    Write_Path_Counts(outfile,pathcnts,times,alltime);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( name : in string; outfile : in file_type;
                nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);
    times : Array_of_Duration(0..integer(topdim));
    pocotime,alltime : duration;

  begin
    tstart(timer);
    times := (times'range => 0.0);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        times(integer(i)) := pocotime;
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
    alltime := Elapsed_User_Time(timer);
    Write_Path_Counts(outfile,pathcnts,times,alltime);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);
    times : Array_of_Duration(0..integer(topdim));
    pocotime,alltime : duration;

  begin
    tstart(timer);
    times := (times'range => 0.0);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        times(integer(i)) := pocotime;
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
    alltime := Elapsed_User_Time(timer);
    Write_Path_Counts(outfile,pathcnts,times,alltime);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);
    times : Array_of_Duration(0..integer(topdim));
    pocotime,alltime : duration;

  begin
    tstart(timer);
    times := (times'range => 0.0);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        times(integer(i)) := pocotime;
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
    alltime := Elapsed_User_Time(timer);
    Write_Path_Counts(outfile,pathcnts,times,alltime);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);
    times : Array_of_Duration(0..integer(topdim));
    pocotime,alltime : duration;

  begin
    tstart(timer);
    times := (times'range => 0.0);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        times(integer(i)) := pocotime;
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
    alltime := Elapsed_User_Time(timer);
    Write_Path_Counts(outfile,pathcnts,times,alltime);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    embsys : Array_of_Poly_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);
    times : Array_of_Duration(0..integer(topdim));
    pocotime,alltime : duration;

  begin
    tstart(timer);
    times := (times'range => 0.0);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        times(integer(i)) := pocotime;
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
    alltime := Elapsed_User_Time(timer);
    Write_Path_Counts(outfile,pathcnts,times,alltime);
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
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    embsys : Array_of_Laur_Sys(0..integer32(topdim));
    pathcnts : Standard_Natural_VecVecs.VecVec(embsys'range);
    times : Array_of_Duration(0..integer(topdim));
    pocotime,alltime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (outfile,sols,integer32(n),integer32(topdim),zerotol,sols0,sols1);
   -- Zero_Singular_Split_Filter
   --   (outfile,sols,integer32(n),integer32(topdim),
   --    zerotol,tolsing,sols0,sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(sols0),Length_Of(sols1));
    Write_Witness_Superset(name,ep,sols0,topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols0); Clear(sols1);
        Down_Continuation
          (outfile,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        times(integer(i)) := pocotime;
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
    alltime := Elapsed_User_Time(timer);
    Write_Path_Counts(outfile,pathcnts,times,alltime);
    new_line(outfile);
    print_times(outfile,timer,"Witness Generate with Cascade of Homotopies");
    new_line(outfile);
    Write_Seed_Number(outfile);
    put_line(outfile,Greeting_Banners.Version);
  end Witness_Generate;

  procedure Witness_Generate
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
   -- Zero_Singular_Split_Filter
   --   (file,sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (file,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (file,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (file,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (file,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (file,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

  procedure Witness_Generate
              ( file : in file_type; nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (file,nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate;

-- OUTPUT WITH CALLBACK PROCEDURE :

  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    Report_Witness_Set
      (embsys(integer32(topdim)).all,esols0(integer32(topdim)),topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Report_Witness_Set(embsys(i-1).all,esols0(i-1),natural32(i-1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate_Callback;

  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    Report_Witness_Set
      (embsys(integer32(topdim)).all,esols0(integer32(topdim)),topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Report_Witness_Set(embsys(i-1).all,esols0(i-1),natural32(i-1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate_Callback;

  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    Report_Witness_Set
      (embsys(integer32(topdim)).all,esols0(integer32(topdim)),topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Report_Witness_Set(embsys(i-1).all,esols0(i-1),natural32(i-1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate_Callback;

  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    Report_Witness_Set
      (embsys(integer32(topdim)).all,esols0(integer32(topdim)),topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Report_Witness_Set(embsys(i-1).all,esols0(i-1),natural32(i-1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate_Callback;

  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Poly_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Poly_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    Report_Witness_Set
      (embsys(integer32(topdim)).all,esols0(integer32(topdim)),topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Report_Witness_Set(embsys(i-1).all,esols0(i-1),natural32(i-1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate_Callback;

  procedure Witness_Generate_Callback
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32;
                zerotol,tolsing : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    timer : Timing_Widget;
    wsols,sols0,sols1 : Solution_List;
    n : constant natural32 := natural32(ep'last)-topdim; -- #variables
    pocotime : duration;

  begin
    tstart(timer);
    embsys(integer32(topdim)) := new Laur_Sys'(ep);
    for i in 0..topdim-1 loop
      embsys(integer32(i)) := new Laur_Sys'(Remove_Embedding1(ep,topdim-i));
    end loop;
    Filter_and_Split_Solutions
      (sols,integer32(n),integer32(topdim),zerotol,
       esols0(integer32(topdim)),sols1);
   -- Zero_Singular_Split_Filter
   --   (sols,integer32(n),integer32(topdim),zerotol,tolsing,
   --    esols0(integer32(topdim)),sols1);
    Update_Path_Counts
      (pathcnts,topdim,Length_Of(sols),Length_Of(esols0(integer32(topdim))),
       Length_Of(sols1));
    Report_Witness_Set
      (embsys(integer32(topdim)).all,esols0(integer32(topdim)),topdim);
    if not Is_Null(sols1) then
      Copy(sols1,wsols);
      for i in reverse integer32(lowdim)+1..integer32(topdim) loop
        Clear(sols1);
        Down_Continuation
          (nt,embsys(i).all,natural32(i),zerotol,tolsing,
           wsols,sols0,sols1,pocotime);
        esols0(i-1) := Remove_Component(sols0); Clear(sols0);
        times(integer(i)) := pocotime;
        Update_Path_Counts
          (pathcnts,natural32(i-1),
           Length_Of(wsols),Length_Of(esols0(i-1)),Length_Of(sols1));
        Report_Witness_Set(embsys(i-1).all,esols0(i-1),natural32(i-1));
        Clear(wsols);
        exit when Is_Null(sols1);
        wsols := Remove_Component(sols1);
      end loop;
    end if;
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Generate_Callback;

end Cascade_Homotopies;
