with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Homotopy_Membership_Filters;

package body Cascade_Membership_Filters is

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop -- removing junk at dimension dim
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop -- removing junk at dimension dim
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Homotopy_Membership_Filters.Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

end Cascade_Membership_Filters;
