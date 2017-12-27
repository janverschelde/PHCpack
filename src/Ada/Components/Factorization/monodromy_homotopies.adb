with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Witness_Sets;
with Sampling_Machine;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with Sample_Point_Lists;                use Sample_Point_Lists;
with DoblDobl_Sample_Lists;             use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Lists;             use QuadDobl_Sample_Lists;
with Rectangular_Sample_Grids;
with DoblDobl_Rectangular_Sample_Grids;
with QuadDobl_Rectangular_Sample_Grids;
with Monodromy_Component_Breakup;

package body Monodromy_Homotopies is

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    hyp : Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Witness_Sets.Slices(eqs,dim);
    sps : constant Standard_Sample_List := Create(pts,hyp);
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    Sampling_Machine.Initialize(eqs);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    Rectangular_Sample_Grids.Set_Polynomial_Type(laurent=>false);
    grid := Rectangular_Sample_Grids.Create1(sps,2);
    Monodromy_Component_Breakup.Factor(dim,nbl,grid,f);
    Standard_Complex_VecVecs.Clear(hyp);
    Sampling_Machine.Clear;
  exception
    when others =>
      put_line("An exception happened in Witness_Factor."); raise;
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    hyp : Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Witness_Sets.Slices(eqs,dim);
    sps : constant Standard_Sample_List := Create(pts,hyp);
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    Sampling_Laurent_Machine.Initialize(eqs);
    Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    Rectangular_Sample_Grids.Set_Polynomial_Type(laurent=>true);
    grid := Rectangular_Sample_Grids.Create1(sps,2);
    Monodromy_Component_Breakup.Laurent_Factor(dim,nbl,grid,f);
    Standard_Complex_VecVecs.Clear(hyp);
    Sampling_Laurent_Machine.Clear;
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    hyp : DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Witness_Sets.Slices(eqs,dim);
    sps : constant DoblDobl_Sample_List := Create(pts,hyp);
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    DoblDobl_Sampling_Machine.Initialize(eqs);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    DoblDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(laurent=>false);
    grid := DoblDobl_Rectangular_Sample_Grids.Create1(sps,2);
    Monodromy_Component_Breakup.Factor(dim,nbl,grid,f);
    DoblDobl_Complex_VecVecs.Clear(hyp);
    DoblDobl_Sampling_Machine.Clear;
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    hyp : DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Witness_Sets.Slices(eqs,dim);
    sps : constant DoblDobl_Sample_List := Create(pts,hyp);
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(eqs);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    DoblDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(laurent=>true);
    grid := DoblDobl_Rectangular_Sample_Grids.Create1(sps,2);
    Monodromy_Component_Breakup.Laurent_Factor(dim,nbl,grid,f);
    DoblDobl_Complex_VecVecs.Clear(hyp);
    DoblDobl_Sampling_Laurent_Machine.Clear;
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    hyp : QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Witness_Sets.Slices(eqs,dim);
    sps : constant QuadDobl_Sample_List := Create(pts,hyp);
    grid : Array_of_QuadDobl_Sample_Lists(0..2);
    copiedeqs : QuadDobl_Complex_Poly_Systems.Poly_Sys(eqs'range);

  begin
    QuadDobl_Sampling_Machine.Initialize(eqs);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    QuadDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(laurent=>false);
    grid := QuadDobl_Rectangular_Sample_Grids.Create1(sps,2);
    Monodromy_Component_Breakup.Factor(dim,nbl,grid,f);
    QuadDobl_Complex_VecVecs.Clear(hyp);
    QuadDobl_Sampling_Machine.Clear;
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    hyp : QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Witness_Sets.Slices(eqs,dim);
    sps : constant QuadDobl_Sample_List := Create(pts,hyp);
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(eqs);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(0);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    QuadDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(laurent=>true);
    grid := QuadDobl_Rectangular_Sample_Grids.Create1(sps,2);
    Monodromy_Component_Breakup.Laurent_Factor(dim,nbl,grid,f);
    QuadDobl_Complex_VecVecs.Clear(hyp);
    QuadDobl_Sampling_Laurent_Machine.Clear;
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration ) is

    level_timer,total_timer : Timing_Widget;

  begin
    tstart(total_timer);
    for k in reverse 1..topdim loop
      tstart(level_timer);
      if not Standard_Complex_Solutions.Is_Null(pts(integer32(k))) then
        Witness_Factor
          (verbose,eqs(integer32(k)).all,pts(integer32(k)),
           k,nbl,tol,f(integer32(k)));
      end if;
      tstop(level_timer);
      times(integer(k)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration ) is

    level_timer,total_timer : Timing_Widget;

  begin
    tstart(total_timer);
    for k in reverse 1..topdim loop
      tstart(level_timer);
      if not Standard_Complex_Solutions.Is_Null(pts(integer32(k))) then
        Witness_Factor
          (verbose,eqs(integer32(k)).all,pts(integer32(k)),
           k,nbl,tol,f(integer32(k)));
      end if;
      tstop(level_timer);
      times(integer(k)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration ) is

    level_timer,total_timer : Timing_Widget;

  begin
    tstart(total_timer);
    for k in reverse 1..topdim loop
      tstart(level_timer);
      if not DoblDobl_Complex_Solutions.Is_Null(pts(integer32(k))) then
        Witness_Factor
          (verbose,eqs(integer32(k)).all,pts(integer32(k)),
           k,nbl,tol,f(integer32(k)));
      end if;
      tstop(level_timer);
      times(integer(k)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration ) is

    level_timer,total_timer : Timing_Widget;

  begin
    tstart(total_timer);
    for k in reverse 1..topdim loop
      tstart(level_timer);
      if not DoblDobl_Complex_Solutions.Is_Null(pts(integer32(k))) then
        Witness_Factor
          (verbose,eqs(integer32(k)).all,pts(integer32(k)),
           k,nbl,tol,f(integer32(k)));
      end if;
      tstop(level_timer);
      times(integer(k)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration ) is

    level_timer,total_timer : Timing_Widget;

  begin
    tstart(total_timer);
    for k in reverse 1..topdim loop
      tstart(level_timer);
      if not QuadDobl_Complex_Solutions.Is_Null(pts(integer32(k))) then
        Witness_Factor
          (verbose,eqs(integer32(k)).all,pts(integer32(k)),
           k,nbl,tol,f(integer32(k)));
      end if;
      tstop(level_timer);
      times(integer(k)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Witness_Factor;

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Array_of_VecVecs;
                times : out Array_of_Duration; alltime : out duration ) is

    level_timer,total_timer : Timing_Widget;

  begin
    tstart(total_timer);
    for k in reverse 1..topdim loop
      tstart(level_timer);
      if not QuadDobl_Complex_Solutions.Is_Null(pts(integer32(k))) then
        Witness_Factor
          (verbose,eqs(integer32(k)).all,pts(integer32(k)),
           k,nbl,tol,f(integer32(k)));
      end if;
      tstop(level_timer);
      times(integer(k)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Witness_Factor;

end Monodromy_Homotopies;
