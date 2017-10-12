with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_VecVecs;
with Witness_Sets;
with Sampling_Machine;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Rectangular_Sample_Grids;
with Monodromy_Component_Breakup;

package body Monodromy_Homotopies is

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    hyp : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    sps : constant Standard_Sample_List := Create(pts,hyp);
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    Sampling_Machine.Initialize(eqs);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    grid := Rectangular_Sample_Grids.Create1(sps,2);
    Monodromy_Component_Breakup.Factor(eqs,dim,grid,f);
    Sampling_Machine.Clear;
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
      Witness_Factor
        (verbose,eqs(integer32(k)).all,pts(integer32(k)),
         k,nbl,tol,f(integer32(k)));
      tstop(level_timer);
      times(integer(k)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Witness_Factor;

end Monodromy_Homotopies;
