with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package body Monodromy_Homotopies is

  procedure Witness_Factor
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim,nbl : in natural32; tol : in double_float;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is
  begin
    null;
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
