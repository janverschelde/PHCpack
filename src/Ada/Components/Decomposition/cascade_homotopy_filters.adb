with Cascade_Homotopies;                 use Cascade_Homotopies;
with Homotopy_Membership_Filters;        use Homotopy_Membership_Filters;
with Monodromy_Homotopies;               use Monodromy_Homotopies;

package body Cascade_Homotopy_Filters is

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    Witness_Generate(outfile,resfile,nt,ep,sols,topdim,lowdim,zerotol);
    tstop(timer);
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is
 
    timer : Timing_Widget;

  begin
    tstart(timer);
    Witness_Generate(outfile,resfile,nt,ep,sols,topdim,lowdim,zerotol);
    tstop(timer);
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    Witness_Generate(outfile,resfile,nt,ep,sols,topdim,lowdim,zerotol);
    tstop(timer);
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    Witness_Generate(outfile,resfile,nt,ep,sols,topdim,lowdim,zerotol);
    tstop(timer);
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    Witness_Generate(outfile,resfile,nt,ep,sols,topdim,lowdim,zerotol);
    tstop(timer);
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is

    timer : Timing_Widget;

  begin
    tstart(timer);
    Witness_Generate(outfile,resfile,nt,ep,sols,topdim,lowdim,zerotol);
    tstop(timer);
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim,lowdim : in natural32; zerotol : in double_float;
                 rcotol,restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration ) is

    timer : Timing_Widget;

  begin
    castms := (castms'range => 0.0);
    filtms := (filtms'range => 0.0);
    tstart(timer);
    Witness_Generate
      (nt,ep,sols,topdim,lowdim,zerotol,embsys,esols0,pathcnts,castms,totcas);
    Filter(false,embsys,esols0,integer32(topdim),rcotol,restol,homtol,
           filtcnts,filtms,totfil);
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration ) is

    timer : Timing_Widget;

  begin
    castms := (castms'range => 0.0);
    filtms := (filtms'range => 0.0);
    tstart(timer);
    Witness_Generate
      (nt,ep,sols,topdim,lowdim,zerotol,embsys,esols0,pathcnts,castms,totcas);
    Filter(false,embsys,esols0,integer32(topdim),rcotol,restol,homtol,
           filtcnts,filtms,totfil);
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration ) is

    timer : Timing_Widget;

  begin
    castms := (castms'range => 0.0);
    filtms := (filtms'range => 0.0);
    tstart(timer);
    Witness_Generate
      (nt,ep,sols,topdim,lowdim,zerotol,embsys,esols0,pathcnts,castms,totcas);
    Filter(false,embsys,esols0,integer32(topdim),rcotol,restol,homtol,
           filtcnts,filtms,totfil);
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration ) is

    timer : Timing_Widget;

  begin
    castms := (castms'range => 0.0);
    filtms := (filtms'range => 0.0);
    tstart(timer);
    Witness_Generate
      (nt,ep,sols,topdim,lowdim,zerotol,embsys,esols0,pathcnts,castms,totcas);
    Filter(false,embsys,esols0,integer32(topdim),rcotol,restol,homtol,
           filtcnts,filtms,totfil);
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration ) is

    timer : Timing_Widget;

  begin
    castms := (castms'range => 0.0);
    filtms := (filtms'range => 0.0);
    tstart(timer);
    Witness_Generate
      (nt,ep,sols,topdim,lowdim,zerotol,embsys,esols0,pathcnts,castms,totcas);
    Filter(false,embsys,esols0,integer32(topdim),rcotol,restol,homtol,
           filtcnts,filtms,totfil);
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim,lowdim : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms : out Array_of_Duration;
                totcas,totfil,alltime : out duration ) is

    timer : Timing_Widget;

  begin
    castms := (castms'range => 0.0);
    filtms := (filtms'range => 0.0);
    tstart(timer);
    Witness_Generate
      (nt,ep,sols,topdim,lowdim,zerotol,embsys,esols0,pathcnts,castms,totcas);
    Filter(false,embsys,esols0,integer32(topdim),rcotol,restol,homtol,
           filtcnts,filtms,totfil);
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Filter;

  procedure Witness_Factor
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim,lowdim,nitfix : in natural32; zerotol : in double_float;
                rcotol,restol,homtol : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                factors : out Standard_Natural_VecVecs.Array_of_VecVecs;
                pathcnts,filtcnts : out Standard_Natural_VecVecs.VecVec;
                castms,filtms,factms : out Array_of_Duration;
                totcas,totfil,totfac,alltime : out duration ) is

    timer : Timing_Widget;

  begin
    castms := (castms'range => 0.0);
    filtms := (filtms'range => 0.0);
    factms := (factms'range => 0.0);
    tstart(timer);
    Witness_Generate
      (nt,ep,sols,topdim,lowdim,zerotol,embsys,esols0,pathcnts,castms,totcas);
    Filter(false,embsys,esols0,integer32(topdim),rcotol,restol,homtol,
           filtcnts,filtms,totfil);
    Witness_Factor
      (false,embsys,esols0,topdim,nitfix,zerotol,factors,factms,totfac);
    tstop(timer);
    alltime := Elapsed_User_Time(timer);
  end Witness_Factor;

end Cascade_Homotopy_Filters;
