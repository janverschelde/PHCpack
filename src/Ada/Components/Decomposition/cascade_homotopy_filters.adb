with Cascade_Homotopies;                 use Cascade_Homotopies;

package body Cascade_Homotopy_Filters is

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( outfile,resfile : in file_type; nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in DoblDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
               ( name : in string; outfile : in file_type;
                 nt : in natural32;
                 ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in QuadDobl_Complex_Solutions.Solution_List;
                 topdim : in natural32; zerotol : in double_float;
                 restol,homtol : in double_float ) is
  begin
    null;
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float;
                restol,homtol : in double_float;
                embsys : out Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is
  begin
    Witness_Generate
      (nt,ep,sols,topdim,zerotol,embsys,esols0,pathcnts,times,alltime);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float;
                restol,homtol : in double_float;
                embsys : out Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out Standard_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is
  begin
    Witness_Generate
      (nt,ep,sols,topdim,zerotol,embsys,esols0,pathcnts,times,alltime);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float;
                restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is
  begin
    Witness_Generate
      (nt,ep,sols,topdim,zerotol,embsys,esols0,pathcnts,times,alltime);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float;
                restol,homtol : in double_float;
                embsys : out DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is
  begin
    Witness_Generate
      (nt,ep,sols,topdim,zerotol,embsys,esols0,pathcnts,times,alltime);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float;
                restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is
  begin
    Witness_Generate
      (nt,ep,sols,topdim,zerotol,embsys,esols0,pathcnts,times,alltime);
  end Witness_Filter;

  procedure Witness_Filter
              ( nt : in natural32;
                ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                topdim : in natural32; zerotol : in double_float;
                restol,homtol : in double_float;
                embsys : out QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                esols0 : out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                pathcnts : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is
  begin
    Witness_Generate
      (nt,ep,sols,topdim,zerotol,embsys,esols0,pathcnts,times,alltime);
  end Witness_Filter;

end Cascade_Homotopy_Filters;
