with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with DoblDobl_Solution_Filters;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with Root_Refining_Parameters;           use Root_Refining_Parameters;
with Multitasking_Root_Refiners;         use Multitasking_Root_Refiners;

package body DoblDobl_BlackBox_Refiners is

-- WITHOUT MULTITASKING :

  procedure Silent_Black_Box_Refine
              ( p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Root_Refiner
        (p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Silent_Black_Box_Refine
              ( p : in Laur_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Root_Refiner(p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,deflate,wout);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Reporting_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Laur_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,false);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Reporting_Black_Box_Refine;

-- WITH MULTITASKING ON LAURENT SYSTEMS :

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    ref_sols : Solution_List;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Mute_Multitasking_Root_Refiner
        (nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      ref_sols := DoblDobl_Solution_Filters.Vanishing_Filter(sols,epsfa);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    ref_sols : Solution_List;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Multitasking_Root_Refiner -- tasks remain silent
        (file,nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      ref_sols := DoblDobl_Solution_Filters.Vanishing_Filter(sols,epsfa);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Reporting_Black_Box_Refine;

-- WITH MULTITASKING ON POLYNOMIAL SYSTEMS :

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    ref_sols : Solution_List;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Mute_Multitasking_Root_Refiner
        (nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      ref_sols := DoblDobl_Solution_Filters.Vanishing_Filter(sols,epsfa);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    ref_sols : Solution_List;

  begin
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Multitasking_Root_Refiner -- tasks remain silent
        (file,nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      ref_sols := DoblDobl_Solution_Filters.Vanishing_Filter(sols,epsfa);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Reporting_Black_Box_Refine;

end DoblDobl_BlackBox_Refiners;
