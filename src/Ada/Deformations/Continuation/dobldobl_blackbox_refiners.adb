with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Solution_Manipulators;
with DoblDobl_Solution_Filters;
with DoblDobl_Solution_Splitters;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with Root_Refining_Parameters;           use Root_Refining_Parameters;
with Multitasking_Root_Refiners;         use Multitasking_Root_Refiners;

package body DoblDobl_BlackBox_Refiners is

-- WITHOUT MULTITASKING :

  procedure Silent_Black_Box_Refine
              ( p : in Poly_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if verbose > 0 then
      put("-> in dobldobl_blackbox_refiners.");
      put_line("Silent_Black_Box_Refine 1 ...");
    end if;
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Root_Refiner
        (p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,deflate,verbose-1);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Silent_Black_Box_Refine
              ( p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if verbose > 0 then
      put("-> in dobldobl_blackbox_refiners.");
      put_line("Silent_Black_Box_Refine 2 ...");
    end if;
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Root_Refiner
        (p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,verbose-1);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Poly_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if verbose > 0 then
      put("-> in dobldobl_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 1 ...");
    end if;
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,
         nb,maxit,deflate,wout,verbose-1);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Reporting_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if verbose > 0 then
      put("-> in dobldobl_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 2 ...");
    end if;
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,false,verbose-1);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Reporting_Black_Box_Refine;

-- WITH MULTITASKING ON LAURENT SYSTEMS :

 -- procedure Silent_Black_Box_Refine
 --             ( nt : in integer32;
 --               p : in Laur_Sys; sols : in out Solution_List ) is
 --
 --   epsxa,epsfa,tolsing : double_float;
 --   maxit,nb : natural32 := 0;
 --   deflate,wout : boolean;
 --   ref_sols : Solution_List;
 --
 -- begin
 --   if Length_Of(sols) > 0 then
 --     DoblDobl_Default_Root_Refining_Parameters
 --       (epsxa,epsfa,tolsing,maxit,deflate,wout);
 --     Mute_Multitasking_Root_Refiner
 --       (nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
 --     ref_sols := DoblDobl_Solution_Filters.Vanishing_Filter(sols,epsfa);
 --     Clear(sols);
 --     sols := ref_sols;
 --   end if;
 -- end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    ref_sols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in dobldobl_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 3 ...");
    end if;
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
                p : in Poly_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    tarsols,vansols,regsols,sinsols,ref_sinsols : Solution_List;
    ddtarget : constant double_double := create(1.0);
    target : constant Complex_Number := Create(ddtarget);

  begin
    if verbose > 0 then
      put("-> in dobldobl_blackbox_refiners.");
      put_line("Silent_Black_Box_Refine 3 ...");
    end if;
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Mute_Multitasking_Root_Refiner
        (nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
      tarsols := DoblDobl_Solution_Filters.On_Target_Filter(sols,target,epsfa);
      vansols := DoblDobl_Solution_Filters.Vanishing_Filter(tarsols,epsfa);
      DoblDobl_Solution_Splitters.Silent_Singular_Filter
        (vansols,tolsing,sinsols,regsols);
      nb := 0;
      if not Is_Null(sinsols) then
        Silent_Root_Refiner
          (p,sinsols,ref_sinsols,epsxa,epsfa,tolsing,
           nb,maxit,deflate,verbose-1);
        Push(ref_sinsols,regsols);
      end if;
      Clear(sols); Clear(vansols); Clear(sinsols); Clear(tarsols);
      sols := regsols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    tarsols,vansols,regsols,sinsols,ref_sinsols : Solution_List;
    ddtarget : constant double_double := create(1.0);
    target : constant Complex_Number := Create(ddtarget);

  begin
    if verbose > 0 then
      put("-> in dobldobl_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 4 ...");
    end if;
    if Length_Of(sols) > 0 then
      DoblDobl_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Multitasking_Root_Refiner -- tasks remain silent
        (file,nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      DoblDobl_Solution_Manipulators.Remove_Imaginary_Target(sols);
      tarsols := DoblDobl_Solution_Filters.On_Target_Filter(sols,target,epsfa);
      vansols := DoblDobl_Solution_Filters.Vanishing_Filter(tarsols,epsfa);
      DoblDobl_Solution_Splitters.Silent_Singular_Filter
        (vansols,tolsing,sinsols,regsols);
      if not Is_Null(sinsols) then
        Reporting_Root_Refiner
          (file,p,sinsols,ref_sinsols,epsxa,epsfa,tolsing,
           nb,maxit,deflate,verbose-1);
        Push(ref_sinsols,regsols);
      end if;
      Clear(sols); Clear(vansols); Clear(sinsols); Clear(tarsols);
      sols := regsols;
    end if;
  end Reporting_Black_Box_Refine;

end DoblDobl_BlackBox_Refiners;
