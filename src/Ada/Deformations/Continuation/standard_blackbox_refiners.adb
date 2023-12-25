with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Solution_Manipulators;
with Standard_Solution_Filters;
with Standard_Solution_Splitters;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Root_Refining_Parameters;           use Root_Refining_Parameters;
with Multitasking_Root_Refiners;         use Multitasking_Root_Refiners;

--with standard_natural_numbers_io;
-- use standard_natural_numbers_io;
--with Standard_Complex_Solutions_io;
-- use Standard_Complex_Solutions_io;

package body Standard_BlackBox_Refiners is

-- WITHOUT MULTITASKING :

  procedure Silent_Black_Box_Refine
              ( p : in Poly_Sys; sols : in out Solution_List;
                deflate : in boolean; verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    default_deflate,wout : boolean;
    vardeflate : boolean := deflate;

  begin
    if verbose > 0 then
      put("-> in standard_blackbox_refiners.");
      put_line("Silent_Black_Box_Refine 1 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,default_deflate,wout);
      Silent_Root_Refiner
        (p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,vardeflate,verbose-1);
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
      put("-> in standard_blackbox_refiners.");
      put_line("Silent_Black_Box_Refine 2 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
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
                deflate : in boolean; verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    default_deflate,wout : boolean;
    vardeflate : boolean := deflate;

  begin
    if verbose > 0 then
      put("-> in standard_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 1 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,default_deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,
         vardeflate,wout,verbose-1);
    end if;
    Clear(sols);
    sols := ref_sols;
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
      put("-> in standard_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 2 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,wout,verbose-1);
    end if;
    Clear(sols);
    sols := ref_sols;
  end Reporting_Black_Box_Refine;

-- WITH MULTITASKING ON LAURENT SYSTEMS :

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    ref_sols : Solution_List;

  begin
    if verbose > 0 then
      put("-> in standard_blackbox_refiners.");
      put_line("Silent_Black_Box_Refine 3 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Mute_Multitasking_Root_Refiner
        (nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      ref_sols := Standard_Solution_Filters.Vanishing_Filter(sols,epsfa);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Silent_Black_Box_Refine;

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
      put("-> in standard_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 3 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Silent_Multitasking_Root_Refiner -- tasks remain silent
        (file,nt,p,sols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      ref_sols := Standard_Solution_Filters.Vanishing_Filter(sols,epsfa);
      Clear(sols);
      sols := ref_sols;
    end if;
  end Reporting_Black_Box_Refine;

-- WITH MULTITASKING ON POLYNOMIAL SYSTEM :

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List;
                deflate : in boolean; verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    default_deflate,wout : boolean;
    vardeflate : boolean := deflate;
    tarsols,vansols,regsols,sinsols,ref_sinsols : Solution_List;
    target : constant Complex_Number := Create(1.0);

  begin
    if verbose > 0 then
      put("-> in standard_blackbox_refiners.");
      put_line("Silent_Black_Box_Refine 4 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,default_deflate,wout);
     -- refine only the vanishing solutions that reached the target
     -- put("The number of solutions : "); put(Length_Of(sols),1); new_line;
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
      tarsols := Standard_Solution_Filters.On_Target_Filter(sols,target,epsfa);
     -- put("The number of solutions on target : ");
     -- put(Length_Of(tarsols),1); new_line;
      vansols := Standard_Solution_Filters.Vanishing_Filter(tarsols,epsfa);
     -- put("The number of vanishing solutions on target : ");
     -- put(Length_Of(vansols),1); new_line;
      if not Is_Null(vansols) then
        if not deflate then
          Mute_Multitasking_Root_Refiner
            (nt,p,vansols,epsxa,epsfa,tolsing,nb,maxit,vardeflate);
          Clear(sols); sols := vansols;
        else
         -- apply deflation to singular solutions, split the list first
          Standard_Solution_Splitters.Silent_Singular_Filter
            (vansols,tolsing,sinsols,regsols);
         -- run Newton's method only on the regular solutions
          if not Is_Null(regsols) then
            Mute_Multitasking_Root_Refiner
              (nt,p,regsols,epsxa,epsfa,tolsing,nb,maxit,vardeflate);
          end if;
         -- apply deflation only on the singular solutions
          if not Is_Null(sinsols) then
            nb := 0;
            Silent_Root_Refiner
              (p,sinsols,ref_sinsols,epsxa,epsfa,tolsing,nb,maxit,vardeflate);
            Push(ref_sinsols,regsols);
          end if;
          Clear(sols); Clear(vansols); Clear(sinsols);
          sols := regsols;
        end if;
      end if;
      Clear(tarsols);
    end if;
  end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List;
                deflate : in boolean; verbose : in integer32 := 0 ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    default_deflate,wout : boolean;
    vardeflate : boolean := deflate;
    tarsols,vansols,regsols,sinsols,ref_sinsols : Solution_List;
    target : constant Complex_Number := Create(1.0);

  begin
    if verbose > 0 then
      put("-> in standard_blackbox_refiners.");
      put_line("Reporting_Black_Box_Refine 4 ...");
    end if;
    if not Is_Null(sols) then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,default_deflate,wout);
     -- refine only the vanishing solutions that reached the target
      Standard_Solution_Manipulators.Remove_Imaginary_Target(sols);
      tarsols := Standard_Solution_Filters.On_Target_Filter(sols,target,epsfa);
      vansols := Standard_Solution_Filters.Vanishing_Filter(tarsols,epsfa);
      if not Is_Null(vansols) then
        if not deflate then
          Silent_Multitasking_Root_Refiner -- tasks remain silent
            (file,nt,p,vansols,epsxa,epsfa,tolsing,nb,maxit,vardeflate);
          Clear(sols); sols := vansols;
        else
         -- apply deflation to singular solutions, split the list first
          Standard_Solution_Splitters.Silent_Singular_Filter
            (vansols,tolsing,sinsols,regsols);
         -- refine the regular solutions just with Newton's method
          if not Is_Null(regsols) then
            Silent_Multitasking_Root_Refiner -- tasks remain silent
              (file,nt,p,regsols,epsxa,epsfa,tolsing,nb,maxit,vardeflate);
          end if;
         -- running Newton's method prior to deflation may not be good
          if not Is_Null(sinsols) then
            nb := 0;
            Reporting_Root_Refiner
              (file,p,sinsols,ref_sinsols,epsxa,epsfa,tolsing,
               nb,maxit,vardeflate,wout);
            Push(ref_sinsols,regsols);
          end if;
          Clear(sols); Clear(vansols); Clear(sinsols);
          sols := regsols;
        end if;
      end if;
      Clear(tarsols);
    end if;
  end Reporting_Black_Box_Refine;

end Standard_BlackBox_Refiners;
