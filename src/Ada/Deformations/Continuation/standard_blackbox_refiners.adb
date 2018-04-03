with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Solution_Filters;
with Standard_Solution_Splitters;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with Root_Refining_Parameters;           use Root_Refining_Parameters;
with Multitasking_Root_Refiners;         use Multitasking_Root_Refiners;

with Standard_Complex_Solutions_io;
 use Standard_Complex_Solutions_io;

package body Standard_BlackBox_Refiners is

-- WITHOUT MULTITASKING :

  procedure Silent_Black_Box_Refine
              ( p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    ref_sols : Solution_List;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;

  begin
    if Length_Of(sols) > 0 then
      Standard_Default_Root_Refining_Parameters
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
      Standard_Default_Root_Refining_Parameters
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
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,deflate,wout);
    end if;
    Clear(sols);
    sols := ref_sols;
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
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
      Reporting_Root_Refiner
        (file,p,sols,ref_sols,epsxa,epsfa,tolsing,nb,maxit,wout);
    end if;
    Clear(sols);
    sols := ref_sols;
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
                p : in Laur_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    ref_sols : Solution_List;

  begin
    if Length_Of(sols) > 0 then
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
                p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    vansols,regsols,sinsols,ref_sinsols : Solution_List;

  begin
    if Length_Of(sols) > 0 then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
     -- refine only the vanishing solutions
      vansols := Standard_Solution_Filters.Vanishing_Filter(sols,epsfa);
     -- Copy(sols,vansols);
     -- apply deflation to singular solutions, split the list first
      Standard_Solution_Splitters.Silent_Singular_Filter
        (vansols,tolsing,sinsols,regsols);
     -- run Newton's method only on the regular solutions
      Mute_Multitasking_Root_Refiner
        (nt,p,regsols,epsxa,epsfa,tolsing,nb,maxit,deflate);
     -- apply deflation only on the singular solutions
      nb := 0;
      Silent_Root_Refiner
        (p,sinsols,ref_sinsols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      Push(ref_sinsols,regsols);
      Clear(sols); Clear(vansols); Clear(sinsols);
      sols := regsols;
    end if;
  end Silent_Black_Box_Refine;

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List ) is

    epsxa,epsfa,tolsing : double_float;
    maxit,nb : natural32 := 0;
    deflate,wout : boolean;
    vansols,regsols,sinsols,ref_sinsols : Solution_List;

  begin
    if Length_Of(sols) > 0 then
      Standard_Default_Root_Refining_Parameters
        (epsxa,epsfa,tolsing,maxit,deflate,wout);
     -- refine only the vanishing solutions
      vansols := Standard_Solution_Filters.Vanishing_Filter(sols,epsfa);
     -- Copy(sols,vansols);
     -- apply deflation to singular solutions, split the list first
      Standard_Solution_Splitters.Silent_Singular_Filter
        (vansols,tolsing,sinsols,regsols);
     -- refine the regular solutions just with Newton's method
      Silent_Multitasking_Root_Refiner -- tasks remain silent
        (file,nt,p,regsols,epsxa,epsfa,tolsing,nb,maxit,deflate);
     -- running Newton's method prior to deflation may not be good
      nb := 0;
      Reporting_Root_Refiner
        (file,p,sinsols,ref_sinsols,epsxa,epsfa,tolsing,
         nb,maxit,deflate,wout);
     -- Silent_Root_Refiner
     --   (p,sinsols,ref_sinsols,epsxa,epsfa,tolsing,nb,maxit,deflate);
      new_line(file);
      put_line(file,"THE REFINED SINGULAR SOLUTIONS :");
      put(file,Length_Of(ref_sinsols),natural32(Head_Of(ref_sinsols).n),
          ref_sinsols);
      Push(ref_sinsols,regsols);
      Clear(sols); Clear(vansols); Clear(sinsols);
      sols := regsols;
    end if;
  end Reporting_Black_Box_Refine;

end Standard_BlackBox_Refiners;
