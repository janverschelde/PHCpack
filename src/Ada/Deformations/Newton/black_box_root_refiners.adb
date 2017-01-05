with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Root_Refiners;             use Standard_Root_Refiners;
with DoblDobl_Root_Refiners;             use DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;             use QuadDobl_Root_Refiners;
with Root_Refining_Parameters;           use Root_Refining_Parameters;

package body Black_Box_Root_Refiners is

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;
    dim : constant integer32 := Head_Of(sols).n;

  begin
    Standard_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    if p'last = dim then
      tstart(timer);
      Reporting_Root_Refiner
        (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
      tstop(timer);
    else
      tstart(timer);
      Reporting_Root_Sharpener
        (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
      tstop(timer);
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Refine_Roots
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;

  begin
    DoblDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    tstart(timer);
    Reporting_Root_Refiner
      (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Refine_Roots
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    epsxa,epsfa,tolsing : double_float;
    maxit,numb : natural32 := 0;
    deflate,wout : boolean;
    refsols : Solution_List;
    timer : Timing_Widget;

  begin
    QuadDobl_Default_Root_Refining_Parameters
      (epsxa,epsfa,tolsing,maxit,deflate,wout);
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    Standard_Put_Root_Refining_Parameters
      (file,epsxa,epsfa,tolsing,maxit,deflate,wout);
    tstart(timer);
    Reporting_Root_Refiner
      (file,p,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
    tstop(timer);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

end Black_Box_Root_Refiners;
