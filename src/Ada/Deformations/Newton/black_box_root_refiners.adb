with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Root_Refiners;             use Standard_Root_Refiners;

package body Black_Box_Root_Refiners is

-- DESCRIPTION :
--   Wraps the root refiners with basic settins for the tolerances.

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    epsxa,epsfa : constant double_float := 1.0E-8;
    tolsing : constant double_float := 1.0E-8;
    maxit : constant natural32 := 3;
    numb : natural32 := 0;
    deflate : boolean := false;
    refsols : Solution_List;
    timer : Timing_Widget;
    dim : constant integer32 := Head_Of(sols).n;

  begin
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS");
    put(file,"  tolerance for error on the root : ");
    put(file,epsxa,2,3,3); new_line(file);
    put(file,"  tolerance for residual          : ");
    put(file,epsfa,2,3,3); new_line(file);
    put(file,"  tolerance for singular roots    : ");
    put(file,tolsing,2,3,3); new_line(file);
    put(file,"  maximum number of iterations    : ");
    put(file,maxit,2); new_line(file);
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

end Black_Box_Root_Refiners;
