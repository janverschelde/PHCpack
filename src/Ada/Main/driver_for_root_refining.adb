with Timing_Package;                     use Timing_Package;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Scaling;                   use Standard_Scaling;
with Projective_Transformations;         use Projective_Transformations;
with Standard_Root_Refiners;             use Standard_Root_Refiners;

procedure Driver_for_Root_Refining
             ( file : in file_type;
               scalp,p : in Poly_Sys; basis : in natural32;
               scalvec : in Link_to_Vector; sols : in out Solution_List;
               verbose : in integer32 := 0 ) is

  numb : natural32;
  epsxa,epsfa : constant double_float := 1.0E-8;
  tolsing : constant double_float := 1.0E-8;
  timer : timing_widget;
  len : constant natural32 := Length_Of(sols);
  deflate : boolean := false;

begin
  if verbose > 0
   then put_line("-> in driver_for_root_refining ...");
  end if;
  if (len /= 0) and then Head_Of(sols).n > p'last
   then Affine_Transformation(sols);
  end if;
  if scalvec /= null then
    put_line(file,"ROOT REFINING ON THE SCALED SYSTEM :");
    tstart(timer);
    numb := 0;
    Reporting_Root_Refiner
      (file,scalp,sols,epsxa,epsfa,tolsing,numb,5,deflate,false);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Root Refining on the Scaled System");
    Scale(basis,scalvec.all,sols);
  end if;
  tstart(timer);
  numb := 0;
  Reporting_Root_Refiner(file,p,sols,epsxa,epsfa,tolsing,numb,5,deflate,false);
  tstop(timer);
  new_line(file);
  print_times(file,timer,"Root Refining");
end Driver_for_Root_Refining;
