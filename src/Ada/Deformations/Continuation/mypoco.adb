with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Drivers_for_Homotopy_Creation;      use Drivers_for_Homotopy_Creation;
with Increment_and_Fix_Continuation;     use Increment_and_Fix_Continuation;
with Drivers_for_Poly_Continuation;      use Drivers_for_Poly_Continuation;
with MyHomotopy;

procedure mypoco is

  procedure Continue ( file : in file_type; sols : in out Solution_List;
                       report : in boolean;
                       target : in Complex_Number ) is

  -- DESCRIPTION :
  --   Instantiates the path-trackers.

    timer : Timing_Widget;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,
                          MyHomotopy.Eval,MyHomotopy.Diff,MyHomotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,
                             MyHomotopy.Eval,MyHomotopy.Diff,MyHomotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,false,target);
     else Sil_Cont(sols,false,target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end Continue;

  procedure Main is

    infile,outfile : file_type;
    found : boolean;
    sols : Solution_List;
    k : positive;
    a,target : Complex_Number;
    oc : natural;
    report : boolean;

  begin
    put_line
      ("Polynomial Continuation with inline evaluators and differentiators.");
    new_line;
    put_line("Reading the name of the file where the start solutions are.");
    Read_Name_and_Open_File(infile);
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found
     then get(infile,sols);
     else Reset(infile);
          get(infile,sols);
    end if;
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    put_line(outfile,"THE START SOLUTIONS : ");
    put(outfile,Length_Of(sols),Head_Of(sols).n,sols); new_line(outfile);
    Default_Homotopy_Settings(k,a,target);
    Menu_for_Homotopy_Settings(outfile,k,a,target);
    MyHomotopy.Init(a,k);
    new_line;
    Driver_for_Continuation_Parameters(outfile);
    new_line;
    Driver_for_Process_io(outfile,oc);
    report := not (oc = 0);
    Continue(outfile,sols,report,target);
  end Main;

begin
  Main;
end mypoco;
