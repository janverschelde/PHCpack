with text_io;                            use text_io;
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

procedure bablvali ( infilename,outfilename : in string ) is

  lp : Link_to_Poly_Sys;
  sysonfile : boolean;

  procedure Read_System ( file : in out file_type; filename : in string ) is
  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      sysonfile := true;
    else
      sysonfile := false;
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   sysonfile := false;
                   lp := null; return;
  end Read_System;

  procedure Refine_Roots
               ( file : in file_type; sols : in out Solution_List ) is

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
    if lp'last = dim then
      tstart(timer);
      Reporting_Root_Refiner
        (file,lp.all,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
      tstop(timer);
    else
      tstart(timer);
      Reporting_Root_Sharpener
        (file,lp.all,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,deflate,false);
      tstop(timer);
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Length_Of(refsols),natural32(Head_Of(refsols).n),refsols);
    new_line(file);
    print_times(file,timer,"Root refining");
  end Refine_Roots;

  procedure Main is

    infile,outfile : file_type;
    nbvar : natural32;
    ans : character;
    sols : Solution_List;
    n : integer32 := 0;
    found : boolean;

  begin
    Read_System(infile,infilename);
    if lp = null then
      new_line;
      put("Is the system on file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put_line("Reading the name of the input file.");
        Read_Name_and_Open_File(infile);
        get(infile,lp);
        sysonfile := true;
      else
        put("Give the dimension : "); get(n);
        lp := new Poly_Sys(1..n);
        put("Give "); put(n,1); put(" "); put(n,1); 
        put_line("-variate polynomials :");
        get(natural32(n),lp.all);
        skip_line;
        sysonfile := false;
      end if;
    end if;
    nbvar := Standard_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    Create_Output_File(outfile,outfilename);
    if lp'last = integer32(nbvar)
     then put(outfile,natural32(lp'last),lp.all);
     else put(outfile,natural32(lp'last),nbvar,lp.all);
    end if;
    if sysonfile then
      Scan_and_Skip(infile,"THE SOLUTIONS",found);
      if found
       then get(infile,sols);
      end if;
      Close(infile);
    else
      found := false;
    end if;
    if not found
     then new_line; Read(sols);
    end if;
    Refine_Roots(outfile,sols);
  end Main;

begin
  Main;
end bablvali;
