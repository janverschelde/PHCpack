with text_io;                           use text_io;

package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.19 released 2016-06-03";

  begin
    return res;
  end Version;

  procedure show_help is
  begin
    put_line("PHC stands for Polynomial Homotopy Continuation,");
    put_line("to numerically solve systems of polynomial equations.");
    new_line;
    put_line("The quickest use of phc is via the option -b,");
    put_line("on the command line as phc -b input output,");
    put_line("where input and output are file names for input and output.");
    put_line("Note that the quickest use may not give the fastest solver.");
    new_line;
    put_line("The input file starts with the number of polynomials,");
    put_line("followed by the polynomials in symbolic form.");
    put_line("Each polynomial ends with a semi-colon, for example:");
    new_line;
    put_line("2 ");
    put_line(" x*y + 2*y^3 - 3.14E-01;");
    put_line(" 3*x^2 + y^2 - 2/3;");
    new_line;
    put_line("In double precision, 2/3 becomes 6.66666666666667E-01.");
    put_line("For double double precision, use phc -b2,");
    put_line("and for quad double precision, use phc -b4.");
    put_line("To compensate for the higher cost of extended precision,");
    put_line("use multithreading, e.g.: phc -t4, to use 4 threads.");
  end show_help;

  procedure help4setseed is
  begin
    put_line("phc -0 fixes the seed in the random number generators.");
  end help4setseed;

  procedure help4eqnbyeqn is
  begin
    put_line("phc -a calls an equation-by-equation solver.");
  end help4eqnbyeqn;

  procedure help4blackbox is
  begin
    put_line("phc -b calls the blackbox solver.");
  end help4blackbox;

  procedure help4components is
  begin
    put_line("phc -c provides a numerical irreducible decomposition.");
  end help4components;

  procedure help4reduction is
  begin
    put_line("phc -d attempts to lower the degrees of the polynomials.");
  end help4reduction;

  procedure help4enumeration is
  begin
    put_line("phc -e gives homotopies for numerical Schubert calculus.");
  end help4enumeration;

  procedure help4feedback is
  begin
    put_line("phc -k computes feedback laws to control linear systems.");
  end help4feedback;

  procedure help4factor is
  begin
    put_line("phc -f factors a solution set into irreducible components.");   
  end help4factor;

  procedure help4goodformat is
  begin
    put_line("phc -g checks whether the format of an input system is good.");
  end help4goodformat;

  procedure help4help is
  begin
    put_line("phc -h or phc --help tries to help...");
  end help4help;

  procedure help4hypersurface is
  begin
    put_line("phc -l to compute a witness set for a hypersurface.");
  end help4hypersurface;

  procedure help4mixvol is
  begin
    put_line("phc -m computes mixed volumes and runs polyhedral homotopies.");
  end help4mixvol;

  procedure help4symbols is
  begin
    put_line("phc -o writes the symbol table after parsing an input system.");
  end help4symbols;

  procedure help4continuation is
  begin
    put_line("phc -p runs continuation with a homotopy in one parameter.");
  end help4continuation;

  procedure help4jumpstart is
  begin
    put_line("phc -q runs continuation with jump starting for huge systems.");
  end help4jumpstart;

  procedure help4rootcounts is
  begin
    put_line("phc -r computes root counts and constructs start systems.");
  end help4rootcounts;

  procedure help4scaling is
  begin
    put_line("phc -s applies equation and variable scaling to a system.");
  end help4scaling;

  procedure help4tasking is
  begin
    put_line("phc -t# computes with # tasks, for shared memory parallelism.");
  end help4tasking;

  procedure help4verification is
  begin
    put_line("phc -v filters, verifies, and refines lists of solutions.");
  end help4verification;

  procedure help4witsetinsect is
  begin
    put_line("phc -w for witness set intersection with diagonal homotopies.");
  end help4witsetinsect;

  procedure help4pythondict is
  begin
    put_line("phc -x converts lists of solutions into Python dictionaries.");
  end help4pythondict;

  procedure help4sampler is
  begin
    put_line("phc -y samples points on a positive dimensional solution set.");
  end help4sampler;

  procedure help4mapleform is
  begin
    put_line("phc -z converts lists of solutions into Maple format.");
  end help4mapleform;

end Greeting_Banners;
