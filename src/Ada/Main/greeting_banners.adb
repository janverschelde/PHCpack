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
    put_line("2");
    put_line(" x*y + 2*y^3 - 3.14E-01;");
    put_line(" 3*(x + y)*x - 2/3;");
    new_line;
    put_line("In double precision, 2/3 becomes 6.66666666666667E-01.");
    put_line("The letters e and E can not be used as names of variables,");
    put_line("as they are used in the scientific notation of numbers.");
    put_line("The letters i and I denote the imaginary unit, sqrt(-1).");
    new_line;
    put_line("For double double precision, use phc -b2,");
    put_line("and for quad double precision, use phc -b4.");
    put_line("To compensate for the higher cost of extended precision,");
    put_line("use multithreading, e.g.: phc -t4, to use 4 threads.");
  end show_help;

  procedure help4setseed is
  begin
    put_line("phc -0 fixes the seed in the random number generators.");
    new_line;
    put_line("Random numbers are used in the construction of a start system.");
    put_line("Even when the user provides a start system in phc -p or -q,");
    put_line("a random 'gamma' constant is generated to avoid singularities");
    put_line("along the solution paths with probability one.");
    put_line("Because of this, different runs on the same inputs may lead to");
    put_line("(1) fluctuations in execution time; and");
    put_line("(2) a different order of the solutions in the list on output.");
    new_line;
    put_line("By default, without the option -0, phc generates the seed");
    put_line("based on the moment of time, at the start of the program.");
    put_line("The generated seed is printed at the end of the output.");
    new_line;
    put_line("For example, to use 17891 as seed, do phc -017891");
    put_line("and this will result in a reproducible execution.");
  end help4setseed;

  procedure help4eqnbyeqn is
  begin
    put_line("phc -a calls an equation-by-equation solver.");
  end help4eqnbyeqn;

  procedure help4blackbox is
  begin
    put_line("phc -b calls the blackbox solvers.");
    new_line;
    put_line("For polynomial systems with as many equations as unknowns,");
    put_line("the blackbox solver operates in four stages:");
    put_line("(1) preprocessing with some form of scaling;");
    put_line("(2) root counting and the construction of a start system");
    put_line("(3) path tracking of the solutions from start to target");
    put_line("(4) postprocessing: filtering and refinement of roots.");
    new_line;
    put_line("The four stages can be called separately with the options:");
    put_line("(1) phc -s for scaling, phc -d for degree reduction;");
    put_line("(2) phc -r for root counting, phc -m for mixed volumes;");
    put_line("(3) phc -p or -q for path tracking;");
    put_line("(4) phc -v for verification and root refinement.");
    new_line;
    put_line("The blackbox solver recognizes several special cases:");
    put_line("(1) one polynomial in one variable;");
    put_line("(2) one system of linear equations;");
    put_line("(3) a system with exactly two monomials in every equation.");
    put_line("For these special, no polynomial continuation is needed.");
    new_line;
    put_line("By default, phc -b solves in hardware double precision.");
    put_line("With phc -b2, double double arithmetic is applied,");
    put_line("and with -b4, path tracking happens in quad double arithmetic.");
    put_line("To offset the extra cost of higher precision arithmetic,");
    put_line("multitasking can be applied, e.g.: as in phc -b4 -t8");
    put_line("to solve a system in quad double precision with 8 tasks.");
    new_line;
    put_line("The focus of phc -b is on isolated solutions and on systems");
    put_line("with as many equations as unknowns.  For more general systems"); 
    put_line("and to compute positive dimensional solution sets, see the");
    put_line("options -a and -c, respectively for an equation-by-equation");
    put_line("solver and for a numerical irreducible decomposition.");
  end help4blackbox;

  procedure help4components is
  begin
    put_line("phc -c provides a numerical irreducible decomposition.");
    new_line;
    put_line("In a numerical irreducible decomposition, positive dimensional");
    put_line("solution sets are represented by a set of generic points that");
    put_line("satisfy the given system and as many linear equations with");
    put_line("random coefficients as the dimension of the solution set.");
    put_line("The number of generic points in that so-called witness set");
    put_line("then equals the degree of the solution set.");
  end help4components;

  procedure help4reduction is
  begin
    put_line("phc -d attempts to lower the degrees of the polynomials.");
  end help4reduction;

  procedure help4enumeration is
  begin
    put_line("phc -e gives homotopies for numerical Schubert calculus.");
    new_line;
    put_line("A Pieri problem asks to compute all p-planes that meet");
    put_line("a finite number of m-planes in (m+p)-space.");
    put_line("If the given m-planes are sufficiently generic, and");
    put_line("if their number equals m times p, then the number of"); 
    put_line("solutions is finite in number.  If this is the case,");
    put_line("then SAGBI and Pieri homotopies are optimal homotopies");
    put_line("to solve a Pieri problem: no solution path diverges.");
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
