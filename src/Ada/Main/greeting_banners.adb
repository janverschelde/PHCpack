with text_io;                           use text_io;

package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.90 released 2024-03-20";

  begin
    return res;
  end Version;

  procedure How_to_Cite is
  begin
    put_line("To cite PHCpack in publications use:");
    new_line;
    put_line("Jan Verschelde:");
    put_line("Algorithm 795: PHCpack: A general-purpose solver for polynomial");
    put_line("systems by homotopy continuation.");
    put_line
      ("ACM Transactions on Mathematical Software, 25(2):251-276, 1999.");
  end How_to_Cite;

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
    new_line;
    put_line("Consider for example the system :");
    put_line("3");
    put_line(" (x-1)*(x^2 - y)*(x-0.5);");
    put_line(" (x-1)*(x^3 - z)*(y-0.5);");
    put_line(" (x-1)*(x*y - z)*(z-0.5);");
    new_line;
    put_line("which contains among its solutions a two dimensional plane,");
    put_line("the twisted cubic, and at least one isolated point.");
    put_line("Running 'phc -a' will produce three files, with suffixes");
    put_line("_w0, _w1, _w2 for the three dimensions in the solution set.");
    new_line;
    put_line("In a numerical irreducible decomposition, a k-dimensional set");
    put_line("is represented by a witness set, which consists of the input");
    put_line("system, augmented with k random linear equations, and the list");
    put_line("of generic points which satisfy the augmented system.");
    new_line;
    put_line("The output of 'phc -a' produces candidate witness sets, which");
    put_line("may still contains points on higher dimensional sets that need");
    put_line("filtering by the homotopy membership test in 'phc -f'.");
    put_line("After filtering, the factorization methods of 'phc -f' can be");
    put_line("applied to partition the witness points along the irreducible");
    put_line("factors in a numerical irreducible decomposition.");
    new_line;
    put_line("The equation-by-equation solver applies a bottom up method,");
    put_line("applying diagonal homotopies in 'phc -c'.");
    put_line("The top down method is available in 'phc -c' as well.");
  end help4eqnbyeqn;

  procedure help4blackbox is
  begin
    put_line("phc -b calls the blackbox solvers.");
    new_line;
    put_line("A typical call is 'phc -b input output' where 'input'");
    put_line("and 'output' are respectively the input and output files.");
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
    put_line("For these cases, no polynomial continuation is needed.");
    new_line;
    put_line("By default, phc -b solves in hardware double precision.");
    put_line("With phc -b2, double double arithmetic is applied,");
    put_line("and with -b4, path tracking happens in quad double arithmetic.");
    put_line("To offset the extra cost of higher precision arithmetic,");
    put_line("multitasking can be applied, e.g.: as in phc -b4 -t8");
    put_line("to solve a system in quad double precision with 8 tasks.");
    put_line("To use all available cores, type phc -b -t, without the #tasks.");
    new_line;
    put_line("The focus of phc -b is on isolated solutions and on systems");
    put_line("with as many equations as unknowns.  For more general systems"); 
    put_line("and to compute positive dimensional solution sets, see the");
    put_line("options -a and -c, respectively for an equation-by-equation");
    put_line("solver and for a numerical irreducible decomposition.");
  end help4blackbox;

  procedure help4compsolve is
  begin
    put_line("phc -B gives a blackbox numerical irreducible decomposition.");
    new_line;
    put_line("A numerical irreducible decomposition represents the solution");
    put_line("set of a polynomial system as a sequence of witness sets.");
    put_line("Each witness set in a numerical irreducible decomposition");
    put_line("represents one irreducible factor of the solution set.");
    new_line;
    put_line("A menu driven approach to compute a numerical irreducible");
    put_line("decompositon is provided via the options -c and -f.");
    put_line("In blackbox mode, the user does not have to know the various");
    put_line("stages in the computation and the tolerance settings.");
  end help4compsolve;

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
    new_line;
    put_line("The first menus offer a cascade of homotopies to compute");
    put_line("witness points on each positive dimensional solution set.");
    put_line("The user is prompted to enter the top dimension, which is");
    put_line("by default the number of variables minus one.  The cascade");
    put_line("runs faster if one would know the correct top dimension.");
    put_line("The output of the cascade needs filtering of points which");
    put_line("may lie on higher dimensionals set and factoring, see '-f'");
    new_line;
    put_line("The second part of the menu provides witness set intersection");
    put_line("with diagonal homotopies, in extrinsic coordinates (doubling)");
    put_line("or in intrinsic coordinates which avoids doubling #variables.");
    put_line("The witness set intersection is wrapped with '-w'.");
    new_line;
    put_line("The third functionality provided is the decomposition of the");
    put_line("solution set of binomial systems.  Solutions sets of a system");
    put_line("with exactly two monomials with nonzero coefficient in every");
    put_line("equation are represented by monomial maps, which is much more");
    put_line("efficient than a witness set data structure.");
  end help4components;

  procedure help4reduction is
  begin
    put_line("phc -d attempts to lower the degrees of the polynomials.");
    new_line;
    put_line("As an example, consider the intersection of two circles:");
    new_line;
    put_line("2");
    put_line(" x^2 + y^2 - 1;");
    put_line(" (x -  0.5)^2 + y^2 - 1;");
    new_line;
    put_line("A simple linear combination of the two polynomials gives:");
    new_line;
    put_line("2");
    put_line(" x^2 + y^2 - 1;");
    put_line(" x -  2.5E-1;");
    new_line;
    put_line("This reduced system has the same solutions, but only two");
    put_line("instead of four solution paths need to be tracked.");
    put_line("Nonlinear reduction attempts to replace higher degree");
    put_line("polynomials in the system by S-polynomials.");
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
    new_line;
    put_line("Pieri homotopies generalize SAGBI homotopies in two ways :");
    put_line("(1) the dimensions of the space of intersection may be");
    put_line("    higher than one, which reduces the number of solutions;");
    put_line("(2) the solutions may be curves of degree q which meet as");
    put_line("    many as m*p + q*(m+p) given m-planes at as many points.");
    new_line;
    put_line("An application of Pieri homotopies is the output pole place");
    put_line("placement problem in linear systems control, see phc -k.");
    new_line;
    put_line("Littlewood-Richardson homotopies solve general Schubert");
    put_line("problems.  On input is a sequence of square matrices.");
    put_line("With each matrix corresponds a bracket of conditions on");
    put_line("p-planes.  Each intersection condition is the dimension of");
    put_line("the intersection of a solution p-plane with a linear space");
    put_line("with generators in one of the matrices in the input sequence.");
  end help4enumeration;

  procedure help4feedback is
  begin
    put_line("phc -k computes feedback laws to control linear systems.");
    new_line;
    put_line("The two command line arguments are respectively the names");
    put_line("of the input and output files.  For examples of the input,");
    put_line("see the Feedback/Data folder of the source code.");
    put_line("The input consists of the (A, B, C) definition of a plant,");
    put_line("and for dynamic output, a number of poles.");
    new_line;
    put_line("The feedback laws in the output file are realized and written");
    put_line("in a format ready for parsing with MATLAB or Octave.");
    new_line;
    put_line("The computation of output feedback laws is an application");
    put_line("of the Pieri homotopies, available via phc -e.");
  end help4feedback;

  procedure help4factor is
  begin
    put_line("phc -f factors a solution set into irreducible components.");   
    new_line;
    put_line("The '-f' also stands for filter.  The first option allows");
    put_line("for example to extract only the real solutions from a file.");
    put_line("Other filtering criteria involve for example the residual,");
    put_line("the estimate of the inverse condition numbers, and a test");
    put_line("whether a coordinate of a solution is zero or not.");
    new_line;
    put_line("The homotopy membership test is needed to process the super");
    put_line("witness sets produced by 'phc -a' or '-c'.  Given a witness");
    put_line("set and a point, the test decides whether the point belongs");
    put_line("to the solution set represented by a witness set.");
    new_line;
    put_line("To factor a pure dimensional solution set, given as a witness");
    put_line("set, two different methods are available.  For small degrees,");
    put_line("a combinatorial factorization method may be fastest.");
    put_line("The other factorization run monodromy loops, using the linear");
    put_line("trace as stop criterion.");
    new_line;
    put_line("Another option of '-f' applies a tropical method to detect");
    put_line("a common factor of two Laurent polynomials.");
  end help4factor;

  procedure help4goodformat is
  begin
    put_line("phc -g checks whether the format of an input system is good.");
    new_line;
    put_line("Use '-g' as 'phc -g input output' where 'input' and 'output'");
    put_line("are the names of input and output files respectively.");
    put_line("If 'output' is omitted, then output can be written to screen.");
    put_line("If both 'input' and 'output' are omitted, then the user");
    put_line("is prompted to type the polynomials in the input system.");
    new_line;
    put_line("The first line of the input must start with the number of");
    put_line("polynomials in the system.  If the number of variables differs");
    put_line("from the number of polynomials, then the number of variables");
    put_line("must occur as the next number on the first line of the input.");
    new_line;
    put_line("Every polynomial must end with a semi colon, the ';' symbol.");
    put_line("Starting the name of variable with ';' is not good anyway,");
    put_line("other symbols which cannot be variable names are 'e' and 'E'");
    put_line("for the scientific notation of numbers.  Moreover, 'i' and");
    put_line("'I' denote the imaginary unit in complex coefficients.");
    new_line;
    put_line("Exponentiation is represented both by '^' and '**'");
    put_line("If negative exponents occur, then the system is parsed as a");
    put_line("Laurent polynomial system, see phc -m for polyhedral methods.");
    put_line("The division '/' may not appear in a monomial, e.g.: 'x/y',");
    put_line("but may be used in a coefficient, such as in '5/7'.");
    new_line;
    put_line("Round brackets are for grouping the real and imaginary parts");
    put_line("of complex coefficients, e.g.: (1.2e-3 + 3/7*I)*x^2*y");
    put_line("or for grouping factors, e.g.: 3.14*(x+y)*(x-1)^4.");
  end help4goodformat;

  procedure help4help is
  begin
    put_line("phc -h or phc --help tries to help...");
    new_line;
    put_line("For help on the options, type phc -h followed by");
    put_line(" -0 : random numbers with fixed seed for repeatable runs");
    put_line(" -a : solving polynomial systems equation-by-equation");
    put_line(" -b : batch or black-box processing, the blackbox solver");
    put_line(" -c : irreducible decomposition for solution components");
    put_line(" -d : linear and nonlinear reduction w.r.t. the total degree");
    put_line(" -e : SAGBI/Pieri/Littlewood-Richardson homotopies");
    put_line(" -f : factor pure dimensional solution set into irreducibles");
    put_line(" -g : checking whether an input system has the right syntax");
    put_line(" -k : realization of dynamic output feedback placing poles");
    put_line(" -l : witness set for hypersurface cutting with random line");
    put_line(" -m : mixed volume computation via lift+prune and MixedVol");
    put_line(" -o : write order of symbols after parsing polynomial system");
    put_line(" -p : polynomial continuation by a homotopy in one parameter");
    put_line(" -q : tracking solution paths with incremental read/write");  
    put_line(" -r : root counting and construction of start systems");
    put_line(" -s : equation and variable scaling on system and solutions");
    put_line(" -t : tasking for tracking paths using multiple threads");
    put_line(" -v : verification, refinement and purification of solutions");
    put_line(" -w : witness set intersection using diagonal homotopies");   
    put_line(" -x : convert solutions from PHCpack into Python dictionary"); 
    put_line(" -y : sample points from an algebraic set, given witness set");
    put_line(" -z : strip phc output solution lists into Maple format");
    new_line;
    put_line("For example, typing 'phc -h -z' or 'phc -z -h' (or '--help'");
    put_line("instead of '-h') displays information about 'phc -z'.");
  end help4help;

  procedure help4hypersurface is
  begin
    put_line("phc -l to compute a witness set for a hypersurface.");
    new_line;
    put_line("A hypersurface defined by a polynomial in many variables is cut");
    put_line("with one general line. The number of points on the hypersurface");
    put_line("and the general line equal the degree of the hypersurface.");
    put_line("This collection of points on the intersection of a hypersurface");
    put_line("and a general line form a witness set.");
    new_line;
    put_line("For example, if the file 'sphere' contains");
    new_line;
    put_line("1 3");
    put_line(" x^2 + y^2 + z^2 - 1;");
    new_line;
    put_line("then typing at the command line 'phc -l sphere output' produces");
    put_line("the file 'sphere_w2' with a witness set for the sphere.");
    new_line;
    put_line("By default, the witness points are computed in double precision");
    put_line("For double double and quad double precision,");
    put_line("use 'phc -l2' and 'phc -l4' respectively.");
    new_line;
    put_line("Related options are '-w' to intersect witness sets and");
    put_line("'-c' for more choices to compute witness sets.");
  end help4hypersurface;

  procedure help4mixvol is
  begin
    put_line("phc -m computes mixed volumes and runs polyhedral homotopies.");
    new_line;
    put_line("For polynomial systems with as many equations as unknowns,");
    put_line("the mixed volume of the Newton polytopes of the polynomials");
    put_line("gives a generically sharp upper bound on the number of");
    put_line("isolated solutions with their coordinates different from zero.");
    put_line("Polyhedral homotopies solve random coefficient start systems,");
    put_line("tracking as many solution paths as the mixed volume.");
    new_line;
    put_line("Methods to compute mixed volumes lift the polytopes,");
    put_line("adding one extra value to the coordinates of each point.");
    put_line("This extra lifting value will be the power of the continuation");
    put_line("parameter in the polyhedral homotopy.");
    new_line;
    put_line("Several lifting strategies are available :");
    put_line(" 0. static : lift and compute cells, the user may provide");
    put_line("    point-by-point lifting values to exploit structure;");
    put_line(" 1. implicit : based on a recursive formula for mixed volumes;");
    put_line(" 2. dynamic : points are considered in a sequence and lifting");
    put_line("    values are generated to preserve already computed cells;");
    put_line(" 3. symmetric : to exploit permutation symmetry,");
    put_line("    points in the same orbit receive the same lifting value;");
    put_line(" 4. MixedVol : a faster mixed volume calculator,");
    put_line("    which is applied in the blackbox solver.");
    put_line(" 5. DEMiCs : dynamic enumeration produces the mixed cells at a");
    put_line("    faster pace than MixedVol; is fastest for many polytopes.");
    put_line("For the Cayley trick to compute the Minkowski polynomial,");
    put_line("the dynamic lifting strategy should be selected.");
    new_line;
    put_line("Mixed volumes exclude solutions with zero coordinates.");
    put_line("With stable mixed volumes we count all affine solutions.");
    new_line;
    put_line("A related option is '-r' to compute bounds on the number of");
    put_line("isolated solutions based on the degrees of the polynomials.");
    put_line("The random coefficient start system constructed by '-m' can");
    put_line("be used in '-p' as start system to solve a specific system,");
    put_line("which has the same Newton polytopes as the start system.");
    put_line("If combined with '-t', followed by the number of tasks,");
    put_line("then multithreading is applied to solve the start system.");
    put_line("To use all available cores, type phc -m -t, without the #tasks.");
  end help4mixvol;

  procedure help4symbols is
  begin
    put_line("phc -o writes the symbol table after parsing an input system.");
    new_line;
    put_line("The order of the variables appears in the order in which the");
    put_line("the variables occur in the polynomial system.  For example,");
    put_line("1 2");
    put_line("  y + x + 1;");
    put_line("as saved in the file '/tmp/ex1' will result in the output");
    put_line("  y x");
    put_line("when running 'phc -o /tmp/ex1' or as 'phc -o /tmp/ex1 out'.");
    new_line;
    put_line("To force a particular order in the symbol table,");
    put_line("One can add to the first polynomial 'x - x + y - y + ...'");
    put_line("While 'x - x + y - y' results in zero, its effect is that it");
    put_line("initializes the symbol table with 'x' and 'y', in that order.");
  end help4symbols;

  procedure help4continuation is
  begin
    put_line("phc -p runs continuation with a homotopy in one parameter.");
    new_line;
    put_line("We have two types of homotopies:");
    put_line("(1) Using an artificial parameter t in a convex combination");
    put_line("    between a target and start system; or");
    put_line("(2) Using a natural parameter, using a variable as parameter,");
    put_line("    a variable which occurs in the given input system.");
    new_line;
    put_line("In an artificial parameter homotopy, the t is incremented");
    put_line("by the predictor and remains fixed in the corrector.");
    put_line("In a natural parameter homotopy, arc length continuation");
    put_line("is applied and the variable used as continuation parameter");
    put_line("remains a variable in the corrector stage.");
    new_line;
    put_line("The start system can be constructed with '-r' or '-m'.");
    put_line("The working precision can be adjusted in the interactive menu");
    put_line("or at the command line via '-p2' for double double,");
    put_line("or via '-p4' for quad double precision.");
  end help4continuation;

  procedure help4jumpstart is
  begin
    put_line("phc -q runs continuation with jump starting for huge systems.");
    new_line;
    put_line("For huge systems, all solutions may not fit in memory.");
    put_line("Jump starting a homotopy means that start solutions will be");
    put_line("computed or retrieved from file when they are needed.");
    new_line;
    put_line("The following homotopies are supported :");
    put_line(" 1. start system based on total degree, see 'phc -r';");
    put_line(" 2. a linear-product start system, see 'phc -r';");
    put_line(" 3. cheater: start system and solutions are on file;");
    put_line(" 4. polyhedral continuation to solve a generic system,");
    put_line("    given a mixed cell configuration computed with 'phc -m';");
    put_line(" 5. diagonal homotopies to intersect sets, see 'phc -c';");
    put_line(" 6. cascade homotopies to compute witness sets, see 'phc -c'.");
    new_line;
    put_line("With 'phc -q' we invert the order of solving the start system");
    put_line("and tracking all paths ('phc -r' followed by 'phc -p'),");
    put_line("as with '-q' we compute or retrieve a start solution only");
    put_line("when needed in the tracking of the next solution path.");
  end help4jumpstart;

  procedure help4adepath is
  begin
    put_line("phc -j tracks solution paths with algorithmic differentiation.");
    new_line;
    put_line("Instead of computing all partial derivatives separately,");
    put_line("the Jacobian matrix can be evaluated more efficiently");
    put_line("using methods of algorithmic differentiation.");
  end help4adepath;

  procedure help4rootcounts is
  begin
    put_line("phc -r computes root counts and constructs start systems.");
    new_line;
    put_line("The root count determines the number of solution paths that");
    put_line("are tracked in a homotopy connecting the input system with");
    put_line("the start system that has as many solutions as the root count.");
    put_line("A homotopy is optimal to solve a given system if the number of");
    put_line("solution paths equals the number of solutions of the system.");
    new_line;
    put_line("For root counts in numerical Schubert calculus, see 'phc -e'.");
    put_line("Mixed volumes are a subset of 'phc -r', but do 'phc -m -h'");
    put_line("for a brief overview of the lifting strategies.");
    new_line;
    put_line("Based on the degrees, we have four different upper bounds :");
    put_line(" 0. the product of the degrees of the polynomials in the");
    put_line("    system corresponds to the classical Bezout bound;");
    put_line(" 1. a multi-homogeneous Bezout number based on a partition");
    put_line("    of the set of variables that occur in the system;");
    put_line(" 2. a partitioned linear-product Bezout number where every");
    put_line("    polynomial can have a different partition;");
    put_line(" 3. a general linear-product Bezout number based on a set");
    put_line("    structure where every set corresponds to a hyperplane.");
    put_line("Upper bounds 1, 2, and 3 give rise to a linear-product start");
    put_line("system.  Via an interactive menu, the user can try several");
    put_line("bounds before settling on a choice of a start system.");
    new_line;
    put_line("There is support for a symmetric set structure.");
    new_line;
    put_line("The start system can be given as input to 'phc -p' or '-q'");
    put_line("to compute all isolated solutions of the input system.");
    put_line("If after selecting the type of start system, the computation");
    put_line("of all start solutions is delayed, then '-q' must be used.");
    put_line("To use '-p' all start solutions must be computed by 'phc -r'.");
  end help4rootcounts;

  procedure help4scaling is
  begin
    put_line("phc -s applies equation and variable scaling to a system.");
    new_line;
    put_line("A badly scaled system has tiny and huge coefficients.");
    put_line("Solutions of a badly scaled system are ill conditioned.");
    new_line;
    put_line("To transform a badly scaled system into a better one,");
    put_line("equation scaling multiplies the coefficients in the same");
    put_line("polynomial by the same factor and variable scaling multiplies");
    put_line("all occurrences of the same variable by the same factor.");
    put_line("To center the coefficients and reduce their variability,");
    put_line("the scaling factors are obtained via an optimization problem.");
    new_line;
    put_line("The output are a scaled system and the scaling factors.");
    put_line("These factors are needed to transform the computed solutions");
    put_line("back into their original coordinates.");
  end help4scaling;

  procedure help4tasking is
  begin
    put_line("phc -t# computes with # tasks, for shared memory parallelism.");
    put_line("To use all available cores, type phc -t, without the #tasks.");
    new_line;
    put_line("Because the solution paths in a polynomial homotopy can be");
    put_line("tracked independently from each other, the problem of tracking");
    put_line("a number of paths is a pleasingly parallel problem.");
    put_line("Already on relatively small systems, with a modest number");
    put_line("of threads, a significant speedup can be obtained.");
    new_line;
    put_line("To apply multithreading to the blackbox solver,");
    put_line("type for example 'phc -b -t8 input output' and 8 tasks will be");
    put_line("used in the path tracking.  The obtained speedup can offset");
    put_line("the overhead caused by the cost of higher precision.");
    put_line("For example, to compensate for the cost of double double");
    put_line("precision arithmetic, do 'phc -b2 -t16 input output' and");
    put_line("if sufficiently many cores are available, one may observe");
    put_line("a 'quality up', obtaining more accurate results faster.");
    new_line;
    put_line("The option -t can be added to phc -p and -m to track solution");
    put_line("path defined by a polynomial homotopy with many threads.");
  end help4tasking;

  procedure help4series is
  begin
    put_line("phc -u uses power series, Pade approximants for path tracking.");
    new_line;
    put_line("A power series solution expresses the coordinates of the");
    put_line("solution as a truncated power series in one variable.");
    put_line("That one variable is typically one of the other variables in");
    put_line("the system.  The system on input defines the polynomials for");
    put_line("a space curve.  Given any point on the space curve, the power");
    put_line("series provide an approximation for points along that curve,");
    put_line("similar to a Taylor series expansion for any function.");
    new_line;
    put_line("The input to phc -u is a system and its isolated solutions.");
    put_line("The format for the input is similar to that of a start system.");
    put_line("If no solutions are available for the system, then launching");
    put_line("the blackbox solver phc -b will append solutions to the input.");
    put_line("The solutions provided should be for a zero value of the one");
    put_line("variable of the power series.");
    new_line;
    put_line("The coefficients of the power series are complex numbers,");
    put_line("by default in standard double precision.  Also double double");
    put_line("and quad double power expansions are available, respectively");
    put_line("by calling phc -u2 and phc -u4.");
    new_line;
    put_line("After reading the input file, phc -u prompts the user for the");
    put_line("index of the series parameter.  This index points to one of");
    put_line("the variables in the input system.  The next prompt is for the");
    put_line("number of Newton updates to the solution.");
    put_line("If the maximum degree of the power series is set to 16,");
    put_line("by the quadratic convergence of Newton's method, starting at");
    put_line("degree 1, five stepts will suffice to compute all coefficients.");
    put_line("For a regular solution, the quadratic convergence can be seen");
    put_line("when answering 'y' to the question on extra diagnostic output.");
    new_line;
    put_line("With power series, Pade approximants are computed to predict");
    put_line("the next point on a solution path in a polynomial homotopy.");
  end help4series;

  procedure help4verification is
  begin
    put_line("phc -v filters, verifies, and refines lists of solutions.");
    new_line;
    put_line("Every solution terminates with a line such as");
    put_line("'== err :  3.374E-09 = rco :  2.326E-06 = res :  3.613E-16 =='");
    put_line("where the quality indicators are as follows:");
    put_line("* err : the magnitude of the last correction of Newton's method");
    put_line("is the forward error on the solution, '3E-9' in the example");
    put_line("above means that we may hope to have 8 correct decimal places;");
    put_line("* rco : an estimate for the inverse condition number of the");
    put_line("Jacobian matrix evaluated at the approximate solution, which"); 
    put_line("measures how sensitive the solution is to changes in the input");
    put_line("coefficients, '2E-6' in the example, indicated that errors on");
    put_line("the input can be multiplied by a factor of 10^6;");
    put_line("* res : the magnitude of the approximate solution evaluated at");
    put_line("the input polynomials is the backward error, i.e.: by much");
    put_line("should the input coefficients change so that the approximate");
    put_line("solution is an exact solution of the changed problem.");
    new_line;
    put_line("The main result of a basic verification is the tally of");
    put_line("good solutions, versus solutions at infinity and/or singular");
    put_line("solutions.  Solution paths may also have ended at points");
    put_line("that are clustered at a regular solution so with '-v' we can");
    put_line("detect some cases of occurrences of path crossing.");
    put_line("Applying Newton's method in a higher precision may resolve");
    put_line("the distinction between regular and singular solutions.");
    put_line("For isolated singular solutions, deflation may restore the");
    put_line("quadratic convergence of Newton's method.");
    new_line;
    put_line("The blackbox version uses default settings for the parameters,");
    put_line("use as phc -v -b or phc -b -v, for double precision.");
    put_line("For double double precision, use as phc -b2 -v or phc -b -v2.");
    put_line("For quad double precision, use as phc -b4 -v or phc -b -v4.");
    put_line("The order of -b and -v at the command line does not matter.");
    new_line;
    put_line("To select solutions from a list subject to certain criteria,");
    put_line("use phc -f.");
  end help4verification;

  procedure help4verbose is
  begin
    put_line("phc -V# provides verbose information for debugging code.");
    new_line;
    put_line("The number # following the -V should be a positive integer,");
    put_line("otherwise the verbose option has no effect.");
    new_line;
    put_line("For example, phc -b -V9 input output, will display extra");
    put_line("internal data, listing the procedures called during the");
    put_line("running of the blackbox solver.");
    new_line;
    put_line("More information is shown for larger verbose levels.");
  end help4verbose;

  procedure help4witsetinsect is
  begin
    put_line("phc -w for witness set intersection with diagonal homotopies.");
    new_line;
    put_line("This option wraps the diagonal homotopies to intersect two");
    put_line("positive dimensional solution sets given as witness sets.");
    new_line;
    put_line("Consider for example a sphere and a cylinder in 3-space.");
    put_line("Their intersection defines a quartic, a curve of degree 4.");
    put_line("If the file sphere_w2 contains a witness set for the sphere");
    put_line("and if a witness set for the cylinder is contained in the");
    put_line("file cylinder_w2, then we can compute their intersection as");
    new_line;
    put_line("phc -w sphere_w2 cylinder_w2 quartic");
    new_line;
    put_line("Note that the symbols in both input files must line up.");
    put_line("If one of the variables would not be needed in the natural");
    put_line("formulation of a set, e.g.: such as in 'x^2 + y^2 - 1' for a");
    put_line("vertical cylinder in the (x, y, z)-space, then the z variable");
    put_line("must be added in an artificial manner, e.g., adding: 'z - z'.");
    new_line;
    put_line("To compute a witness set for a hypersurface, use phc -l.");
    put_line("For more choices to intersect witness sets, see phc -c.");
  end help4witsetinsect;

  procedure help4pythondict is
  begin
    put_line("phc -x converts lists of solutions into Python dictionaries.");
    new_line;
    put_line("To work with solution lists in Python scripts, phc -x converts");
    put_line("a solution list in PHCpack format to a list of dictionaries.");
    put_line("Given a Python list of dictionaries, phc -x returns a list of");
    put_line("solutions in PHCpack format.");
    new_line;
    put_line("The Python format of a list of solutions is a list, enclosed by");
    put_line("square brackets.  Every solution is a dictionary.  The names of");
    put_line("the variables are keys to the values of the coordinates.");
    new_line;
    put_line("At the command line, type 'phc -x input output' where 'input'");
    put_line("contains a list of solutions in PHCpack or in Python format.");
    put_line("The result of the conversion (from PHCpack into Python format");
    put_line("or from Python to PHCpack format, depending on the input type)");
    put_line("is written to the file 'output'.");
    put_line("If 'output' is omitted, then the result is written to screen.");
  end help4pythondict;

  procedure help4sampler is
  begin
    put_line("phc -y samples points on a positive dimensional solution set.");
    new_line;
    put_line("The points on a positive dimensional solution set are fixed by");
    put_line("the position of hyperplanes that define a linear space of the");
    put_line("dimension equal to the co-dimension of the solution set.");
    put_line("For example, in 3-space, a 2-dimensional set is cut with");
    put_line("a line and a 1-dimensional set is cut with a plane.");
    new_line;
    put_line("If 'sphere_w2' contains a witness set for a sphere, then");
    put_line("typing at the command prompt 'phc -y sphere_w2 new_sphere'");
    put_line("will make a new witness set in the file 'new_sphere_w2'.");
    new_line;
    put_line("Use 'phc -c' to compute a first witness set,");
    put_line("or 'phc -l' for the particular case of a hypersurface.");
  end help4sampler;

  procedure help4mapleform is
  begin
    put_line("phc -z converts lists of solutions into Maple format.");
    new_line;
    put_line("The Maple format of a solution list is a list of lists,");
    put_line("where each list is enclosed by square brackets.");
    put_line("The items in a list are separated by commas.");
    put_line("A solution is represented as a list of equations.");
    new_line;
    put_line("At the command line, type 'phc -z input output' where 'input'");
    put_line("contains a list of solutions, in PHCpack or in Maple format.");
    put_line("The result of the conversion (from PHCpack into Maple format");
    put_line("or from Maple to PHCpack format, depending on the input format)");
    put_line("is written to the file 'output'.");
    put_line("If 'output' is omitted, then the result is written to screen.");
    new_line;
    put_line("A related option is phc -z -p or phc -p -z to extract the");
    put_line("start system and solution from the phc -b output file."); 
  end help4mapleform;

  procedure help4getstart is
  begin
    put_line("phc -z -p returns the start system from the phc -b output.");
    new_line;
    put_line("The output file of the blackbox solver contains the start");
    put_line("system and start solutions used to solve a system.");
    put_line("If a system with the same structure but other coefficients");
    put_line("needs to be solved, then the start system computed by phc -b");
    put_line("can be used as the start system for phc -p or phc -u.");
    new_line;
    put_line("At the command line, type 'phc -z -p input output'");
    put_line("where 'input' is the name of the output file of 'phc -b'");
    put_line("If 'output' is omitted, then 'phc -z -p' writes to screen.");
  end help4getstart;

end Greeting_Banners;
