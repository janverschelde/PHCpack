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

end Greeting_Banners;
