with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with Random_Conditioned_Evaluations;    use Random_Conditioned_Evaluations;
with VarbPrec_Polynomial_Evaluations;   use VarbPrec_Polynomial_Evaluations;

procedure ts_vmpeval is

-- DESCRIPTION :
--   Tests condition number computation of polynomial evaluation.
--   The problem is to evaluate a polynomial in several variables
--   at some point.

  procedure Standard_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : Standard_Complex_Polynomials.Poly;
    x : Standard_Complex_Vectors.Vector(1..integer32(n));
    rco : double_float;

  begin
    Random_Conditioned_Evaluation_Problem(n,d,m,c,cffsz,pntsz,close,p,x);
    put_line("The generated random polynomial :");
    put_line(p); new_line;
    rco := Inverse_Condition_Number(p,x);
    put("The inverse condition number :"); put(rco,3); new_line;
  end Standard_Test;

  procedure DoblDobl_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : DoblDobl_Complex_Polynomials.Poly;
    x : DoblDobl_Complex_Vectors.Vector(1..integer32(n));
    rco : double_double;

  begin
    Random_Conditioned_Evaluation_Problem(n,d,m,c,cffsz,pntsz,close,p,x);
    put_line("The generated random polynomial :");
    put_line(p); new_line;
    rco := Inverse_Condition_Number(p,x);
    put("The inverse condition number : "); put(rco,3); new_line;
  end DoblDobl_Test;

  procedure QuadDobl_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : QuadDobl_Complex_Polynomials.Poly;
    x : QuadDobl_Complex_Vectors.Vector(1..integer32(n));
    rco : quad_double;

  begin
    Random_Conditioned_Evaluation_Problem(n,d,m,c,cffsz,pntsz,close,p,x);
    put_line("The generated random polynomial :");
    put_line(p); new_line;
    rco := Inverse_Condition_Number(p,x);
    put("The inverse condition number : "); put(rco,3); new_line;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of a random polynomial.

    n,d,m : natural32 := 0;
    coeffsize,pointsize,closeroot,cond : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put_line("Condition numbers for polynomial evaluation ...");
    new_line;
    put_line("First part, parameters of the random polynomial :");
    put("Give number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    new_line;
    put_line("Second part, factors in numerical condition of evaluation :");
    put("Give magnitude of the coefficients : "); get(coeffsize);
    put("Give magnitude of the coordinates of the point : "); get(pointsize);
    put("Give closeness to a root : "); get(closeroot);
    cond := coeffsize*(pointsize**integer(d))/closeroot;
    put("Predicted condition number : "); put(cond,3); new_line;
    new_line;
    put_line("Choose the level of precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test(n,d,m,0,coeffsize,pointsize,closeroot);
      when '1' => DoblDobl_Test(n,d,m,0,coeffsize,pointsize,closeroot);
      when '2' => QuadDobl_Test(n,d,m,0,coeffsize,pointsize,closeroot);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_vmpeval;
