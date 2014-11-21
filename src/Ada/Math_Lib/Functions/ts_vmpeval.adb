with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Mathematical_Functions;   use Multprec_Mathematical_Functions;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Multprec_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Multprec_Complex_Matrices;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;   use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Strings;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;   use DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laurentials_io;   use DoblDobl_Complex_Laurentials_io;
with DoblDobl_Complex_Laur_Strings;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Functions;   use QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laurentials_io;   use QuadDobl_Complex_Laurentials_io;
with QuadDobl_Complex_Laur_Strings;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;   use Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laurentials_io;   use Multprec_Complex_Laurentials_io;
with Multprec_Complex_Laur_Strings;
with Random_Conditioned_Matrices;       use Random_Conditioned_Matrices;
with Random_Conditioned_Evaluations;    use Random_Conditioned_Evaluations;
with Varbprec_Complex_Linear_Solvers;   use Varbprec_Complex_Linear_Solvers; 
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
    g : Standard_Complex_Vectors.Vector(1..integer32(n))
      := Standard_Random_Vectors.Random_Vector(1,integer32(n));
    rco : double_float;

    use Standard_Complex_Numbers;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,cffsz,pntsz,close,g,p,x);
    put_line("The generated random polynomial :");
    put_line(p); new_line;
    put_line("Its random gradient vector :");
    for k in g'range loop
      declare
        dpk : Standard_Complex_Polynomials.Poly
            := Standard_Complex_Polynomials.Diff(p,k);
        val : Complex_Number := Eval(dpk,x);
      begin
        put("at k = "); put(k,1); put_line(" :");
        put(val); new_line;
        put(g(k)); new_line;
        Standard_Complex_Polynomials.Clear(dpk);
      end;
    end loop;
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
    g : DoblDobl_Complex_Vectors.Vector(1..integer32(n))
      := DoblDobl_Random_Vectors.Random_Vector(1,integer32(n));
    rco : double_double;

    use DoblDobl_Complex_Numbers;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,cffsz,pntsz,close,g,p,x);
    put_line("The generated random polynomial :");
    put_line(p); new_line;
    put_line("Its random gradient vector :");
    for k in g'range loop
      declare
        dpk : DoblDobl_Complex_Polynomials.Poly
            := DoblDobl_Complex_Polynomials.Diff(p,k);
        val : Complex_Number := Eval(dpk,x);
      begin
        put("at k = "); put(k,1); put_line(" :");
        put(val); new_line;
        put(g(k)); new_line;
        DoblDobl_Complex_Polynomials.Clear(dpk);
      end;
    end loop;
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
    g : QuadDobl_Complex_Vectors.Vector(1..integer32(n))
      := QuadDobl_Random_Vectors.Random_Vector(1,integer32(n));
    rco : quad_double;

    use QuadDobl_Complex_Numbers;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,cffsz,pntsz,close,g,p,x);
    put_line("The generated random polynomial :");
    put_line(p); new_line;
    put_line("Its random gradient vector :");
    for k in g'range loop
      declare
        dpk : QuadDobl_Complex_Polynomials.Poly
            := QuadDobl_Complex_Polynomials.Diff(p,k);
        val : Complex_Number := Eval(dpk,x);
      begin
        put("at k = "); put(k,1); put_line(" :");
        put(val); new_line;
        put(g(k)); new_line;
        QuadDobl_Complex_Polynomials.Clear(dpk);
      end;
    end loop;
    rco := Inverse_Condition_Number(p,x);
    put("The inverse condition number : "); put(rco,3); new_line;
  end QuadDobl_Test;

  procedure Multprec_Test
              ( n,d,m,c,sz : in natural32;
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
  --   sz       the size of the numbers;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : Multprec_Complex_Polynomials.Poly;
    x : Multprec_Complex_Vectors.Vector(1..integer32(n));
    g : Multprec_Complex_Vectors.Vector(1..integer32(n))
      := Multprec_Random_Vectors.Random_Vector(1,integer32(n),sz);
    rco : Floating_Number;

    use Multprec_Complex_Numbers;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,sz,cffsz,pntsz,close,g,p,x);
    put_line("The generated random polynomial :");
    put_line(p); new_line;
    put_line("Its random gradient vector :");
    for k in g'range loop
      declare
        dpk : Multprec_Complex_Polynomials.Poly
            := Multprec_Complex_Polynomials.Diff(p,k);
        val : Complex_Number := Eval(dpk,x);
      begin
        put("at k = "); put(k,1); put_line(" :");
        put(val); new_line;
        put(g(k)); new_line;
        Multprec_Complex_Polynomials.Clear(dpk);
        Clear(val);
      end;
    end loop;
    rco := Inverse_Condition_Number(p,x);
    put("The inverse condition number : "); put(rco,3); new_line;
  end Multprec_Test;

  procedure Test_Random_Polynomial is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of a random polynomial.

    n,d,m,deci,size : natural32 := 0;
    coeffsize,pointsize,closeroot,cond : double_float := 0.0;
    ans : character;

  begin
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
    put_line("  2. quad double precision;");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Standard_Test(n,d,m,0,coeffsize,pointsize,closeroot);
      when '1' => DoblDobl_Test(n,d,m,0,coeffsize,pointsize,closeroot);
      when '2' => QuadDobl_Test(n,d,m,0,coeffsize,pointsize,closeroot);
      when '3' =>
        new_line;
        deci := 2*natural32(log10(cond));
        put("The number of decimal places : "); put(deci,1); new_line;
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        Multprec_Test(n,d,m,0,size,coeffsize,pointsize,closeroot);
      when others => null;
    end case;
  end Test_Random_Polynomial;

  procedure Standard_Test_Taschini_Laurential is

  -- DESCRIPTION :
  --   Parses a particular polynomial used to motivate interval arithmetic,
  --   using standard double precision arithmetic.
  --   The example is taken from the paper of Stefano Taschini:
  --   "Interval Arithmetic: Python Implementation and Applications"
  --   Proceedings of the 7th Python in Science Conference (SciPy 2008).
  --   While we get the final correct approximation when we use 36
  --   decimal places using mpmath of SymPy, the choice of 36 is not obvious.

    use Standard_Complex_Laurentials;

    n : constant natural32 := 2;
    s : constant string
      := "(333.75 - x**2)*y**6 + x**2*(11*x**2*y**2 - 121*y**4 - 2) "
       & "+ 5.5*y**8 + (1/2)*x*y^-1;";
    p : Poly := Standard_Complex_Laur_Strings.Parse(n,s);
    x : Standard_Complex_Vectors.Vector(1..integer32(n));
    rco,absfz,denrco : double_float;
    val : Standard_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put_line("The test polynomial : ");
    put_line(s);
    put_line("The parsed polynomial :"); put(p); new_line;
    x(1) := Standard_Complex_Numbers.Create(77617.0);
    x(2) := Standard_Complex_Numbers.Create(33096.0);
    for i in x'range loop
      put("x("); put(i,1); put(") : "); put(x(i)); new_line;
    end loop;
    rco := Inverse_Condition_Number(p,x);
    put("The inverse of the condition number : "); put(rco,3); new_line;
    Evaluate_with_Inverse_Condition(p,x,absfz,denrco,rco,val);
    put("Condition number without residual : "); put(denrco,3); new_line;
    put("The value : "); put(val); new_line;
  end Standard_Test_Taschini_Laurential;

  procedure DoblDobl_Test_Taschini_Laurential is

  -- DESCRIPTION :
  --   Parses a particular polynomial used to motivate interval arithmetic,
  --   using double double precision arithmetic.
  --   The example is taken from the paper of Stefano Taschini:
  --   "Interval Arithmetic: Python Implementation and Applications"
  --   Proceedings of the 7th Python in Science Conference (SciPy 2008).
  --   While we get the final correct approximation when we use 36
  --   decimal places using mpmath of SymPy, the choice of 36 is not obvious.

    use DoblDobl_Complex_Laurentials;

    n : natural32 := 2;
    s : constant string
      := "(333.75 - x**2)*y**6 + x**2*(11*x**2*y**2 - 121*y**4 - 2) "
       & "+ 5.5*y**8 + (1/2)*x*y^-1;";
    p : Poly := DoblDobl_Complex_Laur_Strings.Parse(n,s);
    x : DoblDobl_Complex_Vectors.Vector(1..integer32(n));
    rco,absfz,denrco : double_double;
    dd_x1 : constant double_double := create(77617.0);
    dd_x2 : constant double_double := create(33096.0);
    val : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put_line("The test polynomial : ");
    put_line(s);
    put_line("The parsed polynomial :"); put(p); new_line;
    x(1) := DoblDobl_Complex_Numbers.Create(dd_x1);
    x(2) := DoblDobl_Complex_Numbers.Create(dd_x2);
    for i in x'range loop
      put("x("); put(i,1); put(") : "); put(x(i)); new_line;
    end loop;
    put_line("The point where to evaluate :"); put_line(x);
    rco := Inverse_Condition_Number(p,x);
    put("The inverse of the condition number : "); put(rco,3); new_line;
    Evaluate_with_Inverse_Condition(p,x,absfz,denrco,rco,val);
    put("Condition number without residual : "); put(denrco,3); new_line;
    put("The value : "); put(val); new_line;
  end DoblDobl_Test_Taschini_Laurential;

  procedure QuadDobl_Test_Taschini_Laurential is

  -- DESCRIPTION :
  --   Parses a particular polynomial used to motivate interval arithmetic,
  --   using double double precision arithmetic.
  --   The example is taken from the paper of Stefano Taschini:
  --   "Interval Arithmetic: Python Implementation and Applications"
  --   Proceedings of the 7th Python in Science Conference (SciPy 2008).
  --   While we get the final correct approximation when we use 36
  --   decimal places using mpmath of SymPy, the choice of 36 is not obvious.

    use QuadDobl_Complex_Laurentials;

    n : natural32 := 2;
    s : constant string
      := "(333.75 - x**2)*y**6 + x**2*(11*x**2*y**2 - 121*y**4 - 2) "
       & "+ 5.5*y**8 + (1/2)*x*y^-1;";
    p : Poly := QuadDobl_Complex_Laur_Strings.Parse(n,s);
    x : QuadDobl_Complex_Vectors.Vector(1..integer32(n));
    rco,absfz,denrco : quad_double;
    qd_x1 : constant quad_double := create(77617.0);
    qd_x2 : constant quad_double := create(33096.0);
    val : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put_line("The test polynomial : ");
    put_line(s);
    put_line("The parsed polynomial :"); put(p); new_line;
    x(1) := QuadDobl_Complex_Numbers.Create(qd_x1);
    x(2) := QuadDobl_Complex_Numbers.Create(qd_x2);
    for i in x'range loop
      put("x("); put(i,1); put(") : "); put(x(i)); new_line;
    end loop;
    put_line("The point where to evaluate :"); put_line(x);
    rco := Inverse_Condition_Number(p,x);
    put("The inverse of the condition number : "); put(rco,3); new_line;
    Evaluate_with_Inverse_Condition(p,x,absfz,denrco,rco,val);
    put("Condition number without residual : "); put(denrco,3); new_line;
    put("The value : "); put(val); new_line;
  end QuadDobl_Test_Taschini_Laurential;

  procedure Multprec_Test_Taschini_Laurential ( dp : in natural32 ) is

  -- DESCRIPTION :
  --   Parses a particular polynomial used to motivate interval arithmetic,
  --   using double double precision arithmetic.
  --   The example is taken from the paper of Stefano Taschini:
  --   "Interval Arithmetic: Python Implementation and Applications"
  --   Proceedings of the 7th Python in Science Conference (SciPy 2008).
  --   While we get the final correct approximation when we use 36
  --   decimal places using mpmath of SymPy, the choice of 36 is not obvious.

    use Multprec_Complex_Laurentials;

    size : natural32 := Multprec_Floating_Numbers.Decimal_to_Size(dp);
    n : natural32 := 2;
    s : constant string
      := "(333.75 - x**2)*y**6 + x**2*(11*x**2*y**2 - 121*y**4 - 2) "
       & "+ 5.5*y**8 + (1/2)*x*y^-1;";
    p : Poly := Multprec_Complex_Laur_Strings.Parse(n,size,s);
    x : Multprec_Complex_Vectors.Vector(1..integer32(n));
    rco,absfz,denrco : Floating_Number;
    qd_x1 : Floating_Number := create(77617.0);
    qd_x2 : Floating_Number := create(33096.0);
    val : Multprec_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put_line("The test polynomial : ");
    put_line(s);
    put_line("The parsed polynomial :"); put(p); new_line;
    x(1) := Multprec_Complex_Numbers.Create(qd_x1);
    x(2) := Multprec_Complex_Numbers.Create(qd_x2);
    for i in x'range loop
      put("x("); put(i,1); put(") : "); put(x(i)); new_line;
    end loop;
    put_line("The point where to evaluate :"); put_line(x);
    rco := Inverse_Condition_Number(p,x);
    put("The inverse of the condition number : "); put(rco,3); new_line;
    Evaluate_with_Inverse_Condition(p,x,absfz,denrco,rco,val);
    put("Condition number without residual : "); put(denrco,3); new_line;
    put("The value : "); put(val); new_line;
  end Multprec_Test_Taschini_Laurential;

  procedure Test_Taschini_Laurential is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables and the precision.

    n : natural32 := 0;
    ans : character;

  begin
   -- new_line;
   -- put("Give the number of variables : "); get(n);
    n := 2;
    Symbol_Table.Init(n);
    new_line;
    put_line("MENU to select the precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Standard_Test_Taschini_Laurential;
      when '1' => DoblDobl_Test_Taschini_Laurential;
      when '2' => QuadDobl_Test_Taschini_Laurential;
      when '3' => 
        declare
          deci : natural32 := 0;
        begin
          new_line;
          put("Give the number of decimal places : "); get(deci);
          Multprec_Test_Taschini_Laurential(deci);
        end;
      when others => null;
    end case;
  end Test_Taschini_Laurential;

  procedure Standard_Check_Values
              ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
                x : in Standard_Complex_Vectors.Vector;
                jm : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes the residual |f(x)| and checks whether the Jacobian
  --   matrix of f evaluated at x equals the matrix jm.

    y : Standard_Complex_Vectors.Vector
      := Standard_Complex_Poly_SysFun.Eval(f,x);
    res : double_float := Max_Norm(y);
    jmf : Standard_Complex_Jaco_Matrices.Jaco_Mat(f'range,x'range)
        := Standard_Complex_Jaco_Matrices.Create(f);
    mtx : Standard_Complex_Matrices.Matrix(f'range,x'range)
        := Standard_Complex_Jaco_Matrices.Eval(jmf,x);
 
  begin
    new_line;
    put_line("Checking the values ...");
    put("-> the residual :"); put(res,3); new_line;
    for i in mtx'range(1) loop
      for j in mtx'range(2) loop
        put("i = "); put(i,1); put("  j = "); put(j,1);
        put_line(" :");
        put(mtx(i,j)); new_line;
        put(jm(i,j)); new_line;
      end loop;
    end loop;
  end Standard_Check_Values;

  procedure DoblDobl_Check_Values
              ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                x : in DoblDobl_Complex_Vectors.Vector;
                jm : in DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes the residual |f(x)| and checks whether the Jacobian
  --   matrix of f evaluated at x equals the matrix jm.

    y : DoblDobl_Complex_Vectors.Vector
      := DoblDobl_Complex_Poly_SysFun.Eval(f,x);
    res : double_double := Max_Norm(y);
    jmf : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(f'range,x'range)
        := DoblDobl_Complex_Jaco_Matrices.Create(f);
    mtx : DoblDobl_Complex_Matrices.Matrix(f'range,x'range)
        := DoblDobl_Complex_Jaco_Matrices.Eval(jmf,x);
 
  begin
    new_line;
    put_line("Checking the values ...");
    put("-> the residual :"); put(res,3); new_line;
    for i in mtx'range(1) loop
      for j in mtx'range(2) loop
        put("i = "); put(i,1); put("  j = "); put(j,1);
        put_line(" :");
        put(mtx(i,j)); new_line;
        put(jm(i,j)); new_line;
      end loop;
    end loop;
  end DoblDobl_Check_Values;

  procedure QuadDobl_Check_Values
              ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                x : in QuadDobl_Complex_Vectors.Vector;
                jm : in QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes the residual |f(x)| and checks whether the Jacobian
  --   matrix of f evaluated at x equals the matrix jm.

    y : QuadDobl_Complex_Vectors.Vector
      := QuadDobl_Complex_Poly_SysFun.Eval(f,x);
    res : quad_double := Max_Norm(y);
    jmf : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(f'range,x'range)
        := QuadDobl_Complex_Jaco_Matrices.Create(f);
    mtx : QuadDobl_Complex_Matrices.Matrix(f'range,x'range)
        := QuadDobl_Complex_Jaco_Matrices.Eval(jmf,x);
 
  begin
    new_line;
    put_line("Checking the values ...");
    put("-> the residual :"); put(res,3); new_line;
    for i in mtx'range(1) loop
      for j in mtx'range(2) loop
        put("i = "); put(i,1); put("  j = "); put(j,1);
        put_line(" :");
        put(mtx(i,j)); new_line;
        put(jm(i,j)); new_line;
      end loop;
    end loop;
  end QuadDobl_Check_Values;

  procedure Multprec_Check_Values
              ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
                x : in Multprec_Complex_Vectors.Vector;
                jm : in Multprec_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes the residual |f(x)| and checks whether the Jacobian
  --   matrix of f evaluated at x equals the matrix jm.

    y : Multprec_Complex_Vectors.Vector
      := Multprec_Complex_Poly_SysFun.Eval(f,x);
    res : Floating_Number := Max_Norm(y);
    jmf : Multprec_Complex_Jaco_Matrices.Jaco_Mat(f'range,x'range)
        := Multprec_Complex_Jaco_Matrices.Create(f);
    mtx : Multprec_Complex_Matrices.Matrix(f'range,x'range)
        := Multprec_Complex_Jaco_Matrices.Eval(jmf,x);
 
  begin
    new_line;
    put_line("Checking the values ...");
    put("-> the residual :"); put(res,3); new_line;
    for i in mtx'range(1) loop
      for j in mtx'range(2) loop
        put("i = "); put(i,1); put("  j = "); put(j,1);
        put_line(" :");
        put(mtx(i,j)); new_line;
        put(jm(i,j)); new_line;
      end loop;
    end loop;
  end Multprec_Check_Values;

  procedure Standard_Jacobian_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+16.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Standard_Complex_Vectors.Vector(p'range);
    jm : Standard_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : double_float;
    loss_jac,loss_eva : integer32;
    ans : character;

  begin
    new_line;
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if rco = 0.0
     then loss_eva := -2**30;
     else loss_eva := integer32(log10(rco));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    put("Check values of the Jacobian matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Check_Values(p,x,jm);
    end if;
  end Standard_Jacobian_Test;

  procedure DoblDobl_Jacobian_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+32.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : DoblDobl_Complex_Vectors.Vector(p'range);
    jm : DoblDobl_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : double_double;
    loss_jac,loss_eva : integer32;
    ans : character;

  begin
    new_line;
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Is_Zero(rco)
     then loss_eva := -2**30;
     else loss_eva := integer32(to_double(log10(rco)));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    put("Check the values of the Jacobian matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Check_Values(p,x,jm);
    end if;
  end DoblDobl_Jacobian_Test;

  procedure QuadDobl_Jacobian_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+64.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : QuadDobl_Complex_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : quad_double;
    loss_jac,loss_eva : integer32;
    ans : character;

  begin
    new_line;
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Is_Zero(rco)
     then loss_eva := -2**30;
     else loss_eva := integer32(to_double(log10(rco)));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    put("Check the values of the Jacobian matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Check_Values(p,x,jm);
    end if;
  end QuadDobl_Jacobian_Test;

  procedure Multprec_Jacobian_Test
              ( n,d,m,c,sz : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed, provided the value sz for the precision
  --   is sufficiently high.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   sz       size of the numbers in the working precision;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Multprec_Complex_Vectors.Vector(p'range);
    jm : Multprec_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    rco : Floating_Number;
    loss_jac,loss_eva : integer32;
    ans : character;

  begin
    new_line;
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,sz,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Equal(rco,0.0) then
      loss_eva := -2**30;
    else
      declare
        mp_log10rco : Floating_Number := log10(rco);
        st_log10rco : double_float := Round(mp_log10rco);
      begin
        loss_eva := integer32(st_log10rco);
      end;
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    new_line;
    put("Check the values of the Jacobian matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Multprec_Check_Values(p,x,jm);
    end if;
  end Multprec_Jacobian_Test;

  procedure Test_Random_Conditioned_Jacobian_Evaluation is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of the numerical problem
  --   to set the condition of the Jacobian matrix and
  --   the condition of the polynomial evaluation problem.

    n,d,m,deci,size : natural32 := 0;
    condjm,coeffsize,pointsize,closeroot,cond : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put_line("First part, condition number of Jacobian matrix :");
    put("Give the condition of the Jacobian matrix : "); get(condjm);
    new_line;
    put_line("Second part, parameters of the random polynomials :");
    put("Give number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    new_line;
    put_line("Third part, factors in numerical condition of evaluation :");
    put("Give magnitude of the coefficients : "); get(coeffsize);
    put("Give magnitude of the coordinates of the point : "); get(pointsize);
    put("Give closeness to a root : "); get(closeroot);
    cond := coeffsize*(pointsize**integer(d))/closeroot;
    put("Predicted condition number : "); put(cond,3); new_line;
    new_line;
    put_line("Choose the level of precision : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => 
        Standard_Jacobian_Test(n,d,m,0,coeffsize,pointsize,closeroot,condjm);
      when '1' =>
        DoblDobl_Jacobian_Test(n,d,m,0,coeffsize,pointsize,closeroot,condjm);
      when '2' =>
        QuadDobl_Jacobian_Test(n,d,m,0,coeffsize,pointsize,closeroot,condjm);
      when '3' =>
        declare
          deci,size : natural32 := 0;
        begin
          new_line;
          put("Give the number of decimal places : "); get(deci);
          size := Decimal_to_Size(deci);
          Multprec_Jacobian_Test
            (n,d,m,0,size,coeffsize,pointsize,closeroot,condjm);
        end;
      when others => null;
    end case;
  end Test_Random_Conditioned_Jacobian_Evaluation;

  procedure Main is

  -- DESCRIPTION :
  --   Offers user choice between random polynomials or
  --   to give a Laurent polynomial.

    ans : character;

  begin
    new_line;
    put_line("Condition numbers for polynomial evaluation ...");
    new_line;
    put_line("MENU to test condition number computation :");
    put_line("  1. generate a random polynomial and select precision;");
    put_line("  2. test specific Laurent polynomial of Stefano Taschini;");
    put_line("  3. conditioned Jacobian matrix and evaluation problem.");
    put("Type 1, 2, or 3 to select the type of test : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Test_Random_Polynomial;
      when '2' => Test_Taschini_Laurential;
      when '3' => Test_Random_Conditioned_Jacobian_Evaluation;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_vmpeval;
