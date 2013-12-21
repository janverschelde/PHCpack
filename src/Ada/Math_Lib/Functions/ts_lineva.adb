with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;

procedure ts_lineva is

-- DESCRIPTION :
--   Benchmark testing on the evaluation of a linear form.

  function Poly_Linear_Form ( form : Vector ) return Poly is

  -- DESCRIPTION :
  --   Returns the representation of the form as a polynomial.
  --   The input parameter form should be a vector of range 0..n,
  --   where n = form'last.

    res : Poly := Null_Poly;
    n : constant integer32 := form'last;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    t.cf := form(0);
    Add(res,t);
    for i in 1..n loop
      t.dg(i) := 1;
      t.cf := form(i);
      Add(res,t);
      t.dg(i) := 0;
    end loop;
    return res;
  end Poly_Linear_Form;

  procedure Timer_Test ( form : in Vector; pfor : in Poly;
                         pfor_eval : in Eval_Poly ) is

    n : constant integer32 := form'last;
    nb : integer32 := 0;
    timer : Timing_Widget;
    x : Vector(1..n);
    y : Complex_Number;

  begin
    new_line;
    put_line("Evaluating the form many many times.");
    new_line;
    put("Give the number of times : "); get(nb);
    tstart(timer);
    for i in 1..nb loop
      x := Random_Vector(1,n);
      y := form(0) + form(1..n)*x;
    end loop;
    tstop(timer);
    print_times(Standard_Output,timer,"evaluating form as vector");
    tstart(timer);
    for i in 1..nb loop
      x := Random_Vector(1,n);
      y := Eval(pfor,x);
    end loop;
    tstop(timer);
    print_times(Standard_Output,timer,"evaluating form as polynomial");
    tstart(timer);
    for i in 1..nb loop
      x := Random_Vector(1,n);
      y := Eval(pfor_eval,x);
    end loop;
    tstop(timer);
    print_times(Standard_Output,timer,"evaluating form as evaluation poly");
  end Timer_Test;

  procedure Evaluate_Linear_Form ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random linear form in n variables and evaluates it.

    form : constant Vector(0..n) := Random_Vector(0,n);
    pfor : constant Poly := Poly_Linear_Form(form);
    pfor_eval : constant Eval_Poly := Create(pfor);
    x : Vector(1..n);
    y : Complex_Number;

  begin
    put("The random linear form : "); put_line(pfor);
    x := Random_Vector(1,n);
    y := form(0) + form(1..n)*x;
    put_line("evaluation of the form at random point");
    put("  value : "); put(y); new_line;
    y := Eval(pfor,x);
    put_line("computed again with polynomial representation");
    put("  value : "); put(y); new_line;
    y := Eval(pfor_eval,x);
    put_line("computed again with evaluation polynomial");
    put("  value : "); put(y); new_line;
    Timer_Test(form,pfor,pfor_eval);
  end Evaluate_Linear_Form;

  procedure Main is

    n : integer32 := 0;

  begin
    new_line;
    put_line("Benchmark testing on the evaluation of a linear form.");
    new_line;
    put("Give the dimension : "); get(n);
    Evaluate_Linear_Form(n);
  end Main;

begin
  Main;
end ts_lineva;
