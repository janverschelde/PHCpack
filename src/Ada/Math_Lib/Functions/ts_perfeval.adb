with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Matrices;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Random_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Random_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Random_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lexicographical_Supports;
with Coefficient_Supported_Polynomials;  use Coefficient_Supported_Polynomials;
with Standard_Polynomial_Flatteners;
with DoblDobl_Polynomial_Flatteners;
with QuadDobl_Polynomial_Flatteners;
with Standard_Gradient_Evaluations;
with Standard_Jacobian_Evaluations;
with DoblDobl_Jacobian_Evaluations;
with QuadDobl_Jacobian_Evaluations;

procedure ts_perfeval is

  procedure Prompt_for_Random_Parameters ( n,d,m : out natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension n, largest degree d,
  --   and number of monomials m (0 for dense) as the parameters
  --   to generate a random polynomial system.
  --   The symbol table is initialized for n variables.

  begin
    new_line;
    put_line("Generating a random polynomial system ...");
    n := 0; put("Give the dimension : "); get(n);
    Symbol_Table.Init(n);
    d := 0; put("Give the maximum degree : "); get(d);
    m := 0; put("Give the number of monomials (0 for dense) : "); get(m);
  end Prompt_for_Random_Parameters;

  procedure Standard_Random_System
              ( p : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension, largest degree, number of
  --   monomials in every polynomial to generate a random system.

    n,d,m : natural32;

  begin
    Prompt_for_Random_Parameters(n,d,m);
    declare
      q : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    begin
      for i in 1..integer32(n) loop
        if m = 0
         then q(i) := Standard_Random_Polynomials.Random_Dense_Poly(n,d,0);
         else q(i) := Standard_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
        end if;
      end loop;
      put_line("The polynomial system : "); put_line(q);
      p := new Standard_Complex_Poly_Systems.Poly_Sys'(q);
    end;
  end Standard_Random_System;

  procedure DoblDobl_Random_System
              ( p : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension, largest degree, number of
  --   monomials in every polynomial to generate a random system.

    n,d,m : natural32;

  begin
    Prompt_for_Random_Parameters(n,d,m);
    declare
      q : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    begin
      for i in 1..integer32(n) loop
        if m = 0
         then q(i) := DoblDobl_Random_Polynomials.Random_Dense_Poly(n,d,0);
         else q(i) := DoblDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
        end if;
      end loop;
      put_line("The polynomial system : "); put_line(q);
      p := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(q);
    end;
  end DoblDobl_Random_System;

  procedure QuadDobl_Random_System
              ( p : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for dimension, largest degree, number of
  --   monomials in every polynomial to generate a random system.

    n,d,m : natural32;

  begin
    Prompt_for_Random_Parameters(n,d,m);
    declare
      q : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    begin
      for i in 1..integer32(n) loop
        if m = 0
         then q(i) := QuadDobl_Random_Polynomials.Random_Dense_Poly(n,d,0);
         else q(i) := QuadDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
        end if;
      end loop;
      put_line("The polynomial system : "); put_line(q);
      p := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(q);
    end;
  end QuadDobl_Random_System;

  procedure Standard_Evaluation
              ( m : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Creates the Jacobian matrix for p and evaluates it and p
  --   at m randomly generated complex vectors.

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;

    timer : Timing_Widget;
    n : constant integer32
      := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    x : Standard_Complex_Vectors.Vector(1..n);
    y : Standard_Complex_Vectors.Vector(p'range);
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    A : Standard_Complex_Matrices.Matrix(jf'range(1),jf'range(2));

  begin
    new_line;
    put_line("evaluating system and Jacobian matrix separately ...");
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      y := Eval(f,x);
      A := Eval(jf,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"double polynomial evaluation");
    Clear(f); Clear(jm); Clear(jf);
  end Standard_Evaluation;

  procedure DoblDobl_Evaluation
              ( m : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Creates the Jacobian matrix for p and evaluates it and p
  --   at m randomly generated complex vectors.

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;

    timer : Timing_Widget;
    n : constant integer32
      := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    x : DoblDobl_Complex_Vectors.Vector(1..n);
    y : DoblDobl_Complex_Vectors.Vector(p'range);
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    A : DoblDobl_Complex_Matrices.Matrix(jf'range(1),jf'range(2));

  begin
    new_line;
    put_line("evaluating system and Jacobian matrix separately ...");
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      y := Eval(f,x);
      A := Eval(jf,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"double double polynomial evaluation");
    Clear(f); Clear(jm); Clear(jf);
  end DoblDobl_Evaluation;

  procedure QuadDobl_Evaluation
              ( m : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Creates the Jacobian matrix for p and evaluates it and p
  --   at m randomly generated complex vectors.

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    timer : Timing_Widget;
    n : constant integer32
      := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    x : QuadDobl_Complex_Vectors.Vector(1..n);
    y : QuadDobl_Complex_Vectors.Vector(p'range);
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    A : QuadDobl_Complex_Matrices.Matrix(jf'range(1),jf'range(2));

  begin
    new_line;
    put_line("evaluating system and Jacobian matrix separately ...");
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      y := Eval(f,x);
      A := Eval(jf,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"quad double polynomial evaluation");
    Clear(f); Clear(jm); Clear(jf);
  end QuadDobl_Evaluation;

  function Add_Derivatives
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system p, added with all partial derivates
  --   ordered row after row from the Jacobian matrix.

    use Standard_Complex_Jaco_Matrices;

    n : constant integer32
      := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    d : constant integer32 := p'last + n*n;
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..d);
    ind : integer32 := 0;

  begin
    for i in p'range loop
      res(i) := p(i);
    end loop;
    ind := p'last;
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        ind := ind + 1;
        res(ind) := jm(i,j);
      end loop;
    end loop;
    return res;
  end Add_Derivatives;

  function Add_Derivatives
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system p, added with all partial derivates
  --   ordered row after row from the Jacobian matrix.

    use DoblDobl_Complex_Jaco_Matrices;

    n : constant integer32
      := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    d : constant integer32 := p'last + n*n;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..d);
    ind : integer32 := 0;

  begin
    for i in p'range loop
      res(i) := p(i);
    end loop;
    ind := p'last;
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        ind := ind + 1;
        res(ind) := jm(i,j);
      end loop;
    end loop;
    return res;
  end Add_Derivatives;

  function Add_Derivatives
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the system p, added with all partial derivates
  --   ordered row after row from the Jacobian matrix.

    use QuadDobl_Complex_Jaco_Matrices;

    n : constant integer32
      := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    d : constant integer32 := p'last + n*n;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..d);
    ind : integer32 := 0;

  begin
    for i in p'range loop
      res(i) := p(i);
    end loop;
    ind := p'last;
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        ind := ind + 1;
        res(ind) := jm(i,j);
      end loop;
    end loop;
    return res;
  end Add_Derivatives;

  procedure Standard_Flattened_Evaluation
              ( m : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Constructs a flattened representation for p and its derivatives
  --   and evaluates at m randomly generated vectors.

    use Standard_Polynomial_Flatteners;

    timer : Timing_Widget;
    n : constant integer32
      := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    x : Standard_Complex_Vectors.Vector(1..n);
    f : constant Standard_Complex_Poly_Systems.Poly_Sys := Add_Derivatives(p);
    y : Standard_Complex_Vectors.Vector(f'range);
    sp : List := Distinct_Supports(f);
    ssp : List := Lexicographical_Supports.Sort(sp);
    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(ssp);
    fv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Factor(ssp);
    cfv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Compress(fv);
    c : Standard_Complex_VecVecs.VecVec(f'range);
    k : Standard_Natural_VecVecs.VecVec(f'range);

  begin
    new_line;
    put("# different monomials in support : "); put(v'last,1); new_line;
    Coefficients_of_Supports(f,v,c,k);
    new_line;
    put_line("running in unfactored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      y := Eval(c,v,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened evaluation unfactored");
    new_line;
    put_line("running in factored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      y := Factored_Eval(c,fv,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened evaluation factored");
    new_line;
    put_line("running in compressed factored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      y := Factored_Compressed_Eval(c,cfv,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened compressed factored");
    Lists_of_Integer_Vectors.Clear(sp);
    Lists_of_Integer_Vectors.Clear(ssp);
  end Standard_Flattened_Evaluation;

  procedure DoblDobl_Flattened_Evaluation
              ( m : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Constructs a flattened representation for p and its derivatives
  --   and evaluates at m randomly generated vectors.

    use DoblDobl_Polynomial_Flatteners;

    timer : Timing_Widget;
    n : constant integer32
      := integer32(DoblDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    x : DoblDobl_Complex_Vectors.Vector(1..n);
    f : constant DoblDobl_Complex_Poly_Systems.Poly_Sys := Add_Derivatives(p);
    y : DoblDobl_Complex_Vectors.Vector(f'range);
    sp : List := Distinct_Supports(f);
    ssp : List := Lexicographical_Supports.Sort(sp);
    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(ssp);
    fv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Factor(ssp);
    cfv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Compress(fv);
    c : DoblDobl_Complex_VecVecs.VecVec(f'range);
    k : Standard_Natural_VecVecs.VecVec(f'range);

  begin
    new_line;
    put("# different monomials in support : "); put(v'last,1); new_line;
    Coefficients_of_Supports(f,v,c,k);
    new_line;
    put_line("running in unfactored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      y := Eval(c,v,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened evaluation unfactored");
    new_line;
    put_line("running in factored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      y := Factored_Eval(c,fv,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened evaluation factored");
    new_line;
    put_line("running in compressed factored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      y := Factored_Compressed_Eval(c,cfv,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened compressed factored");
    Lists_of_Integer_Vectors.Clear(sp);
    Lists_of_Integer_Vectors.Clear(ssp);
  end DoblDobl_Flattened_Evaluation;

  procedure QuadDobl_Flattened_Evaluation
              ( m : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Constructs a flattened representation for p and its derivatives
  --   and evaluates at m randomly generated vectors.

    use QuadDobl_Polynomial_Flatteners;

    timer : Timing_Widget;
    n : constant integer32
      := integer32(QuadDobl_Complex_Polynomials.Number_of_Unknowns(p(p'first)));
    x : QuadDobl_Complex_Vectors.Vector(1..n);
    f : constant QuadDobl_Complex_Poly_Systems.Poly_Sys := Add_Derivatives(p);
    y : QuadDobl_Complex_Vectors.Vector(f'range);
    sp : List := Distinct_Supports(f);
    ssp : List := Lexicographical_Supports.Sort(sp);
    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(ssp);
    fv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Factor(ssp);
    cfv : constant Standard_Integer_VecVecs.VecVec
       := Lexicographical_Supports.Compress(fv);
    c : QuadDobl_Complex_VecVecs.VecVec(f'range);
    k : Standard_Natural_VecVecs.VecVec(f'range);

  begin
    new_line;
    put("# different monomials in support : "); put(v'last,1); new_line;
    Coefficients_of_Supports(f,v,c,k);
    new_line;
    put_line("running in unfactored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      y := Eval(c,v,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened evaluation unfactored");
    new_line;
    put_line("running in factored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      y := Factored_Eval(c,fv,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened evaluation factored");
    new_line;
    put_line("running in compressed factored representation ...");
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      y := Factored_Compressed_Eval(c,cfv,k,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"flattened compressed factored");
    Lists_of_Integer_Vectors.Clear(sp);
    Lists_of_Integer_Vectors.Clear(ssp);
  end QuadDobl_Flattened_Evaluation;

  procedure Run_Standard_Jacobian_Evaluation
              ( m : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                e : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Generates m different random points and evaluates the system
  --   and its Jacobian matrix at those m points.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(e);
    v_n : constant Standard_Natural_VecVecs.VecVec
        := Standard_Jacobian_Evaluations.Integer_to_Natural(v);
    f,b : Standard_Natural_VecVecs.VecVec(v'range);
    c : Standard_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    n : constant integer32 := v(v'first)'last;
    x : Standard_Complex_Vectors.Vector(1..n);
    z : Standard_Complex_Vectors.Vector(p'range);
    A : Standard_Complex_Matrices.Matrix(p'range,1..n);
    wrk : Standard_Complex_VecVecs.VecVec(1..v'last);
    timer : Timing_Widget;

    use Standard_Jacobian_Evaluations;

  begin
    Standard_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    Split_Common_Factors(v_n,f,b);
    for k in wrk'range loop
      wrk(k) := new Standard_Complex_Vectors.Vector(0..n);
    end loop;
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      EvalDiff(f,b,c,k,x,z,wrk,A);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"evaluation & differentiation");
  end Run_Standard_Jacobian_Evaluation;

  procedure Run_DoblDobl_Jacobian_Evaluation
              ( m : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                e : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Generates m different random points and evaluates the system
  --   and its Jacobian matrix at those m points.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(e);
    c : DoblDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    n : constant integer32 := v(v'first)'last;
    x : DoblDobl_Complex_Vectors.Vector(1..n);
    z : DoblDobl_Complex_Vectors.Vector(p'range);
    A : DoblDobl_Complex_Matrices.Matrix(p'range,1..n);
    timer : Timing_Widget;

    use DoblDobl_Jacobian_Evaluations;

  begin
    DoblDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      DoblDobl_Jacobian_Evaluation(v,c,k,x,z,A);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"evaluation & differentiation");
  end Run_DoblDobl_Jacobian_Evaluation;

  procedure Run_QuadDobl_Jacobian_Evaluation
              ( m : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                e : in Lists_of_Integer_Vectors.List ) is

  -- DESCRIPTION :
  --   Generates m different random points and evaluates the system
  --   and its Jacobian matrix at those m points.

    v : constant Standard_Integer_VecVecs.VecVec
      := Lists_of_Integer_Vectors.Shallow_Create(e);
    c : QuadDobl_Complex_VecVecs.VecVec(p'range);
    k : Standard_Natural_VecVecs.VecVec(p'range);
    n : constant integer32 := v(v'first)'last;
    x : QuadDobl_Complex_Vectors.Vector(1..n);
    z : QuadDobl_Complex_Vectors.Vector(p'range);
    A : QuadDobl_Complex_Matrices.Matrix(p'range,1..n);
    timer : Timing_Widget;

    use QuadDobl_Jacobian_Evaluations;

  begin
    QuadDobl_Polynomial_Flatteners.Coefficients_of_Supports(p,v,c,k);
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      QuadDobl_Jacobian_Evaluation(v,c,k,x,z,A);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"evaluation & differentiation");
  end Run_QuadDobl_Jacobian_Evaluation;

  procedure Standard_Jacobian_Evaluation
              ( m : in natural32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Applies the algorithms to evaluate and differentiate the
  --   polynomial system p as many as m times for timings.

    s,e : Lists_of_Integer_Vectors.List;

  begin
    s := Standard_Polynomial_Flatteners.Distinct_Supports(p);
    e := Lexicographical_Supports.Sort(s);
    Run_Standard_Jacobian_Evaluation(m,p,e);
  end Standard_Jacobian_Evaluation;

  procedure DoblDobl_Jacobian_Evaluation
              ( m : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Applies the algorithms to evaluate and differentiate the
  --   polynomial system p as many as m times for timings.

    s,e : Lists_of_Integer_Vectors.List;

  begin
    s := DoblDobl_Polynomial_Flatteners.Distinct_Supports(p);
    e := Lexicographical_Supports.Sort(s);
    Run_DoblDobl_Jacobian_Evaluation(m,p,e);
  end DoblDobl_Jacobian_Evaluation;

  procedure QuadDobl_Jacobian_Evaluation
              ( m : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Applies the algorithms to evaluate and differentiate the
  --   polynomial system p as many as m times for timings.

    s,e : Lists_of_Integer_Vectors.List;

  begin
    s := QuadDobl_Polynomial_Flatteners.Distinct_Supports(p);
    e := Lexicographical_Supports.Sort(s);
    Run_QuadDobl_Jacobian_Evaluation(m,p,e);
  end QuadDobl_Jacobian_Evaluation;

  procedure Standard_Performance_Test is

  -- DESCRIPTION :
  --   Runs a performance test on polynomial evaluation
  --   in standard double precision arithmetic.

    use Standard_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    m : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Generate a random system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      Standard_Random_System(lp);
    else
      new_line;
      put_line("reading a polynomial system ...");
      get(lp);
      new_line;
      put("Read "); put(lp'last,1); put_line(" polynomials ...");
    end if;
    new_line;
    put("Give number of evaluations : "); get(m);
    Standard_Evaluation(m,lp.all);
    Standard_Flattened_Evaluation(m,lp.all);
    Standard_Jacobian_Evaluation(m,lp.all);
  end Standard_Performance_Test;

  procedure DoblDobl_Performance_Test is

  -- DESCRIPTION :
  --   Runs a performance test on polynomial evaluation
  --   in double double precision arithmetic.

    use DoblDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    m : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Generate a random system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      DoblDobl_Random_System(lp);
    else
      new_line;
      put_line("reading a polynomial system ...");
      get(lp);
      new_line;
      put("Read "); put(lp'last,1); put_line(" polynomials ...");
    end if;
    new_line;
    put("Give number of evaluations : "); get(m);
    DoblDobl_Evaluation(m,lp.all);
    DoblDobl_Flattened_Evaluation(m,lp.all);
    DoblDobl_Jacobian_Evaluation(m,lp.all);
  end DoblDobl_Performance_Test;

  procedure QuadDobl_Performance_Test is

  -- DESCRIPTION :
  --   Runs a performance test on polynomial evaluation
  --   in quad double precision arithmetic.

    use QuadDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    m : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Generate a random system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      QuadDobl_Random_System(lp);
    else
      new_line;
      put_line("reading a polynomial system ...");
      get(lp);
      new_line;
      put("Read "); put(lp'last,1); put_line(" polynomials ...");
    end if;
    new_line;
    put("Give number of evaluations : "); get(m);
    QuadDobl_Evaluation(m,lp.all);
    QuadDobl_Flattened_Evaluation(m,lp.all);
    QuadDobl_Jacobian_Evaluation(m,lp.all);
  end QuadDobl_Performance_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test performance of polynomial and Jacobian evaluation:");
    put_line("  1. standard doubles to evaluate polynomials and derivatives;");
    put_line("  2. double doubles to evaluate polynomials and derivatives;");
    put_line("  3. quad doubles to evaluate polynomials and derivatives.");
    put("Type 1, 2, or 3 to make your choice : ");
    Ask_Alternative(ans,"123");
    case ans is
       when '1' => Standard_Performance_Test;
       when '2' => DoblDobl_Performance_Test;
       when '3' => QuadDobl_Performance_Test;
       when others => null;
    end case;
  end Main;

begin
  Main;
end ts_perfeval;
