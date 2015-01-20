with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;      use Standard_Floating_VecVecs_io;
with Standard_Floating_Matrices;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;
with Double_Double_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;       use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Matrices;
with Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;       use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Matrices;
with Standard_Random_VecVecs;
with DoblDobl_Random_VecVecs;
with QuadDobl_Random_VecVecs;
with Standard_Floating_GramSchmidt;
with Standard_Complex_GramSchmidt;
with DoblDobl_Complex_GramSchmidt;
with QuadDobl_Complex_GramSchmidt;
with Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;
with Test_LU_Decompositions;            use Test_LU_Decompositions;

procedure ts_mgs is

-- DESCRIPTION :
--   Development of modified Gram-Schmidt orthonormalization.

  procedure Show_Performance
              ( m,n,f : in integer32; timer : in Timing_Widget ) is

  -- DESCRIPTION :
  --   The number of flops of the modified Gram-Schmidt method is 2*m*n^2.
  --   Given m, n, and the timer, this procedure shows the number of
  --   flops per second.

    package duration_io is new text_io.fixed_io(duration);

    flops : constant integer32 := 2*m*n*n*f;
    user_cpu_time : constant duration := Elapsed_User_Time(timer);
    seconds : constant double_float := double_float(user_cpu_time);
    performance : constant integer32 := integer32(double_float(flops)/seconds);

  begin
    new_line;
    put("number of operations : "); put(flops,1); new_line;
    put("user CPU time : "); duration_io.put(user_cpu_time); new_line;
    put("#flops per second : "); put(performance,1); new_line;
  end Show_Performance;

  procedure Show_Accuracy_Report ( orterr,eqserr : in double_float ) is

  -- DESCRIPTION :
  --   Shows the largest deviation in orthonormality and decomposition,
  --   given respectively in orterr and eqserr.

  begin
    put_line("Accuracy level report :");
    put("  largest error in orthonormality :"); put(orterr,3); new_line;
    put("  largest error in decomposition  :"); put(eqserr,3); new_line;
  end Show_Accuracy_Report;

  procedure Standard_Floating_Random_Test ( n,m,g : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method.
 
    use Standard_Floating_GramSchmidt;

    v : constant Standard_Floating_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(m),natural32(g));
    q,r : Standard_Floating_VecVecs.VecVec(1..m);
    tol : constant double_float := 1.0E-8;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new Standard_Floating_Vectors.Vector'(v(i).all);
      r(i) := new Standard_Floating_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := 0.0;
      end loop;
    end loop;
    QR(n,m,q,r);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end Standard_Floating_Random_Test;

  procedure Standard_Complex_Random_Test ( n,m,g : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method.
 
    use Standard_Complex_GramSchmidt;

    v : constant Standard_Complex_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(m),natural32(g));
    q,r : Standard_Complex_VecVecs.VecVec(1..m);
    tol : constant double_float := 1.0E-8;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new Standard_Complex_Vectors.Vector'(v(i).all);
      r(i) := new Standard_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(0.0);
      end loop;
    end loop;
    QR(n,m,q,r);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end Standard_Complex_Random_Test;

  procedure DoblDobl_Complex_Random_Test ( n,m,g : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method.
 
    use DoblDobl_Complex_GramSchmidt;

    v : constant DoblDobl_Complex_VecVecs.VecVec(1..m)
      := DoblDobl_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(m),natural32(g));
    q,r : DoblDobl_Complex_VecVecs.VecVec(1..m);
    tol : constant double_float := 1.0E-16;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new DoblDobl_Complex_Vectors.Vector'(v(i).all);
      r(i) := new DoblDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    QR(n,m,q,r);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end DoblDobl_Complex_Random_Test;

  procedure QuadDobl_Complex_Random_Test ( n,m,g : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method.
 
    use QuadDobl_Complex_GramSchmidt;

    v : constant QuadDobl_Complex_VecVecs.VecVec(1..m)
      := QuadDobl_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(m),natural32(g));
    q,r : QuadDobl_Complex_VecVecs.VecVec(1..m);
    tol : constant double_float := 1.0E-32;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new QuadDobl_Complex_Vectors.Vector'(v(i).all);
      r(i) := new QuadDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    QR(n,m,q,r);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end QuadDobl_Complex_Random_Test;

  procedure Standard_Floating_Performance_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    use Standard_Floating_GramSchmidt;

    v : constant Standard_Floating_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : Standard_Floating_VecVecs.VecVec(1..m);
    A,B : Standard_Floating_Matrices.Matrix(1..n,1..m);
    qq : Standard_Floating_Matrices.Matrix(1..n,1..m);
    rr : Standard_Floating_Matrices.Matrix(1..m,1..m);
    f : integer32 := 0;
    timer : Timing_Widget;

  begin
    for i in 1..m loop
      q(i) := new Standard_Floating_Vectors.Vector(1..n);
      r(i) := new Standard_Floating_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := 0.0;
      end loop;
    end loop;
    put("give number of QR decompositions : "); get(f);
    tstart(timer);
    for k in 1..f loop
      for j in 1..m loop
        for i in 1..n loop
          q(j)(i) := v(j)(i);
        end loop;
      end loop;
      QR(n,m,q,r);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"QR decomposition on vectors");
    Show_Performance(m,n,f,timer);
  --  for j in 1..m loop
  --    for i in 1..n loop
  --      A(i,j) := v(j)(i);
  --    end loop; 
  --  end loop;
  --  tstart(timer);
  --  for k in 1..f loop
  --    for i in 1..n loop
  --      for j in 1..m loop
  --        B(i,j) := A(i,j);
  --      end loop;
  --    end loop;
  --    QR(n,m,B,qq,rr);
  --  end loop;
  --  tstop(timer);
  --  new_line;
  --  print_times(standard_output,timer,"QR decomposition on matrices");
  end Standard_Floating_Performance_Test;

  procedure Standard_Complex_Performance_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    use Standard_Complex_GramSchmidt;

    v : constant Standard_Complex_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : Standard_Complex_VecVecs.VecVec(1..m);
    f : integer32 := 0;
    timer : Timing_Widget;

  begin
    for i in 1..m loop
      q(i) := new Standard_Complex_Vectors.Vector(1..n);
      r(i) := new Standard_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(0.0);
      end loop;
    end loop;
    put("give number of QR decompositions : "); get(f);
    tstart(timer);
    for k in 1..f loop
      for j in 1..m loop
        for i in 1..n loop
          q(j)(i) := v(j)(i);
        end loop;
      end loop;
      QR(n,m,q,r);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"QR decomposition");
  end Standard_Complex_Performance_Test;

  procedure DoblDobl_Complex_Performance_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    use DoblDobl_Complex_GramSchmidt;

    v : constant DoblDobl_Complex_VecVecs.VecVec(1..m)
      := DoblDobl_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : DoblDobl_Complex_VecVecs.VecVec(1..m);
    f : integer32 := 0;
    timer : Timing_Widget;

  begin
    for i in 1..m loop
      q(i) := new DoblDobl_Complex_Vectors.Vector(1..n);
      r(i) := new DoblDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    put("give number of QR decompositions : "); get(f);
    tstart(timer);
    for k in 1..f loop
      for j in 1..m loop
        for i in 1..n loop
          q(j)(i) := v(j)(i);
        end loop;
      end loop;
      QR(n,m,q,r);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"QR decomposition");
  end DoblDobl_Complex_Performance_Test;

  procedure QuadDobl_Complex_Performance_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    use QuadDobl_Complex_GramSchmidt;

    v : constant QuadDobl_Complex_VecVecs.VecVec(1..m)
      := QuadDobl_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : QuadDobl_Complex_VecVecs.VecVec(1..m);
    f : integer32 := 0;
    timer : Timing_Widget;

  begin
    for i in 1..m loop
      q(i) := new QuadDobl_Complex_Vectors.Vector(1..n);
      r(i) := new QuadDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    put("give number of QR decompositions : "); get(f);
    tstart(timer);
    for k in 1..f loop
      for j in 1..m loop
        for i in 1..n loop
          q(j)(i) := v(j)(i);
        end loop;
      end loop;
      QR(n,m,q,r);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"QR decomposition");
  end QuadDobl_Complex_Performance_Test;

  function Minimum ( e : in Standard_Floating_Vectors.Vector )
                   return double_float is

  -- DESCRIPTION :
  --   Returns the minimum value of the elements in e.

    res : double_float := e(e'first);

  begin
    for i in e'first+1..e'last loop
      if e(i) < res
       then res := e(i);
      end if;
    end loop;
    return res;
  end Minimum;

  function Maximum ( e : in Standard_Floating_Vectors.Vector )
                   return double_float is

  -- DESCRIPTION :
  --   Returns the largest value of the elements in e.

    res : double_float := e(e'first);

  begin
    for i in e'first+1..e'last loop
      if e(i) > res
       then res := e(i);
      end if;
    end loop;
    return res;
  end Maximum;

  procedure Show_Statistics 
              ( g : in integer32;
                e : in Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shows the minima and maxima for the errors of number 
  --   of magnitude g.

    min : constant double_float := Minimum(e);
    log10min : constant double_float := LOG10(min);
    max : constant double_float := Maximum(e);
    log10max : constant double_float := LOG10(max);

  begin
    -- put("  min : "); put(min,3);
    put(g,2); put(" : ");
    put("  Log10(min) : "); put(LOG10min,3,1,0);
    -- put("  max : "); put(max,3);
    put("  Log10(max) : "); put(LOG10max,3,1,0);
    put("  D = "); put(LOG10max - LOG10min,3,1,0);
    new_line;
  end Show_Statistics;

  procedure Standard_Complex_Initialize_Data
               ( n,m,g : in integer32;
                 v,q,r : out Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Generates m vectors of length n and of magnitude g in v,
  --   copies v to q and initializes r.

  -- REQUIRED : v'range = 1..m = q'range, r'range = 1..n.

  begin
    v := Standard_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(m),natural32(g));
    for i in 1..m loop
      q(i) := new Standard_Complex_Vectors.Vector'(v(i).all);
      r(i) := new Standard_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
  end Standard_Complex_Initialize_Data;

  procedure Standard_Complex_Initialize_Data
               ( n,g : in integer32;
                 A : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Generates n vectors of length n and of magnitude g in v,
  --   copies v to A.

  -- REQUIRED : A'range(1) = A'range(2) = 1..n.

    v : Standard_Complex_VecVecs.VecVec(1..n);

  begin
    v := Standard_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(n),natural32(g));
    for j in 1..n loop
      for i in 1..n loop
         A(i,j) := v(j)(i);
      end loop;
    end loop;
  end Standard_Complex_Initialize_Data;

  procedure DoblDobl_Complex_Initialize_Data
               ( n,m,g : in integer32;
                 v,q,r : out DoblDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Generates m vectors of length n and of magnitude g in v,
  --   copies v to w and initializes q and r.

  -- REQUIRED : v'range = 1..m = w'range = q'range, r'range = 1..n.

  begin
    v := DoblDobl_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(m),natural32(g));
    for i in 1..m loop
      q(i) := new DoblDobl_Complex_Vectors.Vector'(v(i).all);
      r(i) := new DoblDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
  end DoblDobl_Complex_Initialize_Data;

  procedure DoblDobl_Complex_Initialize_Data
               ( n,g : in integer32;
                 A : out DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Generates n vectors of length n and of magnitude g in v,
  --   copies v to A.

  -- REQUIRED : A'range(1) = A'range(2) = 1..n.

   v : DoblDobl_Complex_VecVecs.VecVec(1..n);

  begin
    v := DoblDobl_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(n),natural32(g));
    for j in 1..n loop
      for i in 1..n loop
        A(i,j) := v(j)(i);
      end loop;
    end loop;
  end DoblDobl_Complex_Initialize_Data;

  procedure QuadDobl_Complex_Initialize_Data
               ( n,m,g : in integer32;
                 v,q,r : out QuadDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Generates m vectors of length n and of magnitude g in v,
  --   copies v to w and initializes q and r.

  -- REQUIRED : v'range = 1..m = q'range, r'range = 1..n.

  begin
    v := QuadDobl_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(m),natural32(g));
    for i in 1..m loop
      q(i) := new QuadDobl_Complex_Vectors.Vector'(v(i).all);
      r(i) := new QuadDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
  end QuadDobl_Complex_Initialize_Data;

  procedure QuadDobl_Complex_Initialize_Data
               ( n,g : in integer32;
                 A : out QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Generates n vectors of length n and of magnitude g in v,
  --   copies v to A.

  -- REQUIRED : A'range(1) = Arange(2) = 1..n.

    v : QuadDobl_Complex_VecVecs.VecVec(1..n);

  begin
    v := QuadDobl_Random_VecVecs.Random_VecVec
           (natural32(n),natural32(n),natural32(g));
    for j in 1..n loop
      for i in 1..n loop
        A(i,j) := v(j)(i);
      end loop;
    end loop;
  end QuadDobl_Complex_Initialize_Data;

  procedure Standard_Complex_Run_MGS_Accuracy_Test
              ( n,m,f,g_start,g_end : in integer32 ) is

  -- DESCRIPTION :
  --   For the range of magnitudes in g_start to g_end,
  --   f samples are taken and the largest error is computed.

    use Standard_Complex_GramSchmidt;

    v,q,r : Standard_Complex_VecVecs.VecVec(1..m);
    tol : constant double_float := 1.0E-32;
    err : Standard_Floating_Vectors.Vector(1..f);
    fail : boolean;

  begin
    for g in g_start..g_end loop
     -- put("Generating "); put(f,1); 
     -- put(" samples of magnitude "); put(g,1);
     -- put_line("...");
      for i in 1..f loop
        Standard_Complex_Initialize_Data(n,m,g,v,q,r);
        QR(n,m,q,r);
        Test_Decomposition(n,m,v,q,r,tol,false,err(i),fail);
       -- put("The largest error : "); put(err(i),3);
       -- put("  and its log : "); put(LOG10(err(i)),3);  new_line;
      end loop;
      Show_Statistics(g,err);
    end loop;
  end Standard_Complex_Run_MGS_Accuracy_Test;

  procedure Standard_Complex_Run_LU_Accuracy_Test
              ( n,f,g_start,g_end : in integer32 ) is

  -- DESCRIPTION :
  --   For the range of magnitudes in g_start to g_end,
  --   f samples are taken and the largest error is computed.

    use Standard_Complex_Linear_Solvers;

    A,LU : Standard_Complex_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    tol : constant double_float := 1.0E-32;
    err : Standard_Floating_Vectors.Vector(1..f);
    info : integer32;
    fail : boolean;

  begin
    for g in g_start..g_end loop
     -- put("Generating "); put(f,1); 
     -- put(" samples of magnitude "); put(g,1);
     -- put_line("...");
      for i in 1..f loop
        Standard_Complex_Initialize_Data(n,g,A);
        LU := A;
        lufac(LU,n,ipvt,info);
        Test_Decomposition(n,A,LU,ipvt,tol,false,err(i),fail);
       -- put("The largest error : "); put(err(i),3);
       -- put("  and its log : "); put(LOG10(err(i)),3);  new_line;
      end loop;
      Show_Statistics(g,err);
    end loop;
  end Standard_Complex_Run_LU_Accuracy_Test;

  procedure DoblDobl_Complex_Run_MGS_Accuracy_Test
              ( n,m,f,g_start,g_end : in integer32 ) is

  -- DESCRIPTION :
  --   For the range of magnitudes in g_start to g_end,
  --   f samples are taken and the largest error is computed.

    use DoblDobl_Complex_GramSchmidt;

    v,q,r : DoblDobl_Complex_VecVecs.VecVec(1..m);
    tol : constant double_float := 1.0E-32;
    err : Standard_Floating_Vectors.Vector(1..f);
    fail : boolean;

  begin
    for g in g_start..g_end loop
     -- put("Generating "); put(f,1); 
     -- put(" samples of magnitude "); put(g,1);
     -- put_line("...");
      for i in 1..f loop
        DoblDobl_Complex_Initialize_Data(n,m,g,v,q,r);
        QR(n,m,q,r);
        Test_Decomposition(n,m,v,q,r,tol,false,err(i),fail);
       -- put("The largest error : "); put(err(i),3);
       -- put("  and its log : "); put(LOG10(err(i)),3);  new_line;
      end loop;
      Show_Statistics(g,err);
    end loop;
  end DoblDobl_Complex_Run_MGS_Accuracy_Test;

  procedure DoblDobl_Complex_Run_LU_Accuracy_Test
              ( n,f,g_start,g_end : in integer32 ) is

  -- DESCRIPTION :
  --   For the range of magnitudes in g_start to g_end,
  --   f samples are taken and the largest error is computed.

    A,LU : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    tol : constant Double_Double_Numbers.double_double
        := Double_Double_Numbers.create(1.0E-32);
    err : Standard_Floating_Vectors.Vector(1..f);
    maxerr : Double_Double_Numbers.double_double;
    info : integer32;
    fail : boolean;

  begin
    for g in g_start..g_end loop
     -- put("Generating "); put(f,1); 
     -- put(" samples of magnitude "); put(g,1);
     -- put_line("...");
      for i in 1..f loop
        DoblDobl_Complex_Initialize_Data(n,g,A);
        LU := A;
        DoblDobl_Complex_Linear_Solvers.lufac(LU,n,ipvt,info);
        Test_Decomposition(n,A,LU,ipvt,tol,false,maxerr,fail);
        err(i) := Double_Double_Numbers.to_double(maxerr);
       -- put("The largest error : "); put(err(i),3);
       -- put("  and its log : "); put(LOG10(err(i)),3);  new_line;
      end loop;
      Show_Statistics(g,err);
    end loop;
  end DoblDobl_Complex_Run_LU_Accuracy_Test;

  procedure QuadDobl_Complex_Run_MGS_Accuracy_Test
              ( n,m,f,g_start,g_end : in integer32 ) is

  -- DESCRIPTION :
  --   For the range of magnitudes in g_start to g_end,
  --   f samples are taken and the largest error is computed.

    use QuadDobl_Complex_GramSchmidt;

    v,q,r : QuadDobl_Complex_VecVecs.VecVec(1..m);
    tol : constant double_float := 1.0E-32;
    err : Standard_Floating_Vectors.Vector(1..f);
    fail : boolean;

  begin
    for g in g_start..g_end loop
     -- put("Generating "); put(f,1); 
     -- put(" samples of magnitude "); put(g,1);
     -- put_line("...");
      for i in 1..f loop
        QuadDobl_Complex_Initialize_Data(n,m,g,v,q,r);
        QR(n,m,q,r);
        Test_Decomposition(n,m,v,q,r,tol,false,err(i),fail);
       -- put("The largest error : "); put(err(i),3);
       -- put("  and its log : "); put(LOG10(err(i)),3);  new_line;
      end loop;
      Show_Statistics(g,err);
    end loop;
  end QuadDobl_Complex_Run_MGS_Accuracy_Test;

  procedure QuadDobl_Complex_Run_LU_Accuracy_Test
              ( n,f,g_start,g_end : in integer32 ) is

  -- DESCRIPTION :
  --   For the range of magnitudes in g_start to g_end,
  --   f samples are taken and the largest error is computed.

    A,LU : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    tol : constant Quad_Double_Numbers.quad_double
        := Quad_Double_Numbers.create(1.0E-32);
    err : Standard_Floating_Vectors.Vector(1..f);
    maxerr : Quad_Double_Numbers.quad_double;
    info : integer32;
    fail : boolean;

  begin
    for g in g_start..g_end loop
     -- put("Generating "); put(f,1); 
     -- put(" samples of magnitude "); put(g,1);
     -- put_line("...");
      for i in 1..f loop
        QuadDobl_Complex_Initialize_Data(n,g,A);
        LU := A;
        QuadDobl_Complex_Linear_Solvers.lufac(LU,n,ipvt,info);
        Test_Decomposition(n,A,LU,ipvt,tol,false,maxerr,fail);
        err(i) := Quad_Double_Numbers.to_double(maxerr);
       -- put("The largest error : "); put(err(i),3);
       -- put("  and its log : "); put(LOG10(err(i)),3);  new_line;
      end loop;
      Show_Statistics(g,err);
    end loop;
  end QuadDobl_Complex_Run_LU_Accuracy_Test;

  procedure Standard_Complex_MGS_Accuracy_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    use Standard_Complex_GramSchmidt;

    f,g_start,g_end : integer32 := 0;

  begin
    put("Give start of magnitude range : "); get(g_start);
    put("Give end of magnitude range : "); get(g_end);
    put("Give number of tests per magnitude : "); get(f);
    Standard_Complex_Run_MGS_Accuracy_Test(n,m,f,g_start,g_end);
  end Standard_Complex_MGS_Accuracy_Test;

  procedure Standard_Complex_LU_Accuracy_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the LU decomposition to time the performance.

    f,g_start,g_end : integer32 := 0;

  begin
    put("Give start of magnitude range : "); get(g_start);
    put("Give end of magnitude range : "); get(g_end);
    put("Give number of tests per magnitude : "); get(f);
    Standard_Complex_Run_LU_Accuracy_Test(n,f,g_start,g_end);
  end Standard_Complex_LU_Accuracy_Test;

  procedure DoblDobl_Complex_MGS_Accuracy_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    use DoblDobl_Complex_GramSchmidt;

    f,g_start,g_end : integer32 := 0;

  begin
    put("Give start of magnitude range : "); get(g_start);
    put("Give end of magnitude range : "); get(g_end);
    put("Give number of tests per magnitude : "); get(f);
    DoblDobl_Complex_Run_MGS_Accuracy_Test(n,m,f,g_start,g_end);
  end DoblDobl_Complex_MGS_Accuracy_Test;

  procedure DoblDobl_Complex_LU_Accuracy_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the LU decomposition to time the performance.

    f,g_start,g_end : integer32 := 0;

  begin
    put("Give start of magnitude range : "); get(g_start);
    put("Give end of magnitude range : "); get(g_end);
    put("Give number of tests per magnitude : "); get(f);
    DoblDobl_Complex_Run_LU_Accuracy_Test(n,f,g_start,g_end);
  end DoblDobl_Complex_LU_Accuracy_Test;

  procedure QuadDobl_Complex_MGS_Accuracy_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    use QuadDobl_Complex_GramSchmidt;

    f,g_start,g_end : integer32 := 0;

  begin
    put("Give start of magnitude range : "); get(g_start);
    put("Give end of magnitude range : "); get(g_end);
    put("Give number of tests per magnitude : "); get(f);
    QuadDobl_Complex_Run_MGS_Accuracy_Test(n,m,f,g_start,g_end);
  end QuadDobl_Complex_MGS_Accuracy_Test;

  procedure QuadDobl_Complex_LU_Accuracy_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a frequency f and generates f many
  --   random vector configurations of m vectors of range 1..n
  --   and applies the QR decomposition to time the performance.

    f,g_start,g_end : integer32 := 0;

  begin
    put("Give start of magnitude range : "); get(g_start);
    put("Give end of magnitude range : "); get(g_end);
    put("Give number of tests per magnitude : "); get(f);
    QuadDobl_Complex_Run_LU_Accuracy_Test(n,f,g_start,g_end);
  end QuadDobl_Complex_LU_Accuracy_Test;

  procedure Standard_Floating_Test_Solver ( n,m : in integer32 ) is

    use Standard_Floating_GramSchmidt;

    v : constant Standard_Floating_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : Standard_Floating_VecVecs.VecVec(1..m);
    x : constant Standard_Floating_Vectors.Vector(1..m) := (1..m => 1.0);
    b : Standard_Floating_Vectors.Vector(1..n) := Matrix_Product(n,m,v,x);
    Qb,Qb2 : Standard_Floating_Vectors.Vector(1..m);
    sol : Standard_Floating_Vectors.Vector(1..m);
    d : double_float;

  begin
    for i in 1..m loop
      q(i) := new Standard_Floating_Vectors.Vector(1..n);
      r(i) := new Standard_Floating_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := 0.0;
      end loop;
    end loop;
    for j in 1..m loop
      for i in 1..n loop
        q(j)(i) := v(j)(i);
      end loop;
    end loop;
    QR(n,m,q,r);
    Qb := Matrix_Projection(n,m,q,b);
    put_line("the vector Q^T*b :"); put_line(Qb);
    Orthogonal_Projection(n,m,q,b,Qb2);
    put_line("the projected vector b : "); put_line(b);
    put_line("the vector Q^T*b :"); put_line(Qb2);
    d := 0.0;
    for i in 1..m loop
      d := d + abs(Qb(i) - Qb2(i));
    end loop;
    put("sum norm of difference :"); put(d); new_line;
    sol := Solve(m,r,Qb2);
    put_line("the solution : "); put_line(sol); 
    d := 0.0;
    for i in 1..m loop
      d := d + abs(x(i) - sol(i));
    end loop;
    put("sum norm of difference :"); put(d); new_line;
  end Standard_Floating_Test_Solver;

  procedure Standard_Complex_Test_Solver ( n,m : in integer32 ) is

    use Standard_Complex_GramSchmidt;

    v : constant Standard_Complex_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : Standard_Complex_VecVecs.VecVec(1..m);
    x : constant Standard_Complex_Vectors.Vector(1..m)
      := (1..m => Create(1.0));
    b : Standard_Complex_Vectors.Vector(1..n) := Matrix_Product(n,m,v,x);
    Qb,Qb2 : Standard_Complex_Vectors.Vector(1..m);
    sol : Standard_Complex_Vectors.Vector(1..m);
    d : double_float;

  begin
    for i in 1..m loop
      q(i) := new Standard_Complex_Vectors.Vector(1..n);
      r(i) := new Standard_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(0.0);
      end loop;
    end loop;
    for j in 1..m loop
      for i in 1..n loop
        q(j)(i) := v(j)(i);
      end loop;
    end loop;
    QR(n,m,q,r);
    Qb := Matrix_Projection(n,m,q,b);
    put_line("the vector Q^T*b :"); put_line(Qb);
    Orthogonal_Projection(n,m,q,b,Qb2);
    put_line("the projected vector b : "); put_line(b);
    put_line("the vector Q^T*b :"); put_line(Qb2);
    d := 0.0;
    for i in 1..m loop
      d := d + AbsVal(Qb(i) - Qb2(i));
    end loop;
    put("sum norm of difference :"); put(d); new_line;
    sol := Solve(m,r,Qb2);
    put_line("the solution : "); put_line(sol); 
    d := 0.0;
    for i in 1..m loop
      d := d + AbsVal(x(i) - sol(i));
    end loop;
    put("sum norm of difference :"); put(d); new_line;
  end Standard_Complex_Test_Solver;

  procedure DoblDobl_Complex_Test_Solver ( n,m : in integer32 ) is

    use Double_Double_Numbers;
    use DoblDobl_Complex_GramSchmidt;

    v : constant DoblDobl_Complex_VecVecs.VecVec(1..m)
      := DoblDobl_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : DoblDobl_Complex_VecVecs.VecVec(1..m);
    x : constant DoblDobl_Complex_Vectors.Vector(1..m)
      := (1..m => Create(integer(1)));
    b : DoblDobl_Complex_Vectors.Vector(1..n) := Matrix_Product(n,m,v,x);
    Qb,Qb2 : DoblDobl_Complex_Vectors.Vector(1..m);
    sol : DoblDobl_Complex_Vectors.Vector(1..m);
    d : double_double;

  begin
    for i in 1..m loop
      q(i) := new DoblDobl_Complex_Vectors.Vector(1..n);
      r(i) := new DoblDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    for j in 1..m loop
      for i in 1..n loop
        q(j)(i) := v(j)(i);
      end loop;
    end loop;
    QR(n,m,q,r);
    Qb := Matrix_Projection(n,m,q,b);
    put_line("the vector Q^T*b :"); put_line(Qb);
    Orthogonal_Projection(n,m,q,b,Qb2);
    put_line("the projected vector b : "); put_line(b);
    put_line("the vector Q^T*b :"); put_line(Qb2);
    d := create(0.0);
    for i in 1..m loop
      d := d + AbsVal(Qb(i) - Qb2(i));
    end loop;
    put("sum norm of difference :"); put(to_double(d)); new_line;
    sol := Solve(m,r,Qb2);
    put_line("the solution : "); put_line(sol); 
    d := create(0.0);
    for i in 1..m loop
      d := d + AbsVal(x(i) - sol(i));
    end loop;
    put("sum norm of difference :"); put(to_double(d)); new_line;
  end DoblDobl_Complex_Test_Solver;

  procedure QuadDobl_Complex_Test_Solver ( n,m : in integer32 ) is

    use Quad_Double_Numbers;
    use QuadDobl_Complex_GramSchmidt;

    v : constant QuadDobl_Complex_VecVecs.VecVec(1..m)
      := QuadDobl_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : QuadDobl_Complex_VecVecs.VecVec(1..m);
    x : constant QuadDobl_Complex_Vectors.Vector(1..m)
      := (1..m => Create(integer(1)));
    b : QuadDobl_Complex_Vectors.Vector(1..n) := Matrix_Product(n,m,v,x);
    Qb,Qb2 : QuadDobl_Complex_Vectors.Vector(1..m);
    sol : QuadDobl_Complex_Vectors.Vector(1..m);
    d : quad_double;

  begin
    for i in 1..m loop
      q(i) := new QuadDobl_Complex_Vectors.Vector(1..n);
      r(i) := new QuadDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    for j in 1..m loop
      for i in 1..n loop
        q(j)(i) := v(j)(i);
      end loop;
    end loop;
    QR(n,m,q,r);
    Qb := Matrix_Projection(n,m,q,b);
    put_line("the vector Q^T*b :"); put_line(Qb);
    Orthogonal_Projection(n,m,q,b,Qb2);
    put_line("the projected vector b : "); put_line(b);
    put_line("the vector Q^T*b :"); put_line(Qb2);
    d := create(0.0);
    for i in 1..m loop
      d := d + AbsVal(Qb(i) - Qb2(i));
    end loop;
    put("sum norm of difference :"); put(to_double(d)); new_line;
    sol := Solve(m,r,Qb2);
    put_line("the solution : "); put_line(sol); 
    d := create(0.0);
    for i in 1..m loop
      d := d + AbsVal(x(i) - sol(i));
    end loop;
    put("sum norm of difference :"); put(to_double(d)); new_line;
  end QuadDobl_Complex_Test_Solver;

  procedure Main is

    n,m,g : integer32 := 0;
    ans : character;
    type_test : integer32 := 0;

  begin
    new_line;
    put_line("MENU to test modified Gram-Schmidt orthonormalization :");
    put_line("  1. interactive test on random vectors in nice range;");
    put_line("  2. performance test on random vectors in nice range;");
    put_line("  3. test accuracy of MGS on random vectors in a wide range;");
    put_line("  4. test accuracy of LU on random vectors in a wide range;");
    put_line("  5. test the solver on a random linear system.");
    put("Type 1, 2, 3, 4, or 5 to choose : "); Ask_Alternative(ans,"12345");
    case ans is
       when '1' => type_test := 1;
       when '2' => type_test := 2;
       when '3' => type_test := 3;
       when '4' => type_test := 4;
       when others => type_test := 5;
    end case;
    new_line;
    put_line("MENU to choose the type of arithmetic :");
    put_line("  1. standard floating point numbers;");
    put_line("  2. standard complex floating point numbers;");
    put_line("  3. double double complex floating point numbers;");
    put_line("  4. quad double complex floating point numbers;");
    put("Type 1, 2, 3, or 4 to choose : "); Ask_Alternative(ans,"1234");
    new_line;
    put("Give dimension : "); get(n);
    if type_test /= 4
     then put("Give number of vectors : "); get(m);
    end if;
    new_line;
    case type_test is
      when 1 =>
        put_line("running one test on random vectors ...");
        case ans is
          when '1' => Standard_Floating_Random_Test(n,m,1);
          when '2' => Standard_Complex_Random_Test(n,m,1);
          when '3' => DoblDobl_Complex_Random_Test(n,m,1);
          when '4' => QuadDobl_Complex_Random_Test(n,m,1);
          when others => null;
        end case;
      when 2 =>
        put_line("running a performance test ...");
        case ans is 
          when '1' => Standard_Floating_Performance_Test(n,m);
          when '2' => Standard_Complex_Performance_Test(n,m);
          when '3' => DoblDobl_Complex_Performance_Test(n,m);
          when '4' => QuadDobl_Complex_Performance_Test(n,m);
          when others => null;
        end case;
      when 3 =>
        put_line("running an accuracy test on MGS ...");
        put("Give the magnitude (0 for range) : "); get(g);
        if g > 0 then
          case ans is
            when '1' => Standard_Floating_Random_Test(n,m,g);
            when '2' => Standard_Complex_Random_Test(n,m,g);
            when '3' => DoblDobl_Complex_Random_Test(n,m,g);
            when '4' => QuadDobl_Complex_Random_Test(n,m,g);
            when others => null;
          end case;
        else
          case ans is
            when '2' => Standard_Complex_MGS_Accuracy_Test(n,m);
            when '3' => DoblDobl_Complex_MGS_Accuracy_Test(n,m);
            when '4' => QuadDobl_Complex_MGS_Accuracy_Test(n,m);
            when others => null;
          end case;
        end if;
      when 4 =>
        put_line("running an accuracy test on LU ...");
        case ans is
          when '2' => Standard_Complex_LU_Accuracy_Test(n);
          when '3' => DoblDobl_Complex_LU_Accuracy_Test(n);
          when '4' => QuadDobl_Complex_LU_Accuracy_Test(n);
          when others => null;
        end case;
      when 5 =>
        put_line("Testing the solver...");
        case ans is
          when '1' => Standard_Floating_Test_Solver(n,m);
          when '2' => Standard_Complex_Test_Solver(n,m);
          when '3' => DoblDobl_Complex_Test_Solver(n,m);
          when '4' => QuadDobl_Complex_Test_Solver(n,m);
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mgs;
