with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Standard_Polynomial_Interpolators;
with Multprec_Polynomial_Interpolators;
with Sampling_Machine;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Rectangular_Sample_Grids;           use Rectangular_Sample_Grids;
with Standard_Stacked_Sample_Grids;
with Multprec_Stacked_Sample_Grids;
--with Drivers_to_Grid_Creators;           use Drivers_to_Grid_Creators;
with Make_Sample_Grids;                  use Make_Sample_Grids;
with Divided_Differences;
with Standard_Divided_Differences;
with QuadDobl_Divided_Differences;
with DoblDobl_Divided_Differences;
with Multprec_Divided_Differences;

procedure ts_divdif is

-- DESCRIPTION :
--   Interactive test facility for the multivariate interpolation
--   with divided differences.

-- EXPERIMENT : COMPARE WITH DIRECT APPROACH :

  procedure Standard_Experiment
               ( file : in file_type;
                 eq : in Standard_Complex_Polynomials.Poly;
                 grid : in Array_of_Standard_Sample_Lists ) is

  -- DESCRIPTION :
  --   Computes the interpolating polynomial by solving a linear system.

    d : constant integer32 := grid'last;
    samples : constant Standard_Complex_VecVecs.VecVec
            := Extract_Samples(grid);
    dense,ip : Standard_Complex_Polynomials.Poly;
    rcond : double_float;

    use Standard_Complex_Polynomials;
    use Standard_Polynomial_Interpolators;

  begin
    put_line(file,"Testing with direct interpolation approach.");
    dense := Create(natural32(d),2,1);
    put_line(file,"The test polynomial : "); put_line(file,dense);
    Interpolate(dense,samples,ip,rcond);
    put(file,"The estimate of the inverse condition number : ");
    put(file,rcond,3); new_line(file);
    put_line(file,"The interpolating polynomial : ");
    put_line(file,ip);
    put(file,"Distance expanded Newton form and exact : ");
    put(file,Distance(dense,eq),3); new_line(file);
    put(file,"Distance direct interpolator and exact  : ");
    put(file,Distance(dense,ip),3); new_line(file);
  end Standard_Experiment;

  procedure Multprec_Experiment
               ( file : in file_type;
                 eq : in Multprec_Complex_Polynomials.Poly;
                 grid : in Array_of_Multprec_Sample_Lists ) is

  -- DESCRIPTION :
  --   Computes the interpolating polynomial by solving a linear system.

    d : constant integer32 := grid'last;
    samples : constant Multprec_Complex_VecVecs.VecVec
            := Extract_Samples(grid);
    dense,ip : Multprec_Complex_Polynomials.Poly;
    rcond : Floating_Number;
    m : natural32;

    use Multprec_Complex_Polynomials;
    use Multprec_Polynomial_Interpolators;

  begin
    put_line(file,"Testing with direct interpolation approach.");
    dense := Create(natural32(d),2,1);
    m := Number_of_Terms(dense)-1;
    put_line(file,"The test polynomial : "); put_line(file,dense);
    Interpolate1(dense,samples(1..integer32(m)),ip,rcond);
    put(file,"The estimate of the inverse condition number : ");
    put(file,rcond,3); new_line(file);
    put_line(file,"The interpolating polynomial : ");
    put_line(file,ip);
    put(file,"Distance expanded Newton form and exact : ");
    put(file,Distance(dense,eq),3); new_line(file);
    put(file,"Distance direct interpolator and exact  : ");
    put(file,Distance(dense,ip),3); new_line(file);
  end Multprec_Experiment;

-- TEST EVALUATION AND SYMBOLIC EXPANSION :

  procedure Test_Standard_Eval
               ( file : in file_type;
                 q : in Standard_Divided_Differences.Newton_Interpolator1;
                 eq : in Standard_Divided_Differences.Newton_Form_Evaluator1;
                 x : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Evaluates both polynomial at x and compares the evaluations.

    use Standard_Complex_Numbers;
    y1 : constant Complex_Number := Standard_Divided_Differences.Eval(q,x);
    y2 : constant Complex_Number := Standard_Divided_Differences.Eval(eq,x);
    diff : constant Complex_Number := y1-y2;

  begin
    put(file,"interpolator : "); put(file,y1); new_line(file);
    put(file,"at evaluator : "); put(file,y2); new_line(file);
    put(file,"difference   : "); put(file,diff); new_line(file);
  end Test_Standard_Eval;

  procedure Test_Multprec_Eval
               ( file : in file_type;
                 q : in Multprec_Divided_Differences.Newton_Interpolator1;
                 eq : in Multprec_Divided_Differences.Newton_Form_Evaluator1;
                 x : in Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Evaluates both polynomial at x and compares the evaluations.

    use Multprec_Complex_Numbers;
    y1 : constant Complex_Number := Multprec_Divided_Differences.Eval(q,x);
    y2 : constant Complex_Number := Multprec_Divided_Differences.Eval(eq,x);
    diff : constant Complex_Number := y1-y2;

  begin
    put(file,"interpolator : "); put(file,y1); new_line(file);
    put(file,"at evaluator : "); put(file,y2); new_line(file);
    put(file,"difference   : "); put(file,diff); new_line(file);
  end Test_Multprec_Eval;

  procedure Shadow_Evaluation
               ( file : in file_type;
                 q : in Divided_Differences.Newton_Interpolator;
                 stgrid : in Array_of_Standard_Sample_Lists;
                 mpgrid : in Array_of_Multprec_Sample_Lists ) is

  -- DESCRIPTION :
  --   Evaluates the Newton interpolator at all points in the grids.

    stptr : Standard_Sample_List;
    mpptr : Multprec_Sample_List;
    steva : Standard_Complex_Numbers.Complex_Number;
    mpeva : Multprec_Complex_Numbers.Complex_Number;

  begin
    for i in stgrid'range loop
      put(file,"Evaluating at sample list "); put(file,i,1);
      put_line(file," :");
      stptr := stgrid(i);        mpptr := mpgrid(i);
      while not Is_Null(stptr) loop
        Divided_Differences.Eval
          (file,q,Sample_Point(Head_Of(stptr)).v,
                  Sample_Point(Head_Of(mpptr)).v,steva,mpeva);
        put(file,"ST eva : "); put(file,steva); new_line(file);
        put(file,"MP eva : "); put(file,mpeva); new_line(file);
        stptr := Tail_Of(stptr); mpptr := Tail_Of(mpptr);
      end loop;
    end loop;
  end Shadow_Evaluation;

  procedure Test_Standard_Symbolic_Expansion
               ( file : in file_type; 
                 q : in Standard_Divided_Differences.Newton_Form_Evaluator1;
                 exq : out Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Creates the symbolic expansion of the table of divided differences.

    use Standard_Divided_Differences;
    timer : Timing_Widget;
    dvd : Standard_Complex_Poly_Systems.Poly_Sys(0..Degree(q));
    expq : Standard_Complex_Polynomials.Poly;

  begin
    tstart(timer);
    dvd := Expand(q);
    expq := Expand(q,dvd);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Symbolic expansion of the Newton evaluator");
    new_line(file);
    put_line(file,"The symbolic expansion of the divided differences : ");
    put_line(file,dvd);
    put_line(file,"The symbolic expansion of the Newton form : ");
    put_line(file,expq);
    exq := expq;
  end Test_Standard_Symbolic_Expansion;

  procedure Test_Multprec_Symbolic_Expansion
               ( file : in file_type; 
                 q : in Multprec_Divided_Differences.Newton_Form_Evaluator1;
                 exq : out Multprec_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Creates the symbolic expansion of the table of divided differences.

    use Multprec_Divided_Differences;
    timer : Timing_Widget;
    dvd : Multprec_Complex_Poly_Systems.Poly_Sys(0..Degree(q));
    expq : Multprec_Complex_Polynomials.Poly;

  begin
    tstart(timer);
    dvd := Expand(q);
    expq := Expand(q,dvd);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Symbolic expansion of the Newton evaluator");
    new_line(file);
    put_line(file,"The symbolic expansion of the divided differences : ");
    put_line(file,dvd);
    put_line(file,"The symbolic expansion of the Newton form : ");
    put_line(file,expq);
    exq := expq;
  end Test_Multprec_Symbolic_Expansion;

  procedure Test_Standard_Newton_Taylor_Form1
               ( file : in file_type; n : in integer32;
                 q : in Standard_Divided_Differences.Newton_Interpolator1;
                 c : in Standard_Complex_Numbers.Complex_Number; 
                 grid : in Array_of_Standard_Sample_Lists;
                 testsps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Test on the Newton-Taylor form of the interpolating polynomial.

    use Standard_Divided_Differences;
    use Standard_Complex_Numbers;

    y : Standard_Complex_Vectors.Vector(0..Degree(q))
      := Random_Vector(0,Degree(q));
    x : Standard_Complex_Vectors.Vector(1..n);
    tq : Newton_Taylor_Form1;
    tmp : Standard_Sample_List;
    ev1,ev2,diff : Complex_Number;
    timer : Timing_Widget;

  begin
    y(0) := c;
    tstart(timer);
    tq := Create(q,y);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of 1-D Newton Taylor form");
    new_line(file);
    put_line(file,"Evaluation of Newton-Taylor form at grid : ");
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        put(file,Eval(tq,Sample_Point(Head_Of(tmp)).v)); new_line(file);
        tmp := Tail_Of(tmp);
      end loop;
      put(file,"---------------------------------------------------");
      put(file," at sample set "); put(file,i,1); new_line(file);
    end loop;
    put_line(file,"Evaluation of Newton-Taylor form at test samples : ");
    tmp := testsps;
    while not Is_Null(tmp) loop
      put(file,Eval(tq,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    for i in 0..5 loop
      if i = 0
       then put_line(file,"Evaluation at the origin : ");
            x := (1..n => Create(0.0));
       else put_line(file,"Evaluation at random point : ");
            x := Random_Vector(1,n);
      end if;
      put(file,"Taylor form  : "); ev1 := Eval(tq,x);
      put(file,ev1); new_line(file);
      put(file,"interpolator : "); ev2 := Eval(q,x);
      put(file,ev2); new_line(file);
      put(file,"difference   : "); diff := ev1 - ev2;
      put(file,diff); new_line(file);
    end loop;
  end Test_Standard_Newton_Taylor_Form1;

  procedure Test_Multprec_Newton_Taylor_Form1
               ( file : in file_type;
                 n : in integer32; size : in natural32;
                 q : in Multprec_Divided_Differences.Newton_Interpolator1;
                 c : in Multprec_Complex_Numbers.Complex_Number; 
                 grid : in Array_of_Multprec_Sample_Lists;
                 testsps : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Test on the Newton-Taylor form of the interpolating polynomial.

    use Multprec_Divided_Differences;
    use Multprec_Complex_Numbers;

    timer : Timing_Widget;
    y : Multprec_Complex_Vectors.Vector(0..Degree(q))
      := Random_Vector(0,Degree(q),size);
    x : Multprec_Complex_Vectors.Vector(1..n);
    tq : Newton_Taylor_Form1;
    tmp : Multprec_Sample_List;
    ev1,ev2,diff : Complex_Number;

  begin
    y(0) := c;
    Set_Size(y(0),size);
    tstart(timer);
    tq := Create(q,y);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of 1-D Newton-Taylor form");
    new_line(file);
    put_line(file,"Evaluation of Newton-Taylor form at grid : ");
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        put(file,Eval(tq,Sample_Point(Head_Of(tmp)).v)); new_line(file);
        tmp := Tail_Of(tmp);
      end loop;
      put(file,"---------------------------------------------------");
      put(file," at sample set "); put(file,i,1); new_line(file);
    end loop;
    put_line(file,"Evaluation of Newton-Taylor form at test samples : ");
    tmp := testsps;
    while not Is_Null(tmp) loop
      put(file,Eval(tq,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    for i in 0..5 loop
      if i = 0 then
        put_line(file,"Evaluation at the origin : ");
        x := (1..n => Create(integer(0)));
      else
        put_line(file,"Evaluation at random point : ");
        x := Random_Vector(1,n,size);
      end if;
      put(file,"Taylor form  : "); ev1 := Eval(tq,x);
      put(file,ev1); new_line(file);
      put(file,"interpolator : "); ev2 := Eval(q,x);
      put(file,ev2); new_line(file);
      put(file,"difference   : "); diff := ev1 - ev2;
      put(file,diff); new_line(file);
      Clear(ev1); Clear(ev2); Clear(diff);
      Multprec_Complex_Vectors.Clear(x);
    end loop;
  end Test_Multprec_Newton_Taylor_Form1;

  procedure Standard_Newton_Interpolation1
               ( file : in file_type; n : natural32;
                 c : in Standard_Complex_Numbers.Complex_Number;
                 grid : in Array_of_Standard_Sample_Lists;
                 testsps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Creates the Newton interpolator from the grid and evaluates
  --   at the grid and at the test points.  Also runs a test on the
  --   Newton Taylor form for a one dimensional component.

    use Standard_Divided_Differences;
    timer : Timing_Widget;
    q : Newton_Interpolator1;
    eq : Newton_Form_Evaluator1;
    exq : Standard_Complex_Polynomials.Poly;
    tmp : Standard_Sample_List;
    x : Standard_Complex_Vectors.Vector(1..integer32(n));

  begin
    tstart(timer);
    q := Create(grid,c);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of Newton interpolator");
    new_line(file);
    tstart(timer);
    eq := Create(q);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of Newton evaluator form");
    new_line(file);
    put_line(file,"Evaluating the Newton interpolator at the grid : ");
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        Test_Standard_Eval(file,q,eq,Sample_Point(Head_Of(tmp)).v);
        tmp := Tail_Of(tmp);
      end loop;
      put(file,"---------------------------------------------------");
      put(file," at sample set "); put(file,i,1); new_line(file);
    end loop;
    put_line(file,"Matrix of residuals after evaluation at grid : ");
    put(file,Errors(q,grid),3);
    put(file,"Maximal residual of evaluation at grid : ");
    put(file,Maximal_Error(q,grid),3); new_line(file);
    put_line(file,"Evaluating the Newton interpolator at the test samples : ");
    tmp := testsps;
    while not Is_Null(tmp) loop
      Test_Standard_Eval(file,q,eq,Sample_Point(Head_Of(tmp)).v);
      tmp := Tail_Of(tmp);
    end loop;
    put_line(file,"Evaluation at 5 random points : ");
    for i in 1..5 loop
      x := Random_Vector(1,integer32(n));
      Test_Standard_Eval(file,q,eq,x);
    end loop;
    Test_Standard_Newton_Taylor_Form1(file,integer32(n),q,c,grid,testsps);
    Test_Standard_Symbolic_Expansion(file,eq,exq);
    Standard_Experiment(file,exq,grid);
  end Standard_Newton_Interpolation1;

  procedure Multprec_Newton_Interpolation1
               ( file : in file_type; n : in integer32; size : in natural32;
                 c : in Multprec_Complex_Numbers.Complex_Number;
                 grid : in Array_of_Multprec_Sample_Lists;
                 testsps : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Creates the Newton interpolator from the grid and evaluates
  --   at the test points.

    use Multprec_Divided_Differences;
    timer : Timing_Widget;
    q : Newton_Interpolator1;
    eq : Newton_Form_Evaluator1;
    exq : Multprec_Complex_Polynomials.Poly;
    tmp : Multprec_Sample_List;
    x : Multprec_Complex_Vectors.Vector(1..n);

  begin
    tstart(timer);
    q := Create(grid,c);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of 1-D Newton interpolator");
    new_line(file);
    tstart(timer);
    eq := Create(q);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of Newton evaluator form");
    new_line(file);
    put_line(file,"Evaluating the Newton interpolator at the grid : ");
    for i in grid'range loop
      tmp := grid(i);
      while not Is_Null(tmp) loop
        Test_Multprec_Eval(file,q,eq,Sample_Point(Head_Of(tmp)).v);
        tmp := Tail_Of(tmp);
      end loop;
      put(file,"---------------------------------------------------");
      put(file," at sample set "); put(file,i,1); new_line(file);
    end loop;
    put_line(file,"Evaluating the Newton interpolator at the test samples : ");
    tmp := testsps;
    while not Is_Null(tmp) loop
      Test_Multprec_Eval(file,q,eq,Sample_Point(Head_Of(tmp)).v);
      tmp := Tail_Of(tmp);
    end loop;
    put_line(file,"Evaluation at 5 random points : ");
    for i in 1..5 loop
      x := Random_Vector(1,n,size);
      Test_Multprec_Eval(file,q,eq,x);
      Multprec_Complex_Vectors.Clear(x);
    end loop;
    Test_Multprec_Newton_Taylor_Form1(file,n,size,q,c,grid,testsps);
    Test_Multprec_Symbolic_Expansion(file,eq,exq);
    Multprec_Experiment(file,exq,grid);
  end Multprec_Newton_Interpolation1;

  procedure Standard_Write_Evaluations
               ( file : in file_type;
                 q : in Standard_Divided_Differences.Newton_Taylor_Form;
                 sps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of q at all samples in sps
  --   on file, each time taking a new line for every evaluation.

    use Standard_Divided_Differences;
    tmp : Standard_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(q,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Write_Evaluations;

  procedure Multprec_Write_Evaluations
               ( file : in file_type;
                 q : in Multprec_Divided_Differences.Newton_Taylor_Form;
                 sps : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of q at all samples in sps
  --   on file, each time taking a new line for every evaluation.

    use Multprec_Complex_Numbers;
    use Multprec_Divided_Differences;
    tmp : Multprec_Sample_List := sps;
    eva : Multprec_Complex_Numbers.Complex_Number;

  begin
    while not Is_Null(tmp) loop
      eva := Eval(q,Sample_Point(Head_Of(tmp)).v);
      put(file,eva); new_line(file);
      Clear(eva);
      tmp := Tail_Of(tmp);
    end loop;
  end Multprec_Write_Evaluations;

  procedure Standard_Test_Grid_Evaluation
               ( file : in file_type;
                 q : in Standard_Divided_Differences.Newton_Taylor_Form;
                 grid : in Standard_Stacked_Sample_Grids.
                           Stacked_Sample_Grid ) is

  -- DESCRIPTION :
  --   Evaluates the Newton-Taylor form of the interpolator at
  --   all grid points.

    use Standard_Divided_Differences;
    use Standard_Stacked_Sample_Grids;

  begin
    if grid.k = 1
     then put_line(file,"Evaluation at dimension 1 in the grid :");
          for i in grid.g'range loop
            put(file,"At sample list "); put(file,i,1); put_line(file," :");
            Standard_Write_Evaluations(file,q,grid.g(i));
          end loop;
     else put(file,"Evaluation at dimension "); put(file,grid.k,1);
          put(file," for degrees "); put(file,grid.d,1);
          put_line(file," down to 1 :");
          for i in reverse 1..grid.d loop
            Standard_Test_Grid_Evaluation(file,q,grid.a(i).all);
          end loop;
          put_line(file,"Evaluation at the last sample : ");
          put(file,Eval(q,Sample_Point(grid.spt).v)); new_line(file);
    end if;
  end Standard_Test_Grid_Evaluation;

  procedure Multprec_Test_Grid_Evaluation
               ( file : in file_type;
                 q : in Multprec_Divided_Differences.Newton_Taylor_Form;
                 grid : in Multprec_Stacked_Sample_Grids.
                           Stacked_Sample_Grid ) is

  -- DESCRIPTION :
  --   Evaluates the Newton-Taylor form of the interpolator at
  --   all grid points.

    use Multprec_Complex_Numbers;
    use Multprec_Divided_Differences;
    use Multprec_Stacked_Sample_Grids;

    eva : Multprec_Complex_Numbers.Complex_Number;

  begin
    if grid.k = 1
     then put_line(file,"Evaluation at dimension 1 in the grid :");
          for i in grid.g'range loop
            put(file,"At sample list "); put(file,i,1); put_line(file," :");
            Multprec_Write_Evaluations(file,q,grid.g(i));
          end loop;
     else put(file,"Evaluation at dimension "); put(file,grid.k,1);
          put(file," for degrees "); put(file,grid.d,1);
          put_line(file," down to 1 :");
          for i in reverse 1..grid.d loop
            Multprec_Test_Grid_Evaluation(file,q,grid.a(i).all);
          end loop;
          put_line(file,"Evaluation at the last sample : ");
          eva := Eval(q,Sample_Point(grid.spt).v);
          put(file,eva); new_line(file);
          Clear(eva);
    end if;
  end Multprec_Test_Grid_Evaluation;

  procedure Standard_Newton_Interpolation
               ( file : in file_type;
                 grid : in Standard_Stacked_Sample_Grids.Stacked_Sample_Grid;
                 c : Standard_Complex_Numbers.Complex_Number;
                 y : Standard_Complex_Vectors.Vector;
                 testsps : in Standard_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the Newton Taylor form of the interpolating polynomial
  --   through the samples in the grid.  Tests whether the interpolator
  --   evaluates to zero at the grid and the test samples.

    use Standard_Stacked_Sample_Grids;
    use Standard_Divided_Differences;
    timer : Timing_Widget;
    q : Newton_Taylor_Form;

  begin
    tstart(timer);
    if output
     then q := Create(file,grid,c,y);
     else q := Create(grid,c,y);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of Newton-Taylor form");
    new_line(file);
    Standard_Test_Grid_Evaluation(file,q,grid);
    put_line(file,"Evaluation of Newton-Taylor form at test samples :");
    Standard_Write_Evaluations(file,q,testsps);
  end Standard_Newton_Interpolation;

  procedure Multprec_Newton_Interpolation
               ( file : in file_type;
                 grid : in Multprec_Stacked_Sample_Grids.Stacked_Sample_Grid;
                 c : Multprec_Complex_Numbers.Complex_Number;
                 y : Multprec_Complex_Vectors.Vector;
                 testsps : in Multprec_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the Newton Taylor form of the interpolating polynomial
  --   through the samples in the grid.  Tests whether the interpolator
  --   evaluates to zero at the grid and the test samples.

    use Multprec_Stacked_Sample_Grids;
    use Multprec_Divided_Differences;
    timer : Timing_Widget;
    q : Newton_Taylor_Form;

  begin
    tstart(timer);
    if output
     then q := Create(file,grid,c,y);
     else q := Create(grid,c,y);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Creation of Newton-Taylor form");
    new_line(file);
    Multprec_Test_Grid_Evaluation(file,q,grid);
    put_line(file,"Evaluation of Newton-Taylor form at test samples :");
    Multprec_Write_Evaluations(file,q,testsps);
  end Multprec_Newton_Interpolation;

-- MAIN PROGRAMS :

  procedure Test_Standard_Newton_Interpolation1
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim : in natural32 ) is

  -- DESCRIPTION :
  --   This tests the one dimensional Newton interpolation.

    n : constant integer32 := p'last;
    sli : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    testsps : Standard_Sample_List;
    len : constant integer32 := integer32(Length_Of(sps));
    grid : Array_of_Standard_Sample_Lists(0..len);
    eps,dst : double_float;
    c : Standard_Complex_Numbers.Complex_Number := Random1;
    ans : character;

  begin
    loop
      put("Normalization constant is "); put(c); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'n');
      if ans = 'y'
       then put("Give complex number : "); get(c);
      end if;
    end loop;
    put_line(file,"The embedded polynomial system :");
    put(file,p);
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    Standard_Test_Samples(file,sps,sli,testsps);
    Standard_Rectangular_Grid_Creator(file,sps,natural32(len),grid,eps,dst);
    Standard_Newton_Interpolation1(file,natural32(n),c,grid,testsps);
    Sampling_Machine.Clear;
  end Test_Standard_Newton_Interpolation1;

  procedure Test_Standard_Newton_Interpolation
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim : in natural32 ) is

  -- DESCRIPTION :
  --   This tests Newton interpolation for any dimension with
  --   standard floating-point numbers.

    use Standard_Stacked_Sample_Grids;

   -- n : constant integer32 := p'last;
    deg : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    gsz : constant natural32 := Grid_Size(dim+1,deg);
    sli : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    testsps : Standard_Sample_List;
    grid : Stacked_Sample_Grid(integer32(dim),integer32(deg));
    c : Standard_Complex_Numbers.Complex_Number := Random1;
    ans : character;
    y : Standard_Complex_Vectors.Vector(0..integer32(deg))
      := Random_Vector(0,integer32(deg));
    output : boolean;

  begin
    put("The number of samples : "); put(gsz,1); new_line;
    loop
      put("Normalization constant is "); put(c); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'n');
      if ans = 'y'
       then put("Give complex number : "); get(c);
      end if;
    end loop;
    put("Do you want intermediate output of interpolation tests ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    new_line;
    y(0) := c;
    put_line(file,"The embedded polynomial system :");
    put(file,p);
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    Standard_Test_Samples(file,sps,sli,testsps);
    Standard_Stacked_Grid_Creator(file,sps,false,grid);
    Standard_Newton_Interpolation(file,grid,c,y,testsps,output);
    Sampling_Machine.Clear;
  end Test_Standard_Newton_Interpolation;

  procedure Test_Multprec_Newton_Interpolation1
               ( file : in file_type;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim,size : in natural32 ) is

  -- DESCRIPTION :
  --   Test on Newton interpolation for curves with multi-precision
  --   floating-point numbers.

    n : constant integer32 := ep'last;
    sli : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(ep,dim);
    sps : Standard_Sample_List := Create(sols,sli);
    len : constant integer32 := integer32(Length_Of(sps));
    testsps : Multprec_Sample_List;
    grid : Array_of_Multprec_Sample_Lists(0..len);
    eps,dst : double_float;
    c : Multprec_Complex_Numbers.Complex_Number := Create(Random1);
    ans : character;

  begin
    loop
      put("Normalization constant is "); put(c); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'n');
      if ans = 'y'
       then put("Give complex number : "); get(c);
      end if;
    end loop;
    put_line(file,"The embedded polynomial system :");
    put(file,ep);
    Sampling_Machine.Initialize(ep,mp,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Multprec_Test_Samples(file,sps,sli,testsps);
    Multprec_Rectangular_Grid_Creator
      (file,sps,natural32(len),size,grid,eps,dst);
    Multprec_Newton_Interpolation1(file,n,size,c,grid,testsps);
    Sampling_Machine.Clear;
  end Test_Multprec_Newton_Interpolation1;

  procedure Test_Multprec_Newton_Interpolation
               ( file : in file_type;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim,size : in natural32 ) is

  -- DESCRIPTION :
  --   Test on Newton interpolation for any dimension with
  --   multi-precision floating-point numbers.

    use Multprec_Stacked_Sample_Grids;

   -- n : constant integer32 := ep'last;
    deg : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    sli : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(ep,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    testsps : Multprec_Sample_List;
    grid : Stacked_Sample_Grid(integer32(dim),integer32(deg));
    c : Multprec_Complex_Numbers.Complex_Number := Create(Random1);
    ans : character;
    y : Multprec_Complex_Vectors.Vector(0..integer32(deg))
      := Random_Vector(0,integer32(deg),size);
    output : boolean;

  begin
    loop
      put("Normalization constant is "); put(c); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'n');
      if ans = 'y'
       then put("Give complex number : "); get(c);
      end if;
    end loop;
    Set_Size(c,size);
    y(0) := c;
    put("Do you want intermediate output of interpolation tests ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    new_line;
    put_line(file,"The embedded polynomial system :");
    put(file,ep);
    Sampling_Machine.Initialize(ep,mp,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Multprec_Test_Samples(file,sps,sli,testsps);
    Multprec_Stacked_Grid_Creator(file,sps,false,size,grid);
    Multprec_Newton_Interpolation(file,grid,c,y,testsps,output);
    Sampling_Machine.Clear;
  end Test_Multprec_Newton_Interpolation;

  procedure Shadow_Test
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim,size : in natural32 ) is

  -- DESCRIPTION :
  --   Calls the operations in "Divided_Differences" to shadow the
  --   multi-precision calculations with standard arithmetic.

    use Divided_Differences;

    sli : constant  Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(ep,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    stgrid : Array_of_Standard_Sample_Lists(0..integer32(Length_Of(sps)));
    mpgrid : Array_of_Multprec_Sample_Lists(0..integer32(Length_Of(sps)));
    q : Newton_Interpolator;
    c : Standard_Complex_Numbers.Complex_Number := Random1;
    ans : character;

  begin
    loop
      put("Normalization constant is "); put(c); new_line;
      put("Do you want to change this value ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans = 'n');
      if ans = 'y'
       then put("Give complex number : "); get(c);
      end if;
    end loop;
    put_line(file,"The embedded polynomial system :");
    put(file,ep);
    Sampling_Machine.Initialize(ep,mp,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Rectangular_Sample_Grids.Create1(sps,Length_Of(sps),size,stgrid,mpgrid);
    q := Create(file,size,stgrid,mpgrid);
    Shadow_Evaluation(file,q,stgrid,mpgrid);
    Sampling_Machine.Clear;
  end Shadow_Test;

  procedure Main is

    outfile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim,deci,size : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Multivariate interpolation with divided differences.");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(lp.all,mp,size,dim);
      if dim = 1 then
        new_line;
        put("Shadow the multi-precision calculations ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          Shadow_Test(outfile,lp.all,mp.all,sols,dim,size);
          -- elsif dim = 1
          --     then Test_Multprec_Newton_Interpolation1
          --            (outfile,lp.all,mp.all,sols,dim,size);
        else
          Test_Multprec_Newton_Interpolation
            (outfile,lp.all,mp.all,sols,dim,size);
        end if;
      else
        Test_Multprec_Newton_Interpolation
          (outfile,lp.all,mp.all,sols,dim,size);
      end if;
    -- elsif dim = 1
    --     then Test_Standard_Newton_Interpolation1(outfile,lp.all,sols,dim);
    else 
      Test_Standard_Newton_Interpolation(outfile,lp.all,sols,dim);
    end if;
  end Main;

begin
  Main;
end ts_divdif;
