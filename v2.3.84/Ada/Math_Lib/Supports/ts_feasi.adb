with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;      use Standard_Integer_Vectors_io;
with Standard_Floating_Numbers;        use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;     use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;          use Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;     use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;    use Standard_Floating_Matrices_io;
with Standard_Random_Vectors;          use Standard_Random_Vectors;
with Standard_Random_Matrices;         use Standard_Random_Matrices;
with Linear_Minimization;              use Linear_Minimization;

procedure ts_feasi is 

-- DESCRIPTION :
--   Test on the correctness of the linear minimization routines.
--   The standard test is to solve a linear system with one common
--   point that satisfies all linear inequalities.

  function Average ( cff : Standard_Floating_Matrices.Matrix )
                   return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the average of the columns in the coefficient matrix.

    res : Standard_Floating_Vectors.Vector(cff'range(1))
        := (cff'range(1) => 0.0);

  begin
    for i in cff'range(1) loop
      for j in cff'range(2) loop
        res(i) := res(i) + cff(i,j);
      end loop;
    end loop;
    for i in res'range loop
      res(i) := res(i)/double_float(cff'last(2));
    end loop;
    return res;
  end Average;

  function Positive_Random ( n : integer32 )
                           return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a positive random vector of range 1..n.

    res : Standard_Floating_Vectors.Vector(1..n)
        := Random_Vector(1,n);

  begin
    for i in res'range loop
      res(i) := abs(res(i));
    end loop;
    return res;
  end Positive_Random;

  function Positive_Random ( nrows,ncols : integer32 )
                           return Standard_Floating_Matrices.Matrix is
 
  -- DESCRIPTION :
  --   Returns a matrix with all entries random positive.

    res : Standard_Floating_Matrices.Matrix(1..nrows,1..ncols)
        := Random_Matrix(natural32(nrows),natural32(ncols));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := abs(res(i,j));
      end loop;
    end loop;
    return res;
  end Positive_Random;

  function Adjusted_Inequalities 
              ( cff : Standard_Floating_Matrices.Matrix;
                sol : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns right hand side vector for the inequalities with
  --   constraints in the coefficient matrix so that the given vector
  --   satisfies all the constraints as equalities.

    res : Standard_Floating_Vectors.Vector(cff'range(2));

  begin
    for j in cff'range(2) loop
      res(j) := 0.0;
      for i in cff'range(1) loop
        res(j) := res(j) + sol(i)*cff(i,j);
      end loop;
    end loop;
    return res;
  end Adjusted_Inequalities;

  procedure Random_Input
               ( nvars,ncons : in integer32;
                 cff : out Standard_Floating_Matrices.Matrix;
                 cost,rhs,sol,init : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Randomly chosen positive numbers are generated such that a positive
  --   vector is the optimal solution to the minimization problem.

  -- ON ENTRY :
  --   nvars     number of variables;
  --   ncons     number of constraints;

  -- ON RETURN :
  --   cff       nvars*ncons random positive coefficient matrix;
  --   cost      cost vector is [1 1 .. 1], vector of range 1..nvars;
  --   rhs       right-hand side vector of range 1..ncons;
  --   sol       optimal solution to the problem;
  --   init      initial, perturbed solution, feasible but not optimal.

    use Standard_Floating_Vectors;

  begin
    cff := Positive_Random(nvars,ncons);
    cost := (1..nvars => 1.0);
    sol := Positive_Random(nvars);
    rhs := Adjusted_Inequalities(cff,sol);
    init := sol + Average(cff);
  end Random_Input;

  procedure User_Input
               ( nvars,ncons : in integer32;
                 cff : out Standard_Floating_Matrices.Matrix;
                 cost,rhs : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   User is prompted to give in the linear minimization problem.

  -- ON ENTRY :
  --   nvars     number of variables;
  --   ncons     number of constraints;

  -- ON RETURN :
  --   cff       nvars*ncons random positive coefficient matrix;
  --   cost      cost vector is [1 1 .. 1], vector of range 1..nvars;
  --   rhs       right-hand side vector of range 1..ncons.

  begin
    put("Give a "); put(nvars,1); put("*"); put(ncons,1);
    put_line("-coefficient matrix :"); get(cff);
    put("Give a "); put(nvars,1);
    put_line("-vector for cost :"); get(cost);
    put("Give a "); put(ncons,1); 
    put_line("-vector for right-hand side : "); get(rhs);
  end User_Input;

  procedure Interactive_Input
               ( nvars,ncons : in integer32;
                 cff : out Standard_Floating_Matrices.Matrix;
                 cost,rhs,init : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   User is prompted to give minimization problem and initial solution.

  -- ON ENTRY :
  --   nvars     number of variables;
  --   ncons     number of constraints;

  -- ON RETURN :
  --   cff       nvars*ncons random positive coefficient matrix;
  --   cost      cost vector is [1 1 .. 1], vector of range 1..nvars;
  --   rhs       right-hand side vector of range 1..ncons;
  --   init      initial, perturbed solution, feasible but not optimal.

  begin
    User_Input(nvars,ncons,cff,cost,rhs);
    put("Give a "); put(nvars,1);
    put_line("-vector for initial solution : "); get(init);
  end Interactive_Input;

  procedure Test_Constraints 
              ( cffmat : in Standard_Floating_Matrices.Matrix;
                rhs,sol : in Standard_Floating_Vectors.Vector;
                tol : in double_float; bug : out boolean ) is

  -- DESCRIPTION :
  --   Tests whether all constraints hold, bug is set to true if
  --   some constraints are violated.

    val : double_float;

  begin
    put_line("evaluation in constraints minus rhs :");
    for j in cffmat'range(2) loop
      val := 0.0;
      for i in cffmat'range(1) loop
        val := val + sol(i)*cffmat(i,j);
      end loop;
      put(val,3); put(" - "); put(rhs(j),3); put(" = ");
      val := val-rhs(j);
      put(val,3);
      if val > -tol then
        put_line(" >= 0.  ok");
        bug := false;
      else
        put_line(" < 0.  BUG!!!");
        bug := true;
      end if;
    end loop;
  end Test_Constraints;

  procedure Minimize_Once
              ( nvars,ncons : in integer32; randinp : in boolean ) is

    cffmat : Standard_Floating_Matrices.Matrix(1..nvars,1..ncons);
    cost : Standard_Floating_Vectors.Vector(1..nvars) := (1..nvars => 1.0);
    x,init : Standard_Floating_Vectors.Vector(1..nvars);
    rhs : Standard_Floating_Vectors.Vector(1..ncons);
    tol : constant double_Float := 10.0**(-8);
    infty,continue : boolean;
    step : integer32 := -1;
    zeros : constant Standard_Integer_Vectors.Vector := (1..nvars => 0);

    procedure Write ( sol : in Standard_Floating_Vectors.Vector;
                      active : in Standard_Integer_Vectors.Vector;
                      cont : out boolean ) is

      val : double_float := 0.0;
      bug : boolean := false;

    begin
      cont := true;
      step := step + 1;
      put("The solution for columns "); put(active);
      put(" at step "); put(step,1); put_line(" :");
      put_line(sol);
      for i in sol'range loop
        val := val + cost(i)*sol(i);
      end loop;
      put("evaluation in cost : "); put(val); new_line;
      Test_Constraints(cffmat,rhs,sol,tol,bug);
      cont := not bug;
    end Write;
    procedure Rep_Minimize is new Reporting_Minimize(Write);

  begin
    if randinp
     then Random_Input(nvars,ncons,cffmat,cost,rhs,x,init);
     else Interactive_Input(nvars,ncons,cffmat,cost,rhs,init);
    end if;
    put_line("The set up of the problem :");
    put_line("matrix of constraints : "); put(cffmat,3);
    put_line("cost vector :"); put(cost,3); new_line;
    put_line("right-hand side vector : "); put(rhs,3); new_line;
    if randinp then
      step := 1000;
      Write(x,zeros,continue);
    end if;
    put_line("Starting the optimization :");
    step := -1;
    Write(init,zeros,continue);
    Rep_Minimize(nvars,ncons,cffmat,cost,rhs,init,tol,infty);
  end Minimize_Once;

  procedure Random_Test_Minimizer
              ( nvars,ncons : in integer32; bug : out boolean ) is

    cffmat : Standard_Floating_Matrices.Matrix(1..nvars,1..ncons);
    cost : Standard_Floating_Vectors.Vector(1..nvars) := (1..nvars => 1.0);
    x,init : Standard_Floating_Vectors.Vector(1..nvars);
    rhs : Standard_Floating_Vectors.Vector(1..ncons);
    tol : constant double_Float := 10.0**(-8);
    infty : boolean;
    first : boolean := true;
    prev_cost : double_float;

    procedure Test ( sol : in Standard_Floating_Vectors.Vector;
                     active : in Standard_Integer_Vectors.Vector;
                     cont : out boolean ) is

    -- DESCRIPTION :
    --   Tests whether the cost has decreased for this new solution
    --   and whether the vector is still a feasible solution.

      val : double_float := 0.0;

    begin
      cont := true; bug := false;
      for j in cffmat'range(2) loop             -- test feasibility
        val := 0.0;
        for i in cffmat'range(1) loop
          val := val + sol(i)*cffmat(i,j);
        end loop;
        val := val-rhs(j);
        if val < -tol then
          bug := true;
          put_line("BUG : SOLUTION NOT FEASIBLE !!!");
          cont := false;
        end if;
      end loop;
      if cont then                                -- test cost decrease
        for i in sol'range loop
          val := val + cost(i)*sol(i);
        end loop;
        if first then
          prev_cost := val; first := false;
        elsif prev_cost < val+tol then
          bug := true;
          put_line("BUG : SOLUTION WITH INCREASING COST !!!");
          cont := false;
        end if;
      end if;
    end Test;
    procedure Test_Minimize is new Reporting_Minimize(Test);

  begin
    Random_Input(nvars,ncons,cffmat,cost,rhs,x,init);
    Test_Minimize(nvars,ncons,cffmat,cost,rhs,init,tol,infty);
  end Random_Test_Minimizer;

  procedure Minimize_Frequently ( nvars,ncons : in integer32 ) is

    bug : boolean := false;
    num : integer32 := 0;

  begin
    new_line;
    put("Give the number of optimization problems : "); get(num);
    new_line;
    put("Starting the test of "); put(num,1); put_line(" minimizations...");
    new_line;
    for i in 1..num loop
      Random_Test_Minimizer(nvars,ncons,bug);
      exit when bug;
      put("Instance "); put(i,1); put_line(" is ok.");
    end loop;
    if bug
     then put_line("BUG FOUND !!!");
     else put("Tested "); put(num,1); put_line(" minimizations without bug.");
    end if;
  end Minimize_Frequently;

  procedure Interactive_Feasibility_Test ( nvars,ncons : in integer32 ) is

    tol : constant double_float := 10.0**(-8);
    cff : Standard_Floating_Matrices.Matrix(1..nvars,1..ncons);
    rhs : Standard_Floating_Vectors.Vector(1..ncons);
    sol : Standard_Floating_Vectors.Vector(1..nvars);
    fail : boolean;

  begin
    put("Give "); put(nvars,1); put("*"); put(ncons,1);
    put_line("-matrix of coefficients :"); get(cff);
    put("Give "); put(ncons,1);
    put_line("-vector with right-hand side :"); get(rhs);
    Feasible(nvars,ncons,cff,rhs,tol,sol,fail);
    if fail
     then put_line("The system is not feasible.");
     else put_line("The system has a solution :");
          put_line(sol);
    end if;
  end Interactive_Feasibility_Test;

  procedure Random_Feasibility_Test
              ( nvars,ncons : in integer32; bug : out boolean ) is

    tol : constant double_float := 10.0**(-8);
    cff : Standard_Floating_Matrices.Matrix(1..nvars,1..ncons)
       -- := Positive_Random(nvars,ncons);
        := Random_Matrix(natural32(nvars),natural32(ncons));
    sol : Standard_Floating_Vectors.Vector(1..nvars)
       -- := Positive_Random(nvars);
        := Random_Vector(1,nvars);
    rhs : Standard_Floating_Vectors.Vector(1..ncons);
    cost : constant Standard_Floating_Vectors.Vector(1..nvars)
         := (1..nvars => 0.0);
    resi : Standard_Floating_Vectors.Vector(1..ncons);
    val : double_float;
    ran : constant double_float := Random;
    feasi : constant boolean := (ran >= 0.0);
    fail : boolean;

  begin
    for i in 1..nvars loop            -- make sure start solution is zero
      for j in 1..nvars loop          -- to enforce pivoting
        if i = j then
          if sol(i) > 0.0
           then cff(i,j) := +1.0;
           else cff(i,j) := -1.0;
          end if;
        else
          cff(i,j) := 0.0;
        end if;
      end loop;
      rhs(i) := 0.0;
    end loop;
    for j in nvars+1..ncons loop
      rhs(j) := 0.0;
      for i in 1..nvars loop
        rhs(j) := rhs(j) + sol(i)*cff(i,j);
      end loop;
     -- rhs(j) := rhs(j) - abs(Random);
    end loop;
    Eval(nvars,ncons,cff,rhs,cost,sol,resi,val);
   -- put_line("The residuals of the generated system : ");
   -- put_line(resi);
    if feasi then
      put("feasible system : ");
    else
      for i in 1..nvars loop
        cff(i,ncons) := -cff(i,1);    -- opposite normal vector
      end loop;
      rhs(1) := abs(rhs(1));
      rhs(ncons) := rhs(1) + 1.0;     -- half spaces do not intersect
      put("infeasible system : ");
    end if;
    Feasible(nvars,ncons,cff,rhs,tol,sol,fail);
    if fail then
      if not feasi
       then put_line("recognized as infeasible"); bug := false;
       else put_line("BUG!!! no solution found"); bug := true;
      end if;
    else
      if feasi
       then put_line("recognized as feasible"); bug := false;
       else put_line("BUG!!!  solution found"); bug := true;
      end if;
    end if;
  end Random_Feasibility_Test;

  procedure Random_Feasibility_Tests ( nvars,ncons : in integer32 ) is

    num : integer32 := 0;
    bug : boolean;

  begin
    put("Give the number of tests : "); get(num);
    for i in 1..num loop
      Random_Feasibility_Test(nvars,ncons,bug);
      exit when bug;
    end loop;
  end Random_Feasibility_Tests;

  procedure Interactive_Basis_Enumeration ( nvars,ncons : in integer32 ) is

    cff : Standard_Floating_Matrices.Matrix(1..nvars,1..ncons);
    cost : Standard_Floating_Vectors.Vector(1..nvars);
    rhs : Standard_Floating_Vectors.Vector(1..ncons);
    tol : constant double_float := 10.0**(-8);
    fail,infty : boolean;

    procedure Write ( binv : in Standard_Floating_Matrices.Matrix;
                      active : in Standard_Integer_Vectors.Vector;
                      sol : in Standard_Floating_Vectors.Vector;
                      continue : out boolean ) is

      val : double_float := 0.0;

    begin
      put(active); 
      put("  (");
      for i in sol'first..sol'last-1 loop
        put(sol(i),3); put(",");
      end loop;
      put(sol(sol'last),3);
      put(" )  value : ");
      for i in cost'range loop
        val := val + cost(i)*sol(i);
      end loop;
      put(val,3);
      new_line;
      continue := true;
    end Write;
    procedure Enumerate is new Enumerate_Feasible_Vertices(Write);

  begin
    User_Input(nvars,ncons,cff,cost,rhs);
    Enumerate(nvars,ncons,cff,cost,rhs,tol,fail,infty);
    if fail
     then put_line("The optimization problem is unfeasible.");
    end if;
    if infty
     then put_line("The optimization problem is unbounded.");
    end if;
  end Interactive_Basis_Enumeration;

  procedure Main is

    nvars,ncons : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Feasibility test on a system of linear inequalities.");
    put_line("   The model is min cost*x, subject to cff*x >= rhs.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. exit this program");
      put_line("  1. minimization on user input");
      put_line("  2. one time minimization on randomly generated data");
      put_line("  3. massive black-box minimization on random data");
      put_line("  4. feasibility test on user input");
      put_line("  5. massive black-box feasibility tests on random data");
      put_line("  6. enumeration of all feasible vertices on user input");
      put("Type 0,1,2,3,4,5, or 6 to select : "); get(ans);
      exit when (ans = '0');
      new_line;
      put("Give the number of variables : ");   get(nvars);
      put("Give the number of constraints : "); get(ncons);
      case ans is
        when '1' => Minimize_Once(nvars,ncons,false);
        when '2' => Minimize_Once(nvars,ncons,true);
        when '3' => Minimize_Frequently(nvars,ncons);
        when '4' => Interactive_Feasibility_Test(nvars,ncons);
        when '5' => Random_Feasibility_Tests(nvars,ncons);
        when '6' => Interactive_Basis_Enumeration(nvars,ncons);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_feasi;
