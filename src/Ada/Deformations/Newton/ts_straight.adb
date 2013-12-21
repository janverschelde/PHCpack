with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_VecMats;          use Standard_Complex_VecMats;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Deflate_Singularities;    use Standard_Deflate_Singularities;
with Standard_Jacobian_Trees;           use Standard_Jacobian_Trees;

procedure ts_straight is

-- DESCRIPTION :
--   Development of a straight-line version of the deflation.

-- GENERAL UTILITIES :

  function Stack ( a,b : Vector ) return Vector is

  -- DESCRIPTION :
  --   The vector of return contains b stacked to a.

  -- REQUIRED : a'first = 1 and b'first = 1.

    res : Vector(1..a'last+b'last);

  begin
    res(a'range) := a;
    for i in b'range loop
      res(a'last+i) := b(i);
    end loop;
    return res;
  end Stack;

  function Create ( jm : Jaco_Mat ) return Array_of_Jaco_Mat is

  -- DESCRIPTION :
  --   Returns the vector of Jacobian matrices, obtained by deriving
  --   the given Jacobian matrix jm with respect to all its columns.

    res : Array_of_Jaco_Mat(jm'range(2));

  begin
    for i in res'range loop
      res(i) := new Jaco_Mat'(Diff(jm,i));
    end loop;
    return res;
  end Create;

  function Create ( jm : Array_of_Jaco_Mat ) return Array_of_Eval_Jaco_Mat is

  -- DESCRIPTION :
  --   Returns evaluable forms of all matrices in jm.

    res : Array_of_Eval_Jaco_Mat(jm'range);

  begin
    for i in res'range loop
      res(i) := new Eval_Jaco_Mat'(Create(jm(i).all));
    end loop;
    return res;
  end Create;

  function Eval ( jm : Array_of_Eval_Jaco_Mat; x : Vector ) return VecMat is

  -- DESCRIPTION :
  --   Returns the vector of matrices, obtained after evaluation of every
  --   Jacobian matrix in jm at x.

    res : VecMat(jm'range);

  begin
    for i in res'range loop
      res(i) := new Matrix'(Eval(jm(i).all,x));
    end loop;
    return res;
  end Eval;

  function "*" ( a : VecMat; b : Vector ) return Matrix is

  -- DESCRIPTION :
  --   Returns the matrix with in its columns a(i)*b, for i in a'range.

  -- REQUIRED : a(i)'range(2) = b'range.

    col : Vector(a(a'first)'range(1));
    res : Matrix(col'range,a'range);

  begin
    for j in a'range loop
      col := a(j).all*b;
      for i in col'range loop
        res(i,j) := col(i);
      end loop;
    end loop;
    return res;
  end "*";

-- STRAIGHT-LINE PROGRAMS :

  function Straight_Evaluate_Deflation
              ( n : integer32; f : Eval_Poly_Sys; jf : Eval_Jaco_Mat;
                b : Matrix; h,x,lambda : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the value of the deflated system at x,
  --   with multipliers in lambda.

  -- ON ENTRY :
  --   n        number of equations in the deflated system,
  --            the range of the vector on return is 1..n;
  --   f        original polynomial system;
  --   jf       Jacobian matrix of the original system;
  --   b        a random matrix for multiplying Jacobian matrix with;
  --   h        a random vector to have unique values for multipliers;
  --   x        values for the original variables;
  --   lambda   values for the multipliers.

    res : Vector(1..n);
    y : constant Vector := Eval(f,x);
    A : constant Matrix(jf'range(1),jf'range(2)) := Eval(jf,x);
    bl : constant Vector := b*lambda;
    Abl : constant Vector := A*bl;
    ind : integer32 := y'last+1;

  begin
    res(y'range) := y;
    for i in Abl'range loop
      res(ind) := Abl(i);
      ind := ind + 1;
    end loop;
    res(n) := Create(-1.0);
    for i in lambda'range loop
      res(n) := res(n) + h(i)*lambda(i);
    end loop;
    return res;
  end Straight_Evaluate_Deflation;

  function Straight_Differentiate_Deflation
              ( rows,cols : integer32; jf : Eval_Jaco_Mat;
                ajf : Array_of_Eval_Jaco_Mat;
                b : Matrix; h,x,lambda : Vector ) return Matrix is

  -- DESCRIPTION :
  --   Returns the value of the Jacobian matrix of the deflated system
  --   evaluated at x with multipliers in lambda.

  -- ON ENTRY :
  --   rows     number of rows in the matrix on return;
  --   cols     number of columns in the matrix on return;
  --   jf       Jacobian matrix of the original polynomial system;
  --   ajf      all 1st derivatives of jf, the so-called "Hessian";
  --   b        random matrix to multiply Jacobian matrices with;
  --   h        a random vector to have unique values for multipliers;
  --   x        values for the original variables;
  --   lambda   values for the multipliers.

    res : Matrix(1..rows,1..cols);
    A : constant Matrix(jf'range(1),jf'range(2)) := Eval(jf,x);
    Ab : constant Matrix(A'range(1),b'range(2)) := A*b;
    dA : constant VecMat(ajf'range) := Eval(ajf,x);
    bl : constant Vector := b*lambda;
    dAbl : constant Matrix(A'range(1),A'range(2)) := dA*bl;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j);
        res(A'last(1)+i,j) := dAbl(i,j);
        res(rows,j) := Create(0.0);
      end loop;
      for j in Ab'range(2) loop
        res(i,A'last(2)+j) := Create(0.0);
        res(A'last(1)+i,A'last(2)+j) := Ab(i,j);
        res(rows,A'last(2)+j) := h(j);
      end loop;
    end loop;
    return res;
  end Straight_Differentiate_Deflation;

  function Straight_Evaluate_Deflation
              ( n1,n2 : integer32; f : Eval_Poly_Sys; jf : Eval_Jaco_Mat;
                cols : integer32; ajf : Array_of_Eval_Jaco_Mat;
                b1,b2 : Matrix; h1,h2,x,lambda1,lambda2 : Vector )
              return Vector is

  -- DESCRIPTION :
  --   Returns the value of the 2nd deflated system at x,
  --   with multipliers in lambda1 and lambda2.

  -- ON ENTRY :
  --   n1       number of equations in the 1st deflated system;
  --   n2       number of equations in the 2nd deflated system;
  --            the range of the vector on return is 1..n2;
  --   f        original polynomial system;
  --   jf       Jacobian matrix of the original system;
  --   cols     number of columns in the Jacobian matrix of 1st deflation; 
  --   ajf      all 1st derivatives of jf, the so-called "Hessian";
  --   b1       matrix for multiplying Jacobian matrix in 1st deflation;
  --   b2       matrix for multiplying Jacobian matrix in 2nd deflation;
  --   h1       for unique values for multipliers in 1st deflation;
  --   h2       for unique values for multipliers in 2nd deflation;
  --   x        values for the original variables;
  --   lambda1  values for the 1st set of multipliers;
  --   lambda2  values for the 2nd set of multipliers.

  -- ON RETURN :
  --   vector of range 1..n2 with the value of the deflated system
  --   at (x,lambda1,lambda2).

    res : Vector(1..n2);
    eva : constant Vector(1..n1)
        := Straight_Evaluate_Deflation(n1,f,jf,b1,h1,x,lambda1);
    mat : Matrix(1..n1,1..cols)
        := Straight_Differentiate_Deflation(n1,cols,jf,ajf,b1,h1,x,lambda1);
    b2l : Vector(b2'range(1)) := b2*lambda2;
    mbl : Vector(1..n1) := mat*b2l;

  begin
    res(1..n1) := eva;
    for i in mbl'range loop
      res(n1+i) := mbl(i);
    end loop;
    res(n2) := Create(-1.0);
    for i in lambda2'range loop
      res(n2) := res(n2) + h2(i)*lambda2(i);
    end loop;
    return res;
  end Straight_Evaluate_Deflation;

  procedure Random_Correctness_Test
              ( n,m : in integer32; s,f : in Eval_Poly_Sys;
                sjf,jf : in Eval_Jaco_Mat; ajf : Array_of_Eval_Jaco_Mat;
                b : in Matrix; h : in Vector ) is

  -- DESCRIPTION :
  --   Generates a random point and compares the results of the
  --   evaluation at the symbolic deflated system with the value
  --   computed by the straight-line program.
  --   The difference between the two should be of machine precision.

  -- ON ENTRY :
  --   n        number of variables in the original system;
  --   m        number of multipliers in the deflation;
  --   s        deflated polynomial system computed symbolically;
  --   f        original polynomial system;
  --   sjf      Jacobi matrix of the symbolically deflated system;
  --   jf       Jacobi matrix of the original system;
  --   ajf      derivatives of the Jacobi matrix of the original system;
  --   b        random matrix used to create s;
  --   h        random vector used to create s.

    x : Vector(1..n) := Random_Vector(1,n);
    lambda : Vector(1..m) := Random_Vector(1,m);
    x_lambda : Vector(1..n+m) := Stack(x,lambda);
    y1 : Vector(s'range) := Eval(s,x_lambda);
    A1 : Matrix(sjf'range(1),sjf'range(2)) := Eval(sjf,x_lambda);
    y2 : Vector(s'range)
       := Straight_Evaluate_Deflation(s'last,f,jf,b,h,x,lambda);
    A2 : Matrix(sjf'range(1),sjf'range(2))
       := Straight_Differentiate_Deflation(s'last,n+m,jf,ajf,b,h,x,lambda);
    dy : Vector(y1'range) := y1-y2;
    ndy : constant double_float := Max_Norm(dy);
    dA : Matrix(A1'range(1),A1'range(2)) := A1-A2;
    ndA : constant double_float := Max_Norm(dA);
    ans : character;

  begin
    put("Do you wish to see all the numbers ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("A random point : "); put_line(x);
      put_line("with random multipliers : "); put_line(lambda);
      put_line("Evaluation of symbolic deflated system :");
      put_line(y1);
      put_line("Result of straight-line evaluation :");
      put_line(y2);
    end if;
    put("Norm of the difference of evaluation : "); put(ndy,3); new_line;
    if ans = 'y' then
      put_line("Differentiation of symbolic deflated system : ");
      put(A1,3);
      put_line("Result of straight-line differentiation : ");
      put(A2,3);
    end if;
    put("Norm of the difference of differentation : "); put(ndA,3); new_line;
  end Random_Correctness_Test;

  procedure Random_Efficiency_Test
              ( k,n,m : in integer32; s,f : in Eval_Poly_Sys;
                sjf,jf : in Eval_Jaco_Mat; ajf : in Array_of_Eval_Jaco_Mat;
                b : in Matrix; h : in Vector ) is

  -- DESCRIPTION :
  --   Generates a random point and compares the efficiency of
  --   the two evaluation methods.

  -- ON ENTRY :
  --   k        number of evaluations to be run;
  --   n        number of variables in the original system;
  --   m        number of multipliers in the deflation;
  --   s        deflated polynomial system computed symbolically;
  --   f        original polynomial system;
  --   sjf      Jacobi matrix of the symbolically deflated system;
  --   jf       Jacobi matrix of the original system;
  --   ajf      derivatives of the Jacobi matrix of the original system;
  --   b        random matrix used to create s;
  --   h        random vector used to create s.

    x : constant Vector(1..n) := Random_Vector(1,n);
    lambda : constant Vector(1..m) := Random_Vector(1,m);
    x_lambda : constant Vector(1..n+m) := Stack(x,lambda);
    y1,y2 : Vector(s'range);
    A1,A2 : Matrix(s'range,1..n+m);
    timer : Timing_Widget;

  begin
    put_line("A random point : "); put_line(x);
    put_line("with random multipliers : "); put_line(lambda);
    new_line;
    put("Evaluating the symbolic system "); put(k,1);
    put_line(" times ...");
    new_line;
    tstart(timer);
    for i in 1..k loop
      y1 := Eval(s,x_lambda);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"symbolic system evaluation");
    new_line;
    put("Evaluating the straight deflation "); put(k,1);
    put_line(" times ...");
    new_line;
    tstart(timer);
    for i in 1..k loop
      y2 := Straight_Evaluate_Deflation(s'last,f,jf,b,h,x,lambda);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"straight-line deflation");
    new_line;
    put("Evaluating the differentiated symbolic system "); put(k,1);
    put_line(" times ...");
    new_line;
    tstart(timer);
    for i in 1..k loop
      A1 := Eval(sjf,x_lambda);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"symbolic system differentiation");
    new_line;
    put("Straight deflation "); put(k,1);
    put_line(" times ...");
    new_line;
    tstart(timer);
    for i in 1..k loop
      A2 := Straight_Differentiate_Deflation(s'last,n+m,jf,ajf,b,h,x,lambda);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"straight-line differentiation");
  end Random_Efficiency_Test;

  procedure Second_Straight_Deflate
              ( p : in Poly_Sys; n,m1,m2 : in integer32;
                b1 : in Matrix; h1 : in Vector;               
                f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                ajf : Array_of_Eval_Jaco_Mat ) is

  -- DESCRIPTION :
  --   Tests the straight-line program for deflating the system p,
  --   which is already a deflated polynomial system.

  -- ON ENTRY :
  --   p        deflated system, using m1 multipliers, b1 and h1;
  --   n        original number of variables;
  --   m1       number of multipliers used in the first deflation;
  --   m2       number of multipliers for the second deflation;
  --   b1       multiplier matrix used in the 1st deflation;
  --   h1       vector to have unique values for 1st multipliers;
  --   f        evaluable form for the original polynomial system;
  --   jf       Jacobian matrix of the original polynomial system;
  --   ajf      Hessian needed to evaluate Jacobian of 1st deflation.

    n1 : constant integer32 := n+m1;
    h2 : Vector(1..m2) := Random_Vector(1,m2);
    b2 : Matrix(1..n1,1..m2) := Random_Matrix(natural32(n1),natural32(m2));
    s : constant Poly_Sys := Deflate(p,natural32(m2),b2,h2);
    sf : Eval_Poly_Sys(s'range) := Create(s);
    x : Vector(1..n) := Random_Vector(1,n);
    l1 : Vector(1..m1) := Random_Vector(1,m1);
    l2 : Vector(1..m2) := Random_Vector(1,m2);
    xl1 : Vector(1..n+m1) := Stack(x,l1);
    xl1l2 : Vector(1..n+m1+m2) := Stack(xl1,l2);
    y1 : Vector(s'range) := Eval(sf,xl1l2);
    y2 : Vector(s'range) 
       := Straight_Evaluate_Deflation 
            (p'last,s'last,f,jf,n1,ajf,b1,b2,h1,h2,x,l1,l2);
    dy : Vector(s'range) := y1-y2;
    ndy : constant double_float := Max_Norm(dy);
    timer : timing_widget;
    k : integer32;

  begin
    put_line("The value at the symbolically deflated system : ");
    put_line(y1);
    put_line("The value of the straight-line program : ");
    put_line(y2);
    put("Max norm of the difference between the values : ");
    put(ndy,3); new_line;
    new_line;
    put("Give number of tests : "); get(k);
    put("evaluating symbolic 2nd deflated system ");
    put(k,1); put_line(" ...");
    new_line;
    tstart(timer);
    for i in 1..k loop
      y1 := Eval(sf,xl1l2);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"symbolic 2nd deflation");
    new_line;
    put("evaluating 2nd straight-line deflation ");
    put(k,1); put_line(" ...");
    new_line;
    tstart(timer);
    for i in 1..k loop
      y2 := Straight_Evaluate_Deflation 
              (p'last,s'last,f,jf,n1,ajf,b1,b2,h1,h2,x,l1,l2);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"2nd straight deflation");
  end Second_Straight_Deflate;

  procedure Straight_Deflate
              ( p : in Poly_Sys; n,m : in integer32 ) is 

  -- DESCRIPTION :
  --   Tests the straight-line programs for deflating the system p
  --   which has n original variables, using m multipliers.

    h : Vector(1..m) := Random_Vector(1,m);
    b : Matrix(1..n,1..m) := Random_Matrix(natural32(n),natural32(m));
    s : constant Poly_Sys := Deflate(p,natural32(m),b,h);
    sjm : Jaco_Mat(s'range,1..n+m) := Create(s);
    sjf : Eval_Jaco_Mat(sjm'range(1),sjm'range(2)) := Create(sjm);
    f : Eval_Poly_Sys(s'range) := Create(s);
    ep : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..n) := Create(p);
    jf : Eval_Jaco_Mat(jm'range(1),jm'range(2)) := Create(jm);
    ajm : Array_of_Jaco_Mat(1..n) := Create(jm);
    ajf : Array_of_Eval_Jaco_Mat(1..n) := Create(ajm);
    k,m2 : integer32 := 0;

  begin
    Random_Correctness_Test(n,m,f,ep,sjf,jf,ajf,b,h);
    put("Give number of evaluations to test efficiency : "); get(k); 
    if k > 0
     then Random_Efficiency_Test(k,n,m,f,ep,sjf,jf,ajf,b,h);
    end if;
    put("current number of variables : "); put(n+m,1); new_line;
    put("current number of equations : "); put(s'last,1); new_line;
    put("Give number of multipliers in second deflation : "); get(m2);
    if m2 > 0
     then Second_Straight_Deflate(s,n,m,m2,b,h,ep,jf,ajf);
    end if;
  end Straight_Deflate;

  procedure Main is

    lp : Link_to_Poly_Sys;
    m,n : integer32 := 0;

  begin
    new_line;
    put_line("Testing dedicated data structures for deflation.");
    get(lp);
    n := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("Give the number of multipliers : "); get(m);
    Straight_Deflate(lp.all,n,m);
  end Main;

begin
  Main;
end ts_straight;
