with text_io,integer_io;                 use text_io,integer_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Equation_by_Equation_Solvers;       use Equation_by_Equation_Solvers;

procedure testfivehom1 is

-- DESCRIPTION :
--   Use nine-point problem eqs. to solve five point problem with 
--   specified ground pivots.
--   Test program for the evaluation and differentiation of one
--   equation arising in the nine-point problem.
--   One-homogenized version

  function Evaluate ( z,p : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value at z for precision point p.

    x  : Complex_Number := z(1);
    xx : Complex_Number := z(2);
    y  : Complex_Number := z(3);
    yy : Complex_Number := z(4);
    mu : Complex_Number := z(5);
    a  : Complex_Number := p(1);
    aa : Complex_Number := p(2);
    b  : Complex_Number := p(3);
    bb : Complex_Number := p(4);
    d  : Complex_Number := p(5);
    dd : Complex_Number := p(6);
    c1 : Complex_Number := (aa-dd)*x;
    c2 : Complex_Number := ( a- d)*xx;
    c3 : Complex_Number := d*(aa*mu-xx) + dd*(a*mu-x) - d*dd*mu;
    d1 : Complex_Number := (bb-dd)*y;
    d2 : Complex_Number := ( b-d)*yy;
    d3 : Complex_Number := d*(bb*mu-yy) + dd*(b*mu-y) - d*dd*mu;
    g1 : Complex_Number := c2*d3-d2*c3;
    g2 : Complex_Number := c3*d1-d3*c1;
    g3 : Complex_Number := c1*d2-d1*c2;
    f  : Complex_Number := g1*g2+g2*g3+g3*g1;

  begin
    return f;
  end Evaluate;

  function Differentiate ( z,p : Vector ) return Vector is

  -- DESCRIPTION :
  --   Returns the row of the Jacobian matrix evaluated at z
  --   for precision point p.

    x  : Complex_Number := z(1);
    xx : Complex_Number := z(2);
    y  : Complex_Number := z(3);
    yy : Complex_Number := z(4);
    mu : Complex_Number := z(5);
    a  : Complex_Number := p(1);
    aa : Complex_Number := p(2);
    b  : Complex_Number := p(3);
    bb : Complex_Number := p(4);
    d  : Complex_Number := p(5);
    dd : Complex_Number := p(6);
    c1 : Complex_Number := (aa-dd)*x;
    c2 : Complex_Number := ( a- d)*xx;
    c3 : Complex_Number := d*(aa*mu-xx) + dd*(a*mu-x) - d*dd*mu;
    c3mu : Complex_Number := d*aa + dd*a - d*dd;
    d1 : Complex_Number := (bb-dd)*y;
    d2 : Complex_Number := ( b-d)*yy;
    d3 : Complex_Number := d*(bb*mu-yy) + dd*(b*mu-y) - d*dd*mu;
    d3mu : Complex_Number := d*bb + dd*b - d*dd;
    g1 : Complex_Number := c2*d3-d2*c3;
    g2 : Complex_Number := c3*d1-d3*c1;
    g3 : Complex_Number := c1*d2-d1*c2;
    dg1,dg2,dg3,df : Vector(1..5);

  begin
    dg1(1) := dd*d2;
    dg1(2) := (a-d)*d3+d*d2;
    dg1(3) := -dd*c2;
    dg1(4) := -d*c2-(b-d)*c3;
    dg1(5) := c2*d3mu - d2*c3mu;
    dg2(1) := -dd*d1-(aa-dd)*d3; 
    dg2(2) := -d*d1; 
    dg2(3) := (bb-dd)*c3+dd*c1;
    dg2(4) := d*c1;
    dg2(5) := c3mu*d1-d3mu*c1; 
    dg3(1) := (aa-dd)*d2;
    dg3(2) := -(a-d)*d1;
    dg3(3) := -(bb-dd)*c2;
    dg3(4) := (b-d)*c1;
    dg3(5) := Create(0.0);
    df := dg1*(g2+g3) + dg2*(g1+g3) + dg3*(g2+g1);
    return df;
  end Differentiate;

  procedure Test_Evaluation_and_Differentiation is

    eps : constant double_float := 1.0E-7;
    jc : Vector(1..5) := (1..5 => Create(0.0));
    z,zz,df : Vector(1..5);
    p : Vector(1..6);
    r1,r2,mu : double_float;
    bC,bD,bP,t,th,CC,DD,PP,bA,bB,x,y,a,b,f,ff : Complex_Number;

  begin
   -- first test evaluation of the function :
    r1 := Random; r2 := Random; bC := Create(r1,r2);
    r1 := Random; r2 := Random; bD := Create(r1,r2);
    r1 := Random; r2 := Random; bP := Create(r1,r2);
    r1 := Random; r2 := Random; t := Create(r1,r2);
    mu := Random;
    th := Random1;
    CC := t + th*bC;
    DD := t + th*bD;
    PP := t + th*bP;
    r1 := Random; r2 := Random;
    bA := (DD+bD)/2.0 + Create(0.0,r1)*(DD-bD);
    bB := (CC+bC)/2.0 + Create(0.0,r2)*(CC-bC);
    x := bD - bP; y := bC - bP; a := bA - bP; b := bB - bP;
    z(1) := x; z(2) := Conjugate(x);
    z(3) := y; z(4) := Conjugate(y);
    z(5) := Create(1.0);
    z := Create(mu)*z;
    p(1) := a; p(2) := Conjugate(a);
    p(3) := b; p(4) := Conjugate(b);
    p(5) := PP - bP; p(6) := Conjugate(p(5));
    put("function value : "); put(Evaluate(z,p)); new_line;
   -- test computation of the Jacobian matrix :
    z := Random_Vector(1,5);
    p := Random_Vector(1,6);
    f := Evaluate(z,p);
    df := Differentiate(z,p);
    for i in 1..5 loop
      zz := z;
      zz(i) := z(i) + eps;
      ff := Evaluate(zz,p);
      jc(i) := (ff - f)/eps;
    end loop;
    put_line("error on the Jacobian matrix :");
    put_line(df-jc);
  end Test_Evaluation_and_Differentiation;

  function Q ( x : Vector ) return Vector is

  -- DESCRIPTION :
  --   This function evaluation a solution at discriminating equations.

    res : Vector(1..7);

  begin
    res(x'range) := x;        -- any component equation to zero  
    res(6) := x(1) - x(3);    -- x = y
    res(7) := x(2) - x(4);   -- xx = yy
    return res;
  end Q;

  procedure Solve_NinePoint_Problem is 

  -- DESCRIPTION :
  --   Generates random precision points and launches 
  --   the equation-by-equation solver on the ninepoint problem.

    file : file_type;
    name : Link_to_String;
    witset : Array_of_Solution_Lists(1..4);
    planes : VecMat(1..4);
    d : Standard_Natural_Vectors.Vector(1..4) := (1..4 => 4);
    p : VecVec(1..4);
    fixed_p : Vector(1..6);

    function f ( k : natural; x : Vector ) return Complex_Number is
    begin
      return Evaluate(x,p(k).all);
    end f;

    function jf ( k : natural; x : Vector ) return Vector is
    begin
      return Differentiate(x,p(k).all);
    end jf;

    procedure Solve is new GQ_Solve_Equation_by_Equation(f,jf,Q);

  begin
    fixed_p := Random_Vector(1,4);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file,name);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    for i in p'range loop
      p(i) := new Vector'(Random_Vector(1,6));
      p(i)(1) := fixed_p(1);
      p(i)(2) := fixed_p(2);
      p(i)(3) := fixed_p(3);
      p(i)(4) := fixed_p(4);
    end loop;
    put_line(file,"The randomly generated precision points :");
    for i in p'range loop
      put(file,"coordinates of precision point "); put(file,i,1);
      put_line(file," :"); put_line(file,p(i).all);
    end loop;
    -- calling is:
    -- solve(output_file_ptr,outputfile_name,Num_eqs,Num_vars,degrees,witnesspts,slicing_planes);
    Solve(file,name.all,4,5,d,witset,planes);
  end Solve_NinePoint_Problem;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for the 5-point problem :");
    put_line("  1. Test evaluation and differentiation routines;");
    put_line("  2. Run equation-by-equation solver on the problem.");
    put("Type 1 or 2 to select : "); Ask_Alternative(ans,"12");
    if ans = '1'
     then Test_Evaluation_and_Differentiation;
     else Solve_NinePoint_Problem;
    end if;
  end Main;

begin
  Main;
end testfivehom1;
