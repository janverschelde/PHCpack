with text_io,integer_io;                 use text_io,integer_io;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_VecVecs;           use Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Intrinsic_Witness_Sets_io;          use Intrinsic_Witness_Sets_io;
with Equation_by_Equation_Solvers;       use Equation_by_Equation_Solvers;

procedure testnine is

-- DESCRIPTION :
--   Test program for the evaluation and differentiation of one
--   equation arising in the nine-point problem.

  function Evaluate ( z,p : Vector ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the value at z for precision point p.

    x  : Complex_Number := z(1);
    xx : Complex_Number := z(2);
    y  : Complex_Number := z(3);
    yy : Complex_Number := z(4);
    a  : Complex_Number := z(5);
    aa : Complex_Number := z(6);
    b  : Complex_Number := z(7);
    bb : Complex_Number := z(8);
    d  : Complex_Number := p(1);
    dd : Complex_Number := p(2);
    c1 : Complex_Number := (aa-dd)*x;
    c2 : Complex_Number := (a-d)*xx;
    c3 : Complex_Number := d*(aa-xx) + dd*(a-x) - d*dd;
    d1 : Complex_Number := (bb-dd)*y;
    d2 : Complex_Number := (b-d)*yy;
    d3 : Complex_Number := d*(bb-yy) + dd*(b-y) - d*dd;
    g1 : Complex_Number := c2*d3-d2*c3;
    g2 : Complex_Number := c3*d1-d3*c1;
    g3 : Complex_Number := c1*d2-d1*c2;
    f : Complex_Number := g1*g2+g2*g3+g3*g1;

  begin
    return f;
  end Evaluate;

  function Evaluate ( z : Vector; p : VecVec ) return Vector is

  -- DESCRIPTION :
  --   Returns the vector of function values at z for all
  --   precision points in p.

    f : Vector(p'range);

  begin
    for k in p'range loop
      f(k) := Evaluate(z,p(k).all);
    end loop;
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
    a  : Complex_Number := z(5);
    aa : Complex_Number := z(6);
    b  : Complex_Number := z(7);
    bb : Complex_Number := z(8);
    d  : Complex_Number := p(1);
    dd : Complex_Number := p(2);
    c1 : Complex_Number := (aa-dd)*x;
    c2 : Complex_Number := (a-d)*xx;
    c3 : Complex_Number := d*(aa-xx) + dd*(a-x) - d*dd;
    d1 : Complex_Number := (bb-dd)*y;
    d2 : Complex_Number := (b-d)*yy;
    d3 : Complex_Number := d*(bb-yy) + dd*(b-y) - d*dd;
    g1 : Complex_Number := c2*d3-d2*c3;
    g2 : Complex_Number := c3*d1-d3*c1;
    g3 : Complex_Number := c1*d2-d1*c2;
    dg1,dg2,dg3,df : Vector(1..8);

  begin
    dg1(1) := dd*d2;
    dg1(2) := (a-d)*d3+d*d2;
    dg1(3) := -dd*c2;
    dg1(4) := -d*c2-(b-d)*c3;
    dg1(5) := xx*d3-dd*d2;
    dg1(6) := -d*d2;
    dg1(7) := dd*c2-yy*c3;
    dg1(8) := d*c2;
    dg2(1) := -dd*d1-(aa-dd)*d3; 
    dg2(2) := -d*d1; 
    dg2(3) := (bb-dd)*c3+dd*c1;
    dg2(4) := d*c1;
    dg2(5) := dd*d1;
    dg2(6) := d*d1-x*d3; 
    dg2(7) := -dd*c1; 
    dg2(8) := y*c3-d*c1; 
    dg3(1) := (aa-dd)*d2;
    dg3(2) := -(a-d)*d1;
    dg3(3) := -(bb-dd)*c2;
    dg3(4) := (b-d)*c1;
    dg3(5) := -xx*d1;
    dg3(6) := x*d2;
    dg3(7) := yy*c1;
    dg3(8) := -y*c2;
    df := dg1*(g2+g3) + dg2*(g1+g3) + dg3*(g2+g1);
    return df;
  end Differentiate;

  procedure Test_Evaluation_and_Differentiation is

    eps : constant double_float := 1.0E-7;
    jc : Vector(1..8) := (1..8 => Create(0.0));
    z,zz,df : Vector(1..8);
    p : Vector(1..2);
    r1,r2 : double_float;
    bC,bD,bP,t,th,CC,DD,PP,bA,bB,x,y,a,b,f,ff : Complex_Number;

  begin
   -- first test evaluation of the function :
    r1 := Random; r2 := Random; bC := Create(r1,r2);
    r1 := Random; r2 := Random; bD := Create(r1,r2);
    r1 := Random; r2 := Random; bP := Create(r1,r2);
    r1 := Random; r2 := Random; t := Create(r1,r2);
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
    z(5) := a; z(6) := Conjugate(a);
    z(7) := b; z(8) := Conjugate(b);
    p(1) := PP - bP; p(2) := Conjugate(p(1));
    put("function value : "); put(Evaluate(z,p)); new_line;
   -- test computation of the Jacobian matrix :
    z := Random_Vector(1,8);
    p := Random_Vector(1,2);
    f := Evaluate(z,p);
    df := Differentiate(z,p);
    for i in 1..8 loop
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

    res : Vector(1..12);

  begin
    res(x'range) := x;        -- any component equation to zero  
    res(9) := x(1) - x(3);    -- x = y
    res(10) := x(2) - x(4);   -- xx = yy
    res(11) := x(5) - x(7);   -- a = b
    res(12) := x(6) - x(8);   -- aa = bb 
    return res;
  end Q;

  procedure Read_Precision_Points ( p : out VecVec ) is

  -- DESCRIPTION :
  --   Prompts the user for the 8 precision points.

    file : file_type;

  begin
    new_line;
    put_line("Reading the file name for the precision points...");
    Read_Name_and_Open_File(file);
    for i in 1..8 loop
      declare
        point : Vector(1..2);
      begin
        get(file,point(1));
        get(file,point(2));
        p(i) := new Vector'(point);
      end;
    end loop;
    close(file);
  end Read_Precision_Points;

  procedure Solve_NinePoint_Problem is 

  -- DESCRIPTION :
  --   Generates random precision points and launches 
  --   the equation-by-equation solver on the ninepoint problem.

    file : file_type;
    name : Link_to_String;
    witset : Array_of_Solution_Lists(1..8);
    planes : VecMat(1..8);
    d : Standard_Natural_Vectors.Vector(1..8) := (1..8 => 7);
    p : VecVec(1..8);
    ans : character;
    want_stone,have_stone : boolean;
    s_dim,k : natural;
    s_sols : Solution_List;
    s_p : Link_to_Matrix;

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
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file,name);
    new_line;
    put("Do you have already a witness set for the first equations ? ");
    Ask_Yes_or_No(ans);
    have_stone := (ans = 'y');
    if have_stone
     then Read_Precision_Points(p);
          Read_Witness_Stone(s_dim,k,s_sols,s_p);
          planes(k) := s_p;
          witset(k) := s_sols;
     else k := 0;
          for i in p'range loop
            p(i) := new Vector'(Random_Vector(1,2));
          end loop;
    end if;
    new_line;
    put("Do you want stepping stones of indetermediate witness sets ? ");
    Ask_Yes_or_No(ans);
    want_stone := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    put_line(file,"The precision points :");
    for i in p'range loop
      put(file,"coordinates of precision point "); put(file,i,1);
      put_line(file," :"); put_line(file,p(i).all);
    end loop;
    Solve(file,name.all,want_stone,8,8,k,d,witset,planes);
  end Solve_NinePoint_Problem;

  procedure Evaluate_Residuals
              ( file : in file_type; p : in VecVec;
                sols : in Solution_List ) is

  -- DESCRIPTION :
  --   Evaluates all solutions for the given precision points
  --   and prints the residual for every solution.

    tmp : Solution_List := sols;
    x,y : Vector(p'range);
    nrm : double_float;

  begin
    for i in 1..Length_Of(sols) loop
      x := Head_Of(tmp).v;
      y := Evaluate(x,p); 
      nrm := Max_Norm(y);
      put(file,"Residual "); put(file,i,1);
      put(file," : "); put(file,nrm,3); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Evaluate_Residuals;

  procedure Read_Witness_Stone
              ( p : out VecVec; sols : out Solution_List ) is

  -- DESCRIPTION :
  --   Prompts the user for the 8 precision points and solutions.

    file : file_type;
    found : boolean;

  begin
    Read_Precision_Points(p);
    new_line;
    loop
      put_line("Reading the file name for the solutions...");
      Read_Name_and_Open_File(file);
      Scan_and_Skip(file,"SOLUTIONS",found);
      if not found
       then put_line("No SOLUTIONS on this file, correct format?");
            put_line("Please try again...");
       else get(file,sols);
      end if;
      close(file);
      exit when found;
    end loop;
  end Read_Witness_Stone;

  procedure Verify_Solution_List is

    file : file_type;
    p : VecVec(1..8);
    sols : Solution_List;

  begin
    Read_Witness_Stone(p,sols);
    new_line;
    put_line("Reading the name of the output file...");
    Read_Name_and_Create_File(file);
    put_line(file,"THE PRECISION POINTS :");
    put_line(file,p);
    if not Is_Null(sols)
     then new_line(file);
          put_line(file,"THE SOLUTIONS :");
          put(file,Length_Of(sols),Head_Of(sols).n,sols);
          Evaluate_Residuals(file,p,sols);
    end if;
  end Verify_Solution_List;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for the 9-point problem :");
    put_line("  1. Test evaluation and differentiation routines;");
    put_line("  2. Run equation-by-equation solver on the problem;");
    put_line("  3. Verify a list of solutions.");
    put("Type 1, 2, or 3 to select : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Test_Evaluation_and_Differentiation;
      when '2' => Solve_NinePoint_Problem;
      when '3' => Verify_Solution_List;
      when others => null;
    end case;
  end Main;

begin
  Main;
end testnine;
