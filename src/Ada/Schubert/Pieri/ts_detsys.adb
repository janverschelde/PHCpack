with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Matrix_Indeterminates;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;
with Brackets;                           use Brackets;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;
with Determinantal_Systems;              use Determinantal_Systems;

procedure ts_detsys is

-- DESCRIPTION :
--   This procedure tests the operations in the package Determinantal_Systems.

  function Random_Sequence ( n,n1,n2 : natural32 ) return VecMat is

  -- DESCRIPTION :
  --   Returns a sequence of n randomly generated n1-by-n2 matrices.

    res : VecMat(1..integer32(n));

  begin
    for i in 1..integer32(n) loop
      res(i) := new Standard_Complex_Matrices.Matrix'(Random_Matrix(n1,n2));
    end loop;
    return res;
  end Random_Sequence;

  function Vector_Rep ( mat : Standard_Complex_Matrices.Matrix )
                      return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the elements in the matrix as one long vector.

    res : Standard_Complex_Vectors.Vector(1..mat'length(1)*mat'length(2));
    cnt : integer32 := 0;

  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        cnt := cnt+1;
        res(cnt) := mat(i,j);
      end loop;
    end loop;
    return res;
  end Vector_Rep;

  procedure Test_Evaluator ( m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Test the evaluating routine in determinantal systems.

    n : constant natural32 := m+p;
    mp : constant natural32 := m*p;
    xmat : constant Standard_Complex_Matrices.Matrix := Random_Matrix(n,p);
    xvec : constant Standard_Complex_Vectors.Vector := Vector_Rep(xmat);
    L : constant Standard_Complex_Matrices.Matrix := Random_Matrix(n,m);
    planes : VecMat(1..integer32(mp)) := Random_Sequence(mp,n,m);
    deteva : Standard_Complex_Vectors.Vector(1..integer32(mp))
           := Eval(planes,xmat);
    detjac : Standard_Complex_Matrices.Matrix
               (1..integer32(mp),1..integer32(n*p))
           := Diff(planes,xmat);
    top : Bracket(1..integer32(p)) := (1..integer32(p) => 1);
    bottom : Bracket(1..integer32(p)) := (1..integer32(p) => n);
    xpm : Standard_Complex_Poly_Matrices.Matrix
            (1..integer32(n),1..integer32(p))
        := Localization_Pattern(n,top,bottom);
    sys : Poly_Sys(1..integer32(mp))
        := Polynomial_Equations(planes,xpm);
    syseva : Standard_Complex_Vectors.Vector(1..integer32(mp))
           := Eval(sys,xvec);
    sysjac : Jaco_Mat(1..integer32(mp),1..integer32(n*p)) := Create(sys);
    jaceva : Standard_Complex_Matrices.Matrix
               (1..integer32(mp),1..integer32(n*p))
           := Eval(sysjac,xvec);
    nb : natural32 := 0;
    timer : Timing_Widget;

  begin
    Matrix_Indeterminates.Initialize_Symbols(n,p);
    put_line("the matrix of indeterminates : "); put(xpm);
    put("Intersecting random "); put(p,1); put("-plane with ");
    put(m*p,1); put(" random "); put(m,1); put_line("-planes.");
    put_line("Determinantal Evaluation : "); put_line(deteva); new_line;
    put_line("Polynomial Evaluation : ");    put_line(syseva); new_line;
    put_line("Determinantal Differentation : "); put(detjac,2); new_line;
    put_line("Polynomial Differentation : ");    put(jaceva,2); new_line;
    put("Give number of evaluations : "); get(nb);
    tstart(timer);
    for i in 1..nb loop
      deteva := Eval(planes,xmat);
      detjac := Diff(planes,xmat);
    end loop;
    tstop(timer);
    print_times(Standard_Output,timer,"determinantal evaluations");
    tstart(timer);
    for i in 1..nb loop
      syseva := Eval(sys,xvec);
      jaceva := Eval(sysjac,xvec);
    end loop;
    tstop(timer);
    print_times(Standard_Output,timer,"polynomial evaluations");
  end Test_Evaluator;

  procedure Main is

    m,p : natural32 := 0;

  begin
    put("Give m : "); get(m);
    put("Give p : "); get(p);
    Test_Evaluator(m,p);
  end Main;

begin
  new_line;
  put_line("Polynomial systems generated from determinantal expansions.");
  new_line;
  Main;
end ts_detsys;
