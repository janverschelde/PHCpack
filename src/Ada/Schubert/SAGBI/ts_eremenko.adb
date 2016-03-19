with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Natural_Matrices;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;
with Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Integer32_Triangulations;  use Standard_Integer32_Triangulations;
with Standard_Dynamic32_Triangulations;  use Standard_Dynamic32_Triangulations;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Triangulations_and_Subdivisions;    use Triangulations_and_Subdivisions;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;
with Exponent_Vectors;                   use Exponent_Vectors;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Standard_Homotopy;
with Continuation_Parameters;
with Matrix_Indeterminates;
with SAGBI_Homotopies;                   use SAGBI_Homotopies;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Osculating_Planes;                  use Osculating_Planes;
with Complex_Osculating_Planes;          use Complex_Osculating_Planes;

procedure ts_eremenko is

-- DESCRIPTION :
--   Test on the set up of the variant of the Shapiro^2 conjectures,
--   communicated by Eremenko.

  procedure Test_Standard_Basis ( n,d : in natural32 ) is

  -- DESCRIPTION :
  --   This procedure compares the real with the complex generated
  --   matrix in the standard basis.

    s : double_float := 0.0;
    cs : Complex_Number;
    rm : Standard_Floating_Matrices.Matrix(1..integer32(n),1..integer32(d));
    cm : Standard_Complex_Matrices.Matrix(1..integer32(n),1..integer32(d));

  begin
    new_line;
    put("Give a real value for s : "); get(s);
    cs := Create(s);
    rm := Standard_Basis(n,d,s);
    cm := Standard_Basis(n,d,cs);
    put_line("The real matrix : "); put(rm,3);
    put_line("The real part of the complex matrix :");
    for i in cm'range(1) loop 
      for j in cm'range(2) loop
        put(REAL_PART(cm(i,j)),3);
      end loop;
      new_line;
    end loop;
  end Test_Standard_Basis;

  procedure Test_Solution ( file : in file_type; n,d : in natural32;
                            sol : in Solution ) is

    a1,a2,a3,b1,b2,b3,c1,c2,c3 : double_float;

  begin
    a1 := IMAG_PART(sol.v(1));
    put(file,"a1 : "); put(file,a1); new_line(file);
    a2 := REAL_PART(sol.v(2));
    put(file,"a2 : "); put(file,a2); new_line(file);
    a3 := IMAG_PART(sol.v(3));
    put(file,"a3 : "); put(file,a3); new_line(file);
    b1 := IMAG_PART(sol.v(4));
    put(file,"b1 : "); put(file,b1); new_line(file);
    b2 := REAL_PART(sol.v(5));
    put(file,"b2 : "); put(file,b2); new_line(file);
    b3 := IMAG_PART(sol.v(6));
    put(file,"b3 : "); put(file,b3); new_line(file);
    c1 := IMAG_PART(sol.v(7));
    put(file,"c1 : "); put(file,c1); new_line(file);
    c2 := REAL_PART(sol.v(8));
    put(file,"c2 : "); put(file,c2); new_line(file);
    c3 := IMAG_PART(sol.v(9));
    put(file,"c3 : "); put(file,c3); new_line(file);
    put(file,"a3 - a1*b2 : ");
    put(file,a3-a1*b2); new_line(file);
    put(file,"b3 - b1*c2 : ");
    put(file,b3-b1*c2); new_line(file);
  end Test_Solution;

  procedure Test_Solutions ( n,d : in natural32 ) is

    sols,tmp : Solution_List;
    file : file_type;

  begin
    Read(sols);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    tmp := sols;
    for i in 1..Length_Of(sols) loop
      put(file,"Testing solution "); put(file,i,1);
      put_line(file," :");
      Test_Solution(file,n,d,Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Solutions;

  function Imaginary_Poles return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector of nine complex numbers on the imaginary axis,
  --   including the origin, symmetrically about the real axis.

    res : Standard_Complex_Vectors.Vector(1..9);
    ranflt : double_float;

  begin
    res(5) := Create(0.0);
    for i in 1..integer32(4) loop
      ranflt := Random;
      res(i) := Create(0.0,ranflt);
      res(i+5) := -res(i);
    end loop;
    return res;
  end Imaginary_Poles;

  function Poles_on_Circle return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns all nine poles on a circle, symmetrically about the
  --   real axis, with the origin as one of the poles.

    res : Standard_Complex_Vectors.Vector(1..9);

  begin
    res(5) := Create(0.0);
    for i in 1..integer32(4) loop
      res(i) := Random1;
      res(i+5) := -res(i);
    end loop;
    return res;
  end Poles_on_Circle;

  function Localization_Map return Standard_Natural_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a matrix of range 1..6 times 1..3 which has the identity
  --   at the top and two's every where else.

    res : Standard_Natural_Matrices.Matrix(1..6,1..3);

  begin
    for j in 1..integer32(3) loop
      for i in 1..integer32(3) loop
        if i = j
         then res(i,j) := 1;
         else res(i,j) := 0;
        end if;
      end loop;
      for i in 4..integer32(6) loop
        res(i,j) := 2;
      end loop;
    end loop;
    return res;
  end Localization_Map;

  function Complex_Random_SAGBI_Homotopy
             ( locmap : Standard_Natural_Matrices.Matrix ) return Poly_Sys is

  -- DESCRIPTION :
  --   Generates a SAGBI homotopy to solve a complex random case.

    res : Poly_Sys(1..9);
    plane : Standard_Complex_Matrices.Matrix(1..6,1..3);
    p : Poly := Lifted_Localized_Laplace_Expansion(locmap);

  begin
    for i in res'range loop
      plane := Random_Orthogonal_Matrix(6,3);
      res(i) := Intersection_Condition(plane,p);
    end loop;
    return res;
  end Complex_Random_SAGBI_Homotopy;

  function Complex_Osculating_SAGBI_Homotopy
             ( locmap : Standard_Natural_Matrices.Matrix;
               poles : Standard_Complex_Vectors.Vector ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the SAGBI homotopy for the given localization map
  --   and choice of complex poles.

    res : Poly_Sys(1..9);
    plane : Standard_Complex_Matrices.Matrix(1..6,1..3);
    p : Poly := Lifted_Localized_Laplace_Expansion(locmap);

  begin
    for i in res'range loop
      plane := Standard_Basis(6,3,poles(i));
      res(i) := Intersection_Condition(plane,p);
    end loop;
    return res;
  end Complex_Osculating_SAGBI_Homotopy;

  procedure Flat_Deformation
               ( sh : in Poly_Sys; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Performs the flat deformation defined by the SAGBI homotopy.
  --   This is the second stage in the continuation.

    use Standard_Complex_Jaco_Matrices;

    sh_eval : Eval_Poly_Sys(sh'range);
    jac_mat : Jaco_Mat(sh'range,sh'first..sh'last+1);
    eva_jac : Eval_Jaco_Mat(jac_mat'range(1),jac_mat'range(2));

    function Eval ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      return Eval(sh_eval,xt);
    end Eval;

    function Diff ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Matrices.Matrix is

      res : Standard_Complex_Matrices.Matrix(x'range,x'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range(1) loop
        for j in res'range(2) loop
          res(i,j) := Eval(eva_jac(i,j),xt);
        end loop;
      end loop;
      return res;
    end Diff;

    function Diff ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                  return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(x'range);
      xt : Standard_Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range loop
        res(i) := Eval(eva_jac(i,eva_jac'last(2)),xt);
      end loop;
      return res;
    end Diff;

    procedure Sil_Cont is new Silent_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    sh_eval := Create(sh);
    jac_mat := Create(sh);
    eva_jac := Create(jac_mat);
    Set_Continuation_Parameter(sols,Create(0.0));
    Sil_Cont(sols,false,target=>Create(1.0));
    Clear(sh_eval);
    Clear(jac_mat);
    Clear(eva_jac);
  end Flat_Deformation;

  procedure Solve_Start_Configuration
              ( file : in file_type;
                locmap : in Standard_Natural_Matrices.Matrix;
                start : out Poly_Sys; sols : out Solution_List ) is

  -- DESCRIPTION :
  --   Applies the SAGBI homotopy to solve one generic complex instance.
  --   There are three stages to this solution process :
  --    1) dynamic lifting to the support of the start system;
  --    2) incremental polyhedral continuation to solve the start system;
  --    3) flat deformation to solve the random complex instance.

    use Standard_Complex_Laur_JacoMats;

    saghom : constant Poly_Sys(1..9) := Complex_Random_SAGBI_Homotopy(locmap);
    startsys : Poly_Sys(1..9);
    n : constant natural32 := Number_of_Unknowns(saghom(saghom'first));
    support,lifted,lifted_last : List;
    t : Triangulation;
    vol : natural32;
    lifted_lq,lq : Laur_Sys(startsys'range);
    mix : Standard_Integer_Vectors.Vector(1..1) := (1..1 => startsys'last);
    lif : Array_of_Lists(1..1);
    h : Eval_Coeff_Laur_Sys(startsys'range);
    c : Standard_Complex_VecVecs.VecVec(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));
    mixsub : Mixed_Subdivision;

  begin
    for i in startsys'range loop
      startsys(i) := Eval(saghom(i),Create(0.0),integer32(n));
      start(i) := Eval(saghom(i),Create(1.0),integer32(n));
    end loop;
    support := Create(startsys(startsys'first));
    Dynamic_Lifting(support,false,false,0,lifted,lifted_last,t);
    vol := Volume(t);
   -- put(file,"The volume is "); put(file,vol,1); new_line(file);
    lif(1..1) := (1..1 => lifted);
    mixsub := Shallow_Create(startsys'last,t);
    lq := Polynomial_to_Laurent_System(startsys);
    lifted_lq := Perform_Lifting(startsys'last,mix,lif,lq);
    h := Create(lq);
    for i in c'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        c(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(lq);
    Create(lq,j,m);
    Mixed_Solve(lifted_lq,lif,h,c,e,j,m,mix,mixsub,sols);
   -- put_line(file,"The solutions after polyhedral continuation : ");
   -- put(file,sols);
    Flat_Deformation(saghom,sols);
   -- put_line(file,"The solutions after flat deformation : ");
   -- put(file,sols);
  end Solve_Start_Configuration;

  function Is_Real ( sol : Solution; tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the imaginary part of every component of the solution
  --   vector is less than the tolerance tol.  Returns false otherwise.

  begin
    for i in sol.v'range loop
      if abs(IMAG_PART(sol.v(i))) > tol
       then return false;
      end if;
    end loop;
    return true;
  end Is_Real;

  function Real_Solutions ( sols : Solution_List ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of real solutions in the list.

    res : natural32 := 0;
    tol : constant double_float := 10.0**(-8);
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      if Is_Real(Head_Of(tmp).all,tol)
       then res := res+1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Real_Solutions;

  procedure Solve_Target_Configuration
              ( file : in file_type; 
                locmap : in Standard_Natural_Matrices.Matrix;
                s : in Standard_Complex_Vectors.Vector;
                start : in Standard_Complex_Poly_Systems.Poly_Sys;
                startsols : in Solution_List ) is

  -- DESCRIPTION :
  --   Applies the cheater's homotopy.

    saghom : Poly_Sys(1..9) := Complex_Osculating_SAGBI_Homotopy(locmap,s);
    target : Poly_Sys(1..9);
    n : constant natural32 := Number_of_Unknowns(saghom(saghom'first));
    a : constant Standard_Complex_Vectors.Vector(1..9) := Random_Vector(1,9);
    b : constant Standard_Complex_Vectors.Vector(1..9) := Random_Vector(1,9);
    sols : Solution_List;

  begin
    for i in target'range loop
      target(i) := Eval(saghom(i),Create(1.0),integer32(n));
    end loop;
    Standard_Homotopy.Create(target,start,2,a,b,true);       -- linear cheater
    Copy(startsols,sols);
    Set_Continuation_Parameter(sols,Create(0.0));
    Continuation_Parameters.Tune(2);
    declare
      procedure Sil_Cont is
        new Silent_Continue(Max_Norm,Standard_Homotopy.Eval,
                            Standard_Homotopy.Diff,Standard_Homotopy.Diff);
    begin
      Sil_Cont(sols,false,target=>Create(1.0));
    end;
   -- put_line(file,"The solution of the target system :");
   -- put(file,sols);
    put(file,"Number of real solutions : ");
    put(file,Real_Solutions(sols),1); new_line(file);
    Standard_Homotopy.Clear;
    Clear(saghom); Clear(target); Clear(sols);
  end Solve_Target_Configuration;

  procedure Massive_Random_Test is

  -- DESCRIPTION :
  --   For (m,p) = (3,3) we generate pole configurations like
  --    1) all on the imaginary axis, symmetrically wrt the real axis;
  --    2) 4 random complex points, origin, symmetrically wrt real axis;
  --    3) keeping 0, -i, +i fixed, symmetrically wrt the real axis; 
  --    4) on a circle with center on the real axis.

    file : file_type;
    cnt : natural32 := 0;
    s : Standard_Complex_Vectors.Vector(1..9);
    locmap : Standard_Natural_Matrices.Matrix(1..6,1..3) := Localization_Map;
    start : Standard_Complex_Poly_Systems.Poly_Sys(1..9);
    startsols : Solution_List;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put_line("Starting a very long loop.  See output file for results...");
    new_line;
    Matrix_Indeterminates.Initialize_Symbols(6,3);
    Matrix_Indeterminates.Reduce_Symbols(locmap);
    Solve_Start_Configuration(file,locmap,start,startsols);
    for i in 1..100000000 loop
      cnt := cnt + 1;
      put(file,"Configuration "); put(file,cnt,1); put_line(file," :");
     -- s := Imaginary_Poles;
      s := Poles_on_Circle;
      put_line(file,s);
      Solve_Target_Configuration(file,locmap,s,start,startsols);
    end loop;
  end Massive_Random_Test;

  procedure Main is

    n,d : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing variant of Shapiro^2 conjectures");
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. Generate and compare with the real case.");
    put_line("  2. Check the solutions.");
    put_line("  3. Perform massive random test for (m,p)=(3,3).");
    put("Type 1, 2 or 3 : "); Ask_Alternative(ans,"123");
    if (ans /= '3') then
      put("Give the dimension of the space : "); get(n);
      put("Give the dimension of the subspace : "); get(d);
    end if;
    case ans is
      when '1' => Test_Standard_Basis(n,d);
      when '2' => Test_Solutions(n,d);
      when '3' => Massive_Random_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_eremenko;
