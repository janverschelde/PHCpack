with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with Quad_Double_Two_Norms;              use Quad_Double_Two_Norms;
with Quad_Double_Linear_Solvers;         use Quad_Double_Linear_Solvers;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Quad_Double_Eigenvalues;            use Quad_Double_Eigenvalues;
with QuadDobl_Quad_Parameters;
with QuadDobl_Quad_Turn_Points_io;       use QuadDobl_Quad_Turn_Points_io;

package body QuadDobl_Quad_Turn_Points is

-- I. STEP SIZE CONTROL and CORRECTOR METHODS :

  procedure Set_Step_Size
               ( h : in out quad_double; flag : in integer32 ) is

    reduction : constant quad_double
              := QuadDobl_Quad_Parameters.reduction_multiplier;
    expansion : constant quad_double
              := QuadDobl_Quad_Parameters.expansion_multiplier;
    threshold : constant natural32
              := QuadDobl_Quad_Parameters.expansion_threshold;
    max_step : constant quad_double := QuadDobl_Quad_Parameters.max_step_size;

  begin
    if flag < 0 then
      h := max_step;
    elsif flag = 0 then
      h := h*reduction;
    elsif flag > integer32(threshold) then
      h := h*expansion;
      if h > max_step
       then h := max_step;
      end if;
    end if;
  end Set_Step_Size;

  procedure Step_Size_Control
               ( h : in out quad_double; flag : in integer32 ) is
  begin
   -- put("flag : "); put(flag,1); put(" step :"); put(h,3);
    Set_Step_Size(h,flag);
   -- put("  ->"); put(h,3); new_line;
  end Step_Size_Control;

  procedure Step_Size_Control
               ( file : in file_type;
                 h : in out quad_double; flag : in integer32 ) is
  begin
    put(file,"flag : "); put(file,flag,1); 
    put(file," step :"); put(file,h,3);
    Set_Step_Size(h,flag);
    put(file,"  ->"); put(file,h,3); new_line(file);
  end Step_Size_Control;

  procedure Step_Size_Control 
               ( h : in out Complex_Number; flag : in integer32 ) is

    rh : quad_double := REAL_PART(h);

  begin
    Step_Size_Control(rh,flag);
    h := Create(rh);
  end Step_Size_Control;

  procedure Step_Size_Control 
               ( file : in file_type;
                 h : in out Complex_Number; flag : in integer32 ) is

    rh : quad_double := REAL_PART(h);

  begin
    Step_Size_Control(file,rh,flag);
    h := Create(rh);
  end Step_Size_Control;

  procedure One_Corrector_Step
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                t : in Quad_Double_Vectors.Vector;
                x,y : in out Quad_Double_Vectors.Vector;
                err,res : out quad_double ) is

    m : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(2))
      := Quad_Double_Jaco_Matrices.Eval(jm,x);
    A : Quad_Double_Matrices.Matrix(x'range,x'range);
    b : Quad_Double_Vectors.Vector(x'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;

    use Quad_Double_Vectors;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -y(i);
    end loop;
    for j in t'range loop
      A(A'last(1),j) := t(j);
    end loop;
    b(b'last(1)) := create(0.0);
    lufac(A,A'last(1),piv,info);
    lusolve(A,A'last(1),piv,b);
    Quad_Double_Vectors.Add(x,b);
    y := Quad_Double_Poly_SysFun.Eval(p,x);
    err := Quad_Double_Two_Norms.Norm2(b);
    res := Quad_Double_Two_Norms.Norm2(y);
  end One_Corrector_Step;

  procedure One_Corrector_Step
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                t : in QuadDobl_Complex_Vectors.Vector;
                x,y : in out QuadDobl_Complex_Vectors.Vector;
                err,res : out quad_double ) is

    m : constant QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
      := QuadDobl_Complex_Jaco_Matrices.Eval(jm,x);
    A : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);
    b : QuadDobl_Complex_Vectors.Vector(x'range);
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;

    use QuadDobl_Complex_Vectors;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -y(i);
    end loop;
    for j in t'range loop
      A(A'last(1),j) := t(j);
    end loop;
    b(b'last(1)) := Create(integer(0));
    lufac(A,A'last(1),piv,info);
    lusolve(A,A'last(1),piv,b);
    QuadDobl_Complex_Vectors.Add(x,b);
    y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    err := Norm2(b);
    res := Norm2(y);
  end One_Corrector_Step;

  procedure One_Corrector_Step
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                x,y : in out Quad_Double_Vectors.Vector;
                A : out Quad_Double_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                err,res : out quad_double ) is

  -- DESCRIPTION :
  --   Auxiliary internal routine to the exported One_Corrector_Step
  --   procedure that returns the LU factored Jacobian matrix to
  --   compute the determinant off, along with the pivoting info.

    m : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(2))
      := Quad_Double_Jaco_Matrices.Eval(jm,x);
    b : Quad_Double_Vectors.Vector(x'range);
    info : integer32;

    use Quad_Double_Vectors;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -y(i);
    end loop;
    for j in A'first(2)..(A'last(2)-1) loop
      A(A'last(1),j) := create(integer(0));
    end loop;
    A(A'last(1),A'last(2)) := create(integer(1));
    b(b'last(1)) := create(integer(0));
    lufac(A,A'last(1),piv,info);
    lusolve(A,A'last(1),piv,b);
    Quad_Double_Vectors.Add(x,b);
    y := Quad_Double_Poly_SysFun.Eval(p,x);
    err := Quad_Double_Two_Norms.Norm2(b);
    res := Quad_Double_Two_Norms.Norm2(y);
  end One_Corrector_Step;

  procedure One_Corrector_Step
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                x,y : in out Quad_Double_Vectors.Vector;
                err,res : out quad_double ) is

    A : Quad_Double_Matrices.Matrix(x'range,x'range);
    piv : Standard_Integer_Vectors.Vector(x'range);

  begin
    One_Corrector_Step(p,jm,x,y,A,piv,err,res);
  end One_Corrector_Step;

  procedure One_Corrector_Step
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                x,y : in out Quad_Double_Vectors.Vector;
                err,res,det : out quad_double ) is

    A : Quad_Double_Matrices.Matrix(x'range,x'range);
    piv : Standard_Integer_Vectors.Vector(x'range);

  begin
    One_Corrector_Step(p,jm,x,y,A,piv,err,res);
    det := Determinant_after_LU(A,piv);
  end One_Corrector_Step;

  procedure One_Corrector_Step
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                x,y : in out QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix;
                piv : out Standard_Integer_Vectors.Vector;
                err,res : out quad_double ) is

  -- DESCRIPTION :
  --   Internal auxiliary routine to the exported One_Corrector_Step
  --   procedure to return in addition the LU factored Jacobian matrix at x,
  --   along with the pivoting information for an eventual determinant.

    m : constant QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
      := QuadDobl_Complex_Jaco_Matrices.Eval(jm,x);
    b : QuadDobl_Complex_Vectors.Vector(x'range);
    info : integer32;

    use QuadDobl_Complex_Vectors;

  begin
    for i in jm'range(1) loop
      for j in jm'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -y(i);
    end loop;
    for j in A'first(2)..(A'last(2)-1) loop
      A(A'last(1),j) := Create(integer(0));
    end loop;
    A(A'last(1),A'last(2)) := Create(integer(1));
    b(b'last(1)) := Create(integer(0));
    lufac(A,A'last(1),piv,info);
    lusolve(A,A'last(1),piv,b);
    QuadDobl_Complex_Vectors.Add(x,b);
    y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    err := Norm2(b);
    res := Norm2(y);
  end One_Corrector_Step;

  procedure One_Corrector_Step
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                x,y : in out QuadDobl_Complex_Vectors.Vector;
                err,res : out quad_double ) is

    A : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);
    piv : Standard_Integer_Vectors.Vector(x'range);

  begin
    One_Corrector_Step(p,jm,x,y,A,piv,err,res);
  end One_Corrector_Step;

  procedure One_Corrector_Step
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                x,y : in out QuadDobl_Complex_Vectors.Vector;
                err,res,det : out quad_double ) is

    A : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);
    piv : Standard_Integer_Vectors.Vector(x'range);

  begin
    One_Corrector_Step(p,jm,x,y,A,piv,err,res);
    det := AbsVal(Determinant_after_LU(A,piv));
  end One_Corrector_Step;

  procedure Interactive_Correct_Solution
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                t : in Quad_Double_Vectors.Vector;
                x : in out Quad_Double_Vectors.Vector ) is

    y : Quad_Double_Vectors.Vector(p'range);
    err,res : quad_double;
    ans : character;

  begin
    put_line("correcting the solution ...");
    y := Quad_Double_Poly_SysFun.Eval(p,x);
    loop
      One_Corrector_Step(p,jm,t,x,y,err,res);
      Write_Corrector_Information(x,y,err,res);
      put("Do one more correction step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    put_line("The current solution "); Write_Vector(x);
  end Interactive_Correct_Solution;

  procedure Interactive_Correct_Solution
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                t : in QuadDobl_Complex_Vectors.Vector;
                x : in out QuadDobl_Complex_Vectors.Vector ) is

    y : QuadDobl_Complex_Vectors.Vector(p'range);
    err,res : quad_double;
    ans : character;

  begin
    put_line("correcting the solution ...");
    y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    loop
      One_Corrector_Step(p,jm,t,x,y,err,res);
      Write_Corrector_Information(x,err,res);
      put("Do one more correction step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    put_line("The current solution "); Write_Vector(x);
  end Interactive_Correct_Solution;

  procedure Write_Corrector_Diagnostics
              ( file : in file_type;
                fail,converge : in boolean; nbrit : in natural32 ) is

  -- DESCRIPTION :
  --   Writes corrector diagnostics.

  begin
    if fail then
      put(file,"  failure : ");
      if not converge
       then put_line(file,"corrector is diverging");
       else put_line(file,"corrector needs more steps");
      end if;
    else
      put(file,"  success : reached accuracy in ");
      put(file,nbrit,1); put_line(file," steps");
    end if;
  end Write_Corrector_Diagnostics;

  procedure Correct_Solution
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                t : in Quad_Double_Vectors.Vector;
                x : in out Quad_Double_Vectors.Vector;
                tol_err,tol_res : in quad_double; fail : out boolean;
                nbrit : out natural32; maxit : in natural32 ) is

    y : Quad_Double_Vectors.Vector(p'range);
    err,prev_err,res,prev_res : quad_double;
    converge : boolean := true;

  begin
    y := Quad_Double_Poly_SysFun.Eval(p,x);
    fail := true; nbrit := 0;
    while converge and fail and (nbrit < maxit) loop
      One_Corrector_Step(p,jm,t,x,y,err,res);
      if ((err < tol_err) and (res < tol_res)) then
        fail := false;
      else
        if nbrit = 0 then
          prev_err := err; prev_res := res;
        else
          if ((err > prev_err) or (res > prev_res)) then
            converge := false;
          end if;
        end if;
      end if;
      nbrit := nbrit + 1;
    end loop;
  end Correct_Solution;

  procedure Correct_Solution
              ( file : in file_type;
                p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                t : in Quad_Double_Vectors.Vector;
                x : in out Quad_Double_Vectors.Vector;
                tol_err,tol_res : in quad_double; fail : out boolean;
                nbrit : out natural32; maxit : in natural32 ) is

    y : Quad_Double_Vectors.Vector(p'range);
    err,prev_err,res,prev_res : quad_double;
    converge : boolean := true;

  begin
    put_line(file,"correcting the solution ...");
    y := Quad_Double_Poly_SysFun.Eval(p,x);
    fail := true; nbrit := 0;
    while converge and fail and (nbrit < maxit) loop
      One_Corrector_Step(p,jm,t,x,y,err,res);
      Write_Corrector_Information(file,x,y,err,res);
      if ((err < tol_err) and (res < tol_res)) then
        fail := false;
      else
        if nbrit = 0 then
          prev_err := err; prev_res := res;
        else
          if ((err > prev_err) or (res > prev_res)) then
            converge := false;
          end if;
        end if;
      end if;
      nbrit := nbrit + 1;
    end loop;
    put_line(file,"The current solution "); Write_Vector(file,x);
    Write_Corrector_Diagnostics(file,fail,converge,nbrit);
  end Correct_Solution;

  procedure Correct_Solution
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                t : in QuadDobl_Complex_Vectors.Vector;
                x : in out QuadDobl_Complex_Vectors.Vector;
                tol_err,tol_res : in quad_double; fail : out boolean;
                nbrit : out natural32; maxit : in natural32 ) is

    y : QuadDobl_Complex_Vectors.Vector(p'range);
    err,prev_err,res,prev_res : quad_double;
    converge : boolean := true;

  begin
    y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    fail := true; nbrit := 0;
    while converge and fail and (nbrit < maxit) loop
      One_Corrector_Step(p,jm,t,x,y,err,res);
      if ((err < tol_err) and (res < tol_res)) then
        fail := false;
      else
        if nbrit = 0 then
          prev_err := err; prev_res := res;
        else
          if ((err > prev_err) or (res > prev_res))
           then converge := false;
          end if;
        end if;
      end if;
      nbrit := nbrit + 1;
    end loop;
  end Correct_Solution;

  procedure Correct_Solution
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                t : in QuadDobl_Complex_Vectors.Vector;
                x : in out QuadDobl_Complex_Vectors.Vector;
                tol_err,tol_res : in quad_double; fail : out boolean;
                nbrit : out natural32; maxit : in natural32 ) is

    y : QuadDobl_Complex_Vectors.Vector(p'range);
    err,prev_err,res,prev_res : quad_double;
    converge : boolean := true;

  begin
    put_line(file,"correcting the solution ...");
    y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    fail := true; nbrit := 0;
    while converge and fail and (nbrit < maxit) loop
      One_Corrector_Step(p,jm,t,x,y,err,res);
      Write_Corrector_Information(file,x,err,res);
      if ((err < tol_err) and (res < tol_res)) then
        fail := false;
      else
        if nbrit = 0 then
          prev_err := err; prev_res := res;
        else
          if ((err > prev_err) or (res > prev_res))
           then converge := false;
          end if;
        end if;
      end if;
      nbrit := nbrit + 1;
    end loop;
    put_line(file,"The current solution "); Write_Vector(file,x);
    Write_Corrector_Diagnostics(file,fail,converge,nbrit);
  end Correct_Solution;

  procedure Target_Correction
               ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                 jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                 target : in quad_double;
                 x : in out Quad_Double_Vectors.Vector;
                 tol_err,tol_res : in quad_double; fail : out boolean;
                 nbrit : out natural32; maxit : in natural32 ) is

     y : Quad_Double_Vectors.Vector(p'range);
     err,res,det : quad_double;

  begin
    fail := true; nbrit := 0;
    x(x'last) := target;
    y := Quad_Double_Poly_SysFun.Eval(p,x);
    for i in 1..maxit loop
      One_Corrector_Step(p,jm,x,y,err,res,det);
      if (err < tol_err) and (res < tol_res)
       then fail := false; nbrit := i;
      end if;
      exit when not fail;
    end loop;
  end Target_Correction;

  procedure Target_Correction
               ( file : in file_type;
                 p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                 jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                 target : in quad_double;
                 x : in out Quad_Double_Vectors.Vector;
                 tol_err,tol_res : in quad_double; fail : out boolean;
                 nbrit : out natural32; maxit : in natural32 ) is

     y : Quad_Double_Vectors.Vector(p'range);
     err,res,det : quad_double;

  begin
    put_line(file,"correcting solution back to target...");
    fail := true; nbrit := 0;
    x(x'last) := target;
    y := Quad_Double_Poly_SysFun.Eval(p,x);
    for i in 1..maxit loop
      One_Corrector_Step(p,jm,x,y,err,res,det);
      Write_Corrector_Information(file,x,y,err,res,det);
      if (err < tol_err) and (res < tol_res)
       then fail := false; nbrit := i;
      end if;
      exit when not fail;
    end loop;
  end Target_Correction;

  procedure Target_Correction
               ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                 target : in quad_double;
                 x : in out QuadDobl_Complex_Vectors.Vector;
                 tol_err,tol_res : in quad_double; fail : out boolean;
                 nbrit : out natural32; maxit : in natural32 ) is

     y : QuadDobl_Complex_Vectors.Vector(p'range);
     err,res,det : quad_double;

  begin
    fail := true; nbrit := 0;
    x(x'last) := Create(target);
    y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    for i in 1..maxit loop
      One_Corrector_Step(p,jm,x,y,err,res,det);
      if (err < tol_err) and (res < tol_res)
       then fail := false; nbrit := i;
      end if;
      exit when not fail;
    end loop;
  end Target_Correction;

  procedure Target_Correction
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                 target : in quad_double;
                 x : in out QuadDobl_Complex_Vectors.Vector;
                 tol_err,tol_res : in quad_double; fail : out boolean;
                 nbrit : out natural32; maxit : in natural32 ) is

     y : QuadDobl_Complex_Vectors.Vector(p'range);
     err,res,det : quad_double;

  begin
    put_line(file,"correcting solution back to target...");
    fail := true; nbrit := 0;
    x(x'last) := Create(target);
    y := QuadDobl_Complex_Poly_SysFun.Eval(p,x);
    for i in 1..maxit loop
      One_Corrector_Step(p,jm,x,y,err,res,det);
      Write_Corrector_Information(file,x,err,res,det);
      if (err < tol_err) and (res < tol_res)
       then fail := false; nbrit := i;
      end if;
      exit when not fail;
    end loop;
  end Target_Correction;

-- II. COMPUTING TANGENTS, DETERMINANTS, and OTHER MONITORING INFO :

  function Inner_Product ( x,y : QuadDobl_Complex_Vectors.Vector )
                         return Complex_Number is

    res : Complex_Number := Create(integer(0));

  begin
    for i in x'range loop
      res := res + x(i)*Conjugate(y(i));
    end loop;
    return res;
  end Inner_Product;

  function Tangent ( jm : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                     x : Quad_Double_Vectors.Vector )
                   return Quad_Double_Vectors.Vector is

    t : Quad_Double_Vectors.Vector(x'range);
    m : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(2))
      := Quad_Double_Jaco_Matrices.Eval(jm,x);
    A : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(1));
    b : Quad_Double_Vectors.Vector(jm'range(1));
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -m(i,m'last(2));
    end loop;
    lufac(A,A'last(1),piv,info);
    lusolve(A,A'last(1),piv,b);
    t(b'range) := b;
    t(t'last) := create(1.0);
    Quad_Double_Two_Norms.Normalize(t);
    return t;
  end Tangent;

  function Tangent ( jm : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                     x : QuadDobl_Complex_Vectors.Vector )
                   return QuadDobl_Complex_Vectors.Vector is

    t : QuadDobl_Complex_Vectors.Vector(x'range);
    m : constant QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
      := QuadDobl_Complex_Jaco_Matrices.Eval(jm,x);
    A : QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(1));
    b : QuadDobl_Complex_Vectors.Vector(jm'range(1));
    piv : Standard_Integer_Vectors.Vector(b'range);
    nrm : quad_double;
    info : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -m(i,m'last(2));
    end loop;
    lufac(A,A'last(1),piv,info);
    lusolve(A,A'last(1),piv,b);
    t(b'range) := b;
    t(t'last) := Create(integer(1));
    nrm := Norm2(t);
    QuadDobl_Complex_Vectors.Mul(t,Create(1.0/nrm));
    return t;
  end Tangent;

  function Determinant_after_LU
              ( A : Quad_Double_Matrices.Matrix;
                piv : Standard_Integer_Vectors.Vector )
              return quad_double is

    res : quad_double := create(1.0);

  begin
    for i in A'range(1) loop
      res := res*A(i,i);
    end loop;
    for i in piv'range loop
      if piv(i) > i
       then res := -res;
      end if;
    end loop;
    return res;
  end Determinant_after_LU;

  function Determinant_after_LU
              ( A : QuadDobl_Complex_Matrices.Matrix;
                piv : Standard_Integer_Vectors.Vector )
              return Complex_Number is

    res : Complex_Number := Create(integer(1));

  begin
    for i in A'range(1) loop
      res := res*A(i,i);
    end loop;
    for i in piv'range loop
      if piv(i) > i
       then res := -res;
      end if;
    end loop;
    return res;
  end Determinant_after_LU;

  procedure Tangent_and_Determinant
              ( jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                x : in Quad_Double_Vectors.Vector;
                t : out Quad_Double_Vectors.Vector;
                d : out quad_double ) is

    m : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(2))
      := Quad_Double_Jaco_Matrices.Eval(jm,x);
    A : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(1));
    b : Quad_Double_Vectors.Vector(jm'range(1));
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -m(i,m'last(2));
    end loop;
    lufac(A,A'last(1),piv,info);
    if info /= 0 then
      d := create(0.0);
    else
      d := Determinant_after_LU(A,piv);
      lusolve(A,A'last(1),piv,b);
      t(b'range) := b;
      t(t'last) := create(integer(1));
      Quad_Double_Two_Norms.Normalize(t);
    end if;
  end Tangent_and_Determinant;

  procedure Tangent_and_Determinant
              ( jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in QuadDobl_Complex_Vectors.Vector;
                t : out QuadDobl_Complex_Vectors.Vector;
                d : out Complex_Number ) is

    m : constant QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(2))
      := QuadDobl_Complex_Jaco_Matrices.Eval(jm,x);
    A : QuadDobl_Complex_Matrices.Matrix(jm'range(1),jm'range(1));
    b : QuadDobl_Complex_Vectors.Vector(jm'range(1));
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;
    nrm : quad_double;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -m(i,m'last(2));
    end loop;
    lufac(A,A'last(1),piv,info);
    if info /= 0 then
      d := Create(integer(0));
    else
      d := Determinant_after_LU(A,piv);
      lusolve(A,A'last(1),piv,b);
      t(b'range) := b;
      t(t'last) := Create(integer(1));
      nrm := Norm2(t);
      QuadDobl_Complex_Vectors.Mul(t,Create(1.0/nrm));
    end if;
  end Tangent_and_Determinant;

  function Maximal_Minors
              ( evjm : Quad_Double_Matrices.Matrix )
              return Quad_Double_Vectors.Vector is

    n : constant integer32 := evjm'last(1) - 1;
    res : Quad_Double_Vectors.Vector(1..n+1);
    sA : Quad_Double_Matrices.Matrix(1..n,1..n);
    sp : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    for i in sA'range(1) loop
      for j in sA'range(2) loop
        sA(i,j) := evjm(i,j);
      end loop;
    end loop;
    lufac(sA,n,sp,info);
    res(n+1) := Determinant_after_LU(sA,sp);
    for k in 1..n loop       -- skip k-th column of evjm
      for j in sA'range(2) loop
        if j < k then
          for i in sA'range(1) loop
            sA(i,j) := evjm(i,j);
          end loop;
        else
          for i in sA'range(1) loop
            sA(i,j) := evjm(i,j+1);
          end loop;
        end if;
      end loop;
      lufac(sA,n,sp,info);
      res(k) := Determinant_after_LU(sa,sp);
    end loop;
    return res;
  end Maximal_Minors;

  procedure Eigenvalues ( evjm : in Quad_Double_Matrices.Matrix;
                          L : out QuadDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := evjm'last(1) - 1;
    sA : Quad_Double_Matrices.Matrix(1..n,1..n);
    ierr : integer32;

  begin
    for i in sA'range(1) loop
      for j in sA'range(2) loop
        sA(i,j) := evjm(i,j);
      end loop;
    end loop;
    Eigenvalues(sA,ierr,L);
  end Eigenvalues;

  procedure Eigenvectors ( evjm : in Quad_Double_Matrices.Matrix;
                            L : out QuadDobl_Complex_Vectors.Vector;
                            v : out QuadDobl_Complex_VecVecs.VecVec ) is

    n : constant integer32 := evjm'last(1) - 1;
    sA : Quad_Double_Matrices.Matrix(1..n,1..n);
    ierr : integer32;

  begin
    for i in sA'range(1) loop
      for j in sA'range(2) loop
        sA(i,j) := evjm(i,j);
      end loop;
    end loop;
    Eigenvectors(sA,ierr,L,v);
  end Eigenvectors;

  procedure Tangent_and_Minors
              ( jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                x : in Quad_Double_Vectors.Vector;
                t,d : out Quad_Double_Vectors.Vector ) is

    m : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(2))
      := Quad_Double_Jaco_Matrices.Eval(jm,x);
    A : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(1));
    b : Quad_Double_Vectors.Vector(jm'range(1));
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -m(i,m'last(2));
    end loop;
    lufac(A,A'last(1),piv,info);
    if info /= 0 then
      for i in d'range loop
        d(i) := create(0.0);
      end loop;
    else
      lusolve(A,A'last(1),piv,b);
      t(b'range) := b;
      t(t'last) := create(1.0);
      Quad_Double_Two_Norms.Normalize(t);
      d := Maximal_Minors(m);
    end if;
  end Tangent_and_Minors;

  procedure Tangent_Minors_and_Eigenvalues
              ( jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                x : in Quad_Double_Vectors.Vector;
                t,d : out Quad_Double_Vectors.Vector;
                L : out QuadDobl_Complex_Vectors.Vector ) is

    m : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(2))
      := Quad_Double_Jaco_Matrices.Eval(jm,x);
    A : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(1));
    b : Quad_Double_Vectors.Vector(jm'range(1));
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -m(i,m'last(2));
    end loop;
    lufac(A,A'last(1),piv,info);
    if info /= 0 then
      for i in d'range loop
        d(i) := create(0.0);
      end loop;
    else
      lusolve(A,A'last(1),piv,b);
      t(b'range) := b;
      t(t'last) := create(1.0);
      Quad_Double_Two_Norms.Normalize(t);
      d := Maximal_Minors(m);
      Eigenvalues(m,L);
    end if;
  end Tangent_Minors_and_Eigenvalues;

  procedure Tangent_Minors_and_Eigenvectors
              ( jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                x : in Quad_Double_Vectors.Vector;
                t,d : out Quad_Double_Vectors.Vector;
                L : out QuadDobl_Complex_Vectors.Vector;
                v : out QuadDobl_Complex_VecVecs.VecVec ) is

    m : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(2))
      := Quad_Double_Jaco_Matrices.Eval(jm,x);
    A : Quad_Double_Matrices.Matrix(jm'range(1),jm'range(1));
    b : Quad_Double_Vectors.Vector(jm'range(1));
    piv : Standard_Integer_Vectors.Vector(b'range);
    info : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        A(i,j) := m(i,j);
      end loop;
      b(i) := -m(i,m'last(2));
    end loop;
    lufac(A,A'last(1),piv,info);
    if info /= 0 then
      for i in d'range loop
        d(i) := create(0.0);
      end loop;
    else
      lusolve(A,A'last(1),piv,b);
      t(b'range) := b;
      t(t'last) := create(1.0);
      Quad_Double_Two_Norms.Normalize(t);
      d := Maximal_Minors(m);
      Eigenvectors(m,L,v);
    end if;
  end Tangent_Minors_and_Eigenvectors;

  procedure Report_Minors_and_Eigenvectors
              ( file : in file_type;
                m : in Quad_Double_Vectors.Vector;
                L : in QuadDobl_Complex_Vectors.Vector;
                v : in QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    put(file,"Minors :");
    for k in m'range loop
      put(file," "); put(file,m(k),1);
    end loop; 
    new_line(file);
    put_line(file,"The eigenvalues are :"); put_line(file,L);
    put_line(file,"The corresponding eigenvectors : "); put(file,v);
  end Report_Minors_and_Eigenvectors;

-- III. COMPUTING QUADRATIC TURNING POINTS :

  procedure Seek_Turn
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                x1,t1,x2,t2 : in out Quad_Double_Vectors.Vector;
                step : in quad_double ) is

    x3,t3 : Quad_Double_Vectors.vector(x1'range);
    det : Quad_Double_Vectors.Vector(p'range);
    h : quad_double := step;
    ans : character;
    nb : natural32;
    fail : boolean;

    use Quad_Double_Vectors;

  begin
    put_line("The solution and its tangent before the turn :");
    Write_Vector_and_its_Tangent(x1,t1);
    put_line("The solution and its tangent after the turn :");
    Write_Vector_and_its_Tangent(x2,t2);
    loop
      h := h/2.0;
      x3 := x1 + h*t1;
     -- Correct_Solution(p,jm,t1,x3);
      Correct_Solution(p,jm,t1,x3,create(1.0E-10),create(1.0E-10),fail,nb,4);
      Tangent_and_Minors(jm,x3,t3,det);
      put_line("The solution and its tangent after bisection :");
      Write_Vector_Tangent_and_Determinants(x3,t3,det);
      put("Continue with bisection ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      if t1*t3 < 0.0
       then x2 := x3; t2 := t3;
       else x1 := x3; t1 := t3;
      end if;
    end loop;
  end Seek_Turn;

  procedure Seek_Turn
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                x1,t1,x2,t2 : in out QuadDobl_Complex_Vectors.Vector;
                step : in Complex_Number ) is

    x3,t3 : QuadDobl_Complex_Vectors.vector(x1'range);
    h : Complex_Number := step;
    orientation : Complex_Number;
    ans : character;
    nb : natural32;
    fail : boolean;
    two : constant quad_double := create(2.0);

    use QuadDobl_Complex_Vectors;

  begin
    put_line("The solution and its tangent before the turn :");
    Write_Vector_and_its_Tangent(x1,t1);
    put_line("The solution and its tangent after the turn :");
    Write_Vector_and_its_Tangent(x2,t2);
    loop
      h := h/two;
      x3 := x1 + h*t1;
     -- Correct_Solution(p,jm,t1,x3);
      Correct_Solution(p,jm,t1,x3,create(1.0E-10),create(1.0E-10),fail,nb,4);
      t3 := Tangent(jm,x3);
      put_line("The solution and its tangent after bisection :");
      Write_Vector_and_its_Tangent(x3,t3);
      put("Continue with bisection ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      orientation := Inner_Product(t1,t3);
      if REAL_PART(orientation) < 0.0
       then x2 := x3; t2 := t3;
       else x1 := x3; t1 := t3;
      end if;
    end loop;
  end Seek_Turn;

  procedure Interactive_Shoot_Turn
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                x1,t1,x2,t2 : in out Quad_Double_Vectors.Vector;
                step : in quad_double ) is

    x3,t3 : Quad_Double_Vectors.vector(x1'range);
    det : Quad_Double_Vectors.Vector(p'range);
    h : quad_double := step;
    ans : character;
    nb : natural32;
    fail : boolean;

    use Quad_Double_Vectors;

  begin
    put_line("The solution and its tangent before the turn :");
    Write_Vector_and_its_Tangent(x1,t1);
    put_line("The solution and its tangent after the turn :");
    Write_Vector_and_its_Tangent(x2,t2);
    loop
      h := h*(t1(t1'last))/(t1(t1'last) + t2(t2'last));
      put("new step size h = "); put(h,3); new_line;
      x3 := x1 + h*t1;
      Correct_Solution
        (Standard_Output,p,jm,t1,x3,create(1.0E-10),create(1.0E-10),fail,nb,4);
      Tangent_and_Minors(jm,x3,t3,det);
      put_line("The solution and its tangent after shooting :");
      Write_Vector_Tangent_and_Determinants(x3,t3,det);
      if t1*t3 < 0.0 then
        x2 := x3; t2 := t3;  put_line("replaced x2");
      else
        x1 := x3; t1 := t3;  put("replaced x1, ");
        h := (x1-x2)*t2;
        put("adjusted step size h = "); put(h,3); new_line;
      end if;
      put("Continue with shooting ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Shoot_Turn;

  procedure Interactive_Shoot_Turn
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                x1,t1,x2,t2 : in out QuadDobl_Complex_Vectors.Vector;
                step : in Complex_Number ) is

    x3,t3 : QuadDobl_Complex_Vectors.vector(x1'range);
    h : Complex_Number := step;
    orientation : Complex_Number;
    ans : character;
    nb : natural32;
    fail : boolean;

    use QuadDobl_Complex_Vectors;

  begin
    put_line("The solution and its tangent before the turn :");
    Write_Vector_and_its_Tangent(x1,t1);
    put_line("The solution and its tangent after the turn :");
    Write_Vector_and_its_Tangent(x2,t2);
    loop
      h := h*(t1(t1'last))/(t1(t1'last) + t2(t2'last));
      put("new step size h = "); put(h,3); new_line;
      x3 := x1 + h*t1;
      Correct_Solution
        (Standard_Output,p,jm,t1,x3,create(1.0E-10),create(1.0E-10),fail,nb,4);
      t3 := Tangent(jm,x3);
      put_line("The solution and its tangent after shooting :");
      Write_Vector_and_its_Tangent(x3,t3);
      orientation := Inner_Product(t1,t3);
      if REAL_PART(orientation) < 0.0 then
        x2 := x3; t2 := t3; put_line("replaced x2");
      else
        x1 := x3; t1 := t3; put("replaced x1, ");
        h := Inner_Product((x1-x2),t2);
        put("adjusted step size h = "); put(h,3); new_line;
      end if;
      put("Continue with shooting ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Shoot_Turn;

  procedure Shoot_Turn
              ( p : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jm : in Quad_Double_Jaco_Matrices.Eval_Jaco_mat;
                x1,t1,x2,t2 : in out Quad_Double_Vectors.Vector;
                step,tol_step : in quad_double; max : in natural ) is

    x3,t3 : Quad_Double_Vectors.vector(x1'range);
    h : quad_double := step;
    nb : natural32;
    fail : boolean;

    use Quad_Double_Vectors;

  begin
    for i in 1..max loop
      h := h*(t1(t1'last))/(t1(t1'last) + t2(t2'last));
      x3 := x1 + h*t1;
      Correct_Solution(p,jm,t1,x3,create(1.0E-10),create(1.0E-10),fail,nb,4);
      t3 := Tangent(jm,x3);
      if t1*t3 < 0.0
       then x2 := x3; t2 := t3;
       else x1 := x3; t1 := t3; h := (x1-x2)*t2;
      end if;
      exit when (h < tol_step);
    end loop;
  end Shoot_Turn;

  procedure Shoot_Turn
              ( p : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jm : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_mat;
                x1,t1,x2,t2 : in out QuadDobl_Complex_Vectors.Vector;
                step,tol_step : in Complex_Number; max : in natural ) is

    x3,t3 : QuadDobl_Complex_Vectors.vector(x1'range);
    h : Complex_Number := step;
    orientation : Complex_Number;
    nb : natural32;
    fail : boolean;

    use QuadDobl_Complex_Vectors;

  begin
    for i in 1..max loop
      h := h*(t1(t1'last))/(t1(t1'last) + t2(t2'last));
      x3 := x1 + h*t1;
      Correct_Solution(p,jm,t1,x3,create(1.0E-10),create(1.0E-10),fail,nb,4);
      t3 := Tangent(jm,x3);
      orientation := Inner_Product(t1,t3);
      if REAL_PART(orientation) < 0.0
       then x2 := x3; t2 := t3;
       else x1 := x3; t1 := t3; h := Inner_Product((x1-x2),t2);
      end if;
      exit when (h < tol_step);
    end loop;
  end Shoot_Turn;

-- IV. PARABOLIC MINIMIZATION OF DETERMINANTS :

  procedure Quadratic_Interpolation
               ( file : in file_type;
                 x,y : in Quad_Double_Vectors.Vector;
                 p,q : out quad_double ) is
  begin
    p := x(1)*x(1)*(y(2)-y(3))
       + x(2)*x(2)*(y(3)-y(1))
       + x(3)*x(3)*(y(1)-y(2));
    q := 2.0*( y(3)*(x(2)-x(1)) + y(2)*(x(1)-x(3)) + y(1)*(x(3)-x(2)) );
  exception
    when others =>
       put_line(file,"Exception happened in Quadratic_Interpolation");
       put(file,"  p = "); put(file,p); new_line(file);
       put(file,"  q = "); put(file,q); new_line(file);
       raise;
  end Quadratic_Interpolation;

  procedure Quadratic_Interpolation
               ( x,y : in Quad_Double_Vectors.Vector;
                 p,q : out quad_double ) is

  begin
    Quadratic_Interpolation(standard_output,x,y,p,q);
  end Quadratic_Interpolation;

  procedure Silent_Monitor_Determinants
               ( x,y : in out Quad_Double_Vectors.Vector;
                 i : in out integer32; t,d : in quad_double;
                 crit : out natural32; z : out quad_double ) is

    p,q : quad_double;

  begin
    if i < x'last then
      i := i + 1;
    else -- i = t'last, shift window
      x(1) := x(2); x(2) := x(3);
      y(1) := y(2); y(2) := y(3);
    end if;
    x(i) := t; y(i) := d;
    if i < x'last then
      if i < x'last - 1 then
        crit := 0;
      elsif (y(1)*y(2) < 0.0) then
        crit := 3;
      else
        crit := 0;
      end if;
    else
      if y(2)*y(3) < 0.0 then
        crit := 3;
        z := (x(2) + x(3))/2.0;
      else
        Quadratic_Interpolation(x,y,p,q); z := p/q;
        if (z >= x(1) and z <= x(3))
         then crit := 4;
         else crit := 0;
        end if;
      end if;
    end if;
  end Silent_Monitor_Determinants;

  procedure Monitor_Determinants
               ( x,y : in out Quad_Double_Vectors.Vector;
                 i : in out integer32; t,d : in quad_double;
                 crit : out natural32; z : out quad_double ) is

    p,q : quad_double;

  begin
    if i < x'last then
      i := i + 1;
    else -- i = t'last, shift window
      x(1) := x(2); x(2) := x(3);
      y(1) := y(2); y(2) := y(3);
    end if;
    x(i) := t; y(i) := d;
    if i < x'last then
      if i < x'last - 1 then
        crit := 0;
      elsif (y(1)*y(2) < 0.0) then
        crit := 3;
      else
        crit := 0;
      end if;
    else
      if y(2)*y(3) < 0.0 then
        crit := 3;
        z := (x(2) + x(3))/2.0;
      else
        Quadratic_Interpolation(x,y,p,q); z := p/q;
        put("t values : "); put(x,3); new_line;
        put("d values : "); put(y,3); new_line;
        put("z = "); put(z,3);
        if (z >= x(1) and z <= x(3))
         then crit := 4;
         else crit := 0;
        end if;
      end if;
      case crit is
        when 3 => put_line(" Determinant sign flipped!  critical");
        when 4 => put_line(" Parabolic minimum inside!  critical");
        when others => put_line("  normal");
      end case;
    end if;
  end Monitor_Determinants;

  procedure Monitor_Determinants
               ( file : in file_type; 
                 x,y : in out Quad_Double_Vectors.Vector;
                 i : in out integer32; t,d : in quad_double;
                 crit : out natural32; z : out quad_double ) is

    p,q : quad_double;

  begin
    if i < x'last then
      i := i + 1;
    else -- i = t'last, shift window
      x(1) := x(2); x(2) := x(3);
      y(1) := y(2); y(2) := y(3);
    end if;
    x(i) := t;
    y(i) := d;
    if i < x'last then
      if i < x'last-1 then
        crit := 0;
      elsif (y(1)*y(2) < 0.0) then
        crit := 3;
      else
        crit := 0;
      end if;
    else
      if y(2)*y(3) < 0.0 then
        crit := 3;
        z := (x(2)+x(3))/2.0;
      else
        Quadratic_Interpolation(file,x,y,p,q); z := p/q;
        put(file,"t values : "); put(file,x,3); new_line(file);
        put(file,"d values : "); put(file,y,3); new_line(file);
        put(file,"z = "); put(file,z,3);
        if (z >= x(1) and z <= x(3))
         then crit := 4;
         else crit := 0;
        end if;
      end if;
      case crit is
        when 3 => put_line(file," Determinant sign flipped!  critical");
        when 4 => put_line(file," Parabolic minimum inside!  critical");
        when others => put_line(file,"  normal");
      end case;
    end if;
  exception
    when others => put_line(file,"exception happened in Monitor: critical");
                   put(file,"  p = "); put(file,p); new_line(file);
                   put(file,"  q = "); put(file,q); new_line(file);
                   crit := 4;
  end Monitor_Determinants;

  procedure Silent_Bisection_Singularity
                 ( t1,t2,d1,d2 : in out quad_double;
                   f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                   jf : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                   x : in out Quad_Double_Vectors.Vector;
                   tol_err,tol_res,tol_det : in quad_double;
                   max : in natural32; fail,critical : out boolean;
                   nit : out natural32 ) is

     y : Quad_Double_Vectors.Vector(f'range);
     err,res,det : quad_double;
     mt : quad_double;

  begin
    critical := true; nit := 0;
    for i in 1..10*max loop
      mt := (t1 + t2)/2.0; x(x'last) := mt;
      y := Quad_Double_Poly_SysFun.Eval(f,x);
      fail := true;
      for i in 1..max loop
        One_Corrector_Step(f,jf,x,y,err,res,det);
        if (err < tol_err) and (res < tol_res)
         then fail := false; nit := i;
        end if;
        exit when not fail;
      end loop;
      if d1*det < 0.0 then
        t2 := mt; d2 := det;
      else
        t1 := mt; d1 := det;
      end if;
    end loop;
  end Silent_Bisection_Singularity;

  procedure Bisection_Singularity
                 ( file : in file_type;
                   t1,t2,d1,d2 : in out quad_double;
                   f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                   jf : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                   x : in out Quad_Double_Vectors.Vector;
                   tol_err,tol_res,tol_det : in quad_double;
                   max : in natural32; fail,critical : out boolean;
                   nit : out natural32 ) is

     y : Quad_Double_Vectors.Vector(f'range);
     err,res,det : quad_double;
     mt : quad_double;

  begin
    put(file,"Bisection to Singularity d1 = ");
    put(file,d1,3); put(file,"  d2 = "); put(file,d2,3); new_line(file);
    put(file,"  t1 = "); put(file,t1);
    put(file,"  t2 = "); put(file,t2); new_line(file);
    critical := true; nit := 0;
    for i in 1..10*max loop
      put(file,"stage "); put(file,i,1); put_line(file," of bisection :");
      mt := (t1 + t2)/2.0; x(x'last) := mt;
      put(file,"new value for t : "); put(file,mt); new_line(file);
      y := Quad_Double_Poly_SysFun.Eval(f,x);
      fail := true;
      for i in 1..max loop
        One_Corrector_Step(f,jf,x,y,err,res,det);
        Write_Corrector_Information(file,x,y,err,res,det);
        if (err < tol_err) and (res < tol_res)
         then fail := false; nit := i;
        end if;
        exit when not fail;
      end loop;
      if d1*det < 0.0 then
        t2 := mt; d2 := det;
      else
        t1 := mt; d1 := det;
      end if;
    end loop;
  end Bisection_Singularity;

  procedure Silent_Parabolic_Minimization
                 ( vt,dt : in Quad_Double_Vectors.Vector;
                   zt : in quad_double;
                   f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                   jf : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                   x : in out Quad_Double_Vectors.Vector;
                   tol_err,tol_res,tol_det : in quad_double;
                   max : in natural32; fail,critical : out boolean;
                   nit : out natural32 ) is

    y : Quad_Double_Vectors.Vector(f'range);
    err,res,det,p,q,z : quad_double;
    nvt : Quad_Double_Vectors.Vector(vt'range) := vt;
    ndt : Quad_Double_Vectors.Vector(vt'range) := dt;

  begin
    z := zt;
    for i in 1..max loop
      y := Quad_Double_Poly_SysFun.Eval(f,x);
      fail := true; nit := max;
      for i in 1..max loop
        One_Corrector_Step(f,jf,x,y,err,res,det);
        if (err < tol_err) and (res < tol_res)
         then fail := false; nit := i;
        end if;
        exit when not fail;
      end loop;
      if z < nvt(2) then -- drop (nvt(3),ndt(3))
        nvt(3) := nvt(2); nvt(2) := z;  
        ndt(3) := ndt(2); ndt(2) := det; 
      else -- drop (nvt(1),ndt(1))
        nvt(1) := nvt(2); nvt(2) := z;
        ndt(1) := ndt(2); ndt(2) := det;
      end if;
      Quadratic_Interpolation(nvt,ndt,p,q); z := p/q;
      if z < nvt(1) or z > nvt(3) then
        exit; -- not critical because outside the interval
      end if;
      x(x'last) := z;
    end loop;
    critical := (abs(det) < tol_det);
  end Silent_Parabolic_Minimization;

  procedure Parabolic_Minimization
                 ( file : in file_type;
                   vt,dt : in Quad_Double_Vectors.Vector;
                   zt : in quad_double;
                   f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                   jf : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                   x : in out Quad_Double_Vectors.Vector;
                   tol_err,tol_res,tol_det : in quad_double;
                   max : in natural32; fail,critical : out boolean;
                   nit : out natural32 ) is

    y : Quad_Double_Vectors.Vector(f'range);
    err,res,det,p,q,z : quad_double;
    nvt : Quad_Double_Vectors.Vector(vt'range) := vt;
    ndt : Quad_Double_Vectors.Vector(vt'range) := dt;

  begin
    put(file,"Parabolic Minimization starts at t = ");
    put(file,zt); put_line(file," :");
    z := zt;
    for i in 1..max loop
      put(file,"stage "); put(file,i,1); put_line(file," of minimization :");
      y := Quad_Double_Poly_SysFun.Eval(f,x);
      fail := true; nit := max;
      for i in 1..max loop
        One_Corrector_Step(f,jf,x,y,err,res,det);
        Write_Corrector_Information(file,x,y,err,res,det);
        if (err < tol_err) and (res < tol_res)
         then fail := false; nit := i;
        end if;
        exit when not fail;
      end loop;
      if z < nvt(2) then -- drop (nvt(3),ndt(3))
        nvt(3) := nvt(2); nvt(2) := z;  
        ndt(3) := ndt(2); ndt(2) := det; 
      else -- drop (nvt(1),ndt(1))
        nvt(1) := nvt(2); nvt(2) := z;
        ndt(1) := ndt(2); ndt(2) := det;
      end if;
      put(file,"t values : "); put(file,nvt,5); new_line(file);
      put(file,"d values : "); put(file,ndt,5); new_line(file);
      Quadratic_Interpolation(file,nvt,ndt,p,q); z := p/q;
      put(file,"new critical value for t : "); put(file,z);
      if z < nvt(1) or z > nvt(3) then
        put_line(file," outside: not critical"); exit;
      else
        put_line(file," inside: stay critical");
      end if;
      x(x'last) := z;
      if abs(det) < 1.0E-14 then
        put_line(file,"determinant reached machine precision level"); exit;
      end if;
    end loop;
    critical := (abs(det) < tol_det);
  exception
    when others => put_line(file,"exception happened in Back");
                   put(file,"  p = "); put(file,p); new_line(file);
                   put(file,"  q = "); put(file,q); new_line(file);
  end Parabolic_Minimization;

  procedure Monitor_Singularity
               ( file : in file_type; output : in boolean;
                 nd : in quad_double;
                 vt,dt : in out Quad_Double_Vectors.Vector;
                 i : in out integer32;
                 f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                 jf : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out Quad_Double_Vectors.Vector;
                 px,pt,dx : in Quad_Double_Vectors.Vector;
                 tol_err,tol_res,tol_det : in quad_double;
                 max : in natural32; fail : out boolean;
                 nit,crtp : out natural32 ) is

    orientation,zt,t1,t2,d1,d2 : quad_double;
    critical : boolean;

    use Quad_Double_Vectors;

  begin
    if output then
      put(file,"Determinant : "); put(file,nd,3);
      put_line(file,"  Tangent vector is"); Write_Tangent(file,dx);
    end if;
    crtp := 0;
    if (abs(nd) < tol_det) then
      crtp := 1;
    else
      orientation := pt*dx;
      if output then
        put(file,"orientation of tangent : ");
        put(file,orientation,2); new_line(file);
      end if;
      if (orientation < 0.0) then
        crtp := 2;
      else
        Monitor_Determinants(file,vt,dt,i,x(x'last),nd,crtp,zt);
        if crtp = 3 then
          if i < vt'last then
            t1 := vt(1); t2 := vt(2); d1 := dt(1); d2 := dt(2);
          else
            t1 := vt(2); t2 := vt(3); d1 := dt(2); d2 := dt(3);
          end if;
          x := px; x(x'last) := t1;
          Bisection_Singularity(file,t1,t2,d1,d2,f,jf,x,
               tol_err,tol_res,tol_det,max,fail,critical,nit);
        elsif crtp = 4 then
          x := px; x(x'last) := zt;
          Parabolic_Minimization(file,vt,dt,zt,f,jf,x,
               tol_err,tol_res,tol_det,max,fail,critical,nit);
        end if;
        if not critical then
          crtp := 0;
        end if;
      end if;
    end if;
  end Monitor_Singularity;

  procedure Silent_Monitor_Singularity
               ( nd : in quad_double;
                 vt,dt : in out Quad_Double_Vectors.Vector;
                 i : in out integer32;
                 f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                 jf : Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out Quad_Double_Vectors.Vector;
                 px,pt,dx : in Quad_Double_Vectors.Vector;
                 tol_err,tol_res,tol_det : in quad_double;
                 max : in natural32; fail : out boolean;
                 nit,crtp : out natural32 ) is

    orientation,zt,t1,t2,d1,d2 : quad_double;
    critical : boolean;

    use Quad_Double_Vectors;

  begin
    crtp := 0;
    if (abs(nd) < tol_det) then
      crtp := 1;
    else
      orientation := pt*dx;
      if (orientation < 0.0) then
        crtp := 2;
      else
        Silent_Monitor_Determinants(vt,dt,i,x(x'last),nd,crtp,zt);
        if crtp = 3 then
          if i < vt'last then
            t1 := vt(1); t2 := vt(2); d1 := dt(1); d2 := dt(2);
          else
            t1 := vt(2); t2 := vt(3); d1 := dt(2); d2 := dt(3);
          end if;
          x := px; x(x'last) := t1;
          Silent_Bisection_Singularity(t1,t2,d1,d2,f,jf,x,
            tol_err,tol_res,tol_det,max,fail,critical,nit);
        elsif crtp = 4 then
          x := px; x(x'last) := zt;
          Silent_Parabolic_Minimization(vt,dt,zt,f,jf,x,
            tol_err,tol_res,tol_det,max,fail,critical,nit);
        end if;
        if not critical then
          crtp := 0;
        end if;
      end if;
    end if;
  end Silent_Monitor_Singularity;

  procedure Monitor_Singularity
               ( file : in file_type; output : in boolean;
                 nd : in Complex_Number;
                 vt,dt : in out Quad_Double_Vectors.Vector;
                 i : in out integer32;
                 f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out QuadDobl_Complex_Vectors.Vector;
                 px,pt,dx : in QuadDobl_Complex_Vectors.Vector;
                 tol_err,tol_res,tol_det : in quad_double;
                 max : in natural32; fail : out boolean;
                 nit,crtp : out natural32 ) is

    t,vnd,zt : quad_double;
    orientation : Complex_Number;
 
  begin
    fail := false; nit := 0;
    if output then
      put(file,"Determinant : "); put(file,nd,3);
      put_line(file,"  Tangent vector is"); Write_Tangent(file,dx);
    end if;
    crtp := 0;
    vnd := AbsVal(nd);
    if vnd < tol_det then
      crtp := 1;
    else 
      orientation := Inner_Product(dx,pt);
      if output then
        put(file,"Orientation : "); put(file,orientation); new_line(file);
      end if;
      if (REAL_PART(orientation) < 0.0) then
        crtp := 2;
      else
        t := REAL_PART(x(x'last));
        Monitor_Determinants(file,vt,dt,i,t,vnd,crtp,zt);
      end if;
    end if;
  end Monitor_Singularity;

  procedure Silent_Monitor_Singularity
               ( nd : in Complex_Number;
                 vt,dt : in out Quad_Double_Vectors.Vector;
                 i : in out integer32;
                 f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out QuadDobl_Complex_Vectors.Vector;
                 px,pt,dx : in QuadDobl_Complex_Vectors.Vector;
                 tol_err,tol_res,tol_det : in quad_double;
                 max : in natural32; fail : out boolean;
                 nit,crtp : out natural32 ) is

    t,vnd,zt : quad_double;
    orientation : Complex_Number;
 
  begin
    fail := false; nit := 0;
    crtp := 0;
    vnd := AbsVal(nd);
    if vnd < tol_det then
      crtp := 1;
    else 
      orientation := Inner_Product(dx,pt);
      if (REAL_PART(orientation) < 0.0) then
        crtp := 2;
      else
        t := REAL_PART(x(x'last));
        Silent_Monitor_Determinants(vt,dt,i,t,vnd,crtp,zt);
      end if;
    end if;
  end Silent_Monitor_Singularity;

end QuadDobl_Quad_Turn_Points;
