with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with Double_Double_Poly_Systems_io;      use Double_Double_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_to_Real_Poly;      use DoblDobl_Complex_to_Real_Poly;
with DoblDobl_Solution_Diagnostics;
with DoblDobl_Quad_Parameters;
with DoblDobl_Quad_Turn_Points;          use DoblDobl_Quad_Turn_Points;
with DoblDobl_Quad_Turn_Points_io;       use DoblDobl_Quad_Turn_Points_io;

package body DoblDobl_Quad_Sweepers is

  function Real_Part ( x : DoblDobl_Complex_Vectors.Vector )
                     return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      res(i) := DoblDobl_Complex_Numbers.REAL_PART(x(i));
    end loop;
    return res;
  end Real_Part;

  procedure Interactive_Real_Sweep
               ( nq,nv : in natural32; eigval : in boolean;
                 p : in Double_Double_Poly_Systems.Poly_Sys;
                 sols : in Solution_List ) is

    x,dx,pt,px,x0 : Double_Double_Vectors.Vector(1..integer32(nv));
    y,det : Double_Double_Vectors.Vector(1..integer32(nq));
    xt,yd : Double_Double_Vectors.Vector(1..3);
    wi : integer32 := 0; 
    L : DoblDobl_Complex_Vectors.Vector(1..integer32(nq)-1);
    v : DoblDobl_Complex_VecVecs.VecVec(1..integer32(nq)-1);
    ep : Double_Double_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Double_Double_Poly_SysFun.Create(p);
    jm : Double_Double_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Double_Double_Jaco_Matrices.Create(p);
    ejm : Double_Double_Jaco_Matrices.Eval_Jaco_Mat(p'range,x'range)
        := Double_Double_Jaco_Matrices.Create(jm);
    step : double_double := DoblDobl_Quad_Parameters.max_step_size;
    step_cnt : integer32;
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    tol_det : constant double_double
            := DoblDobl_Quad_Parameters.determinant_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    ans : character;
    nb,crtp : natural32;
    fail,critical : boolean := false;
    nd : double_double;
    ls : Link_to_Solution;
    tmp : Solution_List := sols;

    use Double_Double_Vectors;
    use Double_Double_Poly_SysFun;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      new_line;
      x := Real_Part(ls.v); -- Read_Initial_Vector(x);
      px := x;
      y := Eval(ep,x);
      put_line("The first evaluation gives"); put_line(y);
      pt := (pt'range => create(0.0));
      pt(pt'last) := create(1.0);
      Step_Size_Control(step,-1); step_cnt := 0;
      loop
        if not fail then
          if eigval then
            Tangent_Minors_and_Eigenvectors(ejm,x,dx,det,L,v);
            nd := det(det'last);
            Report_Minors_and_Eigenvectors(standard_output,det,L,v);
          else
            Tangent_and_Determinant(ejm,x,dx,nd);
          end if;
          Monitor_Singularity
            (standard_output,true,nd,xt,yd,wi,ep,ejm,x,px,pt,dx,
            tol_err,tol_res,tol_det,maxits,fail,nb,crtp);
          critical := (crtp > 0);
          exit when critical;
        end if;
        x0 := x; x := x + step*dx;
        Correct_Solution(ep,ejm,dx,x,tol_err,tol_res,fail,nb,maxits);
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y'); 
        if not fail
         then px := x0; pt := dx; step_cnt := step_cnt + 1;
         else x := x0;            step_cnt := 0;
        end if;
        Step_Size_Control(step,step_cnt);
      end loop;
      if crtp = 0 then
        put_line("no critical point detected");
      elsif crtp = 2 then
        put_line("ORIENTATION OF TANGENT FLIPPED");
        Interactive_Shoot_Turn(ep,ejm,px,pt,x,dx,step);
      else
        put_line("CRITICAL POINT DETECTED");
      end if;
      put_line("The final solution :"); Write_Vector(x);
      tmp := Tail_Of(tmp);
    end loop;
    Double_Double_Poly_SysFun.Clear(ep);
    Double_Double_Jaco_Matrices.Clear(jm);
    Double_Double_Jaco_Matrices.Clear(ejm);
  end Interactive_Real_Sweep;

  procedure Append ( first,last : in out Solution_List;
                     x : in Double_Double_Vectors.Vector ) is

    sol : Solution(x'last);

  begin
    sol.t := Create(x(x'last));
    sol.m := 1;
    for i in x'range loop
      sol.v(i) := Create(x(i));
    end loop;
    sol.err := create(0.0);
    sol.rco := create(1.0);
    sol.res := create(0.0);
    Append(first,last,sol);
  end Append;

  procedure Interactive_Real_Sweep
               ( nq,nv : in natural32; eigval : in boolean;
                 p : in Double_Double_Poly_Systems.Poly_Sys;
                 sols : in Solution_List; repsols : in out Solution_List ) is

    x,dx,pt,px,x0 : Double_Double_Vectors.Vector(1..integer32(nv));
    y,det : Double_Double_Vectors.Vector(1..integer32(nq));
    xt,yd : Double_Double_Vectors.Vector(1..3);
    wi : integer32 := 0;
    L : DoblDobl_Complex_Vectors.Vector(1..integer32(nq-1));
    v : DoblDobl_Complex_VecVecs.VecVec(1..integer32(nq-1));
    ep : Double_Double_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Double_Double_Poly_SysFun.Create(p);
    jm : Double_Double_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Double_Double_Jaco_Matrices.Create(p);
    ejm : Double_Double_Jaco_Matrices.Eval_Jaco_Mat(p'range,x'range)
        := Double_Double_Jaco_Matrices.Create(jm);
    step : double_double := DoblDobl_Quad_Parameters.max_step_size;
    step_cnt : integer32;
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    tol_det : constant double_double
            := DoblDobl_Quad_Parameters.determinant_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    ans : character;
    nb,crtp : natural32;
    fail,critical : boolean := false;
    nd : double_double;
    ls : Link_to_Solution;
    tmp : Solution_List := sols;
    repsols_last : Solution_List;

    use Double_Double_Vectors;
    use Double_Double_Poly_SysFun;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      new_line;
      x := Real_Part(ls.v); -- Read_Initial_Vector(x);
      px := x;
      y := Eval(ep,x);
      put_line("The first evaluation gives"); put_line(y);
      pt := (pt'range => create(0.0));
      pt(pt'last) := create(1.0);
      Step_Size_Control(step,-1); step_cnt := 0;
      Append(repsols,repsols_last,x);
      loop
        if not fail then
          if eigval then
            Tangent_Minors_and_Eigenvectors(ejm,x,dx,det,L,v);
            nd := det(det'last);
            Report_Minors_and_Eigenvectors(standard_output,det,L,v);
          else
            Tangent_and_Determinant(ejm,x,dx,nd);
          end if;
          Monitor_Singularity
            (standard_output,true,nd,xt,yd,wi,ep,ejm,x,px,pt,dx,
             tol_err,tol_res,tol_det,maxits,fail,nb,crtp);
          critical := (crtp > 0);
          exit when critical;
        end if;
        x0 := x; x := x + step*dx;
        Correct_Solution(ep,ejm,dx,x,tol_err,tol_res,fail,nb,maxits);
        put("Record solution ? (y/n) "); Ask_Yes_or_No(ans);
        if ans = 'y'
         then Append(repsols,repsols_last,x);
        end if;
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y'); 
        if not fail
         then px := x0; pt := dx; step_cnt := step_cnt + 1;
         else x := x0;  dx := pt; step_cnt := 0;
        end if;
        Step_Size_Control(step,step_cnt);
      end loop;
      if crtp = 0 then
        put_line("no critical point detected");
      elsif crtp = 2 then
        put_line("ORIENTATION OF TANGENT FLIPPED");
        Interactive_Shoot_Turn(ep,ejm,px,pt,x,dx,step);
      else
        put_line("CRITICAL POINT DETECTED");
      end if;
      put_line("The final solution :"); Write_Vector(x);
      Append(repsols,repsols_last,x);
      tmp := Tail_Of(tmp);
    end loop;
    Double_Double_Poly_SysFun.Clear(ep);
    Double_Double_Jaco_Matrices.Clear(jm);
    Double_Double_Jaco_Matrices.Clear(ejm);
  end Interactive_Real_Sweep;

  procedure Silent_Real_Sweep
               ( eigval : in boolean;
                 nq,nv : in natural32; t : in double_double;
                 f : in Double_Double_Poly_SysFun.Eval_Poly_Sys;
                 jf : in Double_Double_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out Double_Double_Vectors.Vector;
                 dx : out Double_Double_Vectors.Vector ) is

    pt,px,x0 : Double_Double_Vectors.Vector(1..integer32(nv));
    det : Double_Double_Vectors.Vector(1..integer32(nq));
    xt,yd : Double_Double_Vectors.Vector(1..3);
    wi : integer32 := 0;
    L : DoblDobl_Complex_Vectors.Vector(1..integer32(nq-1));
    v : DoblDobl_Complex_VecVecs.VecVec(1..integer32(nq-1));
    max : constant natural32 := DoblDobl_Quad_Parameters.max_predictor_steps;
    step : double_double := DoblDobl_Quad_Parameters.max_step_size;
    step_cnt : integer32;
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    tol_det : constant double_double
            := DoblDobl_Quad_Parameters.determinant_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    nb,crtp : natural32;
    fail : boolean := false;
    critical : boolean := false;
    nd : double_double;

    use Double_Double_Vectors;
    use Double_Double_Poly_SysFun;

  begin
    px := x;
    pt := (pt'range => create(0.0));
    pt(pt'last) := create(1.0);
    Step_Size_Control(step,-1); step_cnt := 0;
    for k in 1..max loop
      if not fail then
        if eigval then
          Tangent_Minors_and_Eigenvectors(jf,x,dx,det,L,v);
          nd := det(det'last);
        else
          Tangent_and_Determinant(jf,x,dx,nd);
        end if;
        Silent_Monitor_Singularity(nd,xt,yd,wi,f,jf,x,px,pt,dx,
          tol_err,tol_res,tol_det,maxits,fail,nb,crtp);
        critical := (crtp > 0);
        exit when critical;
      end if;
      x0 := x; x := x + step*dx;
      Correct_Solution(f,jf,dx,x,tol_err,tol_res,fail,nb,maxits);
      exit when (x(x'last) >= t);
      if not fail
       then px := x0; pt := dx; step_cnt := step_cnt + 1;
       else x := x0;  dx := pt; step_cnt := 0;
      end if;
      Step_Size_Control(step,step_cnt);
    end loop;
    if crtp = 0 then
      if x(x'last) > t
       then Target_Correction(f,jf,t,x,tol_err,tol_res,fail,nb,maxits);
      end if;
    elsif crtp = 2 then
      Shoot_Turn(f,jf,px,pt,x,dx,step,create(1.0E-10),4);
    end if;
  end Silent_Real_Sweep;

  procedure Start_Real_Sweep
               ( eigval : in boolean;
                 nq,nv : in natural32; t : in double_double;
                 f : in Double_Double_Poly_SysFun.Eval_Poly_Sys;
                 jf : in Double_Double_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out Double_Double_Vectors.Vector;
                 dx : out Double_Double_Vectors.Vector ) is

    pt,px,x0 : Double_Double_Vectors.Vector(1..integer32(nv));
    det : Double_Double_Vectors.Vector(1..integer32(nq));
    xt,yd : Double_Double_Vectors.Vector(1..3);
    wi : integer32 := 0;
    L : DoblDobl_Complex_Vectors.Vector(1..integer32(nq-1));
    v : DoblDobl_Complex_VecVecs.VecVec(1..integer32(nq-1));
    max : constant natural32 := DoblDobl_Quad_Parameters.max_predictor_steps;
    step : double_double := DoblDobl_Quad_Parameters.max_step_size;
    step_cnt : integer32;
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    tol_det : constant double_double
            := DoblDobl_Quad_Parameters.determinant_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    nb,crtp : natural32;
    fail : boolean := false;
    critical : boolean := false;
    nd : double_double;

    use Double_Double_Vectors;
    use Double_Double_Poly_SysFun;

  begin
    px := x;
    pt := (pt'range => create(0.0));
    pt(pt'last) := create(1.0);
    Step_Size_Control(step,-1); step_cnt := 0;
    for k in 1..max loop
      if not fail then
        if eigval then
          Tangent_Minors_and_Eigenvectors(jf,x,dx,det,L,v);
          nd := det(det'last);
          Report_Minors_and_Eigenvectors(standard_output,det,L,v);
        else
          Tangent_and_Determinant(jf,x,dx,nd);
        end if;
        Monitor_Singularity(standard_output,true,nd,xt,yd,wi,f,jf,x,px,pt,dx,
                            tol_err,tol_res,tol_det,maxits,fail,nb,crtp);
        critical := (crtp > 0);
        exit when critical;
      end if;
      x0 := x; x := x + step*dx;
      Correct_Solution(f,jf,dx,x,tol_err,tol_res,fail,nb,maxits);
      exit when (x(x'last) >= t);
      if not fail
       then px := x0; pt := dx; step_cnt := step_cnt + 1;
       else x := x0;  dx := pt; step_cnt := 0;
      end if;
      Step_Size_Control(step,step_cnt);
    end loop;
    if crtp = 0 then
      put_line("no critical point detected");
      if x(x'last) > t
       then Target_Correction(f,jf,t,x,tol_err,tol_res,fail,nb,maxits);
      end if;
    elsif crtp = 2 then
      Shoot_Turn(f,jf,px,pt,x,dx,step,create(1.0E-10),4);
    else
      put_line("CRITICAL POINT DETECTED");
    end if;
  end Start_Real_Sweep;

  procedure Start_Real_Sweep
               ( file : in file_type; output,eigval : in boolean;
                 nq,nv : in natural32; t : in double_double;
                 f : in Double_Double_Poly_SysFun.Eval_Poly_Sys;
                 jf : in Double_Double_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out Double_Double_Vectors.Vector;
                 dx : out Double_Double_Vectors.Vector ) is

    pt,px,x0 : Double_Double_Vectors.Vector(1..integer32(nv));
    y,det : Double_Double_Vectors.Vector(1..integer32(nq));
    xt,yd : Double_Double_Vectors.Vector(1..3);
    wi : integer32 := 0;
    L : DoblDobl_Complex_Vectors.Vector(1..integer32(nq)-1);
    v : DoblDobl_Complex_VecVecs.VecVec(1..integer32(nq)-1);
    max : constant natural32 := DoblDobl_Quad_Parameters.max_predictor_steps;
    step : double_double := DoblDobl_Quad_Parameters.max_step_size;
    step_cnt : integer32;
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    tol_det : constant double_double
            := DoblDobl_Quad_Parameters.determinant_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    nb,crtp : natural32;
    fail : boolean := false;
    critical : boolean := false;
    nd : double_double;

    use Double_Double_Vectors;
    use Double_Double_Poly_SysFun;

  begin
    if output then
      y := Eval(f,x);
      new_line(file);
      put_line(file,"The first evaluation gives"); put_line(file,y);
      Step_Size_Control(file,step,-1);
    else
      Step_Size_Control(step,-1);
    end if;
    px := x; step_cnt := 0;
    pt := (pt'range => create(0.0));
    pt(pt'last) := create(1.0);
    for k in 1..max loop
      if not fail then
        if eigval then
          Tangent_Minors_and_Eigenvectors(jf,x,dx,det,L,v);
          nd := det(det'last);
          Report_Minors_and_Eigenvectors(file,det,L,v);
        else
          Tangent_and_Determinant(jf,x,dx,nd);
        end if;
        Monitor_Singularity(file,output,nd,xt,yd,wi,f,jf,x,px,pt,dx,
                            tol_err,tol_res,tol_det,maxits,fail,nb,crtp);
        critical := (crtp > 0);
      end if;
      exit when critical;
      x0 := x; x := x + step*dx;
      if output
       then Correct_Solution(file,f,jf,dx,x,tol_err,tol_res,fail,nb,maxits);
       else Correct_Solution(f,jf,dx,x,tol_err,tol_res,fail,nb,maxits);
      end if;
      exit when (x(x'last) >= t);
      if not fail
       then px := x0; pt := dx; step_cnt := step_cnt + 1;
       else x := x0;  dx := pt; step_cnt := 0;
      end if;
      if output
       then Step_Size_Control(file,step,step_cnt);
       else Step_Size_Control(step,step_cnt);
      end if;
    end loop;
    if crtp = 0 then
      put_line(file,"no critical point detected");
      if x(x'last) > t
       then Target_Correction(file,f,jf,t,x,tol_err,tol_res,fail,nb,maxits);
      end if;
    elsif crtp = 2 then
      put_line(file,"ORIENTATION OF TANGENT FLIPPED");
      Shoot_Turn(f,jf,px,pt,x,dx,step,create(1.0E-10),4);
    else
      put_line(file,"CRITICAL POINT DETECTED");
    end if;
  end Start_Real_Sweep;

  procedure Interactive_Complex_Sweep
               ( nq,nv : in natural32;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in Solution_List ) is

    x,dx,pt,px,x0 : DoblDobl_Complex_Vectors.Vector(1..integer32(nv));
    ep : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := DoblDobl_Complex_Poly_SysFun.Create(p);
    y : DoblDobl_Complex_Vectors.Vector(1..integer32(nq));
    jm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(p);
    ejm : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(p'range,x'range)
        := DoblDobl_Complex_Jaco_Matrices.Create(jm);
    step : Complex_Number := Create(DoblDobl_Quad_Parameters.max_step_size);
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    step_cnt : integer32 := 0;
    orientation : Complex_Number;
    ans : character;
    fail : boolean := false;
    nb : natural32;
    ls : Link_to_Solution;
    tmp : Solution_List := sols;

    use DoblDobl_Complex_Vectors;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      new_line;
      x := ls.v; -- Read_Initial_Vector(x);
      px := x; x0 := x;
      y := DoblDobl_Complex_Poly_SysFun.Eval(ep,x);
      put_line("The first evaluation gives"); put_line(y);
      pt := (pt'range => Create(integer(0)));
      pt(pt'last) := Create(integer(1));
      Step_Size_Control(step,-1); step_cnt := 0;
      loop
        if not fail
         then dx := Tangent(ejm,x);
        end if;
        put_line("The tangent vector is"); Write_Tangent(dx);
        orientation := Inner_Product(dx,pt);
        put("Orientation : "); put(orientation); new_line;
        exit when (REAL_PART(orientation) < 0.0);
        x0 := x; x := x + step*dx;
        Correct_Solution(ep,ejm,dx,x,tol_err,tol_res,fail,nb,maxits);
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        if not fail
         then px := x0; pt := dx; step_cnt := step_cnt + 1;
         else x := x0;  dx := pt; step_cnt := 0;
        end if;
        Step_Size_Control(step,step_cnt);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    if REAL_PART(orientation) < 0.0 then
      put_line("ORIENTATION OF TANGENT FLIPPED");
      Interactive_Shoot_Turn(ep,ejm,px,pt,x,dx,step);
    end if;
    put_line("The final solution :"); Write_Vector(x);
    DoblDobl_Complex_Poly_SysFun.Clear(ep);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(ejm);
  end Interactive_Complex_Sweep;

  procedure Start_Complex_Sweep
               ( nq,nv : in natural32; t : in double_double;
                 f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out DoblDobl_Complex_Vectors.Vector;
                 dx : out DoblDobl_Complex_Vectors.Vector ) is

    pt,px,x0 : DoblDobl_Complex_Vectors.Vector(1..integer32(nv));
    max : constant natural32 := DoblDobl_Quad_Parameters.max_predictor_steps;
    step : Complex_Number := Create(DoblDobl_Quad_Parameters.max_step_size);
    step_cnt : integer32 := 0;
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    orientation : Complex_Number;
    fail : boolean := false;
    nb : natural32;

    use DoblDobl_Complex_Vectors;

  begin
    px := x; x0 := x;
    pt := (pt'range => Create(integer(0)));
    pt(pt'last) := Create(integer(1));
    Step_Size_Control(step,-1); step_cnt := 0;
    for k in 1..max loop
      if not fail
       then dx := Tangent(jf,x);
      end if;
      orientation := Inner_Product(dx,pt);
      exit when (REAL_PART(orientation) < 0.0);
      x0 := x; x := x + step*dx;
      Correct_Solution(f,jf,dx,x,tol_err,tol_res,fail,nb,maxits);
      exit when (REAL_PART(x(x'last)) > t);
      if not fail
       then px := x0; pt := dx; step_cnt := step_cnt + 1;
       else x := x0;  dx := pt; step_cnt := 0;
      end if;
      Step_Size_Control(step,step_cnt);
    end loop;
    if REAL_PART(orientation) < 0.0 then
      declare
        dd_tol_step : constant double_double := create(1.0E-10);
        tol_step : constant Complex_Number := Create(dd_tol_step);
      begin
        Shoot_Turn(f,jf,px,pt,x,dx,step,tol_step,4);
      end;
    end if;
  end Start_Complex_Sweep;

  procedure Start_Complex_Sweep
               ( file : in file_type; output : in boolean;
                 nq,nv : in natural32; t : in double_double;
                 f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                 x : in out DoblDobl_Complex_Vectors.Vector;
                 dx : out DoblDobl_Complex_Vectors.Vector ) is

    pt,px,x0 : DoblDobl_Complex_Vectors.Vector(1..integer32(nv));
    y : DoblDobl_Complex_Vectors.Vector(1..integer32(nq));
    xt,yd : Double_Double_Vectors.Vector(1..3);
    wi : integer32 := 0;
    max : constant natural32 := DoblDobl_Quad_Parameters.max_predictor_steps;
    step : Complex_Number := Create(DoblDobl_Quad_Parameters.max_step_size);
    step_cnt : integer32 := 0;
    tol_err : constant double_double
            := DoblDobl_Quad_Parameters.increment_tolerance;
    tol_res : constant double_double
            := DoblDobl_Quad_Parameters.residual_tolerance;
    tol_det : constant double_double
            := DoblDobl_Quad_Parameters.determinant_tolerance;
    maxits : constant natural32
           := DoblDobl_Quad_Parameters.max_corrector_steps;
    orientation : Complex_Number;
    fail : boolean := false;
    nb,crtp : natural32;
    nd : Complex_Number;

    use DoblDobl_Complex_Vectors;

  begin
    px := x; x0 := x;
    y := DoblDobl_Complex_Poly_SysFun.Eval(f,x);
    pt := (pt'range => Create(integer(0)));
    pt(pt'last) := Create(integer(1));
    if output  then
      new_line(file);
      put_line(file,"The first evaluation gives"); put_line(file,y);
      Step_Size_Control(file,step,-1); 
    else
      Step_Size_Control(step,-1); 
    end if;
    step_cnt := 0;
    for k in 1..max loop
      if not fail
       then Tangent_and_Determinant(jf,x,dx,nd);
      end if;
      orientation := Inner_Product(dx,pt);
      if output then
        put_line(file,"The tangent vector is"); Write_Tangent(file,dx);
        put(file,"Orientation : "); put(file,orientation); new_line(file);
      end if;
      if not fail then
        Monitor_Singularity(file,output,nd,xt,yd,wi,f,jf,x,px,pt,dx,
                            tol_err,tol_res,tol_det,maxits,fail,nb,crtp);
      end if;
      exit when (REAL_PART(orientation) < 0.0);
      x0 := x; x := x + step*dx;
      if output
       then Correct_Solution(file,f,jf,dx,x,tol_err,tol_res,fail,nb,maxits);
       else Correct_Solution(f,jf,dx,x,tol_err,tol_res,fail,nb,maxits);
      end if;
      exit when (REAL_PART(x(x'last)) > t);
      if not fail
       then px := x0; pt := dx; step_cnt := step_cnt + 1;
       else x := x0;  dx := pt; step_cnt := 0;
      end if;
      if output
       then Step_Size_Control(file,step,step_cnt);
       else Step_Size_Control(step,step_cnt);
      end if;
    end loop;
    if crtp = 0 then
      put_line(file,"no critical point detected");
      if REAL_PART(x(x'last)) > t
       then Target_Correction(file,f,jf,t,x,tol_err,tol_res,fail,nb,maxits);
      end if;
    elsif crtp = 2 then
      declare
        dd_tol_step : constant double_double := create(1.0E-10);
        tol_step : constant Complex_Number := create(dd_tol_step);
      begin
        put_line(file,"ORIENTATION OF TANGENT FLIPPED");
        Shoot_Turn(f,jf,px,pt,x,dx,step,tol_step,4);
      end;
    else
      put_line(file,"CRITICAL POINT DETECTED");
    end if;
  end Start_Complex_Sweep;

-- DRIVER ROUTINES :

  procedure Ask_Sweep_Parameters ( target : out double_double ) is
  begin
    new_line;
    target := create(1.0);
    put("Give end target value for t : "); get(target);
  end Ask_Sweep_Parameters;

  procedure Ask_Sweep_Parameters
              ( target : out double_double; output : out boolean ) is

    ans : character;

  begin
    Ask_Sweep_Parameters(target);
    put("Intermediate output of corrector ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y'); 
  end Ask_Sweep_Parameters;

  procedure Ask_Sweep_Parameters
              ( target : out double_double; output,eigval : out boolean ) is

    ans : character;

  begin
    Ask_Sweep_Parameters(target,output);
    put("Monitor eigenvalues ? (y/n) ");
    Ask_Yes_or_No(ans);
    eigval := (ans = 'y'); 
  end Ask_Sweep_Parameters;

  procedure Run_Sweep
              ( nq,nv : in natural32;
                p : in Double_Double_Poly_Systems.Poly_Sys;
                sols : in Solution_List ) is

    file : file_type;
    timer : Timing_Widget;
    t : double_double;
    rx,rdx : Double_Double_Vectors.Vector(1..integer32(nv));
    cx,cdx : DoblDobl_Complex_Vectors.Vector(1..integer32(nv));
    q : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range)
      := Convert_Real_to_Complex(p);
    ep : Double_Double_Poly_SysFun.Eval_Poly_Sys(p'range)
       := Double_Double_Poly_SysFun.Create(p);
    eq : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(q'range)
       := DoblDobl_Complex_Poly_SysFun.Create(q);
    jpm : Double_Double_Jaco_Matrices.Jaco_Mat(p'range,1..integer32(nv))
        := Double_Double_Jaco_Matrices.Create(p);
    jqm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(q'range,1..integer32(nv))
        := DoblDobl_Complex_Jaco_Matrices.Create(q);
    ejpm : Double_Double_Jaco_Matrices.Eval_Jaco_Mat
             (p'range,1..integer32(nv))
         := Double_Double_Jaco_Matrices.Create(jpm);
    ejqm : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat
             (q'range,1..integer32(nv))
         := DoblDobl_Complex_Jaco_Matrices.Create(jqm);
    otp,evl : boolean;
    ls : Link_to_Solution;
    tmp : Solution_List := sols;
    tol : constant double_double := create(1.0E-8);
    len : constant natural32 := Length_Of(sols);

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    put(file,nq,1); put(file," "); put(file,nv,1); new_line(file);
    put(file,p);
    Ask_Sweep_Parameters(t,otp,evl);
    new_line;
    DoblDobl_Quad_Parameters.Tune;
    new_line;
    put("running a real sweep on "); put(len,1);
    put_line(" solutions ...");
    new_line;
    tstart(timer);
    for i in 1..len loop
      ls := Head_Of(tmp);
      new_line(file);
      if not DoblDobl_Solution_Diagnostics.Is_Real(ls.all,hi_part(tol)) then
        cx := ls.v;
        put(file,"Starting a complex sweep at solution ");
        put(file,i,1); put_line(file," :");
        Write_Vector(file,cx); new_line(file);
        Start_Complex_Sweep(file,otp,nq,nv,t,eq,ejqm,cx,cdx);
        put_line(file,"The solution and its tangent at the end");
        Write_Vector_and_its_Tangent(file,cx,cdx);
        ls.v := cx; ls.t := cx(cx'last);
      else
        rx := Real_Part(ls.v);
        put(file,"Starting a real sweep at solution ");
        put(file,i,1); put_line(file," :");
        Write_Vector(file,rx); new_line(file);
        Start_Real_Sweep(file,otp,evl,nq,nv,t,ep,ejpm,rx,rdx);
        put_line(file,"The solution and its tangent at the end");
        Write_Vector_and_its_Tangent(file,rx,rdx);
        ls.t := Create(rx(rx'last));
        for i in rx'range loop
          ls.v(i) := Create(rx(i)); 
        end loop;
      end if;
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    DoblDobl_Complex_Poly_Systems.Clear(q);
    Double_Double_Poly_SysFun.Clear(ep);
    DoblDobl_Complex_Poly_SysFun.Clear(eq);
    Double_Double_Jaco_Matrices.Clear(jpm);
    DoblDobl_Complex_Jaco_Matrices.Clear(jqm);
    Double_Double_Jaco_Matrices.Clear(ejpm);
    DoblDobl_Complex_Jaco_Matrices.Clear(ejqm);
    Write_Sweep_Summary(file,sols);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"running real sweep");
  end Run_Sweep;

  procedure Run_Sweep
              ( nq,nv : in natural32;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in Solution_List ) is

    file : file_type;
    timer : Timing_Widget;
    otp : boolean;
    t : double_double;
    x,dx : DoblDobl_Complex_Vectors.Vector(1..integer32(nv));
    ep : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(p'range)
       := DoblDobl_Complex_Poly_SysFun.Create(p);
    jm : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := DoblDobl_Complex_Jaco_Matrices.Create(p);
    ejm : DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat(p'range,x'range)
        := DoblDobl_Complex_Jaco_Matrices.Create(jm);
    ls : Link_to_Solution;
    tmp : Solution_List := sols;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
   -- put(file,nq,nv,p);
    put(file,p);
    Ask_Sweep_Parameters(t,otp);
    new_line;
    DoblDobl_Quad_Parameters.Tune;
    new_line;
    put("running a complex sweep on "); put(Length_Of(sols),1);
    put_line(" solutions ...");
    new_line;
    tstart(timer);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      x := ls.v; -- Read_Initial_Vector(x);
      new_line(file);
      put_line(file,"Starting the sweep at the point");
      Write_Vector(file,x);
      new_line(file);
      Start_Complex_Sweep(file,otp,nq,nv,t,ep,ejm,x,dx);
      put_line(file,"The solution and its tangent at the end");
      Write_Vector_and_its_Tangent(file,x,dx);
      ls.v := x; ls.t := x(x'last);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
    DoblDobl_Complex_Poly_SysFun.Clear(ep);
    DoblDobl_Complex_Jaco_Matrices.Clear(jm);
    DoblDobl_Complex_Jaco_Matrices.Clear(ejm);
    Write_Sweep_Summary(file,sols);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"running a complex sweep");
  end Run_Sweep;

end DoblDobl_Quad_Sweepers;
