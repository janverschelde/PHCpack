with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Standard_IncFix_Continuation;      use Standard_IncFix_Continuation;
--with Standard_Affine_Planes;            use Standard_Affine_Planes;
--with Standard_Affine_Solutions;         use Standard_Affine_Solutions;
with Standard_Point_Coordinates;        use Standard_Point_Coordinates;

package body Affine_Sampling_Machine is

  function Eval ( p : Eval_Poly_Sys; b : Vector; v : VecVec; c : Vector )
                return Vector is

    x : constant Vector := Affine_Expand(c,b,v);
    y : constant Vector := Eval(p,x);

  begin
    return y;
  end Eval;

  function Diff ( jm : Eval_Jaco_Mat; b : Vector; v : VecVec; c : Vector )
                return Matrix is

    res : Matrix(jm'range(1),v'range);
    x : constant Vector := Affine_Expand(c,b,v);
    eva : constant Matrix(jm'range(1),jm'range(2)) := Eval(jm,x);

  begin
    for i in res'range(1) loop
      for j in v'range loop
        res(i,j) := Create(0.0);
        for k in eva'range(2) loop
          res(i,j) := res(i,j) + eva(i,k)*v(j)(k);
        end loop;
      end loop;
    end loop;
    return res;
  end Diff;

  procedure Moving_Plane
              ( start_b,target_b : in Vector; start_v,target_v : in VecVec;
                t : in Complex_Number; b : out Vector; v : in out VecVec ) is

    one_min_t : constant Complex_Number := 1.0 - t;

  begin
    for i in b'range loop
      b(i) := one_min_t*start_b(i) + t*target_b(i);
    end loop;
    for i in v'range loop
      if v(i) = null
       then v(i) := new Standard_Complex_Vectors.Vector(b'range);
      end if;
      for j in v(i)'range loop
        v(i)(j) := one_min_t*start_v(i)(j) + t*target_v(i)(j);
      end loop;
    end loop;
  end Moving_Plane;

  procedure Silent_LU_Newton_Refiner
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; c : in out Vector;
                err,rco,res : out double_float;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    k : constant integer32 := c'length;
    y : Vector(p'range);
    jac : Matrix(jm'range(1),c'range);
    ipvt : Standard_Integer_Vectors.Vector(c'range);
    dc : Vector(c'range);

  begin
    y := Eval(p,b,v,c);
    for i in 1..maxit loop
      jac := Diff(jm,b,v,c);
      lufco(jac,k,ipvt,rco);
      dc := -y;
      lusolve(jac,k,ipvt,dc);
      c := c + dc;
      y := Eval(p,b,v,c);
      err := Max_Norm(dc);
      res := Max_Norm(y);
      exit when (err < epsxa) or (res < epsfa);
    end loop;
  end Silent_LU_Newton_Refiner;

  procedure Reporting_LU_Newton_Refiner
              ( file : in file_type;
                p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; c : in out Vector;
                err,rco,res : out double_float;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    k : constant integer32 := c'length;
    y : Vector(p'range);
    jac : Matrix(jm'range(1),c'range);
    ipvt : Standard_Integer_Vectors.Vector(c'range);
    dc : Vector(c'range);

  begin
    y := Eval(p,b,v,c);
    for i in 1..maxit loop
      jac := Diff(jm,b,v,c);
      lufco(jac,k,ipvt,rco);
      dc := -y;
      lusolve(jac,k,ipvt,dc);
      c := c + dc;
      y := Eval(p,b,v,c);
      err := Max_Norm(dc);
      res := Max_Norm(y);
      put(file,"  err : "); put(file,err,3);
      put(file,"  rco : "); put(file,rco,3);
      put(file,"  res : "); put(file,res,3); new_line(file);
      exit when (err < epsxa) or (res < epsfa);
    end loop;
    put_line(file,"After Newton refinement : "); put_line(file,y);
  end Reporting_LU_Newton_Refiner;

  procedure Silent_SV_Newton_Refiner
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; c : in out Vector;
                err,rco,res : out double_float;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    n : constant integer32 := jm'length(1);
    k : constant integer32 := c'length;
    mm : constant integer32 := Min0(n+1,k);
    y : Vector(p'range);
    jac : Matrix(jm'range(1),c'range);
    u : Matrix(jac'range(1),jac'range(1));
    vv : Matrix(c'range,c'range);
    s : Vector(1..mm);
    e,dc : Vector(c'range);
    info : integer32;

  begin
    y := Eval(p,b,v,c);
    for i in 1..maxit loop
      jac := Diff(jm,b,v,c);
      SVD(jac,n,k,s,e,u,vv,11,info);
      dc := Solve(u,vv,s,-y);
      c := c + dc;
      y := Eval(p,b,v,c);
      err := Max_Norm(dc);
      rco := Inverse_Condition_Number(s);
      res := Max_Norm(y);
      exit when (err < epsxa) or (res < epsfa);
    end loop;
  end Silent_SV_Newton_Refiner;

  procedure Reporting_SV_Newton_Refiner
              ( file : in file_type;
                p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; c : in out Vector;
                err,rco,res : out double_float;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    n : constant integer32 := jm'length(1);
    k : constant integer32 := c'length;
    mm : constant integer32 := Min0(n+1,k);
    y : Vector(p'range);
    jac : Matrix(jm'range(1),c'range);
    u : Matrix(jac'range(1),jac'range(1));
    vv : Matrix(c'range,c'range);
    s : Vector(1..mm);
    e,dc : Vector(c'range);
    info : integer32;
    tol_rank : double_float := 1.0E-15;
    tol_drop : constant double_float := 1.0E-4;
    rnk : integer32;

  begin
    y := Eval(p,b,v,c);
    rnk := y'length;
    for i in 1..maxit loop
      jac := Diff(jm,b,v,c);
      SVD(jac,n,k,s,e,u,vv,11,info);
      put_line(file,"The singular values : ");
      put_line(file,s);
      rnk := Rank(s,tol_rank);
      put(file,"The numerical rank is "); 
      put(file,rnk,1); new_line(file);
      if s(rnk) < s(rnk-1)*tol_drop then
        put_line(file,"Executing a rank drop...");
        tol_rank := AbsVal(s(rnk))*10.0;
        rco := 0.0;
      else
        put_line(file,"Observed no rank drop.");
        if rnk = y'length(1)
         then rco := Inverse_Condition_Number(s);
        end if;
      end if;
      dc := Solve(u,vv,s,-y,tol_rank);
      c := c + dc;
      y := Eval(p,b,v,c);
      err := Max_Norm(dc);
      res := Max_Norm(y);
      put(file,"  err : "); put(file,err,3);
      put(file,"  rco : "); put(file,rco,3);
      put(file,"  res : "); put(file,res,3); new_line(file);
      exit when (err < epsxa) or (res < epsfa) or (rco = 0.0);
    end loop;
    put_line(file,"After Newton refinement : "); put_line(file,y);
  end Reporting_SV_Newton_Refiner;

  procedure Solution_Evaluation
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                sol,b : Vector; v : in VecVec ) is

    c : Vector(v'range) := Project(sol,b,v);
    y : Vector(p'range) := Eval(p,b,v,c);
    t : constant Vector(v'range) := (v'range => Create(0.0));
    err,rco,res : double_float;

  begin
    put_line("The solution evaluated at computed basis : "); put_line(y);
    Reporting_LU_Newton_Refiner
      (Standard_Output,p,jm,b,v,c,err,rco,res,1.0E-16,1.0E-14,4);
    y := Eval(p,sol,v,t);
    put_line("The solution evaluated at its own basis : "); put_line(y);
  end Solution_Evaluation;

  procedure Silent_Affine_Sampler
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(p,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Sil_Cont(sols,false,target=>Create(1.0));
  end Silent_Affine_Sampler;

  procedure Reporting_Affine_Sampler
              ( file : in file_type;
                p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(p,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Rep_Cont(file,sols,false,target=>Create(1.0));
  end Reporting_Affine_Sampler;

  procedure Silent_LU_Newton_Refiner
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in out Solution_List;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Silent_LU_Newton_Refiner
          (p,jm,b,v,ls.v,ls.err,ls.rco,ls.res,epsxa,epsfa,maxit);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Silent_LU_Newton_Refiner;

  procedure Reporting_LU_Newton_Refiner
              ( file : in file_type;
                p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in out Solution_List;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Reporting_LU_Newton_Refiner
          (file,p,jm,b,v,ls.v,ls.err,ls.rco,ls.res,epsxa,epsfa,maxit);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Reporting_LU_Newton_Refiner;

  procedure Silent_SV_Newton_Refiner
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in out Solution_List;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Silent_SV_Newton_Refiner
          (p,jm,b,v,ls.v,ls.err,ls.rco,ls.res,epsxa,epsfa,maxit);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Silent_SV_Newton_Refiner;

  procedure Reporting_SV_Newton_Refiner
              ( file : in file_type;
                p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in out Solution_List;
                epsxa,epsfa : in double_float; maxit : in natural32 ) is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Reporting_SV_Newton_Refiner
          (file,p,jm,b,v,ls.v,ls.err,ls.rco,ls.res,epsxa,epsfa,maxit);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Reporting_SV_Newton_Refiner;

end Affine_Sampling_Machine;
