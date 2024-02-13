with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Natural64_VecVecs;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Singular_Values;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with QuadDobl_Solution_Diagnostics;      use QuadDobl_Solution_Diagnostics;
with QuadDobl_Condition_Tables;
with QuadDobl_Condition_Report;
with QuadDobl_Multiple_Solutions;
with QuadDobl_Mixed_Residuals;
with Handle_Underflow_Gracefully;
with Monomial_Hashing;
with QuadDobl_Jacobian_Trees;
with QuadDobl_Deflation_Trees_io;
with QuadDobl_Deflation_Methods;

package body QuadDobl_Root_Refiners is

-- ROOT ACCOUNTING :

  procedure Write_Info
              ( file : in file_type; zero : in Solution;
                initres : in quad_double;
                i,numb,nbdef : in natural32;
                fail,infty : in boolean ) is
  begin
    put(file,"solution "); put(file,i,1); put(file," : ");
    put(file,"   start residual : "); put(file,initres,3);
    if nbdef = 0
     then put(file,"   #iterations : "); put(file,numb,1);
     else put(file,"   #deflations : "); put(file,nbdef,1);
    end if;
    if infty then
      put_line(file,"   at infinity");
    elsif fail then
      put_line(file,"   failure");
    else
      put_line(file,"   success");
    end if;
    put(file,zero);
  end Write_Info;

  procedure Write_Type
              ( file : in file_type; ls : in Link_to_Solution;
                fail,infty : in boolean;
                tolsing : in double_float;
                nbfail,nbinfty : in out natural32;
                nbreal,nbcomp,nbreg,nbsing : in out natural32 ) is
  begin
    if infty then
      put_line(file," at infinity =="); nbinfty := nbinfty + 1;
    elsif fail then
      put_line(file," no solution =="); nbfail := nbfail + 1;
      ls.m := 0;
    else
      if Is_Real(ls.all,1.0E-13)
       then put(file," real ");    nbreal := nbreal + 1;
       else put(file," complex "); nbcomp := nbcomp + 1;
      end if;
      if ls.rco < tolsing
       then put_line(file,"singular =="); nbsing := nbsing + 1;
       else put_line(file,"regular ==");  nbreg := nbreg + 1;
      end if;  
    end if;
  end Write_Type;

  procedure Multiplicity
              ( h1,h2 : in QuadDobl_Complex_Vectors.Vector;
                pl : in out Point_List; ls : in Link_to_Solution;
                nb : in natural32; sols : in out Solution_List;
                fail,infty,deflate : in boolean;
                tolsing,tolclus : in double_float ) is
  begin
    if infty then
      null;
    elsif fail then
      ls.m := 0;
    elsif ls.rco < tolsing or deflate then
      if ls.m <= 1 then -- do not change input multiplicity field
        declare   -- to determine multiplicity count clustered solutions
         -- m : constant natural32 := Multiplicity(ls.all,sols,tolclus);
          m : natural32;
        begin
          QuadDobl_Condition_Report.Multiplicity
            (ls.all,nb,sols,tolclus,h1,h2,pl,m);
          if ((m = 1) and (not deflate))
           then ls.m := 0;
           else ls.m := integer32(m);
          end if;
        end;
      end if;
    else  -- ls.rco > tolsing, check for clustering
      declare
       -- nb2 : constant natural32 := Is_Clustered(ls.all,nb,sols,tolclus);
        nb2 : natural32;
      begin
        QuadDobl_Condition_Report.Is_Clustered
          (ls.all,nb,sols,tolclus,h1,h2,pl,nb2);
        if nb2 /= nb then
          ls.m := -integer32(nb2);
          Change_Multiplicity(sols,nb2,-integer32(nb));
        end if;
      end;	   
    end if;
  end Multiplicity;

  procedure Write_Global_Info
              ( file : in file_type; tot,nbfail,nbinfty,
                nbreal,nbcomp,nbreg,nbsing,nbclus : in natural32 ) is
  begin
    Standard_Complex_Solutions_io.put_bar(file);
    put(file,"A list of "); put(file,tot,1);
    put_line(file," solutions has been refined :");
    put(file,"Number of regular solutions     : "); put(file,nbreg,1);
    put_line(file,".");
    put(file,"Number of singular solutions    : "); put(file,nbsing,1);
    put_line(file,".");
    put(file,"Number of real solutions        : "); put(file,nbreal,1);
    put_line(file,".");
    put(file,"Number of complex solutions     : "); put(file,nbcomp,1);
    put_line(file,".");
    put(file,"Number of clustered solutions   : "); put(file,nbclus,1);
    put_line(file,".");
    put(file,"Number of solutions at infinity : "); put(file,nbinfty,1);
    put_line(file,".");
    put(file,"Number of failures              : "); put(file,nbfail,1);
    put_line(file,".");
    Standard_Complex_Solutions_io.put_bar(file);
  end Write_Global_Info;

-- ONE NEWTON STEP :

  procedure Write_Diagnostics
              ( file : in file_type; step : natural32;
                err,rco,res : in quad_double ) is
  begin
    put(file,"Step "); put(file,step,4); put(file," : ");
    put(file," |errxa| : "); put(file,err,3);
    put(file," est rco : "); put(file,rco,3);
    put(file," |errfa| : "); put(file,res,3); new_line(file);
  end Write_Diagnostics;

  procedure QuadDobl_SVD_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,x'range) := eval(jf,x);
    n : constant integer32 := f'last;
    p : constant integer32 := x'last;
    dx : QuadDobl_Complex_Vectors.Vector(1..p);
    mm : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    sv : QuadDobl_Complex_Vectors.Vector(1..mm);
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_SVD_Newton_Step 1 ...");
    end if;
    QuadDobl_Complex_Singular_Values.SVD(A,n,p,sv,e,u,v,job,info);
    rco := QuadDobl_Complex_Singular_Values.Inverse_Condition_Number(sv);
    QuadDobl_Complex_Vectors.Min(y);
    dx := QuadDobl_Complex_Singular_Values.Solve(u,v,sv,y);
    QuadDobl_Complex_Vectors.Add(x,dx);
    err := Max_Norm(dx);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_SVD_Newton_Step;

  procedure QuadDobl_SVD_Newton_Step
              ( f,abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,x'range) := eval(jf,x);
    n : constant integer32 := f'last;
    p : constant integer32 := x'last;
    dx : QuadDobl_Complex_Vectors.Vector(1..p);
    mm : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    sv : QuadDobl_Complex_Vectors.Vector(1..mm);
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);
    ay : QuadDobl_Complex_Vectors.Vector(abh'range);
    qd_one : constant quad_double := create(1.0);
    one : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Complex_Numbers.Create(qd_one);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_SVD_Newton_Step 2 ...");
    end if;
    QuadDobl_Complex_Singular_Values.SVD(A,n,p,sv,e,u,v,job,info);
    rco := QuadDobl_Complex_Singular_Values.Inverse_Condition_Number(sv);
    QuadDobl_Complex_Vectors.Min(y);
    dx := QuadDobl_Complex_Singular_Values.Solve(u,v,sv,y);
    QuadDobl_Complex_Vectors.Add(x,dx);
    err := Max_Norm(dx);
    y := eval(f,x); -- res := Max_Norm(y);
    xt(x'range) := QuadDobl_Mixed_Residuals.AbsVal(x);
    xt(xt'last) := one;
    ay := Eval(abh,xt);
    res := QuadDobl_Mixed_Residuals.Residual(y,ay);
  end QuadDobl_SVD_Newton_Step;

  procedure QuadDobl_SVD_Newton_Step
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,x'range) := eval(jf,x);
    n : constant integer32 := f'last;
    p : constant integer32 := x'last;
    dx : QuadDobl_Complex_Vectors.Vector(1..p);
    mm : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    sv : QuadDobl_Complex_Vectors.Vector(1..mm);
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_SVD_Newton_Step 3 ...");
    end if;
    QuadDobl_Complex_Singular_Values.SVD(A,n,p,sv,e,u,v,job,info);
    rco := QuadDobl_Complex_Singular_Values.Inverse_Condition_Number(sv);
    QuadDobl_Complex_Vectors.Min(y);
    dx := QuadDobl_Complex_Singular_Values.Solve(u,v,sv,y);
    QuadDobl_Complex_Vectors.Add(x,dx);
    err := Max_Norm(dx);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_SVD_Newton_Step;

  procedure QuadDobl_LU_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_LU_Newton_Step 1 ...");
    end if;
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_LU_Newton_Step;

  procedure QuadDobl_LU_Newton_Step
              ( f,abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);
    xt : QuadDobl_Complex_Vectors.Vector(x'first..x'last+1);
    ay : QuadDobl_Complex_Vectors.Vector(abh'range);
    qd_one : constant quad_double := create(1.0);
    one : constant QuadDobl_Complex_Numbers.Complex_Number
        := QuadDobl_Complex_Numbers.Create(qd_one);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_LU_Newton_Step 2 ...");
    end if;
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x); -- res := Max_Norm(y);
    xt(x'range) := QuadDobl_Mixed_Residuals.AbsVal(x);
    xt(xt'last) := one;
    ay := Eval(abh,xt);
    res := QuadDobl_Mixed_Residuals.Residual(y,ay);
  end QuadDobl_LU_Newton_Step;

  procedure QuadDobl_LU_Newton_Step
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys; 
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_LU_Newton_Step 3 ...");
    end if;
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_LU_Newton_Step;

  procedure QuadDobl_SVD_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Vectors.Vector;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Jacobian_Circuits;

    y : QuadDobl_Complex_Vectors.Vector(f'range);
    A : Matrix(f'range,x'range);
    n : constant integer32 := f'last;
    p : constant integer32 := x'last;
    dx : QuadDobl_Complex_Vectors.Vector(1..p);
    mm : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    sv : QuadDobl_Complex_Vectors.Vector(1..mm);
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_SVD_Newton_Step 4 ...");
    end if;
    EvalDiff(jf,x,wrk,y,A);
    QuadDobl_Complex_Singular_Values.SVD(A,n,p,sv,e,u,v,job,info);
    rco := QuadDobl_Complex_Singular_Values.Inverse_Condition_Number(sv);
    QuadDobl_Complex_Vectors.Min(y);
    dx := QuadDobl_Complex_Singular_Values.Solve(u,v,sv,y);
    QuadDobl_Complex_Vectors.Add(x,dx);
    err := Max_Norm(dx);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_SVD_Newton_Step;

  procedure QuadDobl_LU_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Vectors.Vector;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Jacobian_Circuits;

    y : QuadDobl_Complex_Vectors.Vector(f'range);
    A : Matrix(f'range,x'range);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : quad_double;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_LU_Newton_Step 4 ...");
    end if;
    EvalDiff(jf,x,wrk,y,A);
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    Anorm := Norm1(A);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_LU_Newton_Step;

-- WRAPPING ONE NEWTON STEP :

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is
  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_Newton_Step 1 ...");
    end if;
    if f'last > x'last
     then QuadDobl_SVD_Newton_Step(f,jf,x,err,rco,res,verbose-1);
     else QuadDobl_LU_Newton_Step(f,jf,x,err,rco,res,verbose-1);
    end if;
  end QuadDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( f,abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is
  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_Newton_Step 2 ...");
    end if;
    if f'last > x'last
     then QuadDobl_SVD_Newton_Step(f,abh,jf,x,err,rco,res,verbose-1);
     else QuadDobl_LU_Newton_Step(f,abh,jf,x,err,rco,res,verbose-1);
    end if;
  end QuadDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is
  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_Newton_Step 3 ...");
    end if;
    if f'last > x'last
     then QuadDobl_SVD_Newton_Step(f,jf,x,err,rco,res,verbose-1);
     else QuadDobl_LU_Newton_Step(f,jf,x,err,rco,res,verbose-1);
    end if;
  end QuadDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Vectors.Vector;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                err,rco,res : out quad_double;
                verbose : in integer32 := 0 ) is
  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("QuadDobl_Newton_Step 4 ...");
    end if;
    if f'last > x'last
     then QuadDobl_SVD_Newton_Step(f,jf,x,wrk,err,rco,res,verbose-1);
     else QuadDobl_LU_Newton_Step(f,jf,x,wrk,err,rco,res,verbose-1);
    end if;
  end QuadDobl_Newton_Step;

-- SEVERAL NEWTON STEPS :

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure Reporting_Newton
              ( file : in file_type;
                f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure Reporting_Newton
              ( file : in file_type;
                f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Solutions.Solution;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,jf,x.v,wrk,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

-- MIXED RESIDUAL VERSIONS :

  procedure Silent_Newton
              ( f,abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,abh,jf,x.v,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure Reporting_Newton
              ( file : in file_type;
                f,abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_float; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      QuadDobl_Newton_Step(f,abh,jf,x.v,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

-- REFINING A LIST OF SOLUTIONS :

  procedure QuadDobl_Root_Refiner
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution ) is
  begin
    for i in 1..5 loop
      QuadDobl_Newton_Step(f,jf,s.v,s.err,s.rco,s.res);
     -- put("err : "); put(s.err,3);
     -- put(" = rco : "); put(s.rco,3);
     -- put(" = res : "); put(s.res,3); new_line;
    end loop;
  end QuadDobl_Root_Refiner;

  procedure QuadDobl_Root_Refiner
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Jacobian_Circuits.Circuit;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in 1..5 loop
      QuadDobl_Newton_Step(f,jf,s.v,wrk,s.err,s.rco,s.res);
     -- put("err : "); put(s.err,3);
     -- put(" = rco : "); put(s.rco,3);
     -- put(" = res : "); put(s.res,3); new_line;
    end loop;
  end QuadDobl_Root_Refiner;

  procedure QuadDobl_Root_Refiner
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                s : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;

    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      QuadDobl_Root_Refiner(f,jf,ls);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(f); Clear(jm); Clear(jf);
  end QuadDobl_Root_Refiner;

  procedure QuadDobl_Root_Refiner
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Jacobian_Circuits;

    f : Eval_Poly_Sys(p'range) := Create(p);
    jf : Circuit := Create(p);
    nm : constant integer32 := integer32(Number_of_Monomials(jf));
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..nm) := WorkSpace(jf);
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      QuadDobl_Root_Refiner(f,jf,ls,wrk);
      tmp := Tail_Of(tmp);
    end loop;
    QuadDobl_Complex_Poly_SysFun.Clear(f);
    QuadDobl_Jacobian_Circuits.Clear(jf);
    QuadDobl_Complex_VecVecs.Clear(wrk);
  end QuadDobl_Root_Refiner;

-- THE MAIN ROOT REFINERS :

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    nb : natural32;
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Silent_Root_Refiner 1 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
        numit := numit + nb;
      else
        fail := true;
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    refs_last : Solution_List;
    nb : natural32;
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Silent_Root_Refiner 2 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
        numit := numit + nb;
        if not fail
         then Append(refs,refs_last,ls.all);
        end if;
      else
        fail := true;
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    nb : natural32;
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Silent_Root_Refiner 3 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
        numit := numit + nb;
      else
        fail := true;
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    refs_last : Solution_List;
    nb : natural32;
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Silent_Root_Refiner 4 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
        numit := numit + nb;
        if not fail
         then Append(refs,refs_last,ls.all);
        end if;
      else
        fail := true;
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    len : constant natural32 := Length_Of(s);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Reporting_Root_Refiner 1 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," "); put(file,nv,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
      Write_Info(file,ls.all,initres,cnt,nb,0,fail,infty);
      Write_Type(file,ls,fail,infty,tolsing,nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.all);
      solsptr := Tail_Of(solsptr);
    end loop;
    Write_Global_Info
      (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    refs_last : Solution_List;
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    len : constant natural32 := Length_Of(s);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Reporting_Root_Refiner 2 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," "); put(file,nv,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
      Write_Info(file,ls.all,initres,cnt,nb,0,fail,infty);
      Write_Type(file,ls,fail,infty,tolsing,nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.all);
      if not fail
       then Append(refs,refs_last,ls.all);
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Write_Global_Info
      (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    len : constant natural32 := Length_Of(s);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Reporting_Root_Refiner 3 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," "); put(file,nv,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
      Write_Info(file,ls.all,initres,cnt,nb,0,fail,infty);
      Write_Type(file,ls,fail,infty,tolsing,nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.all);
      solsptr := Tail_Of(solsptr);
    end loop;
    Write_Global_Info
      (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    
    nv : constant integer32 := Head_Of(s).n;
    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nv) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nv) := Create(jm);
    refs_last : Solution_List;
    nb : natural32;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..60)
                      := QuadDobl_Condition_Tables.Create(60); 
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nv);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    len : constant natural32 := Length_Of(s);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Reporting_Root_Refiner 4 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," "); put(file,nv,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if wout
         then Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
         else Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,false,tolsing,epsxa);
      Write_Info(file,ls.all,initres,cnt,nb,0,fail,infty);
      Write_Type(file,ls,fail,infty,tolsing,nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.all);
      if not fail
       then Append(refs,refs_last,ls.all);
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Write_Global_Info
      (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    Clear(pl); Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

-- REFINEMENT WITH DEFLATION :

  procedure Silent_Deflate
              ( max : in natural32;
                f : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out Link_to_Solution; order : in integer32;
                tolrnk : in double_float;
                nd : in QuadDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                nv,nq,R1 : in out Standard_Natural_Vectors.Vector;
                nbitr,nbdef : out natural32; fail : out boolean;
                verbose : in integer32 := 0 ) is

  -- DECRIPTION :
  --   Performs deflation on a system without intermediate output.

    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Deflation_Methods;

    k : integer32 := 0;
    rank,m : natural32;
    done : boolean := false;
    B : QuadDobl_Complex_VecMats.VecMat(1..order);
    h : QuadDobl_Complex_VecVecs.VecVec(1..order);
    x : QuadDobl_Complex_VecVecs.VecVec(0..order);
    rsd : quad_double;

  begin
    if verbose > 0
     then put_line("-> in quaddobl_root_refiners.Silent_Deflate ...");
    end if;
    nbitr := 0;
    nbdef := 0;
    loop
      if k = 0 then
        Apply_Newton(max,f,jf,ls,tolrnk,rank);
        done := (integer32(rank) = ls.n);
        if not done then
          x(0) := new QuadDobl_Complex_Vectors.Vector'(ls.v);
        end if;
        nbitr := 1;
      else
        Apply_Newton(max,f,jf,nd,monkeys,k,nv,nq,
                     R1,B,h,x,ls.err,ls.rco,ls.res,tolrnk,rank);
        ls.v := x(0).all;
        done := (rank = nv(k));
      end if;
      exit when done or (k >= order);
      k := k + 1; m := rank + 1;
      Add_Deflation_Data(k,m,nv,nq,R1,B,h);
      Add_Multipliers(f,jf,nd,monkeys,k,nv(0..k),
                      nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
      if rsd > 0.1
       then done := true; -- do not report failure if deflation fails!
      end if;
      exit when (rsd > 0.1); -- deflation fails
    end loop;
    QuadDobl_Complex_Vectors.Clear(x(0));
    for i in 1..k loop
      QuadDobl_Complex_Matrices.Clear(B(i)); 
      QuadDobl_Complex_Vectors.Clear(h(i));
      QuadDobl_Complex_Vectors.Clear(x(i));
    end loop;
    fail := not done;
    nbdef := natural32(k);
  exception
   -- when constraint_error -- same underflow as in Reporting_Newton may occur
    when others -- 
       => -- put_line("exception raised in reporting deflate");
          Handle_Underflow_Gracefully.Underflow_to_Zero(ls.v);
          ls.rco := Handle_Underflow_Gracefully.Underflow_to_Zero(ls.rco);
          fail := not done;
  end Silent_Deflate;

  procedure Reporting_Deflate
              ( file : in file_type; output : in boolean; 
                max : in natural32;
                f : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                ls : in out Link_to_Solution; order : in integer32;
                tolrnk : in double_float;
                nd : in QuadDobl_Jacobian_Trees.Link_to_Eval_Node;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                nv,nq,R1 : in out Standard_Natural_Vectors.Vector;
                nbitr,nbdef : out natural32; fail : out boolean;
                verbose : in integer32 := 0 ) is

  -- DECRIPTION :
  --   Performs deflation on a system with intermediate output.

    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Deflation_Methods;
    use QuadDobl_Deflation_Trees_io;

    k : integer32 := 0;
    rank,m : natural32;
    done : boolean := false;
    B : QuadDobl_Complex_VecMats.VecMat(1..order);
    h : QuadDobl_Complex_VecVecs.VecVec(1..order);
    x : QuadDobl_Complex_VecVecs.VecVec(0..order);
    rsd : quad_double;

  begin
    if verbose > 0
     then put_line("-> in quaddobl_root_refiners.Reporting_Deflate ...");
    end if;
    nbitr := 0;
    nbdef := 0;
    loop
      if k = 0 then
        Apply_Newton(file,output,max,f,jf,ls,tolrnk,rank);
        done := (integer32(rank) = ls.n);
        if not done then
          x(0) := new QuadDobl_Complex_Vectors.Vector'(ls.v);
        end if;
        nbitr := 1;
      else
        Apply_Newton(file,output,max,f,jf,nd,monkeys,k,nv,nq,
                     R1,B,h,x,ls.err,ls.rco,ls.res,tolrnk,rank);
        ls.v := x(0).all;
        done := (rank = nv(k));
      end if;
      exit when done or (k >= order);
      k := k + 1; m := rank + 1;
      Add_Multiplier_Symbols(natural32(k),m);
      Add_Deflation_Data(k,m,nv,nq,R1,B,h);
      Add_Multipliers(file,output,f,jf,nd,monkeys,k,nv(0..k),
                      nq(0..k),R1(1..k),B(1..k),h(1..k),x(0..k),rsd);
      if rsd > 0.1
       then done := true; -- do not report failure if deflation fails!
      end if;
      exit when (rsd > 0.1); -- deflation fails
    end loop;
    QuadDobl_Complex_Vectors.Clear(x(0));
    for i in 1..k loop
      QuadDobl_Complex_Matrices.Clear(B(i)); 
      QuadDobl_Complex_Vectors.Clear(h(i));
      QuadDobl_Complex_Vectors.Clear(x(i));
    end loop;
    fail := not done;
    nbdef := natural32(k);
  exception
   -- when constraint_error -- same underflow as in Reporting_Newton may occur
    when others -- 
       => -- put_line("exception raised in reporting deflate");
          Handle_Underflow_Gracefully.Underflow_to_Zero(ls.v);
          ls.rco := Handle_Underflow_Gracefully.Underflow_to_Zero(ls.rco);
          fail := not done;
  end Reporting_Deflate;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Multiple_Solutions;
    use Monomial_Hashing;
    use QuadDobl_Jacobian_Trees;
    
    nbequ : constant integer32 := p'last;
    nbvar : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nbvar) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nbvar) := Create(jm);
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(nbvar);
    merge : boolean := false; -- to merge clustered solutions
    nb,numb,nbdef : natural32 := 0;
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nbvar);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Silent_Root_Refiner 5 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    if deflate then
      declare
      begin
        nv(0) := natural32(nbvar); nq(0) := natural32(nbequ); R1(1) := 0;
        Create_Remember_Derivatives(jm,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if deflate then
          backup := ls.all;
          Silent_Deflate(max,f,jf,ls,order,tolrnk,nd,monkeys,
                         nv,nq,R1,numb,nbdef,fail,verbose-1);
         -- reinstate Newton after deflation
          if (ls.err > epsxa) and (ls.res > epsfa) then
            Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            if fail and backup.res < ls.res then
              ls.all := backup;
              Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
          end if;
        else
          Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,deflate,tolsing,epsxa);
      if not fail and deflate
       then merge := merge and (ls.m > 1);
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      QuadDobl_Jacobian_Trees.Clear(nd);
      if merge then
        Merge_Multiple_Solutions(s,tolsing);
      end if;
    else
      Clear(jm); -- otherwise crash after Clear(nd)
    end if;
    Clear(pl); Clear(f); Clear(jf);
  end Silent_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Multiple_Solutions;
    use Monomial_Hashing;
    use QuadDobl_Jacobian_Trees;
    
    nbequ : constant integer32 := p'last;
    nbvar : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nbvar) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nbvar) := Create(jm);
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(nbvar);
    merge : boolean := false; -- to merge clustered solutions
    refs_last : Solution_List;
    nb,numb,nbdef : natural32 := 0;
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nbvar);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Silent_Root_Refiner 6 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    if deflate then
      declare
      begin
        nv(0) := natural32(nbvar); nq(0) := natural32(nbequ); R1(1) := 0;
        Create_Remember_Derivatives(jm,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if deflate then
          backup := ls.all;
          Silent_Deflate(max,f,jf,ls,order,tolrnk,nd,monkeys,
                         nv,nq,R1,numb,nbdef,fail,verbose-1);
         -- reinstate Newton after deflation
          if (ls.err > epsxa) and (ls.res > epsfa) then
            Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            if fail and backup.res < ls.res then
              ls.all := backup;
              Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
          end if;
        else
          Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,deflate,tolsing,epsxa);
      if not fail then
        Append(refs,refs_last,ls.all);
        if deflate
         then merge := merge and (ls.m > 1);
        end if;
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      QuadDobl_Jacobian_Trees.Clear(nd);
      if merge then
        Merge_Multiple_Solutions(s,tolsing);
        Merge_Multiple_Solutions(refs,tolsing);
      end if;
    else
      Clear(jm); -- otherwise crash after Clear(nd)
    end if;
    Clear(pl); Clear(f); Clear(jf);
  end Silent_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Multiple_Solutions;
    use Monomial_Hashing;
    use QuadDobl_Jacobian_Trees;
    
    nbequ : constant integer32 := p'last;
    nbvar : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nbvar) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nbvar) := Create(jm);
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(nbvar);
    merge : boolean := false; -- to merge clustered solutions
    nb,numb,nbdef : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..30)
                      := QuadDobl_Condition_Tables.Create(30); 
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nbvar);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    len : constant natural32 := Length_Of(s);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Reporting_Root_Refiner 5 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    if deflate then
      declare
      begin
        nv(0) := natural32(nbvar); nq(0) := natural32(nbequ); R1(1) := 0;
        Create_Remember_Derivatives(jm,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," "); put(file,nbvar,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if deflate then
          backup := ls.all;
          Reporting_Deflate(file,wout,max,f,jf,ls,order,tolrnk,nd,monkeys,
                            nv,nq,R1,numb,nbdef,fail,verbose-1);
         -- reinstate Newton after deflation
          if (ls.err > epsxa) and (ls.res > epsfa) then
            if wout
             then Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
             else Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
            if fail and backup.res < ls.res then
              ls.all := backup;
              Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
          end if;
        elsif wout then
          Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        else
          Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,deflate,tolsing,epsxa);
      Write_Info(file,ls.all,initres,cnt,nb,nbdef,fail,infty);
      Write_Type(file,ls,fail,infty,tolsing,nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.all);
      if not fail and deflate
       then merge := merge and (ls.m > 1);
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Write_Global_Info
      (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      QuadDobl_Jacobian_Trees.Clear(nd);
      if merge then
        Merge_Multiple_Solutions(s,tolsing);
      end if;
    else
      Clear(jm); -- otherwise crash after Clear(nd)
    end if;
    Clear(pl); Clear(f); Clear(jf);
  end Reporting_Root_Refiner;

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s,refs : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Multiple_Solutions;
    use Monomial_Hashing;
    use QuadDobl_Jacobian_Trees;
    
    nbequ : constant integer32 := p'last;
    nbvar : constant integer32 := Head_Of(s).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nbvar) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nbvar) := Create(jm);
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(nbvar);
    merge : boolean := false; -- to merge clustered solutions
    refs_last : Solution_List;
    nb,numb,nbdef : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..30)
                      := QuadDobl_Condition_Tables.Create(30); 
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nbvar);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := s;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    len : constant natural32 := Length_Of(s);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Reporting_Root_Refiner 6 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    if deflate then
      declare
      begin
        nv(0) := natural32(nbvar); nq(0) := natural32(nbequ); R1(1) := 0;
        Create_Remember_Derivatives(jm,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," "); put(file,nbvar,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      ls.res := Sum_Norm(Eval(f,ls.v));
      initres := ls.res;
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if deflate then
          backup := ls.all;
          Reporting_Deflate(file,wout,max,f,jf,ls,order,tolrnk,nd,monkeys,
                            nv,nq,R1,numb,nbdef,fail,verbose-1);
         -- reinstate Newton after deflation
          if (ls.err > epsxa) and (ls.res > epsfa) then
            if wout
             then Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
             else Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
            if fail and backup.res < ls.res then
              ls.all := backup;
              Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
          end if;
        elsif wout then
          Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        else
          Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,s,fail,infty,deflate,tolsing,epsxa);
      Write_Info(file,ls.all,initres,cnt,nb,nbdef,fail,infty);
      Write_Type(file,ls,fail,infty,tolsing,nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.all);
      if not fail then
        Append(refs,refs_last,ls.all);
        if deflate
         then merge := merge and (ls.m > 1);
        end if;
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Write_Global_Info
      (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      QuadDobl_Jacobian_Trees.Clear(nd);
      if merge then
        Merge_Multiple_Solutions(s,tolsing);
        Merge_Multiple_Solutions(refs,tolsing);
      end if;
    else
      Clear(jm); -- otherwise crash after Clear(nd)
    end if;
    Clear(pl); Clear(f); Clear(jf);
  end Reporting_Root_Refiner;

-- REFINEMENT with mixed residuals :

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 abh : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                 sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural32; max : in natural32;
                 deflate : in out boolean; wout : in boolean;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Multiple_Solutions;
    use Monomial_Hashing;
    use QuadDobl_Jacobian_Trees;
    
    nbequ : constant integer32 := p'last;
    nbvar : constant integer32 := Head_Of(sols).n;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,1..nbvar) := Create(p);
    jf : Eval_Jaco_Mat(p'range,1..nbvar) := Create(jm);
    tolrnk : constant double_float := 100.0*tolsing;
    order : constant integer32 := 3;
    nv : Standard_Natural_Vectors.Vector(0..order);
    nq : Standard_Natural_Vectors.Vector(0..order);
    R1 : Standard_Natural_Vectors.Vector(1..order);
    monkeys : Standard_Natural64_VecVecs.VecVec(1..order);
    nd : Link_to_Eval_Node;
    backup : Solution(nbvar);
    merge : boolean := false; -- to merge clustered solutions
    nb,numb,nbdef : natural32 := 0;
    nbfail,nbinfty,nbreg,nbsing,nbclus,nbreal,nbcomp : natural32 := 0;
    t_err,t_rco,t_res : Standard_Natural_Vectors.Vector(0..30)
                      := QuadDobl_Condition_Tables.Create(30); 
    fail,infty : boolean;
    seed : integer32 := 1234567;
    h1,h2 : QuadDobl_Complex_Vectors.Vector(1..nbvar);
    pl : Point_List;
    initres : quad_double;
    solsptr : Solution_List := sols;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    len : constant natural32 := Length_Of(sols);

  begin
    if verbose > 0 then
      put("-> in quaddobl_root_refiners.");
      put_line("Reporting_Root_Refiner 7 ...");
    end if;
    QuadDobl_Random_Vectors.Random_Vector(seed,h1);
    QuadDobl_Random_Vectors.Random_Vector(seed,h2);
    if deflate then
      declare
      begin
        nv(0) := natural32(nbvar); nq(0) := natural32(nbequ); R1(1) := 0;
        Create_Remember_Derivatives(jm,order,nd);
        monkeys := Monomial_Keys(natural32(order),nv(0));
      exception
        when others =>
          put_line("The system is too large to apply deflation.");
          deflate := false;
      end;
    end if;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,len,1); put(file," "); put(file,nbvar,1); new_line(file);
    Standard_Complex_Solutions_io.put_bar(file);
    while not Is_Null(solsptr) loop
      ls := Head_Of(solsptr);
      nb := 0; cnt := cnt + 1;
      initres := ls.res; -- ls.res := Sum_Norm(Eval(f,ls.v));
      infty := At_Infinity(ls.all,false,1.0E+8);
      if not infty and ls.res < 0.1 and ls.err < 0.1 then
        if deflate then
          backup := ls.all;
          Reporting_Deflate(file,wout,max,f,jf,ls,order,tolrnk,nd,monkeys,
                            nv,nq,R1,numb,nbdef,fail,verbose-1);
         -- reinstate Newton after deflation
          if (ls.err > epsxa) and (ls.res > epsfa) then
            if wout then
              Reporting_Newton(file,f,abh,jf,ls.all,epsxa,epsfa,nb,max,fail);
            else
              Silent_Newton(f,abh,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
            if fail and backup.res < ls.res then
              ls.all := backup;
              Silent_Newton(f,abh,jf,ls.all,epsxa,epsfa,nb,max,fail);
            end if;
          end if;
        elsif wout then
          Reporting_Newton(file,f,abh,jf,ls.all,epsxa,epsfa,nb,max,fail);
        else
          Silent_Newton(f,abh,jf,ls.all,epsxa,epsfa,nb,max,fail);
        end if;
        numit := numit + nb;
      else
        fail := true;
      end if;
      Multiplicity(h1,h2,pl,ls,cnt,sols,fail,infty,deflate,tolsing,epsxa);
      Write_Info(file,ls.all,initres,cnt,nb,nbdef,fail,infty);
      Write_Type(file,ls,fail,infty,tolsing,nbfail,nbinfty,
                 nbreal,nbcomp,nbreg,nbsing);
      QuadDobl_Condition_Tables.Update_Corrector(t_err,ls.all);
      QuadDobl_Condition_Tables.Update_Condition(t_rco,ls.all);
      QuadDobl_Condition_Tables.Update_Residuals(t_res,ls.all);
      if not fail and deflate
       then merge := merge and (ls.m > 1);
      end if;
      solsptr := Tail_Of(solsptr);
    end loop;
    Write_Global_Info
      (file,len,nbfail,nbinfty,nbreal,nbcomp,nbreg,nbsing,nbclus);
    QuadDobl_Condition_Tables.Write_Tables(file,t_err,t_res,t_rco);
    if deflate then
      Standard_Natural64_VecVecs.Clear(monkeys);
      QuadDobl_Jacobian_Trees.Clear(nd);
      if merge then
        Merge_Multiple_Solutions(sols,tolsing);
      end if;
    else
      Clear(jm); -- otherwise crash after Clear(nd)
    end if;
    Clear(pl); Clear(f); Clear(jf);
  end Reporting_Root_Refiner;

end QuadDobl_Root_Refiners;
