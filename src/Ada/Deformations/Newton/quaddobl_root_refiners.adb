with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;

package body QuadDobl_Root_Refiners is

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);

  begin
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys; 
                jf : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                err,rco,res : out quad_double ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;

    y : QuadDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant quad_double := Norm1(A);

  begin
    QuadDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    QuadDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end QuadDobl_Newton_Step;

  procedure QuadDobl_Newton_Step
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys; 
                jf : in QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Vectors.Vector;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                err,rco,res : out quad_double ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Jacobian_Circuits;

    y : QuadDobl_Complex_Vectors.Vector(f'range);
    A : Matrix(f'range,f'range);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : quad_double;

  begin
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
  end QuadDobl_Newton_Step;

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in quad_double; numit : in out natural32;
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

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in quad_double; numit : in out natural32;
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

  procedure Silent_Newton
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  QuadDobl_Jacobian_Circuits.Circuit;
                x : in out QuadDobl_Complex_Solutions.Solution;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in quad_double; numit : in out natural32;
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
    use QuadDobl_Complex_Solutions;

    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
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
    use QuadDobl_Complex_Solutions;

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

  procedure Silent_Root_Refiner
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out QuadDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in quad_double;
                 numit : in out natural32; max : in natural32 ) is

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    tmp : Solution_List := s;
    ls : Link_to_Solution;
    nb : natural32;
    fail : boolean;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      nb := 0;
      Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
      numit := numit + nb;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(f); Clear(jm); Clear(jf);
  end Silent_Root_Refiner;

end QuadDobl_Root_Refiners;
