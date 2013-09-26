with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;

package body DoblDobl_Root_Refiners is

  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;

    y : DoblDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant double_double := Norm1(A);

  begin
    DoblDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    DoblDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end DoblDobl_Newton_Step;

  procedure DoblDobl_Newton_Step
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Vectors.Vector;
                err,rco,res : out double_double ) is

    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;

    y : DoblDobl_Complex_Vectors.Vector(f'range) := eval(f,x);
    A : Matrix(f'range,f'range) := eval(jf,x);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : constant double_double := Norm1(A);

  begin
    DoblDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    estco(A,A'last(1),ipvt,Anorm,rco);
    lusolve(A,A'last(1),ipvt,y);
    DoblDobl_Complex_Vectors.Add(x,y);
    err := Max_Norm(y);
    y := eval(f,x);
    res := Max_Norm(y);
  end DoblDobl_Newton_Step;

  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      DoblDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      DoblDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure DoblDobl_Root_Refiner
              ( f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution ) is
  begin
    for i in 1..3 loop
      DoblDobl_Newton_Step(f,jf,s.v,s.err,s.rco,s.res);
     -- put("err : "); put(s.err,3);
     -- put(" = rco : "); put(s.rco,3);
     -- put(" = res : "); put(s.res,3); new_line;
    end loop;
  end DoblDobl_Root_Refiner;

  procedure DoblDobl_Root_Refiner
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;

    f : Eval_Laur_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      DoblDobl_Root_Refiner(f,jf,ls);
      tmp := Tail_Of(tmp);
    end loop;
    DoblDobl_Complex_Laur_SysFun.Clear(f);
    DoblDobl_Complex_Laur_JacoMats.Clear(jm);
    DoblDobl_Complex_Laur_JacoMats.Clear(jf);
  end DoblDobl_Root_Refiner;

  procedure Silent_Root_Refiner
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in double_double;
                 numit : in out natural32; max : in natural32 ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    
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

end DoblDobl_Root_Refiners;
