with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Matrices;          use DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;

package body DoblDobl_Root_Refiners is

  procedure Write_Diagnostics
              ( file : in file_type; step : natural32;
                err,rco,res : in double_double ) is

  -- DESCRIPTION :
  --   Writes diagnostics in one line about the Newton step to file.
  --   This procedure defines the formatting in the Reporting
  --   versions of a sequence of Newton steps.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   step     the current step number;
  --   err      forward error, magnitude of the update;
  --   rco      estimate for the inverse of the condition number;
  --   res      backward error, residual.

  begin
    put(file,"Step "); put(file,step,4); put(file," : ");
    put(file," |errxa| : "); put(file,err,3);
    put(file," est rco : "); put(file,rco,3);
    put(file," |errfa| : "); put(file,res,3); new_line(file);
  end Write_Diagnostics;

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
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Vectors.Vector;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                err,rco,res : out double_double ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Jacobian_Circuits;

    y : DoblDobl_Complex_Vectors.Vector(f'range);
    A : Matrix(f'range,f'range);
    ipvt : Standard_Integer_Vectors.Vector(A'range(2));
    info : integer32;
    Anorm : double_double;

  begin
    EvalDiff(jf,x,wrk,y,A);
    DoblDobl_Complex_Vectors.Min(y);
    lufac(A,A'last(1),ipvt,info);
    Anorm := Norm1(A);
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

  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      DoblDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

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

  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jf : in  DoblDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                x : in out DoblDobl_Complex_Solutions.Solution;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      DoblDobl_Newton_Step(f,jf,x.v,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

  procedure Silent_Newton
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Solutions.Solution;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      DoblDobl_Newton_Step(f,jf,x.v,wrk,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Silent_Newton;

  procedure Reporting_Newton
              ( file : in file_type;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in  DoblDobl_Jacobian_Circuits.Circuit;
                x : in out DoblDobl_Complex_Solutions.Solution;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                epsxa,epsfa : in double_double; numit : in out natural32;
                max : in natural32; fail : out boolean ) is
  begin
    fail := true;
    while numit < max loop
      numit := numit + 1;
      DoblDobl_Newton_Step(f,jf,x.v,wrk,x.err,x.rco,x.res);
      Write_Diagnostics(file,numit,x.err,x.rco,x.res);
      if (x.err < epsxa) or (x.res < epsfa)
       then fail := false; exit;
      end if;
    end loop;
  end Reporting_Newton;

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
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Jacobian_Circuits.Circuit;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec ) is
  begin
    for i in 1..3 loop
      DoblDobl_Newton_Step(f,jf,s.v,wrk,s.err,s.rco,s.res);
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

  procedure DoblDobl_Root_Refiner
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Jacobian_Circuits;
    use DoblDobl_Complex_Solutions;

    f : Eval_Poly_Sys(p'range) := Create(p);
    jf : Circuit := Create(p);
    nm : constant integer32 := integer32(Number_of_Monomials(jf));
    wrk : DoblDobl_Complex_VecVecs.VecVec(1..nm) := WorkSpace(jf);
    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      DoblDobl_Root_Refiner(f,jf,ls,wrk);
      tmp := Tail_Of(tmp);
    end loop;
    DoblDobl_Complex_Poly_SysFun.Clear(f);
    DoblDobl_Jacobian_Circuits.Clear(jf);
    DoblDobl_Complex_VecVecs.Clear(wrk);
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

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 s : in out DoblDobl_Complex_Solutions.Solution_List;
                 epsxa,epsfa : in double_double;
                 numit : in out natural32; max : in natural32;
                 wout : in boolean ) is

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
      if wout
       then Reporting_Newton(file,f,jf,ls.all,epsxa,epsfa,nb,max,fail);
       else Silent_Newton(f,jf,ls.all,epsxa,epsfa,nb,max,fail);
      end if;
      numit := numit + nb;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(f); Clear(jm); Clear(jf);
  end Reporting_Root_Refiner;

end DoblDobl_Root_Refiners;
