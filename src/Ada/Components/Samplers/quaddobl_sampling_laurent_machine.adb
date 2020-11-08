with Numbers_io;                         use Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Matrices;          use QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Continuation_Parameters;
with Continuation_Parameters_io;
with QuadDobl_IncFix_Continuation;       use QuadDobl_IncFix_Continuation;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with QuadDobl_Root_Refiners;
with Planes_and_Polynomials;             use Planes_and_Polynomials;

package body QuadDobl_Sampling_Laurent_Machine is

-- INTERNAL DATA :

  qdsys : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
  qdsys_eval : QuadDobl_Complex_Laur_SysFun.Link_to_Eval_Laur_Sys;
  qdjac_eval : QuadDobl_Complex_Laur_JacoMats.Link_to_Eval_Jaco_Mat;

-- PARAMETERS FOR ROOT REFINER :

  epsxa,epsfa,tolsing : double_float;
  maxit : natural32;

-- AUXILIARY ROUTINES :

  function Number_of_Equations return integer32 is

  -- DESCRIPTION :
  --   Returns the number of equations in the embedded system.

  begin
    return qdsys'last;
  end Number_of_Equations;

  function QuadDobl_System return QuadDobl_Complex_Laur_Systems.Laur_Sys is

  -- DESCRIPTION :
  --   Returns the embedded polynomial system.

  begin
    return qdsys.all;
  end QuadDobl_System;

  function Embed ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                   h : in QuadDobl_Complex_VecVecs.VecVec ) 
                 return QuadDobl_Complex_Laur_Systems.Laur_Sys is

  -- DESCRIPTION :
  --   Returns the embedding of p with other slices in h.

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);
    k : constant integer32 := h'last;

  begin
    res(p'first..p'last-k) := p(p'first..p'last-k);
    for i in 1..k loop
      res(p'last-k+i) := Hyperplane(h(i).all);
    end loop;
    return res;
  end Embed;

-- BUILD EVALUATOR and DIFFERENTIATOR for HOMOTOPY :

  function Evaluator
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys is

  -- DESCRIPTION :
  --   Only creates new evaluators for the last k polynomials,
  --   the rest is copied from the global data qdsys_eval.

    res : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(p'range);
    use QuadDobl_Complex_Laur_Functions;

  begin
    for i in p'first..p'last-k loop
      res(i) := qdsys_eval(i);
    end loop;
    for i in p'last-k+1..p'last loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Evaluator;

  function Jacobi_Evaluator
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat is

  -- DESCRIPTION :
  --   Only creates new differentiators for the last k polynomials,
  --   the rest is copied from the global variable qdjac_eval.

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Functions;

    res : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat(p'range,p'range);
    jac : Poly;

  begin
    for i in p'first..p'last-k loop
      for j in p'range loop
        res(i,j) := qdjac_eval(i,j);
      end loop;
    end loop;
    for i in p'last-k+1..p'last loop
      for j in p'range loop
        jac := Diff(p(i),j);
        res(i,j) := Create(jac);
        Clear(jac);
      end loop;
    end loop;
    return res;
  end Jacobi_Evaluator;

  function Equal ( x,y : QuadDobl_Complex_Vectors.Vector;
                   tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the two vectors are the same within
  --   the given tolerance in tol.

    difference : quad_double;

  begin
    for i in x'range loop
      difference := QuadDobl_Complex_Numbers.AbsVal(x(i) - y(i));
      if difference > tol
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

-- HOMOTOPY CONTINUATION :

  procedure Silent_Homotopy_Continuation
              ( target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in QuadDobl_Complex_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output,
  --   with given gamma constants for the homotopy.

  procedure Silent_Homotopy_Continuation
              ( target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output,
  --   chosing random gamma constants in the homotopy.

    a : constant QuadDobl_Complex_Vectors.Vector(target_eval'range)
      := Random_Vector(1,target_eval'last);

  begin
    Silent_Homotopy_Continuation(target_eval,target_jaco,dim,a,sols);
  end Silent_Homotopy_Continuation;

  procedure Silent_Homotopy_Continuation
              ( target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in QuadDobl_Complex_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output,
  --   with given gamma constants in the homotopy.

    one : constant quad_double := Quad_Double_Numbers.Create(1.0);
    tt : constant Complex_Number := Create(one);
    ddzero : constant quad_double := Quad_Double_Numbers.Create(0.0);

    function Hom_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return QuadDobl_Complex_Vectors.Vector is

      res : QuadDobl_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        res(i) := gamma(i)*Eval(qdsys_eval(i),x);
      end loop;
      for i in qdsys_eval'last-dim+1..qdsys_eval'last loop
        res(i) := gamma(i)*res(i)*mint + t*Eval(target_eval(i),x);
      end loop;
      return res;
    end Hom_Eval;

    function Dummy ( x : QuadDobl_Complex_Vectors.Vector; t : Complex_Number )
                   return QuadDobl_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant QuadDobl_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        for j in qdsys_eval'range loop
          res(i,j) := gamma(i)*Eval(qdjac_eval(i,j),x);
        end loop;
      end loop;
      for i in qdsys_eval'last-dim+1..qdjac_eval'last loop
        for j in qdsys_eval'range loop
          res(i,j) := gamma(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
        end loop;
      end loop;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Silent_Continue(Max_Norm,Hom_Eval,Dummy,Diff_Eval);

  begin
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(sols,Create(ddzero));
    Cont(sols,target=>tt);
  end Silent_Homotopy_Continuation;

  procedure Silent_Homotopy_Continuation_with_Stop
              ( target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                x : in QuadDobl_Complex_Vectors.Vector; tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output.

    a : constant QuadDobl_Complex_Vectors.Vector(target_eval'range)
      := Random_Vector(1,target_eval'last);
    one : constant quad_double := Quad_Double_Numbers.Create(1.0);
    tt : constant Complex_Number := Create(one);
    ddzero : constant quad_double := Quad_Double_Numbers.Create(0.0);

    function Stop_Test ( s : QuadDobl_Complex_Solutions.Solution ) 
                       return boolean is
    begin
      return Equal(s.v,x,tol);
    end Stop_Test;

    function Hom_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return QuadDobl_Complex_Vectors.Vector is

      res : QuadDobl_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        res(i) := a(i)*Eval(qdsys_eval(i),x);
      end loop;
      for i in qdsys_eval'last-dim+1..qdsys_eval'last loop
        res(i) := a(i)*res(i)*mint + t*Eval(target_eval(i),x);
      end loop;
      return res;
    end Hom_Eval;

    function Dummy ( x : QuadDobl_Complex_Vectors.Vector; t : Complex_Number )
                   return QuadDobl_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant QuadDobl_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        for j in qdsys_eval'range loop
          res(i,j) := a(i)*Eval(qdjac_eval(i,j),x);
        end loop;
      end loop;
      for i in qdsys_eval'last-dim+1..qdjac_eval'last loop
        for j in qdsys_eval'range loop
          res(i,j) := a(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
        end loop;
      end loop;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Silent_Continue_with_Stop
            (Max_Norm,Hom_Eval,Dummy,Diff_Eval,Stop_Test);

  begin
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(sols,Create(ddzero));
    Cont(sols,target=>tt);
  end Silent_Homotopy_Continuation_with_Stop;

  procedure Reporting_Homotopy_Continuation
              ( file : in file_type;
                target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in QuadDobl_Complex_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output,
  --   with given gamma constants for the homotopy.

  procedure Reporting_Homotopy_Continuation
              ( file : in file_type;
                target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output,
  --   after generating random gamma constant for the homotopy.

    a : constant QuadDobl_Complex_Vectors.Vector(1..target_eval'last)
      := Random_Vector(1,target_eval'last);

  begin
    Reporting_Homotopy_Continuation(file,target_eval,target_jaco,dim,a,sols);
  end Reporting_Homotopy_Continuation;

  procedure Reporting_Homotopy_Continuation
              ( file : in file_type;
                target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in QuadDobl_Complex_Vectors.Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output.

    one : constant quad_double := Quad_Double_Numbers.Create(1.0);
    ddzero : constant quad_double := Quad_Double_Numbers.Create(0.0);
    tt : constant Complex_Number := Create(one);

    function Hom_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return QuadDobl_Complex_Vectors.Vector is

      res : QuadDobl_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        res(i) := gamma(i)*Eval(qdsys_eval(i),x);
      end loop;
      for i in qdsys_eval'last-dim+1..qdsys_eval'last loop
        res(i) := gamma(i)*res(i)*mint + t*Eval(target_eval(i),x);
      end loop;
      return res;
    end Hom_Eval;

    function Dummy ( x : QuadDobl_Complex_Vectors.Vector; t : Complex_Number )
                   return QuadDobl_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant QuadDobl_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        for j in qdsys_eval'range loop
          res(i,j) := gamma(i)*Eval(qdjac_eval(i,j),x);
        end loop;
      end loop;
      for i in qdsys_eval'last-dim+1..qdjac_eval'last loop
        for j in qdsys_eval'range loop
          res(i,j) := gamma(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
        end loop;
      end loop;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Reporting_Continue(Max_Norm,Hom_Eval,Dummy,Diff_Eval);

  begin
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(sols,create(ddzero));
    Cont(file,sols,target=>tt);
  end Reporting_Homotopy_Continuation;

  procedure Reporting_Homotopy_Continuation_with_Stop
              ( file : in file_type;
                target_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                x : in QuadDobl_Complex_Vectors.Vector; tol : in double_float;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output.

    a : constant QuadDobl_Complex_Vectors.Vector(1..target_eval'last)
      := Random_Vector(1,target_eval'last);
    one : constant quad_double := Quad_Double_Numbers.Create(1.0);
    ddzero : constant quad_double := Quad_Double_Numbers.Create(0.0);
    tt : constant Complex_Number := Create(one);

    function Stop_Test ( s : QuadDobl_Complex_Solutions.Solution )
                       return boolean is
    begin
      return Equal(s.v,x,tol);
    end Stop_Test;

    function Hom_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return QuadDobl_Complex_Vectors.Vector is

      res : QuadDobl_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        res(i) := a(i)*Eval(qdsys_eval(i),x);
      end loop;
      for i in qdsys_eval'last-dim+1..qdsys_eval'last loop
        res(i) := a(i)*res(i)*mint + t*Eval(target_eval(i),x);
      end loop;
      return res;
    end Hom_Eval;

    function Dummy ( x : QuadDobl_Complex_Vectors.Vector; t : Complex_Number )
                   return QuadDobl_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant QuadDobl_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : QuadDobl_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use QuadDobl_Complex_Laur_Functions;

    begin
      for i in qdsys_eval'range loop
        for j in qdsys_eval'range loop
          res(i,j) := a(i)*Eval(qdjac_eval(i,j),x);
        end loop;
      end loop;
      for i in qdsys_eval'last-dim+1..qdjac_eval'last loop
        for j in qdsys_eval'range loop
          res(i,j) := a(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
        end loop;
      end loop;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Reporting_Continue_with_Stop
            (Max_Norm,Hom_Eval,Dummy,Diff_Eval,Stop_Test);

  begin
    QuadDobl_Complex_Solutions.Set_Continuation_Parameter(sols,Create(ddzero));
    Cont(file,sols,target=>tt);
  end Reporting_Homotopy_Continuation_with_Stop;

-- CORRECTORS :

  procedure Silent_Correct
              ( p_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jac_eval : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                iter : out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Applies the silent version of Newton's method to the solutions.

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Root_Refiners;

    numit : natural32 := 0;
    tmp : Solution_List := sols;

  begin
    fail := false;
    while not Is_Null(tmp) loop
      numit := 0;
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Silent_Newton(p_eval,jac_eval,ls.all,epsxa,epsfa,numit,maxit,fail);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    iter := numit;
  end Silent_Correct;

  procedure Reporting_Correct
              ( file : in file_type;
                p_eval : in QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys;
                jac_eval : in QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                iter : out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Applies the reporting version of Newton's method to the solutions.

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Root_Refiners;

    tmp : Solution_List := sols;
    numit : natural32 := 0;

  begin
    fail := false;
    while not Is_Null(tmp) loop
      numit := 0;
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Silent_Newton(p_eval,jac_eval,ls.all,epsxa,epsfa,numit,maxit,fail);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    iter := numit;
  end Reporting_Correct;

-- SAMPLING DIAGNOSTICS :

  function Satisfies ( sol : in QuadDobl_Complex_Solutions.Solution ) 
                     return boolean is

  -- DESCRIPTION :
  --   Returns true if the solution satisfies the quality test with
  --   either error or residual smaller than the tolerance and the
  --   rco less than tolerance for singular roots.

  begin
    if sol.rco < tolsing then
      return false;
    elsif sol.err < epsxa then
      return true;
    elsif sol.res < epsfa then
      return true;
    else
      return false;
    end if;
  end Satisfies;

  function Satisfies ( file : in file_type;
                       sol : in QuadDobl_Complex_Solutions.Solution )
                     return boolean is

  -- DESCRIPTION :
  --   Writes the diagnostics of the solution on file before the test.

    res : boolean;

  begin
    put(file,"  err : "); put(file,sol.err,3);
    put(file,"  rco : "); put(file,sol.rco,3);
    put(file,"  res : "); put(file,sol.res,3);
    res := Satisfies(sol);
    if res
     then put_line(file,"  success");
     else put_line(file,"  failure");
    end if;
    return res;
  end Satisfies;

-- THE TARGET OPERATIONS :

  procedure Initialize ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;

    jac : Jaco_Mat(ep'range,ep'range);
    jac_eval : Eval_Jaco_Mat(ep'range,ep'range);
    epcopy : QuadDobl_Complex_Laur_Systems.Laur_Sys(ep'range);

  begin
    QuadDobl_Complex_Laur_Systems.Copy(ep,epcopy);
    jac := Create(epcopy);
    jac_eval := Create(jac);
    qdsys := new QuadDobl_Complex_Laur_Systems.Laur_Sys'(epcopy);
    qdsys_eval := new Eval_Laur_Sys'(Create(epcopy));
    qdjac_eval := new Eval_Jaco_Mat'(jac_eval);
    Clear(jac);
  end Initialize;

  function Embedded_System return QuadDobl_Complex_Laur_Systems.Laur_Sys is
  begin
    return qdsys.all;
  end Embedded_System;

  procedure Default_Tune_Sampler ( level : in natural32 ) is
  begin
    Continuation_Parameters.Tune(level); -- ,64); -- too restrictive !
  end Default_Tune_Sampler;

  procedure Default_Tune_Sampler
              ( file : in file_type; level : in natural32 ) is
  begin
    Continuation_Parameters.Tune(level); -- ,64);
    Continuation_Parameters_io.put;
  end Default_Tune_Sampler;

  procedure Interactive_Tune_Sampler is
  begin
    Driver_for_Continuation_Parameters;
  end Interactive_Tune_Sampler;

  procedure Interactive_Tune_Sampler ( file : in file_type ) is
  begin
    Driver_for_Continuation_Parameters(file);
  end Interactive_Tune_Sampler;

  procedure Write_Refiner_Settings ( file : in file_type ) is
  begin
    new_line(file);
    put_line(file,"Current settings for the root refiner : ");
    put(file,"  1. tolerance on error on root   : ");
    put(file,epsxa,3); new_line(file);
    put(file,"  2. tolerance on residual        : ");
    put(file,epsfa,3); new_line(file);
    put(file,"  3. tolerance on singular rcond  : ");
    put(file,tolsing,3); new_line(file);
    put(file,"  4. maximal number of iterations : ");
    put(file,maxit,1); new_line(file);
  end Write_Refiner_Settings;

  procedure Default_Tune_Refiner is
  begin
    epsxa := 1.0E-48;
    epsfa := 1.0E-48;
    tolsing := 1.0E-32;
    maxit := 4;
  end Default_Tune_Refiner;

  procedure Default_Tune_Refiner ( file : in file_type ) is
  begin
    Default_Tune_Refiner;
    Write_Refiner_Settings(file);
  end Default_Tune_Refiner;

  procedure Interactive_Tune_Refiner is
 
    ans : character;

  begin
    Default_Tune_Refiner;
    loop
      Write_Refiner_Settings(Standard_Output);
      put("Type 0 to exit, 1,2,3 or 4 to change : ");
      Ask_Alternative(ans,"01234");
      exit when (ans = '0');
      case ans is
        when '1' => put("Give new tolerance on error on root : ");
                    Read_Double_Float(epsxa);
        when '2' => put("Give new tolerance on residual : ");
                    Read_Double_Float(epsfa);
        when '3' => put("Give new tolerance on singular rcond : ");
                    Read_Double_Float(tolsing);
        when '4' => put("Give new maximal number of iterations : ");
                    Read_Natural(maxit);
        when others => null;
      end case;
    end loop;
  end Interactive_Tune_Refiner;

  procedure Interactive_Tune_Refiner ( file : in file_type ) is
  begin
    Interactive_Tune_Refiner;
    Write_Refiner_Settings(file);
  end Interactive_Tune_Refiner;

-- MODIFIER :

  procedure Change_Slices ( hyp : QuadDobl_Complex_VecVecs.VecVec ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Functions;

    jac : Poly;

  begin
    for i in hyp'range loop
      Clear(qdsys(qdsys'last-hyp'last+i));
      qdsys(qdsys'last-hyp'last+i) := Hyperplane(hyp(i).all);
      Clear(qdsys_eval(qdsys_eval'last-hyp'last+i));
      qdsys_eval(qdsys_eval'last-hyp'last+i)
        := Create(qdsys(qdsys'last-hyp'last+i));
      for j in qdjac_eval'range(2) loop
        jac := Diff(qdsys(qdsys'last-hyp'last+i),j);
        Clear(qdjac_eval(qdjac_eval'last(1)-hyp'last+i,j));
        qdjac_eval(qdjac_eval'last(1)-hyp'last+i,j) := Create(jac);
        Clear(jac);
      end loop;
    end loop;
  end Change_Slices;

-- SAMPLERS :

  procedure Sample ( startsol : in QuadDobl_Complex_Solutions.Solution;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     newsol : out QuadDobl_Complex_Solutions.Solution ) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter,cnt_reruns : natural32;
    fail : boolean;

  begin
    cnt_reruns := 0;
    loop
      QuadDobl_Complex_Solutions.Add(sols,startsol);
      Silent_Homotopy_Continuation(target_eval,target_jaco,dim,sols);
      Silent_Correct(target_eval,target_jaco,sols,iter,fail);
     -- exit when (QuadDobl_Complex_Solutions.Head_Of(sols).rco > tolsing);
      exit when Satisfies(QuadDobl_Complex_Solutions.Head_Of(sols).all);
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
      QuadDobl_Complex_Solutions.Shallow_Clear(sols);
    end loop;
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsol := QuadDobl_Complex_Solutions.Head_Of(sols).all;
    QuadDobl_Complex_Solutions.Shallow_Clear(sols);
  end Sample;

  procedure Sample ( file : in file_type; full_output : in boolean;
                     startsol : in QuadDobl_Complex_Solutions.Solution;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     newsol : out QuadDobl_Complex_Solutions.Solution ) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter,cnt_reruns : natural32;
    fail : boolean;

  begin
    cnt_reruns := 0;
    loop
      QuadDobl_Complex_Solutions.Add(sols,startsol);
      if full_output then
        Reporting_Homotopy_Continuation(file,target_eval,target_jaco,dim,sols);
        Reporting_Correct(file,target_eval,target_jaco,sols,iter,fail);
      else
        Silent_Homotopy_Continuation(target_eval,target_jaco,dim,sols);
        Silent_Correct(target_eval,target_jaco,sols,iter,fail);
      end if;
     -- exit when (QuadDobl_Complex_Solutions.Head_Of(sols).rco > tolsing);
      exit when Satisfies(file,QuadDobl_Complex_Solutions.Head_Of(sols).all);
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
      QuadDobl_Complex_Solutions.Shallow_Clear(sols);
    end loop;
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsol := QuadDobl_Complex_Solutions.Head_Of(sols).all;
    QuadDobl_Complex_Solutions.Shallow_Clear(sols);
  end Sample;

  procedure Sample ( startsols : in QuadDobl_Complex_Solutions.Solution_List;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     newsols : out QuadDobl_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    QuadDobl_Complex_Solutions.Copy(startsols,sols);
    Silent_Homotopy_Continuation(target_eval,target_jaco,dim,sols);
    Silent_Correct(target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample ( startsols : in QuadDobl_Complex_Solutions.Solution_List;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     gamma : in QuadDobl_Complex_Vectors.Vector;
                     newsols : out QuadDobl_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    QuadDobl_Complex_Solutions.Copy(startsols,sols);
    Silent_Homotopy_Continuation(target_eval,target_jaco,dim,gamma,sols);
    Silent_Correct(target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample ( file : in file_type;
                     startsols : in QuadDobl_Complex_Solutions.Solution_List;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     newsols : out QuadDobl_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    QuadDobl_Complex_Solutions.Copy(startsols,sols);
    Reporting_Homotopy_Continuation(file,target_eval,target_jaco,dim,sols);
    Reporting_Correct(file,target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample ( file : in file_type;
                     startsols : in QuadDobl_Complex_Solutions.Solution_List;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     gamma : in QuadDobl_Complex_Vectors.Vector;
                     newsols : out QuadDobl_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    QuadDobl_Complex_Solutions.Copy(startsols,sols);
    Reporting_Homotopy_Continuation
      (file,target_eval,target_jaco,dim,gamma,sols);
    Reporting_Correct(file,target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample_with_Stop
                   ( startsols : in QuadDobl_Complex_Solutions.Solution_List;
                     x : in QuadDobl_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     newsols : out QuadDobl_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    QuadDobl_Complex_Solutions.Copy(startsols,sols);
    Silent_Homotopy_Continuation_with_Stop
      (target_eval,target_jaco,dim,x,tol,sols);
    Silent_Correct(target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample_with_Stop;

  procedure Sample_with_Stop
                   ( file : in file_type;
                     startsols : in QuadDobl_Complex_Solutions.Solution_List;
                     x : in QuadDobl_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in QuadDobl_Complex_VecVecs.VecVec;
                     newsols : out QuadDobl_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(QuadDobl_System,newhyp);
    target_eval : QuadDobl_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Evaluator(target,dim);
    target_jaco : QuadDobl_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Jacobi_Evaluator(target,dim);
    sols : QuadDobl_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    QuadDobl_Complex_Solutions.Copy(startsols,sols);
    Reporting_Homotopy_Continuation_with_Stop
      (file,target_eval,target_jaco,dim,x,tol,sols);
    Reporting_Correct(file,target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      QuadDobl_Complex_Laurentials.Clear(target(i));
      QuadDobl_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        QuadDobl_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample_with_Stop;

-- DEALLOCATION :

  procedure Clear is
  begin
   -- QuadDobl_Complex_Laur_Systems.Shallow_Clear(qdsys); --  data sharing
    QuadDobl_Complex_Laur_Systems.Clear(qdsys); -- no data sharing!
    QuadDobl_Complex_Laur_SysFun.Clear(qdsys_eval);
    QuadDobl_Complex_Laur_JacoMats.Clear(qdjac_eval);
  end Clear;

end QuadDobl_Sampling_Laurent_Machine;
