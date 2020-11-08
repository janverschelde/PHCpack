with Numbers_io;                         use Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Standard_Integer_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Standard_Complex_Laurentials;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_JacoMats;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laur_Functions;
with Multprec_Complex_Laur_SysFun;
with Multprec_Complex_Laur_JacoMats;
with Continuation_Parameters;
with Continuation_Parameters_io;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;
with Main_Poly_Continuation;             use Main_Poly_Continuation;
with Standard_Root_Refiners;
with Multprec_Root_Refiners;
with Planes_and_Polynomials;             use Planes_and_Polynomials;

package body Sampling_Laurent_Machine is

-- INTERNAL DATA :

  stansys : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
  stansys_eval : Standard_Complex_Laur_SysFun.Link_to_Eval_Laur_Sys;
  stanjac_eval : Standard_Complex_Laur_JacoMats.Link_to_Eval_Jaco_Mat;
  orgsys,multsys : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;
  multsys_eval : Multprec_Complex_Laur_SysFun.Link_to_Eval_Laur_Sys;
  multjac_eval : Multprec_Complex_Laur_JacoMats.Link_to_Eval_Jaco_Mat;

  restricted : boolean;

  rststansys : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
  rststansys_eval : Standard_Complex_Laur_SysFun.Link_to_Eval_Laur_Sys;
  rststanjac_eval : Standard_Complex_Laur_JacoMats.Link_to_Eval_Jaco_Mat;
  rstmultsys : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;
  rstmultsys_eval : Multprec_Complex_Laur_SysFun.Link_to_Eval_Laur_Sys;
  rstmultjac_eval : Multprec_Complex_Laur_JacoMats.Link_to_Eval_Jaco_Mat;

-- PARAMETERS FOR ROOT REFINER :

  epsxa,epsfa,tolsing : double_float;
  maxit : natural32;
  mpepsxa,mpepsfa,mptolsing : Floating_Number;
  mpsize,mpmaxit : natural32;

-- AUXILIARY ROUTINES :

  function Number_of_Equations return integer32 is

  -- DESCRIPTION :
  --   Returns the number of equations in the embedded system.

  begin
    if restricted
     then return rststansys'last;
     else return stansys'last;
    end if;
  end Number_of_Equations;

  function Standard_System return Standard_Complex_Laur_Systems.Laur_Sys is

  -- DESCRIPTION :
  --   Returns the embedded polynomial system, wrt restricted.

  begin
    if restricted
     then return rststansys.all;
     else return stansys.all;
    end if;
  end Standard_System;

  function Embed ( t : Multprec_Complex_Laurentials.Term; k : integer32 )
                 return Multprec_Complex_Laurentials.Term is

  -- DESCRIPTION :
  --   Extends the variable range of t with k variables.

    res : Multprec_Complex_Laurentials.Term;

  begin
    Multprec_Complex_Numbers.Copy(t.cf,res.cf);
    res.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last+k);
    res.dg(t.dg'range) := t.dg.all;
    for i in t.dg'last+1..res.dg'last loop
      res.dg(i) := 0;
    end loop;
    return res;
  end Embed;

  procedure Embedded_Adds
              ( mp : in out Multprec_Complex_Laurentials.Poly;
                sp : in Standard_Complex_Laurentials.Poly;
                n,k : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Adds c_i*z_i, for i in n+1..n+k, where c_i is a random coefficient.

    t : Multprec_Complex_Laurentials.Term;
    spcf : Standard_Complex_Numbers.Complex_Number;

  begin
    t.dg := new Standard_Integer_Vectors.Vector'(1..n+k => 0);
    for i in 1..k loop
      t.dg(n+i) := 1;
      spcf := Standard_Complex_Laurentials.Coeff
                (sp,Standard_Complex_Laurentials.Degrees(t.dg));
      t.cf := Multprec_Complex_Number_Tools.Create(spcf);
      Set_Size(t.cf,size);
      Multprec_Complex_Laurentials.Add(mp,t);
      t.dg(n+i) := 0;
      Multprec_Complex_Numbers.Clear(t.cf);
    end loop;
    Multprec_Complex_Laurentials.Clear(t.dg);
  end Embedded_Adds;

  function Embed ( mp : Multprec_Complex_Laurentials.Poly;
                   sp : Standard_Complex_Laurentials.Poly;
                   n,k : integer32; size : natural32 )
                 return Multprec_Complex_Laurentials.Poly is

  -- DESCRIPTION :
  --   Adds k new variables to every equation.

    res : Multprec_Complex_Laurentials.Poly;

    procedure Embed_Term ( t : Multprec_Complex_Laurentials.Term;
                           cont : out boolean ) is

      et : Multprec_Complex_Laurentials.Term := Embed(t,k);

    begin
      Multprec_Complex_Laurentials.Add(res,et);
      Multprec_Complex_Laurentials.Clear(et);
      cont := true;
    end Embed_Term;
    procedure Embed_Terms is
      new Multprec_Complex_Laurentials.Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(mp);
    Embedded_Adds(res,sp,n,k,size);
    return res;
  end Embed;

  function Embed ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                   h : in Standard_Complex_VecVecs.VecVec ) 
                 return Standard_Complex_Laur_Systems.Laur_Sys is

  -- DESCRIPTION :
  --   Returns the embedding of p with other slices in h.

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);
    k : constant integer32 := h'last;

  begin
    res(p'first..p'last-k) := p(p'first..p'last-k);
    for i in 1..k loop
      res(p'last-k+i) := Hyperplane(h(i).all);
    end loop;
    return res;
  end Embed;

-- BUILD EVALUATOR and DIFFERENTIATOR for HOMOTOPY :

  function Standard_Evaluator
             ( p : Standard_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return Standard_Complex_Laur_SysFun.Eval_Laur_Sys is

  -- DESCRIPTION :
  --   Only creates new evaluators for the last k polynomials,
  --   the rest is copied from the global data stansys_eval.

    res : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(p'range);
    use Standard_Complex_Laur_Functions;

  begin
    if restricted then
      for i in p'first..p'last-k loop
        res(i) := rststansys_eval(i);
      end loop;
    else
      for i in p'first..p'last-k loop
        res(i) := stansys_eval(i);
      end loop;
    end if;
    for i in p'last-k+1..p'last loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Standard_Evaluator;

  function Multprec_Evaluator
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return Multprec_Complex_Laur_SysFun.Eval_Laur_Sys is

  -- DESCRIPTION :
  --   This is the multi-precision analogue of the function above.

    res : Multprec_Complex_Laur_SysFun.Eval_Laur_Sys(p'range);
    use Multprec_Complex_Laur_Functions;

  begin
    if restricted then
      for i in p'first..p'last-k loop
        res(i) := rstmultsys_eval(i);
      end loop;
    else
      for i in p'first..p'last-k loop
        res(i) := multsys_eval(i);
      end loop;
    end if;
    for i in p'last-k+1..p'last loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Multprec_Evaluator;

  function Standard_Jacobi_Evaluator
             ( p : Standard_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat is

  -- DESCRIPTION :
  --   Only creates new differentiators for the last k polynomials,
  --   the rest is copied from the global variable stanjac_eval.

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Functions;

    res : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat(p'range,p'range);

  begin
    if restricted then
      for i in p'first..p'last-k loop
        for j in p'range loop
          res(i,j) := rststanjac_eval(i,j);
        end loop;
      end loop;
    else
      for i in p'first..p'last-k loop
        for j in p'range loop
          res(i,j) := stanjac_eval(i,j);
        end loop;
      end loop;
    end if;
    for i in p'last-k+1..p'last loop
      for j in p'range loop
        declare
          jac : Poly;
        begin
          jac := Diff(p(i),j);
          res(i,j) := Create(jac);
          Clear(jac);
        end;
      end loop;
    end loop;
    return res;
  end Standard_Jacobi_Evaluator;

  function Multprec_Jacobi_Evaluator
             ( p : Multprec_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat is

    use Multprec_Complex_Laurentials;
    use Multprec_Complex_Laur_Functions;

    res : Multprec_Complex_Laur_JacoMats.Eval_Jaco_Mat(p'range,p'range);
    jac : Poly;

  begin
    if restricted then
      for i in p'first..p'last-k loop
        for j in p'range loop
          res(i,j) := rstmultjac_eval(i,j);
        end loop;
      end loop;
    else
      for i in p'first..p'last-k loop
        for j in p'range loop
          res(i,j) := multjac_eval(i,j);
        end loop;
      end loop;
    end if;
    for i in p'last-k+1..p'last loop
      for j in p'range loop
        jac := Diff(p(i),j);
        res(i,j) := Create(jac);
        Clear(jac);
      end loop;
    end loop;
    return res;
  end Multprec_Jacobi_Evaluator;

-- HOMOTOPY CONTINUATION :

  procedure Silent_Homotopy_Continuation
              ( target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in Standard_Complex_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output,
  --   with given gamma constants for the homotopy.

  procedure Silent_Homotopy_Continuation
              ( target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output,
  --   chosing random gamma constants in the homotopy.

    a : constant Standard_Complex_Vectors.Vector(target_eval'range)
      := Random_Vector(1,target_eval'last);

  begin
    Silent_Homotopy_Continuation(target_eval,target_jaco,dim,a,sols);
  end Silent_Homotopy_Continuation;

  procedure Silent_Homotopy_Continuation
              ( target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in Standard_Complex_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output,
  --   with given gamma constants in the homotopy.

    tt : constant Complex_Number := Create(1.0);

    function Hom_Eval ( x : Standard_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          res(i) := gamma(i)*Eval(rststansys_eval(i),x);
        end loop;
        for i in rststansys_eval'last-dim+1..rststansys_eval'last loop
          res(i) := gamma(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      else
        for i in stansys_eval'range loop
          res(i) := gamma(i)*Eval(stansys_eval(i),x);
        end loop;
        for i in stansys_eval'last-dim+1..stansys_eval'last loop
          res(i) := gamma(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      end if;
      return res;
    end Hom_Eval;

    function Dummy ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                   return Standard_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant Standard_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : Standard_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          for j in rststansys_eval'range loop
            res(i,j) := gamma(i)*Eval(rststanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in rststansys_eval'last-dim+1..rststanjac_eval'last loop
          for j in rststansys_eval'range loop
            res(i,j) := gamma(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      else
        for i in stansys_eval'range loop
          for j in stansys_eval'range loop
            res(i,j) := gamma(i)*Eval(stanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in stansys_eval'last-dim+1..stanjac_eval'last loop
          for j in stansys_eval'range loop
            res(i,j) := gamma(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      end if;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Silent_Continue(Max_Norm,Hom_Eval,Dummy,Diff_Eval);

  begin
    Standard_Complex_Solutions.Set_Continuation_Parameter(sols,Create(0.0));
    Cont(sols,false,target=>tt);
  end Silent_Homotopy_Continuation;

  procedure Silent_Homotopy_Continuation_with_Stop
              ( target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                x : in Standard_Complex_Vectors.Vector; tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation without intermediate output.

    a : constant Standard_Complex_Vectors.Vector(target_eval'range)
      := Random_Vector(1,target_eval'last);
    tt : constant Complex_Number := Create(1.0);

    function Stop_Test ( s : Standard_Complex_Solutions.Solution ) 
                       return boolean is
    begin
      return Equal(s.v,x,tol);
    end Stop_Test;

    function Hom_Eval ( x : Standard_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          res(i) := a(i)*Eval(rststansys_eval(i),x);
        end loop;
        for i in rststansys_eval'last-dim+1..rststansys_eval'last loop
          res(i) := a(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      else
        for i in stansys_eval'range loop
          res(i) := a(i)*Eval(stansys_eval(i),x);
        end loop;
        for i in stansys_eval'last-dim+1..stansys_eval'last loop
          res(i) := a(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      end if;
      return res;
    end Hom_Eval;

    function Dummy ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                   return Standard_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant Standard_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : Standard_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          for j in rststansys_eval'range loop
            res(i,j) := a(i)*Eval(rststanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in rststansys_eval'last-dim+1..rststanjac_eval'last loop
          for j in rststansys_eval'range loop
            res(i,j) := a(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      else
        for i in stansys_eval'range loop
          for j in stansys_eval'range loop
            res(i,j) := a(i)*Eval(stanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in stansys_eval'last-dim+1..stanjac_eval'last loop
          for j in stansys_eval'range loop
            res(i,j) := a(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      end if;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Silent_Continue_with_Stop
            (Max_Norm,Hom_Eval,Dummy,Diff_Eval,Stop_Test);

  begin
    Standard_Complex_Solutions.Set_Continuation_Parameter(sols,Create(0.0));
    Cont(sols,false,target=>tt);
  end Silent_Homotopy_Continuation_with_Stop;

  procedure Reporting_Homotopy_Continuation
              ( file : in file_type;
                target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in Standard_Complex_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output,
  --   with given gamma constants for the homotopy.

  procedure Reporting_Homotopy_Continuation
              ( file : in file_type;
                target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output,
  --   after generating random gamma constant for the homotopy.

    a : constant Standard_Complex_Vectors.Vector(1..target_eval'last)
      := Random_Vector(1,target_eval'last);

  begin
    Reporting_Homotopy_Continuation(file,target_eval,target_jaco,dim,a,sols);
  end Reporting_Homotopy_Continuation;

  procedure Reporting_Homotopy_Continuation
              ( file : in file_type;
                target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                gamma : in Standard_Complex_Vectors.Vector;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output,
  --   generating 

    tt : constant Complex_Number := Create(1.0);

    function Hom_Eval ( x : Standard_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          res(i) := gamma(i)*Eval(rststansys_eval(i),x);
        end loop;
        for i in rststansys_eval'last-dim+1..rststansys_eval'last loop
          res(i) := gamma(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      else
        for i in stansys_eval'range loop
          res(i) := gamma(i)*Eval(stansys_eval(i),x);
        end loop;
        for i in stansys_eval'last-dim+1..stansys_eval'last loop
          res(i) := gamma(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      end if;
      return res;
    end Hom_Eval;

    function Dummy ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                   return Standard_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant Standard_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : Standard_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          for j in rststansys_eval'range loop
            res(i,j) := gamma(i)*Eval(rststanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in rststansys_eval'last-dim+1..rststanjac_eval'last loop
          for j in rststansys_eval'range loop
            res(i,j) := gamma(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      else
        for i in stansys_eval'range loop
          for j in stansys_eval'range loop
            res(i,j) := gamma(i)*Eval(stanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in stansys_eval'last-dim+1..stanjac_eval'last loop
          for j in stansys_eval'range loop
            res(i,j) := gamma(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      end if;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Reporting_Continue(Max_Norm,Hom_Eval,Dummy,Diff_Eval);

  begin
    Standard_Complex_Solutions.Set_Continuation_Parameter(sols,Create(0.0));
    Cont(file,sols,false,target=>tt);
  end Reporting_Homotopy_Continuation;

  procedure Reporting_Homotopy_Continuation_with_Stop
              ( file : in file_type;
                target_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                target_jaco : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                dim : in integer32;
                x : in Standard_Complex_Vectors.Vector; tol : in double_float;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Performs homotopy continuation with intermediate output.

    a : constant Standard_Complex_Vectors.Vector(1..target_eval'last)
      := Random_Vector(1,target_eval'last);
    tt : constant Complex_Number := Create(1.0);

    function Stop_Test ( s : Standard_Complex_Solutions.Solution )
                       return boolean is
    begin
      return Equal(s.v,x,tol);
    end Stop_Test;

    function Hom_Eval ( x : Standard_Complex_Vectors.Vector;
                        t : Complex_Number )
                      return Standard_Complex_Vectors.Vector is

      res : Standard_Complex_Vectors.Vector(x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          res(i) := a(i)*Eval(rststansys_eval(i),x);
        end loop;
        for i in rststansys_eval'last-dim+1..rststansys_eval'last loop
          res(i) := a(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      else
        for i in stansys_eval'range loop
          res(i) := a(i)*Eval(stansys_eval(i),x);
        end loop;
        for i in stansys_eval'last-dim+1..stansys_eval'last loop
          res(i) := a(i)*res(i)*mint + t*Eval(target_eval(i),x);
        end loop;
      end if;
      return res;
    end Hom_Eval;

    function Dummy ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
                   return Standard_Complex_Vectors.Vector is

    -- NOTE : this procedure will not be used with Tune(0) and Tune(2).

      res : constant Standard_Complex_Vectors.Vector(x'range) := x;

    begin
      return res;
    end Dummy;

    function Diff_Eval ( x : Standard_Complex_Vectors.Vector;
                         t : Complex_Number ) return Matrix is

      res : Matrix(x'range,x'range);
      mint : constant Complex_Number := Create(1.0) - t;

      use Standard_Complex_Laur_Functions;

    begin
      if restricted then
        for i in rststansys_eval'range loop
          for j in rststansys_eval'range loop
            res(i,j) := a(i)*Eval(rststanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in rststansys_eval'last-dim+1..rststanjac_eval'last loop
          for j in rststansys_eval'range loop
            res(i,j) := a(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      else
        for i in stansys_eval'range loop
          for j in stansys_eval'range loop
            res(i,j) := a(i)*Eval(stanjac_eval(i,j),x);
          end loop;
        end loop;
        for i in stansys_eval'last-dim+1..stanjac_eval'last loop
          for j in stansys_eval'range loop
            res(i,j) := a(i)*res(i,j)*mint + t*Eval(target_jaco(i,j),x);
          end loop;
        end loop;
      end if;
      return res;
    end Diff_Eval;

    procedure Cont is
      new Reporting_Continue_with_Stop
            (Max_Norm,Hom_Eval,Dummy,Diff_Eval,Stop_Test);

  begin
    Standard_Complex_Solutions.Set_Continuation_Parameter(sols,Create(0.0));
    Cont(file,sols,false,target=>tt);
  end Reporting_Homotopy_Continuation_with_Stop;

-- CORRECTORS :

  procedure Standard_Silent_Correct
              ( p_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                jac_eval : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                iter : out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Applies the silent version of Newton's method to the solutions.

    use Standard_Complex_Solutions;
    use Standard_Root_Refiners;

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
  end Standard_Silent_Correct;

  procedure Standard_Reporting_Correct
              ( file : in file_type;
                p_eval : in Standard_Complex_Laur_SysFun.Eval_Laur_Sys;
                jac_eval : in Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat;
                sols : in out Standard_Complex_Solutions.Solution_List;
                iter : out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Applies the reporting version of Newton's method to the solutions.

    use Standard_Complex_Solutions;
    use Standard_Root_Refiners;

    tmp : Solution_List := sols;
    numit : natural32 := 0;

  begin
    fail := false;
    while not Is_Null(tmp) loop
      numit := 0;
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        Reporting_Newton
          (file,p_eval,jac_eval,ls.all,epsxa,epsfa,numit,maxit,fail);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    iter := numit;
  end Standard_Reporting_Correct;

  procedure Multprec_Silent_Correct
              ( sol : in out Multprec_Complex_Solutions.Solution;
                hyp : in Multprec_Complex_VecVecs.VecVec;
                iter : out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Applies the silent version of Newton's method to the solution.

    use Multprec_Complex_Laur_SysFun;
    use Multprec_Complex_Laur_JacoMats;
    use Multprec_Complex_Solutions;
    use Multprec_Root_Refiners;

    dim : constant integer32 := hyp'last;
    numit : natural32 := 0;
    neq : constant integer32 := Number_of_Equations;
    mps : Multprec_Complex_Laur_Systems.Laur_Sys(1..neq);
    p_eval : Eval_Laur_Sys(1..neq);
    jac_eval : Eval_Jaco_Mat(1..neq,1..neq);

  begin
    fail := false;
    for i in hyp'range loop
      mps(mps'last-dim+i) := Hyperplane(hyp(i).all);
    end loop;
    p_eval := Multprec_Evaluator(mps,dim);
    jac_eval := Multprec_Jacobi_Evaluator(mps,dim);
    Silent_Newton(p_eval,jac_eval,sol,mpepsxa,mpepsfa,numit,mpmaxit,fail);
    for i in neq-dim+1..neq loop
      Multprec_Complex_Laurentials.Clear(mps(i));
      Multprec_Complex_Laur_Functions.Clear(p_eval(i));
      for j in 1..neq loop
        Multprec_Complex_Laur_Functions.Clear(jac_eval(i,j));
      end loop;
    end loop;
    iter := numit;
  end Multprec_Silent_Correct;

  procedure Multprec_Reporting_Correct
              ( file : in file_type;
                sol : in out Multprec_Complex_Solutions.Solution;
                hyp : in Multprec_Complex_VecVecs.VecVec;
                iter : out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Applies the reporting version of Newton's method to the solution.

    use Multprec_Complex_Laur_SysFun;
    use Multprec_Complex_Laur_JacoMats;
    use Multprec_Complex_Solutions;
    use Multprec_Root_Refiners;

    dim : constant integer32 := hyp'last;
    numit : natural32 := 0;
    neq : constant integer32 := Number_of_Equations;
    mps : Multprec_Complex_Laur_Systems.Laur_Sys(1..neq);
    p_eval : Eval_Laur_Sys(1..neq);
    jac_eval : Eval_Jaco_Mat(1..neq,1..neq);

  begin
    fail := false;
    for i in hyp'range loop
      mps(mps'last-dim+i) := Hyperplane(hyp(i).all);
    end loop;
    p_eval := Multprec_Evaluator(mps,dim);
    jac_eval := Multprec_Jacobi_Evaluator(mps,dim);
    Reporting_Newton
      (file,p_eval,jac_eval,sol,mpepsxa,mpepsfa,numit,mpmaxit,fail);
    for i in neq-dim+1..neq loop
      Multprec_Complex_Laurentials.Clear(mps(i));
      Multprec_Complex_Laur_Functions.Clear(p_eval(i));
      for j in 1..neq loop
        Multprec_Complex_Laur_Functions.Clear(jac_eval(i,j));
      end loop;
    end loop;
    iter := numit;
  end Multprec_Reporting_Correct;

-- SAMPLING DIAGNOSTICS :

  function Satisfies ( sol : in Standard_Complex_Solutions.Solution ) 
                     return boolean is

  -- DESCRIPTION :
  --   Returns true if the solution satisfies the quality test with
  --   either error or residual smaller than the tolerance and the
  --   rco less than tolerance for singular roots.

  begin
    if sol.rco < tolsing then
      return false;
    elsif sol.err < epsxa/sol.rco then
      return true;
    elsif sol.res < epsfa/sol.rco then
      return true;
    else
      return false;
    end if;
  end Satisfies;

  function Satisfies ( file : in file_type;
                       sol : in Standard_Complex_Solutions.Solution )
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

  procedure Write_Summary ( file : in file_type; iter : in natural32;
                            sol : in Multprec_Complex_Solutions.Solution ) is

  -- DESCIPTION :
  --   Writes the (err,rco,res) of the solution on file.

  begin
    put(file,"  err : "); put(file,sol.err,3);
    put(file,"  rco : "); put(file,sol.rco,3);
    put(file,"  res : "); put(file,sol.res,3);
    put(file,"  #it : "); put(file,iter,1);  new_line(file);
  end Write_Summary;

-- THE TARGET OPERATIONS :

  procedure Initialize ( ep : in Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;

    epcopy : Standard_Complex_Laur_Systems.Laur_Sys(ep'range);
    jac : Jaco_Mat(ep'range,ep'range);
    jac_eval : Eval_Jaco_Mat(ep'range,ep'range);

  begin
    restricted := false;
    Standard_Complex_Laur_Systems.Copy(ep,epcopy);
    jac := Create(epcopy);
    jac_eval := Create(jac);
    stansys := new Standard_Complex_Laur_Systems.Laur_Sys'(epcopy);
    stansys_eval := new Eval_Laur_Sys'(Create(epcopy));
    stanjac_eval := new Eval_Jaco_Mat'(jac_eval);
    Clear(jac);
  end Initialize;

  function Embedded_System return Standard_Complex_Laur_Systems.Laur_Sys is
  begin
    return stansys.all;
  end Embedded_System;

  function Original_System return Multprec_Complex_Laur_Systems.Laur_Sys is
  begin
    return orgsys.all;
  end Original_System;

  procedure Initialize_Restricted
                 ( ep : in Standard_Complex_Laur_Systems.Laur_Sys ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;

    epcopy : Standard_Complex_Laur_Systems.Laur_Sys(ep'range);
    jac : Jaco_Mat(ep'range,ep'range);
    jac_eval : Eval_Jaco_Mat(ep'range,ep'range);

  begin
    restricted := true;
    Standard_Complex_Laur_Systems.Copy(ep,epcopy);
    jac := Create(epcopy);
    jac_eval := Create(jac);
    rststansys := new Standard_Complex_Laur_Systems.Laur_Sys'(epcopy);
    rststansys_eval := new Eval_Laur_Sys'(Create(epcopy));
    rststanjac_eval := new Eval_Jaco_Mat'(jac_eval);
    Clear(jac);
  end Initialize_Restricted;

  procedure Initialize ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                         mp : in Multprec_Complex_Laur_Systems.Laur_Sys;
                         k : in integer32; size : in natural32 ) is

    use Multprec_Complex_Laurentials;
    use Multprec_Complex_Laur_Functions;
    use Multprec_Complex_Laur_SysFun;
    use Multprec_Complex_Laur_JacoMats;

    epcopy : Standard_Complex_Laur_Systems.Laur_Sys(ep'range);
    mps : Multprec_Complex_Laur_Systems.Laur_Sys(ep'range);
    mpjac : Jaco_Mat(ep'range,ep'range);
    stjac : Standard_Complex_Laur_JacoMats.Jaco_Mat(ep'range,ep'range);
    stjac_eval
      : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat(ep'range,ep'range);

  begin
    restricted := false;
    Standard_Complex_Laur_Systems.Copy(ep,epcopy);
    stjac := Standard_Complex_Laur_JacoMats.Create(epcopy);
    stjac_eval := Standard_Complex_Laur_JacoMats.Create(stjac);
    stansys := new Standard_Complex_Laur_Systems.Laur_Sys'(epcopy);
    stansys_eval
      := new Standard_Complex_Laur_SysFun.Eval_Laur_Sys'
               (Standard_Complex_Laur_SysFun.Create(epcopy));
    stanjac_eval
      := new Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat'(stjac_eval);
    Standard_Complex_Laur_JacoMats.Clear(stjac);
    orgsys := new Multprec_Complex_Laur_Systems.Laur_Sys'(mp);
    for i in ep'first..ep'last-k loop
      mps(i) := Embed(mp(i),ep(i),mp'last,k,size);
    end loop;
    for i in ep'last-k+1..ep'last loop
      mps(i) := Convert(ep(i));
      Set_Size(mps(i),size);
    end loop;
    multsys := new Multprec_Complex_Laur_Systems.Laur_Sys'(mps);
    multsys_eval := new Eval_Laur_Sys'(Create(mps));
    mpjac := Create(mps);
    multjac_eval := new Eval_Jaco_Mat'(Create(mpjac));
    Clear(mpjac);
    mpsize := size;
  end Initialize;

  procedure Initialize_Restricted 
                   ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                     mp : in Multprec_Complex_Laur_Systems.Laur_Sys;
                     k : in integer32; size : in natural32 ) is

    use Multprec_Complex_Laurentials;
    use Multprec_Complex_Laur_Functions;
    use Multprec_Complex_Laur_SysFun;
    use Multprec_Complex_Laur_JacoMats;

    epcopy : Standard_Complex_Laur_Systems.Laur_Sys(ep'range);
    mps : Multprec_Complex_Laur_Systems.Laur_Sys(ep'range);
    mpjac : Jaco_Mat(ep'range,ep'range);
    stjac : Standard_Complex_Laur_JacoMats.Jaco_Mat(ep'range,ep'range);
    stjac_eval
      : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat(ep'range,ep'range);

  begin
    restricted := true;
    Standard_Complex_Laur_Systems.Copy(ep,epcopy);
    stjac := Standard_Complex_Laur_JacoMats.Create(ep);
    stjac_eval := Standard_Complex_Laur_JacoMats.Create(stjac);
    rststansys := new Standard_Complex_Laur_Systems.Laur_Sys'(epcopy);
    rststansys_eval
      := new Standard_Complex_Laur_SysFun.Eval_Laur_Sys'
               (Standard_Complex_Laur_SysFun.Create(epcopy));
    rststanjac_eval
      := new Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat'(stjac_eval);
    Standard_Complex_Laur_JacoMats.Clear(stjac);
    for i in ep'first..ep'last-k loop
      mps(i) := Embed(mp(i),ep(i),mp'last,k,size);
    end loop;
    for i in ep'last-k+1..ep'last loop
      mps(i) := Convert(ep(i));
      Set_Size(mps(i),size);
    end loop;
    rstmultsys := new Multprec_Complex_Laur_Systems.Laur_Sys'(mps);
    rstmultsys_eval := new Eval_Laur_Sys'(Create(mps));
    mpjac := Create(mps);
    rstmultjac_eval := new Eval_Jaco_Mat'(Create(mpjac));
    Clear(mpjac);
    mpsize := size;
  end Initialize_Restricted;

  procedure Default_Tune_Sampler ( level : in natural32 ) is
  begin
    Continuation_Parameters.Tune(level);
  end Default_Tune_Sampler;

  procedure Default_Tune_Sampler
              ( file : in file_type; level : in natural32 ) is
  begin
    Continuation_Parameters.Tune(level);
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

  procedure Write_Standard_Refiner_Settings ( file : in file_type ) is
  begin
    new_line(file);
    put_line(file,"Current settings for standard root refiner : ");
    put(file,"  1. tolerance on error on root   : ");
    put(file,epsxa,3); new_line(file);
    put(file,"  2. tolerance on residual        : ");
    put(file,epsfa,3); new_line(file);
    put(file,"  3. tolerance on singular rcond  : ");
    put(file,tolsing,3); new_line(file);
    put(file,"  4. maximal number of iterations : ");
    put(file,maxit,1); new_line(file);
  end Write_Standard_Refiner_Settings;

  procedure Write_Multprec_Refiner_Settings ( file : in file_type ) is
  begin
    new_line(file);
    put_line(file,"Current settings for multi-precision root refiner : ");
    put(file,"  1. number of decimal places     : ");
    put(file,Size_to_Decimal(mpsize),1); new_line(file);   
    put(file,"  2. tolerance on error on root   : ");
    put(file,mpepsxa,3); new_line(file);
    put(file,"  3. tolerance on residual        : ");
    put(file,mpepsfa,3); new_line(file);
    put(file,"  4. tolerance on singular rcond  : ");
    put(file,mptolsing,3); new_line(file);
    put(file,"  5. maximal number of iterations : ");
    put(file,mpmaxit,1); new_line(file);
  end Write_Multprec_Refiner_Settings;

  procedure Default_Tune_Refiner is
  begin
    epsxa := 1.0E-14;
    epsfa := 1.0E-14;
    tolsing := 1.0E-08;
    maxit := 4;
  end Default_Tune_Refiner;

  procedure Default_Tune_Refiner ( file : in file_type ) is
  begin
    Default_Tune_Refiner;
    Write_Standard_Refiner_Settings(file);
  end Default_Tune_Refiner;

  procedure Default_Tune_Refiner ( size : in natural32 ) is

    deci : constant natural32 := Size_to_Decimal(size);
    tol : constant double_float := 10.0**(-8);
    epserr : constant double_float := 10.0**integer(-deci+8);
    epsres : constant double_float := 10.0**integer(-deci+8);

  begin
    mpsize := size;
    mpmaxit := 4;
    mptolsing := Create(tol);
    mpepsxa := Create(epserr);
    mpepsfa := Create(epsres);
  end Default_Tune_Refiner;

  procedure Default_Tune_Refiner
              ( file : in file_type; size : in natural32 ) is
  begin
    Default_Tune_Refiner(size);
    Write_Multprec_Refiner_Settings(file);
  end Default_Tune_Refiner;

  procedure Interactive_Tune_Refiner is
 
    ans : character;

  begin
    Default_Tune_Refiner;
    loop
      Write_Standard_Refiner_Settings(Standard_Output);
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
    Write_Standard_Refiner_Settings(file);
  end Interactive_Tune_Refiner;

  procedure Interactive_Tune_Refiner ( size : in natural32 ) is

    deci : natural32;
    ans : character;

  begin
    Default_Tune_Refiner(size);
    loop
      Write_Multprec_Refiner_Settings(Standard_Output);
      put("Type 0 to exit, 1,2,3,4 or 5 to change : ");
      Ask_Alternative(ans,"012345");
      exit when (ans = '0');
      case ans is
        when '1' => put("Give new number of decimal places : ");
                    Read_Natural(deci);
                    mpsize := Decimal_to_Size(deci);
                    Default_Tune_Refiner(mpsize);
        when '2' => put("Give new tolerance on error on root : ");
                    get(mpepsxa);
        when '3' => put("Give new tolerance on residual : ");
                    get(mpepsfa);
        when '4' => put("Give new tolerance on singular rcond : ");
                    get(mptolsing);
        when '5' => put("Give new maximal number of iterations : ");
                    Read_Natural(mpmaxit);
        when others => null;
      end case;
    end loop;
  end Interactive_Tune_Refiner;

  procedure Interactive_Tune_Refiner
                ( file : in file_type; size : in natural32 ) is
  begin
    Interactive_Tune_Refiner(size);
    Write_Multprec_Refiner_Settings(file);
  end Interactive_Tune_Refiner;

-- MODIFIER :

  procedure Change_Slices ( hyp : Standard_Complex_VecVecs.VecVec ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Functions;

    jac : Poly;

  begin
    for i in hyp'range loop
      Clear(stansys(stansys'last-hyp'last+i));
      stansys(stansys'last-hyp'last+i) := Hyperplane(hyp(i).all);
      Clear(stansys_eval(stansys_eval'last-hyp'last+i));
      stansys_eval(stansys_eval'last-hyp'last+i)
        := Create(stansys(stansys'last-hyp'last+i));
      for j in stanjac_eval'range(2) loop
        jac := Diff(stansys(stansys'last-hyp'last+i),j);
        Clear(stanjac_eval(stanjac_eval'last(1)-hyp'last+i,j));
        stanjac_eval(stanjac_eval'last(1)-hyp'last+i,j) := Create(jac);
        Clear(jac);
      end loop;
    end loop;
  end Change_Slices;

-- SAMPLERS :

  procedure Sample ( startsol : in Standard_Complex_Solutions.Solution;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsol : out Standard_Complex_Solutions.Solution ) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter,cnt_reruns : natural32;
    fail : boolean;

  begin
    cnt_reruns := 0;
    loop
      Standard_Complex_Solutions.Add(sols,startsol);
      Silent_Homotopy_Continuation(target_eval,target_jaco,dim,sols);
      Standard_Silent_Correct(target_eval,target_jaco,sols,iter,fail);
     -- exit when (Standard_Complex_Solutions.Head_Of(sols).rco > tolsing);
      exit when Satisfies(Standard_Complex_Solutions.Head_Of(sols).all);
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
      Standard_Complex_Solutions.Shallow_Clear(sols);
    end loop;
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsol := Standard_Complex_Solutions.Head_Of(sols).all;
    Standard_Complex_Solutions.Shallow_Clear(sols);
  end Sample;

  procedure Sample ( file : in file_type; full_output : in boolean;
                     startsol : in Standard_Complex_Solutions.Solution;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsol : out Standard_Complex_Solutions.Solution ) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter,cnt_reruns : natural32;
    fail : boolean;

  begin
    cnt_reruns := 0;
    loop
      Standard_Complex_Solutions.Add(sols,startsol);
      if full_output
       then Reporting_Homotopy_Continuation
              (file,target_eval,target_jaco,dim,sols);
            Standard_Reporting_Correct
              (file,target_eval,target_jaco,sols,iter,fail);
       else Silent_Homotopy_Continuation(target_eval,target_jaco,dim,sols);
            Standard_Silent_Correct(target_eval,target_jaco,sols,iter,fail);
      end if;
     -- exit when (Standard_Complex_Solutions.Head_Of(sols).rco > tolsing);
      exit when Satisfies(file,Standard_Complex_Solutions.Head_Of(sols).all);
      cnt_reruns := cnt_reruns + 1;
      exit when (cnt_reruns > Continuation_Parameters.max_reruns);
      Standard_Complex_Solutions.Shallow_Clear(sols);
    end loop;
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsol := Standard_Complex_Solutions.Head_Of(sols).all;
    Standard_Complex_Solutions.Shallow_Clear(sols);
  end Sample;

  procedure Sample ( startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    Standard_Complex_Solutions.Copy(startsols,sols);
    Silent_Homotopy_Continuation(target_eval,target_jaco,dim,sols);
    Standard_Silent_Correct(target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample ( startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     gamma : in Standard_Complex_Vectors.Vector;
                     newsols : out Standard_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    Standard_Complex_Solutions.Copy(startsols,sols);
    Silent_Homotopy_Continuation(target_eval,target_jaco,dim,gamma,sols);
    Standard_Silent_Correct(target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample ( file : in file_type;
                     startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    Standard_Complex_Solutions.Copy(startsols,sols);
    Reporting_Homotopy_Continuation(file,target_eval,target_jaco,dim,sols);
    Standard_Reporting_Correct(file,target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample ( file : in file_type;
                     startsols : in Standard_Complex_Solutions.Solution_List;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     gamma : in Standard_Complex_Vectors.Vector;
                     newsols : out Standard_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    Standard_Complex_Solutions.Copy(startsols,sols);
    Reporting_Homotopy_Continuation
      (file,target_eval,target_jaco,dim,gamma,sols);
    Standard_Reporting_Correct(file,target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample;

  procedure Sample_with_Stop
                   ( startsols : in Standard_Complex_Solutions.Solution_List;
                     x : in Standard_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    Standard_Complex_Solutions.Copy(startsols,sols);
    Silent_Homotopy_Continuation_with_Stop
      (target_eval,target_jaco,dim,x,tol,sols);
    Standard_Silent_Correct(target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample_with_Stop;

  procedure Sample_with_Stop
                   ( file : in file_type;
                     startsols : in Standard_Complex_Solutions.Solution_List;
                     x : in Standard_Complex_Vectors.Vector;
                     tol : in double_float;
                     newhyp : in Standard_Complex_VecVecs.VecVec;
                     newsols : out Standard_Complex_Solutions.Solution_List) is

    dim : constant integer32 := newhyp'last;
    target : Standard_Complex_Laur_Systems.Laur_Sys(1..Number_of_Equations)
           := Embed(Standard_System,newhyp);
    target_eval : Standard_Complex_Laur_SysFun.Eval_Laur_Sys(target'range)
                := Standard_Evaluator(target,dim);
    target_jaco : Standard_Complex_Laur_JacoMats.Eval_Jaco_Mat
                    (target'range,target'range)
                := Standard_Jacobi_Evaluator(target,dim);
    sols : Standard_Complex_Solutions.Solution_List;
    iter : natural32;
    fail : boolean;

  begin
    Standard_Complex_Solutions.Copy(startsols,sols);
    Reporting_Homotopy_Continuation_with_Stop
      (file,target_eval,target_jaco,dim,x,tol,sols);
    Standard_Reporting_Correct(file,target_eval,target_jaco,sols,iter,fail);
    for i in target'last-dim+1..target'last loop       -- mind the sharing !
      Standard_Complex_Laurentials.Clear(target(i));
      Standard_Complex_Laur_Functions.Clear(target_eval(i));
      for j in target'range loop
        Standard_Complex_Laur_Functions.Clear(target_jaco(i,j));
      end loop;
    end loop;
    newsols := sols;
  end Sample_with_Stop;

-- REFINERS :

  procedure Refine ( stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution;
                     mphyp : out Multprec_Complex_VecVecs.VecVec ) is

    iter : natural32;
    fail : boolean;

  begin
    mpsol := Multprec_Complex_Solutions.Create(stsol);
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    for i in sthyp'range loop
      mphyp(i) := new Multprec_Complex_Vectors.Vector'(Create(sthyp(i).all));
      Set_Size(mphyp(i).all,mpsize);
    end loop;
    Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
  end Refine;

  procedure Refine ( file : in file_type; full_output : in boolean;
                     stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution;
                     mphyp : out Multprec_Complex_VecVecs.VecVec ) is

    iter : natural32;
    fail : boolean;

  begin
    mpsol := Multprec_Complex_Solutions.Create(stsol);
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    for i in sthyp'range loop
      mphyp(i) := new Multprec_Complex_Vectors.Vector'(Create(sthyp(i).all));
      Set_Size(mphyp(i).all,mpsize);
    end loop;
    if full_output
     then Multprec_Reporting_Correct(file,mpsol,mphyp,iter,fail);
     else Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
    end if;
    Write_Summary(file,iter,mpsol);
  end Refine;

  procedure Refine ( mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in out Multprec_Complex_VecVecs.VecVec ) is

    iter : natural32;
    fail : boolean;

  begin
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    for i in mphyp'range loop
      Set_Size(mphyp(i).all,mpsize);
    end loop;
    Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
  end Refine;

  procedure Refine ( file : in file_type; full_output : in boolean;
                     mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in out Multprec_Complex_VecVecs.VecVec ) is

    iter : natural32;
    fail : boolean;

  begin
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    for i in mphyp'range loop
      Set_Size(mphyp(i).all,mpsize);
    end loop;
    if full_output
     then Multprec_Reporting_Correct(file,mpsol,mphyp,iter,fail);
     else Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
    end if;
    Write_Summary(file,iter,mpsol);
  end Refine;

  procedure Refine_on_Slices
                   ( stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution ) is

    iter : natural32;
    fail : boolean;

  begin
    mpsol := Multprec_Complex_Solutions.Create(stsol);
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
  end Refine_on_Slices;

  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     stsol : in Standard_Complex_Solutions.Solution;
                     sthyp : in Standard_Complex_VecVecs.VecVec;
                     mphyp : in Multprec_Complex_VecVecs.VecVec;
                     mpsol : out Multprec_Complex_Solutions.Solution ) is

    iter : natural32;
    fail : boolean;

  begin
    mpsol := Multprec_Complex_Solutions.Create(stsol);
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    if full_output
     then Multprec_Reporting_Correct(file,mpsol,mphyp,iter,fail);
     else Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
    end if;
    Write_Summary(file,iter,mpsol);
  end Refine_on_Slices;

  procedure Refine_on_Slices
                   ( mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in Multprec_Complex_VecVecs.VecVec ) is

    iter : natural32;
    fail : boolean;

  begin
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
  end Refine_on_Slices;

  procedure Refine_on_Slices
                   ( file : in file_type; full_output : in boolean;
                     mpsol : in out Multprec_Complex_Solutions.Solution;
                     mphyp : in Multprec_Complex_VecVecs.VecVec ) is

    iter : natural32;
    fail : boolean;

  begin
    Multprec_Complex_Solutions.Set_Size(mpsol,mpsize);
    if full_output
     then Multprec_Reporting_Correct(file,mpsol,mphyp,iter,fail);
     else Multprec_Silent_Correct(mpsol,mphyp,iter,fail);
    end if;
    Write_Summary(file,iter,mpsol);
  end Refine_on_Slices;

-- DEALLOCATION :

  procedure Clear is
  begin
   -- Standard_Complex_Laur_Systems.Shallow_Clear(stansys); -- data sharing
    Standard_Complex_Laur_Systems.Clear(stansys); -- no sharing!
    Multprec_Complex_Laur_Systems.Shallow_Clear(multsys);
    Multprec_Complex_Laur_Systems.Shallow_Clear(orgsys);
    Standard_Complex_Laur_SysFun.Clear(stansys_eval);
    Standard_Complex_Laur_JacoMats.Clear(stanjac_eval);
    Multprec_Complex_Laur_SysFun.Clear(multsys_eval);
    Multprec_Complex_Laur_JacoMats.Clear(multjac_eval);
  end Clear;

  procedure Clear_Restricted is
  begin
    restricted := false;
   -- Standard_Complex_Laur_Systems.Shallow_Clear(rststansys); -- data sharing
    Standard_Complex_Laur_Systems.Clear(rststansys); -- no sharing!
    Multprec_Complex_Laur_Systems.Shallow_Clear(rstmultsys);
    Standard_Complex_Laur_SysFun.Clear(rststansys_eval);
    Standard_Complex_Laur_JacoMats.Clear(rststanjac_eval);
    Multprec_Complex_Laur_SysFun.Clear(rstmultsys_eval);
    Multprec_Complex_Laur_JacoMats.Clear(rstmultjac_eval);
  end Clear_Restricted;

end Sampling_Laurent_Machine;
