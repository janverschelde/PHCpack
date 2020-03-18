with unchecked_deallocation;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Complex_Singular_Values;
with Standard_Mixed_Residuals;
with Standard_Rational_Approximations;
with Newton_Convolutions;
with Newton_Power_Convolutions;
with Convergence_Radius_Estimates;
with Jacobian_Convolution_Circuits;
with Hessian_Convolution_Circuits;
with Homotopy_Pade_Approximants;
with Series_and_Predictors;

package body Standard_Predictor_Convolutions is

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) return LU_Predictor is

    res : LU_Predictor(sol'last,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(0.0);
    knum : constant Standard_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant Standard_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

  begin
    res.sol := Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new Standard_Complex_Vectors.Vector'(1..neq => zero);
    for k in sol'range loop
      res.numcff(k) := new Standard_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new Standard_Complex_Vectors.Vector'(kden);
    end loop;
    return res;
  end Create;

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 )
                  return SVD_Predictor is

    dim : constant integer32 := sol'last;
    res : SVD_Predictor(neq,dim,dim+1,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(0.0);
    knum : constant Standard_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant Standard_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

  begin
    res.sol := Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new Standard_Complex_Vectors.Vector'(1..neq => zero);
    res.ewrk := new Standard_Complex_Vectors.Vector'(1..dim => zero);
    res.dx := Allocate_Coefficients(dim,deg);
    res.xd := Linearized_Allocation(dim,deg);
    for k in sol'range loop
      res.numcff(k) := new Standard_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new Standard_Complex_Vectors.Vector'(kden);
    end loop;
    res.svl := (res.svl'range => zero);
    return res;
  end Create;

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) 
                  return Link_to_LU_Predictor is

    prd : constant LU_Predictor(sol'last,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg);
    res : constant Link_to_LU_Predictor := new LU_Predictor'(prd);

  begin
    return res;
  end Create;

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) 
                  return Link_to_SVD_Predictor is

    dim : constant integer32 := sol'last;
    prd : constant SVD_Predictor(neq,dim,dim+1,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg);
    res : constant Link_to_SVD_Predictor := new SVD_Predictor'(prd);

  begin
    return res;
  end Create;

  procedure Newton_Fabry_Report 
              ( file : in file_type;
                nbrit : in integer32; absdx : in double_float;
                fail : in boolean;
                z : in Standard_Complex_Numbers.Complex_Number;
                r,err,step : in double_float;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                output : in boolean ) is
  begin
    put(file,"#iterations : "); put(file,nbrit,1);
    put(file,"  |dx| :"); put(file,absdx,3); new_line(file);
    if fail then
      put_line(file,"Predictor failed!");
    else
      put(file,"z : "); put(file,z); 
      put(file,"  error estimate :"); put(file,err,3); new_line(file);
      put(file,"estimated radius :"); put(file,r,3);
    end if;
    put(file,"  pole step :"); put(file,step,3); new_line(file);
    if output then
      for k in numcff'range loop
        put(file,"Numerator coefficients at "); put(file,k,1);
        put_line(file," :"); put_line(file,numcff(k));
        put(file,"Denominator coefficients at "); put(file,k,1);
        put_line(file," :"); put_line(file,dencff(k));
      end loop;
    end if;
  end Newton_Fabry_Report;

  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Link_to_System; prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean ) is

    use Standard_Rational_Approximations;
    use Newton_Power_Convolutions;

    info : integer32;

  begin
    nbrit := 0;
    if output then
      LU_Newton_Steps
        (file,hom,prd.sol,maxit,nbrit,tol,absdx,fail,
         info,prd.newtpiv,prd.wrk,false);
      Convergence_Radius_Estimates.Fabry(file,prd.sol,z,rad,err,fail,2);
    else
      LU_Newton_Steps
        (hom,prd.sol,maxit,nbrit,tol,absdx,fail,
         info,prd.newtpiv,prd.wrk,false,false);
      Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    end if;
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Link_to_System; prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean ) is

    use Standard_Rational_Approximations;
    use Newton_Power_Convolutions;

    info : integer32;

  begin
    nbrit := 0;
    if output then
      SVD_Newton_Steps
        (file,hom,prd.sol,prd.dx,prd.xd,maxit,nbrit,tol,absdx,
         fail, prd.svl,prd.U,prd.V,info,rcond,prd.ewrk,prd.wrk,false);
      Convergence_Radius_Estimates.Fabry(file,prd.sol,z,rad,err,fail,2);
    else
      SVD_Newton_Steps
        (hom,prd.sol,prd.dx,prd.xd,maxit,nbrit,tol,absdx,fail,
         prd.svl,prd.U,prd.V,info,rcond,prd.ewrk,prd.wrk,false,false);
      Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    end if;
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

  procedure Second
              ( hom : in Link_to_System; svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := svh.dim;
    p : constant integer32 := svh.dim;
    job : constant integer32 := 11;
    info : integer32;

    use Standard_Complex_Singular_Values;

  begin
    for k in hom.crc'range loop
      svh.H := Hessian_Convolution_Circuits.Hessian(hom.crc(k),sol);
      SVD(svh.H,n,p,svh.svl,svh.ewrk,svh.U,svh.V,job,info);
      svh.vals(k) := svh.svl(1);
    end loop;
  end Second;

  function Distance ( svh : in Link_to_SVD_Hessians ) return double_float is

    sigma1 : constant double_float
           := Standard_Complex_Numbers.REAL_PART(svh.vals(0));
    accsum,acc,nrm : double_float := 0.0;

  begin
    for k in 1..svh.dim loop
      acc := Standard_Complex_Numbers.REAL_PART(svh.vals(k));
      accsum := accsum + acc*acc; 
    end loop;
    nrm := Standard_Mathematical_Functions.SQRT(accsum);
    return (2.0*sigma1)/nrm;
  end Distance;

  procedure Hesse_Pade
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Standard_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true ) is

    info : integer32;

    use Standard_Complex_Singular_Values;

  begin -- with LU, the system is square so svh.H work space works
    svh.H := Jacobian_Convolution_Circuits.Jacobian(hom.crc,sol);
    SVD(svh.H,svh.dim,svh.dim,svh.svl,svh.ewrk,svh.U,svh.V,11,info);
    svh.vals(0) := svh.svl(svh.dim);
    Standard_Predictor_Convolutions.Second(hom,svh,sol);
    if verbose
     then put_line(file,"All singular values : "); put_line(file,svh.vals);
    end if;
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := Standard_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,beta2,eta,nrm);
    if verbose then
      put(file,"eta :"); put(file,eta,3);
      put(file,"  nrm :"); put(file,nrm,3);
      put(file,"  curv_step :"); put(file,step,3); new_line(file);
    end if;
  end Hesse_Pade;

  procedure Hesse_Pade
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Standard_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true ) is
  begin
    svh.vals(0) := prd.svl(prd.dim);
    Standard_Predictor_Convolutions.Second(hom,svh,sol);
    if verbose
     then put_line(file,"All singular values : "); put_line(file,svh.vals);
    end if;
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := Standard_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,beta2,eta,nrm);
    if verbose then
      put(file,"eta :"); put(file,eta,3);
      put(file,"  nrm :"); put(file,nrm,3);
      put(file,"  curv_step :"); put(file,step,3); new_line(file);
    end if;
  end Hesse_Pade;

  procedure Predictor_Feedback
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; alpha : in double_float;
                eva,radsol : in out Standard_Complex_Vectors.Vector;
                res,absres : in out Standard_Complex_Vectors.Vector;
                nrm,mixres : out double_float; nbfail : out integer32;
                verbose : in boolean := true ) is

    z : Standard_Complex_Numbers.Complex_Number;

  begin
    nbfail := 0;
    loop
      if verbose then
        put(file,"step in the predictor feedback loop :");
        put(file,step,3); new_line(file);
      end if;
      Standard_Rational_Approximations.Evaluate(numcff,dencff,step,eva);
      z := Standard_Complex_Numbers.Create(step);
      res := Standard_Speelpenning_Convolutions.Eval(hom.crc,eva,z);
      nrm := Standard_Complex_Vector_Norms.Max_Norm(res);
      radsol := Standard_Mixed_Residuals.AbsVal(eva);
      absres := Standard_Speelpenning_Convolutions.Eval(abh.crc,radsol,z);
      mixres := Standard_Mixed_Residuals.Mixed_Residual(res,absres);
      if verbose then
        put_line(file,"Evaluation of the predicted solution : ");
        put_line(file,res);
        put(file,"The predictor residual :"); put(file,nrm,3);
        put(file,"  mixres :"); put(file,mixres,3);
      end if;
      if mixres < alpha then
        if verbose
         then put_line(file,"  okay");
        end if;
        exit;
      else
        if verbose
         then put(file," >"); put(file,alpha,3); new_line(file);
        end if;
        step := step/2.0; nbfail := nbfail + 1;
      end if;
    end loop;
  end Predictor_Feedback;

-- DESTRUCTORS :

  procedure Clear ( p : in out Link_to_LU_Predictor ) is

    procedure free is
      new unchecked_deallocation(LU_Predictor,Link_to_LU_Predictor);

  begin
    if p /= null then
      Standard_Complex_VecVecs.Clear(p.sol);
      Standard_Complex_Vectors.Clear(p.wrk);
      Standard_Complex_VecVecs.Clear(p.numcff);
      Standard_Complex_VecVecs.Clear(p.dencff);
      free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Link_to_SVD_Predictor ) is

    procedure free is
      new unchecked_deallocation(SVD_Predictor,Link_to_SVD_Predictor);

  begin
    if p /= null then
      Standard_Complex_VecVecs.Clear(p.sol);
      Standard_Complex_Vectors.Clear(p.wrk);
      Standard_Complex_Vectors.Clear(p.ewrk);
      Standard_Complex_VecVecs.Clear(p.dx);
      Standard_Complex_VecVecs.Clear(p.xd);
      Standard_Complex_VecVecs.Clear(p.numcff);
      Standard_Complex_VecVecs.Clear(p.dencff);
      free(p);
    end if;
  end Clear;

  procedure Clear ( h : in out Link_to_SVD_Hessians ) is

    procedure free is
      new unchecked_deallocation(SVD_Hessians,Link_to_SVD_Hessians);

  begin
    if h /= null
     then free(h);
    end if;
  end Clear;

end Standard_Predictor_Convolutions;
