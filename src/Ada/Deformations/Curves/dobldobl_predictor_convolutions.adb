with unchecked_deallocation;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Mathematical_Functions;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Singular_Values;
with DoblDobl_Mixed_Residuals;
with DoblDobl_Rational_Approximations;
with DoblDobl_Newton_Convolutions;
with DoblDobl_Newton_Convolution_Steps;
with Convergence_Radius_Estimates;
with Jacobian_Convolution_Circuits;
with Hessian_Convolution_Circuits;
with Homotopy_Pade_Approximants;
with Series_and_Predictors;
with Three_Way_Minima;

package body DoblDobl_Predictor_Convolutions is

-- CONSTRUCTORS :

  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) return LU_Predictor is

    res : LU_Predictor(sol'last,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(integer(0));
    knum : constant DoblDobl_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant DoblDobl_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

  begin
    res.sol := DoblDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new DoblDobl_Complex_Vectors.Vector'(1..neq => zero);
    for k in sol'range loop
      res.numcff(k) := new DoblDobl_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new DoblDobl_Complex_Vectors.Vector'(kden);
    end loop;
    return res;
  end Create;

  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 )
                  return SVD_Predictor is

    dim : constant integer32 := sol'last;
    res : SVD_Predictor(neq,dim,dim+1,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(integer(0));
    knum : constant DoblDobl_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant DoblDobl_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

  begin
    res.sol := DoblDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new DoblDobl_Complex_Vectors.Vector'(1..neq => zero);
    res.ewrk := new DoblDobl_Complex_Vectors.Vector'(1..dim => zero);
    res.dx := Allocate_Coefficients(dim,deg);
    res.xd := Linearized_Allocation(dim,deg);
    for k in sol'range loop
      res.numcff(k) := new DoblDobl_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new DoblDobl_Complex_Vectors.Vector'(kden);
    end loop;
    res.svl := (res.svl'range => zero);
    return res;
  end Create;

  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) 
                  return Link_to_LU_Predictor is

    prd : constant LU_Predictor(sol'last,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg);
    res : constant Link_to_LU_Predictor := new LU_Predictor'(prd);

  begin
    return res;
  end Create;

  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32 ) 
                  return Link_to_SVD_Predictor is

    dim : constant integer32 := sol'last;
    prd : constant SVD_Predictor(neq,dim,dim+1,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg);
    res : constant Link_to_SVD_Predictor := new SVD_Predictor'(prd);

  begin
    return res;
  end Create;

  function Create ( p : Link_to_LU_Predictor ) return Predictor is

    res : Predictor(LU);

  begin
    res.ludata := p;
    return res;
  end Create;

  function Create ( p : Link_to_SVD_Predictor ) return Predictor is

    res : Predictor(SVD);

  begin
    res.svdata := p;
    return res;
  end Create;

  function Create ( sol : DoblDobl_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type ) return Predictor is
  begin
    case kind is
      when LU =>
        declare
          res : Predictor(LU);
        begin
          res.ludata := Create(sol,neq,deg,numdeg,dendeg);
          return res;
        end;
      when SVD =>
        declare
          res : Predictor(SVD);
        begin
          res.svdata := Create(sol,neq,deg,numdeg,dendeg);
          return res;
        end;
    end case;
  end Create;

  function Create ( nbr : integer32;
                    sol : DoblDobl_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type ) return Predictor_Array is

    res : Predictor_Array(1..nbr);

  begin
    for k in 1..nbr loop
      res(k) := Create(sol,neq,deg,numdeg,dendeg,kind);
    end loop;
    return res;
  end Create;

  function Create ( nbr,dim,neq : integer32 )
                  return Predictor_Vectors_Array is

    res : Predictor_Vectors_Array(1..nbr);

  begin
    for k in 1..nbr loop
      res(k) := new Predictor_Vectors(dim,neq);
    end loop;
    return res;
  end Create;

  function Create ( dim : integer32 ) return Link_to_SVD_Hessians is

    res : Link_to_SVD_Hessians;
    hss : SVD_Hessians(dim,dim+1);
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));

  begin
    hss.vals := (hss.vals'range => zero);
    res := new SVD_Hessians'(hss);
    return res;
  end Create;

  function Create ( nbr,dim : integer32 ) return SVD_Hessians_Array is

    res : SVD_Hessians_Array(1..nbr);

  begin
    for k in 1..nbr loop
      res(k) := Create(dim);
    end loop;
    return res;
  end Create;

-- AUXILIARY PREDICTOR PROCEDURES FOR SETUP :

  procedure Set_Lead_Coefficients
              ( p : in Predictor;
                s : in DoblDobl_Complex_Vectors.Vector ) is

    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));

  begin
    case p.kind is
      when LU =>
        for k in p.ludata.sol'range loop
          lnk := p.ludata.sol(k);
          lnk(0) := s(k);
          for i in 1..lnk'last loop
            lnk(i) := zero;
          end loop;
        end loop;
      when SVD =>
        for k in p.svdata.sol'range loop
          lnk := p.svdata.sol(k);
          lnk(0) := s(k);
          for i in 1..lnk'last loop
            lnk(i) := zero;
          end loop;
        end loop;
    end case;
  end Set_Lead_Coefficients;

  procedure Newton_Fabry_Report 
              ( file : in file_type;
                nbrit : in integer32; absdx : in double_double;
                fail : in boolean;
                z : in DoblDobl_Complex_Numbers.Complex_Number;
                r,err,step : in double_double;
                numcff,dencff : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean ) is
  begin
    put(file,"#iterations : "); put(file,nbrit,1);
    put(file,"  |dx| : "); put(file,absdx,3); new_line(file);
    if fail then
      put_line(file,"Predictor failed!");
    else
      put(file,"z : "); put(file,z); 
      put(file,"  error estimate : "); put(file,err,3); new_line(file);
      put(file,"estimated radius : "); put(file,r,3);
    end if;
    put(file,"  pole step : "); put(file,step,3); new_line(file);
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
              ( hom : in Link_to_System; prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_double ) is

    use DoblDobl_Rational_Approximations;
    use DoblDobl_Newton_Convolution_Steps;

    info : integer32;

  begin
    nbrit := 0;
    LU_Newton_Steps
      (hom,prd.sol,maxit,nbrit,tol,absdx,fail,
       info,prd.newtpiv,prd.wrk,false,false);
    Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Link_to_System; prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_double; output : in boolean ) is

    use DoblDobl_Rational_Approximations;
    use DoblDobl_Newton_Convolution_Steps;

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
              ( hom : in Link_to_System; prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_double ) is

    use DoblDobl_Rational_Approximations;
    use DoblDobl_Newton_Convolution_Steps;

    info : integer32;

  begin
    nbrit := 0;
    SVD_Newton_Steps
      (hom,prd.sol,prd.dx,prd.xd,maxit,nbrit,tol,absdx,fail,
       prd.svl,prd.U,prd.V,info,rcond,prd.ewrk,prd.wrk,false,false);
    Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Link_to_System; prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_double;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_double; output : in boolean ) is

    use DoblDobl_Rational_Approximations;
    use DoblDobl_Newton_Convolution_Steps;

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
                sol : in DoblDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := svh.dim;
    p : constant integer32 := svh.dim;
    info : integer32;

    use DoblDobl_Complex_Singular_Values;

  begin
    for k in hom.crc'range loop
      svh.H := Hessian_Convolution_Circuits.Hessian(hom.crc(k),sol);
      SVD(svh.H,n,p,svh.svl,svh.ewrk,svh.U,svh.V,0,info,svh.work);
      svh.vals(k) := svh.svl(1);
    end loop;
  end Second;

  function Distance ( svh : in Link_to_SVD_Hessians ) return double_double is

    sigma1 : constant double_double
           := DoblDobl_Complex_Numbers.REAL_PART(svh.vals(0));
    accsum,acc,nrm : double_double := create(0.0);

  begin
    for k in 1..svh.dim loop
      acc := DoblDobl_Complex_Numbers.REAL_PART(svh.vals(k));
      accsum := accsum + acc*acc; 
    end loop;
    nrm := DoblDobl_Mathematical_Functions.SQRT(accsum);
    return (2.0*sigma1)/nrm;
  end Distance;

  procedure Hesse_Pade
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in DoblDobl_Complex_Vectors.Vector;
                res : out DoblDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_double ) is

    info : integer32;
    dd_beta2 : constant double_double := create(beta2);

    use DoblDobl_Complex_Singular_Values;

  begin -- with LU, the system is square so svh.H work space works
    svh.H := Jacobian_Convolution_Circuits.Jacobian(hom.crc,sol);
    SVD(svh.H,svh.dim,svh.dim,svh.svl,svh.ewrk,svh.U,svh.V,0,info,svh.work);
    svh.vals(0) := svh.svl(svh.dim);
    Second(hom,svh,sol);
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := DoblDobl_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,dd_beta2,eta,nrm);
  end Hesse_Pade;

  procedure Hesse_Pade
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in DoblDobl_Complex_Vectors.Vector;
                res : out DoblDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_double;
                verbose : in boolean := true ) is

    info : integer32;
    dd_beta2 : constant double_double := create(beta2);

    use DoblDobl_Complex_Singular_Values;

  begin -- with LU, the system is square so svh.H work space works
    svh.H := Jacobian_Convolution_Circuits.Jacobian(hom.crc,sol);
    SVD(svh.H,svh.dim,svh.dim,svh.svl,svh.ewrk,svh.U,svh.V,0,info,svh.work);
    svh.vals(0) := svh.svl(svh.dim);
    Second(hom,svh,sol);
    if verbose
     then put_line(file,"All singular values : "); put_line(file,svh.vals);
    end if;
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := DoblDobl_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,dd_beta2,eta,nrm);
    if verbose then
      put(file,"eta : "); put(file,eta,3);
      put(file,"  nrm : "); put(file,nrm,3);
      put(file,"  curv_step : "); put(file,step,3); new_line(file);
    end if;
  end Hesse_Pade;

  procedure Hesse_Pade
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in DoblDobl_Complex_Vectors.Vector;
                res : out DoblDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_double ) is

    dd_beta2 : constant double_double := create(beta2);

  begin
    svh.vals(0) := prd.svl(prd.dim);
    Second(hom,svh,sol);
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := DoblDobl_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,dd_beta2,eta,nrm);
  end Hesse_Pade;

  procedure Hesse_Pade
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                sol : in DoblDobl_Complex_Vectors.Vector;
                res : out DoblDobl_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_double;
                verbose : in boolean := true ) is

    dd_beta2 : constant double_double := create(beta2);

  begin
    svh.vals(0) := prd.svl(prd.dim);
    Second(hom,svh,sol);
    if verbose
     then put_line(file,"All singular values : "); put_line(file,svh.vals);
    end if;
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := DoblDobl_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,dd_beta2,eta,nrm);
    if verbose then
      put(file,"eta : "); put(file,eta,3);
      put(file,"  nrm : "); put(file,nrm,3);
      put(file,"  curv_step : "); put(file,step,3); new_line(file);
    end if;
  end Hesse_Pade;

  procedure Predictor_Feedback
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in DoblDobl_Complex_VecVecs.VecVec;
                step : in out double_double; minstep,alpha : in double_float;
                nrm,mixres : out double_double; nbfail : out integer32 ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    nbfail := 0;
    loop
      DoblDobl_Rational_Approximations.Evaluate(numcff,dencff,step,psv.sol);
      z := DoblDobl_Complex_Numbers.Create(step);
      psv.res := Eval(hom.crc,psv.sol,z);
      nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(psv.res);
      psv.radsol := DoblDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol,z);
      mixres := DoblDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if mixres < alpha then
        exit;
      else
        step := step/2.0; nbfail := nbfail + 1;
        exit when (step < minstep);
      end if;
    end loop;
  end Predictor_Feedback;

  procedure Predictor_Feedback
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in DoblDobl_Complex_VecVecs.VecVec;
                step : in out double_double; minstep,alpha : in double_float;
                nrm,mixres : out double_double; nbfail : out integer32;
                verbose : in boolean := true ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    nbfail := 0;
    loop
      if verbose then
        put(file,"step in the feedback loop : ");
        put(file,step,3); new_line(file);
      end if;
      DoblDobl_Rational_Approximations.Evaluate(numcff,dencff,step,psv.sol);
      z := DoblDobl_Complex_Numbers.Create(step);
      psv.res := Eval(hom.crc,psv.sol,z);
      nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(psv.res);
      psv.radsol := DoblDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol,z);
      mixres := DoblDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put_line(file,"Evaluation of the predicted solution : ");
        put_line(file,psv.res);
        put(file,"The predictor residual : "); put(file,nrm,3);
        put(file,"  mixres : "); put(file,mixres,3);
      end if;
      if mixres < alpha then
        if verbose
         then put_line(file,"  okay");
        end if;
        exit;
      else
        if verbose
         then put(file," > "); put(file,alpha,3); new_line(file);
        end if;
        step := step/2.0; nbfail := nbfail + 1;
        exit when (step < minstep);
      end if;
    end loop;
  end Predictor_Feedback;

-- MAIN PREDICTOR PROCEDURES :

  procedure LU_Prediction
              ( hom,abh : in Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_double;
                fail : out boolean; step : out double_double;
                nbpole,nbhess,nbmaxm : in out natural32 ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,mixres : double_double;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;
    dd_maxstep : constant double_double := create(maxstep);

    use Three_Way_Minima;

  begin
    Newton_Fabry(hom,prd,maxit,tol,nbrit,absdx,fail,z,r,err);
    if fail
     then pole_step := dd_maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Hesse_Pade(hom,prd,svh,psv.sol,psv.res,beta2,eta,nrm,curv_step);
    Minimum(pole_step,curv_step,dd_maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(hom,abh,psv,prd.numcff,prd.dencff,
      step,minstep,alpha,nrm,mixres,nbfail);
    fail := (mixres > alpha);
  end LU_Prediction;

  procedure LU_Prediction
              ( file : in file_type; hom,abh : in Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_double;
                fail : out boolean; step : out double_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,mixres : double_double;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;
    dd_maxstep : constant double_double := create(maxstep);

    use Three_Way_Minima;

  begin
    Newton_Fabry(file,hom,prd,maxit,tol,nbrit,absdx,fail,z,r,err,output);
    if fail
     then pole_step := dd_maxstep;
     else pole_step := beta1*r;
    end if;
    if verbose then
      Newton_Fabry_Report(file,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Hesse_Pade(file,hom,prd,svh,psv.sol,psv.res,beta2,
               eta,nrm,curv_step,verbose);
    Minimum(pole_step,curv_step,dd_maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(file,hom,abh,psv,prd.numcff,prd.dencff,
      step,minstep,alpha,nrm,mixres,nbfail,verbose);
    fail := (mixres > alpha);
  end LU_Prediction;

  procedure SVD_Prediction
              ( file : in file_type; hom,abh : in Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_double;
                fail : out boolean; step : out double_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,mixres : double_double;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;
    dd_maxstep : constant double_double := create(maxstep);

    use Three_Way_Minima;

  begin
    Newton_Fabry(file,hom,prd,maxit,tol,nbrit,absdx,rcond,fail,z,r,err,output);
    if fail
     then pole_step := dd_maxstep;
     else pole_step := beta1*r;
    end if;
    if verbose then
      Newton_Fabry_Report(file,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Hesse_Pade(file,hom,prd,svh,psv.sol,psv.res,beta2,
               eta,nrm,curv_step,verbose);
    Minimum(pole_step,curv_step,dd_maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(file,hom,abh,psv,prd.numcff,prd.dencff,
      step,minstep,alpha,nrm,mixres,nbfail,verbose);
    fail := (mixres > alpha);
  end SVD_Prediction;

  procedure SVD_Prediction
              ( hom,abh : in Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_double;
                fail : out boolean; step : out double_double;
                nbpole,nbhess,nbmaxm : in out natural32 ) is

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,mixres : double_double;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;
    dd_maxstep : constant double_double := create(maxstep);

    use Three_Way_Minima;

  begin
    Newton_Fabry(hom,prd,maxit,tol,nbrit,absdx,rcond,fail,z,r,err);
    if fail
     then pole_step := dd_maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Hesse_Pade(hom,prd,svh,psv.sol,psv.res,beta2,eta,nrm,curv_step);
    Minimum(pole_step,curv_step,dd_maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(hom,abh,psv,prd.numcff,prd.dencff,
      step,minstep,alpha,nrm,mixres,nbfail);
    fail := (mixres > alpha);
  end SVD_Prediction;

-- DESTRUCTORS :

  procedure Clear ( p : in out Link_to_LU_Predictor ) is

    procedure free is
      new unchecked_deallocation(LU_Predictor,Link_to_LU_Predictor);

  begin
    if p /= null then
      DoblDobl_Complex_VecVecs.Clear(p.sol);
      DoblDobl_Complex_VecVecs.Clear(p.numcff);
      DoblDobl_Complex_VecVecs.Clear(p.dencff);
      free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Link_to_SVD_Predictor ) is

    procedure free is
      new unchecked_deallocation(SVD_Predictor,Link_to_SVD_Predictor);

  begin
    if p /= null then
      DoblDobl_Complex_VecVecs.Clear(p.sol);
      DoblDobl_Complex_Vectors.Clear(p.wrk);
      DoblDobl_Complex_Vectors.Clear(p.ewrk);
      DoblDobl_Complex_VecVecs.Clear(p.dx);
      DoblDobl_Complex_VecVecs.Clear(p.xd);
      DoblDobl_Complex_VecVecs.Clear(p.numcff);
      DoblDobl_Complex_VecVecs.Clear(p.dencff);
      free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Predictor ) is
  begin
    case p.kind is
      when LU  => Clear(p.ludata);
      when SVD => Clear(p.svdata);
    end case;
  end Clear;

  procedure Clear ( p : in out Predictor_Array ) is
  begin
    for i in p'range loop
      case p(i).kind is
        when LU => Clear(p(i).ludata);
        when SVD => Clear(p(i).svdata);
      end case;
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Predictor_Vectors ) is

    procedure free is
      new unchecked_deallocation(Predictor_Vectors,Link_to_Predictor_Vectors);

  begin
    if p /= null
     then free(p);
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

  procedure Clear ( p : in out LU_Predictor_Array ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out SVD_Predictor_Array ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Predictor_Vectors_Array ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( h : in out SVD_Hessians_Array ) is
  begin
    for i in h'range loop
      Clear(h(i));
    end loop;
  end Clear;

end DoblDobl_Predictor_Convolutions;
