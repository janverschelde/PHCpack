with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Solution_Drops;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with Series_and_Homotopies;
with Test_Series_Predictors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Three_Way_Minima;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;

procedure ts_padepcnv is

-- DESCRIPTION :
--   Development of the Pade predictor on convolution circuits.

  procedure Standard_LU_Prediction
              ( hom,abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Standard_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep : in double_float;
                fail : out boolean; nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coeffiecients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations;
  --   verbose  flag for intermediate numerical output.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status;
  --   nbpole   updated number of times pole step was minimal;
  --   nbhess   updated number of times curve step was minimal;
  --   nbmaxm   updated number of times maximum step size was minimal.

    use Standard_Speelpenning_Convolutions;
    use Standard_Predictor_Convolutions;

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    sol,radsol : Standard_Complex_Vectors.Vector(1..prd.dim);
    res,absres : Standard_Complex_Vectors.Vector(hom.crc'range);
    nbrit,nbfail : integer32;

  begin
    Newton_Fabry(standard_output,hom,prd,maxit,tol,nbrit,absdx,fail,
                 z,r,err,output);
    pole_step := beta1*r;
    if verbose then
      Newton_Fabry_Report(standard_output,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    Hesse_Pade(standard_output,hom,prd,svh,sol,res,beta2,
               eta,nrm,curv_step,verbose);
    Three_Way_Minima.Minimum
      (pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Predictor_Feedback(standard_output,hom,abh,prd.numcff,prd.dencff,step,
      alpha,sol,radsol,res,absres,nrm,mixres,nbfail,verbose);
  end Standard_LU_Prediction;

  procedure Standard_SVD_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Standard_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep : in double_float;
                fail : out boolean; nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations;
  --   verbose  flag for intermediate numerical output.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status;
  --   nbpole   updated number of times pole step was minimal;
  --   nbhess   updated number of times curve step was minimal;
  --   nbmaxm   updated number of times maximum step size was minimal.

    use Standard_Speelpenning_Convolutions;
    use Standard_Predictor_Convolutions;

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    sol,radsol : Standard_Complex_Vectors.Vector(1..prd.dim);
    res,absres : Standard_Complex_Vectors.Vector(hom.crc'range);
    nbrit,nbfail : integer32;

  begin
    Newton_Fabry(standard_output,hom,prd,maxit,tol,nbrit,absdx,rcond,fail,
                 z,r,err,output);
    pole_step := beta1*r;
    if verbose then
      Newton_Fabry_Report(standard_output,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    Hesse_Pade(standard_output,hom,prd,svh,sol,res,beta2,
               eta,nrm,curv_step,verbose);
    Three_Way_Minima.Minimum
      (pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Predictor_Feedback(standard_output,hom,abh,prd.numcff,prd.dencff,step,
      alpha,sol,radsol,res,absres,nrm,mixres,nbfail,verbose);
  end Standard_SVD_Prediction;

  procedure DoblDobl_LU_Prediction
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep : in double_float;
		fail : out boolean; nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double double precision.

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor on the pole radius;
  --   beta2    multiplication factor on the curvature step;
  --   maxstep  the maximum step size;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations;
  --   verbose  flag for intermediate numerical output.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the curve step was minimal;
  --   nbmaxm   updated number of times the maximum step size was minimal.

    use DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Predictor_Convolutions;

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,step,mixres : double_double;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    sol,radsol : DoblDobl_Complex_Vectors.Vector(1..prd.dim);
    res,absres : DoblDobl_Complex_Vectors.Vector(hom.crc'range);
    nbrit,nbfail : integer32;
    dd_maxstep : constant double_double := create(maxstep);

  begin
    Newton_Fabry(standard_output,hom,prd,maxit,tol,nbrit,absdx,fail,
                 z,r,err,output);
    pole_step := beta1*r;
    if verbose then
      Newton_Fabry_Report(standard_output,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    Hesse_Pade(standard_output,hom,prd,svh,sol,res,beta2,
               eta,nrm,curv_step,verbose);
    Three_Way_Minima.Minimum
      (pole_step,curv_step,dd_maxstep,step,nbpole,nbhess,nbmaxm);
    Predictor_Feedback(standard_output,hom,abh,prd.numcff,prd.dencff,step,
      alpha,sol,radsol,res,absres,nrm,mixres,nbfail,verbose);
  end DoblDobl_LU_Prediction;

  procedure DoblDobl_SVD_Prediction
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep : in double_float;
		fail : out boolean; nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate extra output during computations;
  --   verbose  flag for intermediate numerical output.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the curve step was minimal;
  --   nbmaxm   updated number of times the maximum step size was minimal.

    use DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Predictor_Convolutions;

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,step,mixres : double_double;
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    sol,radsol : DoblDobl_Complex_Vectors.Vector(1..prd.dim);
    res,absres : DoblDobl_Complex_Vectors.Vector(hom.crc'range);
    nbrit,nbfail : integer32;
    dd_maxstep : constant double_double := create(maxstep);

  begin
    Newton_Fabry(standard_output,hom,prd,maxit,tol,nbrit,absdx,rcond,fail,
                 z,r,err,output);
    pole_step := beta1*r;
    if verbose then
      Newton_Fabry_Report(standard_output,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    Hesse_Pade(standard_output,hom,prd,svh,sol,res,beta2,
               eta,nrm,curv_step,verbose);
    Three_Way_Minima.Minimum
      (pole_step,curv_step,dd_maxstep,step,nbpole,nbhess,nbmaxm);
    Predictor_Feedback(standard_output,hom,abh,prd.numcff,prd.dencff,step,
      alpha,sol,radsol,res,absres,nrm,mixres,nbfail,verbose);
  end DoblDobl_SVD_Prediction;

  procedure QuadDobl_LU_Prediction
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep : in double_float;
		fail : out boolean; nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in quad double precision.

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate data output during computations;
  --   verbose  flag for intermediate numerical output.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the curve step was minimal;
  --   nbmaxm   updated number of times the maximum step size was minimal.

    use QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Predictor_Convolutions;

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,step,mixres : quad_double;
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    sol,radsol : QuadDobl_Complex_Vectors.Vector(1..prd.dim);
    res,absres : QuadDobl_Complex_Vectors.Vector(hom.crc'range);
    nbrit,nbfail : integer32;
    qd_maxstep : constant quad_double := create(maxstep);

  begin
    Newton_Fabry(standard_output,hom,prd,maxit,tol,nbrit,absdx,fail,
                 z,r,err,output);
    pole_step := beta1*r;
    if verbose then
      Newton_Fabry_Report(standard_output,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    Hesse_Pade(standard_output,hom,prd,svh,sol,res,beta2,
               eta,nrm,curv_step,verbose);
    Three_Way_Minima.Minimum
      (pole_step,curv_step,qd_maxstep,step,nbpole,nbhess,nbmaxm);
    Predictor_Feedback(standard_output,hom,abh,prd.numcff,prd.dencff,step,
      alpha,sol,radsol,res,absres,nrm,mixres,nbfail,verbose);
  end QuadDobl_LU_Prediction;

  procedure QuadDobl_SVD_Prediction
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep : in double_float;
                fail : out boolean; nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in quad double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   abh      circuits with radii as coefficients, for mixed residuals;
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   beta2    multiplication factor for the curvature step;
  --   maxstep  the maximum step size;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the curve step was minimal;
  --   nbmaxm   number of times the maximum step size was minimal;
  --   output   flag to indicate extra output during computations;
  --   verbose  flag for intermediate numerical output.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the curve step was minimal;
  --   nbmaxm   updated number of times the maximum step size was minimal.

    use QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Predictor_Convolutions;

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,step,mixres : quad_double;
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    sol,radsol : QuadDobl_Complex_Vectors.Vector(1..prd.dim);
    res,absres : QuadDobl_Complex_Vectors.Vector(hom.crc'range);
    nbrit,nbfail : integer32;
    qd_maxstep : constant quad_double := create(maxstep);

  begin
    Newton_Fabry(standard_output,hom,prd,maxit,tol,nbrit,absdx,rcond,fail,
      z,r,err,output);
    pole_step := beta1*r;
    if verbose then
      Newton_Fabry_Report(standard_output,nbrit,absdx,fail,z,r,err,
        pole_step,prd.numcff,prd.dencff,output);
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    Hesse_Pade(standard_output,hom,prd,svh,sol,res,beta2,
               eta,nrm,curv_step,verbose);
    Three_Way_Minima.Minimum
      (pole_step,curv_step,qd_maxstep,step,nbpole,nbhess,nbmaxm);
    Predictor_Feedback(standard_output,hom,abh,prd.numcff,prd.dencff,step,
      alpha,sol,radsol,res,absres,nrm,mixres,nbfail,verbose);
  end QuadDobl_SVD_Prediction;

  procedure Standard_Run_Prediction
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List;
                deg,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy chom,
  --   starting at the solutions in sols, in double precision.

    use Standard_Complex_Solutions;
    use Standard_Predictor_Convolutions;
  
    neq : constant integer32 := chom.crc'last;
    ls : Link_to_Solution := Head_Of(sols);
    dim : constant integer32 := ls.v'last;
    tmp : Solution_List := sols;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-12;
    fail,otp,usesvd : boolean;
    ans : character;
    alpha : constant double_float := 1.0E-3;
    beta1 : constant double_float := 5.0E-1;
    beta2 : constant double_float := 5.0E-3;
    maxstep : constant double_float := 1.0E-1;
    prd : Predictor;
    hss : SVD_Hessians(dim,dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    nbpole,nbhess,nbmaxm : natural32 := 0;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    if usesvd
     then prd := Create(ls.v,neq,deg,numdeg,dendeg,SVD);
     else prd := Create(ls.v,neq,deg,numdeg,dendeg,LU);
    end if;
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      Set_Lead_Coefficients(prd,ls.v);
      hss.vals := (hss.vals'range => Standard_Complex_Numbers.Create(0.0));
      if usesvd then
        Standard_SVD_Prediction(chom,abh,prd.svdata,svh,maxit,tol,alpha,
          beta1,beta2,maxstep,fail,nbpole,nbhess,nbmaxm,otp,true);
      else
        Standard_LU_Prediction(chom,abh,prd.ludata,svh,maxit,tol,alpha,
          beta1,beta2,maxstep,fail,nbpole,nbhess,nbmaxm,otp,true);
      end if;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    Clear(prd); Clear(svh);
  end Standard_Run_Prediction;

  procedure DoblDobl_Run_Prediction
              ( chom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                deg,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy chom,
  --   starting at the solutions in sols, in double double precision.

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Predictor_Convolutions;
  
    neq : constant integer32 := chom.crc'last;
    ls : Link_to_Solution := Head_Of(sols);
    dim : constant integer32 := ls.v'last;
    tmp : Solution_List := sols;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-24;
    fail,otp,usesvd : boolean;
    ans : character;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    alpha : constant double_float := 1.0E-3;
    beta1 : constant double_float := 5.0E-1;
    beta2 : constant double_float := 5.0E-3;
    maxstep : constant double_float := 1.0E-1;
    prd : Predictor;
    hss : SVD_Hessians(dim,dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    nbpole,nbhess,nbmaxm : natural32 := 0;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    if usesvd
     then prd := Create(ls.v,neq,deg,numdeg,dendeg,SVD);
     else prd := Create(ls.v,neq,deg,numdeg,dendeg,LU);
    end if;
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      Set_Lead_Coefficients(prd,ls.v);
      hss.vals := (hss.vals'range => zero);
      if usesvd then
        DoblDobl_SVD_Prediction(chom,abh,prd.svdata,svh,maxit,tol,alpha,
          beta1,beta2,maxstep,fail,nbpole,nbhess,nbmaxm,otp,true);
      else
        DoblDobl_LU_Prediction(chom,abh,prd.ludata,svh,maxit,tol,alpha,
          beta1,beta2,maxstep,fail,nbpole,nbhess,nbmaxm,otp,true);
      end if;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    Clear(prd); Clear(svh);
  end DoblDobl_Run_Prediction;

  procedure QuadDobl_Run_Prediction
              ( chom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                deg,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy chom,
  --   starting at the solutions in sols, in quad double precision.

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Predictor_Convolutions;
  
    neq : constant integer32 := chom.crc'last;
    ls : Link_to_Solution := Head_Of(sols);
    dim : constant integer32 := ls.v'last;
    tmp : Solution_List := sols;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-48;
    fail,otp,usesvd : boolean;
    ans : character;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    alpha : constant double_float := 1.0E-3;
    beta1 : constant double_float := 5.0E-1;
    beta2 : constant double_float := 5.0E-3;
    maxstep : constant double_float := 1.0E-1;
    prd : Predictor;
    hss : SVD_Hessians(dim,dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    nbpole,nbhess,nbmaxm : natural32 := 0;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    if usesvd
     then prd := Create(ls.v,neq,deg,numdeg,dendeg,SVD);
     else prd := Create(ls.v,neq,deg,numdeg,dendeg,LU);
    end if;
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      Set_Lead_Coefficients(prd,ls.v);
      hss.vals := (hss.vals'range => zero);
      if usesvd then
        QuadDobl_SVD_Prediction(chom,abh,prd.svdata,svh,maxit,tol,alpha,
          beta1,beta2,maxstep,fail,nbpole,nbhess,nbmaxm,otp,true);
      else
        QuadDobl_LU_Prediction(chom,abh,prd.ludata,svh,maxit,tol,alpha,
          beta1,beta2,maxstep,fail,nbpole,nbhess,nbmaxm,otp,true);
      end if;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    Clear(prd); Clear(svh);
  end QuadDobl_Run_Prediction;

  procedure Standard_Test_Prediction
              ( nq,idxpar,numdeg,dendeg : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The Standard_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    use Standard_Complex_Solutions;

    hom : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
        := Standard_Homotopy.Homotopy_System;
    serhom : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idxpar);
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    deg : constant integer32 := numdeg + dendeg + 2;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);
    ans : character;

    use Standard_Complex_Vectors;

  begin
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    cnvhom := Make_Convolution_System(serhom,natural32(deg));
    abshom := Residual_Convolution_System(cnvhom);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    new_line;
    put("Check all start solutions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Checking the start solutions ...");
      for k in 1..Length_Of(sols) loop
        ls := Head_Of(tmp);
        declare
          y : constant Standard_Complex_Vectors.Vector
            := Standard_Speelpenning_Convolutions.Eval(cnvhom.crc,ls.v,zero);
        begin
          put("Value at solution "); put(k,1); put_line(" :");
          put_line(y);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    Standard_Run_Prediction(cnvhom,abshom,sols,deg,numdeg,dendeg);
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction
              ( nq,idxpar,numdeg,dendeg : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The DoblDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    use DoblDobl_Complex_Solutions;

    hom : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
        := DoblDobl_Homotopy.Homotopy_System;
    serhom : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idxpar);
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    deg : constant integer32 := numdeg + dendeg + 2;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    ans : character;

  begin
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    cnvhom := Make_Convolution_System(serhom,natural32(deg));
    abshom := Residual_Convolution_System(cnvhom);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    new_line;
    put("Check all start solutions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Checking the start solutions ...");
      for k in 1..Length_Of(sols) loop
        ls := Head_Of(tmp);
        declare
          y : constant DoblDobl_Complex_Vectors.Vector
            := DoblDobl_Speelpenning_Convolutions.Eval(cnvhom.crc,ls.v,zero);
        begin
          put("Value at solution "); put(k,1); put_line(" :");
          put_line(y);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    DoblDobl_Run_Prediction(cnvhom,abshom,sols,deg,numdeg,dendeg);
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction
              ( nq,idxpar,numdeg,dendeg : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The QuadDobl_Homotopy is initialized with nq equations
  --   and sols contains the solutions of the start system.
  --   The parameter idxpar is the index to the continuation parameter.

    use QuadDobl_Complex_Solutions;

    hom : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
        := QuadDobl_Homotopy.Homotopy_System;
    serhom : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
           := Series_and_Homotopies.Create(hom,idxpar);
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    deg : constant integer32 := numdeg + dendeg + 2;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    ans : character;

  begin
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    cnvhom := Make_Convolution_System(serhom,natural32(deg));
    abshom := Residual_Convolution_System(cnvhom);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    new_line;
    put("Check all start solutions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Checking the start solutions ...");
      for k in 1..Length_Of(sols) loop
        ls := Head_Of(tmp);
        declare
          y : constant QuadDobl_Complex_Vectors.Vector
            := QuadDobl_Speelpenning_Convolutions.Eval(cnvhom.crc,ls.v,zero);
        begin
          put("Value at solution "); put(k,1); put_line(" :");
          put_line(y);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    QuadDobl_Run_Prediction(cnvhom,abshom,sols,deg,numdeg,dendeg);
  end QuadDobl_Test_Prediction;

  procedure Standard_Main ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    nbeq,idxpar : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.Standard_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      Standard_Test_Prediction(nbeq,nbeq+1,numdeg,dendeg,sols);
    else
      declare
        dropsols : constant Standard_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        Standard_Test_Prediction(nbeq,idxpar,numdeg,dendeg,dropsols);
      end;
    end if;
  end Standard_Main;

  procedure DoblDobl_Main ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq,idxpar : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.DoblDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      DoblDobl_Test_Prediction(nbeq,nbeq+1,numdeg,dendeg,sols);
    else
      declare
        dropsols : constant DoblDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        DoblDobl_Test_Prediction(nbeq,idxpar,numdeg,dendeg,dropsols);
      end;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq,idxpar : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.QuadDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      QuadDobl_Test_Prediction(nbeq,nbeq+1,numdeg,dendeg,sols);
    else
      declare
        dropsols : constant QuadDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        QuadDobl_Test_Prediction(nbeq,idxpar,numdeg,dendeg,dropsols);
      end;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision,
  --   the degrees of the Pade approximants, and then launches
  --   the test corresponding to the selected precision.

    precision : character;
    numdeg,dendeg : integer32 := 0;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    new_line;
    put("Give the degree of the Pade numerator : "); get(numdeg);
    put("Give the degree of the Pade denominator : "); get(dendeg);
    case precision is
      when '0' => Standard_Main(numdeg,dendeg);
      when '1' => DoblDobl_Main(numdeg,dendeg);
      when '2' => QuadDobl_Main(numdeg,dendeg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_padepcnv;
