with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;
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
with Jacobian_Convolution_Circuits;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;

procedure ts_padepcnv is

-- DESCRIPTION :
--   Development of the Pade predictor on convolution circuits.

  procedure Standard_LU_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Standard_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1 : in double_float;
                fail : out boolean; output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with LU,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status.

    use Standard_Complex_Singular_Values;
    use Standard_Rational_Approximations;
    use Standard_Speelpenning_Convolutions;
    use Standard_Predictor_Convolutions;

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,eta,pole_step,nrm : double_float;
    eva : Standard_Complex_Vectors.Vector(1..prd.dim);
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    sol : Standard_Complex_Vectors.Vector(1..prd.dim);
    res : Standard_Complex_Vectors.Vector(hom.crc'range);
    info,nbrit : integer32;

  begin
    Predict(hom,prd,maxit,tol,nbrit,absdx,fail,z,r,err,output);
    put("#iterations : "); put(nbrit,1);
    put("  |dx| :"); put(absdx,3); new_line;
    if fail then
      put_line("Predictor failed!");
    else
      put("z : "); put(z); 
      put("  error estimate :"); put(err,3); new_line;
      put("estimated radius :"); put(r,3); new_line;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
   -- with LU, the system should be square, so the svh work space works
    svh.H := Jacobian_Convolution_Circuits.Jacobian(hom.crc,sol);
    SVD(svh.H,svh.dim,svh.dim,svh.svl,svh.ewrk,svh.U,svh.V,11,info);
    svh.vals(0) := svh.svl(svh.dim);
    Second(hom,svh,sol);
    put_line("All singular values : "); put_line(svh.vals);
    eta := Standard_Predictor_Convolutions.Standard_Distance(svh);
    put("  eta :"); put(eta,3); new_line;
    pole_step := beta1*r;
    Evaluate(prd.numcff,prd.dencff,pole_step,eva);
    z := Standard_Complex_Numbers.Create(pole_step);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : "); put_line(res);
    nrm := Standard_Complex_Vector_Norms.Max_Norm(res);
    put("The predictor residual :"); put(nrm,3);
    if nrm < alpha
     then put_line("  okay");
     else put(" >"); put(alpha,3); new_line;
    end if;
  end Standard_LU_Prediction;

  procedure Standard_SVD_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Standard_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1 : in double_float;
                fail : out boolean; output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status.

    use Standard_Rational_Approximations;
    use Standard_Speelpenning_Convolutions;
    use Standard_Predictor_Convolutions;

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,eta,pole_step,nrm : double_float;
    eva : Standard_Complex_Vectors.Vector(1..prd.dim);
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    sol : Standard_Complex_Vectors.Vector(1..prd.dim);
    res : Standard_Complex_Vectors.Vector(hom.crc'range);
    nbrit : integer32;

  begin
    Predict(hom,prd,maxit,tol,nbrit,absdx,rcond,fail,z,r,err,output);
    put("#iterations : "); put(nbrit,1);
    put("  |dx| :"); put(absdx,3);
    put("  rcond :"); put(rcond,3); new_line;
    put_line("The singular values : "); put_line(prd.svl);
    if fail then
      put_line("Predictor failed!");
    else
      put("z : "); put(z); 
      put("  error estimate :"); put(err,3); new_line;
      put("estimated radius :"); put(r,3); new_line;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    svh.vals(0) := prd.svl(prd.dim);
    Second(hom,svh,sol);
    put_line("All singular values : "); put_line(svh.vals);
    eta := Standard_Predictor_Convolutions.Standard_Distance(svh);
    put("  eta :"); put(eta,3); new_line;
    pole_step := beta1*r;
    Evaluate(prd.numcff,prd.dencff,pole_step,eva);
    z := Standard_Complex_Numbers.Create(pole_step);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : "); put_line(res);
    nrm := Standard_Complex_Vector_Norms.Max_Norm(res);
    put("The predictor residual :"); put(nrm,3);
    if nrm < alpha
     then put_line("  okay");
     else put(" > "); put(alpha,3); new_line;
    end if;
  end Standard_SVD_Prediction;

  procedure DoblDobl_LU_Prediction
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1 : in double_float;
                fail : out boolean; output : in boolean ) is

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
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor on the pole step;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status.

    use DoblDobl_Complex_Singular_Values;
    use DoblDobl_Rational_Approximations;
    use DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Predictor_Convolutions;

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,eta,pole_step,nrm : double_double;
    eva : DoblDobl_Complex_Vectors.Vector(1..prd.dim);
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    sol : DoblDobl_Complex_Vectors.Vector(1..prd.dim);
    res : DoblDobl_Complex_Vectors.Vector(hom.crc'range);
    info,nbrit : integer32;

  begin
    Predict(hom,prd,maxit,tol,nbrit,absdx,fail,z,r,err,output);
    put("#iterations : "); put(nbrit,1);
    put("  |dx| : "); put(absdx,3); new_line;
    if fail then
      put_line("Predictor failed!");
    else
      put("z : "); put(z); 
      put("  error estimate : "); put(err,3); new_line;
      put("estimated radius : "); put(r,3); new_line;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
   -- with LU, the system should be square, so the svh work space works
    svh.H := Jacobian_Convolution_Circuits.Jacobian(hom.crc,sol);
    SVD(svh.H,svh.dim,svh.dim,svh.svl,svh.ewrk,svh.U,svh.V,11,info);
    svh.vals(0) := svh.svl(svh.dim);
    Second(hom,svh,sol);
    put_line("All singular values : "); put_line(svh.vals);
    eta := DoblDobl_Predictor_Convolutions.DoblDobl_Distance(svh);
    put("  eta : "); put(eta,3); new_line;
    pole_step := beta1*r;
    Evaluate(prd.numcff,prd.dencff,pole_step,eva);
    z := DoblDobl_Complex_Numbers.Create(pole_step);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : "); put_line(res);
    nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(res);
    put("The predictor residual :"); put(nrm,3);
    if nrm < alpha
     then put_line("  okay");
     else put(" > "); put(alpha,3); new_line;
    end if;
  end DoblDobl_LU_Prediction;

  procedure DoblDobl_SVD_Prediction
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1 : in double_float;
                fail : out boolean; output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status.

    use DoblDobl_Rational_Approximations;
    use DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Predictor_Convolutions;

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,eta,pole_step,nrm : double_double;
    eva : DoblDobl_Complex_Vectors.Vector(1..prd.dim);
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    sol : DoblDobl_Complex_Vectors.Vector(1..prd.dim);
    res : DoblDobl_Complex_Vectors.Vector(hom.crc'range);
    nbrit : integer32;

  begin
    Predict(hom,prd,maxit,tol,nbrit,absdx,rcond,fail,z,r,err,output);
    put("#iterations : "); put(nbrit,1);
    put("  |dx| : "); put(absdx,3);
    put("  rcond : "); put(rcond,3); new_line;
    put_line("The singular values : "); put_line(prd.svl);
    if fail then
      put_line("Predictor failed!");
    else
      put("z : "); put(z); 
      put("  error estimate : "); put(err,3); new_line;
      put("estimated radius : "); put(r,3); new_line;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    svh.vals(0) := prd.svl(prd.dim);
    Second(hom,svh,sol);
    put_line("All singular values : "); put_line(svh.vals);
    eta := DoblDobl_Predictor_Convolutions.DoblDobl_Distance(svh);
    put("  eta : "); put(eta,3); new_line;
    pole_step := beta1*r;
    Evaluate(prd.numcff,prd.dencff,pole_step,eva);
    z := DoblDobl_Complex_Numbers.Create(pole_step);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : "); put_line(res);
    nrm := DoblDobl_Complex_Vector_Norms.Max_Norm(res);
    put("The predictor residual :"); put(nrm,3);
    if nrm < alpha
     then put_line("  okay");
     else put(" > "); put(alpha,3); new_line;
    end if;
  end DoblDobl_SVD_Prediction;

  procedure QuadDobl_LU_Prediction
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_LU_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1 : in double_float;
                fail : out boolean; output : in boolean ) is

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
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status.

    use QuadDobl_Complex_Singular_Values;
    use QuadDobl_Rational_Approximations;
    use QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Predictor_Convolutions;

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,eta,pole_step,nrm : quad_double;
    eva : QuadDobl_Complex_Vectors.Vector(1..prd.dim);
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    sol : QuadDobl_Complex_Vectors.Vector(1..prd.dim);
    res : QuadDobl_Complex_Vectors.Vector(hom.crc'range);
    info,nbrit : integer32;

  begin
    Predict(hom,prd,maxit,tol,nbrit,absdx,fail,z,r,err,output);
    put("#iterations : "); put(nbrit,1);
    put("  |dx| : "); put(absdx,3); new_line;
    if fail then
      put_line("Predictor failed!");
    else
      put("z : "); put(z); 
      put("  error estimate : "); put(err,3); new_line;
      put("estimated radius : "); put(r,3); new_line;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
   -- with LU, the system should be square, so the svh work space works
    svh.H := Jacobian_Convolution_Circuits.Jacobian(hom.crc,sol);
    SVD(svh.H,svh.dim,svh.dim,svh.svl,svh.ewrk,svh.U,svh.V,11,info);
    svh.vals(0) := svh.svl(svh.dim);
    Second(hom,svh,sol);
    put_line("All singular values : "); put_line(svh.vals);
    eta := QuadDobl_Predictor_Convolutions.QuadDobl_Distance(svh);
    put("  eta : "); put(eta,3); new_line;
    pole_step := beta1*r;
    Evaluate(prd.numcff,prd.dencff,pole_step,eva);
    z := QuadDobl_Complex_Numbers.Create(pole_Step);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : "); put_line(res);
    nrm := QuadDobl_Complex_Vector_Norms.Max_Norm(res);
    put("The predictor residual :"); put(nrm,3);
    if nrm < alpha
     then put_line("  okay");
     else put(" > "); put(alpha,3); new_line;
    end if;
  end QuadDobl_LU_Prediction;

  procedure QuadDobl_SVD_Prediction
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Predictor;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                maxit : in integer32; tol : in double_float;
                alpha,beta1 : in double_float;
                fail : out boolean; output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd with SVD,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in quad double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   svh      data for the curvature estimation;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   alpha    tolerance on the predictor residual;
  --   beta1    multiplication factor for the pole radius;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   svh      contains largest singular values of all Hessians;
  --   fail     indicates failure status.

    use QuadDobl_Rational_Approximations;
    use QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Predictor_Convolutions;

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,eta,pole_step,nrm : quad_double;
    eva : QuadDobl_Complex_Vectors.Vector(1..prd.dim);
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    sol : QuadDobl_Complex_Vectors.Vector(1..prd.dim);
    res : QuadDobl_Complex_Vectors.Vector(hom.crc'range);
    nbrit : integer32;

  begin
    Predict(hom,prd,maxit,tol,nbrit,absdx,rcond,fail,z,r,err,output);
    put("#iterations : "); put(nbrit,1);
    put("  |dx| : "); put(absdx,3);
    put("  rcond : "); put(rcond,3); new_line;
    put_line("The singular values : "); put_line(prd.svl);
    if fail then
      put_line("Predictor failed!");
    else
      put("z : "); put(z); 
      put("  error estimate : "); put(err,3); new_line;
      put("estimated radius : "); put(r,3); new_line;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); sol(k) := lnk(0);
    end loop;
    svh.vals(0) := prd.svl(prd.dim);
    Second(hom,svh,sol);
    put_line("All singular values : "); put_line(svh.vals);
    eta := QuadDobl_Predictor_Convolutions.QuadDobl_Distance(svh);
    put("  eta : "); put(eta,3); new_line;
    pole_step := beta1*r;
    Evaluate(prd.numcff,prd.dencff,pole_step,eva);
    z := QuadDobl_Complex_Numbers.Create(pole_step);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : "); put_line(res);
    nrm := QuadDobl_Complex_Vector_Norms.Max_Norm(res);
    put("The predictor residual :"); put(nrm,3);
    if nrm < alpha
     then put_line("  okay");
     else put(" > "); put(alpha,3); new_line;
    end if;
  end QuadDobl_SVD_Prediction;

  procedure Standard_Run_Prediction
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List;
                deg,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy chom,
  --   starting at the solutions in sols, in double precision.

    use Standard_Complex_Solutions;
    use Standard_Predictor_Convolutions;
  
    neq : constant integer32 := chom.crc'last;
    dim : integer32;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-12;
    fail,otp,usesvd : boolean;
    ans : character;
    alpha : constant double_float := 1.0E-3;
    beta1 : constant double_float := 5.0E-1;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp); dim := ls.v'last;
      declare
        hss : SVD_Hessians(dim,dim+1);
        svh : Link_to_SVD_Hessians;
      begin
        hss.vals := (hss.vals'range => Standard_Complex_Numbers.Create(0.0));
        svh := new SVD_Hessians'(hss);
        if usesvd then
          declare
            prd : Link_to_SVD_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
          begin
            Standard_SVD_Prediction
              (chom,prd,svh,maxit,tol,alpha,beta1,fail,otp);
            Clear(prd);
          end;
        else
          declare
            prd : Link_to_LU_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
          begin
            Standard_LU_Prediction
              (chom,prd,svh,maxit,tol,alpha,beta1,fail,otp);
            Clear(prd);
          end;
        end if;
        Clear(svh);
      end;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Run_Prediction;

  procedure DoblDobl_Run_Prediction
              ( chom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                deg,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy chom,
  --   starting at the solutions in sols, in double double precision.

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Predictor_Convolutions;
  
    neq : constant integer32 := chom.crc'last;
    dim : integer32;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-24;
    fail,otp,usesvd : boolean;
    ans : character;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    alpha : constant double_float := 1.0E-3;
    beta1 : constant double_float := 5.0E-1;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp); dim := ls.v'last;
      declare
        hss : SVD_Hessians(dim,dim+1);
        svh : Link_to_SVD_Hessians;
      begin
        hss.vals := (hss.vals'range => zero);
        svh := new SVD_Hessians'(hss);
        if usesvd then
          declare
            prd : Link_to_SVD_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
          begin
            DoblDobl_SVD_Prediction
              (chom,prd,svh,maxit,tol,alpha,beta1,fail,otp);
            Clear(prd);
          end;
        else
          declare
            prd : Link_to_LU_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
          begin
            DoblDobl_LU_Prediction
              (chom,prd,svh,maxit,tol,alpha,beta1,fail,otp);
            Clear(prd);
          end;
        end if;
        Clear(svh);
      end;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end DoblDobl_Run_Prediction;

  procedure QuadDobl_Run_Prediction
              ( chom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                deg,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy chom,
  --   starting at the solutions in sols, in quad double precision.

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Predictor_Convolutions;
  
    neq : constant integer32 := chom.crc'last;
    dim : integer32;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-48;
    fail,otp,usesvd : boolean;
    ans : character;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    alpha : constant double_float := 1.0E-3;
    beta1 : constant double_float := 5.0E-1;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp); dim := ls.v'last;
      declare
        hss : SVD_Hessians(dim,dim+1);
        svh : Link_to_SVD_Hessians;
      begin
        hss.vals := (hss.vals'range => zero);
        svh := new SVD_Hessians'(hss);
        if usesvd then
          declare
            prd : Link_to_SVD_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
          begin
            QuadDobl_SVD_Prediction
              (chom,prd,svh,maxit,tol,alpha,beta1,fail,otp);
            Clear(prd);
          end;
        else
          declare
            prd : Link_to_LU_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
          begin
            QuadDobl_LU_Prediction
              (chom,prd,svh,maxit,tol,alpha,beta1,fail,otp);
            Clear(prd);
          end;
        end if;
        Clear(svh);
      end;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end QuadDobl_Run_Prediction;

  procedure Standard_Test_Prediction
              ( nq,idxpar : in integer32;
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
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    deg,numdeg,dendeg : integer32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    put("Give the degree of the series : "); get(deg);
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    cnvhom := Make_Convolution_System(serhom,natural32(deg));
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
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
    Standard_Run_Prediction(cnvhom,sols,deg,numdeg,dendeg);
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction
              ( nq,idxpar : in integer32;
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
    cnvhom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    deg,numdeg,dendeg : integer32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));

  begin
    put("Give the degree of the series : "); get(deg);
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    cnvhom := Make_Convolution_System(serhom,natural32(deg));
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
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
    DoblDobl_Run_Prediction(cnvhom,sols,deg,numdeg,dendeg);
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction
              ( nq,idxpar : in integer32;
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
    cnvhom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    deg,numdeg,dendeg : integer32 := 0;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));

  begin
    put("Give the degree of the series : "); get(deg);
    put("Give the degree of the numerator : "); get(numdeg);
    put("Give the degree of the denominator : "); get(dendeg);
    Complex_Series_and_Polynomials.Set_Degree(serhom,deg);
    cnvhom := Make_Convolution_System(serhom,natural32(deg));
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
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
    QuadDobl_Run_Prediction(cnvhom,sols,deg,numdeg,dendeg);
  end QuadDobl_Test_Prediction;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    nbeq,idxpar : integer32;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.Standard_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      Standard_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant Standard_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        Standard_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    nbeq,idxpar : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.DoblDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      DoblDobl_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant DoblDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        DoblDobl_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    nbeq,idxpar : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Test_Series_Predictors.QuadDobl_Homotopy_Reader(nbeq,idxpar,sols);
    new_line;
    if idxpar = 0 then
      QuadDobl_Test_Prediction(nbeq,nbeq+1,sols);
    else
      declare
        dropsols : constant QuadDobl_Complex_Solutions.Solution_List
                 := Solution_Drops.Drop(sols,natural32(idxpar));
      begin
        QuadDobl_Test_Prediction(nbeq,idxpar,dropsols);
      end;
    end if;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    case precision is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_padepcnv;
