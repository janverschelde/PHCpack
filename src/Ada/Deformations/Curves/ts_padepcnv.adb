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
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
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
with Newton_Power_Convolutions;          use Newton_Power_Convolutions;
with Convergence_Radius_Estimates;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;

procedure ts_padepcnv is

-- DESCRIPTION :
--   Development of the Pade predictor on convolution circuits.

  procedure Standard_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Standard_Predictor_Convolutions.Link_to_Predictor;
                maxit : in integer32; tol : in double_float;
                usesvd : in boolean;
                fail : out boolean; output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   usesvd   flag to indicate if the SVD has to be used;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   fail     indicates failure status.

    use Standard_Rational_Approximations;
    use Standard_Speelpenning_Convolutions;
    use Standard_Predictor_Convolutions;

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond : double_float;
    eva : Standard_Complex_Vectors.Vector(1..prd.dim);
    res : Standard_Complex_Vectors.Vector(hom.crc'range);
    info,nbrit : integer32;
    dim : constant integer32 := prd.dim;
    deg : constant integer32 := prd.deg;
    neq : constant integer32 := hom.crc'last;
    dx : Standard_Complex_VecVecs.VecVec(1..dim);
    xd : Standard_Complex_VecVecs.VecVec(0..deg);
    svl : Standard_Complex_Vectors.Vector(1..dim+1);
    U : Standard_Complex_Matrices.Matrix(1..neq,1..neq);
    V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    ewrk : constant Standard_Complex_Vectors.Link_to_Vector
         := new Standard_Complex_Vectors.Vector(1..dim);

  begin
    if output then
      if usesvd then
        dx := Standard_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
        xd := Standard_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
        SVD_Newton_Steps
          (standard_output,hom,prd.sol,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,prd.wrk,false);
      else
        LU_Newton_Steps
          (standard_output,hom,prd.sol,maxit,nbrit,tol,absdx,fail,
           info,prd.newtpiv,prd.wrk,false);
      end if;
      Convergence_Radius_Estimates.Fabry(prd.sol,z,r,err,fail);
    else
      if usesvd then
        dx := Standard_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
        xd := Standard_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
        SVD_Newton_Steps
          (hom,prd.sol,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,prd.wrk,false,false);
      else
        LU_Newton_Steps
          (hom,prd.sol,maxit,nbrit,tol,absdx,fail,
           info,prd.newtpiv,prd.wrk,false,false);
      end if;
      Convergence_Radius_Estimates.Fabry(prd.sol,z,r,err,fail,false);
    end if;
    if not fail then
      put("z : "); put(z); 
      put("  error estimate :"); put(err,3); new_line;
      put("estimated radius :"); put(r,3); new_line;
    end if;
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
    Evaluate(prd.numcff,prd.dencff,r/2.0,eva);
    z := Standard_Complex_Numbers.Create(r/2.0);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : ");
    put_line(res);
  end Standard_Prediction;

  procedure DoblDobl_Prediction
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in DoblDobl_Predictor_Convolutions.Link_to_Predictor;
                maxit : in integer32; tol : in double_float;
                usesvd : in boolean;
                fail : out boolean; output : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double double precision.

  -- DESCRIPTION :
  --   Runs Newton's method on the power series in prd,
  --   applies Fabry's theorem and constructs Pade approximants
  --   to predict the next solution, in double precision.

  -- ON ENTRY :
  --   hom      homotopy convolution circuit system
  --   prd      predictor data for LU Newton and Pade approximants;
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   usesvd   flag to indicate if the SVD has to be used;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   fail     indicates failure status.

    use DoblDobl_Rational_Approximations;
    use DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Predictor_Convolutions;

    z : DoblDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond : double_double;
    eva : DoblDobl_Complex_Vectors.Vector(1..prd.dim);
    res : DoblDobl_Complex_Vectors.Vector(hom.crc'range);
    info,nbrit : integer32;
    dim : constant integer32 := prd.dim;
    deg : constant integer32 := prd.deg;
    neq : constant integer32 := hom.crc'last;
    dx : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    xd : DoblDobl_Complex_VecVecs.VecVec(0..deg);
    svl : DoblDobl_Complex_Vectors.Vector(1..dim+1);
    U : DoblDobl_Complex_Matrices.Matrix(1..neq,1..neq);
    V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    ewrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := new DoblDobl_Complex_Vectors.Vector(1..dim);

  begin
    if output then
      if usesvd then
        dx := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
        xd := DoblDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
        SVD_Newton_Steps
          (standard_output,hom,prd.sol,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,prd.wrk,false);
      else
        LU_Newton_Steps
          (standard_output,hom,prd.sol,maxit,nbrit,tol,absdx,fail,
           info,prd.newtpiv,prd.wrk,false);
      end if;
      Convergence_Radius_Estimates.Fabry(prd.sol,z,r,err,fail);
    else
      if usesvd then
        dx := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
        xd := DoblDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
        SVD_Newton_Steps
          (hom,prd.sol,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,prd.wrk,false,false);
      else
        LU_Newton_Steps
          (hom,prd.sol,maxit,nbrit,tol,absdx,fail,
           info,prd.newtpiv,prd.wrk,false,false);
      end if;
      Convergence_Radius_Estimates.Fabry(prd.sol,z,r,err,fail,false);
    end if;
    if not fail then
      put("z : "); put(z); 
      put("  error estimate :"); put(err,3); new_line;
      put("estimated radius :"); put(r,3); new_line;
    end if;
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
    Evaluate(prd.numcff,prd.dencff,r/2.0,eva);
    z := DoblDobl_Complex_Numbers.Create(r/2.0);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : ");
    put_line(res);
  end DoblDobl_Prediction;

  procedure QuadDobl_Prediction
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                prd : in QuadDobl_Predictor_Convolutions.Link_to_Predictor;
                maxit : in integer32; tol : in double_float;
                usesvd : in boolean;
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
  --   maxit    maximum number of iterations in Newton's method;
  --   tol      tolerance on the correction term;
  --   usesvd   flag to indicate if the SVD has to be used;
  --   output   flag to indicate extra output during computations.

  -- ON RETURN :
  --   prd      contains solution series and Pade approximants;
  --   fail     indicates failure status.

    use QuadDobl_Rational_Approximations;
    use QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Predictor_Convolutions;

    z : QuadDobl_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond : quad_double;
    eva : QuadDobl_Complex_Vectors.Vector(1..prd.dim);
    res : QuadDobl_Complex_Vectors.Vector(hom.crc'range);
    info,nbrit : integer32;
    dim : constant integer32 := prd.dim;
    deg : constant integer32 := prd.deg;
    neq : constant integer32 := hom.crc'last;
    dx : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    xd : QuadDobl_Complex_VecVecs.VecVec(0..deg);
    svl : QuadDobl_Complex_Vectors.Vector(1..dim+1);
    U : QuadDobl_Complex_Matrices.Matrix(1..neq,1..neq);
    V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    ewrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
         := new QuadDobl_Complex_Vectors.Vector(1..dim);

  begin
    if output then
      if usesvd then
        dx := QuadDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
        xd := QuadDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
        SVD_Newton_Steps
          (standard_output,hom,prd.sol,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,prd.wrk,false);
      else
        LU_Newton_Steps
          (standard_output,hom,prd.sol,maxit,nbrit,tol,absdx,fail,
           info,prd.newtpiv,prd.wrk,false);
      end if;
      Convergence_Radius_Estimates.Fabry(prd.sol,z,r,err,fail);
    else
      if usesvd then
        dx := QuadDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
        xd := QuadDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
        SVD_Newton_Steps
          (hom,prd.sol,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,prd.wrk,false,false);
      else
        LU_Newton_Steps
          (hom,prd.sol,maxit,nbrit,tol,absdx,fail,
           info,prd.newtpiv,prd.wrk,false,false);
      end if;
      Convergence_Radius_Estimates.Fabry(prd.sol,z,r,err,fail,false);
    end if;
    if not fail then
      put("z : "); put(z); 
      put("  error estimate : "); put(err,3); new_line;
      put("estimated radius : "); put(r,3); new_line;
    end if;
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
    Evaluate(prd.numcff,prd.dencff,r/2.0,eva);
    z := QuadDobl_Complex_Numbers.Create(r/2.0);
    res := Eval(hom.crc,eva,z);
    put_line("Evaluation of the predicted solution : ");
    put_line(res);
  end QuadDobl_Prediction;

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
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-12;
    fail,otp,usesvd : boolean;
    ans : character;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        prd : Link_to_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
      begin
        Standard_Prediction(chom,prd,maxit,tol,usesvd,fail,otp);
        Clear(prd);
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
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-24;
    fail,otp,usesvd : boolean;
    ans : character;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        prd : Link_to_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
      begin
        DoblDobl_Prediction(chom,prd,maxit,tol,usesvd,fail,otp);
        Clear(prd);
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
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    maxit : integer32 := 0; 
    tol : constant double_float := 1.0E-48;
    fail,otp,usesvd : boolean;
    ans : character;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    put("Output during Newton's method ? (y/n) "); Ask_Yes_or_No(ans);
    otp := (ans = 'y');
    put("Use SVD ? (y/n) "); Ask_Yes_or_No(ans); usesvd := (ans = 'y');
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        prd : Link_to_Predictor := Create(ls.v,neq,deg,numdeg,dendeg);
      begin
        QuadDobl_Prediction(chom,prd,maxit,tol,usesvd,fail,otp);
        Clear(prd);
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
