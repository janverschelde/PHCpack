with unchecked_deallocation;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Vector_Splitters;
with Standard_Matrix_Splitters;
with Standard_Complex_Singular_Values;
with Standard_Inlined_Singular_Values;
with Standard_Mixed_Residuals;
with Standard_Rational_Approximations;
with Standard_Newton_Convolutions;
with Standard_Newton_Convolution_Steps;
with Staggered_Newton_Convolutions;
with Convergence_Radius_Estimates;
with Jacobian_Convolution_Circuits;
with Hessian_Convolution_Circuits;
with Homotopy_Pade_Approximants;
with Series_and_Predictors;
with Three_Way_Minima;

package body Standard_Predictor_Convolutions is

-- CONSTRUCTORS :

  function Create ( sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    inlined : boolean := true ) return LU_Predictor is

    res : LU_Predictor(sol'last,deg,numdeg,dendeg);
    zero : constant Complex_Number := Create(0.0);
    knum : constant Standard_Complex_Vectors.Vector(0..numdeg)
         := (0..numdeg => zero);
    kden : constant Standard_Complex_Vectors.Vector(0..dendeg)
         := (0..dendeg => zero);

  begin
    res.sol := Standard_Newton_Convolutions.Series_Coefficients(sol,deg);
    res.wrk := new Standard_Complex_Vectors.Vector'(1..neq => zero);
    for k in sol'range loop
      res.numcff(k) := new Standard_Complex_Vectors.Vector'(knum);
      res.dencff(k) := new Standard_Complex_Vectors.Vector'(kden);
    end loop;
    if inlined then
      Standard_Floating_VecVecVecs.Allocate(res.rv,1,deg,1,neq,1,neq);
      Standard_Floating_VecVecVecs.Allocate(res.iv,1,deg,1,neq,1,neq);
      res.rc := Standard_Vector_Splitters.Allocate(neq,neq,1,1);
      res.ic := Standard_Vector_Splitters.Allocate(neq,neq,1,1);
      res.rb := Standard_Vector_Splitters.Allocate(deg,neq,0,1);
      res.ib := Standard_Vector_Splitters.Allocate(deg,neq,0,1);
      res.ry := new Standard_Floating_Vectors.Vector'(1..neq => 0.0);
      res.iy := new Standard_Floating_Vectors.Vector'(1..neq => 0.0);
    end if;
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

    use Standard_Speelpenning_Convolutions;

  begin
    res.sol := Standard_Newton_Convolutions.Series_Coefficients(sol,deg);
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
                    neq,deg,numdeg,dendeg : integer32;
                    inlined : boolean := true ) 
                  return Link_to_LU_Predictor is

    prd : constant LU_Predictor(sol'last,deg,numdeg,dendeg)
        := Create(sol,neq,deg,numdeg,dendeg,inlined);
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

  function Create ( sol : Standard_Complex_Vectors.Vector;
	            neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type;
                    inlined : boolean := true ) return Predictor is
  begin
    case kind is
      when LU =>
        declare
          res : Predictor(LU);
        begin
          res.ludata := Create(sol,neq,deg,numdeg,dendeg,inlined);
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
                    sol : Standard_Complex_Vectors.Vector;
                    neq,deg,numdeg,dendeg : integer32;
                    kind : Predictor_Type;
                    inlined : boolean := true ) return Predictor_Array is

    res : Predictor_Array(1..nbr);

  begin
    for k in 1..nbr loop
      res(k) := Create(sol,neq,deg,numdeg,dendeg,kind,inlined);
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

  begin
    hss.vals := (hss.vals'range => Standard_Complex_Numbers.Create(0.0));
    hss.first := true;
    hss.sr := new Standard_Floating_Vectors.Vector'(1..dim+1 => 0.0);
    hss.si := new Standard_Floating_Vectors.Vector'(1..dim+1 => 0.0);
    hss.er := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    hss.ei := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    hss.wr := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    hss.wi := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    hss.xrv := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    hss.xiv := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
   -- job = 0, so no allocations for urv, uiv, vrv, and viv
   -- hss.urv := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
   -- hss.urv := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
   -- hss.viv := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
   -- hss.viv := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
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
                s : in Standard_Complex_Vectors.Vector ) is

    lnk : Standard_Complex_Vectors.Link_to_Vector;

  begin
    case p.kind is
      when LU =>
        for k in p.ludata.sol'range loop
          lnk := p.ludata.sol(k);
          lnk(0) := s(k);
          for i in 1..lnk'last loop
            lnk(i) := Standard_Complex_Numbers.Create(0.0);
          end loop;
        end loop;
      when SVD =>
        for k in p.svdata.sol'range loop
          lnk := p.svdata.sol(k);
          lnk(0) := s(k);
          for i in 1..lnk'last loop
            lnk(i) := Standard_Complex_Numbers.Create(0.0);
          end loop;
        end loop;
    end case;
  end Set_Lead_Coefficients;

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

-- NEWTON-FABRY ON COEFFICIENT CONVOLUTIONS :

  procedure Newton_Fabry
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;

    info : integer32;
    idxtol : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 1 ...");
    end if;
    nbrit := 0;
   -- Standard_Newton_Convolution_Steps.LU_Newton_Steps
   --  LU_Newton_Steps was first replaced by a staggered version
   -- Staggered_Newton_Convolutions.LU_Newton_Steps
   --   (hom,prd.sol,rx,ix,maxit,nbrit,tol,absdx,fail,
   --    info,prd.newtpiv,prd.wrk,false,false);
   --  and then by an inlined version
   -- Staggered_Newton_Convolutions.Inlined_LU_Newton_Steps
   --   (hom,prd.sol,rx,ix,maxit,nbrit,tol,absdx,fail,info,prd.newtpiv,
   --    prd.rc,prd.ic,prd.rv,prd.iv,prd.rb,prd.ib,prd.ry,prd.iy,
   --    false,false);
   -- and finally with an indexed version ...
    Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
      (hom,prd.sol,rx,ix,maxit,nbrit,tol,idxtol,absdx,fail,info,prd.newtpiv,
       prd.rc,prd.ic,prd.rv,prd.iv,prd.rb,prd.ib,prd.ry,prd.iy,
       false,false);
    Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;

    info : integer32;
    idxtol : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 2 ...");
    end if;
    nbrit := 0;
    if output then
     -- Standard_Newton_Convolution_Steps.LU_Newton_Steps
     -- Staggered_Newton_Convolutions.LU_Newton_Steps
     --   (file,hom,prd.sol,rx,ix,maxit,nbrit,tol,absdx,fail,
     --    info,prd.newtpiv,prd.wrk,false);
     -- Staggered_Newton_Convolutions.Inlined_LU_Newton_Steps
     --   (file,hom,prd.sol,rx,ix,maxit,nbrit,tol,absdx,fail,info,
     --    prd.newtpiv,prd.rc,prd.ic,prd.rv,prd.iv,prd.rb,prd.ib,prd.ry,prd.iy,
     --    true);
      Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
        (file,hom,prd.sol,rx,ix,maxit,nbrit,tol,idxtol,absdx,fail,info,
         prd.newtpiv,prd.rc,prd.ic,prd.rv,prd.iv,prd.rb,prd.ib,prd.ry,prd.iy,
         true);
      Convergence_Radius_Estimates.Fabry(file,prd.sol,z,rad,err,fail,2);
    else
     -- Standard_Newton_Convolution_Steps.LU_Newton_Steps
     -- Staggered_Newton_Convolutions.LU_Newton_Steps
     --   (hom,prd.sol,rx,ix,maxit,nbrit,tol,absdx,fail,
     --    info,prd.newtpiv,prd.wrk,false,false);
     -- Staggered_Newton_Convolutions.Inlined_LU_Newton_Steps
     --   (hom,prd.sol,rx,ix,maxit,nbrit,tol,absdx,fail,info,
     --    prd.newtpiv,prd.rc,prd.ic,prd.rv,prd.iv,prd.rb,prd.ib,prd.ry,prd.iy,
     --    false,false);
      Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
        (hom,prd.sol,rx,ix,maxit,nbrit,tol,idxtol,absdx,fail,info,
         prd.newtpiv,prd.rc,prd.ic,prd.rv,prd.iv,prd.rb,prd.ib,prd.ry,prd.iy,
         false,false);
      Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    end if;
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

  procedure Newton_Fabry
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;

    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 3 ...");
    end if;
    nbrit := 0;
   -- Standard_Newton_Convolution_Steps.SVD_Newton_Steps
    Staggered_Newton_Convolutions.SVD_Newton_Steps
      (hom,prd.sol,prd.dx,prd.xd,rx,ix,maxit,nbrit,tol,absdx,fail,
       prd.svl,prd.U,prd.V,info,rcond,prd.ewrk,prd.wrk,false,false);
    Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

  procedure Newton_Fabry
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;

    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 4 ...");
    end if;
    nbrit := 0;
    if output then
     -- Standard_Newton_Convolution_Steps.SVD_Newton_Steps
      Staggered_Newton_Convolutions.SVD_Newton_Steps
        (file,hom,prd.sol,prd.dx,prd.xd,rx,ix,maxit,nbrit,tol,absdx,
         fail,prd.svl,prd.U,prd.V,info,rcond,prd.ewrk,prd.wrk,false);
      Convergence_Radius_Estimates.Fabry(file,prd.sol,z,rad,err,fail,2);
    else
     -- Standard_Newton_Convolution_Steps.SVD_Newton_Steps
      Staggered_Newton_Convolutions.SVD_Newton_Steps
        (hom,prd.sol,prd.dx,prd.xd,rx,ix,maxit,nbrit,tol,absdx,fail,
         prd.svl,prd.U,prd.V,info,rcond,prd.ewrk,prd.wrk,false,false);
      Convergence_Radius_Estimates.Fabry(prd.sol,z,rad,err,fail,2,false);
    end if;
    Pade_Vector(prd.numdeg,prd.dendeg,prd.sol,prd.numcff,prd.dencff,
                prd.mat,prd.rhs,prd.padepiv,info,false);
  end Newton_Fabry;

-- NEWTON-FABRY ON COMPLEX CONVOLUTION CIRCUITS :

  procedure Newton_Fabry
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;
    use Standard_Newton_Convolution_Steps;

    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 5 ...");
    end if;
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
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;
    use Standard_Newton_Convolution_Steps;

    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 6 ...");
    end if;
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
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;
    use Standard_Newton_Convolution_Steps;

    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 7 ...");
    end if;
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
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                maxit : in integer32; tol : in double_float;
                nbrit : out integer32; absdx,rcond : out double_float;
                fail : out boolean; z : out Complex_Number;
                rad,err : out double_float; output : in boolean;
                vrblvl : in integer32 := 0 ) is

    use Standard_Rational_Approximations;
    use Standard_Newton_Convolution_Steps;

    info : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Newton_Fabry 8 ...");
    end if;
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

-- HESSE-PADE ON COEFFICIENT CIRCUITS :

  procedure Cached_Singular_Values
              ( cfs : in Standard_Coefficient_Circuits.Link_to_System;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec ) is

    use Standard_Coefficient_Circuits;

    info : integer32;
    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;
    mat : Standard_Complex_Matrices.Link_to_Matrix;
    hrprow,hiprow : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if svh.first then
      Singular_Values(cfs,xr,xi,vh,svh.U,svh.V,svh.ewrk,svls);
      svh.first := false;
    else
      Power_Table(cfs.mxe,xr,xi,cfs.rpwt,cfs.ipwt);
     -- EvalDiff2(cfs.crc,xr,xi,cfs.ryd,cfs.iyd,cfs.rpwt,cfs.ipwt,cfs.jrc.all,
     --           cfs.jic.all,cfs.hrp.all,cfs.hip.all,cfs.fx,cfs.jm,vh);
      for i in cfs.crc'range loop
        if cfs.crc(i).pdg <= 2 then -- no Hessian evaluation needed
         -- code copied from Standard_Coefficient_Circuits.EvalDiff
          Speel(cfs.crc(i).all,xr,xi,cfs.ryd,cfs.iyd,cfs.rpwt,cfs.ipwt);
          cfs.fx(i) := Create(cfs.ryd(0),cfs.iyd(0));
          for j in cfs.jm'range(2) loop
            cfs.jm(i,j) := Create(cfs.ryd(j),cfs.iyd(j));
            rlnk := cfs.jrc(j);    ilnk := cfs.jic(j);
            rlnk(i) := cfs.ryd(j); ilnk(i) := cfs.iyd(j);
            cfs.ryd(j) := 0.0;     cfs.iyd(j) := 0.0;
          end loop;
        else -- for pdg > 2, compute the Hessian matrix
         -- code copied from Standard_Coefficient_Circuits.EvalDiff2
          Speel(cfs.crc(i).all,xr,xi,cfs.ryd,cfs.iyd,cfs.rpwt,cfs.ipwt,
                cfs.hrp.all,cfs.hip.all);
          cfs.fx(i) := Create(cfs.ryd(0),cfs.iyd(0));
          for j in cfs.jm'range(2) loop
            cfs.jm(i,j) := Create(cfs.ryd(j),cfs.iyd(j));
            rlnk := cfs.jrc(j);    ilnk := cfs.jic(j);
            rlnk(i) := cfs.ryd(j); ilnk(i) := cfs.iyd(j);
            cfs.ryd(j) := 0.0;     cfs.iyd(j) := 0.0;
          end loop;
          mat := vh(i);
          for j in cfs.hrp'range loop
            hrprow := cfs.hrp(j); hiprow := cfs.hip(j);
            for k in hrprow'range loop
              mat(j,k) := Create(hrprow(k),hiprow(k));
            end loop;
          end loop;
        end if;
      end loop;
     -- Standard_Complex_Singular_Values.SVD
     --   (cfs.jm,cfs.dim,cfs.dim,svls(0).all,svh.ewrk,svh.U,svh.V,0,info);
      Standard_Matrix_Splitters.Complex_Parts(cfs.jm,svh.xrv,svh.xiv);
      Standard_Inlined_Singular_Values.SVD
        (svh.xrv,svh.xiv,svh.dim,svh.dim,svh.sr,svh.si,svh.er,svh.ei,
         svh.urv,svh.uiv,svh.vrv,svh.viv,0,info,svh.wr,svh.wi);
      for k in vh'range loop
        if cfs.crc(k).pdg > 2 then
         -- Standard_Complex_Singular_Values.SVD
         --   (vh(k).all,cfs.dim,cfs.dim,svls(k).all,
         --    svh.ewrk,svh.U,svh.V,0,info);
          Standard_Matrix_Splitters.Complex_Parts(vh(k).all,svh.xrv,svh.xiv);
          Standard_Inlined_Singular_Values.SVD
            (svh.xrv,svh.xiv,svh.dim,svh.dim,svh.sr,svh.si,svh.er,svh.ei,
             svh.urv,svh.uiv,svh.vrv,svh.viv,0,info,svh.wr,svh.wi);
        end if;
      end loop;
    end if;
  end Cached_Singular_Values;

  procedure Hesse_Pade
              ( cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 1 ...");
    end if;
    Cached_Singular_Values(cfs,svh,xr,xi,vh,svls);
    svh.vals(0) := svls(0)(svh.dim); -- smallest singular value of Jacobian
    for k in 1..svls'last loop
      svh.vals(k) := svls(k)(1);  -- largest singular value of k-th Hessian
    end loop;
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := Standard_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,beta2,eta,nrm);
  end Hesse_Pade;

  procedure Hesse_Pade
              ( file : in file_type;
                cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 2 ...");
    end if;
    Cached_Singular_Values(cfs,svh,xr,xi,vh,svls);
    svh.vals(0) := svls(0)(svh.dim); -- smallest singular value of Jacobian
    for k in 1..svls'last loop
      svh.vals(k) := svls(k)(1);  -- largest singular value of k-th Hessian
    end loop;
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
              ( cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 3 ...");
    end if;
    Cached_Singular_Values(cfs,svh,xr,xi,vh,svls);
    svh.vals(0) := svls(0)(svh.dim); -- smallest singular value of Jacobian
    for k in 1..svls'last loop
      svh.vals(k) := svls(k)(1);  -- largest singular value of k-th Hessian
    end loop;
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := Standard_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,beta2,eta,nrm);
  end Hesse_Pade;

  procedure Hesse_Pade
              ( file : in file_type;
                cfs : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 4 ...");
    end if;
    Cached_Singular_Values(cfs,svh,xr,xi,vh,svls);
    svh.vals(0) := svls(0)(svh.dim); -- smallest singular value of Jacobian
    for k in 1..svls'last loop
      svh.vals(k) := svls(k)(1);  -- largest singular value of k-th Hessian
    end loop;
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

-- HESSE-PADE ON COMPLEX CONVOLUTION CIRCUITS :

  procedure Second
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := svh.dim;
    p : constant integer32 := svh.dim;
    info : integer32;

    use Standard_Complex_Singular_Values;

  begin
    for k in hom.crc'range loop
      svh.H := Hessian_Convolution_Circuits.Hessian(hom.crc(k),sol);
      SVD(svh.H,n,p,svh.svl,svh.ewrk,svh.U,svh.V,0,info,svh.work);
      svh.vals(k) := svh.svl(1);
    end loop;
  end Second;

  procedure Hesse_Pade
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 ) is

    info : integer32;

    use Standard_Complex_Singular_Values;

  begin -- with LU, the system is square so svh.H work space works
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 5 ...");
    end if;
    svh.H := Jacobian_Convolution_Circuits.Jacobian(hom.crc,sol);
    SVD(svh.H,svh.dim,svh.dim,svh.svl,svh.ewrk,svh.U,svh.V,0,info,svh.work);
    svh.vals(0) := svh.svl(svh.dim);
    Second(hom,svh,sol);
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := Standard_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,beta2,eta,nrm);
  end Hesse_Pade;

  procedure Hesse_Pade
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is

    info : integer32;

    use Standard_Complex_Singular_Values;

  begin -- with LU, the system is square so svh.H work space works
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 6 ...");
    end if;
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
    nrm := Standard_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,beta2,eta,nrm);
    if verbose then
      put(file,"eta :"); put(file,eta,3);
      put(file,"  nrm :"); put(file,nrm,3);
      put(file,"  curv_step :"); put(file,step,3); new_line(file);
    end if;
  end Hesse_Pade;

  procedure Hesse_Pade
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 7 ...");
    end if;
    svh.vals(0) := prd.svl(prd.dim);
    Second(hom,svh,sol);
    eta := Distance(svh);
    Homotopy_Pade_Approximants.Solution_Error
      (prd.sol,prd.numcff,prd.dencff,res);
    nrm := Standard_Complex_Vector_Norms.Norm2(res);
    step := Series_and_Predictors.Step_Distance(prd.deg,beta2,eta,nrm);
  end Hesse_Pade;

  procedure Hesse_Pade
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor;
                svh : in Link_to_SVD_Hessians;
                sol : in Standard_Complex_Vectors.Vector;
                res : out Standard_Complex_Vectors.Vector;
                beta2 : in double_float; eta,nrm,step : out double_float;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.Hesse_Pade 8 ...");
    end if;
    svh.vals(0) := prd.svl(prd.dim);
    Second(hom,svh,sol);
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

-- PREDICTOR FEEDBACK LOOPS ON COEFFICIENT SYSTEMS :

  procedure EvalCoeff
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                t : in double_float ) is

    c : Standard_Coefficient_Circuits.Link_to_Circuit;

  begin
    for k in hom.crc'range loop
      c := cfh.crc(k);
      Standard_Coefficient_Convolutions.EvalCoeff
        (hom.crc(k).all,t,c.rcst,c.icst,c.rcf,c.icf);
    end loop;
  end EvalCoeff;

  procedure AbsVal ( cfh : in Standard_Coefficient_Circuits.Link_to_System ) is

    c : Standard_Coefficient_Circuits.Link_to_Circuit;

  begin
    for k in cfh.crc'range loop
      c := cfh.crc(k);
      Standard_Coefficient_Circuits.AbsVal(c.all);
    end loop;
  end AbsVal;

  procedure EvalCffRad
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                t : in double_float ) is

    c,a : Standard_Coefficient_Circuits.Link_to_Circuit;

    use Standard_Mathematical_Functions;

  begin
    EvalCoeff(hom,cfh,t);
    for k in abh.crc'range loop 
      c := cfh.crc(k); -- take coefficients of cfh on input
      a := abh.crc(k); -- assign to the coefficients of abh
      a.rcst := SQRT(c.rcst**2 + c.icst**2);
      a.icst := 0.0;
      for j in 1..c.nbr loop
        a.rcf(j) := SQRT(c.rcf(j)**2 + c.icf(j)**2);
        a.icf(j) := 0.0;
      end loop;
    end loop;
  end EvalCffRad;

  procedure Predictor_Feedback
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in standard_predictor_convolutions.");
      put_line("Predictor_Feedback 1 ...");
    end if;
    nbfail := 0;
    loop
      Standard_Rational_Approximations.Evaluate(numcff,dencff,step,psv.sol);
      EvalCoeff(hom,cfh,step);
      Standard_Vector_Splitters.Complex_Parts(psv.sol,xr,xi);
      Standard_Coefficient_Circuits.Eval(cfh,xr,xi);
      psv.res := cfh.fx;
      nrm := Standard_Complex_Vector_Norms.Max_Norm(psv.res);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      AbsVal(cfh);
      Standard_Vector_Splitters.Complex_Parts(psv.radsol,xr,xi);
      Standard_Coefficient_Circuits.Eval(cfh,xr,xi);
      psv.radres := cfh.fx;
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
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
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in standard_predictor_convolutions.");
      put_line("Predictor_Feedback 2 ...");
    end if;
    nbfail := 0;
    loop
      if verbose then
        put(file,"step in the predictor feedback loop :");
        put(file,step,3); new_line(file);
      end if;
      Standard_Rational_Approximations.Evaluate(numcff,dencff,step,psv.sol);
      EvalCoeff(hom,cfh,step);
      Standard_Vector_Splitters.Complex_Parts(psv.sol,xr,xi);
      Standard_Coefficient_Circuits.Eval(cfh,xr,xi);
      psv.res := cfh.fx;
      nrm := Standard_Complex_Vector_Norms.Max_Norm(psv.res);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      AbsVal(cfh);
      Standard_Vector_Splitters.Complex_Parts(psv.radsol,xr,xi);
      Standard_Coefficient_Circuits.Eval(cfh,xr,xi);
      psv.radres := cfh.fx;
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put_line(file,"Evaluation of the predicted solution : ");
        put_line(file,psv.res);
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
        exit when (step < minstep);
      end if;
    end loop;
  end Predictor_Feedback;

-- PREDICTOR FEEDBACK LOOPS ON CONVOLUTION SYSTEMS :

  procedure Predictor_Feedback
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;

    use Standard_Speelpenning_Convolutions;

  begin
    if vrblvl > 0 then
      put("-> in standard_predictor_convolutions.");
      put_line("Predictor_Feedback 3 ...");
    end if;
    nbfail := 0;
    loop
      Standard_Rational_Approximations.Evaluate(numcff,dencff,step,psv.sol);
      z := Standard_Complex_Numbers.Create(step);
      psv.res := Eval(hom.crc,psv.sol,z);
      nrm := Standard_Complex_Vector_Norms.Max_Norm(psv.res);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol,z);
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
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
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Predictor_Vectors;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                step : in out double_float; minstep,alpha : in double_float;
                nrm,mixres : out double_float; nbfail : out integer32;
                verbose : in boolean := true; vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;

    use Standard_Speelpenning_Convolutions;

  begin
    if vrblvl > 0 then
      put("-> in standard_predictor_convolutions.");
      put_line("Predictor_Feedback 4 ...");
    end if;
    nbfail := 0;
    loop
      if verbose then
        put(file,"step in the predictor feedback loop :");
        put(file,step,3); new_line(file);
      end if;
      Standard_Rational_Approximations.Evaluate(numcff,dencff,step,psv.sol);
      z := Standard_Complex_Numbers.Create(step);
      psv.res := Eval(hom.crc,psv.sol,z);
      nrm := Standard_Complex_Vector_Norms.Max_Norm(psv.res);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol,z);
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put_line(file,"Evaluation of the predicted solution : ");
        put_line(file,psv.res);
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
        exit when (step < minstep);
      end if;
    end loop;
  end Predictor_Feedback;

-- MAIN PREDICTOR PROCEDURES ON COEFFICIENT CONVOLUTION CIRCUITS :

  procedure LU_Prediction
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.LU_Prediction 1 ...");
    end if;
    Newton_Fabry(hom,prd,rx,ix,maxit,tol,nbrit,absdx,fail,z,r,err,vrblvl-1);
    if fail
     then pole_step := maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Standard_Vector_Splitters.Complex_Parts(psv.sol,xr,xi);
    Hesse_Pade(cfh,prd,svh,xr,xi,vh,svls,psv.res,
               beta2,eta,nrm,curv_step,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(hom,cfh,xr,xi,psv,prd.numcff,prd.dencff,
                       step,minstep,alpha,nrm,mixres,nbfail,vrblvl-1);
    fail := (mixres > alpha);
  end LU_Prediction;

  procedure LU_Prediction
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.LU_Prediction 2 ...");
    end if;
    Newton_Fabry(file,hom,prd,rx,ix,maxit,tol,nbrit,absdx,fail,z,r,err,
                 output,vrblvl-1);
    if fail
     then pole_step := maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Standard_Vector_Splitters.Complex_Parts(psv.sol,xr,xi);
    Hesse_Pade(file,cfh,prd,svh,xr,xi,vh,svls,psv.res,beta2,
               eta,nrm,curv_step,verbose,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(file,hom,cfh,xr,xi,psv,prd.numcff,prd.dencff,
                       step,minstep,alpha,nrm,mixres,nbfail,verbose,vrblvl-1);
    fail := (mixres > alpha);
  end LU_Prediction;

  procedure SVD_Prediction
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.SVD_Prediction 1 ...");
    end if;
    Newton_Fabry(hom,prd,rx,ix,maxit,tol,nbrit,absdx,rcond,fail,
                 z,r,err,vrblvl-1);
    if fail
     then pole_step := maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Standard_Vector_Splitters.Complex_Parts(psv.sol,xr,xi);
    Hesse_Pade(cfh,prd,svh,xr,xi,vh,svls,psv.res,
               beta2,eta,nrm,curv_step,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(hom,cfh,xr,xi,psv,prd.numcff,prd.dencff,
                       step,minstep,alpha,nrm,mixres,nbfail,vrblvl-1);
    fail := (mixres > alpha);
  end SVD_Prediction;

  procedure SVD_Prediction
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh : in Standard_Coefficient_Circuits.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.SVD_Prediction 2 ...");
    end if;
    Newton_Fabry(file,hom,prd,rx,ix,maxit,tol,nbrit,absdx,rcond,
                 fail,z,r,err,output,vrblvl-1);
    if fail
     then pole_step := maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Standard_Vector_Splitters.Complex_Parts(psv.sol,xr,xi);
    Hesse_Pade(file,cfh,prd,svh,xr,xi,vh,svls,psv.res,beta2,
               eta,nrm,curv_step,verbose,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(file,hom,cfh,xr,xi,psv,prd.numcff,prd.dencff,
                       step,minstep,alpha,nrm,mixres,nbfail,verbose,vrblvl-1);
    fail := (mixres > alpha);
  end SVD_Prediction;

-- MAIN PREDICTOR PROCEDURES ON COMPLEX CONVOLUTION CIRCUITS :

  procedure LU_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.LU_Prediction 3 ...");
    end if;
    Newton_Fabry(hom,prd,maxit,tol,nbrit,absdx,fail,z,r,err,vrblvl-1);
    if fail
     then pole_step := maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Hesse_Pade(hom,prd,svh,psv.sol,psv.res,beta2,eta,nrm,curv_step,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(hom,abh,psv,prd.numcff,prd.dencff,
                       step,minstep,alpha,nrm,mixres,nbfail,vrblvl-1);
    fail := (mixres > alpha);
  end LU_Prediction;

  procedure LU_Prediction
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_LU_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.LU_Prediction 4 ...");
    end if;
    Newton_Fabry(file,hom,prd,maxit,tol,nbrit,absdx,fail,
                 z,r,err,output,vrblvl-1);
    if fail
     then pole_step := maxstep;
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
               eta,nrm,curv_step,verbose,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(file,hom,abh,psv,prd.numcff,prd.dencff,
                       step,minstep,alpha,nrm,mixres,nbfail,verbose,vrblvl-1);
    fail := (mixres > alpha);
  end LU_Prediction;

  procedure SVD_Prediction
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.SVD_Prediction 3 ...");
    end if;
    Newton_Fabry(hom,prd,maxit,tol,nbrit,absdx,rcond,fail,z,r,err,vrblvl-1);
    if fail
     then pole_step := maxstep;
     else pole_step := beta1*r;
    end if;
    for k in prd.sol'range loop
      lnk := prd.sol(k); psv.sol(k) := lnk(0);
    end loop;
    Hesse_Pade(hom,prd,svh,psv.sol,psv.res,beta2,eta,nrm,curv_step,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(hom,abh,psv,prd.numcff,prd.dencff,
      step,minstep,alpha,nrm,mixres,nbfail,vrblvl-1);
    fail := (mixres > alpha);
  end SVD_Prediction;

  procedure SVD_Prediction
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                prd : in Link_to_SVD_Predictor; svh : in Link_to_SVD_Hessians;
                psv : in out Predictor_Vectors;
                maxit : in integer32; tol : in double_float;
                alpha,beta1,beta2,maxstep,minstep,endt : in double_float;
                acct : in out double_float;
                fail : out boolean; step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                output : in boolean := false;
                verbose : in boolean := false;
                vrblvl : in integer32 := 0 ) is

    z : Standard_Complex_Numbers.Complex_Number;
    r,err,absdx,rcond,pole_step,eta,nrm,curv_step,mixres : double_float;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    nbrit,nbfail : integer32;

    use Three_Way_Minima;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_predictor_convolutions.SVD_Prediction 4 ...");
    end if;
    Newton_Fabry
      (file,hom,prd,maxit,tol,nbrit,absdx,rcond,fail,z,r,err,output,vrblvl-1);
    if fail
     then pole_step := maxstep;
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
               eta,nrm,curv_step,verbose,vrblvl-1);
    Minimum(pole_step,curv_step,maxstep,step,nbpole,nbhess,nbmaxm);
    Bounded_Update(acct,step,endt,minstep);
    Predictor_Feedback(file,hom,abh,psv,prd.numcff,prd.dencff,
      step,minstep,alpha,nrm,mixres,nbfail,verbose);
    fail := (mixres > alpha);
  end SVD_Prediction;

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
      Standard_Floating_VecVecs.Deep_Clear(p.rc);
      Standard_Floating_VecVecs.Deep_Clear(p.ic);
      Standard_Floating_VecVecVecs.Clear(p.rv);
      Standard_Floating_VecVecVecs.Clear(p.iv);
      Standard_Floating_VecVecs.Deep_Clear(p.rb);
      Standard_Floating_VecVecs.Deep_Clear(p.ib);
      Standard_Floating_Vectors.Clear(p.ry);
      Standard_Floating_Vectors.Clear(p.iy);
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

  procedure Clear ( p : in out Predictor ) is
  begin
    case p.kind is
      when LU  => Clear(p.ludata);
      when SVD => Clear(p.svdata);
    end case;
  end Clear;

  procedure Clear ( p : in out Link_to_Predictor_Vectors ) is

    procedure free is
      new unchecked_deallocation(Predictor_Vectors,Link_to_Predictor_Vectors);

  begin
    if p /= null
     then free(p);
    end if;
  end Clear;

  procedure Clear ( h : in out SVD_Hessians ) is
  begin
    Standard_Floating_Vectors.Clear(h.sr);
    Standard_Floating_Vectors.Clear(h.si);
    Standard_Floating_Vectors.Clear(h.er);
    Standard_Floating_Vectors.Clear(h.ei);
    Standard_Floating_Vectors.Clear(h.wr);
    Standard_Floating_Vectors.Clear(h.wi);
    Standard_Floating_VecVecs.Deep_Clear(h.xrv);
    Standard_Floating_VecVecs.Deep_Clear(h.xiv);
   -- not allocated because job = 0
   -- Standard_Floating_VecVecs.Deep_Clear(h.urv);
   -- Standard_Floating_VecVecs.Deep_Clear(h.uiv);
   -- Standard_Floating_VecVecs.Deep_Clear(h.vrv);
   -- Standard_Floating_VecVecs.Deep_Clear(h.viv);
  end Clear;

  procedure Clear ( h : in out Link_to_SVD_Hessians ) is

    procedure free is
      new unchecked_deallocation(SVD_Hessians,Link_to_SVD_Hessians);

  begin
    if h /= null then
      Clear(h.all);
      free(h);
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

  procedure Clear ( p : in out Predictor_Array ) is
  begin
    for i in p'range loop
      case p(i).kind is
        when LU => Clear(p(i).ludata);
        when SVD => Clear(p(i).svdata);
      end case;
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

end Standard_Predictor_Convolutions;
