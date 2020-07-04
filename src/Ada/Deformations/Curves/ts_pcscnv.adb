with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with Standard_Vector_Splitters;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with Standard_Coefficient_Convolutions;
with Standard_Circuit_Makers;
with Standard_Convolution_Splitters;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Homotopy_Continuation_Parameters;
with Shift_Convolution_Circuits;
with Shift_Coefficient_Convolutions;
with Standard_Coefficient_Storage;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Predictor_Corrector_Loops;          use Predictor_Corrector_Loops;
with Predictor_Corrector_Trackers;       use Predictor_Corrector_Trackers;
with Track_Path_Convolutions;

procedure ts_pcscnv is

-- DESCRIPTION :
--   Development of one predictor-corrector-shift step with
--   a homotopy system of convolution circuits.

  procedure Step_Track
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                nbpole,nbhess,nbmaxm : out natural32;
                fail : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Tracks one path step by step, interactively, prompting the user each
  --   time before moving on to the next step, in double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters;
  --   maxit    maximum number of steps in Newton's method on power series;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   wrk      work space vector for power series coefficients
  --            during the shifting of the coefficients.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

    endt : constant double_float := 1.0;
    acct,step,mixres : double_float := 0.0;
    ans : character;
    cnt,nbrit : integer32 := 0;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0;
    loop
      cnt := cnt + 1;
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,maxit,mhom,idz,prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,
        nbpole,nbhess,nbmaxm,fail,verbose);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
      put("At step "); put(cnt,1);
      put(", t :"); put(acct,3); put(", step :"); put(step,3); new_line;
      put("Do the next step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("t :"); put(acct,3); put_line(" :");
    end loop;
  end Step_Track;

  procedure Step_Track
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                cfs,abh : in Standard_Coefficient_Circuits.Link_to_System;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                vh : in Standard_Complex_VecMats.VecMat;
                svls : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                pwt : in Standard_Floating_Vectors.Link_to_Vector;
                nbpole,nbhess,nbmaxm : out natural32;
                fail : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Tracks one path step by step, interactively, prompting the user each
  --   time before moving on to the next step, in double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   rcf      real parts of all coefficients in the original homotopy;
  --   icf      imaginary parts of all coefficients in the original homotopy;
  --   cfs      coefficient circuits for the corrector loop;
  --   abh      radii as coefficients of cfs for mixed residuals;
  --   pars     values for the tolerances and parameters;
  --   maxit    maximum number of steps in Newton's method on power series;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   rx       work space for the real parts of the series coefficients;
  --   ix       work space for the imaginary parts of the series coefficients;
  --   xr       work space allocated for the real parts of solution vectors;
  --   xi       work space allocated for the imag parts of solution vectors;
  --   vh       space allocated for dim matrices, dim = dimension,
  --            all matrices have 1..dim for range(1) and range(2);
  --   svls     svls(0) contains the singular values of s.jm, and
  --            svls(k) contains the singular values of vh(k),
  --            for k in vh'range.
  --   pwt      work space the powers of the value used in the shift,
  --            pwt'range = 0..hom.deg.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

    endt : constant double_float := 1.0;
    acct,step,mixres : double_float := 0.0;
    ans : character;
    cnt,nbrit : integer32 := 0;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0;
    loop
      cnt := cnt + 1;
      Predictor_Corrector_Loop(standard_output,hom,cfs,abh,
        pars,maxit,mhom,idz,prd,psv,svh,rx,ix,xr,xi,vh,svls,ipvt,
        endt,acct,step,mixres,nbrit,nbpole,nbhess,nbmaxm,fail,verbose);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      Shift_Coefficient_Convolutions.Powers_of_Shift(pwt,-acct);
      Shift_Coefficient_Convolutions.Map(rcf,icf,hom,pwt);
      Standard_Predictor_Convolutions.EvalCffRad(hom,cfs,abh,0.0);
      put("At step "); put(cnt,1);
      put(", t :"); put(acct,3); put(", step :"); put(step,3); new_line;
      put("Do the next step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("t :"); put(acct,3); put_line(" :");
    end loop;
  end Step_Track;

  procedure Step_Track
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                nbpole,nbhess,nbmaxm : out natural32;
                fail : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Tracks one path step by step, interactively, prompting the user each
  --   time before moving on to the next step, in double double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters;
  --   maxit    maximum number of steps in Newton's method on power series;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   wrk      work space vector for power series coefficients
  --            during the shifting of the coefficients.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

    endt : constant double_float := 1.0;
    acct,step,mixres : double_double := Create(0.0);
    ans : character;
    cnt,nbrit : integer32 := 0;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0;
    loop
      cnt := cnt + 1;
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,maxit,mhom,idz,prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,
        nbpole,nbhess,nbmaxm,fail,verbose);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
      put("At step "); put(cnt,1);
      put(", t : "); put(acct,3); put(", step : "); put(step,3); new_line;
      put("Do the next step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("t : "); put(acct,3); put_line(" :");
    end loop;
  end Step_Track;

  procedure Step_Track
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                nbpole,nbhess,nbmaxm : out natural32;
                fail : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Tracks one path step by step, interactively, prompting the user each
  --   time before moving on to the next step, in quad double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters;
  --   maxit    maximum number of steps in Newton's method on power series;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   wrk      work space vector for power series coefficients
  --            during the shifting of the coefficients.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

    endt : constant double_float := 1.0;
    acct,step,mixres : quad_double := Create(0.0);
    ans : character;
    cnt,nbrit : integer32 := 0;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0;
    loop
      cnt := cnt + 1;
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,maxit,mhom,idz,prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,
        nbpole,nbhess,nbmaxm,fail,verbose);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
      put("At step "); put(cnt,1);
      put(", t : "); put(acct,3); put(", step : "); put(step,3); new_line;
      put("Do the next step ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("t : "); put(acct,3); put_line(" :");
    end loop;
  end Step_Track;

  procedure Standard_Run_Loops
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Runs predictor-corrector-shift loops in double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to.

    use Standard_Complex_Solutions,Standard_Predictor_Convolutions;

    maxit : constant integer32 := 4; -- max #steps in Newton on Power Series
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail,stepwise : boolean;
    ans : character;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : Standard_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : Standard_Complex_VecVecs.Link_to_VecVec;
    wrk : constant Standard_Complex_Vectors.Link_to_Vector
        := Standard_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
    t,mixres,minstpz,maxstpz : double_float := 0.0;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    put("Interactive step-by-step run ? (y/n) "); Ask_Yes_or_No(ans);
    stepwise := (ans = 'y');
    loop
      ls := Head_Of(solsptr); psv.sol := ls.v; t := 0.0;
      if stepwise then
        Step_Track(hom,abh,homlead,abhlead,pars,maxit,mhom,idz,prd,psv,svh,
          dx,ipvt,wrk,nbpole,nbhess,nbmaxm,fail,true);
      else   
        Track_One_Path(standard_output,hom,abh,homlead,abhlead,pars,maxit,
          mhom,idz,prd,psv,svh,dx,ipvt,wrk,t,mixres,tnbrit,nbpole,nbhess,
          nbmaxm,nbsteps,minstpz,maxstpz,fail,true);
      end if;
      ls.v := psv.sol; ls.res := mixres;
      ls.t := Standard_Complex_Numbers.Create(t); Set_Head(solsptr,ls);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
    Clear(svh);
    Standard_Complex_VecVecs.Deep_Clear(homlead);
    Standard_Complex_VecVecs.Deep_Clear(abhlead);
    Standard_Speelpenning_Convolutions.Clear(homcff);
  end Standard_Run_Loops;

  procedure Standard_Run_Loops
              ( hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfs,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Runs predictor-corrector-shift loops in double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   cfs      coefficient circuits for the corrector;
  --   abh      radii as coefficients of cfs, for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to.

    use Standard_Complex_Solutions,Standard_Predictor_Convolutions;

    maxit : constant integer32 := 4; -- max #steps in Newton on Power Series
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    dim : constant integer32 := hom.dim;
    deg : constant integer32 := hom.deg;
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(dim);
    solsptr : Solution_List := sols;
    tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail,stepwise : boolean;
    ans : character;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    t,mixres,minstpz,maxstpz : double_float := 0.0;
    rx : Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ix : Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    vxr : constant Standard_Floating_Vectors.Vector(ls.v'range)
        := (ls.v'range => 0.0);
    vxi : constant Standard_Floating_Vectors.Vector(ls.v'range)
        := (ls.v'range => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxr);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(vxi);
    vh : Standard_Complex_VecMats.VecMat(1..hom.neq)
       := Standard_Complex_Circuits.Allocate(hom.neq,dim);
    svls : Standard_Complex_VecVecs.VecVec(0..dim)
         := Standard_Vector_Splitters.Allocate(hom.neq,dim+1,0,1);
    pwt : Standard_Floating_Vectors.Link_to_Vector;
    rcfhom,icfhom : Standard_Floating_VecVecVecs.Link_to_VecVecVec;

  begin
    Standard_Coefficient_Storage.Allocate_and_Store(hom.crc,rcfhom,icfhom);
    pwt := Standard_Vector_Splitters.Allocate_Floating_Coefficients(deg);
    put("Interactive step-by-step run ? (y/n) "); Ask_Yes_or_No(ans);
    stepwise := (ans = 'y');
    loop
      ls := Head_Of(solsptr); psv.sol := ls.v; t := 0.0;
      if stepwise then
        Step_Track
          (hom,rcfhom,icfhom,cfs,abh,pars,maxit,mhom,idz,prd,psv,svh,
           rx,ix,xr,xi,vh,svls,ipvt,pwt,nbpole,nbhess,nbmaxm,fail,true);
      else   
        Track_One_Path(standard_output,hom,rcfhom,icfhom,cfs,abh,pars,maxit,
          mhom,idz,prd,psv,svh,rx,ix,xr,xi,vh,svls,ipvt,pwt,t,mixres,tnbrit,
          nbpole,nbhess,nbmaxm,nbsteps,minstpz,maxstpz,fail,true);
      end if;
      ls.v := psv.sol; ls.res := mixres;
      ls.t := Standard_Complex_Numbers.Create(t); Set_Head(solsptr,ls);
      put_line("The solution vector :"); put_line(ls.v);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Standard_Coefficient_Storage.Restore(rcfhom,icfhom,hom.crc);
    end loop;
    Clear(svh);
    Standard_Floating_VecVecs.Deep_Clear(rx);
    Standard_Floating_VecVecs.Deep_Clear(ix);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
    Standard_Floating_Vectors.Clear(pwt);
    Standard_Complex_VecMats.Clear(vh);
    Standard_Complex_VecVecs.Clear(svls);
  end Standard_Run_Loops;

  procedure DoblDobl_Run_Loops
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Runs predictor-corrector-shift loops in double double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to.

    use DoblDobl_Complex_Solutions,DoblDobl_Predictor_Convolutions;

    maxit : constant integer32 := 4; -- max #steps in Newton on Power Series
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail,stepwise : boolean;
    ans : character;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : DoblDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
    t,mixres : double_double;
    minstpz,maxstpz : double_float;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    put("Interactive step-by-step run ? (y/n) "); Ask_Yes_or_No(ans);
    stepwise := (ans = 'y');
    loop
      ls := Head_Of(solsptr); psv.sol := ls.v; t := Create(0.0);
      if stepwise then
        Step_Track(hom,abh,homlead,abhlead,pars,maxit,mhom,idz,prd,psv,svh,
          dx,ipvt,wrk,nbpole,nbhess,nbmaxm,fail,true);
      else   
        Track_One_Path(standard_output,hom,abh,homlead,abhlead,pars,maxit,
          mhom,idz,prd,psv,svh,dx,ipvt,wrk,t,mixres,tnbrit,nbpole,nbhess,
          nbmaxm,nbsteps,minstpz,maxstpz,fail,true);
      end if;
      ls.v := psv.sol; ls.res := mixres;
      ls.t := DoblDobl_Complex_Numbers.Create(t); Set_Head(solsptr,ls);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
    Clear(svh);
    DoblDobl_Complex_VecVecs.Deep_Clear(homlead);
    DoblDobl_Complex_VecVecs.Deep_Clear(abhlead);
  end DoblDobl_Run_Loops;

  procedure QuadDobl_Run_Loops
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --  Runs predictor-corrector-shift loops in quad double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   mhom     0 if affine coordinates are used,
  --            1 for 1-homogeneous coordinates,
  --            m, for m > 1, for multi-homogenization;
  --   idz      the index representation of the partition of the variables,
  --            idz(k) returns a value between 1 and m,
  --            depending on which set the k-th variable belongs to.

    use QuadDobl_Complex_Solutions,QuadDobl_Predictor_Convolutions;

    maxit : constant integer32 := 4; -- max #steps in Newton on Power Series
    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    ls : Link_to_Solution := Head_Of(sols);
    prd : Predictor := Create(ls.v,hom.neq,hom.deg,numdeg,dendeg,SVD);
    psv : Predictor_Vectors(hom.dim,hom.neq);
    svh : Link_to_SVD_Hessians := Create(hom.dim);
    solsptr : Solution_List := sols;
    tnbrit,nbpole,nbhess,nbmaxm,nbsteps : natural32 := 0;
    fail,stepwise : boolean;
    ans : character;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : QuadDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := QuadDobl_Speelpenning_Convolutions.Allocate_Coefficients(hom.deg);
    homcff : QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
    t,mixres : quad_double;
    minstpz,maxstpz : double_float;

  begin
    Allocate_Coefficients(hom.crc,homcff);
    Store_Coefficients(hom.crc,homcff);
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    put("Interactive step-by-step run ? (y/n) "); Ask_Yes_or_No(ans);
    stepwise := (ans = 'y');
    loop
      ls := Head_Of(solsptr); psv.sol := ls.v; t := create(0.0);
      if stepwise then
        Step_Track(hom,abh,homlead,abhlead,pars,maxit,mhom,idz,prd,psv,svh,
          dx,ipvt,wrk,nbpole,nbhess,nbmaxm,fail,true);
      else   
        Track_One_Path(standard_output,hom,abh,homlead,abhlead,pars,maxit,
          mhom,idz,prd,psv,svh,dx,ipvt,wrk,t,mixres,tnbrit,nbpole,nbhess,
          nbmaxm,nbsteps,minstpz,maxstpz,fail,true);
      end if;
      ls.v := psv.sol; ls.res := mixres;
      ls.t := QuadDobl_Complex_Numbers.Create(t); Set_Head(solsptr,ls);
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Restore_Leading_Coefficients(abhlead,abh.crc);
      Restore_Coefficients(homcff,hom.crc);
    end loop;
    Clear(svh);
    QuadDobl_Complex_VecVecs.Deep_Clear(homlead);
    QuadDobl_Complex_VecVecs.Deep_Clear(abhlead);
  end QuadDobl_Run_Loops;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double precision,
  --   and then launches path tracker tests.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial,stepbystep,oncffcirc : boolean;
    ans : character;
    file : file_type;

  begin
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    new_line;
    put("Step-by-step runs ? (y/n) "); Ask_Yes_or_No(ans);
    stepbystep := (ans = 'y');
    new_line;
    put("Run on coefficient circuits ? (y/n) "); Ask_Yes_or_No(ans);
    oncffcirc := (ans = 'y');
    if not oncffcirc then
      if stepbystep then
        Standard_Run_Loops(cnvhom,abshom,sols,pars,integer32(mhom),idz);
      else
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(file);
        Track_Path_Convolutions.Track
          (file,cnvhom,abshom,sols,pars,integer32(mhom),idz,artificial);
      end if;
    else
      declare
        cffhom : Standard_Coefficient_Convolutions.Link_to_System;
        cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
      begin
        cffhom := Standard_Convolution_Splitters.Split(cnvhom);
        cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
        abh := Standard_Coefficient_Circuits.Copy(cfs);
        Standard_Coefficient_Circuits.AbsVal(abh);
        if stepbystep then
          Standard_Run_Loops(cffhom,cfs,abh,sols,pars,integer32(mhom),idz);
        else
          new_line;
          put_line("Reading the name of the output file ...");
          Read_Name_and_Create_File(file);
          Track_Path_Convolutions.Track
            (file,cffhom,cfs,abh,sols,pars,integer32(mhom),idz,artificial);
        end if;
      end;
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double double precision,
  --   and then launches path tracker tests.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    ans : character;
    file : file_type;

  begin
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    new_line;
    put("Step-by-step runs ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      DoblDobl_Run_Loops(cnvhom,abshom,sols,pars,integer32(mhom),idz);
    else
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
      Track_Path_Convolutions.Track
        (file,cnvhom,abshom,sols,pars,integer32(mhom),idz,artificial);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in quad double precision,
  --   and then launches path tracker tests.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    ans : character;
    file : file_type;

  begin
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    new_line;
    put("Step-by-step runs ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      QuadDobl_Run_Loops(cnvhom,abshom,sols,pars,integer32(mhom),idz);
    else
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file);
      Track_Path_Convolutions.Track
        (file,cnvhom,abshom,sols,pars,integer32(mhom),idz,artificial);
    end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_pcscnv;
