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
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;    use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;    use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;       use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;       use QuadDobl_Complex_Solutions_io;
with Partitions_of_Sets_of_Unknowns;      use Partitions_of_Sets_of_Unknowns;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Homotopy;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy;
with QuadDobl_Homotopy_Convolutions_io;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Shift_Convolution_Circuits;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
with Corrector_Convolutions;             use Corrector_Convolutions;
with Predictor_Corrector_Loops;          use Predictor_Corrector_Loops;
with Series_Path_Trackers;

procedure ts_pcscnv is

-- DESCRIPTION :
--   Development of one predictor-corrector-shift step with
--   a homotopy system of convolution circuits.

  procedure Step_Track
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                maxit : in integer32; hcrd : in boolean;
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
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used;
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
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0;
    loop
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,maxit,hcrd,prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,
        nbpole,nbhess,nbmaxm,fail,verbose);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
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
                maxit : in integer32; hcrd : in boolean;
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
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used;
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
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0;
    loop
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,maxit,hcrd,prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,
        nbpole,nbhess,nbmaxm,fail,verbose);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
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
                maxit : in integer32; hcrd : in boolean;
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
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used;
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
    nbrit : integer32;

  begin
    nbpole := 0; nbhess := 0; nbmaxm := 0;
    loop
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,maxit,hcrd,prd,psv,svh,dx,ipvt,endt,acct,step,mixres,nbrit,
        nbpole,nbhess,nbmaxm,fail,verbose);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      Shift_Convolution_Circuits.Shift(hom,wrk,-step);
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
                hcrd : in boolean ) is

  -- DESCRIPTION :
  --   Runs predictor-corrector-shift loops in double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used.

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
        Step_Track(hom,abh,homlead,abhlead,pars,maxit,hcrd,prd,psv,svh,
          dx,ipvt,wrk,nbpole,nbhess,nbmaxm,fail,true);
      else   
        Track_One_Path(standard_output,hom,abh,homlead,abhlead,pars,maxit,
          hcrd,prd,psv,svh,dx,ipvt,wrk,t,mixres,tnbrit,nbpole,nbhess,nbmaxm,
          nbsteps,minstpz,maxstpz,fail,true);
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

  procedure DoblDobl_Run_Loops
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                hcrd : in boolean ) is

  -- DESCRIPTION :
  --   Runs predictor-corrector-shift loops in double double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used.

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
        Step_Track(hom,abh,homlead,abhlead,pars,maxit,hcrd,prd,psv,svh,
          dx,ipvt,wrk,nbpole,nbhess,nbmaxm,fail,true);
      else   
        Track_One_Path(standard_output,hom,abh,homlead,abhlead,pars,maxit,
          hcrd,prd,psv,svh,dx,ipvt,wrk,t,mixres,tnbrit,nbpole,nbhess,nbmaxm,
          nbsteps,minstpz,maxstpz,fail,true);
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
                hcrd : in boolean ) is

  -- DESCRIPTION :
  --  Runs predictor-corrector-shift loops in quad double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used.

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
        Step_Track(hom,abh,homlead,abhlead,pars,maxit,hcrd,prd,psv,svh,
          dx,ipvt,wrk,nbpole,nbhess,nbmaxm,fail,true);
      else   
        Track_One_Path(standard_output,hom,abh,homlead,abhlead,pars,maxit,
          hcrd,prd,psv,svh,dx,ipvt,wrk,t,mixres,tnbrit,nbpole,nbhess,nbmaxm,
          nbsteps,minstpz,maxstpz,fail,true);
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

  procedure Standard_Track_All
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                hcrd,arth : in boolean ) is

  -- DESCRIPTION :
  --   Tracks all paths defined by the homotopy in hom,
  --   starting at solutions in sols, in double precision.

  -- REQUIRED : the homotopy is square.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used;
  --   arth     true if the homotopy is an artificial-parameter one,
  --            false otherwise.

  -- ON RETURN :
  --   sols     solutions at the end of the paths.

    file : file_type;
    verbose : boolean;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    if not arth then
      put(file,natural32(hom.neq),natural32(hom.neq+1),
               Standard_Homotopy.Homotopy_System);
    else
      declare
        p : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Homotopy.Target_System;
        q : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
             natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(file);
    Track_All_Paths(file,hom,abh,sols,pars,hcrd,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
             natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
  end Standard_Track_All;

  procedure DoblDobl_Track_All
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                hcrd,arth : in boolean ) is

  -- DESCRIPTION :
  --   Tracks all paths defined by the homotopy in hom,
  --   starting at solutions in sols, in double double precision.

  -- REQUIRED : the homotopy is square.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used;
  --   arth     true if the homotopy is an artificial-parameter one,
  --            false otherwise.

  -- ON RETURN :
  --   sols     solutions at the end of the path.

    file : file_type;
    verbose : boolean;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    if not arth then
      put(file,natural32(hom.neq),natural32(hom.neq+1),
               DoblDobl_Homotopy.Homotopy_System);
    else
      declare
        p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := DoblDobl_Homotopy.Target_System;
        q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := DoblDobl_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
             natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(file);
    Track_All_Paths(file,hom,abh,sols,pars,hcrd,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
             natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end DoblDobl_Track_All;

  procedure QuadDobl_Track_All
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                hcrd,arth : in boolean ) is

  -- DESCRIPTION :
  --   Tracks all paths defined by the homotopy in hom,
  --   starting at solutions in sols, in quad double precision.

  -- REQUIRED : the homotopy is square.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   hcrd     true if the homotopy is 1-homogeneous with one linear
  --            equation as its last equation which is then updated in
  --            the scaling of the solution coordinates,
  --            false if affine coordinates are used;
  --   arth     true if the homotopy is an artificial-parameter one,
  --            false otherwise.

  -- ON RETURN :
  --   sols     solutions at the end of the path.

    file : file_type;
    verbose : boolean;
    ans : character;

  begin
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    if not arth then
      put(file,natural32(hom.neq),natural32(hom.neq+1),
               QuadDobl_Homotopy.Homotopy_System);
    else
      declare
        p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := QuadDobl_Homotopy.Target_System;
        q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := QuadDobl_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
             natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
    new_line(file);
    Track_All_Paths(file,hom,abh,sols,pars,hcrd,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
             natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
  end QuadDobl_Track_All;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double precision.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    ans : character;
    artificial : constant boolean  -- flag for artificial-parameter homotopy
               := Series_Path_Trackers.Prompt_for_Artificial;
    homogeneous : boolean; -- flag for homogeneous coordinates
    mhom : natural32 := 0;
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;

  begin
    if not artificial
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    Standard_Homotopy_Convolutions_io.get
      (deg,artificial,pars.gamma,cnvhom,sols,idxpar,mhom,z,idz);
    if mhom > 1 then
      put_line("General m-homogenization is not (yet) supported.");
    else
      homogeneous := (mhom = 1);
      abshom := Residual_Convolution_System(cnvhom);
      if artificial
       then pars.gamma := Standard_Homotopy.Accessibility_Constant;
      end if;
      new_line;
      put("Step-by-step runs ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then Standard_Run_Loops(cnvhom,abshom,sols,pars,homogeneous);
       else Standard_Track_All(cnvhom,abshom,sols,pars,homogeneous,artificial);
      end if;
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double double precisin.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    ddgamma : DoblDobl_Complex_Numbers.Complex_Number;
    ans : character;
    artificial : constant boolean  -- flag for artificial-parameter homotopy
               := Series_Path_Trackers.Prompt_for_Artificial;
    homogeneous : boolean; -- flag for homogeneous coordinates
    mhom : natural32 := 0;
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;

    use DoblDobl_Complex_Numbers_cv;
  
  begin
    if not artificial
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    DoblDobl_Homotopy_Convolutions_io.get
      (deg,artificial,pars.gamma,cnvhom,sols,idxpar,mhom,z,idz);
    if mhom > 1 then
      put_line("General m-homogenization is not (yet) supported.");
    else
      homogeneous := (mhom = 1);
      abshom := Residual_Convolution_System(cnvhom);
      if artificial then
        ddgamma := DoblDobl_Homotopy.Accessibility_Constant;
        pars.gamma := DoblDobl_Complex_to_Standard(ddgamma);
      end if;
      new_line;
      put("Step-by-step runs ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then DoblDobl_Run_Loops(cnvhom,abshom,sols,pars,homogeneous);
       else DoblDobl_Track_All(cnvhom,abshom,sols,pars,homogeneous,artificial);
      end if;
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in quad double precision.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    qdgamma : QuadDobl_Complex_Numbers.Complex_Number;
    ans : character;
    artificial : constant boolean  -- flag for artificial-parameter homotopy
               := Series_Path_Trackers.Prompt_for_Artificial;
    homogeneous : boolean; -- flag for homogeneous coordinates
    mhom : natural32 := 0;
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;

    use QuadDobl_Complex_Numbers_cv;

  begin
    if not artificial
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    QuadDobl_Homotopy_Convolutions_io.get
      (deg,artificial,pars.gamma,cnvhom,sols,idxpar,mhom,z,idz);
    if mhom > 1 then
      put_line("General m-homogenization is not (yet) supported.");
    else
      homogeneous := (mhom = 1);
      abshom := Residual_Convolution_System(cnvhom);
      if artificial then
        qdgamma := QuadDobl_Homotopy.Accessibility_Constant;
        pars.gamma := QuadDobl_Complex_to_Standard(qdgamma);
      end if;
      new_line;
      put("Step-by-step runs ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then QuadDobl_Run_Loops(cnvhom,abshom,sols,pars,homogeneous);
       else QuadDobl_Track_All(cnvhom,abshom,sols,pars,homogeneous,artificial);
      end if;
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
