with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
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
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
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
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
with Corrector_Convolutions;             use Corrector_Convolutions;

procedure ts_pcscnv is

-- DESCRIPTION :
--   Development of one predictor-corrector-shift step with
--   a homotopy system of convolution circuits.

  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in Standard_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                prd : in out Standard_Predictor_Convolutions.Predictor;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                svh : in Standard_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                step : out double_float;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Does one predictor-corrector step in double precision.

  -- ON ENTRY :
  --   file     to write the extra output to;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the Hessian step was minimal;
  --   nbmaxm   number of times the maximum step was minimal;
  --   verbose  flag for extra output.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution
  --   ipvt     pivoting information for the LU Newton steps;
  --   step     the step size;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

    use Standard_Predictor_Convolutions;

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    deg : constant integer32 := numdeg + dendeg + 2;
    maxit : constant integer32 := 4;
    info,nbrit : integer32 := 0;
    mixres : double_float;

  begin
    prd := Create(psv.sol,hom.neq,deg,numdeg,dendeg,SVD);
    Set_Lead_Coefficients(prd,psv.sol);
    svh.vals := (svh.vals'range => Standard_Complex_Numbers.Create(0.0));
    SVD_Prediction(file,hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,fail,step,
      nbpole,nbhess,nbmaxm,false,verbose);
    if verbose then
      if fail
       then put(file,"Predictor failed to reach tolerance");
       else put(file,"Predictor reached tolerance");
      end if;
      put(file,pars.alpha,3); put_line(file,".");
    end if;
    loop
      Step_Coefficient(hom,step);
      Step_Coefficient(abh,step);
      LU_Newton_Steps(file,hom,abh,psv,maxit,nbrit,pars.tolres,mixres,
                      dx,ipvt,info,fail,verbose);
      exit when not fail;   
      step := step/2.0;
      put(file,"Reduced step size to"); put(file,step,3); put_line(file,".");
      exit when (step < pars.minsize);
      Standard_Rational_Approximations.Evaluate
        (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      Restore_Leading_Coefficients(homlead,hom.crc);
      Restore_Leading_Coefficients(abhlead,abh.crc);
    end loop;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                prd : in out DoblDobl_Predictor_Convolutions.Predictor;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in DoblDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                step : out double_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Does one predictor-corrector step in double double precision.

  -- ON ENTRY :
  --   file     to write the extra output to;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the Hessian step was minimal;
  --   nbmaxm   number of times the maximum step was minimal;
  --   verbose  flag for extra output.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   step     the step size;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

    use DoblDobl_Predictor_Convolutions;

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    deg : constant integer32 := numdeg + dendeg + 2;
    maxit : constant integer32 := 4;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    info,nbrit : integer32 := 0;
    mixres : double_double;

  begin
    prd := Create(psv.sol,hom.neq,deg,numdeg,dendeg,SVD);
    Set_Lead_Coefficients(prd,psv.sol);
    svh.vals := (svh.vals'range => zero);
    SVD_Prediction(file,hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,fail,step,
      nbpole,nbhess,nbmaxm,false,verbose);
    if verbose then
      if fail
       then put(file,"Predictor failed to reach tolerance");
       else put(file,"Predictor reached tolerance");
      end if;
      put(file,pars.alpha,3); put_line(file,".");
    end if;
    loop
      Step_Coefficient(hom,step);
      Step_Coefficient(abh,step);
      LU_Newton_Steps(file,hom,abh,psv,maxit,nbrit,pars.tolres,mixres,
                      dx,ipvt,info,fail,verbose);
      exit when not fail;
      step := step/2.0;
      put(file,"Reduced step size to "); put(file,step,3); put_line(file,".");
      exit when (step < pars.minsize);
      DoblDobl_Rational_Approximations.Evaluate
        (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      Restore_Leading_Coefficients(homlead,hom.crc);
      Restore_Leading_Coefficients(abhlead,abh.crc);
    end loop;
  end Predictor_Corrector_Loop;

  procedure Predictor_Corrector_Loop
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                homlead,abhlead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                prd : in out QuadDobl_Predictor_Convolutions.Predictor;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                svh : in QuadDobl_Predictor_Convolutions.Link_to_SVD_Hessians;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                step : out quad_double;
                nbpole,nbhess,nbmaxm : in out natural32;
                fail : out boolean; verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Does one predictor-corrector step in quad double precision.

  -- ON ENTRY :
  --   file     to write the extra output to;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   homlead  leading coefficients for the circuits in hom;
  --   abhlead  leading coefficients for the circuits in abh;
  --   pars     values for the tolerances and parameters;
  --   prd      work space for the Newton-Fabry-Pade predictor;
  --   psv      work space vectors for the predictor,
  --            psv.sol contains a start solution;
  --   svh      work space for Hessian convolutions;
  --   nbpole   number of times the pole step was minimal;
  --   nbhess   number of times the Hessian step was minimal;
  --   nbmaxm   number of times the maximum step was minimal;
  --   verbose  flag for extra output.

  -- ON RETURN :
  --   psv.sol  the corrected solution;
  --   dx       last update to the solution;
  --   ipvt     pivoting information for the LU Newton steps;
  --   step     the step size;
  --   nbpole   updated number of times the pole step was minimal;
  --   nbhess   updated number of times the Hessian step was minimal;
  --   nbmaxm   updated number of times the maximum step was minimal;
  --   fail     true if the prescribed tolerance was not reached,
  --            false otherwise.

    use QuadDobl_Predictor_Convolutions;

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    deg : constant integer32 := numdeg + dendeg + 2;
    maxit : constant integer32 := 4;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    info,nbrit : integer32 := 0;
    mixres : quad_double;

  begin
    prd := Create(psv.sol,hom.neq,deg,numdeg,dendeg,SVD);
    Set_Lead_Coefficients(prd,psv.sol);
    svh.vals := (svh.vals'range => zero);
    SVD_Prediction(file,hom,abh,prd.svdata,svh,psv,maxit,pars.tolres,
      pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,pars.minsize,fail,step,
      nbpole,nbhess,nbmaxm,false,verbose);
    if verbose then
      if fail
       then put(file,"Predictor failed to reach tolerance");
       else put(file,"Predictor reached tolerance");
      end if;
      put(file,pars.alpha,3); put_line(file,".");
    end if;
    loop
      Step_Coefficient(hom,step);
      Step_Coefficient(abh,step);
      LU_Newton_Steps(file,hom,abh,psv,maxit,nbrit,pars.tolres,mixres,
                      dx,ipvt,info,fail,verbose);
      exit when not fail;
      step := step/2.0;
      put(file,"Reduced step size to "); put(file,step,3); put_line(file,".");
      exit when (step < pars.minsize);
      QuadDobl_Rational_Approximations.Evaluate
        (prd.svdata.numcff,prd.svdata.dencff,step,psv.sol);
      Restore_Leading_Coefficients(homlead,hom.crc);
      Restore_Leading_Coefficients(abhlead,abh.crc);
    end loop;
  end Predictor_Corrector_Loop;

  procedure Standard_Run
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

  -- DESCRIPTION :
  --   Does one predictor-corrector-shift step in double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters.

    use Standard_Complex_Solutions,Standard_Predictor_Convolutions;

    prd : Predictor;
    psv : Predictor_Vectors(hom.dim,hom.neq);
    hss : SVD_Hessians(hom.dim,hom.dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    solsptr : Solution_List := sols;
    ls : Link_to_Solution;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    step : double_float;
    fail : boolean;
    ans : character;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : Standard_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    loop
      ls := Head_Of(solsptr); psv.sol := ls.v;
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,prd,psv,svh,dx,ipvt,step,nbpole,nbhess,nbmaxm,fail,true);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Restore_Leading_Coefficients(homlead,hom.crc);
      Restore_Leading_Coefficients(abhlead,abh.crc);
    end loop;
    Clear(svh);
    Standard_Complex_VecVecs.Deep_Clear(homlead);
    Standard_Complex_VecVecs.Deep_Clear(abhlead);
  end Standard_Run;

  procedure DoblDobl_Run
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

  -- DESCRIPTION :
  --   Does one predictor-corrector-shift step in double double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters.

    use DoblDobl_Complex_Solutions,DoblDobl_Predictor_Convolutions;

    prd : Predictor;
    psv : Predictor_Vectors(hom.dim,hom.neq);
    hss : SVD_Hessians(hom.dim,hom.dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    solsptr : Solution_List := sols;
    ls : Link_to_Solution;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    step : double_double;
    fail : boolean;
    ans : character;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : DoblDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : DoblDobl_Complex_VecVecs.Link_to_VecVec;

  begin
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    loop
      ls := Head_Of(solsptr); psv.sol := ls.v;
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,prd,psv,svh,dx,ipvt,step,nbpole,nbhess,nbmaxm,fail,true);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Restore_Leading_Coefficients(homlead,hom.crc);
      Restore_Leading_Coefficients(abhlead,abh.crc);
    end loop;
    Clear(svh);
    DoblDobl_Complex_VecVecs.Deep_Clear(homlead);
    DoblDobl_Complex_VecVecs.Deep_Clear(abhlead);
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters ) is

  -- DESCRIPTION :
  --   Does one predictor-corrector-shift step in quad double precision.

  -- ON ENTRY :
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters.

    use QuadDobl_Complex_Solutions,QuadDobl_Predictor_Convolutions;

    prd : Predictor;
    psv : Predictor_Vectors(hom.dim,hom.neq);
    hss : SVD_Hessians(hom.dim,hom.dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    solsptr : Solution_List := sols;
    ls : Link_to_Solution;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    step : quad_double;
    fail : boolean;
    ans : character;
    ipvt : Standard_Integer_Vectors.Vector(1..hom.dim);
    dx : QuadDobl_Complex_Vectors.Vector(1..hom.dim);
    homlead,abhlead : QuadDobl_Complex_VecVecs.Link_to_VecVec;

  begin
    Allocate_Leading_Coefficients(hom.crc,homlead);
    Allocate_Leading_Coefficients(abh.crc,abhlead);
    Store_Leading_Coefficients(hom.crc,homlead);
    Store_Leading_Coefficients(abh.crc,abhlead);
    loop
      ls := Head_Of(solsptr); psv.sol := ls.v;
      Predictor_Corrector_Loop(standard_output,hom,abh,homlead,abhlead,
        pars,prd,psv,svh,dx,ipvt,step,nbpole,nbhess,nbmaxm,fail,true);
      if fail
       then put_line("Predictor-Corrector loop failed.");
       else put_line("Predictor-Corrector loop succeeded.");
      end if;
      solsptr := Tail_Of(solsptr);
      exit when Is_Null(solsptr);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Restore_Leading_Coefficients(homlead,hom.crc);
      Restore_Leading_Coefficients(abhlead,abh.crc);
    end loop;
    Clear(svh);
    QuadDobl_Complex_VecVecs.Deep_Clear(homlead);
    QuadDobl_Complex_VecVecs.Deep_Clear(abhlead);
  end QuadDobl_Run;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double precision.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;

  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    Standard_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    pars.gamma := Standard_Homotopy.Accessibility_Constant;
    Standard_Run(cnvhom,abshom,sols,pars);
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

    use DoblDobl_Complex_Numbers_cv;
  
  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    DoblDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    ddgamma := DoblDobl_Homotopy.Accessibility_Constant;
    pars.gamma := DoblDobl_Complex_to_Standard(ddgamma);
    DoblDobl_Run(cnvhom,abshom,sols,pars);
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

    use QuadDobl_Complex_Numbers_cv;

  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    QuadDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    qdgamma := QuadDobl_Homotopy.Accessibility_Constant;
    pars.gamma := QuadDobl_Complex_to_Standard(qdgamma);
    QuadDobl_Run(cnvhom,abshom,sols,pars);
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
