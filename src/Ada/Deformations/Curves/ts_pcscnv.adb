with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
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
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;

procedure ts_pcscnv is

-- DESCRIPTION :
--   Development of one predictor-corrector-shift step with
--   a homotopy system of convolution circuits.

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

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    deg : constant integer32 := numdeg + dendeg + 2;
    maxit : constant integer32 := 4;
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

  begin
    loop
      ls := Head_Of(solsptr);
      prd := Create(ls.v,hom.neq,deg,numdeg,dendeg,SVD);
      Set_Lead_Coefficients(prd,ls.v);
      hss.vals := (hss.vals'range => Standard_Complex_Numbers.Create(0.0));
      SVD_Prediction(standard_output,hom,abh,prd.svdata,svh,psv,maxit,
        pars.tolres,pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,fail,step,
        nbpole,nbhess,nbmaxm,false,true);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      solsptr := Tail_Of(solsptr);
    end loop;
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

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    deg : constant integer32 := numdeg + dendeg + 2;
    maxit : constant integer32 := 4;
    prd : Predictor;
    psv : Predictor_Vectors(hom.dim,hom.neq);
    hss : SVD_Hessians(hom.dim,hom.dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    solsptr : Solution_List := sols;
    ls : Link_to_Solution;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    step : double_double;
    fail : boolean;
    zero : constant DoblDobl_Complex_Numbers.Complex_Number
         := DoblDobl_Complex_Numbers.Create(integer(0));
    ans : character;

  begin
    loop
      ls := Head_Of(solsptr);
      prd := Create(ls.v,hom.neq,deg,numdeg,dendeg,SVD);
      Set_Lead_Coefficients(prd,ls.v);
      hss.vals := (hss.vals'range => zero);
      SVD_Prediction(standard_output,hom,abh,prd.svdata,svh,psv,maxit,
        pars.tolres,pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,fail,step,
        nbpole,nbhess,nbmaxm,false,true);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      solsptr := Tail_Of(solsptr);
    end loop;
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

    numdeg : constant integer32 := integer32(pars.numdeg);
    dendeg : constant integer32 := integer32(pars.dendeg);
    deg : constant integer32 := integer32(pars.numdeg + pars.dendeg + 2);
    maxit : constant integer32 := 4;
    prd : Predictor;
    psv : Predictor_Vectors(hom.dim,hom.neq);
    hss : SVD_Hessians(hom.dim,hom.dim+1);
    svh : Link_to_SVD_Hessians := new SVD_Hessians'(hss);
    solsptr : Solution_List := sols;
    ls : Link_to_Solution;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    step : quad_double;
    fail : boolean;
    zero : constant QuadDobl_Complex_Numbers.Complex_Number
         := QuadDobl_Complex_Numbers.Create(integer(0));
    ans : character;

  begin
    loop
      ls := Head_Of(solsptr);
      prd := Create(ls.v,hom.neq,deg,numdeg,dendeg,SVD);
      Set_Lead_Coefficients(prd,ls.v);
      hss.vals := (hss.vals'range => zero);
      SVD_Prediction(standard_output,hom,abh,prd.svdata,svh,psv,maxit,
        pars.tolres,pars.alpha,pars.pbeta,pars.cbeta,pars.maxsize,fail,step,
        nbpole,nbhess,nbmaxm,false,true);
      new_line;
      put("Move to the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      solsptr := Tail_Of(solsptr);
    end loop;
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
