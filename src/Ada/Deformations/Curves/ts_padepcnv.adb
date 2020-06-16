with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with Standard_Vector_Splitters;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Complex_Circuits;
with Standard_Coefficient_Circuits;
with Standard_Coefficient_Convolutions;
with Standard_Convolution_Splitters;
with Standard_Circuit_Makers;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy_Convolutions_io;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Standard_Predictor_Convolutions;
with DoblDobl_Predictor_Convolutions;
with QuadDobl_Predictor_Convolutions;
With Test_Predictor_Convolutions;

procedure ts_padepcnv is

-- DESCRIPTION :
--   Development of the Pade predictor on convolution circuits.

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
    minstep : constant double_float := 1.0E-6;
    endt : constant double_float := 1.0;
    prd : Predictor;
    psv : Predictor_Vectors(dim,neq);
    hss : SVD_Hessians(dim,dim+1);
    svh : Link_to_SVD_Hessians;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    acct,step : double_float := 0.0;

  begin
    hss.vals := (hss.vals'range => Standard_Complex_Numbers.Create(0.0));
    svh := new SVD_Hessians'(hss);
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
        SVD_Prediction(standard_output,chom,abh,prd.svdata,svh,psv,maxit,
          tol,alpha,beta1,beta2,maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
      else
        LU_Prediction(standard_output,chom,abh,prd.ludata,svh,psv,maxit,
          tol,alpha,beta1,beta2,maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
      end if;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    Clear(prd); Clear(svh);
  end Standard_Run_Prediction;

  procedure Standard_Run_Prediction
              ( chom : in Standard_Coefficient_Convolutions.Link_to_System;
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
    minstep : constant double_float := 1.0E-6;
    endt : constant double_float := 1.0;
    prd : Predictor;
    psv : Predictor_Vectors(dim,neq);
    hss : SVD_Hessians(dim,dim+1);
    svh : Link_to_SVD_Hessians;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    acct,step : double_float := 0.0;
    cfh : Standard_Coefficient_Circuits.Link_to_System
        := Standard_Circuit_Makers.Make_Coefficient_System(chom);
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    xrv : constant Standard_Floating_Vectors.Vector(1..dim)
        := (1..dim => 0.0);
    xiv : constant Standard_Floating_Vectors.Vector(1..dim)
        := (1..dim => 0.0);
    xr : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(xrv);
    xi : Standard_Floating_Vectors.Link_to_Vector
       := new Standard_Floating_Vectors.Vector'(xiv);
    vh : Standard_Complex_VecMats.VecMat(1..neq)
       := Standard_Complex_Circuits.Allocate(neq,dim);
    svls : Standard_Complex_VecVecs.VecVec(0..neq)
         := Standard_Vector_Splitters.Allocate(neq,dim+1,0,1);

  begin
    hss.vals := (hss.vals'range => Standard_Complex_Numbers.Create(0.0));
    svh := new SVD_Hessians'(hss);
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
        SVD_Prediction(standard_output,chom,cfh,prd.svdata,svh,
          rx,ix,xr,xi,vh,svls,psv,maxit,tol,alpha,beta1,beta2,
          maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
      else
        LU_Prediction(standard_output,chom,cfh,prd.ludata,svh,
          rx,ix,xr,xi,vh,svls,psv,maxit,tol,alpha,beta1,beta2,
          maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
      end if;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    Clear(prd); Clear(svh);
    Standard_Floating_Vectors.Clear(xr);
    Standard_Floating_Vectors.Clear(xi);
    Standard_Complex_VecMats.Clear(vh);
    Standard_Complex_VecVecs.Clear(svls);
    Standard_Coefficient_Circuits.Clear(cfh);
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
    minstep : constant double_float := 1.0E-6;
    endt : constant double_float := 1.0;
    prd : Predictor;
    psv : Predictor_Vectors(dim,neq);
    hss : SVD_Hessians(dim,dim+1);
    svh : Link_to_SVD_Hessians;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    acct,step : double_double := create(0.0);

  begin
    hss.vals := (hss.vals'range => zero);
    svh := new SVD_Hessians'(hss);
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
        SVD_Prediction(standard_output,chom,abh,prd.svdata,svh,psv,maxit,
          tol,alpha,beta1,beta2,maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
      else
        LU_Prediction(standard_output,chom,abh,prd.ludata,svh,psv,maxit,
          tol,alpha,beta1,beta2,maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
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
    minstep : constant double_float := 1.0E-6;
    endt : constant double_float := 1.0;
    prd : Predictor;
    psv : Predictor_Vectors(dim,neq);
    hss : SVD_Hessians(dim,dim+1);
    svh : Link_to_SVD_Hessians;
    nbpole,nbhess,nbmaxm : natural32 := 0;
    acct,step : quad_double := create(0.0);

  begin
    hss.vals := (hss.vals'range => zero);
    svh := new SVD_Hessians'(hss);
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
        SVD_Prediction(standard_output,chom,abh,prd.svdata,svh,psv,maxit,
          tol,alpha,beta1,beta2,maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
      else
        LU_Prediction(standard_output,chom,abh,prd.ludata,svh,psv,maxit,
          tol,alpha,beta1,beta2,maxstep,minstep,endt,acct,fail,step,
          nbpole,nbhess,nbmaxm,otp,true);
      end if;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    Clear(prd); Clear(svh);
  end QuadDobl_Run_Prediction;

  procedure Standard_Test_Prediction ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    cffhom : Standard_Coefficient_Convolutions.Link_to_System;
    deg : constant integer32 := numdeg + dendeg + 2;
    idxpar : integer32;
    ans : character;

  begin
    Standard_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    new_line;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.Standard_Check_Solutions(cnvhom,sols);
    end if;
    new_line;
    put("Test on coefficient convolution circuits ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'n' then
      Standard_Run_Prediction(cnvhom,abshom,sols,deg,numdeg,dendeg);
    else
      cffhom := Standard_Convolution_Splitters.Split(cnvhom);
      Standard_Run_Prediction(cffhom,sols,deg,numdeg,dendeg);
    end if;
  end Standard_Test_Prediction;

  procedure DoblDobl_Test_Prediction ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    deg : constant integer32 := numdeg + dendeg + 2;
    idxpar : integer32;
    ans : character;

  begin
    DoblDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    new_line;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.DoblDobl_Check_Solutions(cnvhom,sols);
    end if;
    DoblDobl_Run_Prediction(cnvhom,abshom,sols,deg,numdeg,dendeg);
  end DoblDobl_Test_Prediction;

  procedure QuadDobl_Test_Prediction ( numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy to test the prediction.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    deg : constant integer32 := numdeg + dendeg + 2;
    idxpar : integer32;
    ans : character;

  begin
    QuadDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    put_line("The exponents in the circuits :");
    for k in cnvhom.crc'range loop
      Standard_Integer_VecVecs_io.put(cnvhom.crc(k).xps);
    end loop;
    new_line;
    put("Check all start solutions ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Test_Predictor_Convolutions.QuadDobl_Check_Solutions(cnvhom,sols);
    end if;
    QuadDobl_Run_Prediction(cnvhom,abshom,sols,deg,numdeg,dendeg);
  end QuadDobl_Test_Prediction;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision,
  --   the degrees of the Pade approximants, and then launches
  --   the test corresponding to the selected precision.

    precision : constant character := Prompt_for_Precision;
    numdeg,dendeg : integer32 := 0;

  begin
    new_line;
    put("Give the degree of the Pade numerator : "); get(numdeg);
    put("Give the degree of the Pade denominator : "); get(dendeg);
    case precision is
      when '0' => Standard_Test_Prediction(numdeg,dendeg);
      when '1' => DoblDobl_Test_Prediction(numdeg,dendeg);
      when '2' => QuadDobl_Test_Prediction(numdeg,dendeg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_padepcnv;
