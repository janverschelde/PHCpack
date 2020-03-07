with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Homotopy;
with Solution_Drops;
with Standard_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with Series_and_Homotopies;
with Test_Series_Predictors;
with Standard_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Newton_Convolutions;
with Newton_Power_Convolutions;          use Newton_Power_Convolutions;
with Convergence_Radius_Estimates;
with Standard_Rational_Approximations;

procedure ts_padepcnv is

-- DESCRIPTION :
--   Development of the Pade predictor on convolution circuits.

  procedure Standard_Prediction
              ( chom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in Standard_Complex_Solutions.Solution_List;
                deg,numdeg,dendeg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method on the homotopy chom,
  --   starting at the solutions in sols.

    use Standard_Complex_Solutions;
    use Standard_Rational_Approximations;
    use Standard_Speelpenning_Convolutions;
  
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    info,maxit,nbrit : integer32 := 0; 
    tol : constant double_float := 1.0E-12;
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(chom.crc'range);
    absdx : double_float;
    fail : boolean;
    ans : character;

  begin
    put("Give the maximum number of iterations : "); get(maxit);
    for k in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      declare
        scf : constant Standard_Complex_VecVecs.VecVec(1..ls.n)
            := Newton_Convolutions.Series_Coefficients(ls.v,deg);
        ipvt : Standard_Integer_Vectors.Vector(1..ls.n);
        z : Standard_Complex_Numbers.Complex_Number;
        r,err : double_float;
        rpiv : Standard_Integer_Vectors.Vector(1..dendeg);
        mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
        rhs : Standard_Complex_Vectors.Vector(1..dendeg);
        knum : constant Standard_Complex_Vectors.Vector(0..numdeg)
             := (0..numdeg => Standard_Complex_Numbers.Create(0.0));
        kden : constant Standard_Complex_Vectors.Vector(0..dendeg)
             := (0..dendeg => Standard_Complex_Numbers.Create(0.0));
        numcff : Standard_Complex_VecVecs.VecVec(1..ls.n);
        dencff : Standard_Complex_VecVecs.VecVec(1..ls.n);
        eva : Standard_Complex_Vectors.Vector(1..ls.n);
        res : Standard_Complex_Vectors.Vector(chom.crc'range);
      begin
        LU_Newton_Steps
          (standard_output,chom,scf,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wrk,false);
        Convergence_Radius_Estimates.Fabry(scf,z,r,err,fail);
        if not fail then
          put("z : "); put(z); 
          put("  error estimate :"); put(err,3); new_line;
          put("estimated radius :"); put(r,3); new_line;
        end if;
        for k in scf'range loop
          numcff(k) := new Standard_Complex_Vectors.Vector'(knum);
          dencff(k) := new Standard_Complex_Vectors.Vector'(kden);
        end loop;
        Pade_Vector(numdeg,dendeg,scf,numcff,dencff,mat,rhs,rpiv,info,false);
        Evaluate(numcff,dencff,r/2.0,eva);
        z := Standard_Complex_Numbers.Create(r/2.0);
        res := Eval(chom.crc,eva,z);
        put_line("Evaluation of the predicted solution : ");
        put_line(res);
      end;
      put("Continue to the next solution ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Complex_Vectors.Clear(wrk);
  end Standard_Prediction;

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
    Standard_Prediction(cnvhom,sols,deg,numdeg,dendeg);
  end Standard_Test_Prediction;

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

begin
  Standard_Main;
end ts_padepcnv;
