with text_io;                            use text_io;
with duration_io;
with Ada.Calendar;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Vector_Splitters;
with Standard_Complex_Circuits;
with Standard_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Standard_Circuit_Makers;
with Convergence_Radius_Estimates;
with Multitasked_Series_Linearization;
with Multitasked_Newton_Convolutions;    use Multitasked_Newton_Convolutions;
with Multitasked_Hessian_Circuits;

procedure ts_mtprdcnv is

-- DESCRIPTION :
--   Development of the multitasked predictor on convolution circuits.

  procedure Newton_Fabry
              ( nbt,maxit : in integer32;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in Standard_Complex_VecVecs.VecVec;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in out Standard_Complex_VecVecs.VecVec;
                radius,raderr : out double_float;
                walltime : out Duration ) is

  -- DESCRIPTION :
  --   Runs the Newton-Fabry method for the convergence radius.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   maxit    the maximum number of iterations;
  --   hom      homotopy system with series coefficients;
  --   sol      solution to start the homotopy;
  --   ipvt     pivoting information to factor the lead matrix
  --            in the linearized solving;
  --   wrk      work space of range 1..nbt of vectors of
  --            the dimension, the number of variables.

  -- ON RETURN :
  --   radius   estimated radius of convergence;
  --   raderr   estimated error on the radius;
  --   walltime is the elapsed wall clock time.

    info,nbrit : integer32 := 0;
    tol : constant double_float := 1.0E-12;
    absdx : double_float;
    fail : boolean;
    z : Standard_Complex_Numbers.Complex_Number;
    tstart,tstop : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_LU_Newton_Steps
      (nbt,hom,sol,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk);
    Convergence_Radius_Estimates.Fabry(sol,z,radius,raderr,fail,0,false);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Newton_Fabry;

  procedure Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                svl : in out Standard_Complex_VecVecs.VecVec;
                walltime : out Duration ) is

  -- DESCRIPTION :
  --   Computes all singular values of the Hessian matrices.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        polynomials defined as a complex circuit system;
  --   x        leading coefficients of a series;
  --   vh       allocated space for all Hessian matrices;
  --   svl      allocate space for the singular values.

  -- ON RETURN :
  --   svl      the singular values of all Hessian matrices;
  --   walltime is the elapsed wall time.

    tstart,tstop : Ada.Calendar.Time;

    use Multitasked_Hessian_Circuits;
    use Ada.Calendar;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_Singular_Values(nbt,s,x,svl,false,false);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Singular_Values;

  procedure Prompt_for_Dimensions
              ( dim,deg,nbr,pwr : in out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user of the dimensions of the random input data.

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the power series : "); get(deg);
    put("Give the number of terms in each circuit : "); get(nbr);
    put("Give the largest power of the variables : "); get(pwr);
  end Prompt_for_Dimensions;           

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for dimensions and generates a random Newton homotopy.
  --   Runs a multitasked Newton-Fabry and singular value computation
  --   of all Hessian matrices.

    dim,deg,nbr,pwr : integer32 := 0;
    hom : Standard_Speelpenning_Convolutions.Link_to_System;
    sol : Standard_Complex_VecVecs.Link_to_VecVec;
    lnkcff,leadsol : Standard_Complex_Vectors.Link_to_Vector;
    crs : Standard_Complex_Circuits.Link_to_System;
    nbt,mxt : integer32 := 0;
    radius,raderr : double_float;
    walltime : duration;

  begin
    new_line;
    put_line("Testing the multitasked predictor convolutions ...");
    Prompt_for_Dimensions(dim,deg,nbr,pwr);
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,hom,sol);
    leadsol := new Standard_Complex_Vectors.Vector(1..dim);
    for k in 1..dim loop
      lnkcff := sol(k);        -- all coefficients of the series
      leadsol(k) := lnkcff(0); -- copy leading coefficient
    end loop;
   -- leadsol := sol(0); -- leading coefficients of the series sol
    crs := Standard_Circuit_Makers.Make_Complex_System(hom);
    put("Give the number of tasks : "); get(nbt);
    put("Give the maximum number of iterations : "); get(mxt);
    new_line;
    put_line("Running Newton-Fabry ...");
    declare
      pvt : Standard_Integer_Vectors.Vector(1..dim);
      wks : Standard_Complex_VecVecs.VecVec(1..nbt)
          := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    begin
      Newton_Fabry(nbt,mxt,hom,sol.all,pvt,wks,radius,raderr,walltime);
      put("estimated pole radius :"); put(radius,3);
      put(", with error :"); put(raderr,3); new_line;
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
      Standard_Complex_VecVecs.Clear(wks);
    end;
    new_line;
    put_line("Computing all singular values ...");
    declare
      svl : Standard_Complex_VecVecs.VecVec(0..crs.neq)
          := Standard_Vector_Splitters.Allocate(crs.neq,crs.dim+1,0,1);
    begin
      Singular_Values(nbt,crs,leadsol,svl,walltime);
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
    end;
    Standard_Speelpenning_Convolutions.Clear(hom);
    Standard_Complex_Circuits.Clear(crs);
  end Main;

begin
  Main;
end ts_mtprdcnv;
