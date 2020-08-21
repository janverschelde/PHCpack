with text_io;                            use text_io;
with duration_io;
with Ada.Calendar;
with Communications_with_User;
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
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with Standard_Complex_Circuits;
with DoblDobl_Complex_Circuits;
with QuadDobl_Complex_Circuits;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Standard_Circuit_Makers;
with DoblDobl_Circuit_Makers;
with QuadDobl_Circuit_Makers;
with Convergence_Radius_Estimates;
with Multitasking;
with Multitasked_Series_Linearization;
with Multitasked_Newton_Convolutions;
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
  --   Runs the Newton-Fabry method for the convergence radius,
  --   in double precision, with multitasking.

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
    use Multitasked_Newton_Convolutions;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_LU_Newton_Steps
      (nbt,hom,sol,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk);
    Convergence_Radius_Estimates.Fabry(sol,z,radius,raderr,fail,0,false);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Newton_Fabry;

  procedure Newton_Fabry
              ( nbt,maxit : in integer32;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in out DoblDobl_Complex_VecVecs.VecVec;
                radius,raderr : out double_double;
                walltime : out Duration ) is

  -- DESCRIPTION :
  --   Runs the Newton-Fabry method for the convergence radius,
  --   in double double precision, with multitasking.

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
    tol : constant double_double := create(1.0E-24);
    absdx : double_double;
    fail : boolean;
    z : DoblDobl_Complex_Numbers.Complex_Number;
    tstart,tstop : Ada.Calendar.Time;

    use Ada.Calendar;
    use Multitasked_Newton_Convolutions;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_LU_Newton_Steps
      (nbt,hom,sol,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk);
    Convergence_Radius_Estimates.Fabry(sol,z,radius,raderr,fail,0,false);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Newton_Fabry;

  procedure Newton_Fabry
              ( nbt,maxit : in integer32;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in out Standard_Integer_Vectors.Vector;
                wrk : in out QuadDobl_Complex_VecVecs.VecVec;
                radius,raderr : out quad_double;
                walltime : out Duration ) is

  -- DESCRIPTION :
  --   Runs the Newton-Fabry method for the convergence radius,
  --   in quad double precision, with multitasking.

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
    tol : constant quad_double := create(1.0E-48);
    absdx : quad_double;
    fail : boolean;
    z : QuadDobl_Complex_Numbers.Complex_Number;
    tstart,tstop : Ada.Calendar.Time;

    use Ada.Calendar;
    use Multitasked_Newton_Convolutions;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_LU_Newton_Steps
      (nbt,hom,sol,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk);
    Convergence_Radius_Estimates.Fabry(sol,z,radius,raderr,fail,0,false);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Newton_Fabry;

  procedure Singular_Values_of_Hessians
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                svl : in out Standard_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out Standard_Complex_VecMats.VecMat;
                e : in out Standard_Complex_VecVecs.VecVec;
                yd : in out Standard_Complex_VecVecs.VecVec;
                eta : out double_float; walltime : out Duration ) is

  -- DESCRIPTION :
  --   Computes all singular values of the Hessian matrices,
  --   in double precision, with multitasking.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        polynomials defined as a complex circuit system;
  --   x        leading coefficients of a series;
  --   vh       allocated space for all Hessian matrices;
  --   svl      allocated space for the singular values.
  --   pwtdone  array of nbt flags to synchronize power table computation,
  --            must all be equal to false on entry;
  --   gradone  array of nbt flags to synchronize gradient computation,
  --            must all be equal to false on entry;
  --   A        vector of range 1..nbt of work space matrices;
  --   U        vector of range 1..nbt of work space matrices;
  --   V        vector of range 1..nbt of work space matrices;
  --   e        vector of range 1..nbt of vector space vectors;
  --   yd       vector of range 1..s.neq of work space for gradients.

  -- ON RETURN :
  --   svl      the singular values of all Hessian matrices;
  --   eta      the eta constant computed from the singular values;
  --   walltime is the elapsed wall time.

    tstart,tstop : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_Hessian_Circuits.Dynamic_Singular_Values
      (nbt,s,x,svl,pwtdone,gradone,A,U,V,e,yd,false);
    eta := Multitasked_Hessian_Circuits.Standard_Distance(svl);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Singular_Values_of_Hessians;

  procedure Singular_Values_of_Hessians
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                svl : in out DoblDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DoblDobl_Complex_VecMats.VecMat;
                e : in out DoblDobl_Complex_VecVecs.VecVec;
                yd : in out DoblDobl_Complex_VecVecs.VecVec;
                eta : out double_double; walltime : out Duration ) is

  -- DESCRIPTION :
  --   Computes all singular values of the Hessian matrices,
  --   in double double precision, with multitasking.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        polynomials defined as a complex circuit system;
  --   x        leading coefficients of a series;
  --   vh       allocated space for all Hessian matrices;
  --   svl      allocate space for the singular values;
  --   pwtdone  array of nbt flags to synchronize power table computation,
  --            must all be equal to false on entry;
  --   gradone  array of nbt flags to synchronize gradient computation,
  --            must all be equal to false on entry;
  --   A        vector of range 1..nbt of work space matrices;
  --   U        vector of range 1..nbt of work space matrices;
  --   V        vector of range 1..nbt of work space matrices;
  --   e        vector of range 1..nbt of vector space vectors;
  --   yd       vector of range 1..s.neq of work space for gradients.

  -- ON RETURN :
  --   svl      the singular values of all Hessian matrices;
  --   eta      the eta constant computed from the singular values;
  --   walltime is the elapsed wall time.

    tstart,tstop : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_Hessian_Circuits.Dynamic_Singular_Values
      (nbt,s,x,svl,pwtdone,gradone,A,U,V,e,yd,false);
    eta := Multitasked_Hessian_Circuits.DoblDobl_Distance(svl);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Singular_Values_of_Hessians;

  procedure Singular_Values_of_Hessians
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                svl : in out QuadDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out QuadDobl_Complex_VecMats.VecMat;
                e : in out QuadDobl_Complex_VecVecs.VecVec;
                yd : in out QuadDobl_Complex_VecVecs.VecVec;
                eta : out quad_double; walltime : out Duration ) is

  -- DESCRIPTION :
  --   Computes all singular values of the Hessian matrices,
  --   in quad double precision, with multitasking.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        polynomials defined as a complex circuit system;
  --   x        leading coefficients of a series;
  --   vh       allocated space for all Hessian matrices;
  --   svl      allocate space for the singular values.
  --   pwtdone  array of nbt flags to synchronize power table computation,
  --            must all be equal to false on entry;
  --   gradone  array of nbt flags to synchronize gradient computation,
  --            must all be equal to false on entry;
  --   A        vector of range 1..nbt of work space matrices;
  --   U        vector of range 1..nbt of work space matrices;
  --   V        vector of range 1..nbt of work space matrices;
  --   e        vector of range 1..nbt of vector space vectors;
  --   yd       vector of range 1..s.neq of work space for gradients.

  -- ON RETURN :
  --   svl      the singular values of all Hessian matrices;
  --   eta      the eta constant computed from the singular values;
  --   walltime is the elapsed wall time.

    tstart,tstop : Ada.Calendar.Time;

    use Ada.Calendar;

  begin
    tstart := Ada.Calendar.Clock;
    Multitasked_Hessian_Circuits.Dynamic_Singular_Values
      (nbt,s,x,svl,pwtdone,gradone,A,U,V,e,yd,false);
    eta := Multitasked_Hessian_Circuits.QuadDobl_Distance(svl);
    tstop := Ada.Calendar.Clock;
    walltime := tstop - tstart;
  end Singular_Values_of_Hessians;

  procedure Standard_Main ( nbt1,nbt2,dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts for dimensions and generates a random Newton homotopy.
  --   Runs a multitasked Newton-Fabry and singular value computation
  --   of all Hessian matrices in double precision.

  -- ON ENTRY :
  --   nbt1     the number of tasks for Newton-Fabry;
  --   nbt2     the number of tasks for SVD of Hessians;
  --   dim      the dimension is the number of equations and variables;
  --   deg      degree of the power series;
  --   nbr      number of terms per equations;
  --   pwr      largest power of the variables.

    hom : Standard_Speelpenning_Convolutions.Link_to_System;
    sol : Standard_Complex_VecVecs.Link_to_VecVec;
    lnkcff,leadsol : Standard_Complex_Vectors.Link_to_Vector;
    crs : Standard_Complex_Circuits.Link_to_System;
    mxt : integer32 := 0;
    radius,raderr,eta : double_float;
    walltime : duration;

  begin
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,hom,sol);
    leadsol := new Standard_Complex_Vectors.Vector(1..dim);
    for k in 1..dim loop
      lnkcff := sol(k);        -- all coefficients of the series
      leadsol(k) := lnkcff(0); -- copy leading coefficient
    end loop;
    crs := Standard_Circuit_Makers.Make_Complex_System(hom);
    put("Give the maximum number of iterations : "); get(mxt);
    new_line;
    put("Running Newton-Fabry with "); put(nbt1,1); put_line(" tasks ...");
    declare
      pvt : Standard_Integer_Vectors.Vector(1..dim);
      wks : Standard_Complex_VecVecs.VecVec(1..nbt1)
          := Multitasked_Series_Linearization.Allocate_Work_Space(nbt1,dim);
    begin
      Newton_Fabry(nbt1,mxt,hom,sol.all,pvt,wks,radius,raderr,walltime);
      put("estimated pole radius :"); put(radius,3);
      put(", with error :"); put(raderr,3); new_line;
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
      Standard_Complex_VecVecs.Clear(wks);
    end;
    new_line;
    put("Computing all singular values with "); put(nbt2,1);
    put_line(" tasks ...");
    declare
      svl : Standard_Complex_VecVecs.VecVec(0..crs.neq)
          := Standard_Vector_Splitters.Allocate(crs.neq,crs.dim+1,0,1);
      pwtdone : Multitasking.boolean_array(1..nbt2) := (1..nbt2 => false);
      gradone : Multitasking.boolean_array(1..nbt2) := (1..nbt2 => false);
      A,U,V : Standard_Complex_VecMats.VecMat(1..nbt2); -- space for SVD
      e : Standard_Complex_VecVecs.VecVec(1..nbt2);
      yd : Standard_Complex_VecVecs.VecVec(1..crs.neq); -- gradient space
    begin
      A := Standard_Complex_Circuits.Allocate(nbt2,crs.dim);
      U := Standard_Complex_Circuits.Allocate(nbt2,crs.dim);
      V := Standard_Complex_Circuits.Allocate(nbt2,crs.dim);
      e := Standard_Vector_Splitters.Allocate(nbt2,crs.dim,1,1);
      yd := Standard_Vector_Splitters.Allocate(crs.neq,crs.dim,1,0);
      Singular_Values_of_Hessians
        (nbt2,crs,leadsol,svl,pwtdone,gradone,A,U,V,e,yd,eta,walltime);
      put("eta :"); put(eta,3); new_line;
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
      Standard_Complex_VecMats.Clear(A);
      Standard_Complex_VecMats.Clear(U);
      Standard_Complex_VecMats.Clear(V);
      Standard_Complex_VecVecs.Clear(e);
      Standard_Complex_VecVecs.Clear(yd);
    end;
    Standard_Speelpenning_Convolutions.Clear(hom);
    Standard_Complex_Circuits.Clear(crs);
  end Standard_Main;

  procedure DoblDobl_Main ( nbt1,nbt2,dim,deg,nbr,pwr : integer32 ) is

  -- DESCRIPTION :
  --   Prompts for dimensions and generates a random Newton homotopy.
  --   Runs a multitasked Newton-Fabry and singular value computation
  --   of all Hessian matrices in double double precision.

  -- ON ENTRY :
  --   nbt1     the number of tasks for Newton-Fabry;
  --   nbt2     the number of tasks for SVD of Hessians;
  --   dim      the dimension is the number of equations and variables;
  --   deg      degree of the power series;
  --   nbr      number of terms per equations;
  --   pwr      largest power of the variables.

    hom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    sol : DoblDobl_Complex_VecVecs.Link_to_VecVec;
    lnkcff,leadsol : DoblDobl_Complex_Vectors.Link_to_Vector;
    crs : DoblDobl_Complex_Circuits.Link_to_System;
    mxt : integer32 := 0;
    radius,raderr,eta : double_double;
    walltime : duration;

  begin
    DoblDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,hom,sol);
    leadsol := new DoblDobl_Complex_Vectors.Vector(1..dim);
    for k in 1..dim loop
      lnkcff := sol(k);        -- all coefficients of the series
      leadsol(k) := lnkcff(0); -- copy leading coefficient
    end loop;
    crs := DoblDobl_Circuit_Makers.Make_Complex_System(hom);
    put("Give the maximum number of iterations : "); get(mxt);
    new_line;
    put("Running Newton-Fabry with "); put(nbt1,1); put_line(" tasks ...");
    declare
      pvt : Standard_Integer_Vectors.Vector(1..dim);
      wks : DoblDobl_Complex_VecVecs.VecVec(1..nbt1)
          := Multitasked_Series_Linearization.Allocate_Work_Space(nbt1,dim);
    begin
      Newton_Fabry(nbt1,mxt,hom,sol.all,pvt,wks,radius,raderr,walltime);
      put("estimated pole radius : "); put(radius,3);
      put(", with error : "); put(raderr,3); new_line;
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
      DoblDobl_Complex_VecVecs.Clear(wks);
    end;
    new_line;
    put("Computing all singular values with "); put(nbt2,1);
    put_line(" tasks ...");
    declare
      svl : DoblDobl_Complex_VecVecs.VecVec(0..crs.neq)
          := DoblDobl_Vector_Splitters.Allocate(crs.neq,crs.dim+1,0,1);
      pwtdone : Multitasking.boolean_array(1..nbt2) := (1..nbt2 => false);
      gradone : Multitasking.boolean_array(1..nbt2) := (1..nbt2 => false);
      A,U,V : DoblDobl_Complex_VecMats.VecMat(1..nbt2); -- space for SVD
      e : DoblDobl_Complex_VecVecs.VecVec(1..nbt2);
      yd : DoblDobl_Complex_VecVecs.VecVec(1..crs.neq); -- gradient space
    begin
      A := DoblDobl_Complex_Circuits.Allocate(nbt2,crs.dim);
      U := DoblDobl_Complex_Circuits.Allocate(nbt2,crs.dim);
      V := DoblDobl_Complex_Circuits.Allocate(nbt2,crs.dim);
      e := DoblDobl_Vector_Splitters.Allocate(nbt2,crs.dim,1,1);
      yd := DoblDobl_Vector_Splitters.Allocate(crs.neq,crs.dim,1,0);
      Singular_Values_of_Hessians
        (nbt2,crs,leadsol,svl,pwtdone,gradone,A,U,V,e,yd,eta,walltime);
      put("eta : "); put(eta,3); new_line;
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
      DoblDobl_Complex_VecMats.Clear(A);
      DoblDobl_Complex_VecMats.Clear(U);
      DoblDobl_Complex_VecMats.Clear(V);
      DoblDobl_Complex_VecVecs.Clear(e);
      DoblDobl_Complex_VecVecs.Clear(yd);
    end;
    DoblDobl_Speelpenning_Convolutions.Clear(hom);
    DoblDobl_Complex_Circuits.Clear(crs);
  end DoblDobl_Main;

  procedure QuadDobl_Main ( nbt1,nbt2,dim,deg,nbr,pwr : integer32 ) is

  -- DESCRIPTION :
  --   Prompts for dimensions and generates a random Newton homotopy.
  --   Runs a multitasked Newton-Fabry and singular value computation
  --   of all Hessian matrices in quad double precision.

  -- ON ENTRY :
  --   nbt1     the number of tasks for Newton-Fabry;
  --   nbt2     the number of tasks for SVD of Hessians;
  --   dim      the dimension is the number of equations and variables;
  --   deg      degree of the power series;
  --   nbr      number of terms per equations;
  --   pwr      largest power of the variables.

    hom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    sol : QuadDobl_Complex_VecVecs.Link_to_VecVec;
    lnkcff,leadsol : QuadDobl_Complex_Vectors.Link_to_Vector;
    crs : QuadDobl_Complex_Circuits.Link_to_System;
    mxt : integer32 := 0;
    radius,raderr,eta : quad_double;
    walltime : duration;

  begin
    QuadDobl_Random_Newton_Homotopy(dim,deg,nbr,pwr,hom,sol);
    leadsol := new QuadDobl_Complex_Vectors.Vector(1..dim);
    for k in 1..dim loop
      lnkcff := sol(k);        -- all coefficients of the series
      leadsol(k) := lnkcff(0); -- copy leading coefficient
    end loop;
    crs := QuadDobl_Circuit_Makers.Make_Complex_System(hom);
    put("Give the maximum number of iterations : "); get(mxt);
    new_line;
    put("Running Newton-Fabry with "); put(nbt1,1); put_line(" tasks...");
    declare
      pvt : Standard_Integer_Vectors.Vector(1..dim);
      wks : QuadDobl_Complex_VecVecs.VecVec(1..nbt1)
          := Multitasked_Series_Linearization.Allocate_Work_Space(nbt1,dim);
    begin
      Newton_Fabry(nbt1,mxt,hom,sol.all,pvt,wks,radius,raderr,walltime);
      put("estimated pole radius : "); put(radius,3);
      put(", with error : "); put(raderr,3); new_line;
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
      QuadDobl_Complex_VecVecs.Clear(wks);
    end;
    new_line;
    put("Computing all singular values with "); put(nbt2,1);
    put_line(" tasks ...");
    declare
      svl : QuadDobl_Complex_VecVecs.VecVec(0..crs.neq)
          := QuadDobl_Vector_Splitters.Allocate(crs.neq,crs.dim+1,0,1);
      pwtdone : Multitasking.boolean_array(1..nbt2) := (1..nbt2 => false);
      gradone : Multitasking.boolean_array(1..nbt2) := (1..nbt2 => false);
      A,U,V : QuadDobl_Complex_VecMats.VecMat(1..nbt2); -- space for SVD
      e : QuadDobl_Complex_VecVecs.VecVec(1..nbt2);
      yd : QuadDobl_Complex_VecVecs.VecVec(1..crs.neq); -- gradient space
    begin
      A := QuadDobl_Complex_Circuits.Allocate(nbt2,crs.dim);
      U := QuadDobl_Complex_Circuits.Allocate(nbt2,crs.dim);
      V := QuadDobl_Complex_Circuits.Allocate(nbt2,crs.dim);
      e := QuadDobl_Vector_Splitters.Allocate(nbt2,crs.dim,1,1);
      yd := QuadDobl_Vector_Splitters.Allocate(crs.neq,crs.dim,1,0);
      Singular_Values_of_Hessians
        (nbt2,crs,leadsol,svl,pwtdone,gradone,A,U,V,e,yd,eta,walltime);
      put("eta : "); put(eta,3); new_line;
      put("Wall clock time : "); duration_io.put(walltime,3,3); new_line;
      QuadDobl_Complex_VecMats.Clear(A);
      QuadDobl_Complex_VecMats.Clear(U);
      QuadDobl_Complex_VecMats.Clear(V);
      QuadDobl_Complex_VecVecs.Clear(e);
      QuadDobl_Complex_VecVecs.Clear(yd);
    end;
    QuadDobl_Speelpenning_Convolutions.Clear(hom);
    QuadDobl_Complex_Circuits.Clear(crs);
  end QuadDobl_Main;

  procedure Prompt_for_Dimensions
              ( dim,deg,nbr,pwr : in out integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user of the dimensions of the random input data.

  -- ON RETURN :
  --   dim      the dimension is the number of equations and variables;
  --   deg      degree of the power series;
  --   nbr      number of terms per equations;
  --   pwr      largest power of the variables.

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the power series : "); get(deg);
    put("Give the number of terms in each circuit : "); get(nbr);
    put("Give the largest power of the variables : "); get(pwr);
  end Prompt_for_Dimensions;           

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision, the dimensions,
  --   and then launches the test.

    prc : character;
    dim,deg,nbr,pwr,nbt1,nbt2 : integer32 := 0;

  begin
    new_line;
    put_line("Testing the multitasked predictor convolutions ...");
    prc := Communications_with_User.Prompt_for_Precision;
    Prompt_for_Dimensions(dim,deg,nbr,pwr);
    new_line;
    put("Give the number of tasks for Newton-Fabry : "); get(nbt1);
    put("Give the number of tasks for SVD of Hessians : "); get(nbt2);
    case prc is
      when '0' => Standard_Main(nbt1,nbt2,dim,deg,nbr,pwr);
      when '1' => DoblDobl_Main(nbt1,nbt2,dim,deg,nbr,pwr);
      when '2' => QuadDobl_Main(nbt1,nbt2,dim,deg,nbr,pwr);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtprdcnv;
