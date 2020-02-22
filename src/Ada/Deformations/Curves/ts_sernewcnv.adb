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
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with Newton_Convolutions;
with Newton_Power_Convolutions;          use Newton_Power_Convolutions;

procedure ts_sernewcnv is

-- DESCRIPTION :
--   Procedure to develop the linearized Newton's method for power series,
--   on convolution circuits.

  procedure Standard_Run
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sol : in Standard_Complex_Solutions.Link_to_Solution;
                deg,maxit : in integer32;
                scale,usesvd,useqrls,lurcond : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double precision.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   sol      a solution;
  --   deg      degree of the power series;
  --   maxit    maximum number of iterations;
  --   scale    if scaling is needed;
  --   usesvd   for singular value decomposition;
  --   useqrls  for least squares after QR decomposition;
  --   lurcond  lu with condition number estimate.

    use Standard_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    neq : constant integer32 := p'last;
    dim : constant integer32 := sol.n;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : constant Standard_Complex_VecVecs.VecVec(1..sol.n)
        := Newton_Convolutions.Series_Coefficients(sol.v,deg);
    info,nbrit : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
    ewrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..dim);
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..neq);
    qraux,w1,w2,w3,w4,w5 : Standard_Complex_Vectors.Vector(1..neq);
    dx : Standard_Complex_VecVecs.VecVec(1..dim);
    xd : Standard_Complex_VecVecs.VecVec(0..deg);
    svl : Standard_Complex_Vectors.Vector(1..dim+1);
    U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    absdx,rcond : double_float;
    tol : constant double_float := 1.0E-14;
    fail : boolean;

  begin
    Add_Parameter_to_Constant(s);
    if useqrls or usesvd then
      dx := Standard_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := Standard_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
      if usesvd then
        SVD_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,wrk,scale);
        put("rcond :"); put(rcond,3); new_line;
      else
        QR_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      end if;
    elsif lurcond then
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale);
      put("rcond :"); put(rcond,3); new_line;
    else
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk,scale);
    end if;
    if fail then
      put("Failed to reach"); put(tol,3);
      put("  absdx : "); put(absdx,3); new_line;
    else
      put("Reached"); put(tol,3);
      put("  absdx : "); put(absdx,3); new_line;
    end if;
    put("after "); put(nbrit,1); put_line(" iterations.");
    Standard_Complex_Vectors.Clear(ewrk);
    Standard_Complex_Vectors.Clear(wrk);
  end Standard_Run;

  procedure DoblDobl_Run
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Link_to_Solution;
                deg,maxit : in integer32;
                scale,usesvd,useqrls,lurcond : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double double precision.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   sol      a solution;
  --   deg      degree of the power series;
  --   maxit    maximum number of iterations;
  --   scale    if scaling is needed;
  --   usesvd   for singular value decomposition;
  --   useqrls  for least squares after QR decomposition;
  --   lurcond  lu with condition number estimate.

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    neq : constant integer32 := p'last;
    dim : constant integer32 := sol.n;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : constant DoblDobl_Complex_VecVecs.VecVec(1..sol.n)
        := Newton_Convolutions.Series_Coefficients(sol.v,deg);
    info,nbrit : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
    ewrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..neq);
    qraux,w1,w2,w3,w4,w5 : DoblDobl_Complex_Vectors.Vector(1..neq);
    dx : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    xd : DoblDobl_Complex_VecVecs.VecVec(0..deg);
    svl : DoblDobl_Complex_Vectors.Vector(1..dim+1);
    U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    absdx,rcond : double_double;
    tol : constant double_float := 1.0E-14;
    fail : boolean;

  begin
    Add_Parameter_to_Constant(s);
    if useqrls or usesvd then
      dx := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := DoblDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
      if usesvd then
        SVD_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,wrk,scale);
        put("rcond : "); put(rcond,3); new_line;
      else
        QR_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      end if;
    elsif lurcond then
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale);
      put("rcond : "); put(rcond,3); new_line;
    else
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk,scale);
    end if;
    if fail then
      put("Failed to reach"); put(tol,3);
      put("  absdx : "); put(absdx,3); new_line;
    else
      put("Reached"); put(tol,3);
      put("  absdx : "); put(absdx,3); new_line;
    end if;
    put("after "); put(nbrit,1); put_line(" iterations.");
    DoblDobl_Complex_Vectors.Clear(ewrk);
    DoblDobl_Complex_Vectors.Clear(wrk);
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sol : in QuadDobl_Complex_Solutions.Link_to_Solution;
                deg,maxit : in integer32;
                scale,usesvd,useqrls,lurcond : in boolean ) is

  -- DESCRIPTION :
  --   Runs Newton's method in quad double precision.

  -- ON ENTRY :
  --   p        a polynomial system;
  --   sol      a solution;
  --   deg      degree of the power series;
  --   maxit    maximum number of iterations;
  --   scale    if scaling is needed;
  --   usesvd   for singular value decomposition;
  --   useqrls  for least squares after QR decomposition;
  --   lurcond  lu with condition number estimate.

    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    neq : constant integer32 := p'last;
    dim : constant integer32 := sol.n;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : constant QuadDobl_Complex_VecVecs.VecVec(1..sol.n)
        := Newton_Convolutions.Series_Coefficients(sol.v,deg);
    info,nbrit : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
    ewrk : QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector(1..dim);
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector(1..neq);
    qraux,w1,w2,w3,w4,w5 : QuadDobl_Complex_Vectors.Vector(1..neq);
    dx : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    xd : QuadDobl_Complex_VecVecs.VecVec(0..deg);
    svl : QuadDobl_Complex_Vectors.Vector(1..dim+1);
    U,V : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    absdx,rcond : quad_double;
    tol : constant double_float := 1.0E-14;
    fail : boolean;

  begin
    Add_Parameter_to_Constant(s);
    if useqrls or usesvd then
      dx := QuadDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := QuadDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
      if usesvd then
        SVD_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,wrk,scale);
        put("rcond : "); put(rcond,3); new_line;
      else
        QR_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      end if;
    elsif lurcond then
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale);
      put("rcond : "); put(rcond,3); new_line;
    else
      LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wrk,scale);
    end if;
    if fail then
      put("Failed to reach"); put(tol,3);
      put("  absdx : "); put(absdx,3); new_line;
    else
      put("Reached"); put(tol,3);
      put("  absdx : "); put(absdx,3); new_line;
    end if;
    put("after "); put(nbrit,1); put_line(" iterations.");
    QuadDobl_Complex_Vectors.Clear(ewrk);
    QuadDobl_Complex_Vectors.Clear(wrk);
  end QuadDobl_Run;

  procedure Prompt_for_Parameters
              ( overdt : in boolean; maxit : out integer32;
                scale,usesvd,useqrls,lurcond : out boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for the parameters of a run.

  -- ON ENTRY :
  --   overdt   if overdetermined or not.

  -- ON RETURN :
  --   maxit    maximum number of iterations;
  --   scale    if scaling is needed;
  --   usesvd   for singular value decomposition;
  --   useqrls  for least squares after QR decomposition;
  --   lurcond  lu with condition number estimate.

    ans : character;

  begin
    new_line;
    maxit := 0;
    put("Give the number of iterations : "); get(maxit);
    put("Apply scaling ? (y/n) "); Ask_Yes_or_No(ans);
    scale := (ans = 'y');
    put("Solve with singular value decomposition ? (y/n) ");
    Ask_Yes_or_No(ans);
    usesvd := (ans = 'y');
    if usesvd then
      useqrls := false; lurcond := false;
    elsif overdt then
      useqrls := true; lurcond := false;
    else
      put("Solve with least squares and QR ? (y/n) ");
      Ask_Yes_or_No(ans); useqrls := (ans = 'y');
      if useqrls then
        lurcond := false;
      else
        put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
        lurcond := (ans = 'y');
      end if;
    end if;
  end Prompt_for_Parameters;

  procedure Standard_Run
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double precision on the first solution
  --   in the list sols of the system p.

    sol : constant Standard_Complex_Solutions.Link_to_Solution
        := Standard_Complex_Solutions.Head_Of(sols);
    maxit : integer32 := 0;
    scale,usesvd,useqrls,needrcond : boolean := false;
    overdet : constant boolean := (p'last > dim);

  begin
    Prompt_for_Parameters(overdet,maxit,scale,usesvd,useqrls,needrcond);
    Standard_Run(p,sol,deg,maxit,scale,usesvd,useqrls,needrcond);
  end Standard_Run;

  procedure DoblDobl_Run
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method in double double precision on the first solution
  --   in the list sols of the system p.

    sol : constant DoblDobl_Complex_Solutions.Link_to_Solution
        := DoblDobl_Complex_Solutions.Head_Of(sols);
    maxit : integer32 := 0;
    scale,usesvd,useqrls,needrcond : boolean;
    overdet : constant boolean := (p'last > dim);

  begin
    Prompt_for_Parameters(overdet,maxit,scale,usesvd,useqrls,needrcond);
    DoblDobl_Run(p,sol,deg,maxit,scale,usesvd,useqrls,needrcond);
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( p : in QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim,deg : in integer32 ) is

  -- DESCRIPTION :
  --   Runs Newton's method in quad double precision on the first solution
  --   in the list sols of the system p.

    sol : constant QuadDobl_Complex_Solutions.Link_to_Solution
        := QuadDobl_Complex_Solutions.Head_Of(sols);
    maxit : integer32 := 0;
    scale,usesvd,useqrls,needrcond : boolean;
    overdet : constant boolean := (p'last > dim);

  begin
    Prompt_for_Parameters(overdet,maxit,scale,usesvd,useqrls,needrcond);
    QuadDobl_Run(p,sol,deg,maxit,scale,usesvd,useqrls,needrcond);
  end QuadDobl_Run;

  procedure Standard_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Given the degree of the power series, prompts the user
  --   for a polynomial system with a parameter,
  --   runs Newton's method in standard double precision.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    nbr : natural32;
    dim : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    nbr := Standard_Complex_Solutions.Length_Of(sols);
    if nbr = 0 then
      put_line("No solutions ?");
    else
      dim := Standard_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbr,1); put(" solutions in dimension ");
      put(dim,1); put_line(".");
      Standard_Run(lp,sols,dim,deg);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Given the degree of the power series, prompts the user
  --   for a polynomial system with a parameter,
  --   runs Newton's method in double double precision.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    nbr : natural32;
    dim : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    nbr := DoblDobl_Complex_Solutions.Length_Of(sols);
    if nbr = 0 then
      put_line("No solutions ?");
    else
      dim := DoblDobl_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbr,1); put(" solutions in dimension ");
      put(dim,1); put_line(".");
      DoblDobl_Run(lp,sols,dim,deg);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Given the degree of the power series, prompts the user
  --   for a polynomial system with a parameter,
  --   runs Newton's method in double double precision.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    nbr : natural32;
    dim : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    nbr := QuadDobl_Complex_Solutions.Length_Of(sols);
    if nbr = 0 then
      put_line("No solutions ?");
    else
      dim := QuadDobl_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbr,1); put(" solutions in dimension ");
      put(dim,1); put_line(".");
      QuadDobl_Run(lp,sols,dim,deg);
    end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree of the power series
  --   and then launches the test.

    deg : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Linearized Newton on power series with convolution circuits.");
    new_line;
    put("Give the degree of the power series : "); get(deg); skip_line;
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test(deg);
      when '1' => DoblDobl_Test(deg);
      when '2' => QuadDobl_Test(deg);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_sernewcnv;
