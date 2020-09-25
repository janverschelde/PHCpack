with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Matrices;
with DoblDobl_System_and_Solutions_io;
with DoblDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with DoblDobl_Newton_Convolutions;
with DoblDobl_Newton_Convolution_Steps;

package body Test_DoblDobl_NewtConvSteps is

  procedure DoblDobl_Run
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sol : in DoblDobl_Complex_Solutions.Link_to_Solution;
                deg,maxit : in integer32;
                scale,usesvd,useqrls,lurcond : in boolean ) is

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    neq : constant integer32 := p'last;
    dim : constant integer32 := sol.n;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : constant DoblDobl_Complex_VecVecs.VecVec(1..sol.n)
        := DoblDobl_Newton_Convolutions.Series_Coefficients(sol.v,deg);
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
    U : DoblDobl_Complex_Matrices.Matrix(1..neq,1..neq);
    V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    absdx,rcond : double_double;
    tol : constant double_float := 1.0E-14;
    fail : boolean;

  begin
    Add_Parameter_to_Constant(s);
    if useqrls or usesvd then
      dx := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := DoblDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
      if usesvd then
        DoblDobl_Newton_Convolution_Steps.SVD_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,wrk,scale);
        put("rcond : "); put(rcond,3); new_line;
      else
        DoblDobl_Newton_Convolution_Steps.QR_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      end if;
    elsif lurcond then
      DoblDobl_Newton_Convolution_Steps.LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale);
      put("rcond : "); put(rcond,3); new_line;
    else
      DoblDobl_Newton_Convolution_Steps.LU_Newton_Steps
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

  procedure Prompt_for_Parameters
              ( overdt : in boolean; maxit : out integer32;
                scale,usesvd,useqrls,lurcond,inlined : out boolean ) is

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
      useqrls := false; lurcond := false; inlined := false;
    elsif overdt then
      useqrls := true; lurcond := false; inlined := false;
    else
      put("Solve with least squares and QR ? (y/n) ");
      Ask_Yes_or_No(ans); useqrls := (ans = 'y');
      if useqrls then
        lurcond := false; inlined := false;
      else
        put("Estimate condition number ? (y/n) "); Ask_Yes_or_No(ans);
        lurcond := (ans = 'y');
        put("Inlined LU factorization ? (y/n) "); Ask_Yes_or_No(ans);
        inlined := (ans = 'y');
      end if;
    end if;
  end Prompt_for_Parameters;

  procedure DoblDobl_Run
              ( p : in DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim,deg : in integer32 ) is

    sol : constant DoblDobl_Complex_Solutions.Link_to_Solution
        := DoblDobl_Complex_Solutions.Head_Of(sols);
    maxit : integer32 := 0;
    scale,usesvd,useqrls,needrcond,inln : boolean;
    overdet : constant boolean := (p'last > dim);

  begin
    Prompt_for_Parameters(overdet,maxit,scale,usesvd,useqrls,needrcond,inln);
    DoblDobl_Run(p,sol,deg,maxit,scale,usesvd,useqrls,needrcond);
  end DoblDobl_Run;

  procedure DoblDobl_Test ( deg : in integer32 ) is

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

  procedure Main is

    deg : integer32 := 0;

  begin
    new_line;
    put_line("Newton on convolution circuits in double double precision.");
    new_line;
    put("Give the degree of the power series : "); get(deg); skip_line;
    DoblDobl_Test(deg);
  end Main;

end Test_DoblDobl_NewtConvSteps;
