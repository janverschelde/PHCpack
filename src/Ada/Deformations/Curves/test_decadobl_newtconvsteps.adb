with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Standard_Integer_Vectors;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_VecVecs;
with DecaDobl_Complex_Matrices;
with DecaDobl_System_and_Solutions_io;
with DecaDobl_Speelpenning_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with DecaDobl_Newton_Convolutions;
with DecaDobl_Newton_Convolution_Steps;

package body Test_DecaDobl_NewtConvSteps is

  procedure DecaDobl_Run
              ( p : in DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sol : in DecaDobl_Complex_Solutions.Link_to_Solution;
                deg,maxit : in integer32;
                scale,usesvd,useqrls,lurcond : in boolean ) is

    use DecaDobl_Speelpenning_Convolutions;

   -- c : constant Circuits(p'range)
   --   := Make_Convolution_Circuits(p.all,natural32(deg));
    c : Circuits(p'range);
    s : Link_to_System;
    neq : constant integer32 := p'last;
    dim : constant integer32 := sol.n;
   -- s : constant Link_to_System := Create(c,dim,deg);
    scf : constant DecaDobl_Complex_VecVecs.VecVec(1..sol.n)
        := DecaDobl_Newton_Convolutions.Series_Coefficients(sol.v,deg);
    info,nbrit : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
    ewrk : DecaDobl_Complex_Vectors.Link_to_Vector
        := new DecaDobl_Complex_Vectors.Vector(1..dim);
    wrk : DecaDobl_Complex_Vectors.Link_to_Vector
        := new DecaDobl_Complex_Vectors.Vector(1..neq);
    qraux,w1,w2,w3,w4,w5 : DecaDobl_Complex_Vectors.Vector(1..neq);
    dx : DecaDobl_Complex_VecVecs.VecVec(1..dim);
    xd : DecaDobl_Complex_VecVecs.VecVec(0..deg);
    svl : DecaDobl_Complex_Vectors.Vector(1..dim+1);
    U : DecaDobl_Complex_Matrices.Matrix(1..neq,1..neq);
    V : DecaDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    absdx,rcond : deca_double;
    tol : constant double_float := 1.0E-14;
    fail : boolean;

  begin
    c := Make_Convolution_Circuits(p.all,natural32(deg));
    s := Create(c,dim,deg);
    Add_Parameter_to_Constant(s);
    if useqrls or usesvd then
      dx := DecaDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := DecaDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
      if usesvd then
        DecaDobl_Newton_Convolution_Steps.SVD_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,wrk,scale);
        put("rcond : "); put(rcond,3); new_line;
      else
        DecaDobl_Newton_Convolution_Steps.QR_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      end if;
    elsif lurcond then
      DecaDobl_Newton_Convolution_Steps.LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale);
      put("rcond : "); put(rcond,3); new_line;
    else
      DecaDobl_Newton_Convolution_Steps.LU_Newton_Steps
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
    DecaDobl_Complex_Vectors.Clear(ewrk);
    DecaDobl_Complex_Vectors.Clear(wrk);
  end DecaDobl_Run;

  procedure Prompt_for_Parameters
              ( overdt : in boolean; maxit : out integer32;
                scale,usesvd,useqrls,lurcond : out boolean ) is

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

  procedure DecaDobl_Run
              ( p : in DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in DecaDobl_Complex_Solutions.Solution_List;
                dim,deg : in integer32 ) is

    sol : constant DecaDobl_Complex_Solutions.Link_to_Solution
        := DecaDobl_Complex_Solutions.Head_Of(sols);
    maxit : integer32 := 0;
    scale,usesvd,useqrls,needrcond : boolean;
    overdet : constant boolean := (p'last > dim);

  begin
    Prompt_for_Parameters(overdet,maxit,scale,usesvd,useqrls,needrcond);
    DecaDobl_Run(p,sol,deg,maxit,scale,usesvd,useqrls,needrcond);
  end DecaDobl_Run;

  procedure DecaDobl_Test ( deg : in integer32 ) is

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DecaDobl_Complex_Solutions.Solution_List;
    nbr : natural32;
    dim : integer32;

  begin
    new_line;
    put_line("Reading a polynomial system with solutions ...");
    DecaDobl_System_and_Solutions_io.get(lp,sols);
    nbr := DecaDobl_Complex_Solutions.Length_Of(sols);
    if nbr = 0 then
      put_line("No solutions ?");
    else
      dim := DecaDobl_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbr,1); put(" solutions in dimension ");
      put(dim,1); put_line(".");
      DecaDobl_Run(lp,sols,dim,deg);
    end if;
  end DecaDobl_Test;

  procedure Main is

    deg : integer32 := 0;

  begin
    new_line;
    put_line("Newton on convolution circuits in deca double precision.");
    new_line;
    put("Give the degree of the power series : "); get(deg); skip_line;
    DecaDobl_Test(deg);
  end Main;

end Test_DecaDobl_NewtConvSteps;
