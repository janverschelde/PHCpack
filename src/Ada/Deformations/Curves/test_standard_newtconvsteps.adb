with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;
with Standard_System_and_Solutions_io;
with Standard_Speelpenning_Convolutions;
with Standard_Vector_Splitters;
with Standard_Convolution_Splitters;
with Standard_Coefficient_Convolutions;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with Standard_Newton_Convolutions;
with Standard_Newton_Convolution_Steps;
with Staggered_Newton_Convolutions;

package body Test_Standard_NewtConvSteps is

  procedure Standard_Coefficient_Run
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sol : in Standard_Complex_Solutions.Link_to_Solution;
                deg,maxit : in integer32;
                scale,usesvd,useqrls,lurcond : in boolean;
                stagdeg,inlined,indexed : in boolean ) is

    c : constant Standard_Speelpenning_Convolutions.Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    neq : constant integer32 := p'last;
    dim : constant integer32 := sol.n;
    s : constant Standard_Speelpenning_Convolutions.Link_to_System
      := Standard_Speelpenning_Convolutions.Create(c,dim,deg);
    cs : Standard_Coefficient_Convolutions.Link_to_System;
    scf : constant Standard_Complex_VecVecs.VecVec(1..sol.n)
        := Standard_Newton_Convolutions.Series_Coefficients(sol.v,deg);
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    info,nbrit,idxtoldx : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..sol.n);
    ewrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..dim);
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..neq);
    qraux,w1,w2,w3,w4,w5 : Standard_Complex_Vectors.Vector(1..neq);
    dx : Standard_Complex_VecVecs.VecVec(1..dim);
    xd : Standard_Complex_VecVecs.VecVec(0..deg);
    svl : Standard_Complex_Vectors.Vector(1..dim+1);
    U : Standard_Complex_Matrices.Matrix(1..neq,1..neq);
    V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    absdx,rcond : double_float;
    tol : constant double_float := 1.0E-12;
    fail : boolean;
    rc,ic,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;

  begin
    Add_Parameter_to_Constant(s);
    cs := Standard_Convolution_Splitters.Split(s);
    put_line("The coefficients of the circuits : ");
    for k in s.crc'range loop
      Standard_Complex_VecVecs_io.put_line(s.crc(k).cff);
      put("The constant of circuit "); put(k,1); put_line(" :");
      Standard_Complex_Vectors_io.put_line(s.crc(k).cst);
    end loop;
    if inlined then -- allocate work space vectors
      Standard_Floating_VecVecVecs.Allocate(rv,1,deg,1,dim,1,dim);
      Standard_Floating_VecVecVecs.Allocate(iv,1,deg,1,dim,1,dim);
      rc := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
      ic := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
      rb := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
      ib := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
      ry := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
      iy := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    end if;
    if useqrls or usesvd then
      dx := Standard_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := Standard_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
      if usesvd then
        if stagdeg then
          Staggered_Newton_Convolutions.SVD_Newton_Steps
            (standard_output,cs,scf,dx,xd,rx,ix,maxit,nbrit,tol,absdx,fail,
             svl,U,V,info,rcond,ewrk,wrk,scale);
        else
          Standard_Newton_Convolution_Steps.SVD_Newton_Steps
            (standard_output,cs,scf,dx,xd,rx,ix,maxit,nbrit,tol,absdx,fail,
             svl,U,V,info,rcond,ewrk,wrk,scale);
        end if;
        put("rcond :"); put(rcond,3); new_line;
      else
        if stagdeg then
          Staggered_Newton_Convolutions.QR_Newton_Steps
            (standard_output,cs,scf,dx,xd,rx,ix,maxit,nbrit,tol,absdx,fail,
             qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
        else
          Standard_Newton_Convolution_Steps.QR_Newton_Steps
            (standard_output,cs,scf,dx,xd,rx,ix,maxit,nbrit,tol,absdx,fail,
             qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
        end if;
      end if;
    elsif lurcond then
      if stagdeg then
        if inlined then
          if indexed then
            Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
              (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,idxtoldx,absdx,
               fail,rcond,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
          else
            Staggered_Newton_Convolutions.Inlined_LU_Newton_Steps
              (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,absdx,fail,rcond,
               ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
          end if;
        else
          Staggered_Newton_Convolutions.LU_Newton_Steps
            (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,absdx,fail,rcond,
             ipvt,wrk,scale);
        end if;
      else
        Standard_Newton_Convolution_Steps.LU_Newton_Steps
          (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,absdx,fail,rcond,
           ipvt,wrk,scale);
      end if;
      put("rcond :"); put(rcond,3); new_line;
    else
      if stagdeg then
        if inlined then
          if indexed then
            Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
              (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,idxtoldx,absdx,
               fail,info,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
          else
            Staggered_Newton_Convolutions.Inlined_LU_Newton_Steps
              (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,absdx,fail,
               info,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
          end if;
        else
          Staggered_Newton_Convolutions.LU_Newton_Steps
            (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,absdx,fail,
             info,ipvt,wrk,scale);
        end if;
      else
        Standard_Newton_Convolution_Steps.LU_Newton_Steps
          (standard_output,cs,scf,rx,ix,maxit,nbrit,tol,absdx,fail,
           info,ipvt,wrk,scale);
      end if;
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
  end Standard_Coefficient_Run;

  procedure Standard_Run
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sol : in Standard_Complex_Solutions.Link_to_Solution;
                deg,maxit : in integer32;
                scale,usesvd,useqrls,lurcond : in boolean ) is

    use Standard_Speelpenning_Convolutions;

    c : constant Circuits(p'range)
      := Make_Convolution_Circuits(p.all,natural32(deg));
    neq : constant integer32 := p'last;
    dim : constant integer32 := sol.n;
    s : constant Link_to_System := Create(c,dim,deg);
    scf : constant Standard_Complex_VecVecs.VecVec(1..sol.n)
        := Standard_Newton_Convolutions.Series_Coefficients(sol.v,deg);
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
    U : Standard_Complex_Matrices.Matrix(1..neq,1..neq);
    V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    absdx,rcond : double_float;
    tol : constant double_float := 1.0E-14;
    fail : boolean;

  begin
    Add_Parameter_to_Constant(s);
    put_line("The coefficients of the circuits : ");
    for k in s.crc'range loop
      Standard_Complex_VecVecs_io.put_line(s.crc(k).cff);
      put("The constant of circuit "); put(k,1); put_line(" :");
      Standard_Complex_Vectors_io.put_line(s.crc(k).cst);
    end loop;
    if useqrls or usesvd then
      dx := Standard_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := Standard_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
      if usesvd then
        Standard_Newton_Convolution_Steps.SVD_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           svl,U,V,info,rcond,ewrk,wrk,scale);
        put("rcond :"); put(rcond,3); new_line;
      else
        Standard_Newton_Convolution_Steps.QR_Newton_Steps
          (standard_output,s,scf,dx,xd,maxit,nbrit,tol,absdx,fail,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      end if;
    elsif lurcond then
      Standard_Newton_Convolution_Steps.LU_Newton_Steps
        (standard_output,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale);
      put("rcond :"); put(rcond,3); new_line;
    else
      Standard_Newton_Convolution_Steps.LU_Newton_Steps
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

  procedure Standard_Run
              ( p : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim,deg : in integer32 ) is

    sol : constant Standard_Complex_Solutions.Link_to_Solution
        := Standard_Complex_Solutions.Head_Of(sols);
    maxit : integer32 := 0;
    scale,usesvd,useqrls,needrcond,staggered,inlined : boolean := false;
    indexed : boolean := false;
    overdet : constant boolean := (p'last > dim);
    ans : character;

  begin
    Prompt_for_Parameters
      (overdet,maxit,scale,usesvd,useqrls,needrcond,inlined);
    new_line;
    put("Apply coefficient convolution circuits ? (y/n) ");
    Ask_Yes_or_No(ans); 
    if ans = 'y' then
      put("Staggered on highest degree ? (y/n) "); Ask_Yes_or_No(ans);
      staggered := (ans = 'y');
      put("Indexed on lower degrees ? (y/n) "); Ask_Yes_or_No(ans);
      indexed := (ans = 'y');
      Standard_Coefficient_Run
        (p,sol,deg,maxit,scale,usesvd,useqrls,needrcond,
         staggered,inlined,indexed);
    else
      Standard_Run(p,sol,deg,maxit,scale,usesvd,useqrls,needrcond);
    end if;
  end Standard_Run;

  procedure Standard_Test ( deg : in integer32 ) is

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

  procedure Main is

    deg : integer32 := 0;

  begin
    new_line;
    put_line("Newton on convolution circuits in double precision.");
    new_line;
    put("Give the degree of the power series : "); get(deg); skip_line;
    Standard_Test(deg);
  end Main;

end Test_Standard_NewtConvSteps;
