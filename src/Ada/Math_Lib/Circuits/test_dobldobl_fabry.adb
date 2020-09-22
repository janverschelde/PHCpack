with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Matrices;
with Standard_Vector_Splitters;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Convolution_Splitters;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with Convergence_Radius_Estimates;
with DoblDobl_Newton_Convolutions;
with Newton_Coefficient_Convolutions;

package body Test_DoblDobl_Fabry is

  procedure DoblDobl_Newton_Steps
              ( s : in DoblDobl_Coefficient_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 ) is

    info : integer32;
    absdx,rcond : double_double;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    ewrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..s.neq);
    qraux,w1,w2,w3,w4,w5 : DoblDobl_Complex_Vectors.Vector(1..s.neq);
    svl : DoblDobl_Complex_Vectors.Vector(1..dim+1);
    U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    dx : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    xd : DoblDobl_Complex_VecVecs.VecVec(0..deg);
    ans : character;
    scale,usesvd,useqrls,needrcond : boolean := false;

  begin
    put("Apply scaling ? (y/n) "); Ask_Yes_or_No(ans);
    scale := (ans = 'y');
    put("Solve with SVD ? (y/n) "); Ask_Yes_or_No(ans);
    usesvd := (ans = 'y');
    if not usesvd then
      put("Apply least squares with QR ? (y/n) "); Ask_Yes_or_No(ans);
      useqrls := (ans = 'y');
      if not useqrls then
        put("Need condition number estimate ? (y/n) "); Ask_Yes_or_No(ans);
        needrcond := (ans = 'y');
      end if;
    end if;
    if useqrls or usesvd then
      dx := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := DoblDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
    end if;
    for k in 1..maxit loop
      put("Step "); put(k,1); put_line(" :");
      if usesvd then
        Newton_Coefficient_Convolutions.SVD_Newton_Step
          (standard_output,s,scf,dx,xd,rhx,ihx,rlx,ilx,absdx,svl,U,V,
           info,rcond,ewrk,wrk,scale);
      elsif useqrls then
        Newton_Coefficient_Convolutions.QR_Newton_Step
          (standard_output,s,scf,dx,xd,rhx,ihx,rlx,ilx,absdx,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      else
        if needrcond then
          Newton_Coefficient_Convolutions.LU_Newton_Step
            (standard_output,s,scf,rhx,ihx,rlx,ilx,absdx,rcond,ipvt,wrk,scale);
          put("  rcond : "); put(rcond,3); new_line;
        else
          Newton_Coefficient_Convolutions.LU_Newton_Step
            (standard_output,s,scf,rhx,ihx,rlx,ilx,absdx,info,ipvt,wrk,scale);
          put("  info : "); put(info,1); new_line;
        end if;
      end if;
      put("absdx : "); put(absdx,3);
      put("  Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    DoblDobl_Complex_Vectors.Clear(ewrk);
    DoblDobl_Complex_Vectors.Clear(wrk);
  end DoblDobl_Newton_Steps;

  procedure DoblDobl_Newton_Steps
              ( s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 ) is

    info : integer32;
    absdx,rcond : double_double;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    ewrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..s.neq);
    qraux,w1,w2,w3,w4,w5 : DoblDobl_Complex_Vectors.Vector(1..s.neq);
    svl : DoblDobl_Complex_Vectors.Vector(1..dim+1);
    U,V : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    dx : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    xd : DoblDobl_Complex_VecVecs.VecVec(0..deg);
    ans : character;
    scale,usesvd,useqrls,needrcond : boolean := false;

  begin
    put("Apply scaling ? (y/n) "); Ask_Yes_or_No(ans);
    scale := (ans = 'y');
    put("Solve with SVD ? (y/n) "); Ask_Yes_or_No(ans);
    usesvd := (ans = 'y');
    if not usesvd then
      put("Apply least squares with QR ? (y/n) "); Ask_Yes_or_No(ans);
      useqrls := (ans = 'y');
      if not useqrls then
        put("Need condition number estimate ? (y/n) "); Ask_Yes_or_No(ans);
        needrcond := (ans = 'y');
      end if;
    end if;
    if useqrls or usesvd then
      dx := DoblDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := DoblDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
    end if;
    for k in 1..maxit loop
      put("Step "); put(k,1); put_line(" :");
      if usesvd then
        DoblDobl_Newton_Convolutions.SVD_Newton_Step
          (standard_output,s,scf,dx,xd,absdx,svl,U,V,
           info,rcond,ewrk,wrk,scale);
      elsif useqrls then
        DoblDobl_Newton_Convolutions.QR_Newton_Step
          (standard_output,s,scf,dx,xd,absdx,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      else
        if needrcond then
          DoblDobl_Newton_Convolutions.LU_Newton_Step
            (standard_output,s,scf,absdx,rcond,ipvt,wrk,scale);
          put("  rcond : "); put(rcond,3); new_line;
        else
          DoblDobl_Newton_Convolutions.LU_Newton_Step
            (standard_output,s,scf,absdx,info,ipvt,wrk,scale);
          put("  info : "); put(info,1); new_line;
        end if;
      end if;
      put("absdx : "); put(absdx,3);
      put("  Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    DoblDobl_Complex_Vectors.Clear(ewrk);
    DoblDobl_Complex_Vectors.Clear(wrk);
  end DoblDobl_Newton_Steps;

  procedure Main is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    sol : DoblDobl_Complex_Solutions.Link_to_Solution;
    dim,deg : integer32 := 0;
    nbr : natural32;
    ans : character;

  begin
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    nbr := DoblDobl_Complex_Solutions.Length_Of(sols);
    sol := DoblDobl_Complex_Solutions.Head_Of(sols);
    dim := sol.n;
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the power series : "); get(deg);
    declare
      c : constant DoblDobl_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : DoblDobl_Speelpenning_Convolutions.Link_to_System
        := DoblDobl_Speelpenning_Convolutions.Create(c,dim,deg);
      cs : DoblDobl_Coefficient_Convolutions.Link_to_System;
      scf : constant DoblDobl_Complex_VecVecs.VecVec(1..sol.n)
          := DoblDobl_Newton_Convolutions.Series_Coefficients(sol.v,deg);
      rhx : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
      ihx : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
      rlx : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
      ilx : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
      z : DoblDobl_Complex_Numbers.Complex_Number;
      r,err : double_double;
      fail : boolean;
    begin
      Add_Continuation_Parameter(c);
      put_line("The coefficients of the convolution circuits :");
      for i in c'range loop
        put_line(c(i).cff);
      end loop;
      put("Run with coefficient convolutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        cs := DoblDobl_Convolution_Splitters.Split(s);
        DoblDobl_Newton_Steps(cs,scf,rhx,ihx,rlx,ilx,lp'last,deg);
      else
        DoblDobl_Newton_Steps(s,scf,lp'last,deg);
      end if;
      DoblDobl_Speelpenning_Convolutions.Clear(s);
      DoblDobl_Coefficient_Convolutions.Clear(cs);
      Convergence_Radius_Estimates.Fabry(scf,z,r,err,fail);
      if not fail then
        put("z : "); put(z); 
        put("  error estimate : "); put(err,3); new_line;
        put("estimated radius : "); put(r,3); new_line;
      end if;
    end;
  end Main;

end Test_DoblDobl_Fabry;
