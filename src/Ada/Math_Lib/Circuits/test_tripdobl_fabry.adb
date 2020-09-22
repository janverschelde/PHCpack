with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with TripDobl_Complex_Numbers;
with TripDobl_Complex_Numbers_io;        use TripDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs_io;        use TripDobl_Complex_VecVecs_io;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_Poly_Systems;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with TripDobl_Complex_Solutions;
with TripDobl_System_and_Solutions_io;
with Convergence_Radius_Estimates;
with TripDobl_Newton_Convolutions;

package body Test_TripDobl_Fabry is

  procedure TripDobl_Newton_Steps
              ( s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 ) is

    info : integer32;
    absdx,rcond : triple_double;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    ewrk : TripDobl_Complex_Vectors.Link_to_Vector
         := new TripDobl_Complex_Vectors.Vector(1..dim);
    wrk : TripDobl_Complex_Vectors.Link_to_Vector
        := new TripDobl_Complex_Vectors.Vector(1..s.neq);
    qraux,w1,w2,w3,w4,w5 : TripDobl_Complex_Vectors.Vector(1..s.neq);
    svl : TripDobl_Complex_Vectors.Vector(1..dim+1);
    U,V : TripDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    dx : TripDobl_Complex_VecVecs.VecVec(1..dim);
    xd : TripDobl_Complex_VecVecs.VecVec(0..deg);
    ans : character;
    scale,usesvd,useqrls,needrcond : boolean;

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
      dx := TripDobl_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := TripDobl_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
    end if;
    for k in 1..maxit loop
      put("Step "); put(k,1); put_line(" :");
      if usesvd then
        TripDobl_Newton_Convolutions.SVD_Newton_Step
          (standard_output,s,scf,dx,xd,absdx,svl,U,V,
           info,rcond,ewrk,wrk,scale);
      elsif useqrls then
        TripDobl_Newton_Convolutions.QR_Newton_Step
          (standard_output,s,scf,dx,xd,absdx,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      else
        if needrcond then
          TripDobl_Newton_Convolutions.LU_Newton_Step
            (standard_output,s,scf,absdx,rcond,ipvt,wrk,scale);
          put("rcond : "); put(rcond,3); new_line;
        else
          TripDobl_Newton_Convolutions.LU_Newton_Step
            (standard_output,s,scf,absdx,info,ipvt,wrk,scale);
          put("info : "); put(info,1); new_line;
        end if;
      end if;
      put("absdx : "); put(absdx,3);
      put("  Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    TripDobl_Complex_Vectors.Clear(ewrk);
    TripDobl_Complex_Vectors.Clear(wrk);
  end TripDobl_Newton_Steps;

  procedure Main is

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : TripDobl_Complex_Solutions.Solution_List;
    sol : TripDobl_Complex_Solutions.Link_to_Solution;
    dim,degree : integer32 := 0;
    nbr : natural32;

    use TripDobl_Speelpenning_Convolutions;

  begin
    TripDobl_System_and_Solutions_io.get(lp,sols);
    nbr := TripDobl_Complex_Solutions.Length_Of(sols);
    sol := TripDobl_Complex_Solutions.Head_Of(sols);
    dim := sol.n;
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the power series : "); get(degree);
    declare
      c : constant Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(degree));
      s : Link_to_System := Create(c,dim,degree);
      scf : constant TripDobl_Complex_VecVecs.VecVec(1..sol.n)
          := TripDobl_Newton_Convolutions.Series_Coefficients(sol.v,degree);
      z : TripDobl_Complex_Numbers.Complex_Number;
      r,err : triple_double;
      fail : boolean;
    begin
      Add_Continuation_Parameter(c);
      put_line("The coefficients of the convolution circuits :");
      for i in c'range loop
        put_line(c(i).cff);
      end loop;
      TripDobl_Newton_Steps(s,scf,lp'last,degree);
      Clear(s);
      Convergence_Radius_Estimates.Fabry(scf,z,r,err,fail);
      if not fail then
        put("z : "); put(z); 
        put("  error estimate : "); put(err,3); new_line;
        put("estimated radius : "); put(r,3); new_line;
      end if;
    end;
  end Main;

end Test_TripDobl_Fabry;
