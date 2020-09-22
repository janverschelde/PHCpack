with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;
with Standard_Vector_Splitters;
with Standard_Complex_Poly_Systems;
with Standard_Convolution_Splitters;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Homotopy_Convolution_Circuits;      use Homotopy_Convolution_Circuits;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Convergence_Radius_Estimates;
with Standard_Newton_Convolutions;
with Newton_Coefficient_Convolutions;

package body Test_Standard_Fabry is

  procedure Standard_Newton_Steps
              ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 ) is

    info : integer32;
    rcond,absdx : double_float;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..s.neq);
    ewrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..dim);
    qraux,w1,w2,w3,w4,w5 : Standard_Complex_Vectors.Vector(1..s.neq);
    svl : Standard_Complex_Vectors.Vector(1..dim+1);
    U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    dx : Standard_Complex_VecVecs.VecVec(1..dim);
    xd : Standard_Complex_VecVecs.VecVec(0..deg);
    ans : character;
    scale,usesvd,useqrls,needrcond : boolean := false;
    staggered,inlined,indexed : boolean := false;
    wrkdeg,wrkidx,idxtoldx : integer32 := 0;
    rc,ic,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;

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
        put("Inlined solver ? (y/n) "); Ask_Yes_or_No(ans);
        inlined := (ans = 'y');
      end if;
    end if;
    put("Staggered degrees ? (y/n) "); Ask_Yes_or_No(ans);
    staggered := (ans = 'y');
    if useqrls or usesvd then
      dx := Standard_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := Standard_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
    end if;
    if inlined then -- allocate work space vectors
      Standard_Floating_VecVecVecs.Allocate(rv,1,deg,1,dim,1,dim);
      Standard_Floating_VecVecVecs.Allocate(iv,1,deg,1,dim,1,dim);
      rc := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
      ic := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
      rb := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
      ib := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
      ry := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
      iy := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
      if staggered then
        put("Indexed ? (y/n) "); Ask_Yes_or_No(ans);
        indexed := (ans = 'y');
      end if;
    end if;
    if staggered
     then wrkdeg := 1; wrkidx := 0;
    end if;
    for k in 1..maxit loop
      put("Step "); put(k,1); put_line(" :");
      if not staggered then
        if usesvd then
          Newton_Coefficient_Convolutions.SVD_Newton_Step
            (standard_output,s,scf,dx,xd,rx,ix,absdx,svl,U,V,
             info,rcond,ewrk,wrk,scale);
        elsif useqrls then
          Newton_Coefficient_Convolutions.QR_Newton_Step
            (standard_output,s,scf,dx,xd,rx,ix,absdx,
             qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
        else
          if needrcond then
            if inlined then
              Newton_Coefficient_Convolutions.Inlined_LU_Newton_Step
                (standard_output,s,scf,rx,ix,absdx,rcond,ipvt,
                 rc,ic,rv,iv,rb,ib,ry,iy,scale);
            else
              Newton_Coefficient_Convolutions.LU_Newton_Step
                (standard_output,s,scf,rx,ix,absdx,rcond,ipvt,wrk,scale);
            end if;
            put("  rcond :"); put(rcond,3); new_line;
          else
            if inlined then
              Newton_Coefficient_Convolutions.Inlined_LU_Newton_Step
                (standard_output,s,scf,rx,ix,absdx,info,ipvt,
                 rc,ic,rv,iv,rb,ib,ry,iy,scale);
            else
              Newton_Coefficient_Convolutions.LU_Newton_Step
                (standard_output,s,scf,rx,ix,absdx,info,ipvt,wrk,scale);
            end if;
            put("  info : "); put(info,1); new_line;
          end if;
        end if;
      else
        if usesvd then
          Newton_Coefficient_Convolutions.SVD_Newton_Step
            (standard_output,wrkdeg,s,scf,dx,xd,rx,ix,absdx,svl,U,V,
             info,rcond,ewrk,wrk,scale);
        elsif useqrls then
          Newton_Coefficient_Convolutions.QR_Newton_Step
            (standard_output,wrkdeg,s,scf,dx,xd,rx,ix,absdx,
             qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
        else
          if needrcond then
            if inlined then
              if indexed then
                put("wrkidx : "); put(wrkidx,1);
                put("  wrkdeg : "); put(wrkdeg,1); new_line;
                Newton_Coefficient_Convolutions.Inlined_LU_Newton_Step
                  (standard_output,wrkidx,wrkdeg,s,scf,rx,ix,1.0E-8,
                   idxtoldx,absdx,rcond,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
              else
                Newton_Coefficient_Convolutions.Inlined_LU_Newton_Step
                  (standard_output,wrkdeg,s,scf,rx,ix,absdx,rcond,ipvt,
                   rc,ic,rv,iv,rb,ib,ry,iy,scale);
              end if;
            else
              Newton_Coefficient_Convolutions.LU_Newton_Step
                (standard_output,wrkdeg,s,scf,rx,ix,absdx,rcond,ipvt,wrk,scale);
            end if;
            put("  rcond :"); put(rcond,3); new_line;
          else
            if inlined then
              if indexed then
                put("wrkidx : "); put(wrkidx,1);
                put("  wrkdeg : "); put(wrkdeg,1); new_line;
                Newton_Coefficient_Convolutions.Inlined_LU_Newton_Step
                  (standard_output,wrkidx,wrkdeg,s,scf,rx,ix,1.0E-8,
                   idxtoldx,absdx,info,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
              else
                Newton_Coefficient_Convolutions.Inlined_LU_Newton_Step
                  (standard_output,wrkdeg,s,scf,rx,ix,absdx,info,ipvt,
                   rc,ic,rv,iv,rb,ib,ry,iy,scale);
              end if;
            else
              Newton_Coefficient_Convolutions.LU_Newton_Step
                (standard_output,wrkdeg,s,scf,rx,ix,absdx,info,ipvt,wrk,scale);
            end if;
            put("  info : "); put(info,1); new_line;
          end if;
        end if;
      end if;
      put("absdx :"); put(absdx,3);
      put("  Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      if staggered then
        if indexed
         then wrkidx := idxtoldx+1;
        end if;
        wrkdeg := 2*wrkdeg;
        if wrkdeg > deg
         then wrkdeg := deg;
        end if;
      end if;
      exit when indexed and then (wrkidx > deg);
    end loop;
    if indexed then
      loop
        put("One extra step ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        Newton_Coefficient_Convolutions.Inlined_LU_Newton_Step
          (standard_output,s,scf,rx,ix,absdx,info,ipvt,
           rc,ic,rv,iv,rb,ib,ry,iy,scale);
      end loop;
    end if;
    Standard_Complex_Vectors.Clear(wrk);
    Standard_Complex_Vectors.Clear(ewrk);
  end Standard_Newton_Steps;

  procedure Standard_Newton_Steps
              ( s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                dim,deg : in integer32; maxit : in natural32 := 100 ) is

    info : integer32;
    rcond,absdx : double_float;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..s.neq);
    ewrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..dim);
    qraux,w1,w2,w3,w4,w5 : Standard_Complex_Vectors.Vector(1..s.neq);
    svl : Standard_Complex_Vectors.Vector(1..dim+1);
    U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    dx : Standard_Complex_VecVecs.VecVec(1..dim);
    xd : Standard_Complex_VecVecs.VecVec(0..deg);
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
      dx := Standard_Speelpenning_Convolutions.Allocate_Coefficients(dim,deg);
      xd := Standard_Speelpenning_Convolutions.Linearized_Allocation(dim,deg);
    end if;
    for k in 1..maxit loop
      put("Step "); put(k,1); put_line(" :");
      if usesvd then
        Standard_Newton_Convolutions.SVD_Newton_Step
          (standard_output,s,scf,dx,xd,absdx,svl,U,V,
           info,rcond,ewrk,wrk,scale);
      elsif useqrls then
        Standard_Newton_Convolutions.QR_Newton_Step
          (standard_output,s,scf,dx,xd,absdx,
           qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,scale);
      else
        if needrcond then
          Standard_Newton_Convolutions.LU_Newton_Step
            (standard_output,s,scf,absdx,rcond,ipvt,wrk,scale);
          put("  rcond :"); put(rcond,3); new_line;
        else
          Standard_Newton_Convolutions.LU_Newton_Step
            (standard_output,s,scf,absdx,info,ipvt,wrk,scale);
          put("  info : "); put(info,1); new_line;
        end if;
      end if;
      put("absdx :"); put(absdx,3);
      put("  Continue ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Complex_Vectors.Clear(wrk);
    Standard_Complex_Vectors.Clear(ewrk);
  end Standard_Newton_Steps;

  procedure Main is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    sol : Standard_Complex_Solutions.Link_to_Solution;
    dim,deg : integer32 := 0;
    nbr : natural32;
    ans : character;

  begin
    Standard_System_and_Solutions_io.get(lp,sols);
    nbr := Standard_Complex_Solutions.Length_Of(sols);
    sol := Standard_Complex_Solutions.Head_Of(sols);
    dim := sol.n;
    new_line;
    put("Read "); put(nbr,1); put(" solutions in dimension ");
    put(dim,1); put_line(".");
    new_line;
    put("Give the degree of the power series : "); get(deg);
    declare
      c : constant Standard_Speelpenning_Convolutions.Circuits(lp'range)
        := Make_Convolution_Circuits(lp.all,natural32(deg));
      s : Standard_Speelpenning_Convolutions.Link_to_System
        := Standard_Speelpenning_Convolutions.Create(c,dim,deg);
      cs : Standard_Coefficient_Convolutions.Link_to_System;
      scf : constant Standard_Complex_VecVecs.VecVec(1..sol.n)
          := Standard_Newton_Convolutions.Series_Coefficients(sol.v,deg);
      rx : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
      ix : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
      z : Standard_Complex_Numbers.Complex_Number;
      r,err : double_float;
      fail : boolean;
    begin
      Add_Continuation_Parameter(c);
      put_line("The coefficients of the convolution circuits :");
      for i in c'range loop
        put_line(c(i).cff);
      end loop;
      new_line;
      put("Run with coefficient convolutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        cs := Standard_Convolution_Splitters.Split(s);
        Standard_Newton_Steps(cs,scf,rx,ix,lp'last,deg);
      else
        Standard_Newton_Steps(s,scf,lp'last,deg);
      end if;
      Standard_Speelpenning_Convolutions.Clear(s);
      Standard_Coefficient_Convolutions.Clear(cs);
      Convergence_Radius_Estimates.Fabry(scf,z,r,err,fail);
      if not fail then
        put("z : "); put(z); 
        put("  error estimate :"); put(err,3); new_line;
        put("estimated radius :"); put(r,3); new_line;
      end if;
    end;
  end Main;

end Test_Standard_Fabry;
