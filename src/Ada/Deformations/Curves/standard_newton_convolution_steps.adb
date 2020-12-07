with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Newton_Convolutions;
with Newton_Coefficient_Convolutions;

package body Standard_Newton_Convolution_Steps is

-- NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE,
-- ON COEFFICIENT CONVOLUTION CIRCUITS :

  procedure LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      Newton_Coefficient_Convolutions.LU_Newton_Step
        (csr,scf,rx,ix,absdx,info,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Newton_Coefficient_Convolutions.LU_Newton_Step
        (file,csr,scf,rx,ix,absdx,info,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

-- NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 3 ...");
    end if;
    for k in 1..nbrit loop
      Standard_Newton_Convolutions.LU_Newton_Step
        (csr,scf,absdx,info,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 4 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Standard_Newton_Convolutions.LU_Newton_Step
        (file,csr,scf,absdx,info,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

-- NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE,
-- ON COEFFICIENT CONVOLUTIONS :

  procedure LU_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 5 ...");
    end if;
    for k in 1..nbrit loop
      Newton_Coefficient_Convolutions.LU_Newton_Step
        (csr,scf,rx,ix,absdx,rcond,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 6 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Newton_Coefficient_Convolutions.LU_Newton_Step
        (file,csr,scf,rx,ix,absdx,rcond,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

-- NEWTON STEPS WITH LU WITH CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 7 ...");
    end if;
    for k in 1..nbrit loop
      Standard_Newton_Convolutions.LU_Newton_Step
        (csr,scf,absdx,rcond,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; rcond : out double_float;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("LU_Newton_Steps 8 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Standard_Newton_Convolutions.LU_Newton_Step
        (file,csr,scf,absdx,rcond,ipvt,wrk,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

-- NEWTON STEPS WITH QR ON COEFFICIENT CONVOLUTIONS :

  procedure QR_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("QR_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      Newton_Coefficient_Convolutions.QR_Newton_Step
        (csr,scf,dx,xd,rx,ix,absdx,qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,
         scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("QR_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Newton_Coefficient_Convolutions.QR_Newton_Step
        (file,csr,scf,dx,xd,rx,ix,absdx,qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,
         scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end QR_Newton_Steps;

-- NEWTON STEPS WITH QR :

  procedure QR_Newton_Steps
              ( csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("QR_Newton_Steps 3 ...");
    end if;
    for k in 1..nbrit loop
      Standard_Newton_Convolutions.QR_Newton_Step
        (csr,scf,dx,xd,absdx,qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,
         scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end QR_Newton_Steps;

  procedure QR_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("QR_Newton_Steps 4 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Standard_Newton_Convolutions.QR_Newton_Step
        (file,csr,scf,dx,xd,absdx,qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,
         scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end QR_Newton_Steps;

-- NEWTON STEPS WITH SVD ON COEFFICIENT CONVOLUTIONS :

  procedure SVD_Newton_Steps
              ( csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("SVD_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      Newton_Coefficient_Convolutions.SVD_Newton_Step
        (csr,scf,dx,xd,rx,ix,absdx,svl,U,V,info,rcond,ewrk,wrkv,
         scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Coefficient_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("SVD_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Newton_Coefficient_Convolutions.SVD_Newton_Step
        (file,csr,scf,dx,xd,rx,ix,absdx,svl,U,V,info,rcond,ewrk,wrkv,
         scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end SVD_Newton_Steps;

-- NEWTON STEPS WITH SVD :

  procedure SVD_Newton_Steps
              ( csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps ...");
      put_line("SVD_Newton_Steps 3 ...");
    end if;
    for k in 1..nbrit loop
      Standard_Newton_Convolutions.SVD_Newton_Step
        (csr,scf,dx,xd,absdx,svl,U,V,info,rcond,ewrk,wrkv,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| ="); put(maxval,3);
        if idx < csr.vy'first
         then put_line(" too large");
         else put(" at index "); put(idx,1); new_line;
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end SVD_Newton_Steps;

  procedure SVD_Newton_Steps
              ( file : in file_type;
                csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean;
                svl : out Standard_Complex_Vectors.Vector;
                U,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv : in Standard_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : double_float;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in standard_newton_convolution_steps.");
      put_line("SVD_Newton_Steps 4 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      Standard_Newton_Convolutions.SVD_Newton_Step
        (file,csr,scf,dx,xd,absdx,svl,U,V,info,rcond,ewrk,wrkv,scale,vrblvl-1);
      Standard_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| ="); put(file,maxval,3);
        if idx < csr.vy'first
         then put_line(file," too large");
         else put(file," at index "); put(file,idx,1); new_line(file);
        end if;
      end if;
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end SVD_Newton_Steps;

end Standard_Newton_Convolution_Steps;
