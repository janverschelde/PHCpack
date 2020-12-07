with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Newton_Convolutions;

package body QuadDobl_Newton_Convolution_Steps is

-- NEWTON STEPS WITH LU WITHOUT CONDITION NUMBER ESTIMATE :

  procedure LU_Newton_Steps
              ( csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("LU_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      QuadDobl_Newton_Convolutions.LU_Newton_Step
        (csr,scf,absdx,info,ipvt,wrk,scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| = "); put(maxval,3);
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
                csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("LU_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      QuadDobl_Newton_Convolutions.LU_Newton_Step
        (file,csr,scf,absdx,info,ipvt,wrk,scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| = "); put(file,maxval,3);
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
              ( csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("LU_Newton_Steps 3 ...");
    end if;
    for k in 1..nbrit loop
      QuadDobl_Newton_Convolutions.LU_Newton_Step
        (csr,scf,absdx,rcond,ipvt,wrk,scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| = "); put(maxval,3);
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
                csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean; rcond : out quad_double;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
		scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("LU_Newton_Steps 4 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      QuadDobl_Newton_Convolutions.LU_Newton_Step
        (file,csr,scf,absdx,rcond,ipvt,wrk,scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| = "); put(file,maxval,3);
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

-- NEWTON STEPS WITH QR :

  procedure QR_Newton_Steps
              ( csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("QR_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      QuadDobl_Newton_Convolutions.QR_Newton_Step
        (csr,scf,dx,xd,absdx,qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,
	 scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| = "); put(maxval,3);
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
                csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_Vectors.Vector;
                info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("QR_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      QuadDobl_Newton_Convolutions.QR_Newton_Step
        (file,csr,scf,dx,xd,absdx,qraux,w1,w2,w3,w4,w5,info,ipvt,wrk,
         scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| = "); put(file,maxval,3);
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

-- NEWTON STEPS WITH SVD :

  procedure SVD_Newton_Steps
              ( csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean;
                svl : out QuadDobl_Complex_Vectors.Vector;
                U,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("SVD_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      QuadDobl_Newton_Convolutions.SVD_Newton_Step
        (csr,scf,dx,xd,absdx,svl,U,V,info,rcond,ewrk,wrkv,scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put("max |dx| = "); put(maxval,3);
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
                csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf,dx,xd : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean;
                svl : out QuadDobl_Complex_Vectors.Vector;
                U,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in QuadDobl_Complex_Vectors.Link_to_Vector;
                scale : in boolean := true; verbose : in boolean := true;
                vrblvl : in integer32 := 0 ) is

    maxval : quad_double;
    idx : integer32;

  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0 then
      put("-> in quaddobl_newton_convolution_steps.");
      put_line("SVD_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      QuadDobl_Newton_Convolutions.SVD_Newton_Step
        (file,csr,scf,dx,xd,absdx,svl,U,V,info,rcond,ewrk,wrkv,scale,vrblvl-1);
      QuadDobl_Newton_Convolutions.MaxIdx(csr.vy,tol,maxval,idx);
      if verbose then
        put(file,"max |dx| = "); put(file,maxval,3);
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

end QuadDobl_Newton_Convolution_Steps;
