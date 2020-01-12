with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Newton_Convolutions;                use Newton_Convolutions;

package body Newton_Power_Convolutions is

  procedure LU_Newton_Steps
              ( csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_float;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is
  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      LU_Newton_Step(csr,scf,absdx,info,ipvt,wrk);
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
                vrblvl : in integer32 := 0 ) is
  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      LU_Newton_Step(file,csr,scf,absdx,info,ipvt,wrk);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( csr : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_double;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is
  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 3 ...");
    end if;
    for k in 1..nbrit loop
      LU_Newton_Step(csr,scf,absdx,info,ipvt,wrk);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out double_double;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is
  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 4 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      LU_Newton_Step(file,csr,scf,absdx,info,ipvt,wrk);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( csr : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; absdx : out quad_double;
                fail : out boolean; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is
  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 5 ...");
    end if;
    for k in 1..nbrit loop
      LU_Newton_Step(csr,scf,absdx,info,ipvt,wrk);
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
                vrblvl : in integer32 := 0 ) is
  begin
    fail := true; nbrit := maxit;
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 6 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      LU_Newton_Step(file,csr,scf,absdx,info,ipvt,wrk);
      if absdx <= tol
       then fail := false; nbrit := k; exit;
      end if;
    end loop;
  end LU_Newton_Steps;

end Newton_Power_Convolutions;
