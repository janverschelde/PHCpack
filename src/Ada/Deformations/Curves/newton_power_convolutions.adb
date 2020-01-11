with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Newton_Convolutions;                use Newton_Convolutions;

package body Newton_Power_Convolutions is

  procedure LU_Newton_Steps
              ( csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                nbrit : in integer32; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 1 ...");
    end if;
    for k in 1..nbrit loop
      LU_Newton_Step(csr,scf,info,ipvt,wrk);
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type; 
                csr : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                nbrit : in integer32; info : out integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in newton_power_convolutions.LU_Newton_Steps 2 ...");
    end if;
    for k in 1..nbrit loop
      put(file,"Step "); put(file,k,1); put_line(file," :");
      LU_Newton_Step(file,csr,scf,info,ipvt,wrk);
    end loop;
  end LU_Newton_Steps;

end Newton_Power_Convolutions;
