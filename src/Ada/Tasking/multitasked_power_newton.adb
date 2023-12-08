with Standard_Integer_Vectors;
with Multitasked_Series_Linearization;
with Multitasked_Newton_Convolutions;    use Multitasked_Newton_Convolutions;

package body Multitasked_Power_Newton is

  procedure Standard_Run
              ( nbt,dim,maxit : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_float; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : Standard_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,tol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    Standard_Complex_VecVecs.Clear(wks);
  end Standard_Run;

  procedure Standard_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                scf : in Standard_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_float; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : Standard_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    Standard_Complex_VecVecs.Clear(wks);
  end Standard_Run;

  procedure DoblDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : DoblDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    ddtol : constant double_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    DoblDobl_Complex_VecVecs.Clear(wks);
  end DoblDobl_Run;

  procedure DoblDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DoblDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out double_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : DoblDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    ddtol : constant double_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,ddtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    DoblDobl_Complex_VecVecs.Clear(wks);
  end DoblDobl_Run;

  procedure TripDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out triple_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : TripDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    tdtol : constant triple_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    TripDobl_Complex_VecVecs.Clear(wks);
  end TripDobl_Run;

  procedure TripDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in TripDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in TripDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out triple_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : TripDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    tdtol : constant triple_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,tdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    TripDobl_Complex_VecVecs.Clear(wks);
  end TripDobl_Run;

  procedure QuadDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out quad_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : QuadDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant quad_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    QuadDobl_Complex_VecVecs.Clear(wks);
  end QuadDobl_Run;

  procedure QuadDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in QuadDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out quad_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : QuadDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant quad_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    QuadDobl_Complex_VecVecs.Clear(wks);
  end QuadDobl_Run;

  procedure PentDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in PentDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out penta_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : PentDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant penta_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    PentDobl_Complex_VecVecs.Clear(wks);
  end PentDobl_Run;

  procedure PentDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in PentDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in PentDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out penta_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : PentDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant penta_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    PentDobl_Complex_VecVecs.Clear(wks);
  end PentDobl_Run;

  procedure OctoDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in OctoDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out octo_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : OctoDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant octo_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    OctoDobl_Complex_VecVecs.Clear(wks);
  end OctoDobl_Run;

  procedure OctoDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in OctoDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in OctoDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out octo_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : OctoDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant octo_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    OctoDobl_Complex_VecVecs.Clear(wks);
  end OctoDobl_Run;

  procedure DecaDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DecaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out deca_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : DecaDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant deca_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    DecaDobl_Complex_VecVecs.Clear(wks);
  end DecaDobl_Run;

  procedure DecaDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in DecaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in DecaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out deca_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : DecaDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant deca_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    DecaDobl_Complex_VecVecs.Clear(wks);
  end DecaDobl_Run;

  procedure HexaDobl_Run
              ( nbt,dim,maxit : in integer32;
                s : in HexaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in HexaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out hexa_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : HexaDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant hexa_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (standard_output,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,
          info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    HexaDobl_Complex_VecVecs.Clear(wks);
  end HexaDobl_Run;

  procedure HexaDobl_Run
              ( file : in file_type; nbt,dim,maxit : in integer32;
                s : in HexaDobl_Speelpenning_Convolutions.Link_to_System;
                scf : in HexaDobl_Complex_VecVecs.VecVec;
                tol : in double_float; estco : in boolean;
                fail : out boolean; info,nbrit : out integer32;
                rcond,absdx : out hexa_double; 
                output : in boolean := false;
                verbose : in boolean := true ) is

    wks : HexaDobl_Complex_VecVecs.VecVec(1..nbt)
        := Multitasked_Series_Linearization.Allocate_Work_Space(nbt,dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    qdtol : constant hexa_double := create(tol);

  begin
    if verbose then
      if estco then
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks,output);
      else
        Multitasked_LU_Newton_Steps
         (file,nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks,output);
      end if;
    else
      if estco then
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,rcond,ipvt,wks);
      else
        Multitasked_LU_Newton_Steps
         (nbt,s,scf,maxit,nbrit,qdtol,absdx,fail,info,ipvt,wks);
      end if;
    end if;
    HexaDobl_Complex_VecVecs.Clear(wks);
  end HexaDobl_Run;

end Multitasked_Power_Newton;
