with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with DoblDobl_Mathematical_Functions;
with QuadDobl_Mathematical_Functions;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_VecMats;
with Multitasking;

package body Multitasked_Hessian_Circuits is

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in res'range loop
      declare
        v : constant Standard_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => Standard_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in res'range loop
      declare
        v : constant DoblDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( neq,dim : integer32; neqstart,dimstart : integer32 )
                    return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(neqstart..neq);

  begin
    for i in res'range loop
      declare
        v : constant QuadDobl_Complex_Vectors.Vector(dimstart..dim)
          := (dimstart..dim => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new QuadDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
		verbose : in boolean := true ) is

    A,U,V : Standard_Complex_VecMats.VecMat(1..nbt);
    e : Standard_Complex_VecVecs.VecVec(1..nbt);
    yd : Standard_Complex_VecVecs.VecVec(1..s.neq);

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd : Standard_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        Standard_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd : Standard_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        Standard_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    A := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    U := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    V := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    e := Allocate(nbt,s.dim,1,1);
    yd := Allocate(s.neq,s.dim,1,0); -- space for all gradients
    Standard_Complex_Circuits.Power_Table(s.mxe,x,s.pwt);
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
    Standard_Complex_VecMats.Clear(A);
    Standard_Complex_VecMats.Clear(U);
    Standard_Complex_VecMats.Clear(V);
    Standard_Complex_VecVecs.Clear(e);
  end Multitasked_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
		verbose : in boolean := true ) is

    A,U,V : DoblDobl_Complex_VecMats.VecMat(1..nbt);
    e : DoblDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : DoblDobl_Complex_VecVecs.VecVec(1..s.neq);

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : DoblDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        DoblDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.

      idx : integer32 := i; 
      pA,pU,pV : DoblDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        DoblDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    A := DoblDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := DoblDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := DoblDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := Allocate(nbt,s.dim,1,1);
    yd := Allocate(s.neq,s.dim,1,0); -- space for gradients
    DoblDobl_Complex_Circuits.Power_Table(s.mxe,x,s.pwt);
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
    DoblDobl_Complex_VecMats.Clear(A);
    DoblDobl_Complex_VecMats.Clear(U);
    DoblDobl_Complex_VecMats.Clear(V);
    DoblDobl_Complex_VecVecs.Clear(e);
  end Multitasked_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
		verbose : in boolean := true ) is

    A,U,V : QuadDobl_Complex_VecMats.VecMat(1..nbt);
    e : QuadDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : QuadDobl_Complex_VecVecs.VecVec(1..s.neq);

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        QuadDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        QuadDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    A := QuadDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := QuadDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := QuadDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := Allocate(nbt,s.dim,1,1);
    yd := Allocate(s.neq,s.dim,1,0); -- space for gradients
    QuadDobl_Complex_Circuits.Power_Table(s.mxe,x,s.pwt);
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
    QuadDobl_Complex_VecMats.Clear(A);
    QuadDobl_Complex_VecMats.Clear(U);
    QuadDobl_Complex_VecMats.Clear(V);
    QuadDobl_Complex_VecVecs.Clear(e);
  end Multitasked_Singular_Values;

  function Standard_Distance
              ( values : Standard_Complex_VecVecs.VecVec )
              return double_float is

    jmsvls : constant Standard_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant double_float
           := Standard_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last));
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    sumvls : double_float := 0.0;
    val,nrm : double_float;

  begin
    for k in values'range loop
      lnk := values(k);
      val := Standard_Complex_Numbers.REAL_PART(lnk(lnk'first));
      sumvls := sumvls + val*val;
    end loop;
    nrm := Standard_Mathematical_Functions.SQRT(sumvls);
    return (2.0*sigma1)/nrm;
  end Standard_Distance;

  function DoblDobl_Distance
              ( values : DoblDobl_Complex_VecVecs.VecVec )
              return double_double is

    jmsvls : constant DoblDobl_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant double_double
           := DoblDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last));
    lnk : DoblDobl_Complex_Vectors.Link_to_Vector;
    sumvls : double_double := create(0.0);
    val,nrm : double_double;

  begin
    for k in values'range loop
      lnk := values(k);
      val := DoblDobl_Complex_Numbers.REAL_PART(lnk(lnk'first));
      sumvls := sumvls + val*val;
    end loop;
    nrm := DoblDobl_Mathematical_Functions.SQRT(sumvls);
    return (2.0*sigma1)/nrm;
  end DoblDobl_Distance;

  function QuadDobl_Distance
              ( values : QuadDobl_Complex_VecVecs.VecVec )
              return quad_double is

    jmsvls : constant QuadDobl_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant quad_double
           := QuadDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last));
    lnk : QuadDobl_Complex_Vectors.Link_to_Vector;
    sumvls : quad_double := create(0.0);
    val,nrm : quad_double;

  begin
    for k in values'range loop
      lnk := values(k);
      val := QuadDobl_Complex_Numbers.REAL_PART(lnk(lnk'first));
      sumvls := sumvls + val*val;
    end loop;
    nrm := QuadDobl_Mathematical_Functions.SQRT(sumvls);
    return (2.0*sigma1)/nrm;
  end QuadDobl_Distance;

end Multitasked_Hessian_Circuits;
