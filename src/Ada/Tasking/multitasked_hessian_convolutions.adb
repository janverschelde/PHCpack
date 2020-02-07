with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with DoblDobl_Mathematical_Functions;
with QuadDobl_Mathematical_Functions;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Hessian_Convolution_Circuits;       use Hessian_Convolution_Circuits;
with Jacobian_Convolution_Circuits;      use Jacobian_Convolution_Circuits;
with Multitasking;

package body Multitasked_Hessian_Convolutions is

  function Allocate ( nbr,dim : integer32 )
                    return Standard_Complex_VecMats.VecMat is

    res : Standard_Complex_VecMats.VecMat(1..nbr);

  begin
    for k in res'range loop
      declare
        mat : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := Standard_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        res(k) := new Standard_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( nbr,dim : integer32 )
                    return DoblDobl_Complex_VecMats.VecMat is

    res : DoblDobl_Complex_VecMats.VecMat(1..nbr);

  begin
    for k in res'range loop
      declare
        mat : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := DoblDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        res(k) := new DoblDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( nbr,dim : integer32 )
                    return QuadDobl_Complex_VecMats.VecMat is

    res : QuadDobl_Complex_VecMats.VecMat(1..nbr);

  begin
    for k in res'range loop
      declare
        mat : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := QuadDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        res(k) := new QuadDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( nbr,dim : integer32 )
                    return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..nbr);

  begin
    for i in res'range loop
      declare
        v : constant Standard_Complex_Vectors.Vector(1..dim)
          := (1..dim => Standard_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( nbr,dim : integer32 )
                    return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..nbr);

  begin
    for i in res'range loop
      declare
        v : constant DoblDobl_Complex_Vectors.Vector(1..dim)
          := (1..dim => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  function Allocate ( nbr,dim : integer32 )
                    return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..nbr);

  begin
    for i in res'range loop
      declare
        v : constant QuadDobl_Complex_Vectors.Vector(1..dim)
          := (1..dim => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        res(i) := new QuadDobl_Complex_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Allocate;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_Vectors.Vector;
                jmsvls : out Standard_Complex_Vectors.Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
		verbose : in boolean := true ) is

    A,U,V : Standard_Complex_VecMats.VecMat(1..nbt);
    e : Standard_Complex_VecVecs.VecVec(1..nbt);

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls : Standard_Complex_Vectors.Link_to_Vector;

    begin
      if i = 1 then
        put_line("the first first tasks compute the Jacobian and its SVD");
        pA := A(1); pU := U(1); pV := V(1); pe := e(1);
        Singular_Values(s.crc,x,pA.all,pU.all,pV.all,pe.all,jmsvls);
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx);
        Singular_Values(s.crc(idx),x,pA.all,pU.all,pV.all,pe.all,vls.all);
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
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls : Standard_Complex_Vectors.Link_to_Vector;

    begin
      if i = 1 then
        pA := A(1); pU := U(1); pV := V(1); pe := e(1);
        Singular_Values(s.crc,x,pA.all,pU.all,pV.all,pe.all,jmsvls);
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx);
        Singular_Values(s.crc(idx),x,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    A := Allocate(nbt,s.dim);
    U := Allocate(nbt,s.dim);
    V := Allocate(nbt,s.dim);
    e := Allocate(nbt,s.dim);
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
                s : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Vector;
                jmsvls : out DoblDobl_Complex_Vectors.Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
		verbose : in boolean := true ) is

    A,U,V : DoblDobl_Complex_VecMats.VecMat(1..nbt);
    e : DoblDobl_Complex_VecVecs.VecVec(1..nbt);

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : DoblDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      if i = 1 then
        put_line("the first first tasks compute the Jacobian and its SVD");
        pA := A(1); pU := U(1); pV := V(1); pe := e(1);
        Singular_Values(s.crc,x,pA.all,pU.all,pV.all,pe.all,jmsvls);
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx);
        Singular_Values(s.crc(idx),x,pA.all,pU.all,pV.all,pe.all,vls.all);
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
      pe,vls : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      if i = 1 then
        pA := A(1); pU := U(1); pV := V(1); pe := e(1);
        Singular_Values(s.crc,x,pA.all,pU.all,pV.all,pe.all,jmsvls);
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx);
        Singular_Values(s.crc(idx),x,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    A := Allocate(nbt,s.dim);
    U := Allocate(nbt,s.dim);
    V := Allocate(nbt,s.dim);
    e := Allocate(nbt,s.dim);
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
                s : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Vector;
                jmsvls : out QuadDobl_Complex_Vectors.Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
		verbose : in boolean := true ) is

    A,U,V : QuadDobl_Complex_VecMats.VecMat(1..nbt);
    e : QuadDobl_Complex_VecVecs.VecVec(1..nbt);

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      if i = 1 then
        put_line("the first first tasks compute the Jacobian and its SVD");
        pA := A(1); pU := U(1); pV := V(1); pe := e(1);
        Singular_Values(s.crc,x,pA.all,pU.all,pV.all,pe.all,jmsvls);
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx);
        Singular_Values(s.crc(idx),x,pA.all,pU.all,pV.all,pe.all,vls.all);
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
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      if i = 1 then
        pA := A(1); pU := U(1); pV := V(1); pe := e(1);
        Singular_Values(s.crc,x,pA.all,pU.all,pV.all,pe.all,jmsvls);
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx);
        Singular_Values(s.crc(idx),x,pA.all,pU.all,pV.all,pe.all,vls.all);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    A := Allocate(nbt,s.dim);
    U := Allocate(nbt,s.dim);
    V := Allocate(nbt,s.dim);
    e := Allocate(nbt,s.dim);
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
              ( jmsvls : Standard_Complex_Vectors.Vector;
                values : Standard_Complex_VecVecs.VecVec )
              return double_float is

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
              ( jmsvls : DoblDobl_Complex_Vectors.Vector;
                values : DoblDobl_Complex_VecVecs.VecVec )
              return double_double is

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
              ( jmsvls : QuadDobl_Complex_Vectors.Vector;
                values : QuadDobl_Complex_VecVecs.VecVec )
              return quad_double is

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

end Multitasked_Hessian_Convolutions;
