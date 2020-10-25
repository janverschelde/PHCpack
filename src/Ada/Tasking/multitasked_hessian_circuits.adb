with text_io;                            use text_io;
with Standard_Mathematical_Functions;
with DoblDobl_Mathematical_Functions;
with TripDobl_Mathematical_Functions;
with QuadDobl_Mathematical_Functions;
with PentDobl_Mathematical_Functions;
with OctoDobl_Mathematical_Functions;
with DecaDobl_Mathematical_Functions;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with Standard_Floating_VecVecs;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with TripDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with PentDobl_Complex_Matrices;
with OctoDobl_Complex_Matrices;
with DecaDobl_Complex_Matrices;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with TripDobl_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with PentDobl_Vector_Splitters;
with OctoDobl_Vector_Splitters;
with DecaDobl_Vector_Splitters;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with TripDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;
with PentDobl_Complex_Singular_Values;
with OctoDobl_Complex_Singular_Values;
with DecaDobl_Complex_Singular_Values;
with Semaphore;

package body Multitasked_Hessian_Circuits is

  procedure Allocate_Hessian_Spaces
              ( dim : in integer32;
                hrp,hip : out Standard_Floating_VecVecVecs.VecVecVec ) is
  begin
    for k in hrp'range loop
      Standard_Coefficient_Circuits.Allocate_Hessian_Space(dim,hrp(k),hip(k));
    end loop;
  end Allocate_Hessian_Spaces;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Coefficient_Circuits.Link_to_System;
                xr : in Standard_Floating_Vectors.Link_to_Vector;
                xi : in Standard_Floating_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    A,U,V : Standard_Complex_VecMats.VecMat(1..nbt);
    e : Standard_Complex_VecVecs.VecVec(1..nbt);
    ryd,iyd : Standard_Floating_VecVecs.VecVec(1..s.neq);
    pwtneeded : boolean := false;
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;
    hrps,hips : Standard_Floating_VecVecVecs.VecVecVec(1..nbt);

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start with the computation of the
    --   power table and then compute the SVD of the Hessian.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls : Standard_Complex_Vectors.Link_to_Vector;
      rlnk,ilnk,rpyd,ipyd : Standard_Floating_Vectors.Link_to_Vector;
      zr,zi,yr,yi,xrk,xik : double_float;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            rlnk := s.rpwt(mdx); ilnk := s.ipwt(mdx);
           -- lnk(1) := x(k)*x(k);
            xrk := xr(mdx); xik := xi(mdx);
            zr := xrk*xrk - xik*xik;
            zi := 2.0*xrk*xik;
            rlnk(1) := zr; ilnk(1) := zi;
            for i in 2..s.mxe(mdx)-1 loop
             -- lnk(i) := lnk(i-1)*x(k);
              yr := zr; yi := zi;
              zr := xrk*yr - xik*yi;
              zi := xrk*yi + xik*yr;
              rlnk(i) := zr; ilnk(i) := zi;
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(n,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); rpyd := ryd(mdx); ipyd := iyd(mdx);
        Standard_Coefficient_Circuits.Singular_Values
          (s.crc(mdx),xr,xi,rpyd,ipyd,s.rpwt,s.ipwt,hrps(i).all,hips(i).all,
           pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := Standard_Complex_Numbers.Create(rpyd(j),ipyd(j));
          rpyd(j) := 0.0; ipyd(j) := 0.0;
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        Standard_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start with the computation of the
    --   power table and then compute the SVD of the Hessian.
    --   Writes no output.

      mdx : integer32; -- local index
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls : Standard_Complex_Vectors.Link_to_Vector;
      rlnk,ilnk,rpyd,ipyd : Standard_Floating_Vectors.Link_to_Vector;
      zr,zi,yr,yi,xrk,xik : double_float;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            rlnk := s.rpwt(mdx); ilnk := s.ipwt(mdx);
           -- lnk(1) := x(k)*x(k);
            xrk := xr(mdx); xik := xi(mdx);
            zr := xrk*xrk - xik*xik;
            zi := 2.0*xrk*xik;
            rlnk(1) := zr; ilnk(1) := zi;
            for i in 2..s.mxe(mdx)-1 loop
             -- lnk(i) := lnk(i-1)*x(k);
              yr := zr; yi := zi;
              zr := xrk*yr - xik*yi;
              zi := xrk*yi + xik*yr;
              rlnk(i) := zr; ilnk(i) := zi;
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(n,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); rpyd := ryd(mdx); ipyd := iyd(mdx);
        Standard_Coefficient_Circuits.Singular_Values
          (s.crc(mdx),xr,xi,rpyd,ipyd,s.rpwt,s.ipwt,hrps(i).all,hips(i).all,
           pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := Standard_Complex_Numbers.Create(rpyd(j),ipyd(j));
          rpyd(j) := 0.0; ipyd(j) := 0.0;
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        Standard_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    A := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    U := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    V := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    e := Standard_Vector_Splitters.Allocate(nbt,s.dim,1,1);
   -- space for all gradients
    ryd := Standard_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    iyd := Standard_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    Allocate_Hessian_Spaces(s.dim,hrps,hips);
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
    Standard_Complex_VecMats.Clear(A);
    Standard_Complex_VecMats.Clear(U);
    Standard_Complex_VecMats.Clear(V);
    Standard_Complex_VecVecs.Clear(e);
    Standard_Floating_VecVecs.Clear(ryd);
    Standard_Floating_VecVecs.Clear(iyd);
    Standard_Floating_VecVecVecs.Clear(hrps);
    Standard_Floating_VecVecVecs.Clear(hips);
  end Multitasked_Singular_Values;

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out Standard_Complex_VecMats.VecMat;
                e : in out Standard_Complex_VecVecs.VecVec;
                yd : in out Standard_Complex_VecVecs.VecVec;
		verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;

    use Standard_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start with the computation of the
    --   power table and then compute the SVD of the Hessian.
    --   The index is i + k*n, for all k starting at 0,
    --   as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : Standard_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.neq loop
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(idx));
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.neq loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        Standard_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        Standard_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : Standard_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.neq loop
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.neq loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        Standard_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        Standard_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Static_Singular_Values;

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out Standard_Complex_VecMats.VecMat;
                e : in out Standard_Complex_VecVecs.VecVec;
                yd : in out Standard_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    use Standard_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start with the computation of the
    --   power table and then compute the SVD of the Hessian.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : Standard_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(n,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        Standard_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        Standard_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start with the computation of the power
    --   table and then compute the SVD of the Hessian.

      mdx : integer32;  -- local index
      pA,pU,pV : Standard_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : Standard_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        Standard_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := Standard_Complex_Numbers.Create(0.0);
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        Standard_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Dynamic_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Complex_Circuits.Link_to_System;
                x : in Standard_Complex_Vectors.Link_to_Vector;
                values : in out Standard_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false ) is

    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    A,U,V : Standard_Complex_VecMats.VecMat(1..nbt); -- space for SVD
    e : Standard_Complex_VecVecs.VecVec(1..nbt);
    yd : Standard_Complex_VecVecs.VecVec(1..s.neq); -- gradient space

  begin
    A := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    U := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    V := Standard_Complex_Circuits.Allocate(nbt,s.dim);
    e := Standard_Vector_Splitters.Allocate(nbt,s.dim,1,1);
    yd := Standard_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    if static then
      Static_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    else
      Dynamic_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    end if;
    Standard_Complex_VecMats.Clear(A);
    Standard_Complex_VecMats.Clear(U);
    Standard_Complex_VecMats.Clear(V);
    Standard_Complex_VecVecs.Clear(e);
    Standard_Complex_VecVecs.Clear(yd);
  end Multitasked_Singular_Values;

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DoblDobl_Complex_VecMats.VecMat;
                e : in out DoblDobl_Complex_VecVecs.VecVec;
                yd : in out DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;

    use DoblDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : DoblDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.neq loop
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(idx));
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        DoblDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := DoblDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        DoblDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.

      idx : integer32 := i; 
      pA,pU,pV : DoblDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        DoblDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := DoblDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        DoblDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Static_Singular_Values;

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DoblDobl_Complex_VecMats.VecMat;
                e : in out DoblDobl_Complex_VecVecs.VecVec;
                yd : in out DoblDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    use DoblDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and then compute the SVD of the Hessian.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : DoblDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        DoblDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := DoblDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        DoblDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.

      mdx : integer32; -- local index
      pA,pU,pV : DoblDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.dim loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        DoblDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := DoblDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        DoblDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Dynamic_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in DoblDobl_Complex_Circuits.Link_to_System;
                x : in DoblDobl_Complex_Vectors.Link_to_Vector;
                values : in out DoblDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false ) is

    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    A,U,V : DoblDobl_Complex_VecMats.VecMat(1..nbt); -- space for SVD
    e : DoblDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : DoblDobl_Complex_VecVecs.VecVec(1..s.neq); -- gradient space

  begin
    A := DoblDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := DoblDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := DoblDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := DoblDobl_Vector_Splitters.Allocate(nbt,s.dim,1,1);
    yd := DoblDobl_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    if static then
      Static_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    else
      Dynamic_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    end if;
    DoblDobl_Complex_VecMats.Clear(A);
    DoblDobl_Complex_VecMats.Clear(U);
    DoblDobl_Complex_VecMats.Clear(V);
    DoblDobl_Complex_VecVecs.Clear(e);
    DoblDobl_Complex_VecVecs.Clear(yd);
  end Multitasked_Singular_Values;

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                values : in out TripDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out TripDobl_Complex_VecMats.VecMat;
                e : in out TripDobl_Complex_VecVecs.VecVec;
                yd : in out TripDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;

    use TripDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : TripDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : TripDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.neq loop
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(idx));
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        TripDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := TripDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        TripDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.

      idx : integer32 := i; 
      pA,pU,pV : TripDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : TripDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        TripDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := TripDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        TripDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Static_Singular_Values;

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                values : in out TripDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out TripDobl_Complex_VecMats.VecMat;
                e : in out TripDobl_Complex_VecVecs.VecVec;
                yd : in out TripDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    use TripDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and then compute the SVD of the Hessian.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : TripDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : TripDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        TripDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := TripDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        TripDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   The first task computes the Jacobian and its singular values.

      mdx : integer32; -- local index
      pA,pU,pV : TripDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : TripDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.dim loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        TripDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := TripDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        TripDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Dynamic_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in TripDobl_Complex_Circuits.Link_to_System;
                x : in TripDobl_Complex_Vectors.Link_to_Vector;
                values : in out TripDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false ) is

    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    A,U,V : TripDobl_Complex_VecMats.VecMat(1..nbt); -- space for SVD
    e : TripDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : TripDobl_Complex_VecVecs.VecVec(1..s.neq); -- gradient space

  begin
    A := TripDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := TripDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := TripDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := TripDobl_Vector_Splitters.Allocate(nbt,s.dim,1,1);
    yd := TripDobl_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    if static then
      Static_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    else
      Dynamic_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    end if;
    TripDobl_Complex_VecMats.Clear(A);
    TripDobl_Complex_VecMats.Clear(U);
    TripDobl_Complex_VecMats.Clear(V);
    TripDobl_Complex_VecVecs.Clear(e);
    TripDobl_Complex_VecVecs.Clear(yd);
  end Multitasked_Singular_Values;

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out QuadDobl_Complex_VecMats.VecMat;
                e : in out QuadDobl_Complex_VecVecs.VecVec;
                yd : in out QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;

    use QuadDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(idx));
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        QuadDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := QuadDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        QuadDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        QuadDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := QuadDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        QuadDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Static_Singular_Values;

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out QuadDobl_Complex_VecMats.VecMat;
                e : in out QuadDobl_Complex_VecVecs.VecVec;
                yd : in out QuadDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    use QuadDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and then compute the SVD of Hessians.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.dim loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        QuadDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := QuadDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        QuadDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and compute SVD of Hessians, without output.

      mdx : integer32; 
      pA,pU,pV : QuadDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        QuadDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := QuadDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        QuadDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Dynamic_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in QuadDobl_Complex_Circuits.Link_to_System;
                x : in QuadDobl_Complex_Vectors.Link_to_Vector;
                values : in out QuadDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false ) is

    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    A,U,V : QuadDobl_Complex_VecMats.VecMat(1..nbt); -- space for SVD
    e : QuadDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : QuadDobl_Complex_VecVecs.VecVec(1..s.neq); -- gradient space

  begin
    A := QuadDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := QuadDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := QuadDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := QuadDobl_Vector_Splitters.Allocate(nbt,s.dim,1,1);
    yd := QuadDobl_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    if static then
      Static_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    else
      Dynamic_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    end if;
    QuadDobl_Complex_VecMats.Clear(A);
    QuadDobl_Complex_VecMats.Clear(U);
    QuadDobl_Complex_VecMats.Clear(V);
    QuadDobl_Complex_VecVecs.Clear(e);
    QuadDobl_Complex_VecVecs.Clear(yd);
  end Multitasked_Singular_Values;

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                values : in out PentDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out PentDobl_Complex_VecMats.VecMat;
                e : in out PentDobl_Complex_VecVecs.VecVec;
                yd : in out PentDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;

    use PentDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : PentDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : PentDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(idx));
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        PentDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := PentDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        PentDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 
      pA,pU,pV : PentDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : PentDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        PentDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := PentDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        PentDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Static_Singular_Values;

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                values : in out PentDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out PentDobl_Complex_VecMats.VecMat;
                e : in out PentDobl_Complex_VecVecs.VecVec;
                yd : in out PentDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    use PentDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and then compute the SVD of Hessians.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : PentDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : PentDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.dim loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        PentDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := PentDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        PentDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and compute SVD of Hessians, without output.

      mdx : integer32; 
      pA,pU,pV : PentDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : PentDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        PentDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := PentDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        PentDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Dynamic_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in PentDobl_Complex_Circuits.Link_to_System;
                x : in PentDobl_Complex_Vectors.Link_to_Vector;
                values : in out PentDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false ) is

    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    A,U,V : PentDobl_Complex_VecMats.VecMat(1..nbt); -- space for SVD
    e : PentDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : PentDobl_Complex_VecVecs.VecVec(1..s.neq); -- gradient space

  begin
    A := PentDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := PentDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := PentDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := PentDobl_Vector_Splitters.Allocate(nbt,s.dim,1,1);
    yd := PentDobl_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    if static then
      Static_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    else
      Dynamic_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    end if;
    PentDobl_Complex_VecMats.Clear(A);
    PentDobl_Complex_VecMats.Clear(U);
    PentDobl_Complex_VecMats.Clear(V);
    PentDobl_Complex_VecVecs.Clear(e);
    PentDobl_Complex_VecVecs.Clear(yd);
  end Multitasked_Singular_Values;

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                values : in out OctoDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out OctoDobl_Complex_VecMats.VecMat;
                e : in out OctoDobl_Complex_VecVecs.VecVec;
                yd : in out OctoDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;

    use OctoDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : OctoDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : OctoDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(idx));
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        OctoDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := OctoDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        OctoDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 
      pA,pU,pV : OctoDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : OctoDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        OctoDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := OctoDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        OctoDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Static_Singular_Values;

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                values : in out OctoDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out OctoDobl_Complex_VecMats.VecMat;
                e : in out OctoDobl_Complex_VecVecs.VecVec;
                yd : in out OctoDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    use OctoDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and then compute the SVD of Hessians.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : OctoDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : OctoDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.dim loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        OctoDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := OctoDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        OctoDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and compute SVD of Hessians, without output.

      mdx : integer32; 
      pA,pU,pV : OctoDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : OctoDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        OctoDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := OctoDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        OctoDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Dynamic_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in OctoDobl_Complex_Circuits.Link_to_System;
                x : in OctoDobl_Complex_Vectors.Link_to_Vector;
                values : in out OctoDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false ) is

    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    A,U,V : OctoDobl_Complex_VecMats.VecMat(1..nbt); -- space for SVD
    e : OctoDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : OctoDobl_Complex_VecVecs.VecVec(1..s.neq); -- gradient space

  begin
    A := OctoDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := OctoDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := OctoDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := OctoDobl_Vector_Splitters.Allocate(nbt,s.dim,1,1);
    yd := OctoDobl_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    if static then
      Static_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    else
      Dynamic_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    end if;
    OctoDobl_Complex_VecMats.Clear(A);
    OctoDobl_Complex_VecMats.Clear(U);
    OctoDobl_Complex_VecMats.Clear(V);
    OctoDobl_Complex_VecVecs.Clear(e);
    OctoDobl_Complex_VecVecs.Clear(yd);
  end Multitasked_Singular_Values;

  procedure Static_Singular_Values
              ( nbt : in integer32;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                values : in out DecaDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DecaDobl_Complex_VecMats.VecMat;
                e : in out DecaDobl_Complex_VecVecs.VecVec;
                yd : in out DecaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;

    use DecaDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.
    --   Writes one line for each job to screen.

      idx : integer32 := i; 
      pA,pU,pV : DecaDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DecaDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(idx));
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(idx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        DecaDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := DecaDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        DecaDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 
      pA,pU,pV : DecaDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DecaDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx <= s.dim loop
          if s.mxe(idx) > 1 then
            xpw := s.pwt(idx);
            xpw(1) := x(idx)*x(idx);
            for k in 2..(s.mxe(idx)-1) loop
              xpw(k) := xpw(k-1)*x(idx);
            end loop;
          end if;
          idx := idx + n;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
        idx := i; -- reset the index for Hessian and SVD computation
      end if;
      while idx <= s.dim loop
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(idx); pyd := yd(idx);
        DecaDobl_Complex_Circuits.Singular_Values
          (s.crc(idx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(idx,j) := pyd(j);
          pyd(j) := DecaDobl_Complex_Numbers.Create(integer(0));
        end loop;
        idx := idx + n;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        DecaDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Static_Singular_Values;

  procedure Dynamic_Singular_Values
              ( nbt : in integer32;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                values : in out DecaDobl_Complex_VecVecs.VecVec;
                pwtdone,gradone : in out Multitasking.boolean_array;
                A,U,V : in out DecaDobl_Complex_VecMats.VecMat;
                e : in out DecaDobl_Complex_VecVecs.VecVec;
                yd : in out DecaDobl_Complex_VecVecs.VecVec;
                verbose : in boolean := false ) is

    pwtneeded : boolean := false;
    info : integer32;
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    use DecaDobl_Complex_Numbers;

    procedure Reporting_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and then compute the SVD of Hessians.
    --   Writes one line for each job to screen.

      mdx : integer32; -- local index
      pA,pU,pV : DecaDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DecaDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.dim loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          put_line("task " & Multitasking.to_string(i)
                           & " computes power table row "
                           & Multitasking.to_string(mdx));
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.dim loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        put_line("task " & Multitasking.to_string(i)
                         & " computes gradient & Hessian "
                         & Multitasking.to_string(mdx));
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        DecaDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := DecaDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        put_line("task " & Multitasking.to_string(i)
                         & " computes SVD of Jacobian");
        DecaDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Reporting_Job;
    procedure report_jobs is new Multitasking.Silent_Workers(Reporting_Job);

    procedure Silent_Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will start the computation of the power table
    --   and compute SVD of Hessians, without output.

      mdx : integer32; 
      pA,pU,pV : DecaDobl_Complex_Matrices.Link_to_Matrix;
      pe,vls,pyd,xpw : DecaDobl_Complex_Vectors.Link_to_Vector;

    begin
      if pwtneeded then
        while idx1 <= s.neq loop
          Semaphore.Request(lck1);
          if idx1 <= s.neq
           then idx1 := idx1 + 1; mdx := idx1;
          end if;
          Semaphore.Release(lck1);
          exit when (mdx > s.neq);
          if s.mxe(mdx) > 1 then
            xpw := s.pwt(mdx);
            xpw(1) := x(mdx)*x(mdx);
            for k in 2..(s.mxe(mdx)-1) loop
              xpw(k) := xpw(k-1)*x(mdx);
            end loop;
          end if;
        end loop;
        pwtdone(i) := true;
       -- make sure all tasks are done with the power table
        while not Multitasking.all_true(nbt,pwtdone) loop
          delay 0.001;
        end loop;
      end if;
      while idx2 <= s.neq loop
        Semaphore.Request(lck2);
        if idx2 <= s.neq
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > s.neq);
        pA := A(i); pU := U(i); pV := V(i); pe := e(i);
        vls := values(mdx); pyd := yd(mdx);
        DecaDobl_Complex_Circuits.Singular_Values
          (s.crc(mdx),x,pyd,s.pwt,pA.all,pU.all,pV.all,pe.all,vls.all);
        for j in s.jm'range(2) loop
          s.jm(mdx,j) := pyd(j);
          pyd(j) := DecaDobl_Complex_Numbers.Create(integer(0));
        end loop;
      end loop;
      gradone(i) := true;
      if i = n then
        while not Multitasking.all_true(nbt,gradone) loop
          delay 0.001;
        end loop;
        DecaDobl_Complex_Singular_Values.SVD
          (s.jm,s.dim,s.dim,values(0).all,pe.all,pU.all,pV.all,11,info);
      end if;
    end Silent_Job;
    procedure silent_jobs is new Multitasking.Silent_Workers(Silent_Job);

  begin
    for k in s.mxe'range loop -- determine if power table is needed
      if s.mxe(k) > 1
       then pwtneeded := true; exit;
      end if;
    end loop;
    if verbose
     then report_jobs(nbt);
     else silent_jobs(nbt);
    end if;
  end Dynamic_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in DecaDobl_Complex_Circuits.Link_to_System;
                x : in DecaDobl_Complex_Vectors.Link_to_Vector;
                values : in out DecaDobl_Complex_VecVecs.VecVec;
                static : in boolean := false;
                verbose : in boolean := false ) is

    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    gradone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    A,U,V : DecaDobl_Complex_VecMats.VecMat(1..nbt); -- space for SVD
    e : DecaDobl_Complex_VecVecs.VecVec(1..nbt);
    yd : DecaDobl_Complex_VecVecs.VecVec(1..s.neq); -- gradient space

  begin
    A := DecaDobl_Complex_Circuits.Allocate(nbt,s.dim);
    U := DecaDobl_Complex_Circuits.Allocate(nbt,s.dim);
    V := DecaDobl_Complex_Circuits.Allocate(nbt,s.dim);
    e := DecaDobl_Vector_Splitters.Allocate(nbt,s.dim,1,1);
    yd := DecaDobl_Vector_Splitters.Allocate(s.neq,s.dim,1,0);
    if static then
      Static_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    else
      Dynamic_Singular_Values
        (nbt,s,x,values,pwtdone,gradone,A,U,V,e,yd,verbose);
    end if;
    DecaDobl_Complex_VecMats.Clear(A);
    DecaDobl_Complex_VecMats.Clear(U);
    DecaDobl_Complex_VecMats.Clear(V);
    DecaDobl_Complex_VecVecs.Clear(e);
    DecaDobl_Complex_VecVecs.Clear(yd);
  end Multitasked_Singular_Values;

  function Standard_Distance
              ( values : Standard_Complex_VecVecs.VecVec )
              return double_float is

    jmsvls : constant Standard_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant double_float -- jmsvls'range = 1..dim+1 = jmsvls'last
           := Standard_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last-1));
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
    sigma1 : constant double_double -- jmsvls'range = 1..dim+1 = jmsvls'last
           := DoblDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last-1));
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

  function TripDobl_Distance
              ( values : TripDobl_Complex_VecVecs.VecVec )
              return triple_double is

    jmsvls : constant TripDobl_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant triple_double -- jmsvls'range = 1..dim+1 = jmsvls'last
           := TripDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last-1));
    lnk : TripDobl_Complex_Vectors.Link_to_Vector;
    sumvls : triple_double := create(0.0);
    val,nrm : triple_double;

  begin
    for k in values'range loop
      lnk := values(k);
      val := TripDobl_Complex_Numbers.REAL_PART(lnk(lnk'first));
      sumvls := sumvls + val*val;
    end loop;
    nrm := TripDobl_Mathematical_Functions.SQRT(sumvls);
    return (2.0*sigma1)/nrm;
  end TripDobl_Distance;

  function QuadDobl_Distance
              ( values : QuadDobl_Complex_VecVecs.VecVec )
              return quad_double is

    jmsvls : constant QuadDobl_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant quad_double -- jmsvls'range = 1..dim+1 = jmsvls'last
           := QuadDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last-1));
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

  function PentDobl_Distance
              ( values : PentDobl_Complex_VecVecs.VecVec )
              return penta_double is

    jmsvls : constant PentDobl_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant penta_double -- jmsvls'range = 1..dim+1 = jmsvls'last
           := PentDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last-1));
    lnk : PentDobl_Complex_Vectors.Link_to_Vector;
    sumvls : penta_double := create(0.0);
    val,nrm : penta_double;

  begin
    for k in values'range loop
      lnk := values(k);
      val := PentDobl_Complex_Numbers.REAL_PART(lnk(lnk'first));
      sumvls := sumvls + val*val;
    end loop;
    nrm := PentDobl_Mathematical_Functions.SQRT(sumvls);
    return (2.0*sigma1)/nrm;
  end PentDobl_Distance;

  function OctoDobl_Distance
              ( values : OctoDobl_Complex_VecVecs.VecVec )
              return octo_double is

    jmsvls : constant OctoDobl_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant octo_double -- jmsvls'range = 1..dim+1 = jmsvls'last
           := OctoDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last-1));
    lnk : OctoDobl_Complex_Vectors.Link_to_Vector;
    sumvls : octo_double := create(0.0);
    val,nrm : octo_double;

  begin
    for k in values'range loop
      lnk := values(k);
      val := OctoDobl_Complex_Numbers.REAL_PART(lnk(lnk'first));
      sumvls := sumvls + val*val;
    end loop;
    nrm := OctoDobl_Mathematical_Functions.SQRT(sumvls);
    return (2.0*sigma1)/nrm;
  end OctoDobl_Distance;

  function DecaDobl_Distance
              ( values : DecaDobl_Complex_VecVecs.VecVec )
              return deca_double is

    jmsvls : constant DecaDobl_Complex_Vectors.Link_to_Vector := values(0);
    sigma1 : constant deca_double -- jmsvls'range = 1..dim+1 = jmsvls'last
           := DecaDobl_Complex_Numbers.REAL_PART(jmsvls(jmsvls'last-1));
    lnk : DecaDobl_Complex_Vectors.Link_to_Vector;
    sumvls : deca_double := create(0.0);
    val,nrm : deca_double;

  begin
    for k in values'range loop
      lnk := values(k);
      val := DecaDobl_Complex_Numbers.REAL_PART(lnk(lnk'first));
      sumvls := sumvls + val*val;
    end loop;
    nrm := DecaDobl_Mathematical_Functions.SQRT(sumvls);
    return (2.0*sigma1)/nrm;
  end DecaDobl_Distance;

end Multitasked_Hessian_Circuits;
