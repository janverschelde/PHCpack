with text_io;                            use text_io;
with Standard_Complex_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with PentDobl_Complex_Vectors;
with PentDobl_Complex_Matrices;
with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Matrices;
with DecaDobl_Complex_Vectors;
with DecaDobl_Complex_Matrices;
with HexaDobl_Complex_Vectors;
with HexaDobl_Complex_Matrices;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with Semaphore,Multitasking;

package body Multitasked_AlgoDiff_Convolutions is

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return Standard_Floating_VecVecVecs.VecVecVec is

    use Standard_Floating_VecVecVecs;
    use Standard_Vector_Splitters;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant Standard_Floating_VecVecs.VecVec(1..dim+1)
            := Allocate_Floating_Coefficients(dim+1,deg);
      begin
        res(i) := new Standard_Floating_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return Standard_Speelpenning_Convolutions.VecVecVec is

    use Standard_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new Standard_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return DoblDobl_Speelpenning_Convolutions.VecVecVec is

    use DoblDobl_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new DoblDobl_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return TripDobl_Speelpenning_Convolutions.VecVecVec is

    use TripDobl_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant TripDobl_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new TripDobl_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return QuadDobl_Speelpenning_Convolutions.VecVecVec is

    use QuadDobl_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new QuadDobl_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return PentDobl_Speelpenning_Convolutions.VecVecVec is

    use PentDobl_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant PentDobl_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new PentDobl_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return OctoDobl_Speelpenning_Convolutions.VecVecVec is

    use OctoDobl_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant OctoDobl_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new OctoDobl_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return DecaDobl_Speelpenning_Convolutions.VecVecVec is

    use DecaDobl_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant DecaDobl_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new DecaDobl_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return HexaDobl_Speelpenning_Convolutions.VecVecVec is

    use HexaDobl_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant HexaDobl_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new HexaDobl_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  procedure Standard_Static_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Coefficient_Convolutions.Circuits;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                output : in boolean := false ) is

  -- DESCRIPTION :
  --   Original procedure to define Standard_Multitasked_EvalDiff
  --   with static load balancing.  Except for the flag "static",
  --   the specification is the same as Standard_Multitasked_EvalDiff.

    use Standard_Floating_VecVecVecs;
    use Standard_Coefficient_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    ryd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    iyd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : Standard_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          rxpw := rpwt(idx); ixpw := ipwt(idx);
          Multiply(rx(idx),ix(idx),rx(idx),ix(idx),rxpw(1),ixpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),rx(idx),ix(idx),rxpw(k),ixpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,rx.all,ix.all,rpwt,ipwt,ryd(i).all,iyd(i).all);
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(rx'last+1); ivright := iydi(ix'last+1);
        for j in rvright'range loop
          vleft := vy(j);
          vleft(idx) := Standard_Complex_Numbers.Create(rvright(j),ivright(j));
          rvright(j) := 0.0; ivright(j) := 0.0;
        end loop;
        for j in 1..rx'last loop
          rvright := rydi(j); ivright := iydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j)               -- column j in vm(k) is the variable
              := Standard_Complex_Numbers.Create(rvright(k),ivright(k));
            rvright(k) := 0.0; ivright(k) := 0.0;
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : Standard_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          rxpw := rpwt(idx); ixpw := ipwt(idx);
          Multiply(rx(idx),ix(idx),rx(idx),ix(idx),rxpw(1),ixpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),rx(idx),ix(idx),rxpw(k),ixpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,rx.all,ix.all,rpwt,ipwt,ryd(i).all,iyd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(rx'last+1); ivright := iydi(ix'last+1);
        for j in rvright'range loop
          vleft := vy(j);
          vleft(idx) := Standard_Complex_Numbers.Create(rvright(j),ivright(j));
          rvright(j) := 0.0; ivright(j) := 0.0;
        end loop;
        for j in 1..rx'last loop
          rvright := rydi(j); ivright := iydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j)               -- column j in vm(k) is the variable
              := Standard_Complex_Numbers.Create(rvright(k),ivright(k));
            rvright(k) := 0.0; ivright(k) := 0.0;
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end Standard_Static_EvalDiff;

  procedure Standard_Dynamic_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Coefficient_Convolutions.Circuits;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                output : in boolean := false ) is

  -- DESCRIPTION :
  --   New procedure to define Standard_Multitasked_EvalDiff
  --   with dynamic load balancing.  Except for the flag "static",
  --   the specification is the same as Standard_Multitasked_EvalDiff.

    use Standard_Floating_VecVecVecs;
    use Standard_Coefficient_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    ryd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    iyd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : Standard_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
      mdx : integer32; -- local index
 
    begin
      while idx1 <= c'last loop
        Semaphore.Request(lck1);
        if idx1 <= c'last
         then idx1 := idx1 + 1; mdx := idx1;
        end if;
        Semaphore.Release(lck1);
        exit when (mdx > c'last);
        if mxe(mdx) > 2 then
          rxpw := rpwt(mdx); ixpw := ipwt(mdx);
          Multiply(rx(mdx),ix(mdx),rx(mdx),ix(mdx),rxpw(1),ixpw(1));
          for k in 2..(mxe(mdx)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),rx(mdx),ix(mdx),rxpw(k),ixpw(k));
          end loop;
        end if;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(n,pwtdone) loop
        delay 0.001;
      end loop;
      while idx2 <= c'last loop
        Semaphore.Request(lck2);
        if idx2 <= c'last
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > c'last);
        EvalDiff(c(mdx).all,rx.all,ix.all,rpwt,ipwt,ryd(i).all,iyd(i).all);
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(rx'last+1); ivright := iydi(ix'last+1);
        for j in rvright'range loop
          vleft := vy(j);
          vleft(mdx) := Standard_Complex_Numbers.Create(rvright(j),ivright(j));
          rvright(j) := 0.0; ivright(j) := 0.0;
        end loop;
        for j in 1..rx'last loop
          rvright := rydi(j); ivright := iydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(mdx,j)              -- column j in vm(k) is the variable
              := Standard_Complex_Numbers.Create(rvright(k),ivright(k));
            rvright(k) := 0.0; ivright(k) := 0.0;
          end loop;
        end loop;

      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : Standard_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
      mdx : integer32; -- local index
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx1 <= c'last loop
        Semaphore.Request(lck1);
        if idx1 > c'last then
          put_line("task " & Multitasking.to_string(i)
                           & " leaves the loop");
        else
          put_line("task " & Multitasking.to_string(i)
                           & " increments idx "
                           & Multitasking.to_string(idx1));
          idx1 := idx1 + 1; mdx := idx1;
          if idx1 <= c'last then
            put_line("task " & Multitasking.to_string(i)
                             & " will do job "
                             & Multitasking.to_string(idx1));
          else
            put_line("task " & Multitasking.to_string(i)
                             & " leaves the loop");
          end if;
        end if;
        Semaphore.Release(lck1);
        exit when (mdx > c'last);
        if mxe(mdx) > 2 then
          rxpw := rpwt(mdx); ixpw := ipwt(mdx);
          Multiply(rx(mdx),ix(mdx),rx(mdx),ix(mdx),rxpw(1),ixpw(1));
          for k in 2..(mxe(idx1)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),rx(mdx),ix(mdx),rxpw(k),ixpw(k));
          end loop;
        end if;
      end loop;
      put_line("task " & Multitasking.to_string(i)
                       & " is done with first loop");
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(n,pwtdone) loop
        delay 0.001;
      end loop;
      while idx2 <= c'last loop
        Semaphore.Request(lck2);
        if idx2 > c'last then
          put_line("task " & Multitasking.to_string(i)
                           & " leaves the loop");
        else
          put_line("task " & Multitasking.to_string(i)
                           & " increments idx "
                           & Multitasking.to_string(idx2));
          idx2 := idx2 + 1; mdx := idx2;
          if idx2 <= c'last then
            put_line("task " & Multitasking.to_string(i)
                             & " will do job "
                             & Multitasking.to_string(idx2));
          else
            put_line("task " & Multitasking.to_string(i)
                             & " leaves the loop");
          end if;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > c'last);
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(mdx));
        EvalDiff(c(mdx).all,rx.all,ix.all,rpwt,ipwt,ryd(i).all,iyd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(mdx));
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(rx'last+1); ivright := iydi(ix'last+1);
        for j in rvright'range loop
          vleft := vy(j);
          vleft(mdx) := Standard_Complex_Numbers.Create(rvright(j),ivright(j));
          rvright(j) := 0.0; ivright(j) := 0.0;
        end loop;
        for j in 1..rx'last loop
          rvright := rydi(j); ivright := iydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(mdx,j)               -- column j in vm(k) is the variable
              := Standard_Complex_Numbers.Create(rvright(k),ivright(k));
            rvright(k) := 0.0; ivright(k) := 0.0;
          end loop;
        end loop;
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
  end Standard_Dynamic_EvalDiff;

  procedure Standard_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Coefficient_Convolutions.Circuits;
                rx,ix : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                static : in boolean := false;
                output : in boolean := false ) is
  begin
    if static
     then Standard_Static_EvalDiff(nbt,c,rx,ix,mxe,rpwt,ipwt,vy,vm,output);
     else Standard_Dynamic_EvalDiff(nbt,c,rx,ix,mxe,rpwt,ipwt,vy,vm,output);
    end if;
  end Standard_Multitasked_EvalDiff;

  procedure DoblDobl_Static_EvalDiff
              ( nbt : in integer32;
                c : in DoblDobl_Coefficient_Convolutions.Circuits;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

  -- DESCRIPTION :
  --   Original procedure to define DoblDobl_Multitasked_EvalDiff
  --   with static load balancing.  Except for the flag "static",
  --   the specification is the same as DoblDobl_Multitasked_EvalDiff.

    use DoblDobl_Vector_Splitters;
    use DoblDobl_Coefficient_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    rhyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    ihyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    rlyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    ilyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      rhydi,ihydi : Standard_Floating_VecVecs.Link_to_VecVec;
      rlydi,ilydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : DoblDobl_Complex_Vectors.Link_to_Vector;
      rhvright,ihvright : Standard_Floating_Vectors.Link_to_Vector;
      rlvright,ilvright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
      rhxpw,ihxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rlxpw,ilxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rdd,idd : double_double;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          rhxpw := rhpwt(idx); ihxpw := ihpwt(idx);
          rlxpw := rlpwt(idx); ilxpw := ilpwt(idx);
          Multiply(rhx(idx),ihx(idx),rlx(idx),ilx(idx),
                   rhx(idx),ihx(idx),rlx(idx),ilx(idx),
                   rhxpw(1),ihxpw(1),rlxpw(1),ilxpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(rhxpw(k-1),ihxpw(k-1),rlxpw(k-1),ilxpw(k-1),
                     rhx(idx),ihx(idx),rlx(idx),ilx(idx),
                     rhxpw(k),ihxpw(k),rlxpw(k),ilxpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,rhx.all,ihx.all,rlx.all,ilx.all,rhpwt,ihpwt,
                 rlpwt,ilpwt,rhyd(i).all,ihyd(i).all,rlyd(i).all,ilyd(i).all);
        rhydi := rhyd(i); ihydi := ihyd(i);
        rlydi := rlyd(i); ilydi := ilyd(i);
        rhvright := rhydi(rhx'last+1); ihvright := ihydi(ihx'last+1);
        rlvright := rlydi(rlx'last+1); ilvright := ilydi(ilx'last+1);
        for j in rhvright'range loop
          vleft := vy(j);
          rdd := Create(rhvright(j),rlvright(j));
          idd := Create(ihvright(j),ilvright(j));
          vleft(idx) := DoblDobl_Complex_Numbers.Create(rdd,idd);
          rhvright(j) := 0.0; ihvright(j) := 0.0;
          rlvright(j) := 0.0; ilvright(j) := 0.0;
        end loop;
        for j in 1..rhx'last loop
          rhvright := rhydi(j); ihvright := ihydi(j);
          rlvright := rlydi(j); ilvright := ilydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rdd := Create(rhvright(k),rlvright(k));
            idd := Create(ihvright(k),ilvright(k));
            mleft(idx,j) := DoblDobl_Complex_Numbers.Create(rdd,idd);
            rhvright(k) := 0.0; ihvright(k) := 0.0;
            rlvright(k) := 0.0; ilvright(k) := 0.0;
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      rhydi,ihydi : Standard_Floating_VecVecs.Link_to_VecVec;
      rlydi,ilydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : DoblDobl_Complex_Vectors.Link_to_Vector;
      rhvright,ihvright : Standard_Floating_Vectors.Link_to_Vector;
      rlvright,ilvright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
      rhxpw,ihxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rlxpw,ilxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rdd,idd : double_double;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          rhxpw := rhpwt(idx); ihxpw := ihpwt(idx);
          rlxpw := rlpwt(idx); ilxpw := ilpwt(idx);
          Multiply(rhx(idx),ihx(idx),rlx(idx),ilx(idx),
                   rhx(idx),ihx(idx),rlx(idx),ilx(idx),
                   rhxpw(1),ihxpw(1),rlxpw(1),ilxpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(rhxpw(k-1),ihxpw(k-1),rlxpw(k-1),ilxpw(k-1),
                     rhx(idx),ihx(idx),rlx(idx),ilx(idx),
                     rhxpw(k),ihxpw(k),rlxpw(k),ilxpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,rhx.all,ihx.all,rlx.all,ilx.all,rhpwt,ihpwt,
                 rlpwt,ilpwt,rhyd(i).all,ihyd(i).all,rlyd(i).all,ilyd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        rhydi := rhyd(i); ihydi := ihyd(i);
        rlydi := rlyd(i); ilydi := ilyd(i);
        rhvright := rhydi(rhx'last+1); ihvright := ihydi(ihx'last+1);
        rlvright := rlydi(rlx'last+1); ilvright := ilydi(ilx'last+1);
        for j in rhvright'range loop
          vleft := vy(j);
          rdd := Create(rhvright(j),rlvright(j));
          idd := Create(ihvright(j),ilvright(j));
          vleft(idx) := DoblDobl_Complex_Numbers.Create(rdd,idd);
          rhvright(j) := 0.0; ihvright(j) := 0.0;
          rlvright(j) := 0.0; ilvright(j) := 0.0;
        end loop;
        for j in 1..rhx'last loop
          rhvright := rhydi(j); ihvright := ihydi(j);
          rlvright := rlydi(j); ilvright := ilydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rdd := Create(rhvright(k),rlvright(k));
            idd := Create(ihvright(k),ilvright(k));
            mleft(idx,j) := DoblDobl_Complex_Numbers.Create(rdd,idd);
            rhvright(k) := 0.0; ihvright(k) := 0.0;
            rlvright(k) := 0.0; ilvright(k) := 0.0;
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end DoblDobl_Static_EvalDiff;

  procedure DoblDobl_Dynamic_EvalDiff
              ( nbt : in integer32;
                c : in DoblDobl_Coefficient_Convolutions.Circuits;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

  -- DESCRIPTION :
  --   New procedure to define DoblDobl_Multitasked_EvalDiff
  --   with dynamic load balancing.  Except for the flag "static",
  --   the specification is the same as DoblDobl_Multitasked_EvalDiff.

    use DoblDobl_Vector_Splitters;
    use DoblDobl_Coefficient_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    rhyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    ihyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    rlyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    ilyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
         := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    idx1,idx2 : integer32 := 0; -- global indices
    lck1,lck2 : Semaphore.Lock;

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      mdx : integer32; -- local index
      rhydi,ihydi : Standard_Floating_VecVecs.Link_to_VecVec;
      rlydi,ilydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : DoblDobl_Complex_Vectors.Link_to_Vector;
      rhvright,ihvright : Standard_Floating_Vectors.Link_to_Vector;
      rlvright,ilvright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
      rhxpw,ihxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rlxpw,ilxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rdd,idd : double_double;
 
    begin
      while idx1 <= c'last loop
        Semaphore.Request(lck1);
        if idx1 <= c'last
         then idx1 := idx1 + 1; mdx := idx1;
        end if;
        Semaphore.Release(lck1);
        exit when (mdx > c'last);
        if mxe(mdx) > 2 then
          rhxpw := rhpwt(mdx); ihxpw := ihpwt(mdx);
          rlxpw := rlpwt(mdx); ilxpw := ilpwt(mdx);
          Multiply(rhx(mdx),ihx(mdx),rlx(mdx),ilx(mdx),
                   rhx(mdx),ihx(mdx),rlx(mdx),ilx(mdx),
                   rhxpw(1),ihxpw(1),rlxpw(1),ilxpw(1));
          for k in 2..(mxe(mdx)-2) loop
            Multiply(rhxpw(k-1),ihxpw(k-1),rlxpw(k-1),ilxpw(k-1),
                     rhx(mdx),ihx(mdx),rlx(mdx),ilx(mdx),
                     rhxpw(k),ihxpw(k),rlxpw(k),ilxpw(k));
          end loop;
        end if;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(n,pwtdone) loop
        delay 0.001;
      end loop;
      while idx2 <= c'last loop
        Semaphore.Request(lck2);
        if idx2 <= c'last
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > c'last);
        EvalDiff(c(mdx).all,rhx.all,ihx.all,rlx.all,ilx.all,rhpwt,ihpwt,
                 rlpwt,ilpwt,rhyd(i).all,ihyd(i).all,rlyd(i).all,ilyd(i).all);
        rhydi := rhyd(i); ihydi := ihyd(i);
        rlydi := rlyd(i); ilydi := ilyd(i);
        rhvright := rhydi(rhx'last+1); ihvright := ihydi(ihx'last+1);
        rlvright := rlydi(rlx'last+1); ilvright := ilydi(ilx'last+1);
        for j in rhvright'range loop
          vleft := vy(j);
          rdd := Create(rhvright(j),rlvright(j));
          idd := Create(ihvright(j),ilvright(j));
          vleft(mdx) := DoblDobl_Complex_Numbers.Create(rdd,idd);
          rhvright(j) := 0.0; ihvright(j) := 0.0;
          rlvright(j) := 0.0; ilvright(j) := 0.0;
        end loop;
        for j in 1..rhx'last loop
          rhvright := rhydi(j); ihvright := ihydi(j);
          rlvright := rlydi(j); ilvright := ilydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row mdx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rdd := Create(rhvright(k),rlvright(k));
            idd := Create(ihvright(k),ilvright(k));
            mleft(mdx,j) := DoblDobl_Complex_Numbers.Create(rdd,idd);
            rhvright(k) := 0.0; ihvright(k) := 0.0;
            rlvright(k) := 0.0; ilvright(k) := 0.0;
          end loop;
        end loop;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      mdx : integer32; -- local index
      rhydi,ihydi : Standard_Floating_VecVecs.Link_to_VecVec;
      rlydi,ilydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : DoblDobl_Complex_Vectors.Link_to_Vector;
      rhvright,ihvright : Standard_Floating_Vectors.Link_to_Vector;
      rlvright,ilvright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
      rhxpw,ihxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rlxpw,ilxpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rdd,idd : double_double;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx1 <= c'last loop
        Semaphore.Request(lck1);
        if idx1 > c'last then
          put_line("task " & Multitasking.to_string(i)
                           & " leaves the loop");
        else
          put_line("task " & Multitasking.to_string(i)
                           & " increments idx "
                           & Multitasking.to_string(idx1));
          idx1 := idx1 + 1; mdx := idx1;
          if idx1 <= c'last then
            put_line("task " & Multitasking.to_string(i)
                             & " will do job "
                             & Multitasking.to_string(idx1));
          else
            put_line("task " & Multitasking.to_string(i)
                             & " leaves the loop");
          end if;
        end if;
        Semaphore.Release(lck1);
        exit when (mdx > c'last);
        if mxe(mdx) > 2 then
          rhxpw := rhpwt(mdx); ihxpw := ihpwt(mdx);
          rlxpw := rlpwt(mdx); ilxpw := ilpwt(mdx);
          Multiply(rhx(mdx),ihx(mdx),rlx(mdx),ilx(mdx),
                   rhx(mdx),ihx(mdx),rlx(mdx),ilx(mdx),
                   rhxpw(1),ihxpw(1),rlxpw(1),ilxpw(1));
          for k in 2..(mxe(mdx)-2) loop
            Multiply(rhxpw(k-1),ihxpw(k-1),rlxpw(k-1),ilxpw(k-1),
                     rhx(mdx),ihx(mdx),rlx(mdx),ilx(mdx),
                     rhxpw(k),ihxpw(k),rlxpw(k),ilxpw(k));
          end loop;
        end if;
      end loop;
      put_line("task " & Multitasking.to_string(i)
                       & " is done with first loop");
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(n,pwtdone) loop
        delay 0.001;
      end loop;
      while idx2 <= c'last loop
        Semaphore.Request(lck2);
        if idx2 > c'last then
          put_line("task " & Multitasking.to_string(i)
                           & " leaves the loop");
        else
          put_line("task " & Multitasking.to_string(i)
                           & " increments idx "
                           & Multitasking.to_string(idx2));
          idx2 := idx2 + 1; mdx := idx2;
          if idx2 <= c'last then
            put_line("task " & Multitasking.to_string(i)
                             & " will do job "
                             & Multitasking.to_string(idx2));
          else
            put_line("task " & Multitasking.to_string(i)
                             & " leaves the loop");
          end if;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > c'last);
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(mdx));
        EvalDiff(c(mdx).all,rhx.all,ihx.all,rlx.all,ilx.all,rhpwt,ihpwt,
                 rlpwt,ilpwt,rhyd(i).all,ihyd(i).all,rlyd(i).all,ilyd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(mdx));
        rhydi := rhyd(i); ihydi := ihyd(i);
        rlydi := rlyd(i); ilydi := ilyd(i);
        rhvright := rhydi(rhx'last+1); ihvright := ihydi(ihx'last+1);
        rlvright := rlydi(rlx'last+1); ilvright := ilydi(ilx'last+1);
        for j in rhvright'range loop
          vleft := vy(j);
          rdd := Create(rhvright(j),rlvright(j));
          idd := Create(ihvright(j),ilvright(j));
          vleft(mdx) := DoblDobl_Complex_Numbers.Create(rdd,idd);
          rhvright(j) := 0.0; ihvright(j) := 0.0;
          rlvright(j) := 0.0; ilvright(j) := 0.0;
        end loop;
        for j in 1..rhx'last loop
          rhvright := rhydi(j); ihvright := ihydi(j);
          rlvright := rlydi(j); ilvright := ilydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row mdx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rdd := Create(rhvright(k),rlvright(k));
            idd := Create(ihvright(k),ilvright(k));
            mleft(mdx,j) := DoblDobl_Complex_Numbers.Create(rdd,idd);
            rhvright(k) := 0.0; ihvright(k) := 0.0;
            rlvright(k) := 0.0; ilvright(k) := 0.0;
          end loop;
        end loop;
      end loop;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
  end DoblDobl_Dynamic_EvalDiff;

  procedure DoblDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in DoblDobl_Coefficient_Convolutions.Circuits;
                rhx,ihx : in Standard_Floating_VecVecs.Link_to_VecVec;
                rlx,ilx : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rhpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                static : in boolean := false;
                output : in boolean := false ) is

  begin
    if static then
      DoblDobl_Static_EvalDiff
        (nbt,c,rhx,ihx,rlx,ilx,mxe,rhpwt,ihpwt,rlpwt,ilpwt,vy,vm,output);
    else
      DoblDobl_Dynamic_EvalDiff
        (nbt,c,rhx,ihx,rlx,ilx,mxe,rhpwt,ihpwt,rlpwt,ilpwt,vy,vm,output);
    end if;
  end DoblDobl_Multitasked_EvalDiff;

  procedure QuadDobl_Static_EvalDiff
              ( nbt : in integer32;
                c : in QuadDobl_Coefficient_Convolutions.Circuits;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

  -- DESCRIPTION :
  --   Original procedure to define QuadDobl_Multitasked_EvalDiff
  --   with static load balancing.  Except for the flag "static",
  --   the specification is the same as QuadDobl_Multitasked_EvalDiff.

    use QuadDobl_Vector_Splitters;
    use QuadDobl_Coefficient_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    degdim : constant integer32 := 4*(deg+1)-1;
    ryd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
        := Allocate_Work_Space(nbt,dim,degdim);
    iyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
        := Allocate_Work_Space(nbt,dim,degdim);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    u : constant Standard_Floating_VecVecs.VecVec(1..nbt)
      := Standard_Vector_Splitters.Allocate_Floating_Coefficients(nbt,3);
    v : constant Standard_Floating_VecVecs.VecVec(1..nbt)
      := Standard_Vector_Splitters.Allocate_Floating_Coefficients(nbt,3);
    w : constant Standard_Floating_VecVecs.VecVec(1..nbt)
      := Standard_Vector_Splitters.Allocate_Floating_Coefficients(nbt,3);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : QuadDobl_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rqd,iqd : quad_double;
      m : integer32;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          rxpw := rpwt(idx); ixpw := ipwt(idx);
          Multiply(xr(idx),xi(idx),xr(idx),xi(idx),rxpw(1),ixpw(1),
                   u(i),v(i),w(i));
          for k in 2..(mxe(idx)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),xr(idx),xi(idx),rxpw(k),ixpw(k),
                     u(i),v(i),w(i));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,xr.all,xi.all,rpwt,ipwt,ryd(i).all,iyd(i).all,
                 u(i),v(i),w(i));
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(xr'last+1); ivright := iydi(xi'last+1);
        m := rvright'first;
        for j in vy'range loop
          vleft := vy(j);
          rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
          iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
          vleft(idx) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
          rvright(m) := 0.0;   ivright(m) := 0.0;
          rvright(m+1) := 0.0; ivright(m+1) := 0.0;
          rvright(m+2) := 0.0; ivright(m+2) := 0.0;
          rvright(m+3) := 0.0; ivright(m+3) := 0.0;
          m := m + 4;
        end loop;
        for j in 1..xr'last loop
          rvright := rydi(j); ivright := iydi(j);
          m := rvright'first;
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
            iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
            mleft(idx,j) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
            rvright(m) := 0.0;   ivright(m) := 0.0;
            rvright(m+1) := 0.0; ivright(m+1) := 0.0;
            rvright(m+2) := 0.0; ivright(m+2) := 0.0;
            rvright(m+3) := 0.0; ivright(m+3) := 0.0;
            m := m + 4;
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : QuadDobl_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rqd,iqd : quad_double;
      m : integer32;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          rxpw := rpwt(idx); ixpw := ipwt(idx);
          Multiply(xr(idx),xi(idx),xr(idx),xi(idx),rxpw(1),ixpw(1),
                   u(i),v(i),w(i));
          for k in 2..(mxe(idx)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),xr(idx),xi(idx),rxpw(k),ixpw(k),
                     u(i),v(i),w(i));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,xr.all,xi.all,rpwt,ipwt,ryd(i).all,iyd(i).all,
                 u(i),v(i),w(i));
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(xr'last+1); ivright := iydi(xi'last+1);
        m := rvright'first;
        for j in vy'range loop
          vleft := vy(j);
          rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
          iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
          vleft(idx) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
          rvright(m) := 0.0;   ivright(m) := 0.0;
          rvright(m+1) := 0.0; ivright(m+1) := 0.0;
          rvright(m+2) := 0.0; ivright(m+2) := 0.0;
          rvright(m+3) := 0.0; ivright(m+3) := 0.0;
          m := m + 4;
        end loop;
        for j in 1..xr'last loop
          rvright := rydi(j); ivright := iydi(j);
          m := rvright'first;
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
            iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
            mleft(idx,j) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
            rvright(m) := 0.0;   ivright(m) := 0.0;
            rvright(m+1) := 0.0; ivright(m+1) := 0.0;
            rvright(m+2) := 0.0; ivright(m+2) := 0.0;
            rvright(m+3) := 0.0; ivright(m+3) := 0.0;
            m := m + 4;
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end QuadDobl_Static_EvalDiff;

  procedure QuadDobl_Dynamic_EvalDiff
              ( nbt : in integer32;
                c : in QuadDobl_Coefficient_Convolutions.Circuits;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

  -- DESCRIPTION :
  --   New procedure to define QuadDobl_Multitasked_EvalDiff
  --   with dynamic load balancing.  Except for the flag "static",
  --   the specification is the same as QuadDobl_Multitasked_EvalDiff.

    use QuadDobl_Vector_Splitters;
    use QuadDobl_Coefficient_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    degdim : constant integer32 := 4*(deg+1)-1;
    ryd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
        := Allocate_Work_Space(nbt,dim,degdim);
    iyd : constant Standard_Floating_VecVecVecs.VecVecVec(1..nbt)
        := Allocate_Work_Space(nbt,dim,degdim);
    u : constant Standard_Floating_VecVecs.VecVec(1..nbt)
      := Standard_Vector_Splitters.Allocate_Floating_Coefficients(nbt,3);
    v : constant Standard_Floating_VecVecs.VecVec(1..nbt)
      := Standard_Vector_Splitters.Allocate_Floating_Coefficients(nbt,3);
    w : constant Standard_Floating_VecVecs.VecVec(1..nbt)
      := Standard_Vector_Splitters.Allocate_Floating_Coefficients(nbt,3);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    idx1,idx2 : integer32 := 0; -- global indices
   -- idx1 : integer32 := 0; -- global indices
   -- idx2 : integer32 := c'last+1; -- only power table
    lck1,lck2 : Semaphore.Lock;

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      mdx : integer32; -- local index
      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : QuadDobl_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rqd,iqd : quad_double;
      m : integer32;
 
    begin
      while idx1 <= c'last loop
        Semaphore.Request(lck1);
        if idx1 <= c'last
         then idx1 := idx1 + 1; mdx := idx1;
        end if;
        Semaphore.Release(lck1);
        exit when (mdx > c'last);
        if mxe(mdx) > 2 then
          rxpw := rpwt(mdx); ixpw := ipwt(mdx);
          Multiply(xr(mdx),xi(mdx),xr(mdx),xi(mdx),rxpw(1),ixpw(1),
                   u(i),v(i),w(i));
          for k in 2..(mxe(mdx)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),xr(mdx),xi(mdx),rxpw(k),ixpw(k),
                     u(i),v(i),w(i));
          end loop;
        end if;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(n,pwtdone) loop
        delay 0.001;
      end loop;
      while idx2 <= c'last loop
        Semaphore.Request(lck2);
        if idx2 <= c'last
         then idx2 := idx2 + 1; mdx := idx2;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > c'last);
        EvalDiff(c(mdx).all,xr.all,xi.all,rpwt,ipwt,ryd(i).all,iyd(i).all,
                 u(i),v(i),w(i));
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(xr'last+1); ivright := iydi(xi'last+1);
        m := rvright'first;
        for j in vy'range loop
          vleft := vy(j);
          rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
          iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
          vleft(mdx) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
          rvright(m) := 0.0;   ivright(m) := 0.0;
          rvright(m+1) := 0.0; ivright(m+1) := 0.0;
          rvright(m+2) := 0.0; ivright(m+2) := 0.0;
          rvright(m+3) := 0.0; ivright(m+3) := 0.0;
          m := m + 4;
        end loop;
        for j in 1..xr'last loop
          rvright := rydi(j); ivright := iydi(j);
          m := rvright'first;
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row mdx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
            iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
            mleft(mdx,j) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
            rvright(m) := 0.0;   ivright(m) := 0.0;
            rvright(m+1) := 0.0; ivright(m+1) := 0.0;
            rvright(m+2) := 0.0; ivright(m+2) := 0.0;
            rvright(m+3) := 0.0; ivright(m+3) := 0.0;
            m := m + 4;
          end loop;
        end loop;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      mdx : integer32; -- local index
      rydi,iydi : Standard_Floating_VecVecs.Link_to_VecVec;
      vleft : QuadDobl_Complex_Vectors.Link_to_Vector;
      rvright,ivright : Standard_Floating_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
      rxpw,ixpw : Standard_Floating_VecVecs.Link_to_VecVec;
      rqd,iqd : quad_double;
      m : integer32;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx1 <= c'last loop
        Semaphore.Request(lck1);
        if idx1 > c'last then
          put_line("task " & Multitasking.to_string(i)
                           & " leaves the loop");
        else
          put_line("task " & Multitasking.to_string(i)
                           & " increments idx "
                           & Multitasking.to_string(idx1));
          idx1 := idx1 + 1; mdx := idx1;
          if idx1 <= c'last then
            put_line("task " & Multitasking.to_string(i)
                             & " will do job "
                             & Multitasking.to_string(idx1));
          else
            put_line("task " & Multitasking.to_string(i)
                             & " leaves the loop");
          end if;
        end if;
        Semaphore.Release(lck1);
        exit when (mdx > c'last);
        if mxe(mdx) > 2 then
          rxpw := rpwt(mdx); ixpw := ipwt(mdx);
          Multiply(xr(mdx),xi(mdx),xr(mdx),xi(mdx),rxpw(1),ixpw(1),
                   u(i),v(i),w(i));
          for k in 2..(mxe(mdx)-2) loop
            Multiply(rxpw(k-1),ixpw(k-1),xr(mdx),xi(mdx),rxpw(k),ixpw(k),
                     u(i),v(i),w(i));
          end loop;
        end if;
      end loop;
      put_line("task " & Multitasking.to_string(i)
                       & " is done with first loop");
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(n,pwtdone) loop
        delay 0.001;
      end loop;
      while idx2 <= c'last loop
        Semaphore.Request(lck2);
        if idx2 > c'last then
          put_line("task " & Multitasking.to_string(i)
                           & " leaves the loop");
        else
          put_line("task " & Multitasking.to_string(i)
                           & " increments idx "
                           & Multitasking.to_string(idx2));
          idx2 := idx2 + 1; mdx := idx2;
          if idx2 <= c'last then
            put_line("task " & Multitasking.to_string(i)
                             & " will do job "
                             & Multitasking.to_string(idx2));
          else
            put_line("task " & Multitasking.to_string(i)
                             & " leaves the loop");
          end if;
        end if;
        Semaphore.Release(lck2);
        exit when (mdx > c'last);
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(mdx));
        EvalDiff(c(mdx).all,xr.all,xi.all,rpwt,ipwt,ryd(i).all,iyd(i).all,
                 u(i),v(i),w(i));
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(mdx));
        rydi := ryd(i); iydi := iyd(i);
        rvright := rydi(xr'last+1); ivright := iydi(xi'last+1);
        m := rvright'first;
        for j in vy'range loop
          vleft := vy(j);
          rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
          iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
          vleft(mdx) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
          rvright(m) := 0.0;   ivright(m) := 0.0;
          rvright(m+1) := 0.0; ivright(m+1) := 0.0;
          rvright(m+2) := 0.0; ivright(m+2) := 0.0;
          rvright(m+3) := 0.0; ivright(m+3) := 0.0;
          m := m + 4;
        end loop;
        for j in 1..xr'last loop
          rvright := rydi(j); ivright := iydi(j);
          m := rvright'first;
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row mdx in vm(k) is the equation
                                       -- column j in vm(k) is the variable
            rqd := Create(rvright(m),rvright(m+1),rvright(m+2),rvright(m+3));
            iqd := Create(ivright(m),ivright(m+1),ivright(m+2),ivright(m+3));
            mleft(mdx,j) := QuadDobl_Complex_Numbers.Create(rqd,iqd);
            rvright(m) := 0.0;   ivright(m) := 0.0;
            rvright(m+1) := 0.0; ivright(m+1) := 0.0;
            rvright(m+2) := 0.0; ivright(m+2) := 0.0;
            rvright(m+3) := 0.0; ivright(m+3) := 0.0;
            m := m + 4;
          end loop;
        end loop;
      end loop;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
  end QuadDobl_Dynamic_EvalDiff;

  procedure QuadDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in QuadDobl_Coefficient_Convolutions.Circuits;
                xr,xi : in Standard_Floating_VecVecs.Link_to_VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                rpwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ipwt : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                static : in boolean := false;
                output : in boolean := false ) is
  begin
    if static
     then QuadDobl_Static_EvalDiff(nbt,c,xr,xi,mxe,rpwt,ipwt,vy,vm,output);
     else QuadDobl_Dynamic_EvalDiff(nbt,c,xr,xi,mxe,rpwt,ipwt,vy,vm,output);
    end if;
  end QuadDobl_Multitasked_EvalDiff;

  procedure Standard_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Speelpenning_Convolutions.Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use Standard_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : Standard_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : Standard_Complex_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
      xpw : Standard_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := Standard_Complex_Numbers.Create(0.0);
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : Standard_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : Standard_Complex_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
      xpw : Standard_Complex_VecVecs.Link_to_VecVec;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset index for evaluation and differentation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := Standard_Complex_Numbers.Create(0.0);
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end Standard_Multitasked_EvalDiff;

  procedure DoblDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use DoblDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : DoblDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : DoblDobl_Complex_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
      xpw : DoblDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := DoblDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := DoblDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : DoblDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : DoblDobl_Complex_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
      xpw : DoblDobl_Complex_VecVecs.Link_to_VecVec;

    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := DoblDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := DoblDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end DoblDobl_Multitasked_EvalDiff;

  procedure TripDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in TripDobl_Speelpenning_Convolutions.Circuits;
                x : in TripDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in TripDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in TripDobl_Complex_VecVecs.VecVec;
                vm : in TripDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use TripDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : TripDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : TripDobl_Complex_Vectors.Link_to_Vector;
      mleft : TripDobl_Complex_Matrices.Link_to_Matrix;
      xpw : TripDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := TripDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := TripDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : TripDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : TripDobl_Complex_Vectors.Link_to_Vector;
      mleft : TripDobl_Complex_Matrices.Link_to_Matrix;
      xpw : TripDobl_Complex_VecVecs.Link_to_VecVec;

    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := TripDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := TripDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end TripDobl_Multitasked_EvalDiff;

  procedure QuadDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use QuadDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : QuadDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : QuadDobl_Complex_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
      xpw : QuadDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := QuadDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := QuadDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : QuadDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : QuadDobl_Complex_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
      xpw : QuadDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := QuadDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := QuadDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end QuadDobl_Multitasked_EvalDiff;

  procedure PentDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in PentDobl_Speelpenning_Convolutions.Circuits;
                x : in PentDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in PentDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in PentDobl_Complex_VecVecs.VecVec;
                vm : in PentDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use PentDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : PentDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : PentDobl_Complex_Vectors.Link_to_Vector;
      mleft : PentDobl_Complex_Matrices.Link_to_Matrix;
      xpw : PentDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := PentDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := PentDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : PentDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : PentDobl_Complex_Vectors.Link_to_Vector;
      mleft : PentDobl_Complex_Matrices.Link_to_Matrix;
      xpw : PentDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := PentDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := PentDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end PentDobl_Multitasked_EvalDiff;

  procedure OctoDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in OctoDobl_Speelpenning_Convolutions.Circuits;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in OctoDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in OctoDobl_Complex_VecVecs.VecVec;
                vm : in OctoDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use OctoDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : OctoDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : OctoDobl_Complex_Vectors.Link_to_Vector;
      mleft : OctoDobl_Complex_Matrices.Link_to_Matrix;
      xpw : OctoDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := OctoDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := OctoDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : OctoDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : OctoDobl_Complex_Vectors.Link_to_Vector;
      mleft : OctoDobl_Complex_Matrices.Link_to_Matrix;
      xpw : OctoDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := OctoDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := OctoDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end OctoDobl_Multitasked_EvalDiff;

  procedure DecaDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in DecaDobl_Speelpenning_Convolutions.Circuits;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in DecaDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in DecaDobl_Complex_VecVecs.VecVec;
                vm : in DecaDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use DecaDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : DecaDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : DecaDobl_Complex_Vectors.Link_to_Vector;
      mleft : DecaDobl_Complex_Matrices.Link_to_Matrix;
      xpw : DecaDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := DecaDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := DecaDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : DecaDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : DecaDobl_Complex_Vectors.Link_to_Vector;
      mleft : DecaDobl_Complex_Matrices.Link_to_Matrix;
      xpw : DecaDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := DecaDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := DecaDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end DecaDobl_Multitasked_EvalDiff;

  procedure HexaDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in HexaDobl_Speelpenning_Convolutions.Circuits;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                mxe : in Standard_Integer_Vectors.Vector;
                pwt : in HexaDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in HexaDobl_Complex_VecVecs.VecVec;
                vm : in HexaDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use HexaDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    pwtdone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    alldone : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   without intermediate output.

      idx : integer32 := i;
      ydi : HexaDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : HexaDobl_Complex_Vectors.Link_to_Vector;
      mleft : HexaDobl_Complex_Matrices.Link_to_Matrix;
      xpw : HexaDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := HexaDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := HexaDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        idx := idx + n;
      end loop;
      alldone(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate circuits
    --   with intermediate output.

      idx : integer32 := i;
      ydi : HexaDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : HexaDobl_Complex_Vectors.Link_to_Vector;
      mleft : HexaDobl_Complex_Matrices.Link_to_Matrix;
      xpw : HexaDobl_Complex_VecVecs.Link_to_VecVec;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop  -- start with power table computation
        if mxe(idx) > 2 then
          xpw := pwt(idx);
          Multiply(x(idx),x(idx),xpw(1));
          for k in 2..(mxe(idx)-2) loop
            Multiply(xpw(k-1),x(idx),xpw(k));
          end loop;
        end if;
        idx := idx + n;
      end loop;
      pwtdone(i) := true;
     -- make sure all tasks are done with the power table
      while not Multitasking.all_true(nbt,pwtdone) loop
        delay 0.001;
      end loop;
      idx := i; -- reset the index for evaluation and differentiation
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vright'range loop
          vleft := vy(j);
          vleft(idx) := vright(j);
          vright(j) := HexaDobl_Complex_Numbers.Create(integer32(0));
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := HexaDobl_Complex_Numbers.Create(integer32(0));
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      alldone(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,alldone) loop
      delay 0.001;
    end loop;
  end HexaDobl_Multitasked_EvalDiff;

end Multitasked_AlgoDiff_Convolutions;
