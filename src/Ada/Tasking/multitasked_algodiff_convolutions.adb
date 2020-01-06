with text_io;                            use text_io;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with Multitasking;

package body Multitasked_AlgoDiff_Convolutions is

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

  procedure Standard_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Speelpenning_Convolutions.Convolution_Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                pwt : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use Standard_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   without intermediate output.

      idx : integer32 := i;
      ydi : Standard_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : Standard_Complex_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
 
    begin
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
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   with intermediate output.

      idx : integer32 := i;
      ydi : Standard_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : Standard_Complex_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
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
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,done) loop
      delay 0.1;
    end loop;
  end Standard_Multitasked_EvalDiff;

  procedure DoblDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in DoblDobl_Speelpenning_Convolutions.Convolution_Circuits;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                pwt : in DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in DoblDobl_Complex_VecVecs.VecVec;
                vm : in DoblDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use DoblDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   without intermediate output.

      idx : integer32 := i;
      ydi : DoblDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : DoblDobl_Complex_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
 
    begin
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
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   with intermediate output.

      idx : integer32 := i;
      ydi : DoblDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : DoblDobl_Complex_Vectors.Link_to_Vector;
      mleft : DoblDobl_Complex_Matrices.Link_to_Matrix;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
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
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,done) loop
      delay 0.1;
    end loop;
  end DoblDobl_Multitasked_EvalDiff;

  procedure QuadDobl_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in QuadDobl_Speelpenning_Convolutions.Convolution_Circuits;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                pwt : in QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in QuadDobl_Complex_VecVecs.VecVec;
                vm : in QuadDobl_Complex_VecMats.VecMat;
                output : in boolean := false ) is

    use QuadDobl_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   without intermediate output.

      idx : integer32 := i;
      ydi : QuadDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : QuadDobl_Complex_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
 
    begin
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
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   with intermediate output.

      idx : integer32 := i;
      ydi : QuadDobl_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : QuadDobl_Complex_Vectors.Link_to_Vector;
      mleft : QuadDobl_Complex_Matrices.Link_to_Matrix;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
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
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,done) loop
      delay 0.1;
    end loop;
  end QuadDobl_Multitasked_EvalDiff;

end Multitasked_AlgoDiff_Convolutions;
