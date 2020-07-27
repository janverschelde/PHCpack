with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
-- with Standard_Complex_VecVecs_io;
-- with DoblDobl_Complex_VecVecs_io;
-- with QuadDobl_Complex_VecVecs_io;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_QR_Least_Squares;  use Standard_Complex_QR_Least_Squares;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with DoblDobl_Complex_QR_Least_Squares;  use DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_QR_Least_Squares;  use QuadDobl_Complex_QR_Least_Squares;
with QuadDobl_Complex_Singular_Values;
with Standard_Series_Matrix_Solvers;
with DoblDobl_Series_Matrix_Solvers;
with QuadDobl_Series_Matrix_Solvers;
with Multitasking;

package body Multitasked_Series_Linearization is

  function Allocate_Work_Space
             ( nbt,dim : integer32 )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(1..nbt);

  begin
    for k in res'range loop
      declare
        cff : constant Standard_Complex_Vectors.Vector(1..dim)
            := (1..dim => Standard_Complex_Numbers.Create(0.0));
      begin
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..nbt);
    zero : constant double_double := create(0.0);

  begin
    for k in res'range loop
      declare
        cff : constant DoblDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => DoblDobl_Complex_Numbers.Create(zero));
      begin
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..nbt);
    zero : constant quad_double := Create(0.0);

  begin
    for k in res'range loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => QuadDobl_Complex_Numbers.Create(zero));
      begin
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  procedure MV_Multiply
             ( dim : in integer32;
               A : in Standard_Complex_Matrices.Link_to_Matrix;
               x,y : in Standard_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use Standard_Complex_Numbers;

  begin
    while Ak <= dim loop
      y(Ak) := A(Ak,1)*x(1);
      AL := 2;
      while AL <= dim loop
        y(Ak) := y(Ak) + A(Ak,AL)*x(AL);
        AL := AL + 1;
      end loop;
      Ak := Ak + 1;
    end loop;
  end MV_Multiply;

  procedure MV_Multiply
             ( dim : in integer32;
               A : in DoblDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use DoblDobl_Complex_Numbers;

  begin
    while Ak <= dim loop
      y(Ak) := A(Ak,1)*x(1);
      AL := 2;
      while AL <= dim loop
        y(Ak) := y(Ak) + A(Ak,AL)*x(AL);
        AL := AL + 1;
      end loop;
      Ak := Ak + 1;
    end loop;
  end MV_Multiply;

  procedure MV_Multiply
             ( dim : in integer32;
               A : in QuadDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use QuadDobl_Complex_Numbers;

  begin
    while Ak <= dim loop
      y(Ak) := A(Ak,1)*x(1);
      AL := 2;
      while AL <= dim loop
        y(Ak) := y(Ak) + A(Ak,AL)*x(AL);
        AL := AL + 1;
      end loop;
      Ak := Ak + 1;
    end loop;
  end MV_Multiply;

  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in Standard_Complex_Matrices.Link_to_Matrix;
               x,y : in Standard_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use Standard_Complex_Numbers;

  begin
    while Ak <= nrows loop
      y(Ak) := A(Ak,1)*x(1);
      AL := 2;
      while AL <= ncols loop
        y(Ak) := y(Ak) + A(Ak,AL)*x(AL);
        AL := AL + 1;
      end loop;
      Ak := Ak + 1;
    end loop;
  end MV_Multiply;

  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in DoblDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use DoblDobl_Complex_Numbers;

  begin
    while Ak <= nrows loop
      y(Ak) := A(Ak,1)*x(1);
      AL := 2;
      while AL <= ncols loop
        y(Ak) := y(Ak) + A(Ak,AL)*x(AL);
        AL := AL + 1;
      end loop;
      Ak := Ak + 1;
    end loop;
  end MV_Multiply;

  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in QuadDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use QuadDobl_Complex_Numbers;

  begin
    while Ak <= nrows loop
      y(Ak) := A(Ak,1)*x(1);
      AL := 2;
      while AL <= ncols loop
        y(Ak) := y(Ak) + A(Ak,AL)*x(AL);
        AL := AL + 1;
      end loop;
      Ak := Ak + 1;
    end loop;
  end MV_Multiply;

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in Standard_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use Standard_Complex_Numbers;

  begin
    while idx <= dim loop
      x(idx) := x(idx) - y(idx);
      idx := idx + 1;
    end loop;
  end V_Subtract;

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use DoblDobl_Complex_Numbers;

  begin
    while idx <= dim loop
      x(idx) := x(idx) - y(idx);
      idx := idx + 1;
    end loop;
  end V_Subtract;

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use QuadDobl_Complex_Numbers;

  begin
    while idx <= dim loop
      x(idx) := x(idx) - y(idx);
      idx := idx + 1;
    end loop;
  end V_Subtract;

  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    dim : constant integer32 := ipvt'last;
    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        MV_Multiply(dim,A(myjob-idx+1),b(idx-1),wrk(i));
        V_Subtract(dim,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(dim,A(myjob-idx+1),b(idx-1),wrk(i));
        V_Subtract(dim,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_lusolve;

  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    dim : constant integer32 := ipvt'last;
    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        MV_Multiply(dim,A(myjob-idx+1),b(idx-1),wrk(i));
        V_Subtract(dim,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(dim,A(myjob-idx+1),b(idx-1),wrk(i));
        V_Subtract(dim,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_lusolve;

  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    dim : constant integer32 := ipvt'last;
    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        MV_Multiply(dim,A(myjob-idx+1),b(idx-1),wrk(i));
        V_Subtract(dim,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(dim,A(myjob-idx+1),b(idx-1),wrk(i));
        V_Subtract(dim,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            lusolve(A(0).all,ipvt'last,ipvt,b(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_lusolve;

  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                qraux : in Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_VecVecs.VecVec;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant Standard_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant Standard_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant Standard_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant Standard_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant Standard_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant Standard_Complex_Vectors.Link_to_Vector := w5(i);

    begin
      while myjob <= b'last loop
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          lw1.all := b(idx).all;
          QRLS(lead.all,nrows,ncols,qraux,
               lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            lw1.all := b(idx).all;
            QRLS(lead.all,nrows,ncols,qraux,
	         lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant Standard_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant Standard_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant Standard_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant Standard_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant Standard_Complex_Vectors.Link_to_Vector := w5(i);

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          lw1.all := b(idx).all;
          QRLS(lead.all,nrows,ncols,qraux,
               lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            lw1.all := b(idx).all;
            QRLS(lead.all,nrows,ncols,qraux,
                 lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_QRLS;

  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : in DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_VecVecs.VecVec;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w5(i);

    begin
      while myjob <= b'last loop
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          lw1.all := b(idx).all;
          QRLS(lead.all,nrows,ncols,qraux,
	       lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            lw1.all := b(idx).all;
            QRLS(lead.all,nrows,ncols,qraux,
	         lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant DoblDobl_Complex_Vectors.Link_to_Vector := w5(i);

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          lw1.all := b(idx).all;
          QRLS(lead.all,nrows,ncols,qraux,
               lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            lw1.all := b(idx).all;
            QRLS(lead.all,nrows,ncols,qraux,
                 lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_QRLS;

  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : in QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_VecVecs.VecVec;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant QuadDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w5(i);

    begin
      while myjob <= b'last loop
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          lw1.all := b(idx).all;
          QRLS(lead.all,nrows,ncols,qraux,
	       lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            lw1.all := b(idx).all;
            QRLS(lead.all,nrows,ncols,qraux,
                 lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant QuadDobl_Complex_Vectors.Link_to_Vector := w5(i);

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          lw1.all := b(idx).all;
          QRLS(lead.all,nrows,ncols,qraux,
               lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            lw1.all := b(idx).all;
            QRLS(lead.all,nrows,ncols,qraux,
                 lw1.all,lw2.all,lw3.all,x(idx).all,lw4.all,lw5.all,110,info);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_QRLS;

  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                S : in Standard_Complex_Vectors.Vector;
                Ut,V : in Standard_Complex_Matrices.Matrix;
                wrk,utb,sub : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant Standard_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use Standard_Complex_Singular_Values;

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          -- x(idx).all := Solve(U,V,S,b(idx).all);
          Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            -- x(idx).all := Solve(U,V,S,b(idx).all);
            Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          -- x(idx).all := Solve(U,V,S,b(idx).all);
          Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            -- x(idx).all := Solve(U,V,S,b(idx).all);
            Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_SVD;

  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : in DoblDobl_Complex_Vectors.Vector;
                Ut,V : in DoblDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use DoblDobl_Complex_Singular_Values;

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          -- x(idx).all := Solve(U,V,S,b(idx).all);
          Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            -- x(idx).all := Solve(U,V,S,b(idx).all);
            Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          -- x(idx).all := Solve(U,V,S,b(idx).all);
          Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            -- x(idx).all := Solve(U,V,S,b(idx).all);
            Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_SVD;

  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                S : in QuadDobl_Complex_Vectors.Vector;
                Ut,V : in QuadDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant QuadDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use QuadDobl_Complex_Singular_Values;

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          -- x(idx).all := Solve(U,V,S,b(idx).all);
          Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            -- x(idx).all := Solve(U,V,S,b(idx).all);
            Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, with intermediate output.

      myjob : integer32 := idx+i-1;

    begin
      while myjob <= b'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " updates b(" 
                         & Multitasking.to_string(myjob) & ")");
        MV_Multiply(nrows,ncols,A(myjob-idx+1),x(idx-1),wrk(i));
        V_Subtract(nrows,b(myjob),wrk(i));
        myjob := myjob + n;
        if myjob = b'last + 1 then
          put_line("Task " & Multitasking.to_string(i)
                           & " solves for x(" 
                           & Multitasking.to_string(idx) & ")");
          -- x(idx).all := Solve(U,V,S,b(idx).all);
          Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
        elsif myjob > b'last then
          if i = 1 and (n > b'last-idx) then
            put_line("Task " & Multitasking.to_string(i)
                             & " solves for x(" 
                             & Multitasking.to_string(idx) & ")");
            -- x(idx).all := Solve(U,V,S,b(idx).all);
            Solve(Ut,V,S,b(idx).all,utb(i).all,sub(i).all,x(idx).all);
	  end if;
        end if;
      end loop;
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks finish
    while not Multitasking.all_true(nbt,done) loop
      delay 0.001;
    end loop;
  end Multitasked_Solve_Next_by_SVD;

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_lusolve(k,nbt,A,b,ipvt,wrk,output);
    end loop;
    if output then
      put_line("Norm of solution components of the multitasked solve by lu :");
     -- Standard_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant double_float
              := Standard_Complex_Vector_Norms.Max_Norm(b(k).all);
        begin
          put("||x("); put(k,1); put(")|| :"); put(nrm,3); new_line;
        end;
      end loop;
    end if;
  end Multitasked_Solve_Loop_by_lusolve;

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_lusolve(k,nbt,A,b,ipvt,wrk,output);
    end loop;
    if output then
      put_line("Norm of solution components of the multitasked solve by lu :");
     -- DoblDobl_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant double_double
              := DoblDobl_Complex_Vector_Norms.Max_Norm(b(k).all);
        begin
          put("||x("); put(k,1); put(")|| : "); put(nrm,3); new_line;
        end;
      end loop;
    end if;
  end Multitasked_Solve_Loop_by_lusolve;

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_lusolve(k,nbt,A,b,ipvt,wrk,output);
    end loop;
    if output then
      put_line("Norm of solution components of the multitasked solve by lu :");
     -- QuadDobl_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant quad_double
              := QuadDobl_Complex_Vector_Norms.Max_Norm(b(k).all);
        begin
          put("||x("); put(k,1); put(")|| : "); put(nrm,3); new_line;
        end;
      end loop;
    end if;
  end Multitasked_Solve_Loop_by_lusolve;

  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                qraux : in Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_VecVecs.VecVec;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_QRLS
        (k,nbt,A,b,x,qraux,w1,w2,w3,w4,w5,wrk,output);
    end loop;
  end Multitasked_Solve_Loop_by_QRLS;

  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : in DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_VecVecs.VecVec;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_QRLS
        (k,nbt,A,b,x,qraux,w1,w2,w3,w4,w5,wrk,output);
    end loop;
  end Multitasked_Solve_Loop_by_QRLS;

  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : in QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_VecVecs.VecVec;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_QRLS
        (k,nbt,A,b,x,qraux,w1,w2,w3,w4,w5,wrk,output);
    end loop;
  end Multitasked_Solve_Loop_by_QRLS;

  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                S : in Standard_Complex_Vectors.Vector;
                Ut,V : in Standard_Complex_Matrices.Matrix;
                wrk,utb,sub : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_SVD(k,nbt,A,b,x,S,Ut,V,wrk,utb,sub,output);
    end loop;
  end Multitasked_Solve_Loop_by_SVD;

  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : in DoblDobl_Complex_Vectors.Vector;
                Ut,V : in DoblDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_SVD(k,nbt,A,b,x,S,Ut,V,wrk,utb,sub,output);
    end loop;
  end Multitasked_Solve_Loop_by_SVD;

  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                S : in QuadDobl_Complex_Vectors.Vector;
                Ut,V : in QuadDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    for k in 1..b'last loop
      if output then
        put("calling multitasked solve next for k = ");
        put(k,1); put_line(" ...");
      end if;
      Multitasked_Solve_Next_by_SVD(k,nbt,A,b,x,S,Ut,V,wrk,utb,sub,output);
    end loop;
  end Multitasked_Solve_Loop_by_SVD;

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    Standard_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufac;

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    DoblDobl_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufac;

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    QuadDobl_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufac;

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    Standard_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if rcond + 1.0 /= 1.0
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufco;

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    one : constant double_double := create(1.0);

  begin
    DoblDobl_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if rcond + one /= one
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufco;

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    one : constant quad_double := create(1.0);

  begin
    QuadDobl_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if rcond + one /= one
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufco;

  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                qraux : out Standard_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use Standard_Series_Matrix_Solvers;

  begin
    Solve_Lead_by_QRLS(A,b,x(0),qraux,
      w1(1).all,w2(1).all,w3(1).all,w4(1).all,w5(1).all,ipvt,info);
    if info = 0 then
      Multitasked_Solve_Loop_by_QRLS
        (nbt,A,b,x,qraux,w1,w2,w3,w4,w5,wrk,output);
    end if;
  end Multitasked_Solve_by_QRLS;

  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use DoblDobl_Series_Matrix_Solvers;

  begin
    Solve_Lead_by_QRLS(A,b,x(0),qraux,
      w1(1).all,w2(1).all,w3(1).all,w4(1).all,w5(1).all,ipvt,info);
    if info = 0 then
      Multitasked_Solve_Loop_by_QRLS
        (nbt,A,b,x,qraux,w1,w2,w3,w4,w5,wrk,output);
    end if;
  end Multitasked_Solve_by_QRLS;

  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out QuadDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use QuadDobl_Series_Matrix_Solvers;

  begin
    Solve_Lead_by_QRLS(A,b,x(0),qraux,
      w1(1).all,w2(1).all,w3(1).all,w4(1).all,w5(1).all,ipvt,info);
    if info = 0 then
      Multitasked_Solve_Loop_by_QRLS
        (nbt,A,b,x,qraux,w1,w2,w3,w4,w5,wrk,output);
    end if;
  end Multitasked_Solve_by_QRLS;

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                x : in Standard_Complex_VecVecs.VecVec;
                S : out Standard_Complex_Vectors.Vector;
                U,Ut,V : out Standard_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_float;
                ewrk : in Standard_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in Standard_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use Standard_Series_Matrix_Solvers;

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if 1.0 + rcond /= 1.0 then
      Ut := Standard_Complex_Singular_Values.Conjugate_Transpose(U);
      Multitasked_Solve_Loop_by_SVD(nbt,A,b,x,S,Ut,V,wrkv,utb,sub,output);
    end if;
  end Multitasked_Solve_by_SVD;

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : out DoblDobl_Complex_Vectors.Vector;
                U,Ut,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in DoblDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use DoblDobl_Series_Matrix_Solvers;

    one : constant double_double := create(integer(1));

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if one + rcond /= one then
      Ut := DoblDobl_Complex_Singular_Values.Conjugate_Transpose(U);
      Multitasked_Solve_Loop_by_SVD(nbt,A,b,x,S,Ut,V,wrkv,utb,sub,output);
    end if;
  end Multitasked_Solve_by_SVD;

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in QuadDobl_Complex_VecMats.VecMat;
                b : in QuadDobl_Complex_VecVecs.VecVec;
                x : in QuadDobl_Complex_VecVecs.VecVec;
                S : out QuadDobl_Complex_Vectors.Vector;
                U,Ut,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                ewrk : in QuadDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in QuadDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use QuadDobl_Series_Matrix_Solvers;

    one : constant quad_double := create(integer(1));

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if one + rcond /= one then
      Ut := QuadDobl_Complex_Singular_Values.Conjugate_Transpose(U);
      Multitasked_Solve_Loop_by_SVD(nbt,A,b,x,S,Ut,V,wrkv,utb,sub,output);
    end if;
  end Multitasked_Solve_by_SVD;

end Multitasked_Series_Linearization;
