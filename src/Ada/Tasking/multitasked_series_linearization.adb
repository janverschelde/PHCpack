with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Octo_Double_Numbers_io;             use Octo_Double_Numbers_io;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with TripDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with PentDobl_Complex_Vector_Norms;
with OctoDobl_Complex_Vector_Norms;
with DecaDobl_Complex_Vector_Norms;
with HexaDobl_Complex_Vector_Norms;
-- with Standard_Complex_VecVecs_io;
-- with DoblDobl_Complex_VecVecs_io;
-- with TripDobl_Complex_VecVecs_io;
-- with QuadDobl_Complex_VecVecs_io;
-- with PentDobl_Complex_VecVecs_io;
-- with OctoDobl_Complex_VecVecs_io;
-- with DecaDobl_Complex_VecVecs_io;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_QR_Least_Squares;  use Standard_Complex_QR_Least_Squares;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with DoblDobl_Complex_QR_Least_Squares;  use DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Singular_Values;
with TripDobl_Complex_Linear_Solvers;    use TripDobl_Complex_Linear_Solvers;
with TripDobl_Complex_QR_Least_Squares;  use TripDobl_Complex_QR_Least_Squares;
with TripDobl_Complex_Singular_Values;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_QR_Least_Squares;  use QuadDobl_Complex_QR_Least_Squares;
with QuadDobl_Complex_Singular_Values;
with PentDobl_Complex_Linear_Solvers;    use PentDobl_Complex_Linear_Solvers;
with PentDobl_Complex_QR_Least_Squares;  use PentDobl_Complex_QR_Least_Squares;
with PentDobl_Complex_Singular_Values;
with OctoDobl_Complex_Linear_Solvers;    use OctoDobl_Complex_Linear_Solvers;
with OctoDobl_Complex_QR_Least_Squares;  use OctoDobl_Complex_QR_Least_Squares;
with OctoDobl_Complex_Singular_Values;
with DecaDobl_Complex_Linear_Solvers;    use DecaDobl_Complex_Linear_Solvers;
with DecaDobl_Complex_QR_Least_Squares;  use DecaDobl_Complex_QR_Least_Squares;
with DecaDobl_Complex_Singular_Values;
with HexaDobl_Complex_Linear_Solvers;    use HexaDobl_Complex_Linear_Solvers;
with HexaDobl_Complex_QR_Least_Squares;  use HexaDobl_Complex_QR_Least_Squares;
with HexaDobl_Complex_Singular_Values;
with Standard_Series_Matrix_Solvers;
with DoblDobl_Series_Matrix_Solvers;
with TripDobl_Series_Matrix_Solvers;
with QuadDobl_Series_Matrix_Solvers;
with PentDobl_Series_Matrix_Solvers;
with OctoDobl_Series_Matrix_Solvers;
with DecaDobl_Series_Matrix_Solvers;
with HexaDobl_Series_Matrix_Solvers;
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
             return TripDobl_Complex_VecVecs.VecVec is

    res : TripDobl_Complex_VecVecs.VecVec(1..nbt);
    zero : constant triple_double := create(0.0);

  begin
    for k in res'range loop
      declare
        cff : constant TripDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => TripDobl_Complex_Numbers.Create(zero));
      begin
        res(k) := new TripDobl_Complex_Vectors.Vector'(cff);
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

  function Allocate_Work_Space
             ( nbt,dim : integer32 )
             return PentDobl_Complex_VecVecs.VecVec is

    res : PentDobl_Complex_VecVecs.VecVec(1..nbt);
    zero : constant penta_double := Create(0.0);

  begin
    for k in res'range loop
      declare
        cff : constant PentDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => PentDobl_Complex_Numbers.Create(zero));
      begin
        res(k) := new PentDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim : integer32 )
             return OctoDobl_Complex_VecVecs.VecVec is

    res : OctoDobl_Complex_VecVecs.VecVec(1..nbt);
    zero : constant octo_double := Create(0.0);

  begin
    for k in res'range loop
      declare
        cff : constant OctoDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => OctoDobl_Complex_Numbers.Create(zero));
      begin
        res(k) := new OctoDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim : integer32 )
             return DecaDobl_Complex_VecVecs.VecVec is

    res : DecaDobl_Complex_VecVecs.VecVec(1..nbt);
    zero : constant deca_double := Create(0.0);

  begin
    for k in res'range loop
      declare
        cff : constant DecaDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => DecaDobl_Complex_Numbers.Create(zero));
      begin
        res(k) := new DecaDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  function Allocate_Work_Space
             ( nbt,dim : integer32 )
             return HexaDobl_Complex_VecVecs.VecVec is

    res : HexaDobl_Complex_VecVecs.VecVec(1..nbt);
    zero : constant hexa_double := Create(0.0);

  begin
    for k in res'range loop
      declare
        cff : constant HexaDobl_Complex_Vectors.Vector(1..dim)
            := (1..dim => HexaDobl_Complex_Numbers.Create(zero));
      begin
        res(k) := new HexaDobl_Complex_Vectors.Vector'(cff);
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
               A : in TripDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in TripDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use TripDobl_Complex_Numbers;

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
             ( dim : in integer32;
               A : in PentDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use PentDobl_Complex_Numbers;

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
               A : in OctoDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use OctoDobl_Complex_Numbers;

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
               A : in DecaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DecaDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use DecaDobl_Complex_Numbers;

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
               A : in HexaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in HexaDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use HexaDobl_Complex_Numbers;

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
               A : in TripDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in TripDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use TripDobl_Complex_Numbers;

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

  procedure MV_Multiply
             ( nrows,ncols : in integer32;
               A : in PentDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use PentDobl_Complex_Numbers;

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
               A : in OctoDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use OctoDobl_Complex_Numbers;

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
               A : in DecaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in DecaDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use DecaDobl_Complex_Numbers;

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
               A : in HexaDobl_Complex_Matrices.Link_to_Matrix;
               x,y : in HexaDobl_Complex_Vectors.Link_to_Vector ) is

    Ak,AL : integer32 := 1;

    use HexaDobl_Complex_Numbers;

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
                x,y : in TripDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use TripDobl_Complex_Numbers;

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

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use PentDobl_Complex_Numbers;

  begin
    while idx <= dim loop
      x(idx) := x(idx) - y(idx);
      idx := idx + 1;
    end loop;
  end V_Subtract;

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in OctoDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use OctoDobl_Complex_Numbers;

  begin
    while idx <= dim loop
      x(idx) := x(idx) - y(idx);
      idx := idx + 1;
    end loop;
  end V_Subtract;

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in DecaDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use DecaDobl_Complex_Numbers;

  begin
    while idx <= dim loop
      x(idx) := x(idx) - y(idx);
      idx := idx + 1;
    end loop;
  end V_Subtract;

  procedure V_Subtract
              ( dim : in integer32;
                x,y : in HexaDobl_Complex_Vectors.Link_to_Vector ) is

    idx : integer32 := 1;

    use HexaDobl_Complex_Numbers;

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
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

  procedure Multitasked_Solve_Next_by_lusolve
              ( idx,nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
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
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
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
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
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
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : in TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_VecVecs.VecVec;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant TripDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant TripDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant TripDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant TripDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant TripDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant TripDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
      lw1 : constant TripDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant TripDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant TripDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant TripDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant TripDobl_Complex_Vectors.Link_to_Vector := w5(i);

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

  procedure Multitasked_Solve_Next_by_QRLS
              ( idx,nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : in PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_VecVecs.VecVec;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant PentDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant PentDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant PentDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant PentDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant PentDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
      lw1 : constant PentDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant PentDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant PentDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant PentDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant PentDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                qraux : in OctoDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out OctoDobl_Complex_VecVecs.VecVec;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant OctoDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
      lw1 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant OctoDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                qraux : in DecaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DecaDobl_Complex_VecVecs.VecVec;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant DecaDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
      lw1 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant DecaDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                qraux : in HexaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out HexaDobl_Complex_VecVecs.VecVec;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant HexaDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will update a right hand side vector,
    --   or solve for component idx, without intermediate output.

      myjob : integer32 := idx+i-1;
      info : integer32;
      lw1 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
      lw1 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w1(i);
      lw2 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w2(i);
      lw3 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w3(i);
      lw4 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w4(i);
      lw5 : constant HexaDobl_Complex_Vectors.Link_to_Vector := w5(i);

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : in TripDobl_Complex_Vectors.Vector;
                Ut,V : in TripDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant TripDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use TripDobl_Complex_Singular_Values;

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

  procedure Multitasked_Solve_Next_by_SVD
              ( idx,nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : in PentDobl_Complex_Vectors.Vector;
                Ut,V : in PentDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use PentDobl_Complex_Singular_Values;

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
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                S : in OctoDobl_Complex_Vectors.Vector;
                Ut,V : in OctoDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant OctoDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use OctoDobl_Complex_Singular_Values;

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
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                S : in DecaDobl_Complex_Vectors.Vector;
                Ut,V : in DecaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant DecaDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use DecaDobl_Complex_Singular_Values;

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
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                S : in HexaDobl_Complex_Vectors.Vector;
                Ut,V : in HexaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);
    lead : constant HexaDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

    use HexaDobl_Complex_Singular_Values;

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
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
     -- TripDobl_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant triple_double
              := TripDobl_Complex_Vector_Norms.Max_Norm(b(k).all);
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

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
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
     -- PentDobl_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant penta_double
              := PentDobl_Complex_Vector_Norms.Max_Norm(b(k).all);
        begin
          put("||x("); put(k,1); put(")|| : "); put(nrm,3); new_line;
        end;
      end loop;
    end if;
  end Multitasked_Solve_Loop_by_lusolve;

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
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
     -- OctoDobl_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant octo_double
              := OctoDobl_Complex_Vector_Norms.Max_Norm(b(k).all);
        begin
          put("||x("); put(k,1); put(")|| : "); put(nrm,3); new_line;
        end;
      end loop;
    end if;
  end Multitasked_Solve_Loop_by_lusolve;

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
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
     -- DecaDobl_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant deca_double
              := DecaDobl_Complex_Vector_Norms.Max_Norm(b(k).all);
        begin
          put("||x("); put(k,1); put(")|| : "); put(nrm,3); new_line;
        end;
      end loop;
    end if;
  end Multitasked_Solve_Loop_by_lusolve;

  procedure Multitasked_Solve_Loop_by_lusolve
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
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
     -- HexaDobl_Complex_VecVecs_io.put_line(b);
      for k in b'range loop
        declare
          nrm : constant hexa_double
              := HexaDobl_Complex_Vector_Norms.Max_Norm(b(k).all);
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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : in TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_VecVecs.VecVec;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
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

  procedure Multitasked_Solve_Loop_by_QRLS
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : in PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_VecVecs.VecVec;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
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
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                qraux : in OctoDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out OctoDobl_Complex_VecVecs.VecVec;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
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
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                qraux : in DecaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DecaDobl_Complex_VecVecs.VecVec;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
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
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                qraux : in HexaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out HexaDobl_Complex_VecVecs.VecVec;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : in TripDobl_Complex_Vectors.Vector;
                Ut,V : in TripDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in TripDobl_Complex_VecVecs.VecVec;
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

  procedure Multitasked_Solve_Loop_by_SVD
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : in PentDobl_Complex_Vectors.Vector;
                Ut,V : in PentDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in PentDobl_Complex_VecVecs.VecVec;
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
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                S : in OctoDobl_Complex_Vectors.Vector;
                Ut,V : in OctoDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in OctoDobl_Complex_VecVecs.VecVec;
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
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                S : in DecaDobl_Complex_Vectors.Vector;
                Ut,V : in DecaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in DecaDobl_Complex_VecVecs.VecVec;
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
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                S : in HexaDobl_Complex_Vectors.Vector;
                Ut,V : in HexaDobl_Complex_Matrices.Matrix;
                wrk,utb,sub : in HexaDobl_Complex_VecVecs.VecVec;
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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    TripDobl_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
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

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    PentDobl_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufac;

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    OctoDobl_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufac;

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    DecaDobl_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufac;

  procedure Multitasked_Solve_by_lufac
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is
  begin
    HexaDobl_Series_Matrix_Solvers.Solve_Lead_by_lufac(A,b,ipvt,info);
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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out triple_double;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    one : constant triple_double := create(1.0);

  begin
    TripDobl_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
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

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out penta_double;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    one : constant penta_double := create(1.0);

  begin
    PentDobl_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if rcond + one /= one
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufco;

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out octo_double;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    one : constant octo_double := create(1.0);

  begin
    OctoDobl_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if rcond + one /= one
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufco;

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out deca_double;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    one : constant deca_double := create(1.0);

  begin
    DecaDobl_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if rcond + one /= one
     then Multitasked_Solve_Loop_by_lusolve(nbt,A,b,ipvt,wrk,output);
    end if;
  end Multitasked_Solve_by_lufco;

  procedure Multitasked_Solve_by_lufco
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out hexa_double;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    one : constant hexa_double := create(1.0);

  begin
    HexaDobl_Series_Matrix_Solvers.Solve_Lead_by_lufco(A,b,ipvt,rcond);
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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                qraux : out TripDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out TripDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use TripDobl_Series_Matrix_Solvers;

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

  procedure Multitasked_Solve_by_QRLS
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : out PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use PentDobl_Series_Matrix_Solvers;

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
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                qraux : out OctoDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out OctoDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use OctoDobl_Series_Matrix_Solvers;

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
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                qraux : out DecaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DecaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use DecaDobl_Series_Matrix_Solvers;

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
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                qraux : out HexaDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out HexaDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use HexaDobl_Series_Matrix_Solvers;

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
                A : in TripDobl_Complex_VecMats.VecMat;
                b : in TripDobl_Complex_VecVecs.VecVec;
                x : in TripDobl_Complex_VecVecs.VecVec;
                S : out TripDobl_Complex_Vectors.Vector;
                U,Ut,V : out TripDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out triple_double;
                ewrk : in TripDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in TripDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use TripDobl_Series_Matrix_Solvers;

    one : constant triple_double := create(integer(1));

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if one + rcond /= one then
      Ut := TripDobl_Complex_Singular_Values.Conjugate_Transpose(U);
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

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : out PentDobl_Complex_Vectors.Vector;
                U,Ut,V : out PentDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out penta_double;
                ewrk : in PentDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in PentDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use PentDobl_Series_Matrix_Solvers;

    one : constant penta_double := create(integer(1));

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if one + rcond /= one then
      Ut := PentDobl_Complex_Singular_Values.Conjugate_Transpose(U);
      Multitasked_Solve_Loop_by_SVD(nbt,A,b,x,S,Ut,V,wrkv,utb,sub,output);
    end if;
  end Multitasked_Solve_by_SVD;

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in OctoDobl_Complex_VecMats.VecMat;
                b : in OctoDobl_Complex_VecVecs.VecVec;
                x : in OctoDobl_Complex_VecVecs.VecVec;
                S : out OctoDobl_Complex_Vectors.Vector;
                U,Ut,V : out OctoDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out octo_double;
                ewrk : in OctoDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in OctoDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use OctoDobl_Series_Matrix_Solvers;

    one : constant octo_double := create(integer(1));

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if one + rcond /= one then
      Ut := OctoDobl_Complex_Singular_Values.Conjugate_Transpose(U);
      Multitasked_Solve_Loop_by_SVD(nbt,A,b,x,S,Ut,V,wrkv,utb,sub,output);
    end if;
  end Multitasked_Solve_by_SVD;

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in DecaDobl_Complex_VecMats.VecMat;
                b : in DecaDobl_Complex_VecVecs.VecVec;
                x : in DecaDobl_Complex_VecVecs.VecVec;
                S : out DecaDobl_Complex_Vectors.Vector;
                U,Ut,V : out DecaDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out deca_double;
                ewrk : in DecaDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in DecaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use DecaDobl_Series_Matrix_Solvers;

    one : constant deca_double := create(integer(1));

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if one + rcond /= one then
      Ut := DecaDobl_Complex_Singular_Values.Conjugate_Transpose(U);
      Multitasked_Solve_Loop_by_SVD(nbt,A,b,x,S,Ut,V,wrkv,utb,sub,output);
    end if;
  end Multitasked_Solve_by_SVD;

  procedure Multitasked_Solve_by_SVD
              ( nbt : in integer32;
                A : in HexaDobl_Complex_VecMats.VecMat;
                b : in HexaDobl_Complex_VecVecs.VecVec;
                x : in HexaDobl_Complex_VecVecs.VecVec;
                S : out HexaDobl_Complex_Vectors.Vector;
                U,Ut,V : out HexaDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out hexa_double;
                ewrk : in HexaDobl_Complex_Vectors.Link_to_Vector;
                wrkv,utb,sub : in HexaDobl_Complex_VecVecs.VecVec;
                output : in boolean := true ) is

    use HexaDobl_Series_Matrix_Solvers;

    one : constant hexa_double := create(integer(1));

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv(1));
    if one + rcond /= one then
      Ut := HexaDobl_Complex_Singular_Values.Conjugate_Transpose(U);
      Multitasked_Solve_Loop_by_SVD(nbt,A,b,x,S,Ut,V,wrkv,utb,sub,output);
    end if;
  end Multitasked_Solve_by_SVD;

end Multitasked_Series_Linearization;
