with text_io;                            use text_io;
with Standard_Integer_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Rational_Approximations;
with DoblDobl_Rational_Approximations;
with QuadDobl_Rational_Approximations;
with Multitasking;

package body Multitasked_Pade_Approximations is

  procedure Standard_Multitasking
              ( nbtasks,numdeg,dendeg : in integer32;
                cff : in Standard_Complex_VecVecs.VecVec;
                numcff,dencff : in Standard_Complex_VecVecs.VecVec;
                t : in double_float;
                eva : out Standard_Complex_Vectors.Vector;
                output : in boolean := false ) is

    use Standard_Rational_Approximations;

    procedure Silent_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will construct a rational approximation
    --   without intermediate output.

      idx : integer32 := i;
      mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
      rhs : Standard_Complex_Vectors.Vector(1..dendeg);
      ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
      info : integer32;
      icff,inum,iden : Standard_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= cff'last loop
        icff := cff(idx); inum := numcff(idx); iden := dencff(idx);
        Pade(numdeg,dendeg,icff.all,inum.all,iden.all,mat,rhs,ipvt,info,false);
        eva(idx) := Evaluate(inum.all,iden.all,t);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will construct a rational approximation
    --   with intermediate output.

      idx : integer32 := i;
      mat : Standard_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
      rhs : Standard_Complex_Vectors.Vector(1..dendeg);
      ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
      info : integer32;
      icff,inum,iden : Standard_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= cff'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " computes component "
                         & Multitasking.to_string(idx));
        icff := cff(idx); inum := numcff(idx); iden := dencff(idx);
        Pade(numdeg,dendeg,icff.all,inum.all,iden.all,mat,rhs,ipvt,info,false);
        eva(idx) := Evaluate(inum.all,iden.all,t);
        idx := idx + n;
      end loop;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Reporting_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
  end Standard_Multitasking;

  procedure DoblDobl_Multitasking
              ( nbtasks,numdeg,dendeg : in integer32;
                cff : in DoblDobl_Complex_VecVecs.VecVec;
                numcff,dencff : in DoblDobl_Complex_VecVecs.VecVec;
                t : in double_double;
                eva : out DoblDobl_Complex_Vectors.Vector;
                output : in boolean := false ) is

    use DoblDobl_Rational_Approximations;

    procedure Silent_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will construct a rational approximation
    --   without intermediate output.

      idx : integer32 := i;
      mat : DoblDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
      rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
      ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
      info : integer32;
      icff,inum,iden : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= cff'last loop
        icff := cff(idx); inum := numcff(idx); iden := dencff(idx);
        Pade(numdeg,dendeg,icff.all,inum.all,iden.all,mat,rhs,ipvt,info,false);
        eva(idx) := Evaluate(inum.all,iden.all,t);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will construct a rational approximation
    --   with intermediate output.

      idx : integer32 := i;
      mat : DoblDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
      rhs : DoblDobl_Complex_Vectors.Vector(1..dendeg);
      ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
      info : integer32;
      icff,inum,iden : DoblDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= cff'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " computes component "
                         & Multitasking.to_string(idx));
        icff := cff(idx); inum := numcff(idx); iden := dencff(idx);
        Pade(numdeg,dendeg,icff.all,inum.all,iden.all,mat,rhs,ipvt,info,false);
        eva(idx) := Evaluate(inum.all,iden.all,t);
        idx := idx + n;
      end loop;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Reporting_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
  end DoblDobl_Multitasking;

  procedure QuadDobl_Multitasking
              ( nbtasks,numdeg,dendeg : in integer32;
                cff : in QuadDobl_Complex_VecVecs.VecVec;
                numcff,dencff : in QuadDobl_Complex_VecVecs.VecVec;
                t : in quad_double;
                eva : out QuadDobl_Complex_Vectors.Vector;
                output : in boolean := false ) is

    use QuadDobl_Rational_Approximations;

    procedure Silent_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will construct a rational approximation
    --   without intermediate output.

      idx : integer32 := i;
      mat : QuadDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
      rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
      ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
      info : integer32;
      icff,inum,iden : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= cff'last loop
        icff := cff(idx); inum := numcff(idx); iden := dencff(idx);
        Pade(numdeg,dendeg,icff.all,inum.all,iden.all,mat,rhs,ipvt,info,false);
        eva(idx) := Evaluate(inum.all,iden.all,t);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will construct a rational approximation
    --   with intermediate output.

      idx : integer32 := i;
      mat : QuadDobl_Complex_Matrices.Matrix(1..dendeg,1..dendeg);
      rhs : QuadDobl_Complex_Vectors.Vector(1..dendeg);
      ipvt : Standard_Integer_Vectors.Vector(1..dendeg);
      info : integer32;
      icff,inum,iden : QuadDobl_Complex_Vectors.Link_to_Vector;

    begin
      while idx <= cff'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " computes component "
                         & Multitasking.to_string(idx));
        icff := cff(idx); inum := numcff(idx); iden := dencff(idx);
        Pade(numdeg,dendeg,icff.all,inum.all,iden.all,mat,rhs,ipvt,info,false);
        eva(idx) := Evaluate(inum.all,iden.all,t);
        idx := idx + n;
      end loop;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Reporting_Workers(Report_Job);

  begin
    if output
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
  end QuadDobl_Multitasking;

end Multitasked_Pade_Approximations;
