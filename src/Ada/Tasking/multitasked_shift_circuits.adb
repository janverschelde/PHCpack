with text_io;                            use text_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Shift_Convolution_Circuits;         use Shift_Convolution_Circuits;
with Multitasking;

package body Multitasked_Shift_Circuits is

  procedure Standard_Multitasking
              ( nbtasks,deg : in integer32;
                c : in out Standard_Speelpenning_Convolutions.Circuits;
                t : in double_float;
                output : in boolean := false ) is

    use Standard_Complex_Numbers,Standard_Speelpenning_Convolutions;

    procedure Silent_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will shift coefficients without intermediate output.

      idx : integer32 := i;
      zero : constant Complex_Number := Create(0.0);
      wrk : constant Standard_Complex_Vectors.Link_to_Vector
          := new Standard_Complex_Vectors.Vector'(0..deg => zero);

    begin
      while idx <= c'last loop
        Shift(c(idx),wrk,t);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will shift coefficientswith intermediate output.

      idx : integer32 := i;
      zero : constant Complex_Number := Create(0.0);
      wrk : constant Standard_Complex_Vectors.Link_to_Vector
          := new Standard_Complex_Vectors.Vector'(0..deg => zero);

    begin
      while idx <= c'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " shifts component "
                         & Multitasking.to_string(idx));
        Shift(c(idx),wrk,t);
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
              ( nbtasks,deg : in integer32;
                c : in out DoblDobl_Speelpenning_Convolutions.Circuits;
                t : in double_double;
                output : in boolean := false ) is

    use DoblDobl_Complex_Numbers,DoblDobl_Speelpenning_Convolutions;

    procedure Silent_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will shift coefficients without intermediate output.

      idx : integer32 := i;
      zero : constant Complex_Number := Create(integer(0));
      wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
          := new DoblDobl_Complex_Vectors.Vector'(0..deg => zero);

    begin
      while idx <= c'last loop
        Shift(c(idx),wrk,t);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will shift coefficientswith intermediate output.

      idx : integer32 := i;
      zero : constant Complex_Number := Create(integer(0));
      wrk : constant DoblDobl_Complex_Vectors.Link_to_Vector
          := new DoblDobl_Complex_Vectors.Vector'(0..deg => zero);

    begin
      while idx <= c'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " shifts component "
                         & Multitasking.to_string(idx));
        Shift(c(idx),wrk,t);
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
              ( nbtasks,deg : in integer32;
                c : in out QuadDobl_Speelpenning_Convolutions.Circuits;
                t : in quad_double;
                output : in boolean := false ) is

    use QuadDobl_Complex_Numbers,QuadDobl_Speelpenning_Convolutions;

    procedure Silent_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will shift coefficients without intermediate output.

      idx : integer32 := i;
      zero : constant Complex_Number := Create(integer(0));
      wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
          := new QuadDobl_Complex_Vectors.Vector'(0..deg => zero);

    begin
      while idx <= c'last loop
        Shift(c(idx),wrk,t);
        idx := idx + n;
      end loop;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will shift coefficientswith intermediate output.

      idx : integer32 := i;
      zero : constant Complex_Number := Create(integer(0));
      wrk : constant QuadDobl_Complex_Vectors.Link_to_Vector
          := new QuadDobl_Complex_Vectors.Vector'(0..deg => zero);

    begin
      while idx <= c'last loop
        put_line("Task " & Multitasking.to_string(i)
                         & " shifts component "
                         & Multitasking.to_string(idx));
        Shift(c(idx),wrk,t);
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

end Multitasked_Shift_Circuits;
