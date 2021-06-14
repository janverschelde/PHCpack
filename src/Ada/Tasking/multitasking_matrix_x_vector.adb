with text_io;                            use text_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multitasking;                       use Multitasking;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;

package body Multitasking_Matrix_x_Vector is

  function Silent_Multiply
              ( n : integer32;
                A : Standard_Complex_Matrices.Matrix;
                v : Standard_Complex_Vectors.Vector ) 
              return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(A'range(1));

    procedure Job ( i,n : integer32 ) is

      use Standard_Complex_Numbers;

    begin
      for r in A'range(1) loop
        if r mod n = i-1 then
          res(r) := Create(0.0);
          for j in A'range(2) loop
            res(r) := res(r) + A(r,j)*v(j);
          end loop;
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Silent_Multiply;

  function Silent_Multiply
              ( n : integer32;
                A : DoblDobl_Complex_Matrices.Matrix;
                v : DoblDobl_Complex_Vectors.Vector ) 
              return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(A'range(1));

    procedure Job ( i,n : integer32 ) is

      use DoblDobl_Complex_Numbers;

    begin
      for r in A'range(1) loop
        if r mod n = i-1 then
          res(r) := Create(integer32(0));
          for j in A'range(2) loop
            res(r) := res(r) + A(r,j)*v(j);
          end loop;
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Silent_Multiply;

  function Silent_Multiply
              ( n : integer32;
                A : QuadDobl_Complex_Matrices.Matrix;
                v : QuadDobl_Complex_Vectors.Vector ) 
              return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(A'range(1));

    procedure Job ( i,n : integer32 ) is

      use QuadDobl_Complex_Numbers;

    begin
      for r in A'range(1) loop
        if r mod n = i-1 then
          res(r) := Create(integer32(0));
          for j in A'range(2) loop
            res(r) := res(r) + A(r,j)*v(j);
          end loop;
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Silent_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Silent_Multiply;

  function Reporting_Multiply
              ( n : integer32;
                A : Standard_Complex_Matrices.Matrix;
                v : Standard_Complex_Vectors.Vector ) 
              return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(A'range(1));

    procedure Job ( i,n : integer32 ) is

      use Standard_Complex_Numbers;

    begin
      put_line("hello from task " & to_string(i));
      for r in A'range(1) loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := Create(0.0);
          for j in A'range(2) loop
            res(r) := res(r) + A(r,j)*v(j);
          end loop;
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Reporting_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Reporting_Multiply;

  function Reporting_Multiply
              ( n : integer32;
                A : DoblDobl_Complex_Matrices.Matrix;
                v : DoblDobl_Complex_Vectors.Vector ) 
              return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(A'range(1));

    procedure Job ( i,n : integer32 ) is

      use DoblDobl_Complex_Numbers;

    begin
      put_line("hello from task " & to_string(i));
      for r in A'range(1) loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := Create(integer32(0));
          for j in A'range(2) loop
            res(r) := res(r) + A(r,j)*v(j);
          end loop;
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Reporting_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Reporting_Multiply;

  function Reporting_Multiply
              ( n : integer32;
                A : QuadDobl_Complex_Matrices.Matrix;
                v : QuadDobl_Complex_Vectors.Vector ) 
              return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(A'range(1));

    procedure Job ( i,n : integer32 ) is

      use QuadDobl_Complex_Numbers;

    begin
      put_line("hello from task " & to_string(i));
      for r in A'range(1) loop
        if r mod n = i-1 then
          put_line("task " & to_string(i) & " computes " & to_string(r));
          res(r) := Create(integer32(0));
          for j in A'range(2) loop
            res(r) := res(r) + A(r,j)*v(j);
          end loop;
        end if;
      end loop;
    end Job;
    procedure do_jobs is new Reporting_Workers(Job);

  begin
    do_jobs(n);
    return res;
  end Reporting_Multiply;

end Multitasking_Matrix_x_Vector;
