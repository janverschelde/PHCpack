with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Hessian_Convolution_Circuits;       use Hessian_Convolution_Circuits;
with Multitasking;

procedure ts_mthessian is

-- DESCRIPTION :
--   Tests the development of the Hessian criterion,
--   in double, double double, and quad double arithmetic,
--   with multitasking for shared memory parallel computers.

  procedure Write_Singular_Values 
              ( values : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Writes the singular values in values, line by line.

    lnk : Standard_Complex_Vectors.Link_to_Vector;
    val : double_float;

  begin
    for k in values'range loop
      lnk := values(k);
      for i in lnk'range loop
        val := Standard_Complex_Numbers.REAL_PART(lnk(i));
        put(val,3);
      end loop;
      new_line;
    end loop;
  end Write_Singular_Values;

  procedure Multitasked_Singular_Values
              ( nbt : in integer32;
                s : in Standard_Speelpenning_Convolutions.Link_to_System;
                x : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Evaluates all Hessians of the circuits in s at x
  --   and computes the singular values with nbt tasks.

  -- ON ENTRY :
  --   nbt      the number of tasks;
  --   s        a system of convolution circuits;
  --   x        coordinates of the point to evaluate the Hessians.

    procedure Job ( i,n : integer32 ) is 

    -- DESCRIPTION :
    --   Task i out of n will evaluate the Hessian with index i + k*n,
    --   for all k starting at 0 as long as i + k*n <= s.dim.

      idx : integer32 := i; 

    begin
      while idx <= s.dim loop
        put_line("task " & Multitasking.to_string(i)
                         & " computes Hessian "
                         & Multitasking.to_string(idx));
        idx := idx + n;
      end loop;
    end Job;
    procedure do_jobs is new Multitasking.Reporting_Workers(Job);

  begin
    do_jobs(nbt);
  end Multitasked_Singular_Values;

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random Newton homotopy in double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

    s : Link_to_System;
    x : Link_to_VecVec;
    lnk : Standard_Complex_Vectors.Link_to_Vector;
    vx : Standard_Complex_Vectors.Vector(1..dim);
    A,U,V : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    e : Standard_Complex_Vectors.Vector(1..dim);
    svl : Standard_Complex_VecVecs.VecVec(1..dim);
    nbt : integer32 := 0;

  begin
    Standard_Random_Newton_Homotopy(dim,deg,nbr,pwr,s,x);
    for i in 1..dim loop
      lnk := x(i);
      vx(i) := lnk(0);
    end loop;
    for i in 1..dim loop
      declare
        vls : Standard_Complex_Vectors.Vector(1..dim);
      begin
        Singular_Values(s.crc(i),vx,A,U,V,e,vls);
        svl(i) := new Standard_Complex_Vectors.Vector'(vls);
      end;
    end loop;
    put_line("All singular values :");
    Write_Singular_Values(svl);
    new_line;
    put("Give the number of tasks : "); get(nbt);
    Multitasked_Singular_Values(nbt,s,vx);
  end Standard_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Launches the test after prompting for the parameters
  --   to generate a random problem.

    dim,deg,nbr,pwr : integer32 := 0;

  begin
    new_line;
    put_line("Testing the Hessian criterion ...");
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the power series : "); get(deg);
    put("Give the number of terms in each circuit : "); get(nbr);
    put("Give the largest power of the variables : "); get(pwr);
    Standard_Random_Test(dim,deg,nbr,pwr);
  end Main;

begin
  Main;
end ts_mthessian;
