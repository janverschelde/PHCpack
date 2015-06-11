with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Standard_Random_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;

procedure ts_perflu is

-- DESCRIPTION :
--   Test on the performance of LU factorization.
--   Run as "perf stat -e r538010 ts_perflu < /tmp/input"
--   where /tmp/input contains the values on input.
--   There is a clear benefit for doubles, e.g.: n = 20, f = 500000,
--   but not for complex numbers and the reason for this is that the
--   auto vectorization in gcc does currently not provide any support
--   for complex numbers.

  procedure Swap ( x,y : in out double_float ) is

  -- DESCRIPTION :
  --   Swaps x with y.

    z : constant double_float := y;

  begin
    y := x;
    x := z;
  end Swap;

  function cabs ( c : Complex_Number ) return double_float is
  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

  procedure vvlufac ( a : in out Standard_Floating_VecVecs.VecVec;
                      n : in integer32;
                      ipvt : out Standard_Integer_Vectors.Vector;
                      info : out integer32 ) is

  -- DESCRIPTION :
  --   LU factorization on vector of vectors data type.

    kp1,ell,nm1 : integer32;
    smax,acc : double_float;
    ak,aj : Standard_Floating_Vectors.Link_to_Vector;

  begin
    info := 0;
    nm1 := integer32(n - 1);
    if nm1 >= 1 then
      for k in 1..nm1 loop
        ak := a(k);
        kp1 := k + 1;
        ell := k;                               -- find the pivot index ell
        smax := abs(ak(k));
        for i in kp1..n loop
          acc := abs(ak(i));
          if acc > smax 
           then ell := i; smax := acc;
          end if;
        end loop;
        ipvt(k) := ell;
        if smax = 0.0 then
          info := k;               -- this column is already triangulated
        else
          if ell /= k                         -- interchange if necessary
           then Swap(ak(ell),ak(k));
          end if;                                  -- compute multipliers
          acc := -1.0/ak(k); 
          for i in kp1..n loop
            ak(i) := ak(i)*acc;
          end loop;
          for j in kp1..n loop                         -- row elimination
            aj := a(j);
            if ell /= k
             then Swap(aj(ell),aj(k));
            end if;
            for i in kp1..n loop
              aj(i) := aj(i) + aj(k)*ak(i);
            end loop;
          end loop;
        end if;
      end loop;
    end if;
    ipvt(n) := n;
    if a(n)(n) = 0.0
     then info := n;
    end if;
  end vvlufac;

  procedure Run_Floating_LU_Factorizations ( n,f : in integer32 ) is

    mat : Standard_Floating_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    for i in 1..f loop
      mat := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
      lufac(mat,n,ipvt,info);
    end loop;
  end Run_Floating_LU_Factorizations;

  procedure Run_Complex_LU_Factorizations ( n,f : in integer32 ) is

    mat : Standard_Complex_Matrices.Matrix(1..n,1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    for i in 1..f loop
      mat := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
      lufac(mat,n,ipvt,info);
    end loop;
  end Run_Complex_LU_Factorizations;

  procedure Run_Floating_vvLU_Factorizations ( n,f : in integer32 ) is

    mat : Standard_Floating_VecVecs.VecVec(1..n);
    matj : Standard_Floating_Vectors.Link_to_Vector;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    for i in 1..n loop
      mat(i) := new Standard_Floating_Vectors.Vector(1..n);
    end loop;
    for i in 1..f loop
      for j in 1..n loop
         matj := mat(j);
         for k in 1..n loop
           matj(k) := Standard_Random_Numbers.Random;
         end loop;
      end loop;
      vvlufac(mat,n,ipvt,info);
    end loop;
  end Run_Floating_vvLU_Factorizations;

  procedure Run_Complex_vvLU_Factorizations ( n,f : in integer32 ) is

    mat : Standard_Complex_VecVecs.VecVec(1..n);
    matj : Standard_Complex_Vectors.Link_to_Vector;
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;

  begin
    for i in 1..n loop
      mat(i) := new Standard_Complex_Vectors.Vector(1..n);
    end loop;
    for i in 1..f loop
      for j in 1..n loop
         matj := mat(j);
         for k in 1..n loop
           matj(k) := Standard_Random_Numbers.Random1;
         end loop;
      end loop;
      lufac(mat,n,ipvt,info);
    end loop;
  end Run_Complex_vvLU_Factorizations;

  procedure Main is

    n,f : integer32 := 0;
    timer : Timing_Widget;

  begin
    put("Give the dimension : "); get(n);
    put("Give the frequency : "); get(f);
    tstart(timer);
    Run_Floating_LU_Factorizations(n,f);
    tstop(timer);
    print_times(standard_output,timer,"float matrix lu fac");
    tstart(timer);
    Run_Complex_LU_Factorizations(n,f);
    tstop(timer);
    print_times(standard_output,timer,"complex matrix lu fac");
    tstart(timer);
    Run_Floating_vvLU_Factorizations(n,f);
    tstop(timer);
    print_times(standard_output,timer,"float vecvec lu fac");
    tstart(timer);
    Run_Complex_vvLU_Factorizations(n,f);
    tstop(timer);
    print_times(standard_output,timer,"complex vecvec lu fac");
  end Main;

begin
  Main;
end ts_perflu;
