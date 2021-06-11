with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_Durand_Kerner;             use Standard_Durand_Kerner;
with DoblDobl_Durand_Kerner;             use DoblDobl_Durand_Kerner;
with QuadDobl_Durand_Kerner;             use QuadDobl_Durand_Kerner;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Rectangular_Sample_Grids;           use Rectangular_Sample_Grids;
with Make_Sample_Grids;                  use Make_Sample_Grids;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Random_Vectors;
with Multprec_Durand_Kerner;             use Multprec_Durand_Kerner;
with Sampling_Machine;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Standard_Power_Traces;              use Standard_Power_Traces;
with DoblDobl_Power_Traces;              use DoblDobl_Power_Traces;
with QuadDobl_Power_Traces;              use QuadDobl_Power_Traces;
with Multprec_Power_Traces;              use Multprec_Power_Traces;
with Standard_Univariate_Interpolators;
with DoblDobl_Univariate_Interpolators;
with QuadDobl_Univariate_Interpolators;
with Standard_Trace_Interpolators;
with Multprec_Trace_Interpolators;

procedure ts_powtrc is

-- DESCRIPTION :
--   This procedure tests the conversion of traces into sums of powers
--   of the roots of a polynomial, and vice versa.
--   Moreover, with the Newton identities, we can interpolate with an
--   optimal number of sample points, i.e.: on a triangular grid.

  function Standard_Coefficients_to_Traces 
             ( c : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the coefficients c into traces t, t(i) = c(d-i),
  --   where i runs from 1 till d, d = c'last.

  -- REQUIRED :
  --   c'range = 0..d, c(i) is coefficient with x^i, c(d) is not used.

    use Standard_Complex_Numbers,Standard_Complex_Vectors;

    d : constant integer32 := c'last;
    t : Vector(1..d);

  begin
    for i in 1..d loop
      t(i) := c(d-i);
    end loop;
    return t;
  end Standard_Coefficients_to_Traces;

  function DoblDobl_Coefficients_to_Traces 
             ( c : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the coefficients c into traces t, t(i) = c(d-i),
  --   where i runs from 1 till d, d = c'last.

  -- REQUIRED :
  --   c'range = 0..d, c(i) is coefficient with x^i, c(d) is not used.

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Vectors;

    d : constant integer32 := c'last;
    t : Vector(1..d);

  begin
    for i in 1..d loop
      t(i) := c(d-i);
    end loop;
    return t;
  end DoblDobl_Coefficients_to_Traces;

  function QuadDobl_Coefficients_to_Traces 
             ( c : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the coefficients c into traces t, t(i) = c(d-i),
  --   where i runs from 1 till d, d = c'last.

  -- REQUIRED :
  --   c'range = 0..d, c(i) is coefficient with x^i, c(d) is not used.

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Vectors;

    d : constant integer32 := c'last;
    t : Vector(1..d);

  begin
    for i in 1..d loop
      t(i) := c(d-i);
    end loop;
    return t;
  end QuadDobl_Coefficients_to_Traces;

  function Multprec_Coefficients_to_Traces 
             ( c : Multprec_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the coefficients c into traces t, t(i) = c(d-i),
  --   where i runs from 1 till d, d = c'last.

  -- REQUIRED :
  --   c'range = 0..d, c(i) is coefficient with x^i, c(d) is not used.

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;

    d : constant integer32 := c'last;
    t : Vector(1..d);

  begin
    for i in 1..d loop
      Copy(c(d-i),t(i));
    end loop;
    return t;
  end Multprec_Coefficients_to_Traces;

  function Standard_Roots ( file : file_type;
                            c : Standard_Complex_Vectors.Vector )
                          return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.

    use Standard_Complex_Numbers;

    d : constant integer32 := c'last;
    z : Standard_Complex_Vectors.Vector(1..d)
      := Standard_Random_Vectors.Random_Vector(1,d);
    maxit : constant natural32 := 10 + natural32(d);
    eps : constant double_float := 1.0E-14;
    res : Standard_Complex_Vectors.Vector(1..d);
    maxres : double_float;
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,z,res,maxit,eps,nb,fail);
   -- put_line(file,"The roots of the polynomial : "); put_line(file,z);
   -- put_line(file,"The residuals at the roots : "); put_line(file,res);
    maxres := AbsVal(res(1));
    for i in 2..d loop
      if AbsVal(res(i)) > maxres
       then maxres := AbsVal(res(i));
      end if;
    end loop;
    put(file,"Maximal residual at the roots : ");
    put(file,maxres,3); new_line(file);
    return z;
  end Standard_Roots;

  function DoblDobl_Roots ( file : file_type;
                            c : DoblDobl_Complex_Vectors.Vector )
                          return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.

    use DoblDobl_Complex_Numbers;

    d : constant integer32 := c'last;
    z : DoblDobl_Complex_Vectors.Vector(1..d)
      := DoblDobl_Random_Vectors.Random_Vector(1,d);
    maxit : constant natural32 := 20 + natural32(d);
    eps : constant double_float := 1.0E-26;
    res : DoblDobl_Complex_Vectors.Vector(1..d);
    maxres : double_double;
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,z,res,maxit,eps,nb,fail);
   -- put_line(file,"The roots of the polynomial : "); put_line(file,z);
   -- put_line(file,"The residuals at the roots : "); put_line(file,res);
    maxres := AbsVal(res(1));
    for i in 2..d loop
      if AbsVal(res(i)) > maxres
       then maxres := AbsVal(res(i));
      end if;
    end loop;
    put(file,"Maximal residual at the roots : ");
    put(file,maxres,3); new_line(file);
    return z;
  end DoblDobl_Roots;

  function QuadDobl_Roots ( file : file_type;
                            c : QuadDobl_Complex_Vectors.Vector )
                          return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.

    use QuadDobl_Complex_Numbers;

    d : constant integer32 := c'last;
    z : QuadDobl_Complex_Vectors.Vector(1..d)
      := QuadDobl_Random_Vectors.Random_Vector(1,d);
    maxit : constant natural32 := 30 + natural32(d);
    eps : constant double_float := 1.0E-58;
    res : QuadDobl_Complex_Vectors.Vector(1..d);
    maxres : quad_double;
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,z,res,maxit,eps,nb,fail);
   -- put_line(file,"The roots of the polynomial : "); put_line(file,z);
   -- put_line(file,"The residuals at the roots : "); put_line(file,res);
    maxres := AbsVal(res(1));
    for i in 2..d loop
      if AbsVal(res(i)) > maxres
       then maxres := AbsVal(res(i));
      end if;
    end loop;
    put(file,"Maximal residual at the roots : ");
    put(file,maxres,3); new_line(file);
    return z;
  end QuadDobl_Roots;

  function Multprec_Roots ( file : file_type;
                            c : Multprec_Complex_Vectors.Vector;
                            size : natural32 )
                          return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Applies Durand-Kerner to compute the roots of the polynomial
  --   defined by the coefficients in c.

    use Multprec_Complex_Numbers;

    d : constant integer32 := c'last;
    deci : constant natural32 := Size_to_Decimal(size);
    z : Multprec_Complex_Vectors.Vector(1..d)
      := Multprec_Random_Vectors.Random_Vector(1,d,size) ;
    maxit : constant natural32 := 10 + natural32(d);
    eps : constant double_float := 10.0**integer(-deci);
    mpeps : Floating_Number := Create(eps);
    res : Multprec_Complex_Vectors.Vector(1..d);
    absres,maxres : Floating_Number;
    nb : natural32;
    fail : boolean;

  begin
    Silent_Durand_Kerner(c,z,res,maxit,mpeps,nb,fail);
   -- put_line(file,"The roots of the polynomial : "); put_line(file,z);
   -- put_line(file,"The residuals at the roots : "); put_line(file,res);
    maxres := AbsVal(res(1));
    for i in 2..d loop
      absres := AbsVal(res(i));
      if absres > maxres then
        Copy(absres,maxres);
        Clear(absres);
      end if;
    end loop;
    put(file,"Maximal residual at the roots : ");
    put(file,maxres,3); new_line(file);
    Clear(maxres); Clear(mpeps);
    Multprec_Complex_Vectors.Clear(res);
    return z;
  end Multprec_Roots;

  function Standard_Power_Sums ( x : Standard_Complex_Vectors.Vector )
                               return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Given all d roots in x of a degree d polynomial,
  --   then all power sums up to degree d of the roots are returned.

    use Standard_Complex_Numbers,Standard_Complex_Vectors;

    sums : Vector(x'range);
    powers : Vector(x'range) := x;   -- contains powers of roots

  begin
    sums(1) := Sum(powers);
    for i in 2..x'last loop
      for j in x'range loop
        powers(j) := powers(j)*x(j);
      end loop;
      sums(i) := Sum(powers);
    end loop;
    return sums;
  end Standard_Power_Sums;

  function DoblDobl_Power_Sums ( x : DoblDobl_Complex_Vectors.Vector )
                               return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Given all d roots in x of a degree d polynomial,
  --   then all power sums up to degree d of the roots are returned.

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Vectors;

    sums : Vector(x'range);
    powers : Vector(x'range) := x;   -- contains powers of roots

  begin
    sums(1) := Sum(powers);
    for i in 2..x'last loop
      for j in x'range loop
        powers(j) := powers(j)*x(j);
      end loop;
      sums(i) := Sum(powers);
    end loop;
    return sums;
  end DoblDobl_Power_Sums;

  function QuadDobl_Power_Sums ( x : QuadDobl_Complex_Vectors.Vector )
                               return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Given all d roots in x of a degree d polynomial,
  --   then all power sums up to degree d of the roots are returned.

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Vectors;

    sums : Vector(x'range);
    powers : Vector(x'range) := x;   -- contains powers of roots

  begin
    sums(1) := Sum(powers);
    for i in 2..x'last loop
      for j in x'range loop
        powers(j) := powers(j)*x(j);
      end loop;
      sums(i) := Sum(powers);
    end loop;
    return sums;
  end QuadDobl_Power_Sums;

  function Multprec_Power_Sums ( x : Multprec_Complex_Vectors.Vector )
                               return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Given all d roots in x of a degree d polynomial,
  --   then all power sums up to degree d of the roots are returned.

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;

    sums,powers : Vector(x'range);

  begin
    Copy(x,powers);
    sums(1) := Sum(powers);
    for i in 2..x'last loop
      for j in x'range loop
        Mul(powers(j),x(j));
      end loop;
      sums(i) := Sum(powers);
    end loop;
    Clear(powers);
    return sums;
  end Multprec_Power_Sums;

  procedure Standard_Test_Power_Traces ( d : in integer32 ) is

  -- DESCRIPTION :
  --   The test consists of the following steps :
  --    1) generate random coefficient of a monic polynomial of degree d;
  --    2) compute power sums from the traces;
  --    3) convert power sums into traces => original coefficients?

    use Standard_Complex_Numbers,Standard_Complex_Vectors;

    c : Vector(0..d);
    s,t,st,x,ps,diff : Vector(1..d);
    normdiff : double_float;
    ans : character;
  
  begin
    c(0..d-1) := Standard_Random_Vectors.Random_Vector(0,d-1);
    c(d) := Create(1.0);
    put_line("The random coefficients of a monic polynomial ");
    put_line(c);
    t := Standard_Coefficients_to_Traces(c);
    s := Traces_to_Power_Sums(t);
    st := Power_Sums_to_Traces(s);
    put_line("The traces computed from the coefficients : ");
    put_line(t);
    put_line("The traces computed from the power sums : ");
    put_line(st);
    diff := t - st;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
    put("Do you wish to continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    x := Standard_Roots(Standard_Output,c);
    ps := Standard_Power_Sums(x);
    put_line("Power sums computed from the traces : ");
    put_line(s);
    put_line("Power sums computed from the roots : ");
    put_line(ps);
    diff := s - ps;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
  end Standard_Test_Power_Traces;

  procedure DoblDobl_Test_Power_Traces ( d : in integer32 ) is

  -- DESCRIPTION :
  --   The test consists of the following steps :
  --    1) generate random coefficient of a monic polynomial of degree d;
  --    2) compute power sums from the traces;
  --    3) convert power sums into traces => original coefficients?

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Vectors;

    c : Vector(0..d);
    s,t,st,x,ps,diff : Vector(1..d);
    one : constant double_double := create(1.0);
    normdiff : double_double;
    ans : character;
  
  begin
    c(0..d-1) := DoblDobl_Random_Vectors.Random_Vector(0,d-1);
    c(d) := Create(one);
    put_line("The random coefficients of a monic polynomial ");
    put_line(c);
    t := DoblDobl_Coefficients_to_Traces(c);
    s := Traces_to_Power_Sums(t);
    st := Power_Sums_to_Traces(s);
    put_line("The traces computed from the coefficients : ");
    put_line(t);
    put_line("The traces computed from the power sums : ");
    put_line(st);
    diff := t - st;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
    put("Do you wish to continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    x := DoblDobl_Roots(Standard_Output,c);
    ps := DoblDobl_Power_Sums(x);
    put_line("Power sums computed from the traces : ");
    put_line(s);
    put_line("Power sums computed from the roots : ");
    put_line(ps);
    diff := s - ps;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
  end DoblDobl_Test_Power_Traces;

  procedure QuadDobl_Test_Power_Traces ( d : in integer32 ) is

  -- DESCRIPTION :
  --   The test consists of the following steps :
  --    1) generate random coefficient of a monic polynomial of degree d;
  --    2) compute power sums from the traces;
  --    3) convert power sums into traces => original coefficients?

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Vectors;

    c : Vector(0..d);
    s,t,st,x,ps,diff : Vector(1..d);
    one : constant quad_double := create(1.0);
    normdiff : quad_double;
    ans : character;
  
  begin
    c(0..d-1) := QuadDobl_Random_Vectors.Random_Vector(0,d-1);
    c(d) := Create(one);
    put_line("The random coefficients of a monic polynomial ");
    put_line(c);
    t := QuadDobl_Coefficients_to_Traces(c);
    s := Traces_to_Power_Sums(t);
    st := Power_Sums_to_Traces(s);
    put_line("The traces computed from the coefficients : ");
    put_line(t);
    put_line("The traces computed from the power sums : ");
    put_line(st);
    diff := t - st;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
    put("Do you wish to continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    x := QuadDobl_Roots(Standard_Output,c);
    ps := QuadDobl_Power_Sums(x);
    put_line("Power sums computed from the traces : ");
    put_line(s);
    put_line("Power sums computed from the roots : ");
    put_line(ps);
    diff := s - ps;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
  end QuadDobl_Test_Power_Traces;

  procedure Multprec_Test_Power_Traces
              ( d : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   The test consists of the following steps :
  --    1) generate random coefficient of a monic polynomial of degree d;
  --    2) compute power sums from the traces;
  --    3) convert power sums into traces => original coefficients?

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;

    c : Vector(0..d);
    s,t,st,x,ps,diff : Vector(1..d);
    normdiff : Floating_Number;
    ans : character;
  
  begin
    c(0..d-1) := Multprec_Random_Vectors.Random_Vector(0,d-1,size);
    c(d) := Create(integer(1));
    put_line("The random coefficients of a monic polynomial ");
    put_line(c);
    t := Multprec_Coefficients_to_Traces(c);
    s := Traces_to_Power_Sums(t);
    st := Power_Sums_to_Traces(s);
    put_line("The traces computed from the coefficients : ");
    put_line(t);
    put_line("The traces computed from the power sums : ");
    put_line(st);
    diff := t - st;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
    Clear(diff); Clear(normdiff);
    put("Do you wish to continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    x := Multprec_Roots(Standard_Output,c,size);
    ps := Multprec_Power_Sums(x);
    put_line("Power sums computed from the traces : ");
    put_line(s);
    put_line("Power sums computed from the roots : ");
    put_line(ps);
    diff := s - ps;
    normdiff := Max_Norm(diff);
    put("Max norm of difference : "); put(normdiff,3); new_line;
    Clear(diff); Clear(normdiff);
    Clear(c); Clear(t); Clear(st); Clear(x); Clear(ps);
  end Multprec_Test_Power_Traces;

  procedure Standard_Complete_Reconstruction
               ( file : in file_type; d : in integer32;
                 pts : in Standard_Complex_Vectors.Vector;
                 yval : in Standard_Complex_Matrices.Matrix;
                 t : in Standard_Trace_Interpolators.Trace_Interpolator1 ) is

  -- DESCRIPTION :
  --   Recreates the matrix of y values from the abscisses in the grid
  --   and the trace interpolator.

    use Standard_Complex_Numbers,Standard_Trace_Interpolators;

    tcf : Standard_Complex_Vectors.Vector(0..d);
    y : Standard_Complex_Vectors.Vector(1..d);
    min,dif,maxdif : double_float;
    ind : integer32;

  begin
    tcf(d) := Create(1.0);
    for i in pts'range loop       -- reconstruct for slice x = pts(i)
      for j in 1..d loop
        declare
          trc : constant Standard_Complex_Vectors.Vector(0..j) := Trace(t,j);
        begin
          tcf(d-j) := Standard_Univariate_Interpolators.Evalc(trc,pts(i));
        end;
      end loop;
      y := Standard_Roots(file,tcf);
      put(file,"Comparing roots at slice ");
      put(file,i,1); put_line(file,": ");
      maxdif := 0.0;
      for j in y'range loop        -- compare with values in the matrix
        put(file,y(j));
        ind := yval'first(1);
        min := AbsVal(y(j) - yval(ind,i));
        for k in ind+1..yval'last(1) loop
          dif := AbsVal(y(j) - yval(k,i));
          if dif < min then
            ind := k;
            min := dif;
          end if;
        end loop;
        put(file,"  difference :"); put(file,min,3);
        put(file," at "); put(file,ind,1); new_line(file);
        if min > maxdif
         then maxdif := min;
        end if;
      end loop;
    end loop;
    put(file,"Maximal difference in reconstruction of grid : ");
    put(file,maxdif,3); new_line(file);
  end Standard_Complete_Reconstruction;

  procedure Standard_Traces_Roots
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   This procedure is used in developing tools to create a rectangular
  --   grid of samples from a triangular one.

    use Standard_Complex_Solutions,Standard_Trace_Interpolators;

    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    deg : constant integer32 := integer32(Length_Of(sols));
    grid : Array_of_Standard_Sample_Lists(0..deg);
    eps,dst : double_float;
    pts : Standard_Complex_Vectors.Vector(0..deg);
    yval : Standard_Complex_Matrices.Matrix(1..deg,0..deg);
    t : Trace_Interpolator1;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    Standard_Rectangular_Grid_Creator(file,sps,natural32(deg),grid,eps,dst);
    yval := Rotate_Samples(natural32(deg),natural32(deg),hyp(1).all,grid);
    put_line(file,"The matrix with rotated samples :");
    put(file,yval,3);
    t := Create(grid);
    put(file,"Maximal residual of trace interpolator at grid : ");
    put(file,Maximal_Error(t,grid),3); new_line(file);
    pts := Abscisses(grid,natural32(deg));
    Standard_Complete_Reconstruction(file,deg,pts,yval,t);
    Sampling_Machine.Clear;
  end Standard_Traces_Roots;

  procedure Standard_Triangular_Grid
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

  -- DESCRIPTION :
  --   This procedure is used in developing tools to create a rectangular
  --   grid of samples from a triangular one.

    use Standard_Complex_Solutions,Standard_Trace_Interpolators;

    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    deg : constant integer32 := integer32(Length_Of(sols));
    grid : Array_of_Standard_Sample_Lists(0..deg);
    eps,dist : double_float;
    outputans : character;
    t : Trace_Interpolator1;
    testsps,tmp : Standard_Sample_List;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    put("Do you want intermediate output during creation ? (y/n) ");
    Ask_Yes_or_No(outputans);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Standard_Triangular_Grid_Creator(file,sps,natural32(deg),grid,eps,dist);
    if outputans = 'y'
     then t := Create_on_Triangle(file,grid);
     else t := Create_on_Triangle(grid);
    end if;
    put_line(file,"Errors on the grid : ");
    Write_Errors(file,t,grid);
    new_line(file);
    Standard_Test_Samples(file,sps,hyp,testsps);
    put_line(file,"Evaluating at the test samples : ");
    tmp := testsps;
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Standard_Triangular_Grid;

  procedure Multprec_Triangular_Grid
               ( file : in file_type;
                 ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim,size : in natural32 ) is

  -- DESCRIPTION :
  --   This procedure is used in developing tools to create a rectangular
  --   grid of samples from a triangular one.

    use Standard_Complex_Solutions,Multprec_Trace_Interpolators;

    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(ep,dim);
    sps : Standard_Sample_List := Create(sols,hyp);
    deg : constant natural32 := Length_Of(sols);
    grid : Array_of_Multprec_Sample_Lists(0..integer32(deg));
    eps,dst : double_float;
    outputans : character;
    t : Trace_Interpolator1;
    testsps,tmp : Multprec_Sample_List;
    eva : Multprec_Complex_Numbers.Complex_Number;

  begin
    Sampling_Machine.Initialize(ep,mp,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    put("Do you want intermediate output during creation ? (y/n) ");
    Ask_Yes_or_No(outputans);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Multprec_Triangular_Grid_Creator(file,sps,deg,size,grid,eps,dst);
    if outputans = 'y'
     then t := Create_on_Triangle(file,grid,size);
     else t := Create_on_Triangle(grid,size);
    end if;
    put_line(file,"Errors on the grid : ");
    Write_Errors(file,t,grid);
    new_line(file);
    Multprec_Test_Samples(file,sps,hyp,testsps);
    put_line(file,"Evaluating at the test samples : ");
    tmp := testsps;
    while not Is_Null(tmp) loop
      eva := Eval(t,Sample_Point(Head_Of(tmp)).v);
      put(file,eva); new_line(file);
      Multprec_Complex_Numbers.Clear(eva);
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Multprec_Triangular_Grid;

  procedure Main is

    d : natural32 := 0;
    ans : character;
    size,deci,dim : natural32 := 0;
    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
     
  begin
    loop
      new_line;
      put_line("MENU for testing the operations in power traces :");
      put_line("  1. conversions between power sums and traces;");
      put_line("  2. test relations between traces and roots on grid;");
      put_line("  3. triangular grid creation with Newton identities.");
      put("Type 1,2 or 3 to select : "); Ask_Alternative(ans,"123");
      if ans = '1' then
        new_line;
        put("Give the degree : "); get(d);
        new_line;
        put_line("MENU for the precision of the computations : ");
        put_line("  0. standard double precision;");
        put_line("  1. double double precision;");
        put_line("  2. quad double precision;");
        put_line("  3. arbitrary multiprecision;");
        put("Type 0, 1, 2, or 3 to select precision : ");
        Ask_Alternative(ans,"0123");
        case ans is
           when '0' => Standard_Test_Power_Traces(integer32(d));
           when '1' => DoblDobl_Test_Power_Traces(integer32(d));
           when '2' => QuadDobl_Test_Power_Traces(integer32(d));
           when '3' =>
             put("Give the number of decimal places : ");
             get(deci);
             size := Decimal_to_Size(deci);
             Multprec_Test_Power_Traces(integer32(d),size);
           when others => null;
        end case;
        put("Do you want more tests ? (y/n) ");
        Ask_Yes_or_No(ans);
      else
        Standard_Read_Embedding(lp,sols,dim);
        new_line;
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(file);
        if ans = '2' then
          Standard_Traces_Roots(file,lp.all,sols,dim);
        else
          new_line;
          put("Give the number of decimal places" 
              & " (<= 16 is standard) : ");
          get(deci);
          if deci > 16 then
            size := Decimal_to_Size(deci);
            Get_Multprec_System(lp.all,mp,size,dim);
            Multprec_Triangular_Grid(file,lp.all,mp.all,sols,dim,size);
          else
            Standard_Triangular_Grid(file,lp.all,sols,dim);
          end if;
        end if;
      end if;
      exit when (ans /= 'y');
    end loop;
  end Main;

begin
  new_line;
  put_line("Newton identities to convert power sums into traces, and back.");
  Main;
end ts_powtrc;
