with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Mathematical_Functions;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Standard_Random_Vectors;
with Standard_Vector_Splitters;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;
with Standard_Matrix_Splitters;
with Standard_Complex_BLAS_Helpers;
with Standard_Inlined_BLAS_Helpers;
with Standard_Complex_Singular_Values;
with Standard_Inlined_Singular_Values;

procedure ts_perfdsvd is

-- DESCRIPTION :
--   Test the development of a better performing SVD implementation,
--   in double precision.

  procedure Test_Euclidean_Norm_of_Vector ( n : in integer32 ) is

  -- DESCRIPTION :
  --    Generates a random n-dimensional vector and compares
  --    the output of the dznrm2 on the complex vector with
  --    the output on the splitted vector.

    x : constant Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    nrm1,nrm2 : double_float;

  begin
    Standard_Vector_Splitters.Split_Complex(x,xre,xim);
    nrm1 := Standard_Complex_BLAS_Helpers.dznrm2(n,x,1,1);
    nrm2 := Standard_Inlined_BLAS_Helpers.dznrm2(n,xre,xim,1,1);
    put("nrm1 : "); put(nrm1); new_line;
    put("nrm2 : "); put(nrm2); new_line;
  end Test_Euclidean_Norm_of_Vector;

  procedure Test_Vector_Scaling ( n : in integer32 ) is

  -- DESCRIPTION :
  --    Generates a random n-dimensional vector and compares
  --    the output of the zscal on the complex vector with
  --    the output on the splitted vector.

    x : Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    z : constant Complex_Number := Standard_Random_Numbers.Random1;
    zre : constant double_float := REAL_PART(z);
    zim : constant double_float := IMAG_PART(z);
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;

  begin
    Standard_Vector_Splitters.Split_Complex(x,xre,xim);
    Standard_Complex_BLAS_Helpers.zscal(n,z,x,1,1);
    put_line("zscal on the complex vector :"); put_line(x);
    Standard_Inlined_BLAS_Helpers.zscal(n,zre,zim,xre,xim,1,1);
    put_line("real parts of zscal on splitted vector :"); put_line(xre);
    put_line("imaginary parts of zscal on splitted vector :"); put_line(xim);
  end Test_Vector_Scaling;

  procedure Test_Euclidean_Norm_of_Column ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   dznrm2 on the complex matrix with the output of the splitted matrix.

    A : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
    nrm1,nrm2 : double_float;

  begin
    rvv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    ivv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(A,rvv,ivv);
    nrm1 := Standard_Complex_BLAS_Helpers.dznrm2(n,A,1,1,1);
    nrm2 := Standard_Inlined_BLAS_Helpers.dznrm2(n,rvv,ivv,1,1,1);
    put("nrm1 : "); put(nrm1); new_line;
    put("nrm2 : "); put(nrm2); new_line;
  end Test_Euclidean_Norm_of_Column;

  procedure Test_Column_Scaling ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   zscal on the complex matrix with the output of the splitted matrix.

    A : Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
    z : constant Complex_Number := Standard_Random_Numbers.Random1;
    zre : constant double_float := REAL_PART(z);
    zim : constant double_float := IMAG_PART(z);

  begin
    rvv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    ivv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(A,rvv,ivv);
    Standard_Complex_BLAS_Helpers.zscal(n,z,A,1,1,1);
    Standard_Inlined_BLAS_Helpers.zscal(n,zre,zim,rvv,ivv,1,1,1);
    put_line("zscal on the complex column : ");
    for i in 1..n loop
      put(A(i,1)); new_line;
    end loop;
    put_line("real parts of zscal on splitted column : ");
    put_line(rvv(1));
    put_line("imaginary parts of zscal on splitted column : ");
    put_line(ivv(1));
  end Test_Column_Scaling;

  procedure Test_Update_Column_with_Vector ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   zaxpy on the complex matrix with the output of the splitted matrix.

    x : constant Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    y : Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
    z : constant Complex_Number := Standard_Random_Numbers.Random1;
    zre : constant double_float := REAL_PART(z);
    zim : constant double_float := IMAG_PART(z);

  begin
    Standard_Vector_Splitters.Split_Complex(x,xre,xim);
    rvv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    ivv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(y,rvv,ivv);
    Standard_Complex_BLAS_Helpers.zaxpy(n,z,x,1,1,y,1,1,1);
    Standard_Inlined_BLAS_Helpers.zaxpy(n,zre,zim,xre,xim,1,1,rvv,ivv,1,1,1);
    put_line("zaxpy on the complex column : ");
    for i in 1..n loop
      put(y(i,1)); new_line;
    end loop;
    put_line("real parts of zaxpy on splitted column : ");
    put_line(rvv(1));
    put_line("imaginary parts of zaxpy on splitted column : ");
    put_line(ivv(1));
  end Test_Update_Column_with_Vector;

  procedure Test_Update_Vector_with_Column ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   zaxpy on the complex matrix with the output of the splitted matrix.

    x : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    y : Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    yre,yim : Standard_Floating_Vectors.Link_to_Vector;
    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
    z : constant Complex_Number := Standard_Random_Numbers.Random1;
    zre : constant double_float := REAL_PART(z);
    zim : constant double_float := IMAG_PART(z);

  begin
    rvv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    ivv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(x,rvv,ivv);
    Standard_Vector_Splitters.Split_Complex(y,yre,yim);
    Standard_Complex_BLAS_Helpers.zaxpy(n,z,x,1,1,1,y,1,1);
    Standard_Inlined_BLAS_Helpers.zaxpy(n,zre,zim,rvv,ivv,1,1,1,yre,yim,1,1);
    put_line("zaxpy on the complex vector : "); put_line(y);
    put_line("real parts of zaxpy on splitted column : ");
    put_line(yre);
    put_line("imaginary parts of zaxpy on splitted column : ");
    put_line(yim);
  end Test_Update_Vector_with_Column;

  procedure Test_Update_Column_with_Column ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   zaxpy on the complex matrix with the output of the splitted matrix.

    x : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    y : Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec;
    yrv,yiv : Standard_Floating_VecVecs.Link_to_VecVec;
    z : constant Complex_Number := Standard_Random_Numbers.Random1;
    zre : constant double_float := REAL_PART(z);
    zim : constant double_float := IMAG_PART(z);

  begin
    xrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    xiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(x,xrv,xiv);
    Standard_Matrix_Splitters.Complex_Parts(y,yrv,yiv);
    Standard_Complex_BLAS_Helpers.zaxpy(n,z,x,1,1,1,y,1,1,1);
    Standard_Inlined_BLAS_Helpers.zaxpy(n,zre,zim,xrv,xiv,1,1,1,yrv,yiv,1,1,1);
    put_line("zaxpy on the complex vector : ");
    for i in 1..n loop
      put(y(i,1)); new_line;
    end loop;
    put_line("real parts of zaxpy on splitted column : ");
    put_line(yrv(1));
    put_line("imaginary parts of zaxpy on splitted column : ");
    put_line(yiv(1));
  end Test_Update_Column_with_Column;

  procedure Test_Dot_Product ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   zdotc on the complex matrix with the output of the splitted matrix.

    x : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    y : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec;
    yrv,yiv : Standard_Floating_VecVecs.Link_to_VecVec;
    z : Complex_Number;
    zre,zim : double_float;

  begin
    xrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    xiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(x,xrv,xiv);
    Standard_Matrix_Splitters.Complex_Parts(y,yrv,yiv);
    z := Standard_Complex_BLAS_Helpers.zdotc(n,x,1,1,1,y,1,1,1);
    Standard_Inlined_BLAS_Helpers.zdotc(n,xrv,xiv,1,1,1,yrv,yiv,1,1,1,zre,zim);
    put_line("zdotc on the complex vector : "); put(z); new_line;
    put_line("zdotc on the splitted vector : ");
    put(zre); put("  "); put(zim); new_line;
  end Test_Dot_Product;

  procedure Test_Rotation ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   zdotc on the complex matrix with the output of the splitted matrix.

    x : Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    y : Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec;
    yrv,yiv : Standard_Floating_VecVecs.Link_to_VecVec;
    angle : constant double_float := Standard_Random_Numbers.Random; 
    pi : constant double_float := Standard_Mathematical_Functions.PI;
    c : constant double_float := Standard_Mathematical_Functions.COS(pi*angle);
    s : constant double_float := Standard_Mathematical_Functions.SIN(pi*angle);

  begin
    xrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    xiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(x,xrv,xiv);
    Standard_Matrix_Splitters.Complex_Parts(y,yrv,yiv);
    Standard_Complex_BLAS_Helpers.zdrot(n,x,1,1,1,y,1,1,1,c,s);
    Standard_Inlined_BLAS_Helpers.zdrot(n,xrv,xiv,1,1,1,yrv,yiv,1,1,1,c,s);
    put_line("zdrot on the first complex column : ");
    for i in 1..n loop
      put(x(i,1)); new_line;
    end loop;
    put_line("real parts of zdrot on the first splitted vector : ");
    put_line(xrv(1));
    put_line("imaginary parts of zdrot on the first splitted vector : ");
    put_line(xiv(1));
    put_line("zdrot on the second complex column : ");
    for i in 1..n loop
      put(y(i,1)); new_line;
    end loop;
    put_line("real parts of zdrot on the second splitted vector : ");
    put_line(yrv(1));
    put_line("imaginary parts of zdrot on the second splitted vector : ");
    put_line(yiv(1));
  end Test_Rotation;

  procedure Test_Swap ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   zswap on the complex matrix with the output of the splitted matrix.

    x : Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    y : Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec;
    yrv,yiv : Standard_Floating_VecVecs.Link_to_VecVec;

  begin
    xrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    xiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    yiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(x,xrv,xiv);
    Standard_Matrix_Splitters.Complex_Parts(y,yrv,yiv);
    Standard_Complex_BLAS_Helpers.zswap(n,x,1,1,1,y,1,1,1);
    Standard_Inlined_BLAS_Helpers.zswap(n,xrv,xiv,1,1,1,yrv,yiv,1,1,1);
    put_line("zswap on the first complex column : ");
    for i in 1..n loop
      put(x(i,1)); new_line;
    end loop;
    put_line("real parts of zswap on the first splitted vector : ");
    put_line(xrv(1));
    put_line("imaginary parts of zswap on the first splitted vector : ");
    put_line(xiv(1));
    put_line("zswap on the second complex column : ");
    for i in 1..n loop
      put(y(i,1)); new_line;
    end loop;
    put_line("real parts of zswap on the second splitted vector : ");
    put_line(yrv(1));
    put_line("imaginary parts of zswap on the second splitted vector : ");
    put_line(yiv(1));
  end Test_Swap;

  procedure Test ( n,p : in integer32;
                   A : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in A is an n-by-p matrix, computes its SVD,
  --   and compares the output with the inlined SVD output.

    wrk : Standard_Complex_Matrices.Matrix(1..n,1..p) := A;
    wrk2 : Standard_Complex_Matrices.Matrix(1..n,1..p) := A;
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    S,S2 : Standard_Complex_Vectors.Vector(1..mm);
    e,e2 : Standard_Complex_Vectors.Vector(1..p);
    U,U2 : Standard_Complex_Matrices.Matrix(1..n,1..n);
    V,V2 : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Standard_Complex_Singular_Values.SVD(wrk,n,p,S,e,U,V,job,info);
    Standard_Inlined_Singular_Values.SVD(wrk2,n,p,S2,e2,U2,V2,job,info);
    put_line("The singular values : "); put_line(S);
    put_line("The recomputed singular values : "); put_line(S2);
    put_line("The matrix U :"); put(U,2);
    put_line("The recomputed matrix U :"); put(U2,2);
    put_line("The matrix V :"); put(V,2);
    put_line("The recomputed matrix V :"); put(V2,2);
  end Test;

  procedure Test_SVD ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p matrix and test the SVD.

    A : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;

  begin
    put("A random "); put(n,1); put("-by-"); put(p,1);
    put_line(" complex matrix :"); put(A,2);
    rvv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    ivv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(A,rvv,ivv);
    put_line("The real parts of the columns :"); put(rvv,2);
    put_line("The imaginary parts of the columns :"); put(ivv,2);
    Test(n,p,A);
  end Test_SVD;

  procedure Time_SVD ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p matrix and prompts the user
  --   for a frequence to test the time of SVD computations.

    timer : Timing_Widget;
    A : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    frq : integer32 := 0;
    wrk,wrk2 : Standard_Complex_Matrices.Matrix(1..n,1..p);
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    S,S2 : Standard_Complex_Vectors.Vector(1..mm);
    e,e2 : Standard_Complex_Vectors.Vector(1..p);
    U,U2 : Standard_Complex_Matrices.Matrix(1..n,1..n);
    V,V2 : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : integer32 := 11;
    info : integer32;
    xrv,xiv : Standard_Floating_VecVecs.Link_to_VecVec;
    urv,uiv : Standard_Floating_VecVecs.Link_to_VecVec;
    vrv,viv : Standard_Floating_VecVecs.Link_to_VecVec;
    sr,si,er,ei,wr,wi : Standard_Floating_Vectors.Link_to_Vector;
    ans : character;

  begin
    put("Give the frequency : "); get(frq);
    put("Compute only singular values ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then job := 0;
    end if;
    tstart(timer);
    for k in 1..frq loop
      wrk := A;
      Standard_Complex_Singular_Values.SVD(wrk,n,p,S,e,U,V,job,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex SVD");
    tstart(timer);
    for k in 1..frq loop
      wrk2 := A;
      Standard_Inlined_Singular_Values.SVD(wrk2,n,p,S2,e2,U2,V2,job,info);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"wrapped SVD");
    xrv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    xiv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    if job /= 0 then
      urv := Standard_Vector_Splitters.Allocate(n,n,1,1);
      uiv := Standard_Vector_Splitters.Allocate(n,n,1,1);
      vrv := Standard_Vector_Splitters.Allocate(p,p,1,1);
      viv := Standard_Vector_Splitters.Allocate(p,p,1,1);
    end if;
    sr := new Standard_Floating_Vectors.Vector'(1..mm => 0.0);
    si := new Standard_Floating_Vectors.Vector'(1..mm => 0.0);
    er := new Standard_Floating_Vectors.Vector'(1..p => 0.0);
    ei := new Standard_Floating_Vectors.Vector'(1..p => 0.0);
    wr := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    wi := new Standard_Floating_Vectors.Vector'(1..n => 0.0);
    tstart(timer);
    for k in 1..frq loop -- only do relevant parts/merge
      Standard_Matrix_Splitters.Complex_Parts(A,xrv,xiv);
      Standard_Inlined_Singular_Values.SVD
        (xrv,xiv,n,p,sr,si,er,ei,urv,uiv,vrv,viv,job,info,wr,wi);
      Standard_Vector_Splitters.Complex_Merge(sr,si,s);
      if job /= 0 then
        Standard_Matrix_Splitters.Complex_Merge(urv,uiv,u);
        Standard_Matrix_Splitters.Complex_Merge(vrv,viv,v);
      end if;
    end loop;
    tstop(timer);
    Standard_Floating_Vectors.Clear(sr);
    Standard_Floating_Vectors.Clear(si);
    Standard_Floating_Vectors.Clear(er);
    Standard_Floating_Vectors.Clear(ei);
    Standard_Floating_Vectors.Clear(wr);
    Standard_Floating_Vectors.Clear(wi);
    Standard_Floating_VecVecs.Deep_Clear(xrv);
    Standard_Floating_VecVecs.Deep_Clear(xiv);
    if job /= 0 then
      Standard_Floating_VecVecs.Deep_Clear(urv);
      Standard_Floating_VecVecs.Deep_Clear(uiv);
      Standard_Floating_VecVecs.Deep_Clear(vrv);
      Standard_Floating_VecVecs.Deep_Clear(viv);
    end if;
    new_line;
    print_times(standard_output,timer,"inlined SVD");
  end Time_SVD;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns
  --   and then generates a random complex matrix
  --   with that number of rows and columns.

    ans : character;

  begin
    new_line;
    put_line("MENU to test inlined SVD operations :");
    put_line("  1. Euclidean norm of a random vector; ");
    put_line("  2. scaling a random vector; ");
    put_line("  3. Euclidean norm of a column of a random matrix; ");
    put_line("  4. scaling a column of a random matrix; ");
    put_line("  5. test column update with multiple of a vector;");
    put_line("  6. test vector update with multiple of a column;");
    put_line("  7. test column update with multiple of a column;");
    put_line("  8. test dot product of two columns;");
    put_line("  9. test rotation of two columns;");
    put_line("  A. test swap of two columns.");
    put_line("  B. SVD of a random matrix.");
    put_line("  C. time the SVD of a random matrix.");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, or C to select a test : ");
    Ask_Alternative(ans,"123456789ABC");
    new_line;
    case ans is
      when '1' | '2' =>
        declare
          dim : integer32 := 0;
        begin
          put("Give the dimension : "); get(dim);
          if ans = '1'
           then Test_Euclidean_Norm_of_Vector(dim);
           else Test_Vector_Scaling(dim);
          end if;
        end;
      when '3' | '4' | '5' | '6' | '7' | '8' | '9' | 'A' | 'B' | 'C' =>
        declare
          nrows,ncols : integer32 := 0;
        begin
          put("Give the number of rows : "); get(nrows);
          put("Give the number of columns : "); get(ncols);
          case ans is
            when '3' => Test_Euclidean_Norm_of_Column(nrows,ncols);
            when '4' => Test_Column_Scaling(nrows,ncols);
            when '5' => Test_Update_Column_with_Vector(nrows,ncols);
            when '6' => Test_Update_Vector_with_Column(nrows,ncols);
            when '7' => Test_Update_Column_with_Column(nrows,ncols);
            when '8' => Test_Dot_Product(nrows,ncols);
            when '9' => Test_Rotation(nrows,ncols);
            when 'A' => Test_Swap(nrows,ncols);
            when 'B' => Test_SVD(nrows,ncols);
            when 'C' => Time_SVD(nrows,ncols);
            when others => null;
          end case;
        end;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_perfdsvd;
