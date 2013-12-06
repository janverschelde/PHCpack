with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;          use QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;      use Standard_Floating_VecVecs_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Random_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;       use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Random_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;       use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Random_VecVecs;
with Standard_Floating_GramSchmidt;
with Standard_Complex_GramSchmidt;
with DoblDobl_Complex_GramSchmidt;
with QuadDobl_Complex_GramSchmidt;

procedure ts_mgsrnk is

-- DESCRIPTION :
--   Development of modified Gram-Schmidt orthonormalization
--   to compute the rank and solve in the least squares sense.

  procedure Show_Accuracy_Report ( orterr,eqserr : in double_float ) is

  -- DESCRIPTION :
  --   Shows the largest deviation in orthonormality and decomposition,
  --   given respectively in orterr and eqserr.

  begin
    put_line("Accuracy level report :");
    put("  largest error in orthonormality :"); put(orterr,3); new_line;
    put("  largest error in decomposition  :"); put(eqserr,3); new_line;
  end Show_Accuracy_Report;

  procedure Projection_of_Space
               ( n,m : in integer32;
                 v,q : in Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Shows the result of Q^T*v, where Q is the matrix spanned by
  --   the m columns of q, which are vectors of range 1..n.
  --   The results of this project should be R, as Q^T A = R,
  --   where R is the upper triangular matrix of the QR decomposition.
  --   This is a check to compare with the outcome of Q^T A with
  --   the result of R after Permute_Upper_Triangular(n,r,pivots).

    qTv : Standard_Floating_Vectors.Vector(1..m);

  begin
    put_line("Computing R via Q^T*v : ");
    for i in 1..m loop
      qTv := Standard_Floating_GramSchmidt.Matrix_Projection(n,m,q,v(i).all);
      put("the projection of vector "); put(i,1); put_line(" :");
      put_line(qTv);
    end loop;
  end Projection_of_Space;

  procedure Standard_Floating_Random_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method
  --   with standard floating-point arithmetic.
 
    use Standard_Floating_GramSchmidt;

    v : constant Standard_Floating_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : Standard_Floating_VecVecs.VecVec(1..m);
    pivots : Standard_Integer_Vectors.Vector(1..m);
    tol : constant double_float := 1.0E-8;
    rank : integer32;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new Standard_Floating_Vectors.Vector'(v(i).all);
      r(i) := new Standard_Floating_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := 0.0;
      end loop;
    end loop;
    QR(n,m,tol,q,r,pivots,rank);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    put("the pivots :"); put(pivots); new_line;
    Projection_of_Space(n,m,v,q);
    Permute_Upper_Triangular(n,r,pivots);
    put_line("The computed R (after permuting) : ");
    for i in 1..m loop
      put("column "); put(i,1); put_line(" of the matrix R :");
      put_line(r(i).all);
    end loop;
   -- Revert_Swaps(n,q,pivots);
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end Standard_Floating_Random_Test;

  procedure Standard_Complex_Random_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method
  --   with standard complex arithmetic.
 
    use Standard_Complex_GramSchmidt;

    v : constant Standard_Complex_VecVecs.VecVec(1..m)
      := Standard_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : Standard_Complex_VecVecs.VecVec(1..m);
    pivots : Standard_Integer_Vectors.Vector(1..m);
    tol : constant double_float := 1.0E-8;
    rank : integer32;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new Standard_Complex_Vectors.Vector'(v(i).all);
      r(i) := new Standard_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(0.0);
      end loop;
    end loop;
    QR(n,m,tol,q,r,pivots,rank);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    put("the pivots :"); put(pivots); new_line;
   -- Projection_of_Space(n,m,v,q);
    Permute_Upper_Triangular(n,r,pivots);
    put_line("The computed R (after permuting) : ");
    for i in 1..m loop
      put("column "); put(i,1); put_line(" of the matrix R :");
      put_line(r(i).all);
    end loop;
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end Standard_Complex_Random_Test;

  procedure DoblDobl_Complex_Random_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method
  --   with double double complex arithmetic.
 
    use DoblDobl_Complex_GramSchmidt;

    v : constant DoblDobl_Complex_VecVecs.VecVec(1..m)
      := DoblDobl_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : DoblDobl_Complex_VecVecs.VecVec(1..m);
    pivots : Standard_Integer_Vectors.Vector(1..m);
    tol : constant double_float := 1.0E-8;
    rank : integer32;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new DoblDobl_Complex_Vectors.Vector'(v(i).all);
      r(i) := new DoblDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    QR(n,m,tol,q,r,pivots,rank);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    put("the pivots :"); put(pivots); new_line;
   -- Projection_of_Space(n,m,v,q);
    Permute_Upper_Triangular(n,r,pivots);
    put_line("The computed R (after permuting) : ");
    for i in 1..m loop
      put("column "); put(i,1); put_line(" of the matrix R :");
      put_line(r(i).all);
    end loop;
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end DoblDobl_Complex_Random_Test;

  procedure QuadDobl_Complex_Random_Test ( n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random vectors of range 1..n
  --   and computes the QR decomposition via the
  --   modified Gram-Schmidt orthonormalization method
  --   with quad double complex arithmetic.
 
    use QuadDobl_Complex_GramSchmidt;

    v : constant QuadDobl_Complex_VecVecs.VecVec(1..m)
      := QuadDobl_Random_VecVecs.Random_VecVec(natural32(n),natural32(m));
    q,r : QuadDobl_Complex_VecVecs.VecVec(1..m);
    pivots : Standard_Integer_Vectors.Vector(1..m);
    tol : constant double_float := 1.0E-8;
    rank : integer32;
    orterr,eqserr : double_float;
    fail : boolean;

  begin
    put_line(v);
    for i in 1..m loop
      q(i) := new QuadDobl_Complex_Vectors.Vector'(v(i).all);
      r(i) := new QuadDobl_Complex_Vectors.Vector(1..m);
      for j in 1..m loop
        r(i)(j) := Create(integer(0));
      end loop;
    end loop;
    QR(n,m,tol,q,r,pivots,rank);
    Test_Orthonormality(n,m,q,tol,true,orterr,fail);
    put("the largest error : "); put(orterr,3); new_line;
    if fail
     then put_line("Orthogonality test failed!");
     else put_line("Passed the orthogonality test.");
    end if;
    put("the pivots :"); put(pivots); new_line;
   -- Projection_of_Space(n,m,v,q);
    Permute_Upper_Triangular(n,r,pivots);
    put_line("The computed R (after permuting) : ");
    for i in 1..m loop
      put("column "); put(i,1); put_line(" of the matrix R :");
      put_line(r(i).all);
    end loop;
    Test_Decomposition(n,m,v,q,r,tol,true,eqserr,fail);
    put("the largest error : "); put(eqserr,3); new_line;
    if fail
     then put_line("Decomposition test failed!");
     else put_line("Passed the decomposition test.");
    end if;
    Show_Accuracy_Report(orterr,eqserr);
  end QuadDobl_Complex_Random_Test;

  procedure Main is

    n,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("MENU for modified Gram-Schmidt method with pivoting :");
    put_line("  1. standard floating point numbers");
    put_line("  2. standard complex floating point numbers");
    put_line("  3. double double complex floating point numbers");
    put_line("  4. quad double complex floating point numbers");
    put("Type 1, 2, 3, or 4 to choose test : "); Ask_Alternative(ans,"1234");
    new_line;
    put("Give dimension : "); get(n);
    put("Give number of vectors : "); get(m);
    new_line;
    case ans is
      when '1' => Standard_Floating_Random_Test(n,m);
      when '2' => Standard_Complex_Random_Test(n,m);
      when '3' => DoblDobl_Complex_Random_Test(n,m);
      when '4' => QuadDobl_Complex_Random_Test(n,m);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mgsrnk;
