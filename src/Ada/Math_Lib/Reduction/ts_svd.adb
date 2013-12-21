with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;     use Multprec_Complex_Number_Tools;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Tools;     use Multprec_Complex_Vector_Tools;
with Multprec_Complex_Matrices;
with Multprec_Complex_Matrices_io;      use Multprec_Complex_Matrices_io;
with Multprec_Random_Vectors;           use Multprec_Random_Vectors;
with Multprec_Random_Matrices;          use Multprec_Random_Matrices;
with Multprec_Complex_Singular_Values;  use Multprec_Complex_Singular_Values;

procedure ts_svd is

-- DESCRIPTION :
--   Interactive testing facility for Singular Value Decomposition,
--   with standard and multiprecision complex arithmetic.

-- ROUTINES USING STANDARD ARITHMETIC :

  function Read_Vector
             ( n : in integer32 ) return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user to enter n complex numbers,
  --   returns the vector which contains the user given numbers.

    res : Standard_Complex_Vectors.Vector(1..n);

  begin
    put("Give "); put(n,1);
    put(" complex numbers for a vector of dimension ");
    put(n,1); put_line(" :");
    get(res);
    return res;
  end Read_Vector;

  function Read_Matrix
             ( n,m : in integer32 ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Prompts the user to enter n*m complex numbers to define
  --   an n-by-m matrix.  Returns the matrix of given numbers.

    res : Standard_Complex_Matrices.Matrix(1..n,1..m);

  begin
    put("Give "); put(n*m,1); put(" complex numbers for ");
    put(n,1); put("x"); put(m,1); put_line(" matrix :");
    get(res);
    return res;
  end Read_Matrix;

  function Is_Identity ( a : Standard_Complex_Matrices.Matrix;
                         tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if a is the identity matrix modulo the tolerance.

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        if i = j then
          if Absval(a(i,j) - 1.0) > tol
           then return false;
          end if;
        else
          if AbsVal(a(i,j)) > tol
           then return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Is_Identity;

  function Is_Orthogonal ( a : Standard_Complex_Matrices.Matrix;
                           tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if a*a^H and a^H*a equal the identity matrix
  --   modulo the given tolerance tol.

    use Standard_Complex_Matrices;

    atr : constant Matrix(a'range(2),a'range(1)) := Conjugate_Transpose(a);
    ata : constant Matrix(a'range(2),a'range(2)) := atr*a;
    aat : constant Matrix(a'range(1),a'range(1)) := a*atr;

  begin
    if not Is_Identity(ata,tol) then
      return false;
    elsif not Is_Identity(aat,tol) then
      return false;
    else
      return true;
    end if;
  end Is_Orthogonal;

  function Is_SVD ( x,u,v : Standard_Complex_Matrices.Matrix;
                    s : Standard_Complex_Vectors.Vector;
                    tol : in double_float ) return boolean is

  -- DESCRIPTION :
  --   Forms the product u'*x*v and checks whether this product is a
  --   diagonal matrix with the singular values on the diagonal.

    use Standard_Complex_Matrices;

    ut : constant Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u);
    utx : constant Matrix(u'range(2),x'range(2)) := ut*x;
    utxv : constant Matrix(u'range(2),v'range(2)) := utx*v;

    procedure Error is
    begin
      put_line("The matrix u'*x*v : "); put(utxv,3);
      put_line("is NOT a singular value decomposition!  BUG!");
    end Error;

  begin
    for i in utxv'range(1) loop
      for j in utxv'range(2) loop
        if i = j then
          if AbsVal(utxv(i,j) - s(i)) > tol
           then Error; return false;
          end if;
        else
          if AbsVal(utxv(i,j)) > tol
           then Error; return false;
          end if;
        end if;
      end loop;
    end loop;
    return true;
  end Is_SVD;

  procedure Test_SVD_Output
              ( x,u,v : in Standard_Complex_Matrices.Matrix;
                s,e : in Standard_Complex_Vectors.Vector;
                info : in integer32 ) is

  -- DESCRIPTION :
  --   Tests whether the matrices u and v are orthogonal and
  --   whether u'*x*v yields a diagonal matrix with the singular
  --   values on the diagonal.  Diagnostics are written to screen.

    ortho,decomp : boolean;
    tol : constant double_float := 1.0E-8;

  begin
    put("info = "); put(info,1); new_line;
    put_line("The e values : "); put_line(e);
    put_line("The singular values : "); put_line(s);
    put_line("The matrix u : "); put(u);
    ortho := Is_Orthogonal(u,tol);
    if ortho
     then put_line("The matrix u is orthogonal.");
     else put_line("The matrix u is NOT orthogonal!  BUG!!!");
    end if;
    put_line("The matrix v : "); put(v);
    ortho := Is_Orthogonal(v,tol);
    if ortho
     then put_line("The matrix v is orthogonal.");
     else put_line("The matrix v is NOT orthogonal!  BUG!!!");
    end if;
    decomp := Is_SVD(x,u,v,s,tol);
    if decomp then
      put_line("The singular values stand u'*x*v test.");
    else
      put_line("The singular values do NOT stand u'*x*v test!");
      put_line("The matrix x :"); put(x,3);
      put_line("The computed matrix u :"); put(u,3);
      put_line("The computed matrix v :"); put(v,3);
    end if;
  end Test_SVD_Output;

  procedure Test_SVD_Solver
               ( a,u,v : in Standard_Complex_Matrices.Matrix;
                 s,b : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the solver using the SVD of a to solve a*x = b,
  --   and prints the residual vector b - a*x.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

    res : Vector(b'range) := b;
    x : Vector(a'range(2));
    inv : constant Matrix(a'range(2),a'range(1)) := Inverse(u,v,s);
    ai : constant Matrix(a'range(1),a'range(1)) := a*inv;
    ia : constant Matrix(a'range(2),a'range(2)) := inv*a;
    tol : constant double_float := 1.0E-8;

  begin
    put("The rank of the matrix : ");
    put(Rank(s),1); new_line;
    put("Inverse of condition number : ");
    put(Inverse_Condition_Number(s),3); new_line;
    if Is_Identity(ai,tol) then
      put_line("The product a*inv(a) is the indentity matrix.");
    else
      put_line("The product a*inv(a) is NOT the indentity matrix!");
      if a'length(1) > a'length(2) then
        put_line("... okay, inv(a) can only be a left inverse.");
      else
        put_line("This may be a genuine bug, look at the product:");
        put(ai,3);
      end if;
    end if;
    if Is_Identity(ia,tol) then
      put_line("The product inv(a)*a is the indentity matrix.");
    else
      put_line("The product inv(a)*a is NOT the indentity matrix!");
      if a'length(1) < a'length(2) then
        put_line("... okay, inv(a) can only be a right inverse.");
      else
        put_line("This may be a genuine bug, look at the product:");
        put(ia,3);
      end if;
    end if;
    x := Solve(u,v,s,b);
    res := b - a*x;
    put_line("The residual vector : ");
    put_line(res);
  end Test_SVD_Solver;

  procedure Test_SVD_on_Given_Matrix ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Computes the SVD of a user given n-by-p matrix.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

    x : Matrix(1..n,1..p) := Read_Matrix(n,p);
    y : constant Matrix(1..n,1..p) := x;
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v: Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    SVD(x,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
  end Test_SVD_on_Given_Matrix;

  procedure Test_SVD_on_Given_System ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Uses SVD to solve a given linear system of n equations in p unknowns.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

    a : Matrix(1..n,1..p) := Read_Matrix(n,p);
    b : constant Vector(1..n) := Read_Vector(n);
    y : constant Matrix(1..n,1..p) := a;
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v: Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    SVD(a,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
    Test_SVD_Solver(y,u,v,s,b);
  end Test_SVD_on_Given_System;

  procedure Test_SVD_on_Random_Matrix ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Computes a SVD on a random n-by-p matrix.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
  
    x : Matrix(1..n,1..p) := Random_Matrix(natural32(n),natural32(p));
    y : constant Matrix(1..n,1..p) := x;
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v: Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    SVD(x,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
  end Test_SVD_on_Random_Matrix;

  procedure Test_SVD_on_Random_System ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

    a : Matrix(1..n,1..p) := Random_Matrix(natural32(n),natural32(p));
    b : constant Vector(1..n) := Random_Vector(1,n);
    y : constant Matrix(1..n,1..p) := a;
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v: Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    SVD(a,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
    Test_SVD_Solver(y,u,v,s,b);
  end Test_SVD_on_Random_System;

-- ROUTINES USING MULTIPRECISION ARITHMETIC :

  procedure Set_Size ( a : in out Multprec_Complex_Matrices.Matrix;
                       size : in natural32 ) is

  -- DESCRIPTION :
  --   Converts the numbers in a to the precision given by size.

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        Set_Size(a(i,j),size);
      end loop;
    end loop;
  end Set_Size;

  function Read_Vector
             ( n : in integer32 ) return Multprec_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Prompts the user to enter n complex numbers to define
  --   an n-vector.  Returns the vector of given numbers.

    res : Multprec_Complex_Vectors.Vector(1..n);

  begin
    put("Give "); put(n,1); put(" complex numbers for a ");
    put(n,1); put_line("-vector :");
    get(res);
    return res;
  end Read_Vector;

  function Read_Matrix
             ( n,m : in integer32 ) return Multprec_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Prompts the user to enter n*m complex numbers to define
  --   an n-by-m matrix.  Returns the matrix of given numbers.

    res : Multprec_Complex_Matrices.Matrix(1..n,1..m);

  begin
    put("Give "); put(n*m,1); put(" complex numbers for ");
    put(n,1); put("x"); put(m,1); put_line(" matrix :");
    get(res);
    return res;
  end Read_Matrix;

  function Is_Identity ( a : Multprec_Complex_Matrices.Matrix;
                         tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the matrix a is the identity matrix,
  --   modulo the given tolerance tol.

    val : Floating_Number;
    rnd : double_float;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        val := AbsVal(a(i,j));
        rnd := Round(val);
        if i = j then
          if abs(rnd - 1.0) > tol
           then return false;
          end if;
        else
          if rnd > tol
           then return false;
          end if;
        end if;
        Clear(val);
      end loop;
    end loop;
    return true;
  end Is_Identity;

  function Is_Orthogonal ( a : Multprec_Complex_Matrices.Matrix;
                           tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if a*a^H and a^H*a equal the identity matrix,
  --   modulo the given tolerance tol.

    use Multprec_Complex_Matrices;

    res : boolean;
    atr : Matrix(a'range(2),a'range(1)) := Conjugate_Transpose(a);
    ata : Matrix(a'range(2),a'range(2)) := atr*a;
    aat : Matrix(a'range(1),a'range(1)) := a*atr;

  begin
   -- put_line("The matrix ata : "); put(ata); 
    if not Is_Identity(ata,tol) then
      res := false;
    elsif not Is_Identity(aat,tol) then
      res := false;
    else
      res := true;
    end if;
    Multprec_Complex_Matrices.Clear(atr);
    Multprec_Complex_Matrices.Clear(ata);
    Multprec_Complex_Matrices.Clear(aat);
    return res;
  end Is_Orthogonal;

  function Is_SVD ( x,u,v : Multprec_Complex_Matrices.Matrix;
                    s : Multprec_Complex_Vectors.Vector;
                    tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if u^H*x*v is a diagonal matrix with the singular
  --   values in s on its diagonal, modulo the tolerance tol.

    use Multprec_Complex_Matrices;

    ut : Matrix(u'range(2),u'range(1)) := Conjugate_Transpose(u); 
    utx : Matrix(u'range(2),x'range(2)) := ut*x;
    utxv : Matrix(u'range(2),v'range(2)) := utx*v;

  begin
    put_line("The matrix u^H*x*v : "); put(utxv);
    for i in utxv'range(1) loop
      for j in utxv'range(2) loop
        if i = j then
         if AbsVal(utxv(i,j) - s(i)) > tol
          then return false;
         end if;
       else
         if AbsVal(utxv(i,j)) > tol
          then return false;
         end if;
        end if;
      end loop;
    end loop;
    Multprec_Complex_Matrices.Clear(ut);
    Multprec_Complex_Matrices.Clear(utx);
    Multprec_Complex_Matrices.Clear(utxv);
    return true;
  end Is_SVD;

  procedure Test_SVD_Output
              ( x,u,v : in Multprec_Complex_Matrices.Matrix;
                s,e : in Multprec_Complex_Vectors.Vector;
                info : in integer32 ) is

  -- DESCRIPTION :
  --   Tests whether the matrices u and v are orthogonal and
  --   whether u^H*x*v yields a diagonal matrix with the singular
  --   values on the diagonal.  Diagnostics are written to screen.

    tol : constant double_float := 1.0E-8;

  begin
    put("info = " ); put(info,1); new_line;
    put_line("The e values : "); put_line(e);
    put_line("The singular values : "); put_line(s);
    put_line("The matrix u : "); put(u);
    if Is_Orthogonal(u,tol)
     then put_line("The matrix u is orthogonal.");
     else put_line("The matrix u is NOT orthogonal!  BUG!!!");
    end if;
    put_line("The matrix v : "); put(v);
    if Is_Orthogonal(v,tol)
     then put_line("The matrix v is orthogonal.");
     else put_line("The matrix v is NOT orthogonal!  BUG!!!");
    end if;
    if Is_SVD(x,u,v,s,tol)
     then put_line("Tested u^H*x*v = s successfully.");
     else put_line("u^H*x*v /= s, found BUG!!!");
    end if;
  end Test_SVD_Output;

  procedure Test_SVD_Solver
              ( a,u,v : in Multprec_Complex_Matrices.Matrix;
                s,b : in Multprec_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the solver to use SVD of a to solve the system a*x = b,
  --   and performs several diagnostic tests.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;

    res : Vector(b'range);
    x : Vector(a'range(2));
    inv : Matrix(a'range(2),a'range(1)) := Inverse(u,v,s);
    ai : Matrix(a'range(1),a'range(1)) := a*inv;
    ia : Matrix(a'range(2),a'range(1)) := inv*a;
    tol : constant double_float := 1.0E-8;

  begin
    put("The rank of the matrix : ");
    put(Rank(s),1); new_line;
    put("The rank of the matrix with respect to tolerance : ");
    put(Rank(s,tol),1); new_line;
    put("Inverse of condition number : ");
    put(Inverse_Condition_Number(s),3); new_line;
    if Is_Identity(ai,tol) then
      put_line("The product a*inv(a) is the indentity matrix.");
    else
      put_line("The product a*inv(a) is NOT the indentity matrix!");
      if a'length(1) > a'length(2) then
        put_line("... okay, inv(a) can only be a left inverse.");
      else
        put_line("This may be a genuine bug, look at the product:");
        put(ai,3);
      end if;
    end if;
    if Is_Identity(ia,tol) then
      put_line("The product inv(a)*a is the indentity matrix.");
    else
      put_line("The product inv(a)*a is NOT the indentity matrix!");
      if a'length(1) < a'length(2) then
        put_line("... okay, inv(a) can only be a right inverse.");
      else
        put_line("This may be a genuine bug, look at the product:");
        put(ia,3);
      end if;
    end if;
    x := Solve(u,v,s,b);
    res := a*x;
    Min(res);
    Add(res,b);
    Multprec_Complex_Matrices.Clear(inv);
    Multprec_Complex_Matrices.Clear(ia);
    Multprec_Complex_Matrices.Clear(ai);
    put_line("The residual vector : "); put_line(res);
  end Test_SVD_Solver;

  procedure Test_SVD_on_Random_Matrix
              ( n,p : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns with number of the given size.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;

    x : Matrix(1..n,1..p) := Random_Matrix(natural32(n),natural32(p),size);
    y : Matrix(1..n,1..p);
    mm : constant integer32 := Multprec_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Copy(x,y);
    put_line("The random matrix : "); put(x);
    SVD(x,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
  end Test_SVD_on_Random_Matrix;

  procedure Test_SVD_on_Random_System
              ( n,p : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns with number of the given size.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;

    a : Matrix(1..n,1..p) := Random_Matrix(natural32(n),natural32(p),size);
    b : Vector(1..n) := Random_Vector(1,n,size);
    y : Matrix(1..n,1..p);
    mm : constant integer32 := Multprec_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Copy(a,y);
    put_line("The random matrix : "); put(a);
    SVD(a,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
    Test_SVD_Solver(y,u,v,s,b);
    Multprec_Complex_Matrices.Clear(a);
    Multprec_Complex_Vectors.Clear(b);
  end Test_SVD_on_Random_System;

  procedure Test_SVD_on_Given_System
             ( n,p : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Uses SVD to solve a randomly generated linear system
  --   of n equations in p unknowns with number of the given size.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;
 
    a : Matrix(1..n,1..p) := Read_Matrix(n,p);
    b : Vector(1..n) := Read_Vector(n);
    y : Matrix(1..n,1..p);
    mm : constant integer32 := Multprec_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Set_Size(a,size);
    Set_Size(b,size);
    Copy(a,y);
    put_line("The random matrix : "); put(a);
    SVD(a,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
    Test_SVD_Solver(y,u,v,s,b);
  end Test_SVD_on_Given_System;

  procedure Test_SVD_on_Given_Matrix
              ( n,p : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Uses SVD to solve a user given linear system
  --   of n equations in p unknowns with number of the given size.

    use Multprec_Complex_Vectors;
    use Multprec_Complex_Matrices;

    x : Matrix(1..n,1..p) := Read_Matrix(n,p);
    y : Matrix(1..n,1..p);
    mm : constant integer32 := Multprec_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Set_Size(x,size);
    Copy(x,y);
    put_line("The given matrix : "); put(x);
    SVD(x,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
  end Test_SVD_on_Given_Matrix;

  procedure Main is

    n,m : integer32 := 0;
    sz : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing singular value decomposition (SVD).");
    new_line;
    put("Give number of rows : "); get(n);
    put("Give number of columns : "); get(m);
    put("Give size of the numbers (0 for standard precision): ");
    get(sz);
    new_line;
    put_line("MENU for testing the operations :");
    put_line("  1. give your own matrix to compute the SVD;");
    put_line("  2. use SVD to solve a given linear system;");
    put_line("  3. let computer generate a random matrix to test the SVD;");
    put_line("  4. use SVD to solve a random linear system a*x = b.");
    put("Type 1, 2, 3, or 4 to select a testing option : ");
    Ask_Alternative(ans,"1234");
    if sz = 0
     then case ans is
            when '1' => Test_SVD_on_Given_Matrix(n,m);
            when '2' => Test_SVD_on_Given_System(n,m);
            when '3' => Test_SVD_on_Random_Matrix(n,m);
            when '4' => Test_SVD_on_Random_System(n,m);
            when others => null;
          end case;
     else case ans is
            when '1' => Test_SVD_on_Given_Matrix(n,m,sz);
            when '2' => Test_SVD_on_Given_System(n,m,sz);
            when '3' => Test_SVD_on_Random_Matrix(n,m,sz);
            when '4' => Test_SVD_on_Random_System(n,m,sz);
            when others => null;
          end case;
    end if;
  end Main;

begin
  Main;
end ts_svd;
