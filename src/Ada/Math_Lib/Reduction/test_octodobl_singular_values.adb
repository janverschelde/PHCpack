with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Octo_Double_Numbers_io;            use Octo_Double_Numbers_io;
with OctoDobl_Complex_Numbers;          use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Vectors_io;       use OctoDobl_Complex_Vectors_io;
with OctoDobl_Complex_Matrices_io;      use OctoDobl_Complex_Matrices_io;
with OctoDobl_Random_Vectors;           use OctoDobl_Random_Vectors;
with OctoDobl_Random_Matrices;          use OctoDobl_Random_Matrices;
with OctoDobl_Complex_Singular_Values;  use OctoDobl_Complex_Singular_Values;

package body Test_OctoDobl_Singular_Values is

  function Read_Vector
             ( n : in integer32 ) return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(1..n);

  begin
    put("Give "); put(n,1);
    put(" complex numbers for a vector of dimension ");
    put(n,1); put_line(" :");
    get(res);
    return res;
  end Read_Vector;

  function Read_Matrix
             ( n,m : in integer32 ) return OctoDobl_Complex_Matrices.Matrix is

    res : OctoDobl_Complex_Matrices.Matrix(1..n,1..m);

  begin
    put("Give "); put(n*m,1); put(" complex numbers for ");
    put(n,1); put("x"); put(m,1); put_line(" matrix :");
    get(res);
    return res;
  end Read_Matrix;

  function Is_Identity ( a : OctoDobl_Complex_Matrices.Matrix;
                         tol : double_float ) return boolean is

    one : constant octo_double := create(1.0);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        if i = j then
          if Absval(a(i,j) - one) > tol
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

  function Is_Orthogonal ( a : OctoDobl_Complex_Matrices.Matrix;
                           tol : double_float ) return boolean is

    use OctoDobl_Complex_Matrices;

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

  function Is_SVD ( x,u,v : OctoDobl_Complex_Matrices.Matrix;
                    s : OctoDobl_Complex_Vectors.Vector;
                    tol : in double_float ) return boolean is

    use OctoDobl_Complex_Matrices;

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
              ( x,u,v : in OctoDobl_Complex_Matrices.Matrix;
                s,e : in OctoDobl_Complex_Vectors.Vector;
                info : in integer32; output : in boolean := true ) is

    ortho,decomp : boolean;
    tol : constant double_float := 1.0E-8;

  begin
    put("info = "); put(info,1); new_line;
    if output
     then put_line("The e values : "); put_line(e);
    end if;
    put_line("The singular values : "); put_line(s);
    if output
     then put_line("The matrix u : "); put(u);
    end if;
    ortho := Is_Orthogonal(u,tol);
    if ortho
     then put_line("The matrix u is orthogonal.");
     else put_line("The matrix u is NOT orthogonal!  BUG!!!");
    end if;
    if output
     then put_line("The matrix v : "); put(v);
    end if;
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
               ( a,u,v : in OctoDobl_Complex_Matrices.Matrix;
                 s,b : in OctoDobl_Complex_Vectors.Vector ) is

    use OctoDobl_Complex_Vectors;
    use OctoDobl_Complex_Matrices;

    res : Vector(b'range) := b;
    x : Vector(a'range(2));
    inv : constant Matrix(a'range(2),a'range(1)) := Inverse(u,v,s);
    ai : constant Matrix(a'range(1),a'range(1)) := a*inv;
    ia : constant Matrix(a'range(2),a'range(2)) := inv*a;
    tol : constant double_float := 1.0E-8;

  begin
    put("The rank of the matrix : "); put(Rank(s),1); new_line;
    put("Inverse of condition number : ");
    put(Inverse_Condition_Number(s),3); new_line;
    if Is_Identity(ai,tol) then
      put_line("The product a*inv(a) is the identity matrix.");
    else
      put_line("The product a*inv(a) is NOT the identity matrix!");
      if a'length(1) > a'length(2) then
        put_line("... okay, inv(a) can only be a left inverse.");
      else
        put_line("This may be a genuine bug, look at the product:");
        put(ai,3);
      end if;
    end if;
    if Is_Identity(ia,tol) then
      put_line("The product inv(a)*a is the identity matrix.");
    else
      put_line("The product inv(a)*a is NOT the identity matrix!");
      if a'length(1) < a'length(2) then
        put_line("... okay, inv(a) can only be a right inverse.");
      else
        put_line("This may be a genuine bug, look at the product:");
        put(ia,3);
      end if;
    end if;
    x := Solve(u,v,s,b);
    res := b - a*x;
    put_line("The residual vector : "); put_line(res);
  end Test_SVD_Solver;

  procedure Test_SVD_on_Given_Matrix ( n,p : in integer32 ) is

    use OctoDobl_Complex_Vectors;
    use OctoDobl_Complex_Matrices;

    x : Matrix(1..n,1..p) := Read_Matrix(n,p);
    y : constant Matrix(1..n,1..p) := x;
    mm : constant integer32 := OctoDobl_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    SVD(x,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
  end Test_SVD_on_Given_Matrix;

  procedure Test_SVD_on_Given_System ( n,p : in integer32 ) is

    use OctoDobl_Complex_Vectors;
    use OctoDobl_Complex_Matrices;

    a : Matrix(1..n,1..p) := Read_Matrix(n,p);
    b : constant Vector(1..n) := Read_Vector(n);
    y : constant Matrix(1..n,1..p) := a;
    mm : constant integer32 := OctoDobl_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    SVD(a,n,p,s,e,u,v,job,info);
    Test_SVD_Output(y,u,v,s,e,info);
    Test_SVD_Solver(y,u,v,s,b);
  end Test_SVD_on_Given_System;

  procedure Test_SVD_on_Random_Matrix ( n,p : in integer32 ) is

    use OctoDobl_Complex_Vectors;
    use OctoDobl_Complex_Matrices;
  
    x,y : Matrix(1..n,1..p);
    mm : constant integer32 := OctoDobl_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;
    ans : character;
    otp : boolean;

  begin
    new_line;
    put("See all vectors and matrices ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    loop
      x := Random_Matrix(natural32(n),natural32(p));
      y := x;
      SVD(x,n,p,s,e,u,v,job,info);
      Test_SVD_Output(y,u,v,s,e,info,otp);
      new_line;
      put("Test another random problem ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_SVD_on_Random_Matrix;

  procedure Test_SVD_on_Random_System ( n,p : in integer32 ) is

    use OctoDobl_Complex_Vectors;
    use OctoDobl_Complex_Matrices;

    a,y : Matrix(1..n,1..p);
    b : Vector(1..n);
    mm : constant integer32 := OctoDobl_Complex_Singular_Values.Min0(n+1,p);
    s : Vector(1..mm);
    e : Vector(1..p);
    u : Matrix(1..n,1..n);
    v : Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;
    ans : character;
    otp : boolean;

  begin
    new_line;
    put("See all vectors and matrices ? (y/n) ");
    Ask_Yes_or_No(ans); otp := (ans = 'y');
    loop
      a := Random_Matrix(natural32(n),natural32(p));
      y := a;
      b := Random_Vector(1,n);
      SVD(a,n,p,s,e,u,v,job,info);
      Test_SVD_Output(y,u,v,s,e,info,otp);
      Test_SVD_Solver(y,u,v,s,b);
      new_line;
      put("Test another random problem ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_SVD_on_Random_System;

  procedure Main is

    n,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing the SVD in octo double precision ...");
    new_line;
    put("Give the number of rows : "); get(n);
    put("Give the number of columns : "); get(m);
    new_line;
    put_line("MENU to test the singular value decomposition :");
    put_line("  1. give your own matrix to compute the SVD;");
    put_line("  2. use SVD to solve a given linear system;");
    put_line("  3. let computer generate a random matrix to test the SVD;");
    put_line("  4. use SVD to solve a random linear system a*x = b.");
    put("Type 1, 2, 3, or 4 to select a testing option : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Test_SVD_on_Given_Matrix(n,m);
      when '2' => Test_SVD_on_Given_System(n,m);
      when '3' => Test_SVD_on_Random_Matrix(n,m);
      when '4' => Test_SVD_on_Random_System(n,m);
      when others => null;
    end case;
  end Main;

end Test_OctoDobl_Singular_Values;
