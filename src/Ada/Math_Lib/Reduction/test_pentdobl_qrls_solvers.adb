with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Penta_Double_Numbers_io;           use Penta_Double_Numbers_io;
with Penta_Double_Vectors_io;           use Penta_Double_Vectors_io;
with Penta_Double_Matrices_io;          use Penta_Double_Matrices_io;
with PentDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with PentDobl_Complex_Vectors_io;        use PentDobl_Complex_Vectors_io;
with PentDobl_Complex_Matrices_io;       use PentDobl_Complex_Matrices_io;
with PentDobl_Random_Vectors;            use PentDobl_Random_Vectors;
with PentDobl_Random_Matrices;           use PentDobl_Random_Matrices;
with Penta_Double_Vector_Norms;         use Penta_Double_Vector_Norms;
with PentDobl_Complex_Vector_Norms;      use PentDobl_Complex_Vector_Norms;
with Penta_Double_QR_Least_Squares;     use Penta_Double_QR_Least_Squares;
with PentDobl_Complex_QR_Least_Squares;  use PentDobl_Complex_QR_Least_Squares;

package body Test_PentDobl_QRLS_Solvers is

  function Extract_Upper_Triangular
                ( a : Penta_Double_Matrices.Matrix )
                return Penta_Double_Matrices.Matrix is

    res : Penta_Double_Matrices.Matrix(a'range(1),a'range(2));

    zero : constant penta_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'first(2)..(i-1) loop
        res(i,j) := zero;
      end loop;
      for j in i..a'last(2) loop
        res(i,j) := a(i,j);
      end loop;
    end loop;
    return res;
  end Extract_Upper_Triangular;

  function Extract_Upper_Triangular
                ( a : PentDobl_Complex_Matrices.Matrix )
                return PentDobl_Complex_Matrices.Matrix is

    use PentDobl_Complex_Numbers;

    res : PentDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));
    zero : constant penta_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'first(2)..(i-1) loop
        res(i,j) := Create(zero);
      end loop;
      for j in i..a'last(2) loop
        res(i,j) := a(i,j);
      end loop;
    end loop;
    return res;
  end Extract_Upper_Triangular;

  function Differences ( a,b : in Penta_Double_Matrices.Matrix )
                       return penta_double is

    sum : penta_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + abs(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in PentDobl_Complex_Matrices.Matrix )
                       return penta_double is

    use PentDobl_Complex_Numbers;

    sum : penta_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + AbsVal(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Orthogonality_Check_Sum
             ( q : Penta_Double_Matrices.Matrix ) return penta_double is

    sum,ip : penta_double;
    zero : constant penta_double := create(0.0);

  begin
    sum := zero;
    for j in q'range(2) loop
      for k in j+1..q'last(2) loop
        ip := zero;
        for i in q'range(1) loop
          ip := ip + q(i,j)*q(i,k);
        end loop;
        sum := sum + abs(ip);
      end loop;
    end loop;
    return sum;
  end Orthogonality_Check_Sum;

  function Orthogonality_Check_Sum 
             ( q : PentDobl_Complex_Matrices.Matrix )
             return penta_double is

    use PentDobl_Complex_Numbers;

    zero : constant penta_double := create(0.0);
    sum : penta_double := zero;
    ip : Complex_Number;

  begin
    for j in q'range(2) loop
      for k in j+1..q'last(2) loop
        ip := Create(zero);
        for i in q'range(1) loop
          ip := ip + Conjugate(q(i,j))*q(i,k);
        end loop;
        sum := sum + AbsVal(ip);
      end loop;
    end loop;
    return sum;
  end Orthogonality_Check_Sum;

  procedure Test_QRD ( a,q,r : in Penta_Double_Matrices.Matrix;
                       output : in boolean ) is

    wrk : Penta_Double_Matrices.Matrix(a'range(1),a'range(2));

    use Penta_Double_Matrices;

  begin
    if output
     then put_line("The upper triangular part R :"); put(r,3);
    end if;
    wrk := q*r;
    if output
     then put_line("q*r :"); put(wrk,3); 
    end if;
    put("Difference in 1-norm between the matrix and q*r : ");
    put(Differences(a,wrk),3); new_line;
    put("Orthogonality check sum : ");
    put(Orthogonality_Check_Sum(q),3); new_line;
  end Test_QRD;

  procedure Test_QRD ( a,q,r : in PentDobl_Complex_Matrices.Matrix;
                       output : in boolean ) is

    wrk : PentDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));

    use PentDobl_Complex_Matrices;

  begin
    if output
     then put_line("The upper triangular part R :"); put(r,3);
    end if;
    wrk := q*r;
    if output
     then put_line("q*r :"); put(wrk,3); 
    end if;
    put("Difference in 1-norm between the matrix and q*r : ");
    put(Differences(a,wrk),3); new_line;
    put("Orthogonality check sum : ");
    put(Orthogonality_Check_Sum(q),3); new_line;
  end Test_QRD;

  procedure PentDobl_Real_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Penta_Double_Matrices.Matrix;
                b : in Penta_Double_Vectors.Vector;
                output : in boolean ) is

    zero : constant penta_double := create(0.0);
    wrk : Penta_Double_Matrices.Matrix(1..n,1..m) := a;
    qraux : Penta_Double_Vectors.Vector(1..m) := (1..m => zero);
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : Penta_Double_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : Penta_Double_Vectors.Vector(1..n);
    info : integer32;
    use Penta_Double_Matrices;
    use Penta_Double_Vectors;

  begin
    if output
     then put_line("The matrix : "); put(a,3);
    end if;
    QRD(wrk,qraux,jpvt,piv);
    if output then
      put_line("The matrix after QR : "); put(wrk,3);
      put_line("The vector qraux : "); put(qraux,3); new_line;
    end if;
    if piv then
      put("The vector jpvt : "); put(jpvt); new_line;
      Permute_Columns(wrk,jpvt);
    end if;
    QRLS(wrk,n,n,m,qraux,b,dum2,dum3,sol,rsd,dum,110,info);
    if piv
     then Permute(sol,jpvt);
    end if;
    if output
     then put_line("The solution : "); put(sol,3); new_line;
    end if;
    dum := b - a*sol;
    if output then
      put_line("right-hand size - matrix*solution : "); 
      put(dum,3); new_line;
    end if;
    put("The norm of residual : "); put(Sum_Norm(dum),3); new_line;
  end PentDobl_Real_LS_Test;          

  procedure PentDobl_Real_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Penta_Double_Matrices.Matrix;
                output : in boolean ) is

    zero : constant penta_double := create(0.0);
    wrk : Penta_Double_Matrices.Matrix(1..n,1..m) := a;
    bas : Penta_Double_Matrices.Matrix(1..n,1..n);
    qraux : Penta_Double_Vectors.Vector(1..m) := (1..m => zero);
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);

  begin
    if output
     then put_line("The matrix : "); put(a,3);
    end if;
    QRD(wrk,qraux,jpvt,piv);
    if output then
      put_line("The matrix after QR : "); put(wrk,3);
      put_line("The vector qraux : "); put(qraux,3); new_line;
    end if;
    if piv then
      put("The vector jpvt : "); put(jpvt); new_line;
      Permute_Columns(wrk,jpvt);
    end if;
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        bas(i,j) := wrk(i,j);
      end loop;
      for j in n+1..m loop
        bas(i,j) := zero;
      end loop;
    end loop;
    Basis(bas,a);
    if output
     then put_line("The orthogonal part Q of QR  :"); put(bas,3);
    end if;
    Test_QRD(a,bas,Extract_Upper_Triangular(wrk),output);
  end PentDobl_Real_QR_Test;

  procedure PentDobl_Interactive_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Penta_Double_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    PentDobl_Real_QR_Test(n,m,piv,a,true);
  end PentDobl_Interactive_Real_QR_Test;

  procedure PentDobl_Interactive_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Penta_Double_Matrices.Matrix(1..n,1..m);
    b : Penta_Double_Vectors.Vector(1..n);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    put("Give right-hand size "); put(n,1);
    put_line("-vector : "); get(b);
    PentDobl_Real_LS_Test(n,m,piv,a,b,true);
  end PentDobl_Interactive_Real_LS_Test;

  procedure PentDobl_Random_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Penta_Double_Matrices.Matrix(1..n,1..m);
    nb : integer32 := 0;
    output : boolean;
    ans : character;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Do you want to see all matrices and vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    tstart(timer);
    for i in 1..nb loop
      a := Random_Matrix(natural32(n),natural32(m));
      PentDobl_Real_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factoriziations on PentDobl random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random PentDobl Real QR Factorizations");
  end PentDobl_Random_Real_QR_Test;

  procedure PentDobl_Random_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Penta_Double_Matrices.Matrix(1..n,1..m);
    b : Penta_Double_Vectors.Vector(1..n);
    nb : integer32 := 0;
    ans : character;
    output : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Do you want to see all matrices and vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y'); 
    tstart(timer);
    for i in 1..nb loop
      a := Random_Matrix(natural32(n),natural32(m));
      b := Random_Vector(1,n);
      PentDobl_Real_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" real least squares on PentDobl random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Testing PentDobl Real Least Squares");
  end PentDobl_Random_Real_LS_Test;

  procedure PentDobl_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : PentDobl_Complex_Matrices.Matrix;
                output : in boolean ) is

    use PentDobl_Complex_Numbers;
    wrk : PentDobl_Complex_Matrices.Matrix(1..n,1..m) := a;
    bas : PentDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant penta_double := create(0.0);
    qraux : PentDobl_Complex_Vectors.Vector(1..m) := (1..m => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);

  begin
    if output
     then put_line("The matrix : "); put(a,3);
    end if;
    QRD(wrk,qraux,jpvt,piv);
    if output
     then put_line("The matrix after QR : "); put(wrk,3);
          put_line("The vector qraux : "); put(qraux,3); new_line;
    end if;
   -- put("The vector jpvt : "); put(jpvt); new_line;
    if not piv then
      for i in wrk'range(1) loop
        for j in wrk'range(2) loop
          bas(i,j) := wrk(i,j);
        end loop;
        for j in n+1..m loop
          bas(i,j) := Create(zero);
        end loop;
      end loop;
      Basis(bas,a);
      if output
       then put_line("The orthogonal part Q of QR  :"); put(bas,3);
      end if;
      Test_QRD(a,bas,Extract_Upper_Triangular(wrk),output);
    end if;
  end PentDobl_Complex_QR_Test;

  procedure PentDobl_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : PentDobl_Complex_Matrices.Matrix;
                b : PentDobl_Complex_Vectors.Vector;
                output : in boolean ) is

    use PentDobl_Complex_Numbers;
    wrk : PentDobl_Complex_Matrices.Matrix(1..n,1..m) := a;
    zero : constant penta_double := create(0.0);
    qraux : PentDobl_Complex_Vectors.Vector(1..m) := (1..m => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : PentDobl_Complex_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : PentDobl_Complex_Vectors.Vector(1..n);
    info : integer32;
    use PentDobl_Complex_Matrices;
    use PentDobl_Complex_Vectors; 

  begin
    if output
     then put_line("The matrix : "); put(a,3);
    end if;
    QRD(wrk,qraux,jpvt,piv);
    if output then
      put_line("The matrix after QR : "); put(wrk,3);
      put_line("The vector qraux : "); put(qraux,3); new_line;
    end if;
   -- put("The vector jpvt : "); put(jpvt); new_line;
    QRLS(wrk,n,m,qraux,b,dum,dum2,sol,rsd,dum3,110,info);
    if output
     then put_line("The solution : "); put(sol,3); new_line;
    end if;
    dum := b - a*sol;
    if output then 
      put_line("right-hand size - matrix*solution : ");
      put(dum,3); new_line;
    end if;
    put("Sum norm of residual : "); put(Sum_Norm(dum),3); new_line;
  end PentDobl_Complex_LS_Test;

  procedure PentDobl_Interactive_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : PentDobl_Complex_Matrices.Matrix(1..n,1..m);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      PentDobl_Complex_QR_Test(n,m,piv,a,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end PentDobl_Interactive_Complex_QR_Test;

  procedure PentDobl_Interactive_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : PentDobl_Complex_Matrices.Matrix(1..n,1..m);
    b : PentDobl_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      put("Give right-hand size "); put(n,1);
      put_line("-vector : "); get(b); 
      PentDobl_Complex_LS_Test(n,m,piv,a,b,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end PentDobl_Interactive_Complex_LS_Test;

  procedure PentDobl_Random_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : PentDobl_Complex_Matrices.Matrix(1..n,1..m);
    nb : integer32 := 0;
    ans : character;
    output : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Do you want to see all matrices and vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    tstart(timer);
    for i in 1..nb loop
      a := Random_Matrix(natural32(n),natural32(m));
      PentDobl_Complex_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factorizations on random penta double complex matrices.");
    new_line;
    print_times(Standard_Output,timer,
                "Random PentDobl Complex QR Factorizations");
  end PentDobl_Random_Complex_QR_Test;

  procedure PentDobl_Random_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : PentDobl_Complex_Matrices.Matrix(1..n,1..m);
    b : PentDobl_Complex_Vectors.Vector(1..n);
    nb : integer32 := 0;
    ans : character;
    output : boolean;
    timer : Timing_Widget;

  begin
    put("Give the number of tests : "); get(nb);
    put("Do you want to see all matrices and vectors ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    tstart(timer);
    for i in 1..nb loop
      a := Random_Matrix(natural32(n),natural32(m));
      b := Random_Vector(1,n);
      PentDobl_Complex_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" least squares on random penta double complex matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random PentDobl Complex Least Squares");
  end PentDobl_Random_Complex_LS_Test;

  procedure Main is

    n,m : integer32 := 0;
    choice : character;
    piv : constant boolean := false;

  begin
    loop
      new_line;
      put_line
        ("MENU to test QR and Least Squares in penta double precision : ");
      put_line("  0. Exit this program.");
      put_line("  1. QR-decomposition on given penta double real matrix.");
      put_line("  2.                                        complex matrix.");
      put_line("  3.                  on random penta double real matrix.");
      put_line("  4.                                         complex matrix.");
      put_line("  5. Least Squares on given penta double real matrix.");
      put_line("  6.                                     complex matrix.");
      put_line("  7.               on random penta double real matrix.");
      put_line("  8.                                      complex matrix.");
      put("Type 1, 2, 3, 4, 5, 6, 7, 8 for a test, or 0 to exit : ");
      Ask_Alternative(choice,"012345678");
      exit when (choice = '0');
      new_line;
      put("Give the number of rows of the matrix : "); get(n);
      put("Give the number of columns of the matrix : "); get(m);
      case choice is
        when '1' => PentDobl_Interactive_Real_QR_Test(n,m,piv);
        when '2' => PentDobl_Interactive_Complex_QR_Test(n,m,piv);
        when '3' => PentDobl_Random_Real_QR_Test(n,m,piv);
        when '4' => PentDobl_Random_Complex_QR_Test(n,m,piv);
        when '5' => PentDobl_Interactive_Real_LS_Test(n,m,piv);
        when '6' => PentDobl_Interactive_Complex_LS_Test(n,m,piv);
        when '7' => PentDobl_Random_Real_LS_Test(n,m,piv);
        when '8' => PentDobl_Random_Complex_LS_Test(n,m,piv);
        when others => null;
      end case;
    end loop;
  end Main;

end Test_PentDobl_QRLS_Solvers;
