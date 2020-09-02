with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Floating_Vectors_io;       use Multprec_Floating_Vectors_io;
with Multprec_Floating_Matrices_io;      use Multprec_Floating_Matrices_io;
with Multprec_Floating_Norms_Equals;     use Multprec_Floating_Norms_Equals;
with Multprec_Floating_QR_Least_Squares; use Multprec_Floating_QR_Least_Squares;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Matrices_io;       use Multprec_Complex_Matrices_io;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Random_Matrices;           use Multprec_Random_Matrices;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_QR_Least_Squares;  use Multprec_Complex_QR_Least_Squares;

package body Test_Multprec_QRLS_Solvers is

  function Extract_Upper_Triangular
                ( a : Multprec_Floating_Matrices.Matrix )
                return Multprec_Floating_Matrices.Matrix is

    res : Multprec_Floating_Matrices.Matrix(a'range(1),a'range(2));

  begin
    for i in a'range(1) loop
      for j in a'first(2)..(i-1) loop
        res(i,j) := Create(0.0);
      end loop;
      for j in i..a'last(2) loop
        Copy(a(i,j),res(i,j));
      end loop;
    end loop;
    return res;
  end Extract_Upper_Triangular;

  function Extract_Upper_Triangular
                ( a : Multprec_Complex_Matrices.Matrix )
                return Multprec_Complex_Matrices.Matrix is

    use Multprec_Complex_Numbers;

    res : Multprec_Complex_Matrices.Matrix(a'range(1),a'range(2));

  begin
    for i in a'range(1) loop
      for j in a'first(2)..(i-1) loop
        res(i,j) := Create(integer(0));
      end loop;
      for j in i..a'last(2) loop
        Copy(a(i,j),res(i,j));
      end loop;
    end loop;
    return res;
  end Extract_Upper_Triangular;

  function Differences ( a,b : in Multprec_Floating_Matrices.Matrix )
                       return Floating_Number is

    sum : Floating_Number := Create(0.0);
    dif : Floating_Number;
    absdif : Floating_Number;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        dif := a(i,j) - b(i,j);
        absdif := AbsVal(dif);
        Add(sum,absdif);
        Clear(dif); Clear(absdif);
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in Multprec_Complex_Matrices.Matrix )
                       return Floating_Number is

    use Multprec_Complex_Numbers;

    sum : Floating_Number := Create(0.0);
    dif : Complex_Number;
    absdif : Floating_Number;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        dif := a(i,j) - b(i,j);
        absdif := AbsVal(dif);
        Add(sum,absdif);
        Clear(dif); Clear(absdif);
      end loop;
    end loop;
    return sum;
  end Differences;

  function Orthogonality_Check_Sum
             ( q : Multprec_Complex_Matrices.Matrix )
             return Floating_Number is

    use Multprec_Complex_Numbers;

    sum : Floating_Number := Create(0.0);
    absip : Floating_Number;
    ip,acc : Complex_Number;

  begin
    for j in q'range(2) loop
      for k in j+1..q'last(2) loop
        ip := Create(integer(0));
        for i in q'range(1) loop
          acc := Conjugate(q(i,j));
          Mul(acc,q(i,k));
          Add(ip,acc);
          Clear(acc);
        end loop;
        absip := AbsVal(ip);
        Add(sum,absip);
        Clear(ip);
        Clear(absip);
      end loop;
    end loop;
    return sum;
  end Orthogonality_Check_Sum;

  function Orthogonality_Check_Sum
             ( q : Multprec_Floating_Matrices.Matrix )
             return Floating_Number is

    sum : Floating_Number := Create(0.0);
    absip,ip,acc : Floating_Number;

  begin
    for j in q'range(2) loop
      for k in j+1..q'last(2) loop
        ip := Create(0.0);
        for i in q'range(1) loop
          acc := q(i,j)*q(i,k);
          Add(ip,acc);
          Clear(acc);
        end loop;
        absip := AbsVal(ip);
        Add(sum,absip);
        Clear(ip);
        Clear(absip);
      end loop;
    end loop;
    return sum;
  end Orthogonality_Check_Sum;

  procedure Test_QRD ( a,q,r : in Multprec_Floating_Matrices.Matrix;
                       output : in boolean ) is

    wrk : Multprec_Floating_Matrices.Matrix(a'range(1),a'range(2));
    use Multprec_Floating_Matrices;
    dif : Floating_Number;

  begin
    if output
     then put_line("The upper triangular part R :"); put(r,3);
    end if;
    wrk := q*r;
    if output
     then put_line("q*r :"); put(wrk,3); 
    end if;
    put("Difference in 1-norm between the matrix and q*r : ");
    dif := Differences(a,wrk);
    put(dif,3,3,3); new_line;
    Clear(dif);
    Multprec_Floating_Matrices.Clear(wrk);
    dif := Orthogonality_Check_Sum(q);
    put("Orthogonality check sum : ");
    put(dif,3,3,3); new_line;
    Clear(dif);
  end Test_QRD;

  procedure Test_QRD ( a,q,r : in Multprec_Complex_Matrices.Matrix;
                       output : in boolean ) is

    wrk : Multprec_Complex_Matrices.Matrix(a'range(1),a'range(2));
    use Multprec_Complex_Matrices;
    dif : Floating_Number;

  begin
    if output
     then put_line("The upper triangular part R :"); put(r,3);
    end if;
    wrk := q*r;
    if output
     then put_line("q*r :"); put(wrk,3); 
    end if;
    put("Difference in 1-norm between the matrix and q*r : ");
    dif := Differences(a,wrk);
    put(dif,3,3,3); new_line;
    Clear(dif);
    Multprec_Complex_Matrices.Clear(wrk);
    dif := Orthogonality_Check_Sum(q);
    put("Orthogonality check sum : ");
    put(dif,3,3,3); new_line;
    Clear(dif);
  end Test_QRD;

  procedure Multprec_Real_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Multprec_Floating_Matrices.Matrix;
                b : in Multprec_Floating_Vectors.Vector;
                output : in boolean ) is

    zero : Floating_Number := create(0.0);
    wrk : Multprec_Floating_Matrices.Matrix(1..n,1..m);
    qraux : Multprec_Floating_Vectors.Vector(1..m);
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : Multprec_Floating_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : Multprec_Floating_Vectors.Vector(1..n);
    info : integer32;
    use Multprec_Floating_Matrices;
    use Multprec_Floating_Vectors;

  begin
    Multprec_Floating_Matrices.Copy(a,wrk);
    for i in 1..m loop
      Copy(zero,qraux(i));
    end loop;
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
    Multprec_Floating_Numbers.Clear(zero);
  end Multprec_Real_LS_Test;          

  procedure Multprec_Real_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Multprec_Floating_Matrices.Matrix;
                output : in boolean ) is

    zero : Floating_Number := create(0.0);
    wrk : Multprec_Floating_Matrices.Matrix(1..n,1..m);
    bas : Multprec_Floating_Matrices.Matrix(1..n,1..n);
    qraux : Multprec_Floating_Vectors.Vector(1..m);
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);

  begin
    for i in 1..m loop
      Copy(zero,qraux(i));
    end loop;
    Multprec_Floating_Matrices.Copy(a,wrk);
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
        Copy(wrk(i,j),bas(i,j));
      end loop;
      for j in n+1..m loop
        Copy(zero,bas(i,j));
      end loop;
    end loop;
    Basis(bas,a);
    if output
     then put_line("The orthogonal part Q of QR  :"); put(bas,3);
    end if;
    Test_QRD(a,bas,Extract_Upper_Triangular(wrk),output);
    Multprec_Floating_Numbers.Clear(zero);
  end Multprec_Real_QR_Test;

  procedure Multprec_Interactive_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Multprec_Floating_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    Multprec_Real_QR_Test(n,m,piv,a,true);
  end Multprec_Interactive_Real_QR_Test;

  procedure Multprec_Interactive_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Multprec_Floating_Matrices.Matrix(1..n,1..m);
    b : Multprec_Floating_Vectors.Vector(1..n);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    put("Give right-hand size "); put(n,1);
    put_line("-vector : "); get(b);
    Multprec_Real_LS_Test(n,m,piv,a,b,true);
  end Multprec_Interactive_Real_LS_Test;

  procedure Multprec_Random_Real_QR_Test
              ( n,m : in integer32; sz : in natural32; piv : in boolean ) is

    a : Multprec_Floating_Matrices.Matrix(1..n,1..m);
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
      a := Random_Matrix(natural32(n),natural32(m),sz);
      Multprec_Real_QR_Test(n,m,piv,a,output);
      Multprec_Floating_Matrices.Clear(a);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factoriziations on multprec random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random Multprec Real QR Factorizations");
  end Multprec_Random_Real_QR_Test;

  procedure Multprec_Random_Real_LS_Test
              ( n,m : in integer32; sz : in natural32; piv : in boolean ) is

    a : Multprec_Floating_Matrices.Matrix(1..n,1..m);
    b : Multprec_Floating_Vectors.Vector(1..n);
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
      a := Random_Matrix(natural32(n),natural32(m),sz);
      b := Random_Vector(1,n,sz);
      Multprec_Real_LS_Test(n,m,piv,a,b,output);
      Multprec_Floating_Matrices.Clear(a);
      Multprec_Floating_Vectors.Clear(b);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" real least squares on quaddobl random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Testing Multprec Real Least Squares");
  end Multprec_Random_Real_LS_Test;

  procedure Multprec_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : Multprec_Complex_Matrices.Matrix;
                output : in boolean ) is

    use Multprec_Complex_Numbers;
    wrk,upp : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    bas : Multprec_Complex_Matrices.Matrix(1..n,1..n);
    qraux : Multprec_Complex_Vectors.Vector(1..m)
          := (1..m => Create(integer(0)));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);

  begin
    if output
     then put_line("The matrix : "); put(a,3);
    end if;
    Multprec_Complex_Matrices.Copy(a,wrk);
    QRD(wrk,qraux,jpvt,piv);
    if output then
      put_line("The matrix after QR : "); put(wrk,3);
      put_line("The vector qraux : "); put(qraux,3); new_line;
    end if;
   -- put("The vector jpvt : "); put(jpvt); new_line;
    if not piv then
      for i in wrk'range(1) loop
        for j in wrk'range(2) loop
          Copy(wrk(i,j),bas(i,j));
        end loop;
        for j in n+1..m loop
          bas(i,j) := Create(integer(0));
        end loop;
      end loop;
      Basis(bas,a);
      if output
       then put_line("The orthogonal part Q of QR  :"); put(bas,3);
      end if;
      upp := Extract_Upper_Triangular(wrk);
      Test_QRD(a,bas,upp,output);
      Multprec_Complex_Matrices.Clear(bas); 
      Multprec_Complex_Matrices.Clear(upp);
    end if;
    Multprec_Complex_Matrices.Clear(wrk);
    Multprec_Complex_Vectors.Clear(qraux);
  end Multprec_Complex_QR_Test;

  procedure Multprec_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : Multprec_Complex_Matrices.Matrix;
                b : Multprec_Complex_Vectors.Vector;
                output : in boolean ) is

    use Multprec_Complex_Numbers;
    wrk : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    qraux : Multprec_Complex_Vectors.Vector(1..m)
          := (1..m => Create(integer(0)));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : Multprec_Complex_Vectors.Vector(1..m);
    rsd,dum1,dum2,dum3,res,eva : Multprec_Complex_Vectors.Vector(1..n);
    resi : Floating_Number;
    info : integer32;
    use Multprec_Complex_Matrices;
    use Multprec_Complex_Vectors; 

  begin
    if output
     then put_line("The matrix : "); put(a,3);
    end if;
    Multprec_Complex_Matrices.Copy(a,wrk);
    QRD(wrk,qraux,jpvt,piv);
    if output then
      put_line("The matrix after QR : "); put(wrk,3);
      put_line("The vector qraux : "); put(qraux,3); new_line;
    end if;
    if piv
     then put("The vector jpvt : "); put(jpvt); new_line;
    end if;
    QRLS(wrk,n,n,m,qraux,b,dum1,dum2,sol,rsd,dum3,110,info);
    Multprec_Complex_Vectors.Clear(dum1);
    Multprec_Complex_Vectors.Clear(dum2);
    Multprec_Complex_Vectors.Clear(dum3);
    if output
     then put_line("The solution : "); put(sol,3); new_line;
    end if;
    eva := a*sol;
    res := b - eva;
    if output then
      put_line("right-hand size - matrix*solution : ");
      put(res,3); new_line;
    end if;
    resi := Sum_Norm(res);
    put("Sum norm of residual : "); put(resi,3); new_line;
    Clear(resi);
    Multprec_Complex_Vectors.Clear(eva);
    Multprec_Complex_Vectors.Clear(res);
    Multprec_Complex_Vectors.Clear(sol);
    Multprec_Complex_Vectors.Clear(rsd);
    Multprec_Complex_Vectors.Clear(qraux);
    Multprec_Complex_Matrices.Clear(wrk);
  end Multprec_Complex_LS_Test;

  procedure Multprec_Interactive_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      Multprec_Complex_QR_Test(n,m,piv,a,true);
      Multprec_Complex_Matrices.Clear(a);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Multprec_Interactive_Complex_QR_Test;

  procedure Multprec_Interactive_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    b : Multprec_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      put("Give right-hand size "); put(n,1);
      put_line("-vector : "); get(b); 
      Multprec_Complex_LS_Test(n,m,piv,a,b,true);
      Multprec_Complex_Matrices.Clear(a);
      Multprec_Complex_Vectors.Clear(b);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Multprec_Interactive_Complex_LS_Test;

  procedure Multprec_Random_Complex_QR_Test
              ( n,m : in integer32; sz : in natural32; piv : in boolean ) is

    a : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    nb : natural32 := 0;
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
      a := Random_Matrix(natural32(n),natural32(m),sz);
      Multprec_Complex_QR_Test(n,m,piv,a,output);
      Multprec_Complex_Matrices.Clear(a);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factorizations on multprec random complex matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random Multprec Complex QR Testing");
  end Multprec_Random_Complex_QR_Test;

  procedure Multprec_Random_Complex_LS_Test
              ( n,m : integer32; sz : in natural32; piv : in boolean ) is

    a : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    b : Multprec_Complex_Vectors.Vector(1..n);
    nb : natural32 := 0;
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
      a := Random_Matrix(natural32(n),natural32(m),sz);
      b := Random_Vector(1,n,sz);
      Multprec_Complex_LS_Test(n,m,piv,a,b,output);
      Multprec_Complex_Matrices.Clear(a);
      Multprec_Complex_Vectors.Clear(b);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" least squares on multprec random complex matrices");
    new_line;
    print_times(Standard_Output,timer,"Random Multprec Complex Least Squares");
  end Multprec_Random_Complex_LS_Test;

  procedure Main is

    n,m : integer32 := 0;
    sz : natural32 := 0;
    choice : character;
    piv : constant boolean := false;

  begin
    loop
      new_line;
      put_line
        ("MENU to test QR and Least Squares in arbitrary multiprecision : ");
      put_line("  0. Exit this program.");
      put_line("  1. QR-decomposition on given multprec floating matrix.");
      put_line("  2.                  on random multprec floating matrix.");
      put_line("  3. Least Squares on given multprec floating matrix.");
      put_line("  4.               on random multprec floating matrix.");
      put_line("  5. QR-decomposition on given multprec complex matrix.");
      put_line("  6.                  on random multprec complex matrix.");
      put_line("  7. Least Squares on given multprec complex matrix.");
      put_line("  8.               on random multprec complex matrix.");
      put("Type 1, 2, 3, 4, 5, 6, 7, 8 for a test, or 0 to exit : ");
      Ask_Alternative(choice,"012345678");
      exit when (choice = '0');
      new_line;
      put("Give the number of rows of the matrix : "); get(n);
      put("Give the number of columns of the matrix : "); get(m);
      put("Give the size of the numbers : "); get(sz);
      case choice is
        when '1' => Multprec_Interactive_Real_QR_Test(n,m,piv);
        when '2' => Multprec_Random_Real_QR_Test(n,m,sz,piv);
        when '3' => Multprec_Interactive_Real_LS_Test(n,m,piv);
        when '4' => Multprec_Random_Real_LS_Test(n,m,sz,piv);
        when '5' => Multprec_Interactive_Complex_QR_Test(n,m,piv);
        when '6' => Multprec_Random_Complex_QR_Test(n,m,sz,piv);
        when '7' => Multprec_Interactive_Complex_LS_Test(n,m,piv);
        when '8' => Multprec_Random_Complex_LS_Test(n,m,sz,piv);
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Multprec_QRLS_Solvers;
