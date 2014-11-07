with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Random_Matrices;           use Standard_Random_Matrices;
with Standard_Floating_Norms_Equals;     use Standard_Floating_Norms_Equals;
with Standard_Floating_QR_Least_Squares; use Standard_Floating_QR_Least_Squares;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_QR_Least_Squares;  use Standard_Complex_QR_Least_Squares;
with Double_Double_Vectors;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;
with Double_Double_Matrices;
with Double_Double_Matrices_io;          use Double_Double_Matrices_io;
with Double_Double_Vector_Norms;         use Double_Double_Vector_Norms;
with Double_Double_QR_Least_Squares;     use Double_Double_QR_Least_Squares;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;
with DoblDobl_Random_Matrices;           use DoblDobl_Random_Matrices;
with DoblDobl_Complex_QR_Least_Squares;  use DoblDobl_Complex_QR_Least_Squares;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with Quad_Double_Matrices;
with Quad_Double_Matrices_io;            use Quad_Double_Matrices_io;
with Quad_Double_Vector_Norms;           use Quad_Double_Vector_Norms;
with Quad_Double_QR_Least_Squares;       use Quad_Double_QR_Least_Squares;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with QuadDobl_Random_Matrices;           use QuadDobl_Random_Matrices;
with QuadDobl_Complex_QR_Least_Squares;  use QuadDobl_Complex_QR_Least_Squares;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;       use Multprec_Floating_Vectors_io;
with Multprec_Floating_Matrices;
with Multprec_Floating_Matrices_io;      use Multprec_Floating_Matrices_io;
with Multprec_Floating_Norms_Equals;     use Multprec_Floating_Norms_Equals;
with Multprec_Floating_QR_Least_Squares; use Multprec_Floating_QR_Least_Squares;
with Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Matrices;
with Multprec_Complex_Matrices_io;       use Multprec_Complex_Matrices_io;
with Multprec_Random_Vectors;            use Multprec_Random_Vectors;
with Multprec_Random_Matrices;           use Multprec_Random_Matrices;
with Multprec_Complex_Norms_Equals;      use Multprec_Complex_Norms_Equals;
with Multprec_Complex_QR_Least_Squares;  use Multprec_Complex_QR_Least_Squares;

procedure ts_qrd is

-- DESCRIPTION :
--   This program tests the implementation of the QR decomposition
--   and least squares approximation.

-- GENERAL TESTS ON QR DECOMPOSITION :

  function Extract_Upper_Triangular
                ( a : Standard_Floating_Matrices.Matrix )
                return Standard_Floating_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

    res : Standard_Floating_Matrices.Matrix(a'range(1),a'range(2));

  begin
    for i in a'range(1) loop
      for j in a'first(2)..(i-1) loop
        res(i,j) := 0.0;
      end loop;
      for j in i..a'last(2) loop
        res(i,j) := a(i,j);
      end loop;
    end loop;
    return res;
  end Extract_Upper_Triangular;

  function Extract_Upper_Triangular
                ( a : Double_Double_Matrices.Matrix )
                return Double_Double_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

    res : Double_Double_Matrices.Matrix(a'range(1),a'range(2));
    zero : constant double_double := create(0.0);

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
                ( a : Quad_Double_Matrices.Matrix )
                return Quad_Double_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

    res : Quad_Double_Matrices.Matrix(a'range(1),a'range(2));
    zero : constant quad_double := create(0.0);

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
                ( a : Standard_Complex_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

    use Standard_Complex_Numbers;
    res : Standard_Complex_Matrices.Matrix(a'range(1),a'range(2));

  begin
    for i in a'range(1) loop
      for j in a'first(2)..(i-1) loop
        res(i,j) := Create(0.0);
      end loop;
      for j in i..a'last(2) loop
        res(i,j) := a(i,j);
      end loop;
    end loop;
    return res;
  end Extract_Upper_Triangular;

  function Extract_Upper_Triangular
                ( a : DoblDobl_Complex_Matrices.Matrix )
                return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

    use DoblDobl_Complex_Numbers;
    res : DoblDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));
    zero : constant double_double := create(0.0);

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

  function Extract_Upper_Triangular
                ( a : QuadDobl_Complex_Matrices.Matrix )
                return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

    use QuadDobl_Complex_Numbers;
    res : QuadDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));
    zero : constant quad_double := create(0.0);

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

  function Extract_Upper_Triangular
                ( a : Multprec_Floating_Matrices.Matrix )
                return Multprec_Floating_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

    use Multprec_Floating_Numbers;
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

  -- DESCRIPTION :
  --   Returns the upper triangular part of the matrix a.

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

  function Differences ( a,b : in Standard_Floating_Matrices.Matrix )
                       return double_float is

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

    sum : double_float := 0.0;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + abs(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in Double_Double_Matrices.Matrix )
                       return double_double is

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

    sum : double_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + abs(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in Quad_Double_Matrices.Matrix )
                       return quad_double is

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

    sum : quad_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + abs(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in Standard_Complex_Matrices.Matrix )
                       return double_float is

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

    use Standard_Complex_Numbers;
    sum : double_float := 0.0;

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + AbsVal(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in DoblDobl_Complex_Matrices.Matrix )
                       return double_double is

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

    use DoblDobl_Complex_Numbers;
    sum : double_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + AbsVal(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in QuadDobl_Complex_Matrices.Matrix )
                       return quad_double is

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

    use QuadDobl_Complex_Numbers;
    sum : quad_double := create(0.0);

  begin
    for i in a'range(1) loop
      for j in a'range(2) loop
        sum := sum + AbsVal(a(i,j)-b(i,j));
      end loop;
    end loop;
    return sum;
  end Differences;

  function Differences ( a,b : in Multprec_Floating_Matrices.Matrix )
                       return Floating_Number is

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

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

  -- DESCRIPTION :
  --   Returns the sum of the differences of all elements |a(i,j)-b(i,j)|.

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
             ( q : Standard_Floating_Matrices.Matrix )
             return double_float is

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q 
  --   with all its following columns.

    sum,ip : double_float;

  begin
    sum := 0.0;
    for j in q'range(2) loop
      for k in j+1..q'last(2) loop
        ip := 0.0;
        for i in q'range(1) loop
          ip := ip + q(i,j)*q(i,k);
        end loop;
        sum := sum + abs(ip);
      end loop;
    end loop;
    return sum;
  end Orthogonality_Check_Sum;

  function Orthogonality_Check_Sum
             ( q : Double_Double_Matrices.Matrix ) return double_double is

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q 
  --   with all its following columns.

    sum,ip : double_double;
    zero : constant double_double := create(0.0);

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
             ( q : Quad_Double_Matrices.Matrix ) return quad_double is

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q 
  --   with all its following columns.

    sum,ip : quad_double;
    zero : constant quad_double := create(0.0);

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
             ( q : Standard_Complex_Matrices.Matrix )
             return double_float is

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q with
  --   its following columns.

    use Standard_Complex_Numbers;
    sum : double_float := 0.0;
    ip : Complex_Number;

  begin
    for j in q'range(2) loop
      for k in j+1..q'last(2) loop
        ip := Create(0.0);
        for i in q'range(1) loop
          ip := ip + Conjugate(q(i,j))*q(i,k);
        end loop;
        sum := sum + AbsVal(ip);
      end loop;
    end loop;
    return sum;
  end Orthogonality_Check_Sum;

  function Orthogonality_Check_Sum 
             ( q : DoblDobl_Complex_Matrices.Matrix )
             return double_double is

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q with
  --   its following columns.

    use DoblDobl_Complex_Numbers;
    zero : constant double_double := create(0.0);
    sum : double_double := zero;
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

  function Orthogonality_Check_Sum 
             ( q : QuadDobl_Complex_Matrices.Matrix )
             return quad_double is

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of any column in q with
  --   its following columns.

    use QuadDobl_Complex_Numbers;
    zero : constant quad_double := create(0.0);
    sum : quad_double := zero;
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

  function Orthogonality_Check_Sum
             ( q : Multprec_Complex_Matrices.Matrix )
             return Floating_Number is

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of a column with all
  --   its following columns.

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

  -- DESCRIPTION :
  --   Tests whether the columns are orthogonal w.r.t. each other,
  --   returns the sum of all inner products of a column with all
  --   its following columns.

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

  procedure Test_QRD ( a,q,r : in Standard_Floating_Matrices.Matrix;
                       output : in boolean ) is

    wrk : Standard_Floating_Matrices.Matrix(a'range(1),a'range(2));
    use Standard_Floating_Matrices;

  begin
    if output
     then put_line("The upper triangular part R :"); put(r,3);
    end if;
    wrk := q*r;
    if output
     then put_line("q*r :"); put(wrk,3); 
    end if;
    put("Difference in 1-norm between the matrix and q*r : ");
    put(Differences(a,wrk),3,3,3); new_line;
    put("Orthogonality check sum : ");
    put(Orthogonality_Check_Sum(q),3,3,3); new_line;
  end Test_QRD;

  procedure Test_QRD ( a,q,r : in Double_Double_Matrices.Matrix;
                       output : in boolean ) is

    wrk : Double_Double_Matrices.Matrix(a'range(1),a'range(2));
    use Double_Double_Matrices;

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

  procedure Test_QRD ( a,q,r : in Quad_Double_Matrices.Matrix;
                       output : in boolean ) is

    wrk : Quad_Double_Matrices.Matrix(a'range(1),a'range(2));
    use Quad_Double_Matrices;

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

  procedure Test_QRD ( a,q,r : in Standard_Complex_Matrices.Matrix;
                       output : in boolean ) is

    wrk : Standard_Complex_Matrices.Matrix(a'range(1),a'range(2));
    use Standard_Complex_Matrices;

  begin
    if output
     then put_line("The upper triangular part R :"); put(r,3);
    end if;
    wrk := q*r;
    if output
     then put_line("q*r :"); put(wrk,3); 
    end if;
    put("Difference in 1-norm between the matrix and q*r : ");
    put(Differences(a,wrk),3,3,3); new_line;
    put("Orthogonality check sum : ");
    put(Orthogonality_Check_Sum(q),3,3,3); new_line;
  end Test_QRD;

  procedure Test_QRD ( a,q,r : in DoblDobl_Complex_Matrices.Matrix;
                       output : in boolean ) is

    wrk : DoblDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));
    use DoblDobl_Complex_Matrices;

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

  procedure Test_QRD ( a,q,r : in QuadDobl_Complex_Matrices.Matrix;
                       output : in boolean ) is

    wrk : QuadDobl_Complex_Matrices.Matrix(a'range(1),a'range(2));
    use QuadDobl_Complex_Matrices;

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

-- STANDARD REAL TEST DRIVERS :

  procedure Standard_Real_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Standard_Floating_Matrices.Matrix;
                b : in Standard_Floating_Vectors.Vector;
                output : in boolean ) is

    wrk : Standard_Floating_Matrices.Matrix(1..n,1..m) := a;
    qraux : Standard_Floating_Vectors.Vector(1..m) := (1..m => 0.0);
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : Standard_Floating_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : Standard_Floating_Vectors.Vector(1..n);
    info : integer32;
    use Standard_Floating_Matrices;
    use Standard_Floating_Vectors;

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
    put("The norm of residual : "); put(Sum_Norm(dum),3,3,3); new_line;
  end Standard_Real_LS_Test;          

  procedure Standard_Real_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Standard_Floating_Matrices.Matrix;
                output : in boolean ) is

    wrk : Standard_Floating_Matrices.Matrix(1..n,1..m) := a;
    bas : Standard_Floating_Matrices.Matrix(1..n,1..n);
    qraux : Standard_Floating_Vectors.Vector(1..m) := (1..m => 0.0);
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
        bas(i,j) := 0.0;
      end loop;
    end loop;
    Basis(bas,a);
    if output
     then put_line("The orthogonal part Q of QR  :"); put(bas,3);
    end if;
    Test_QRD(a,bas,Extract_Upper_Triangular(wrk),output);
  end Standard_Real_QR_Test;

  procedure Standard_Interactive_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Floating_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    Standard_Real_QR_Test(n,m,piv,a,true);
  end Standard_Interactive_Real_QR_Test;

  procedure Standard_Interactive_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Floating_Matrices.Matrix(1..n,1..m);
    b : Standard_Floating_Vectors.Vector(1..n);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    put("Give right-hand size "); put(n,1);
    put_line("-vector : "); get(b);
    Standard_Real_LS_Test(n,m,piv,a,b,true);
  end Standard_Interactive_Real_LS_Test;

  procedure Standard_Random_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Floating_Matrices.Matrix(1..n,1..m);
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
      Standard_Real_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factoriziations on standard random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random Standard Real QR Factorizations");
  end Standard_Random_Real_QR_Test;

  procedure Standard_Random_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Floating_Matrices.Matrix(1..n,1..m);
    b : Standard_Floating_Vectors.Vector(1..n);
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
      Standard_Real_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" real least squares on standard random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Testing Standard Real Least Squares");
  end Standard_Random_Real_LS_Test;

-- DOBLDOBL REAL TEST DRIVERS :

  procedure DoblDobl_Real_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Double_Double_Matrices.Matrix;
                b : in Double_Double_Vectors.Vector;
                output : in boolean ) is

    zero : constant double_double := create(0.0);
    wrk : Double_Double_Matrices.Matrix(1..n,1..m) := a;
    qraux : Double_Double_Vectors.Vector(1..m) := (1..m => zero);
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : Double_Double_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : Double_Double_Vectors.Vector(1..n);
    info : integer32;
    use Double_Double_Matrices;
    use Double_Double_Vectors;

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
  end DoblDobl_Real_LS_Test;          

  procedure DoblDobl_Real_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Double_Double_Matrices.Matrix;
                output : in boolean ) is

    zero : constant double_double := create(0.0);
    wrk : Double_Double_Matrices.Matrix(1..n,1..m) := a;
    bas : Double_Double_Matrices.Matrix(1..n,1..n);
    qraux : Double_Double_Vectors.Vector(1..m) := (1..m => zero);
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
  end DoblDobl_Real_QR_Test;

  procedure DoblDobl_Interactive_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Double_Double_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    DoblDobl_Real_QR_Test(n,m,piv,a,true);
  end DoblDobl_Interactive_Real_QR_Test;

  procedure DoblDobl_Interactive_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Double_Double_Matrices.Matrix(1..n,1..m);
    b : Double_Double_Vectors.Vector(1..n);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    put("Give right-hand size "); put(n,1);
    put_line("-vector : "); get(b);
    DoblDobl_Real_LS_Test(n,m,piv,a,b,true);
  end DoblDobl_Interactive_Real_LS_Test;

  procedure DoblDobl_Random_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Double_Double_Matrices.Matrix(1..n,1..m);
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
      DoblDobl_Real_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factoriziations on dobldobl random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random DoblDobl Real QR Factorizations");
  end DoblDobl_Random_Real_QR_Test;

  procedure DoblDobl_Random_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Double_Double_Matrices.Matrix(1..n,1..m);
    b : Double_Double_Vectors.Vector(1..n);
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
      DoblDobl_Real_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" real least squares on dobldobl random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Testing DoblDobl Real Least Squares");
  end DoblDobl_Random_Real_LS_Test;

-- QUADDOBL REAL TEST DRIVERS :

  procedure QuadDobl_Real_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Quad_Double_Matrices.Matrix;
                b : in Quad_Double_Vectors.Vector;
                output : in boolean ) is

    zero : constant quad_double := create(0.0);
    wrk : Quad_Double_Matrices.Matrix(1..n,1..m) := a;
    qraux : Quad_Double_Vectors.Vector(1..m) := (1..m => zero);
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : Quad_Double_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : Quad_Double_Vectors.Vector(1..n);
    info : integer32;
    use Quad_Double_Matrices;
    use Quad_Double_Vectors;

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
  end QuadDobl_Real_LS_Test;          

  procedure QuadDobl_Real_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : in Quad_Double_Matrices.Matrix;
                output : in boolean ) is

    zero : constant quad_double := create(0.0);
    wrk : Quad_Double_Matrices.Matrix(1..n,1..m) := a;
    bas : Quad_Double_Matrices.Matrix(1..n,1..n);
    qraux : Quad_Double_Vectors.Vector(1..m) := (1..m => zero);
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
  end QuadDobl_Real_QR_Test;

  procedure QuadDobl_Interactive_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Quad_Double_Matrices.Matrix(1..n,1..m);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    QuadDobl_Real_QR_Test(n,m,piv,a,true);
  end QuadDobl_Interactive_Real_QR_Test;

  procedure QuadDobl_Interactive_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Quad_Double_Matrices.Matrix(1..n,1..m);
    b : Quad_Double_Vectors.Vector(1..n);

  begin
    put("Give a "); put(n,1); put("x"); put(m,1);   
    put_line(" matrix : "); get(a);
    put("Give right-hand size "); put(n,1);
    put_line("-vector : "); get(b);
    QuadDobl_Real_LS_Test(n,m,piv,a,b,true);
  end QuadDobl_Interactive_Real_LS_Test;

  procedure QuadDobl_Random_Real_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Quad_Double_Matrices.Matrix(1..n,1..m);
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
      QuadDobl_Real_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factoriziations on quaddobl random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random QuadDobl Real QR Factorizations");
  end QuadDobl_Random_Real_QR_Test;

  procedure QuadDobl_Random_Real_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Quad_Double_Matrices.Matrix(1..n,1..m);
    b : Quad_Double_Vectors.Vector(1..n);
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
      QuadDobl_Real_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" real least squares on quaddobl random real matrices.");
    new_line;
    print_times(Standard_Output,timer,"Testing QuadDobl Real Least Squares");
  end QuadDobl_Random_Real_LS_Test;

-- MULTPREC REAL TEST DRIVERS :

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

-- STANDARD COMPLEX TEST DRIVERS :

  procedure Standard_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : Standard_Complex_Matrices.Matrix;
                output : in boolean ) is

    use Standard_Complex_Numbers;
    wrk : Standard_Complex_Matrices.Matrix(1..n,1..m) := a;
    bas : Standard_Complex_Matrices.Matrix(1..n,1..n);
    qraux : Standard_Complex_Vectors.Vector(1..m) := (1..m => Create(0.0));
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
          bas(i,j) := Create(0.0);
        end loop;
      end loop;
      Basis(bas,a);
      if output
       then put_line("The orthogonal part Q of QR  :"); put(bas,3);
      end if;
      Test_QRD(a,bas,Extract_Upper_Triangular(wrk),output);
    end if;
  end Standard_Complex_QR_Test;

  procedure Standard_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : Standard_Complex_Matrices.Matrix;
                b : Standard_Complex_Vectors.Vector;
                output : in boolean ) is

    use Standard_Complex_Numbers;
    wrk : Standard_Complex_Matrices.Matrix(1..n,1..m) := a;
    qraux : Standard_Complex_Vectors.Vector(1..m) := (1..m => Create(0.0));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : Standard_Complex_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : Standard_Complex_Vectors.Vector(1..n);
    info : integer32;
    use Standard_Complex_Matrices;
    use Standard_Complex_Vectors; 

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
    put("Sum norm of residual : "); put(Sum_Norm(dum),3,3,3); new_line;
  end Standard_Complex_LS_Test;

  procedure Standard_Interactive_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Complex_Matrices.Matrix(1..n,1..m);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      Standard_Complex_QR_Test(n,m,piv,a,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Interactive_Complex_QR_Test;

  procedure Standard_Interactive_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Complex_Matrices.Matrix(1..n,1..m);
    b : Standard_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      put("Give right-hand size "); put(n,1);
      put_line("-vector : "); get(b); 
      Standard_Complex_LS_Test(n,m,piv,a,b,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Interactive_Complex_LS_Test;

  procedure Standard_Random_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Complex_Matrices.Matrix(1..n,1..m);
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
      Standard_Complex_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factorizations on random standard complex matrices.");
    new_line;
    print_times(Standard_Output,timer,
                "Random Standard Complex QR Factorizations");
  end Standard_Random_Complex_QR_Test;

  procedure Standard_Random_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : Standard_Complex_Matrices.Matrix(1..n,1..m);
    b : Standard_Complex_Vectors.Vector(1..n);
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
      Standard_Complex_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" least squares on random standard complex matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random Standard Complex Least Squares");
  end Standard_Random_Complex_LS_Test;

-- DOBLDOBL COMPLEX TEST DRIVERS :

  procedure DoblDobl_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : DoblDobl_Complex_Matrices.Matrix;
                output : in boolean ) is

    use DoblDobl_Complex_Numbers;
    wrk : DoblDobl_Complex_Matrices.Matrix(1..n,1..m) := a;
    bas : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant double_double := create(0.0);
    qraux : DoblDobl_Complex_Vectors.Vector(1..m) := (1..m => Create(zero));
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
  end DoblDobl_Complex_QR_Test;

  procedure DoblDobl_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : DoblDobl_Complex_Matrices.Matrix;
                b : DoblDobl_Complex_Vectors.Vector;
                output : in boolean ) is

    use DoblDobl_Complex_Numbers;
    wrk : DoblDobl_Complex_Matrices.Matrix(1..n,1..m) := a;
    zero : constant double_double := create(0.0);
    qraux : DoblDobl_Complex_Vectors.Vector(1..m) := (1..m => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : DoblDobl_Complex_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : DoblDobl_Complex_Vectors.Vector(1..n);
    info : integer32;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Vectors; 

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
  end DoblDobl_Complex_LS_Test;

  procedure DoblDobl_Interactive_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : DoblDobl_Complex_Matrices.Matrix(1..n,1..m);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      DoblDobl_Complex_QR_Test(n,m,piv,a,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Interactive_Complex_QR_Test;

  procedure DoblDobl_Interactive_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : DoblDobl_Complex_Matrices.Matrix(1..n,1..m);
    b : DoblDobl_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      put("Give right-hand size "); put(n,1);
      put_line("-vector : "); get(b); 
      DoblDobl_Complex_LS_Test(n,m,piv,a,b,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Interactive_Complex_LS_Test;

  procedure DoblDobl_Random_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : DoblDobl_Complex_Matrices.Matrix(1..n,1..m);
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
      DoblDobl_Complex_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factorizations on random double double complex matrices.");
    new_line;
    print_times(Standard_Output,timer,
                "Random DoblDobl Complex QR Factorizations");
  end DoblDobl_Random_Complex_QR_Test;

  procedure DoblDobl_Random_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : DoblDobl_Complex_Matrices.Matrix(1..n,1..m);
    b : DoblDobl_Complex_Vectors.Vector(1..n);
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
      DoblDobl_Complex_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" least squares on random double double complex matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random DoblDobl Complex Least Squares");
  end DoblDobl_Random_Complex_LS_Test;

-- QUADDOBL COMPLEX TEST DRIVERS :

  procedure QuadDobl_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean;
                a : QuadDobl_Complex_Matrices.Matrix;
                output : in boolean ) is

    use QuadDobl_Complex_Numbers;
    wrk : QuadDobl_Complex_Matrices.Matrix(1..n,1..m) := a;
    bas : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant quad_double := create(0.0);
    qraux : QuadDobl_Complex_Vectors.Vector(1..m) := (1..m => Create(zero));
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
  end QuadDobl_Complex_QR_Test;

  procedure QuadDobl_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean;
                a : QuadDobl_Complex_Matrices.Matrix;
                b : QuadDobl_Complex_Vectors.Vector;
                output : in boolean ) is

    use QuadDobl_Complex_Numbers;
    wrk : QuadDobl_Complex_Matrices.Matrix(1..n,1..m) := a;
    zero : constant quad_double := create(0.0);
    qraux : QuadDobl_Complex_Vectors.Vector(1..m) := (1..m => Create(zero));
    jpvt : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    sol : QuadDobl_Complex_Vectors.Vector(1..m);
    rsd,dum,dum2,dum3 : QuadDobl_Complex_Vectors.Vector(1..n);
    info : integer32;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Vectors; 

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
  end QuadDobl_Complex_LS_Test;

  procedure QuadDobl_Interactive_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      QuadDobl_Complex_QR_Test(n,m,piv,a,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Interactive_Complex_QR_Test;

  procedure QuadDobl_Interactive_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);
    b : QuadDobl_Complex_Vectors.Vector(1..n);
    ans : character;

  begin
    loop
      put("Give a "); put(n,1); put("x"); put(m,1);
      put_line(" matrix : "); get(a);
      put("Give right-hand size "); put(n,1);
      put_line("-vector : "); get(b); 
      QuadDobl_Complex_LS_Test(n,m,piv,a,b,true);
      put("Do you want more tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Interactive_Complex_LS_Test;

  procedure QuadDobl_Random_Complex_QR_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);
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
      QuadDobl_Complex_QR_Test(n,m,piv,a,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" QR factorizations on random quad double complex matrices.");
    new_line;
    print_times(Standard_Output,timer,
                "Random QuadDobl Complex QR Factorizations");
  end QuadDobl_Random_Complex_QR_Test;

  procedure QuadDobl_Random_Complex_LS_Test
              ( n,m : in integer32; piv : in boolean ) is

    a : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);
    b : QuadDobl_Complex_Vectors.Vector(1..n);
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
      QuadDobl_Complex_LS_Test(n,m,piv,a,b,output);
    end loop;
    tstop(timer);
    put("Tested "); put(nb,1);
    put_line(" least squares on random quad double complex matrices.");
    new_line;
    print_times(Standard_Output,timer,"Random QuadDobl Complex Least Squares");
  end QuadDobl_Random_Complex_LS_Test;

-- MULTPREC COMPLEX TEST DRIVERS :

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

-- MAIN PROGRAM :

  procedure Main is

    n,m : integer32 := 0;
    sz : natural32 := 0;
    choice : character;
    piv : constant boolean := false;

  begin
    new_line;
    put_line("Test on the QR decomposition and Least Squares.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. Exit this program.");
      put_line("*** QR and Least Squares in standard double precision ***");
      put_line("  1. QR-decomposition on given standard floating matrix.");
      put_line("  2.                           standard complex matrix.");
      put_line("  3.                  on random standard floating matrix.");
      put_line("  4.                            standard complex matrix.");
      put_line("  5. Least Squares on given standard floating matrix.");
      put_line("  6.                        standard complex matrix.");
      put_line("  7.               on random standard floating matrix.");
      put_line("  8.                         standard complex matrix.");
      put_line("*** QR and Least Squares in double double precision ***");
      put_line("  9. QR-decomposition on given dobldobl real matrix.");
      put_line("  A.                                    complex matrix.");
      put_line("  B.                  on random dobldobl real matrix.");
      put_line("  C.                                     complex matrix.");
      put_line("  D. Least Squares on given dobldobl real matrix.");
      put_line("  E.                        dobldobl complex matrix.");
      put_line("  F.               on random dobldobl real matrix.");
      put_line("  G.                                  complex matrix.");
      put_line("*** QR and Least Squares in quad double precision ***");
      put_line("  H. QR-decomposition on given quaddobl real matrix.");
      put_line("  I.                                    complex matrix.");
      put_line("  J.                  on random quaddobl real matrix.");
      put_line("  K.                                     complex matrix.");
      put_line("  L. Least Squares on given quaddobl real matrix.");
      put_line("  M.                                 complex matrix.");
      put_line("  N.               on random quaddobl real matrix.");
      put_line("  O.                                  complex matrix.");
      put_line("*** QR and Least Squares in arbitrary multiprecision ***");
      put_line("  P. QR-decomposition on given multprec floating matrix.");
      put_line("  Q.                  on random multprec floating matrix.");
      put_line("  R. Least Squares on given multprec floating matrix.");
      put_line("  S.               on random multprec floating matrix.");
      put_line("  T. QR-decomposition on given multprec complex matrix.");
      put_line("  U.                  on random multprec complex matrix.");
      put_line("  V. Least Squares on given multprec complex matrix.");
      put_line("  W.               on random multprec complex matrix.");
      put("Make your choice (0, 1, .. , A, B, .., W) : ");
      Ask_Alternative(choice,"0123456789ABCDEFGHIJKLMNOPQRSTUVW");
      exit when (choice = '0');
      new_line;
      put("Give the number of rows of the matrix : "); get(n);
      put("Give the number of columns of the matrix : "); get(m);
      if choice = 'P' or choice = 'Q' or choice = 'R' or choice = 'S'
        or choice = 'T' or choice = 'U' or choice = 'V' or choice = 'W'
       then put("Give the size of the numbers : "); get(sz);
      end if;
      case choice is
       -- QR and Least Squares in standard double arithmetic
        when '1' => Standard_Interactive_Real_QR_Test(n,m,piv);
        when '2' => Standard_Interactive_Complex_QR_Test(n,m,piv);
        when '3' => Standard_Random_Real_QR_Test(n,m,piv);
        when '4' => Standard_Random_Complex_QR_Test(n,m,piv);
        when '5' => Standard_Interactive_Real_LS_Test(n,m,piv);
        when '6' => Standard_Interactive_Complex_LS_Test(n,m,piv);
        when '7' => Standard_Random_Real_LS_Test(n,m,piv);
        when '8' => Standard_Random_Complex_LS_Test(n,m,piv);
       -- QR and Least Squares in double double arithmetic
        when '9' => DoblDobl_Interactive_Real_QR_Test(n,m,piv);
        when 'A' => DoblDobl_Interactive_Complex_QR_Test(n,m,piv);
        when 'B' => DoblDobl_Random_Real_QR_Test(n,m,piv);
        when 'C' => DoblDobl_Random_Complex_QR_Test(n,m,piv);
        when 'D' => DoblDobl_Interactive_Real_LS_Test(n,m,piv);
        when 'E' => DoblDobl_Interactive_Complex_LS_Test(n,m,piv);
        when 'F' => DoblDobl_Random_Real_LS_Test(n,m,piv);
        when 'G' => DoblDobl_Random_Complex_LS_Test(n,m,piv);
       -- QR and Least Squares in quad double arithmetic
        when 'H' => QuadDobl_Interactive_Real_QR_Test(n,m,piv);
        when 'I' => QuadDobl_Interactive_Complex_QR_Test(n,m,piv);
        when 'J' => QuadDobl_Random_Real_QR_Test(n,m,piv);
        when 'K' => QuadDobl_Random_Complex_QR_Test(n,m,piv);
        when 'L' => QuadDobl_Interactive_Real_LS_Test(n,m,piv);
        when 'M' => QuadDobl_Interactive_Complex_LS_Test(n,m,piv);
        when 'N' => QuadDobl_Random_Real_LS_Test(n,m,piv);
        when 'O' => QuadDobl_Random_Complex_LS_Test(n,m,piv);
       -- QR and Least Squares in arbitrary precision arithmetic
        when 'P' => Multprec_Interactive_Real_QR_Test(n,m,piv);
        when 'Q' => Multprec_Random_Real_QR_Test(n,m,sz,piv);
        when 'R' => Multprec_Interactive_Real_LS_Test(n,m,piv);
        when 'S' => Multprec_Random_Real_LS_Test(n,m,sz,piv);
        when 'T' => Multprec_Interactive_Complex_QR_Test(n,m,piv);
        when 'U' => Multprec_Random_Complex_QR_Test(n,m,sz,piv);
        when 'V' => Multprec_Interactive_Complex_LS_Test(n,m,piv);
        when 'W' => Multprec_Random_Complex_LS_Test(n,m,sz,piv);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_qrd;
