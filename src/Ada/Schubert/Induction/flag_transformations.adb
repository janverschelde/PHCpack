with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Standard_Matrix_Inversion;
with DoblDobl_Matrix_Inversion;
with QuadDobl_Matrix_Inversion;
with Setup_Flag_Homotopies;

package body Flag_Transformations is

-- DEFINING A LINEAR SYSTEM :

  function Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

    use Standard_Complex_Numbers;

    dim : constant integer32 := 2*n*n;
    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    offset : integer32 := 0;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for k in 0..(n-1) loop
      for i in 1..n loop
        for j in 1..n loop
          res(k*n+i,k*n+j) := f1(j,i);
          res(n*(n+k)+i,k*n+j) := f2(j,i);
        end loop;
      end loop;
    end loop;
    for i in 1..n-1 loop
      for j in 1..n-i loop
        for k in 0..(n-1) loop
         -- put("i = "); put(i,1);
         -- put("  j = "); put(j,1); put("  k = "); put(k,1);
         -- put("  row : "); put(k*n+i+j,1);
         -- put("  col : "); put(offset+n*n+j,1); new_line;
          res(k*n+i+j,offset+n*n+j) := -g1(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i;
    end loop;
   -- new_line;
    for i in 1..n loop
      for j in 1..n-i+1 loop
        for k in 0..(n-1) loop
         -- put("i = "); put(i,1);
         -- put("  j = "); put(j,1); put("  k = "); put(k,1);
         -- put("  row : "); put(n*(n+k)+i+j-1,1);
         -- put("  col : "); put(offset+n*n+j,1); new_line;
          res(n*(n+k)+i+j-1,offset+n*n+j) := -g2(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i + 1;
    end loop;
    return res;
  end Coefficient_Matrix;

  function Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := 2*n*n;
    res : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    offset : integer32 := 0;
    zero : constant double_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(zero);
      end loop;
    end loop;
    for k in 0..(n-1) loop
      for i in 1..n loop
        for j in 1..n loop
          res(k*n+i,k*n+j) := f1(j,i);
          res(n*(n+k)+i,k*n+j) := f2(j,i);
        end loop;
      end loop;
    end loop;
    for i in 1..n-1 loop
      for j in 1..n-i loop
        for k in 0..(n-1) loop
         -- put("i = "); put(i,1);
         -- put("  j = "); put(j,1); put("  k = "); put(k,1);
         -- put("  row : "); put(k*n+i+j,1);
         -- put("  col : "); put(offset+n*n+j,1); new_line;
          res(k*n+i+j,offset+n*n+j) := -g1(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i;
    end loop;
   -- new_line;
    for i in 1..n loop
      for j in 1..n-i+1 loop
        for k in 0..(n-1) loop
         -- put("i = "); put(i,1);
         -- put("  j = "); put(j,1); put("  k = "); put(k,1);
         -- put("  row : "); put(n*(n+k)+i+j-1,1);
         -- put("  col : "); put(offset+n*n+j,1); new_line;
          res(n*(n+k)+i+j-1,offset+n*n+j) := -g2(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i + 1;
    end loop;
    return res;
  end Coefficient_Matrix;

  function Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Matrices.Matrix is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := 2*n*n;
    res : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    offset : integer32 := 0;
    zero : constant quad_double := create(0.0);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(zero);
      end loop;
    end loop;
    for k in 0..(n-1) loop
      for i in 1..n loop
        for j in 1..n loop
          res(k*n+i,k*n+j) := f1(j,i);
          res(n*(n+k)+i,k*n+j) := f2(j,i);
        end loop;
      end loop;
    end loop;
    for i in 1..n-1 loop
      for j in 1..n-i loop
        for k in 0..(n-1) loop
         -- put("i = "); put(i,1);
         -- put("  j = "); put(j,1); put("  k = "); put(k,1);
         -- put("  row : "); put(k*n+i+j,1);
         -- put("  col : "); put(offset+n*n+j,1); new_line;
          res(k*n+i+j,offset+n*n+j) := -g1(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i;
    end loop;
   -- new_line;
    for i in 1..n loop
      for j in 1..n-i+1 loop
        for k in 0..(n-1) loop
         -- put("i = "); put(i,1);
         -- put("  j = "); put(j,1); put("  k = "); put(k,1);
         -- put("  row : "); put(n*(n+k)+i+j-1,1);
         -- put("  col : "); put(offset+n*n+j,1); new_line;
          res(n*(n+k)+i+j-1,offset+n*n+j) := -g2(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i + 1;
    end loop;
    return res;
  end Coefficient_Matrix;

  function Right_Hand_Side
             ( n : integer32; g : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;

    dim : constant integer32 := 2*n*n;
    res : Standard_Complex_Vectors.Vector(1..dim);
    ind : integer32 := 0;

  begin
    for i in g'range(1) loop
      for j in g'range(2) loop
        ind := ind + 1;
        res(ind) := g(i,j);
      end loop;
    end loop;
    for i in ind+1..res'last loop
      res(i) := Create(0.0);
    end loop;
    return res;
  end Right_Hand_Side;

  function Right_Hand_Side
             ( n : integer32; g : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Vectors.Vector is

    use DoblDobl_Complex_Numbers;

    dim : constant integer32 := 2*n*n;
    res : DoblDobl_Complex_Vectors.Vector(1..dim);
    ind : integer32 := 0;
    zero : constant double_double := create(0.0);

  begin
    for i in g'range(1) loop
      for j in g'range(2) loop
        ind := ind + 1;
        res(ind) := g(i,j);
      end loop;
    end loop;
    for i in ind+1..res'last loop
      res(i) := Create(zero);
    end loop;
    return res;
  end Right_Hand_Side;

  function Right_Hand_Side
             ( n : integer32; g : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Vectors.Vector is

    use QuadDobl_Complex_Numbers;

    dim : constant integer32 := 2*n*n;
    res : QuadDobl_Complex_Vectors.Vector(1..dim);
    ind : integer32 := 0;
    zero : constant quad_double := create(0.0);

  begin
    for i in g'range(1) loop
      for j in g'range(2) loop
        ind := ind + 1;
        res(ind) := g(i,j);
      end loop;
    end loop;
    for i in ind+1..res'last loop
      res(i) := Create(zero);
    end loop;
    return res;
  end Right_Hand_Side;

-- PROCESSING THE SOLUTION :

  procedure Extract_Matrices
              ( n : in integer32; sol : in Standard_Complex_Vectors.Vector;
                A,T1,T2 : out Standard_Complex_Matrices.Matrix ) is

    use Standard_Complex_Numbers;

    idx : integer32 := 0;

  begin
    for i in 1..n loop
      for j in 1..n loop
        idx := idx + 1;
        A(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T1(i,j) := Create(0.0);
      end loop;
      T1(i,i) := Create(1.0);
      for j in (i+1)..n loop
        idx := idx + 1;
        T1(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T2(i,j) := Create(0.0);
      end loop;
      for j in i..n loop
        idx := idx + 1;
        T2(i,j) := sol(idx);
      end loop;
    end loop;
  end Extract_Matrices;

  procedure Extract_Matrices
              ( n : in integer32; sol : in DoblDobl_Complex_Vectors.Vector;
                A,T1,T2 : out DoblDobl_Complex_Matrices.Matrix ) is

    use DoblDobl_Complex_Numbers;

    zero : constant double_double := create(0.0);
    one : constant double_double := create(1.0);
    idx : integer32 := 0;

  begin
    for i in 1..n loop
      for j in 1..n loop
        idx := idx + 1;
        A(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T1(i,j) := Create(zero);
      end loop;
      T1(i,i) := Create(one);
      for j in (i+1)..n loop
        idx := idx + 1;
        T1(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T2(i,j) := Create(zero);
      end loop;
      for j in i..n loop
        idx := idx + 1;
        T2(i,j) := sol(idx);
      end loop;
    end loop;
  end Extract_Matrices;

  procedure Extract_Matrices
              ( n : in integer32; sol : in QuadDobl_Complex_Vectors.Vector;
                A,T1,T2 : out QuadDobl_Complex_Matrices.Matrix ) is

    use QuadDobl_Complex_Numbers;

    zero : constant quad_double := create(0.0);
    one : constant quad_double := create(1.0);
    idx : integer32 := 0;

  begin
    for i in 1..n loop
      for j in 1..n loop
        idx := idx + 1;
        A(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T1(i,j) := Create(zero);
      end loop;
      T1(i,i) := Create(one);
      for j in (i+1)..n loop
        idx := idx + 1;
        T1(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T2(i,j) := Create(zero);
      end loop;
      for j in i..n loop
        idx := idx + 1;
        T2(i,j) := sol(idx);
      end loop;
    end loop;
  end Extract_Matrices;

-- THE MAIN TRANSFORMATION :

  procedure Transform
              ( n : in integer32; 
                f1,f2,g1,g2 : in Standard_Complex_Matrices.Matrix;
                A,T1,T2 : out Standard_Complex_Matrices.Matrix ) is

    dim : constant integer32 := 2*n*n;
    mat : Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Coefficient_Matrix(n,f1,f2,g1,g2);
    rhs : Standard_Complex_Vectors.Vector(1..dim)
        := Right_Hand_Side(n,g1);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;

  begin
    lufac(mat,dim,ipvt,info);
    lusolve(mat,dim,ipvt,rhs);
    Extract_Matrices(n,rhs,A,T1,T2);
  end Transform;

  procedure Transform
              ( n : in integer32; 
                f1,f2,g1,g2 : in DoblDobl_Complex_Matrices.Matrix;
                A,T1,T2 : out DoblDobl_Complex_Matrices.Matrix ) is

    dim : constant integer32 := 2*n*n;
    mat : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim)
        := Coefficient_Matrix(n,f1,f2,g1,g2);
    rhs : DoblDobl_Complex_Vectors.Vector(1..dim)
        := Right_Hand_Side(n,g1);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;

  begin
    lufac(mat,dim,ipvt,info);
    lusolve(mat,dim,ipvt,rhs);
    Extract_Matrices(n,rhs,A,T1,T2);
  end Transform;

  procedure Transform
              ( n : in integer32; 
                f1,f2,g1,g2 : in QuadDobl_Complex_Matrices.Matrix;
                A,T1,T2 : out QuadDobl_Complex_Matrices.Matrix ) is

    dim : constant integer32 := 2*n*n;
    mat : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim)
        := Coefficient_Matrix(n,f1,f2,g1,g2);
    rhs : QuadDobl_Complex_Vectors.Vector(1..dim)
        := Right_Hand_Side(n,g1);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;

  begin
    lufac(mat,dim,ipvt,info);
    lusolve(mat,dim,ipvt,rhs);
    Extract_Matrices(n,rhs,A,T1,T2);
  end Transform;

  function Residual
             ( f1,f2,g1,g2,A,T1,T2 : Standard_Complex_Matrices.Matrix )
             return double_float is

    res : double_float := 0.0;

    use Standard_Complex_Numbers;
    use Standard_Complex_Matrices;

    Af1 : constant Standard_Complex_Matrices.Matrix := A*f1;
    Af2 : constant Standard_Complex_Matrices.Matrix := A*f2;
    g1T1 : constant Standard_Complex_Matrices.Matrix := g1*T1;
    g2T2 : constant Standard_Complex_Matrices.Matrix := g2*T2;
    d1 : constant Standard_Complex_Matrices.Matrix := Af1 - g1T1;
    d2 : constant Standard_Complex_Matrices.Matrix := Af2 - g2T2;
    val : double_float;

  begin
    for i in d1'range(1) loop
      for j in d1'range(2) loop
        val := AbsVal(d1(i,j));
        res := res + val;
      end loop;
    end loop;
    for i in d2'range(1) loop
      for j in d2'range(2) loop
        val := AbsVal(d2(i,j));
        res := res + val;
      end loop;
    end loop;
    return res;
  end Residual;

  function Residual
             ( f1,f2,g1,g2,A,T1,T2 : DoblDobl_Complex_Matrices.Matrix )
             return double_double is

    res : double_double := create(0.0);

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Matrices;

    Af1 : constant DoblDobl_Complex_Matrices.Matrix := A*f1;
    Af2 : constant DoblDobl_Complex_Matrices.Matrix := A*f2;
    g1T1 : constant DoblDobl_Complex_Matrices.Matrix := g1*T1;
    g2T2 : constant DoblDobl_Complex_Matrices.Matrix := g2*T2;
    d1 : constant DoblDobl_Complex_Matrices.Matrix := Af1 - g1T1;
    d2 : constant DoblDobl_Complex_Matrices.Matrix := Af2 - g2T2;
    val : double_double;

  begin
    for i in d1'range(1) loop
      for j in d1'range(2) loop
        val := AbsVal(d1(i,j));
        res := res + val;
      end loop;
    end loop;
    for i in d2'range(1) loop
      for j in d2'range(2) loop
        val := AbsVal(d2(i,j));
        res := res + val;
      end loop;
    end loop;
    return res;
  end Residual;

  function Residual
             ( f1,f2,g1,g2,A,T1,T2 : QuadDobl_Complex_Matrices.Matrix )
             return quad_double is

    res : quad_double := create(0.0);

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Matrices;

    Af1 : constant QuadDobl_Complex_Matrices.Matrix := A*f1;
    Af2 : constant QuadDobl_Complex_Matrices.Matrix := A*f2;
    g1T1 : constant QuadDobl_Complex_Matrices.Matrix := g1*T1;
    g2T2 : constant QuadDobl_Complex_Matrices.Matrix := g2*T2;
    d1 : constant QuadDobl_Complex_Matrices.Matrix := Af1 - g1T1;
    d2 : constant QuadDobl_Complex_Matrices.Matrix := Af2 - g2T2;
    val : quad_double;

  begin
    for i in d1'range(1) loop
      for j in d1'range(2) loop
        val := AbsVal(d1(i,j));
        res := res + val;
      end loop;
    end loop;
    for i in d2'range(1) loop
      for j in d2'range(2) loop
        val := AbsVal(d2(i,j));
        res := res + val;
      end loop;
    end loop;
    return res;
  end Residual;

-- GENERAL WRAPPERS :

  function Move_to_Generic_Flag
              ( n : integer32; G : Standard_Complex_Matrices.Matrix )
              return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..n,1..n) := G;
    moved : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
          := Setup_Flag_Homotopies.Moved_Flag(n);
    idemat : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           :=  Setup_Flag_Homotopies.Identity(n);
    A,T1,T2,invT1 : Standard_Complex_Matrices.Matrix(1..n,1..n);

    use Standard_Complex_Matrices;

  begin
    Transform(n,moved,idemat,idemat,res,A,T1,T2);
    invT1 := Standard_Matrix_Inversion.Inverse(T1);
    res := invT1*res;
    return res;
  end Move_to_Generic_Flag;

  function Move_to_Generic_Flag
              ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

    ranflag : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
            := Setup_Flag_Homotopies.Random_Flag(n);
    res : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
        := Move_to_Generic_Flag(n,ranflag);

  begin
    return res;
  end Move_to_Generic_Flag;

  procedure Move_to_Generic_Flag
              ( n : in integer32;
                G : in Standard_Complex_Matrices.Matrix;
                F : out Standard_Complex_Matrices.Matrix;
                rsd : out double_float ) is

    moved : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
          := Setup_Flag_Homotopies.Moved_Flag(n);
    idemat : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           :=  Setup_Flag_Homotopies.Identity(n);
    ranflag : constant Standard_Complex_Matrices.Matrix(1..n,1..n) := G;
    A,T1,T2,invT1 : Standard_Complex_Matrices.Matrix(1..n,1..n);

    use Standard_Complex_Matrices;

  begin
    Transform(n,moved,idemat,idemat,ranflag,A,T1,T2);
    invT1 := Standard_Matrix_Inversion.Inverse(T1);
    F := invT1*ranflag;
    rsd := Residual(moved,idemat,idemat,ranflag,A,T1,T2);
  end Move_to_Generic_Flag;

  procedure Move_to_Generic_Flag
              ( n : in integer32;
                F : out Standard_Complex_Matrices.Matrix;
                rsd : out double_float ) is

    ranflag : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
            := Setup_Flag_Homotopies.Random_Flag(n);

  begin
    Move_to_Generic_Flag(n,ranflag,F,rsd);
  end Move_to_Generic_Flag;

-- FOR APPLICATION TO RESOLVE SCHUBERT PROBLEMS :

  procedure Transform_Sequence_with_Flag
              ( n,i : in integer32;
                flags : in out Standard_Complex_VecMats.VecMat;
                A,invA,sT : out Standard_Complex_Matrices.Matrix ) is

    moved : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
          := Setup_Flag_Homotopies.Moved_Flag(n);
    idemat : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
           :=  Setup_Flag_Homotopies.Identity(n);
    flag,T2 : Standard_Complex_Matrices.Matrix(1..n,1..n);

    use Standard_Complex_Matrices;

  begin
    flag := flags(i).all;
    Transform(n,moved,idemat,idemat,flag,A,sT,T2);
    invA := Standard_Matrix_Inversion.Inverse(A);
    for j in i+1..flags'last loop
      flags(j).all := invA*flags(j).all;
    end loop;
  end Transform_Sequence_with_Flag;

  procedure Transform_Sequence_with_Flag
              ( n,i : in integer32;
                flags : in out DoblDobl_Complex_VecMats.VecMat;
                A,invA,sT : out DoblDobl_Complex_Matrices.Matrix ) is

    moved : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
          := Setup_Flag_Homotopies.Moved_Flag(n);
    idemat : constant DoblDobl_Complex_Matrices.Matrix(1..n,1..n)
           :=  Setup_Flag_Homotopies.Identity(n);
    flag,T2 : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);

    use DoblDobl_Complex_Matrices;

  begin
    flag := flags(i).all;
    Transform(n,moved,idemat,idemat,flag,A,sT,T2);
    invA := DoblDobl_Matrix_Inversion.Inverse(A);
    for j in i+1..flags'last loop
      flags(j).all := invA*flags(j).all;
    end loop;
  end Transform_Sequence_with_Flag;

  procedure Transform_Sequence_with_Flag
              ( n,i : in integer32;
                flags : in out QuadDobl_Complex_VecMats.VecMat;
                A,invA,sT : out QuadDobl_Complex_Matrices.Matrix ) is

    moved : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
          := Setup_Flag_Homotopies.Moved_Flag(n);
    idemat : constant QuadDobl_Complex_Matrices.Matrix(1..n,1..n)
           :=  Setup_Flag_Homotopies.Identity(n);
    flag,T2 : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);

    use QuadDobl_Complex_Matrices;

  begin
    flag := flags(i).all;
    Transform(n,moved,idemat,idemat,flag,A,sT,T2);
    invA := QuadDobl_Matrix_Inversion.Inverse(A);
    for j in i+1..flags'last loop
      flags(j).all := invA*flags(j).all;
    end loop;
  end Transform_Sequence_with_Flag;

  procedure Create ( n : in integer32;
                     flags : in Standard_Complex_VecMats.VecMat;
                     stack : out Standard_Stack_of_Flags;
                     A,invA,sT : out Standard_Complex_VecMats.VecMat ) is

    AA,invAA,T1 : Standard_Complex_Matrices.Matrix(1..n,1..n);
   -- work : Standard_Complex_VecMats.VecMat(flags'range);

  begin
    for i in flags'first..flags'last-1 loop
      declare
        work : Standard_Complex_VecMats.VecMat(flags'range);
        flag : Standard_Complex_Matrices.Matrix(1..n,1..n);
        ptr2flags : Standard_Complex_VecMats.Link_to_VecMat;
      begin
        for j in work'range loop
          if i = flags'first then  -- use original flags at the start
            Standard_Complex_Matrices.Copy(flags(j).all,flag);
          else -- work with transformed flags in later stages
            ptr2flags := stack(i-1);
            Standard_Complex_Matrices.Copy(ptr2flags(j).all,flag);
          end if;
          work(j) := new Standard_Complex_Matrices.Matrix'(flag);
        end loop;
        Transform_Sequence_with_Flag(n,i,work,AA,invAA,T1);
        stack(i) := new Standard_Complex_VecMats.VecMat'(work);
        A(i) := new Standard_Complex_Matrices.Matrix'(AA);
        invA(i) := new Standard_Complex_Matrices.Matrix'(invAA);
        st(i) := new Standard_Complex_Matrices.Matrix'(T1);
      end;
    end loop;
  end Create;

  procedure Create ( n : in integer32;
                     flags : in DoblDobl_Complex_VecMats.VecMat;
                     stack : out DoblDobl_Stack_of_Flags;
                     A,invA,sT : out DoblDobl_Complex_VecMats.VecMat ) is

    AA,invAA,T1 : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
   -- work : DoblDobl_Complex_VecMats.VecMat(flags'range);

  begin
    for i in flags'first..flags'last-1 loop
      declare
        work : DoblDobl_Complex_VecMats.VecMat(flags'range);
        flag : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
        ptr2flags : DoblDobl_Complex_VecMats.Link_to_VecMat;
      begin
        for j in work'range loop
          if i = flags'first then  -- use original flags at the start
            DoblDobl_Complex_Matrices.Copy(flags(j).all,flag);
          else -- work with transformed flags in later stages
            ptr2flags := stack(i-1);
            DoblDobl_Complex_Matrices.Copy(ptr2flags(j).all,flag);
          end if;
          work(j) := new DoblDobl_Complex_Matrices.Matrix'(flag);
        end loop;
        Transform_Sequence_with_Flag(n,i,work,AA,invAA,T1);
        stack(i) := new DoblDobl_Complex_VecMats.VecMat'(work);
        A(i) := new DoblDobl_Complex_Matrices.Matrix'(AA);
        invA(i) := new DoblDobl_Complex_Matrices.Matrix'(invAA);
        st(i) := new DoblDobl_Complex_Matrices.Matrix'(T1);
      end;
    end loop;
  end Create;

  procedure Create ( n : in integer32;
                     flags : in QuadDobl_Complex_VecMats.VecMat;
                     stack : out QuadDobl_Stack_of_Flags;
                     A,invA,sT : out QuadDobl_Complex_VecMats.VecMat ) is

    AA,invAA,T1 : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
   -- work : QuadDobl_Complex_VecMats.VecMat(flags'range);

  begin
    for i in flags'first..flags'last-1 loop
      declare
        work : QuadDobl_Complex_VecMats.VecMat(flags'range);
        flag : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
        ptr2flags : QuadDobl_Complex_VecMats.Link_to_VecMat;
      begin
        for j in work'range loop
          if i = flags'first then  -- use original flags at the start
            QuadDobl_Complex_Matrices.Copy(flags(j).all,flag);
          else -- work with transformed flags in later stages
            ptr2flags := stack(i-1);
            QuadDobl_Complex_Matrices.Copy(ptr2flags(j).all,flag);
          end if;
          work(j) := new QuadDobl_Complex_Matrices.Matrix'(flag);
        end loop;
        Transform_Sequence_with_Flag(n,i,work,AA,invAA,T1);
        stack(i) := new QuadDobl_Complex_VecMats.VecMat'(work);
        A(i) := new QuadDobl_Complex_Matrices.Matrix'(AA);
        invA(i) := new QuadDobl_Complex_Matrices.Matrix'(invAA);
        st(i) := new QuadDobl_Complex_Matrices.Matrix'(T1);
      end;
    end loop;
  end Create;

  procedure Clear ( s : in out Standard_Stack_of_Flags ) is
  begin
    for i in s'range loop
      Standard_Complex_VecMats.Deep_Clear(s(i));
    end loop;
  end Clear;

  procedure Clear ( s : in out DoblDobl_Stack_of_Flags ) is
  begin
    for i in s'range loop
      DoblDobl_Complex_VecMats.Deep_Clear(s(i));
    end loop;
  end Clear;

  procedure Clear ( s : in out QuadDobl_Stack_of_Flags ) is
  begin
    for i in s'range loop
      QuadDobl_Complex_VecMats.Deep_Clear(s(i));
    end loop;
  end Clear;

end Flag_Transformations;
