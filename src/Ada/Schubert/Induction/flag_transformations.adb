with Standard_Complex_Numbers;            use Standard_Complex_Numbers;

package body Flag_Transformations is

  function Flag_Transformation_Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

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
  end Flag_Transformation_Coefficient_Matrix;

  function Flag_Transformation_Right_Hand_Side
             ( n : integer32; g : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector is

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
  end Flag_Transformation_Right_Hand_Side;

  procedure Extract_Matrices
              ( n : in integer32; sol : in Standard_Complex_Vectors.Vector;
                A,T1,T2 : out Standard_Complex_Matrices.Matrix ) is

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

end Flag_Transformations;
