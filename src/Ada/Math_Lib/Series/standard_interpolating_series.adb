with text_io;                             use text_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_Singular_Values;

package body Standard_Interpolating_Series is

  function Eval ( v : Standard_Dense_Vector_Series.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector is

    vec : Standard_Complex_Vectors.Link_to_Vector := v.cff(0);
    res : Standard_Complex_Vectors.Vector(vec'range) := vec.all;
    pwt : Complex_Number := Create(1.0);
 
    use Standard_Complex_Vectors;

  begin
    for k in 1..v.deg loop
      pwt := pwt*t;
      res := res + pwt*v.cff(k).all;
    end loop;
    return res;
  end Eval;

  function Eval ( m : Standard_Dense_Matrix_Series.Matrix;
                  t : Complex_Number )
                return Standard_Complex_Matrices.Matrix is

    mat : Standard_Complex_Matrices.Link_to_Matrix := m.cff(0);
    res : Standard_Complex_Matrices.Matrix(mat'range(1),mat'range(2))
        := mat.all;
    pwt : Complex_Number := Create(1.0);

  begin
    for k in 1..m.deg loop
      pwt := pwt*t;
      mat := m.cff(k);
      for i in mat'range(1) loop
        for j in mat'range(2) loop
          res(i,j) := res(i,j) + pwt*mat(i,j);
        end loop;
      end loop;
    end loop;
    return res;
  end Eval;

  function Full_Rank
             ( m : Standard_Dense_Matrix_Series.Matrix;
               d : integer32; verbose : boolean := true ) return boolean is

    t : constant Complex_Number := Standard_Random_Numbers.Random1;
    mat : Standard_Complex_Matrices.Link_to_Matrix := m.cff(0);
    val : Standard_Complex_Matrices.Matrix(mat'range(1),mat'range(2))
        := mat.all;
    pwt : Complex_Number := Create(1.0);
    n : constant integer32 := val'last(1);
    p : constant integer32 := val'last(2);
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Standard_Complex_Vectors.Vector(1..mm);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Standard_Complex_Matrices.Matrix(1..n,1..n);
    v : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    rnk,info : integer32;
    tol : constant double_float := 1.0E-8;

  begin
    for k in 1..d loop
      pwt := pwt*t;
      mat := m.cff(k);
      for i in mat'range(1) loop
        for j in mat'range(2) loop
          val(i,j) := val(i,j) + pwt*mat(i,j);
        end loop;
      end loop;
    end loop;
    Standard_Complex_Singular_Values.SVD(val,n,p,s,e,u,v,job,info);
    rnk := Standard_Complex_Singular_Values.Rank(s,tol);
    if verbose then
      put_line("The singular values : "); put_line(s);
      put("The numerical rank : "); put(rnk,1); new_line;
    end if;
    return (rnk = n);
  end Full_Rank;

  function Full_Rank
             ( m : Standard_Dense_Matrix_Series.Matrix;
               verbose : boolean := true ) return integer32 is
  begin
    for d in 0..m.deg loop
      if Full_Rank(m,d)
       then return d;
      end if;
    end loop;
    return -1;
  end Full_Rank;

end Standard_Interpolating_Series;
