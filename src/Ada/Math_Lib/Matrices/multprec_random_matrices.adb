with Standard_Random_Numbers;
with Multprec_Integer_Numbers;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;

package body Multprec_Random_Matrices is

  function Random_Matrix ( n,m : natural32; low,upp : integer32 )
                         return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..integer32(n),1..integer32(m));
    rnd : integer32;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        rnd := Standard_Random_Numbers.Random(low,upp);
        res(i,j) := Multprec_Integer_Numbers.Create(rnd);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Integer64_Matrices.Matrix is

    res : Multprec_Integer64_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Floating_Matrices.Matrix is

    res : Multprec_Floating_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

end Multprec_Random_Matrices;
