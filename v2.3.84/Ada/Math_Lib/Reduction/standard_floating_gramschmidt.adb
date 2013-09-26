with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Floating_Vector_Norms;    use Standard_Floating_Vector_Norms;

package body Standard_Floating_GramSchmidt is

-- DESIGN DECISIONS :
--   The modified Gram-Schmidt orthonormalization method is a vector
--   oriented method: data structures are vectors of vectors, not matrices.
--   No memory allocation or deallocation is done in the QR routine.
--   The inner products are inlined.
--   For performance, double indexing q(k)(i) must be avoided:
--   working with qk and vj as below doubled the execution speed!

  procedure QR ( n,m : in integer32;
                 q,r : in out Standard_Floating_VecVecs.VecVec ) is

    w,t : double_float;
    qj,qk,rj,rk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in 1..m loop
      qk := q(k);
      w := 0.0;
      for i in 1..n loop
        w := w + qk(i)**2;
      end loop;
      w := sqrt(w);
      rk := r(k);
      rk(k) := w;
      for i in 1..n loop
        qk(i) := qk(i)/w;
      end loop;
      for j in k+1..m loop
        qj := q(j);
        rj := r(j);
        w := 0.0;
        for i in 1..n loop
          w := w + qk(i)*qj(i);
        end loop;
        rj(k) := w;
        for i in 1..n loop
          t := qj(i) - w*qk(i);
          qj(i) := t;
        end loop;
      end loop;
    end loop;
  end QR;

  procedure Select_Pivot_Column
               ( k,m : in integer32;
                 v : in Standard_Floating_VecVecs.VecVec;
                 ind : out integer32; max : out double_float ) is

    nrm : double_float;

  begin
    max := Max_Norm(v(k).all);
    ind := k;
    for i in k+1..m loop
      nrm := Max_Norm(v(i).all);
      if nrm > max
       then max := nrm; ind := i;
      end if;
    end loop;
  end Select_Pivot_Column;

  procedure Swap ( v : in out Standard_Floating_VecVecs.VecVec;
                   i,j : in integer32 ) is

    w : constant Standard_Floating_Vectors.Link_to_Vector := v(i);

  begin
    v(i) := v(j);
    v(j) := w;
  end Swap;

  procedure Permute_Upper_Triangular
               ( n : in integer32;
                 r : in out Standard_Floating_VecVecs.VecVec;
                 p : in Standard_Integer_Vectors.Vector ) is

  -- NOTE :
  --   To "undo" the effect of the permutation on the upper triangular
  --   matrix R, we must read the permutation backwards and first undo
  --   the most recent permutation.

    tmp : double_float;

  begin
    for i in reverse p'range loop    -- do last permutations first
      if p(i) > 0 then
        if p(i) /= i then
          for j in i..n loop         -- permute only after ith row
            tmp := r(i)(j);
            r(i)(j) := r(p(i))(j);
            r(p(i))(j) := tmp;
          end loop;
        end if;
      end if;
    end loop;
  end Permute_Upper_Triangular;

  procedure Revert_Swaps
               ( n : in integer32;
                 q : in out Standard_Floating_VecVecs.VecVec;
                 p : in Standard_Integer_Vectors.Vector ) is
  begin
    for i in reverse p'range loop    -- do last permutations first
      if p(i) > 0 then
        if p(i) /= i then
          Swap(q,i,p(i));
        end if;
      end if;
    end loop;
  end Revert_Swaps;

  procedure QR ( n,m : in integer32; tol : in double_float;
                 q,r : in out Standard_Floating_VecVecs.VecVec;
                 pivots : out Standard_Integer_Vectors.Vector;
                 rank : out integer32 ) is

    w,t,max : double_float;
    qk,qj,rj,rk : Standard_Floating_Vectors.Link_to_Vector;
    ind : integer32;

  begin
    pivots := (1..m => 0);
    rank := 0;
    for k in 1..m loop
      Select_Pivot_Column(k,m,q,ind,max);
      if max < tol
       then rank := k; return;
      end if;
      pivots(k) := ind;
      if k /= ind
       then Swap(q,k,ind);
      end if;
      qk := q(k);
      w := 0.0;
      for i in 1..n loop
        w := w + qk(i)**2;
      end loop;
      w := sqrt(w);
      rk := r(k);
      rk(k) := w;
      for i in 1..n loop
        qk(i) := qk(i)/w;
      end loop;
      for j in k+1..m loop
        qj := q(j);
        rj := r(j);
        w := 0.0;
        for i in 1..n loop
          w := w + qk(i)*qj(i);
        end loop;
        rj(k) := w;
        for i in 1..n loop
          t := qj(i) - w*qk(i);
          qj(i) := t;
        end loop;
      end loop;
    end loop;
    rank := m;
  end QR;

  procedure QR ( n,m : in integer32;
                 v : in out Standard_Floating_Matrices.Matrix;
                 q,r : out Standard_Floating_Matrices.Matrix ) is

    w : double_float;

  begin
    for k in 1..m loop
      w := 0.0;
      for i in 1..n loop
        w := w + v(k,i)**2;
      end loop;
      w := sqrt(w);
      r(k,k) := w;
      for i in 1..n loop
        q(k,i) := v(k,i)/w;
      end loop;
      for j in k+1..m loop
        w := 0.0;
        for i in 1..n loop
          w := w + q(k,i)*v(j,i);
        end loop;
        r(j,k) := w;
        for i in 1..n loop
          v(j,i) := v(j,i) - w*q(k,i);
        end loop;
      end loop;
    end loop;
  end QR;

  function maximum ( x,y : double_float ) return double_float is
  begin
    if x > y
     then return x;
     else return y;
    end if;
  end maximum;

  procedure Test_Orthonormality
              ( n,m : in integer32;
                q : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float; output : in boolean;
                maxerr : out double_float; fail : out boolean ) is

    w,e : double_float;

  begin
    maxerr := 0.0;
    fail := false;
    for k in 1..m loop
      w := 0.0;
      for i in 1..n loop 
        w := w + q(k)(i)*q(k)(i);
      end loop;
      e := abs(w - 1.0);
      if output then
        put("norm^2 of q("); put(k,1); put(") : "); put(w);
        if e > tol
         then put_line("  fail!");
         else put_line("  okay.");
        end if;
      end if;
      maxerr := maximum(maxerr,e);
      fail := (fail or (e > tol)); 
      for j in k+1..m loop
        w := 0.0;
        for i in 1..n loop
          w := w + q(k)(i)*q(j)(i);
        end loop;
        e := abs(w);
        if output then
          put("  |q("); put(k,1); put(")*q("); put(j,1);
          put(")| = "); put(w);
          if e > tol
           then put_line("  fail!");
           else put_line("  okay.");
          end if;
        end if;
        maxerr := maximum(maxerr,e);
        fail := (fail or (abs(w) > tol));
      end loop;
    end loop;
  end Test_Orthonormality;

  procedure Test_Decomposition
              ( n,m : in integer32;
                v,q,r : in Standard_Floating_VecVecs.VecVec;
                tol : in double_float; output : in boolean;
                maxerr : out double_float; fail : out boolean ) is

    w,e : double_float;

  begin
    maxerr := 0.0;
    fail := false;
    for i in 1..n loop
      for j in 1..m loop
        w := 0.0;
        for k in 1..m loop
          w := w + q(k)(i)*r(j)(k);
        end loop;
        e := abs(v(j)(i) - w);
        if output then
          put("A("); put(i,1); put(","); put(j,1); put(") : ");
          put(v(j)(i)); put(" = "); put(w); put(" ?");
          if e > tol
           then put_line("  fail!");
           else put_line("  okay.");
          end if;
        end if;
        maxerr := maximum(maxerr,e);
        fail := (fail or (e > tol));
      end loop;
    end loop;
  end Test_Decomposition;

  function Matrix_Product
              ( n,m : integer32;
                v : Standard_Floating_VecVecs.VecVec;
                x : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..n) := (1..n => 0.0);
    vj : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for j in 1..m loop   -- column oriented version
      vj := v(j);
      for i in 1..n loop
        res(i) := res(i) + vj(i)*x(j);
      end loop;
    end loop;
    return res;
  end Matrix_Product;

  function Matrix_Projection
              ( n,m : integer32;
                q : Standard_Floating_VecVecs.VecVec;
                b : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..m);
    w : double_float;

  begin
    for k in 1..m loop
      w := 0.0;
      for i in 1..n loop
        w := w + q(k)(i)*b(i);
      end loop;
      res(k) := w;
    end loop;
    return res;
  end Matrix_Projection; 

  procedure Orthogonal_Projection
              ( n,m : in integer32;
                q : in Standard_Floating_VecVecs.VecVec;
                b : in out Standard_Floating_Vectors.Vector;
                qb : out Standard_Floating_Vectors.Vector ) is

    w : double_float;
    qk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in 1..m loop
      qk := q(k);
      w := 0.0;
      for i in 1..n loop
        w := w + qk(i)*b(i);
      end loop;
      for i in 1..n loop
        b(i) := b(i) - w*qk(i);
      end loop;
      qb(k) := w;
    end loop;
  end Orthogonal_Projection; 

  function Solve ( m : integer32;
                   r : Standard_Floating_VecVecs.VecVec;
                   b : Standard_Floating_Vectors.Vector )
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..m) := b;
    rk : Standard_Floating_Vectors.Link_to_Vector;
    w : double_float;

  begin
    for k in reverse 1..m loop
      rk := r(k);
      w := res(k)/rk(k);
      for i in 1..(k-1) loop
        res(i) := res(i) - rk(i)*w;
      end loop;
      res(k) := w;
    end loop;
    return res;
  end Solve;

end Standard_Floating_GramSchmidt;
