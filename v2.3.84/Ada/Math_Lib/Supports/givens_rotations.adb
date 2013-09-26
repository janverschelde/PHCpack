with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Givens_Rotations is

  procedure Givens_Factors ( v: in Vector; i,j : in integer32;
                             cosi,sine : out double_float ) is
    norm : double_float;
 
  begin
    norm := SQRT(v(i)*v(i) + v(j)*v(j));
    cosi := v(i)/norm;  sine := v(j)/norm;
  end Givens_Factors; 

  procedure Givens_Factors ( mat : in Matrix; i,j : in integer32;
                             cosi,sine : out double_float ) is
    norm : double_float;

  begin
    norm := SQRT(mat(i,i)*mat(i,i) + mat(j,i)*mat(j,i));
    cosi := mat(i,i)/norm;  sine := mat(j,i)/norm;
  end Givens_Factors;

  procedure Givens_Rotation ( v : in out Vector; i,j : in integer32;
                              cosi,sine : in double_float ) is
    temp : double_float;

  begin
    temp := v(i);
    v(i) := cosi*v(i) + sine*v(j);
    v(j) := -sine*temp + cosi*v(j);
  end Givens_Rotation;

  procedure Givens_Rotation ( mat : in out Matrix; i,j : in integer32;
                              cosi,sine : in double_float ) is
    temp : double_float;

  begin
   -- for k in mat'range(2) loop -- certain angle preservation
    for k in i..mat'last(2) loop -- only angle preservation if upper triangular
      temp := mat(i,k);
      mat(i,k) := cosi*mat(i,k) + sine*mat(j,k);
      mat(j,k) := -sine*temp + cosi*mat(j,k);
    end loop;
  end Givens_Rotation;

  procedure Givens_Rotation ( mat : in out Matrix; lastcol,i,j : in integer32;
                              cosi,sine : in double_float ) is
    temp : double_float;

  begin
   -- for k in mat'first(2)..lastcol loop -- certain angle preservation
    for k in i..lastcol loop -- only angle preservation if upper triangular
      temp := mat(i,k);
      mat(i,k) := cosi*mat(i,k) + sine*mat(j,k);
      mat(j,k) := -sine*temp + cosi*mat(j,k);
    end loop;
  end Givens_Rotation;

  procedure Givens_Rotation ( v : in out Vector; i,j : in integer32 ) is

    cosi,sine : double_float;

  begin
    Givens_Factors(v,i,j,cosi,sine);
    Givens_Rotation(v,i,j,cosi,sine);
  end Givens_Rotation;

  procedure Givens_Rotation ( mat : in out Matrix; i,j : in integer32 ) is

    cosi,sine : double_float;

  begin
    Givens_Factors(mat,i,j,cosi,sine);
    Givens_Rotation(mat,i,j,cosi,sine);
  end Givens_Rotation;

  procedure Givens_Rotation ( mat : in out Matrix;
                              lastcol,i,j : in integer32 ) is

    cosi,sine : double_float;

  begin
    Givens_Factors(mat,i,j,cosi,sine);
    Givens_Rotation(mat,lastcol,i,j,cosi,sine);
  end Givens_Rotation;

  procedure Upper_Triangulate
              ( row : in integer32; mat : in out Matrix;
                tol : in double_float;
                ipvt : in out Standard_Integer_Vectors.Vector;
                pivot : out integer32 ) is

    piv,tpi : integer32 := 0;
    tmp : double_float;

  begin
    for j in row..mat'last(2) loop  -- search pivot
      if abs(mat(row,j)) > tol
       then piv := j;
      end if;
      exit when (piv /= 0);
    end loop;
    if piv /= 0 then                    -- zero row
      if piv /= row then                  -- interchange columns
        for k in mat'range(1) loop
          tmp := mat(k,row); mat(k,row) := mat(k,piv); mat(k,piv) := tmp;
        end loop;
        tpi := ipvt(row); ipvt(row) := ipvt(piv); ipvt(piv) := tpi;
      end if;
      for k in row+1..mat'last(1) loop  -- perform Givens rotations
        if abs(mat(k,row)) > tol
         then Givens_Rotation(mat,row,k);
        end if;
      end loop;
    end if;
    pivot := piv;
  end Upper_Triangulate;

  procedure Upper_Triangulate ( mat : in out Matrix; tol : in double_float ) is

    pivot : integer32;
    temp : double_float;

  begin
    for i in mat'range(1) loop     -- mat(mat'first(1)..i-1,mat'first(2)..i-1)
      pivot := 0;                  -- is already upper triangular
      for j in i..mat'last(2) loop -- search pivot
        if abs(mat(i,j)) > tol
         then pivot := j;
        end if;
        exit when (pivot /= 0);
      end loop;
      exit when (pivot = 0);              -- zero row
      if pivot /= i then                  -- interchange columns
        for k in mat'range(1) loop
          temp := mat(k,i); mat(k,i) := mat(k,pivot); mat(k,pivot) := temp;
        end loop;
      end if;
      for k in i+1..mat'last(1) loop  -- perform Givens rotations
        if abs(mat(k,i)) > tol
         then Givens_Rotation(mat,i,k);
        end if;
      end loop;                    -- mat(mat'first(1)..i,mat'first(2)..i)
      exit when i > mat'last(2);   -- is upper triangular
    end loop;
  end Upper_Triangulate;

  procedure Upper_Triangulate ( mat : in out Matrix; tol : in double_float;
                                ipvt : out Standard_Integer_Vectors.Vector ) is

    pivot,tpi : integer32;
    pivots : Standard_Integer_Vectors.Vector(mat'range(2));
    temp : double_float;

  begin
    for i in pivots'range loop     -- initialize vector of pivots
      pivots(i) := i;
    end loop;
    for i in mat'range(1) loop     -- mat(mat'first(1)..i-1,mat'first(2)..i-1)
      pivot := 0;                  -- is already upper triangular
      for j in i..mat'last(2) loop -- search pivot
        if abs(mat(i,j)) > tol
         then pivot := j;
        end if;
        exit when (pivot /= 0);
      end loop;
      exit when (pivot = 0);          -- zero row
      if pivot /= i then              -- interchange columns
        for k in mat'range(1) loop
          temp := mat(k,i); mat(k,i) := mat(k,pivot); mat(k,pivot) := temp;
        end loop;
        tpi := pivots(i); pivots(i) := pivots(pivot); pivots(pivot) := tpi;
      end if;
      for k in i+1..mat'last(1) loop  -- perform Givens rotations
        if abs(mat(k,i)) > tol
         then Givens_Rotation(mat,i,k);
        end if;
      end loop;                    -- mat(mat'first(1)..i,mat'first(2)..i)
      exit when i > mat'last(2);   -- is upper triangular
    end loop;
    ipvt := pivots;
  end Upper_Triangulate;

  procedure Upper_Triangulate ( mat : in out Matrix; rhs : in out Vector;
                                tol : in double_float ) is

    pivot : integer32;
    temp,cosi,sine : double_float;

  begin
    for i in mat'range(1) loop     -- mat(mat'first(1)..i-1,mat'first(2)..i-1)
      pivot := 0;                  -- is already upper triangular
      for j in i..mat'last(2) loop -- search pivot
        if abs(mat(i,j)) > tol
         then pivot := j;
        end if;
        exit when (pivot /= 0);
      end loop;
      exit when (pivot = 0);          -- zero row
      if pivot /= i then              -- interchange columns
        for k in mat'range(1) loop
          temp := mat(k,i); mat(k,i) := mat(k,pivot); mat(k,pivot) := temp;
        end loop;
      end if;
      for k in i+1..mat'last(1) loop  -- perform Givens rotations
        if abs(mat(k,i)) > tol then
          Givens_Factors(mat,i,k,cosi,sine);
          Givens_Rotation(mat,i,k,cosi,sine);
          Givens_Rotation(rhs,i,k,cosi,sine);
        end if;
      end loop;                    -- mat(mat'first(1)..i,mat'first(2)..i)
      exit when i > mat'last(2);   -- is upper triangular
    end loop;
  end Upper_Triangulate;

  procedure Upper_Triangulate ( mat : in out Matrix; rhs : in out Vector;
                                tol : in double_float;
                                ipvt : out Standard_Integer_Vectors.Vector ) is

    pivot,tpi : integer32;
    pivots : Standard_Integer_Vectors.Vector(mat'range(2));
    temp,cosi,sine : double_float;

  begin
    for i in pivots'range loop     -- initialize vector of pivots
      pivots(i) := i;
    end loop;
    for i in mat'range(1) loop     -- mat(mat'first(1)..i-1,mat'first(2)..i-1)
      pivot := 0;                  -- is already upper triangular
      for j in i..mat'last(2) loop -- search pivot
        if abs(mat(i,j)) > tol
         then pivot := j;
        end if;
        exit when (pivot /= 0);
      end loop;
      exit when (pivot = 0);          -- zero row
      if pivot /= i then              -- interchange columns
        for k in mat'range(1) loop
          temp := mat(k,i); mat(k,i) := mat(k,pivot); mat(k,pivot) := temp;
        end loop;
        tpi := pivots(i); pivots(i) := pivots(pivot); pivots(pivot) := tpi;
      end if;
      for k in i+1..mat'last(1) loop  -- perform Givens rotations
        if abs(mat(k,i)) > tol then
          Givens_Factors(mat,i,k,cosi,sine);
          Givens_Rotation(mat,i,k,cosi,sine);
          Givens_Rotation(rhs,i,k,cosi,sine);
        end if;
      end loop;                    -- mat(mat'first(1)..i,mat'first(2)..i)
      exit when i > mat'last(2);   -- is upper triangular
    end loop;
    ipvt := pivots;
  end Upper_Triangulate;

  procedure Solve ( mat : in Matrix; rhs : in Vector; tol : in double_float;
                    x : out Vector ) is

    rank : integer32 := 0;
    sol : vector(mat'range(1)) := (mat'range(1) => 0.0);

  begin
    for i in mat'range(1) loop  -- determine the rank of the matrix
      if abs(mat(i,i)) > tol
       then rank := rank + 1;
      end if;
      exit when ((i >= mat'last(2)) or (abs(mat(i,i)) < tol));
    end loop;
    for i in reverse mat'first(1)..rank loop -- back substitution
      for j in i+1..rank loop
        sol(i) := sol(i) + mat(i,j)*sol(j);
      end loop;
      sol(i) := rhs(i) - sol(i);
      if abs(mat(i,i)) > tol then
        sol(i) := sol(i)/mat(i,i);
      elsif abs(sol(i)) < tol then
        sol(i) := 1.0;
      else 
        return;
      end if;
    end loop;                 -- invariant : sol(i..rank) computed
    x := sol;
  end Solve;

  procedure Solve ( mat : in Matrix; rhs : in Vector; tol : in double_float;
                    rank : in integer32; x : out Vector ) is

    sol : vector(mat'range(1)) := (mat'range(1) => 0.0);

  begin
    for i in reverse mat'first(1)..rank loop -- back substitution
      for j in i+1..rank loop
        sol(i) := sol(i) + mat(i,j)*sol(j);
      end loop;
      sol(i) := rhs(i) - sol(i);
      if abs(mat(i,i)) > tol then
        sol(i) := sol(i)/mat(i,i);
      elsif abs(sol(i)) < tol then
        sol(i) := 1.0;
      else
        return;
      end if;
    end loop;                 -- invariant : sol(i..rank) computed
    x := sol;
  end Solve;

end Givens_Rotations;
