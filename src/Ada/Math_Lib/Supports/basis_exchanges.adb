with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;

package body Basis_Exchanges is

-- AUXILIARY TO ENUMERATE :

  function Is_In ( v : Standard_Integer_Vectors.Vector; k : integer32 ) 
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if there exits an index i such that v(i) = k.

  begin
    for i in v'range loop
      if v(i) = k
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

-- AUXILIARIES TO UPDATE :

  function Multipliers ( matcol : Standard_Floating_Vectors.Vector;
                         pivot : integer32 )
                       return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the multipliers that are needed to update the inverse
  --   of the basis matrix from the column in the coefficient matrix.
  --   The vector on entry is binv*mat(col,*).

  -- REQUIRED : abs(matcol(pivot)) > tol.

    res : Standard_Floating_Vectors.Vector(matcol'range);

  begin
    for j in res'range loop
      if j = pivot
       then res(j) := 1.0/matcol(pivot);
       else res(j) := -matcol(j)/matcol(pivot);
      end if;
    end loop;
    return res;
  end Multipliers;

  procedure Multiply ( binv : in out Standard_Floating_Matrices.Matrix;
                       mult : in Standard_Floating_Vectors.Vector;
                       pivot : in integer32 ) is

  -- DESCRIPTION :
  --   Performs the multiplicative update on the inverse of the basis.

  -- REQUIRED : mult'range = binv'range(1).

    new_col : Standard_Floating_Vectors.Vector(binv'range(2));

  begin
    for j in binv'range(2) loop
      for i in mult'range loop
        if i = pivot
         then new_col(i) := mult(i)*binv(pivot,j);
         else new_col(i) := mult(i)*binv(pivot,j) + binv(i,j);
        end if;
      end loop;
      for i in binv'range(1) loop
        binv(i,j) := new_col(i);
      end loop;
    end loop;
  end Multiply;

-- TARGET PROCEDURES :

  procedure Initial_Basis
               ( n,m : in integer32;
                 mat : in Standard_Floating_Matrices.Matrix;
                 tol : in double_float;
                 binv : out Standard_Floating_Matrices.Matrix;
                 cols : out Standard_Integer_Vectors.Vector;
                 fail : out boolean ) is

    invbas : Standard_Floating_Matrices.Matrix(1..n,1..n);
    matcols,used : Standard_Integer_Vectors.Vector(1..n);
    cffmat : Standard_Floating_Matrices.Matrix(1..n,1..m+n);
    todo : boolean;

  begin
    for i in 1..n loop                        -- initialization
      matcols(i) := i;
      for j in 1..n loop
        if i = j
         then invbas(i,j) := 1.0; cffmat(i,j) := 1.0;
         else invbas(i,j) := 0.0; cffmat(i,j) := 0.0;
        end if;
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..m loop
        cffmat(i,j+n) := mat(i,j);
      end loop;
    end loop;
    fail := false;
    used := (used'range => 0);
    for i in 1..n loop                        -- replace ith element in basis
      todo := true;
      for j in n+1..m+n loop
        if not Is_In(matcols,j) and not Is_In(used,j) then
          matcols(i) := j;
          Update(n,m+n,invbas,cffmat,matcols,i,j,tol,todo);
          if not todo
           then used(i) := j;
          end if;
        end if;
        exit when (used(i) > 0);
      end loop;
      exit when todo;
    end loop;
    binv := invbas;
    for i in 1..n loop
      cols(i) := matcols(i) - n;
    end loop;
    fail := todo;
  end Initial_Basis;

  function Solve ( n : integer32;
                   binv : Standard_Floating_Matrices.Matrix;
                   rhs : Standard_Floating_Vectors.Vector;
                   cols : Standard_Integer_Vectors.Vector )
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(1..n);

  begin
    for i in 1..n loop
      res(i) := 0.0;
      for j in 1..n loop
        res(i) := res(i) + binv(j,i)*rhs(cols(j));
      end loop;
    end loop;
    return res;
  end Solve;

  procedure Update ( n,m : in integer32;
                     binv : in out Standard_Floating_Matrices.Matrix;
                     mat : in Standard_Floating_Matrices.Matrix;
                     cols : in Standard_Integer_Vectors.Vector;
                     binv_row,mat_col : in integer32; tol : in double_float;
                     fail : out boolean ) is

    binvmat,mult : Standard_Floating_Vectors.Vector(1..n);

  begin
    for i in 1..n loop                  -- compute binv*mat(*,col)
      binvmat(i) := 0.0;
      for j in 1..n loop
        binvmat(i) := binvmat(i) + binv(i,j)*mat(j,mat_col);
      end loop;
    end loop;
    if abs(binvmat(binv_row)) < tol     -- check pivot
     then -- Column_Basis(n,mat,cols,binv,fail);
          fail := true;
     else fail := false;
          mult := Multipliers(binvmat,binv_row);
          Multiply(binv,mult,binv_row);
    end if;
  end Update;

-- RESTARTING PROCEDURES :

  procedure Column_Basis
                ( n : in integer32;
                  mat : in Standard_Floating_Matrices.Matrix;
                  cols : in Standard_Integer_Vectors.Vector;
                  binv : out Standard_Floating_Matrices.Matrix;
                  fail : out boolean ) is

    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sys : Standard_Floating_Matrices.Matrix(1..n,1..n);
    rhs : Standard_Floating_Vectors.Vector(1..n);

  begin
    for i in 1..n loop                      -- ith constraint is
      for j in 1..n loop                    -- the ith column of mat
        sys(i,j) := mat(i,cols(j));         -- taking cols into account
      end loop;
    end loop;
    lufac(sys,n,ipvt,info);
    if info /= 0 then
      fail := true;
    else
      fail := false;
      for i in 1..n loop
        rhs := (1..n => 0.0);
        rhs(i) := 1.0;
        lusolve(sys,n,ipvt,rhs);
        for j in 1..n loop
          binv(j,i) := rhs(j);
        end loop;
      end loop;
    end if;
  end Column_Basis;

  procedure Column_Solve
                ( n : in integer32;
                  mat : in Standard_Floating_Matrices.Matrix;
                  cols : in Standard_Integer_Vectors.Vector;
                  rhs : in Standard_Floating_Vectors.Vector;
                  sol : out Standard_Floating_Vectors.Vector;
                  fail : out boolean ) is

    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    sys : Standard_Floating_Matrices.Matrix(1..n,1..n);

  begin
    for i in 1..n loop                      -- ith constraint is
      for j in 1..n loop                    -- the ith column of mat
        sys(i,j) := mat(j,cols(i));         -- taking cols into account
      end loop;
    end loop;
    lufac(sys,n,ipvt,info);
    if info /= 0 then
      fail := true;
    else
      fail := false;
      for i in 1..n loop
        sol(i) := rhs(cols(i));
      end loop;
      lusolve(sys,n,ipvt,sol);
    end if;
  end Column_Solve;

-- ENUMERATOR :

  procedure Enumerate_Basis_Inverses
                ( n,m : in integer32;
                  mat,binv : in Standard_Floating_Matrices.Matrix;
                  mat_cols : in Standard_Integer_Vectors.Vector;
                  tol : in double_float ) is

    continue : boolean := true;

    procedure Enumerate ( bas : in Standard_Floating_Matrices.Matrix;
                          act : in Standard_Integer_Vectors.Vector;
                          start,level : in integer32 ) is

    -- DESCRIPTION :
    --   The search for a passive index to come in is started from start.

      new_bas : Standard_Floating_Matrices.Matrix(1..n,1..n);
      new_act : Standard_Integer_Vectors.Vector(act'range) := act;
      ind : integer32;
      fail : boolean;

    begin
      for i in start..m loop
        if not Is_In(act,i) then
          for j in 1..n loop                -- exchange i for binv(j)
            new_bas := bas;
            ind := act(j);
            new_act(j) := i;
            Update(n,m,new_bas,mat,new_act,j,i,tol,fail);
            if not fail then
              Report(new_bas,new_act,level,continue);
              if continue
               then Enumerate(new_bas,new_act,i+1,level+1);
              end if;
              new_act(j) := ind;
            end if;
            exit when not continue;
          end loop;
        end if;
        exit when not continue;
      end loop;
    end Enumerate;

  begin
    Enumerate(binv,mat_cols,mat_cols'last,0);
  end Enumerate_Basis_Inverses;

end Basis_Exchanges;
