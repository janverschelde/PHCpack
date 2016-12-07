with text_io;                             use text_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Double_Double_Numbers;               use Double_Double_Numbers;
with DoblDobl_Random_Numbers;
with DoblDobl_Random_Vectors;
with Standard_Integer_Vectors;
with Double_Double_Vectors_io;            use Double_Double_Vectors_io;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs_io;         use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Matrices_io;        use DoblDobl_Complex_Matrices_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Linear_Solvers;
with DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Singular_Values;

package body DoblDobl_Interpolating_Series is

  function Eval ( v : DoblDobl_Dense_Vector_Series.Vector;
                  t : Complex_Number )
                return DoblDobl_Complex_Vectors.Vector is

    vec : DoblDobl_Complex_Vectors.Link_to_Vector := v.cff(0);
    res : DoblDobl_Complex_Vectors.Vector(vec'range) := vec.all;
    pwt : Complex_Number := Create(integer32(1));
 
    use DoblDobl_Complex_Vectors;

  begin
    for k in 1..v.deg loop
      pwt := pwt*t;
      res := res + pwt*v.cff(k).all;
    end loop;
    return res;
  end Eval;

  function Eval ( m : DoblDobl_Dense_Matrix_Series.Matrix;
                  t : Complex_Number )
                return DoblDobl_Complex_Matrices.Matrix is

    mat : DoblDobl_Complex_Matrices.Link_to_Matrix := m.cff(0);
    res : DoblDobl_Complex_Matrices.Matrix(mat'range(1),mat'range(2))
        := mat.all;
    pwt : Complex_Number := Create(integer32(1));

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

  function Rank ( A : DoblDobl_Complex_Matrices.Matrix;
                  tol : double_float ) return integer32 is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    mm : constant integer32 := DoblDobl_Complex_Singular_Values.Min0(n+1,p);
    s : DoblDobl_Complex_Vectors.Vector(1..mm);
    e : DoblDobl_Complex_Vectors.Vector(1..p);
    u : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : DoblDobl_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;
    x : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2)) := A;

  begin
    DoblDobl_Complex_Singular_Values.SVD(x,n,p,s,e,u,v,job,info);
    return DoblDobl_Complex_Singular_Values.Rank(s,tol);
  end Rank;

  function Full_Rank
             ( m : DoblDobl_Dense_Matrix_Series.Matrix;
               d : integer32; verbose : boolean := true ) return boolean is

    t : constant Complex_Number := DoblDobl_Random_Numbers.Random1;
    mat : DoblDobl_Complex_Matrices.Link_to_Matrix := m.cff(0);
    val : DoblDobl_Complex_Matrices.Matrix(mat'range(1),mat'range(2))
        := mat.all;
    pwt : Complex_Number := Create(integer32(1));
    n : constant integer32 := val'last(1);
    p : constant integer32 := val'last(2);
    mm : constant integer32 := DoblDobl_Complex_Singular_Values.Min0(n+1,p);
    s : DoblDobl_Complex_Vectors.Vector(1..mm);
    e : DoblDobl_Complex_Vectors.Vector(1..p);
    u : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : DoblDobl_Complex_Matrices.Matrix(1..p,1..p);
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
    DoblDobl_Complex_Singular_Values.SVD(val,n,p,s,e,u,v,job,info);
    rnk := DoblDobl_Complex_Singular_Values.Rank(s,tol);
    if verbose then
      put_line("The singular values : "); put_line(s);
      put("The numerical rank : "); put(rnk,1); new_line;
    end if;
    return (rnk = n);
  end Full_Rank;

  function Full_Rank
             ( m : DoblDobl_Dense_Matrix_Series.Matrix;
               verbose : boolean := true ) return integer32 is
  begin
    for d in 0..m.deg loop
      if Full_Rank(m,d)
       then return d;
      end if;
    end loop;
    return -1;
  end Full_Rank;

  function Sample ( v : DoblDobl_Dense_Vector_Series.Vector;
                    t : DoblDobl_Complex_Vectors.Vector )
                  return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(t'range);
    lv : DoblDobl_Complex_Vectors.Link_to_Vector := v.cff(0);

  begin
    for i in t'range loop
      declare
        e : constant DoblDobl_Complex_Vectors.Vector(lv'range)
          := Eval(v,t(i));
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector'(e);
      end;
    end loop;
    return res;
  end Sample;

  function Sample ( m : DoblDobl_Dense_Matrix_Series.Matrix;
                    t : DoblDobl_Complex_Vectors.Vector )
                  return DoblDobl_Complex_VecMats.VecMat is

    res : DoblDobl_Complex_VecMats.VecMat(t'range);
    lm : DoblDobl_Complex_Matrices.Link_to_Matrix := m.cff(0);

  begin
    for i in t'range loop
      declare
        e : constant DoblDobl_Complex_Matrices.Matrix(lm'range(1),lm'range(2))
          := Eval(m,t(i));
      begin
        res(i) := new DoblDobl_Complex_Matrices.Matrix'(e);
      end;
    end loop;
    return res;
  end Sample;

  function Solve_Linear_Systems
             ( m : DoblDobl_Complex_VecMats.VecMat;
               v : DoblDobl_Complex_VecVecs.VecVec )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(m'range);
    dim : constant integer32 := v(0)'last;
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    sol : DoblDobl_Complex_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in v'range loop
      sol := v(i).all;
      wrk := m(i).all;
      DoblDobl_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
      DoblDobl_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,sol);
      res(i) := new DoblDobl_Complex_Vectors.Vector'(sol);
    end loop;
    return res;
  end Solve_Linear_Systems;

  function Residuals
             ( m : DoblDobl_Complex_VecMats.VecMat;
               v,x : DoblDobl_Complex_VecVecs.VecVec )
             return Double_Double_Vectors.Vector is

    res : Double_Double_Vectors.Vector(v'range);
    dim : constant integer32 := v(0)'last;
    lm : DoblDobl_Complex_Matrices.Link_to_Matrix;
    lv,lx : DoblDobl_Complex_Vectors.Link_to_Vector;
    wrk : DoblDobl_Complex_Vectors.Vector(1..dim);

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

  begin
    for i in res'range loop
      lm := m(i);
      lv := v(i);
      lx := x(i);
      wrk := lv.all - lm.all*lx.all;
      res(i) := DoblDobl_Complex_Vector_Norms.Max_Norm(wrk);
    end loop;
    return res;
  end Residuals;

  function Transpose ( x : DoblDobl_Complex_VecVecs.VecVec )
                     return DoblDobl_Complex_VecVecs.VecVec is

    deg : constant integer32 := x'last;
    dim : constant integer32 := x(0)'last;
    res : DoblDobl_Complex_VecVecs.VecVec(1..dim);

  begin
    for i in res'range loop
      res(i) := new DoblDobl_Complex_Vectors.Vector(1..deg+1);
    end loop;
    for i in x'range loop
      for j in x(i)'range loop
        res(j)(i+1) := x(i)(j);
      end loop;
    end loop;
    return res;
  end Transpose;

  function Vandermonde_Matrix
             ( t : DoblDobl_Complex_Vectors.Vector )
             return DoblDobl_Complex_Matrices.Matrix is

    dim : constant integer32 := t'last - t'first + 1;
    res : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    idx : integer32 := 0;

  begin
    for i in t'range loop
      idx := idx + 1;
      res(idx,1) := Create(integer32(1));
      for j in 2..dim loop
        res(idx,j) := res(idx,j-1)*t(i);
      end loop;
    end loop;
    return res;
  end Vandermonde_Matrix;

  function Solve_Interpolation_Systems
             ( v : DoblDobl_Complex_Matrices.Matrix;
               f : DoblDobl_Complex_VecVecs.VecVec )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(f'range);
    dim : constant integer32 := v'last(1);
    wrk : DoblDobl_Complex_Matrices.Matrix(v'range(1),v'range(2)) := v;
    sol : DoblDobl_Complex_Vectors.Vector(1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;

  begin
    DoblDobl_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
    for i in f'range loop
      sol := f(i).all;
      DoblDobl_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,sol);
      res(i) := new DoblDobl_Complex_Vectors.Vector'(sol);
    end loop;
    return res;
  end Solve_Interpolation_Systems;

  function Construct ( x : DoblDobl_Complex_VecVecs.VecVec )
                     return DoblDobl_Dense_Vector_Series.Vector is

    res : DoblDobl_Dense_Vector_Series.Vector;
    dim : constant integer32 := x'last;
    lx0 : constant DoblDobl_Complex_Vectors.Link_to_Vector := x(x'first); 
    deg : constant integer32 := lx0'last - lx0'first;

  begin
    res.deg := deg;
    for i in 0..res.deg loop
      res.cff(i) := new DoblDobl_Complex_Vectors.Vector(1..dim);
    end loop;
    for i in x'range loop
      for j in x(i)'range loop
        res.cff(j-1)(i) := x(i)(j);
      end loop;
    end loop;
    return res;
  end Construct;

  function Interpolate
             ( mat : DoblDobl_Dense_Matrix_Series.Matrix;
               rhs : DoblDobl_Dense_Vector_Series.Vector;
               verbose : boolean := true )
             return DoblDobl_Dense_Vector_Series.Vector is

    res : DoblDobl_Dense_Vector_Series.Vector;
    dim : constant integer32 := rhs.cff(0)'last;
    t : constant DoblDobl_Complex_Vectors.Vector(0..mat.deg)
      := DoblDobl_Random_Vectors.Random_Vector(0,mat.deg);
    m : DoblDobl_Complex_VecMats.VecMat(t'range) := Sample(mat,t);
    v : DoblDobl_Complex_VecVecs.VecVec(t'range) := Sample(rhs,t);
    x : DoblDobl_Complex_VecVecs.VecVec(t'range);
    r : Double_Double_Vectors.Vector(t'range);
    xt : DoblDobl_Complex_VecVecs.VecVec(1..dim);
    vdm : DoblDobl_Complex_Matrices.Matrix(1..mat.deg+1,1..mat.deg+1);
    cff : DoblDobl_Complex_VecVecs.VecVec(1..dim);

  begin
    if verbose then
      put_line("The sample points :");
      put_line(t);
    end if;
    x := Solve_Linear_Systems(m,v);
    r := Residuals(m,v,x);
    if verbose then
      put_line("The solutions to the interpolated linear systems :");
      put_line(x);
      put_line("The residuals of the solved sampled linear systems :");
      put_line(r);
    end if;
    xt := Transpose(x);
    if verbose then
      put_line("The transposed solution vectors :");
      put_line(xt);
    end if;
    vdm := Vandermonde_Matrix(t);
    cff := Solve_Interpolation_Systems(vdm,xt);
    if verbose then
      put_line("The coefficients computed via interpolation :");
      put_line(cff);
    end if;
    res := Construct(cff);
    return res;
  end Interpolate;

  function factorial ( k : integer32 ) return Complex_Number is

    fac : integer32 := 1;
    dd_fac : double_double;
    res : Complex_Number;

  begin
    for i in 2..k loop
      fac := i*fac;
    end loop;
    dd_fac := Double_Double_Numbers.Create(fac);
    res := Create(dd_fac);
    return res;
  end factorial;

  function Diff ( m : DoblDobl_Complex_VecMats.VecMat;
                  t : Complex_Number; pow,ord : integer32 )
                return DoblDobl_Complex_Matrices.Matrix is

    lm0 : DoblDobl_Complex_Matrices.Link_to_Matrix := m(0);
    dim : constant integer32 := lm0'last(1);
    res : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    npw : constant natural := natural(pow);
    idx : integer32;
    pwt,fac : Complex_Number;

    use DoblDobl_Complex_Matrices;

  begin
    if ord = 0 then      -- evaluate
      res := m(0).all;
      if pow = 0 then
        pwt := Create(integer32(1));
      else
        pwt := t**npw;
        res := pwt*res;
      end if;
      for i in 1..m'last loop
        pwt := pwt*t;
        res := res + pwt*m(i).all;
      end loop;
    elsif pow = 0 then   -- differentiate, start at index ord
      res := factorial(ord)*m(ord).all;
      pwt := Create(integer32(1));
      for j in ord+1..m'last loop
        pwt := pwt*t;
        fac := factorial(ord)/factorial(j-ord)*pwt;
        res := res + fac*m(j).all;
      end loop;
    else  -- differentiate with multiplication of power of t
      idx := ord-pow;
      if idx < 0
       then idx := 0;
      end if;
      res := m(idx).all;
      if pow <= ord then
        pwt := Create(integer32(1));
        res := factorial(pow+idx)*m(idx).all;
      else
        pwt := t**natural(pow-ord);
        fac := factorial(pow+idx)/factorial(pow+idx-ord);
        res := fac*pwt*m(idx).all;
      end if;
      for i in idx+1..m'last loop
        pwt := pwt*t;
        fac := factorial(pow+i)/factorial(pow+i-ord);
        res := res + fac*pwt*m(i).all;
      end loop;
    end if;
    return res;
  end Diff;

  function Hermite_Matrix
             ( m : DoblDobl_Complex_VecMats.VecMat;
               t : Complex_Number )
             return DoblDobl_Complex_Matrices.Matrix is

    lmt : constant DoblDobl_Complex_Matrices.Link_to_Matrix := m(0);
    adim : constant integer32 := lmt'last(1);
    rdim : constant integer32 := adim*(m'last+1);
    res : DoblDobl_Complex_Matrices.Matrix(1..rdim,1..rdim);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..adim,1..adim);

  begin
    for col in m'range loop      -- columns are the powers of t
      for row in m'range loop    -- rows are the differentiation orders
        wrk := Diff(m,t,col,row);
        for i in wrk'range(1) loop
          for j in wrk'range(2) loop
            res(row*adim+i,col*adim+j) := wrk(i,j);
          end loop;
        end loop;
      end loop;
    end loop;
    return res;
  end Hermite_Matrix;

  function Hermite_Vector
             ( v : DoblDobl_Complex_VecVecs.VecVec;
               t : Complex_Number )
             return DoblDobl_Complex_Vectors.Vector is

    lv0 : constant DoblDobl_Complex_Vectors.Link_to_Vector := v(0);
    adim : constant integer32 := lv0'last; 
    rdim : constant integer32 := adim*(v'last+1);
    res : DoblDobl_Complex_Vectors.Vector(1..rdim);
    wrk : DoblDobl_Complex_Vectors.Vector(lv0'range);
    pwt : Complex_Number := Create(integer32(1));
    fac : Complex_Number;

    use DoblDobl_Complex_Vectors;

  begin
    wrk := v(0).all;
    for i in 1..v'last loop         -- evaluate at t
      pwt := pwt*t;
      wrk := wrk + pwt*v(i).all;
    end loop;
    for i in 1..adim loop
      res(i) := wrk(i);
    end loop;
    for i in 1..v'last loop         -- the i-th derivative
      wrk := factorial(i)*v(i).all;
      pwt := Create(integer32(1));
      for j in i+1..v'last loop
        pwt := pwt*t;
        fac := factorial(i)/factorial(j-i)*pwt;
        wrk := wrk + fac*v(j).all;
      end loop;
      for j in 1..adim loop
        res(i*adim+j) := wrk(j);
      end loop;
    end loop;
    return res;
  end Hermite_Vector;

  function Hermite_Interpolate
             ( mat : DoblDobl_Dense_Matrix_Series.Matrix;
               rhs : DoblDobl_Dense_Vector_Series.Vector;
               t : Complex_Number; verbose : boolean := true )
             return DoblDobl_Dense_Vector_Series.Vector is

    res : DoblDobl_Dense_Vector_Series.Vector;
    deg : constant integer32 := mat.deg;
    dim : constant integer32 := mat.cff(0)'last(1);
    bigdim : constant integer32 := (deg+1)*dim;
    A : DoblDobl_Complex_Matrices.Matrix(1..bigdim,1..bigdim)
      := Hermite_Matrix(mat.cff(0..deg),t);
    b : DoblDobl_Complex_Vectors.Vector(1..bigdim)
      := Hermite_Vector(rhs.cff(0..deg),t);
    info,rnk : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..bigdim);
    wrk : DoblDobl_Complex_Vectors.Vector(1..dim);

  begin
    if verbose then
      put_line("The coefficient matrix :"); put(A,3);
      put_line("The right hand side vector :"); put_line(b);
      rnk := Rank(A,1.0e-8);
      put("The rank of the Hermite matrix : "); put(rnk,1);
      if bigdim > rnk then
        put(" < "); put(bigdim,1); put_line("  rank deficient.");
      else
        put(" = "); put(bigdim,1); put_line("  okay.");
      end if;
    end if;
    DoblDobl_Complex_Linear_Solvers.lufac(A,bigdim,ipvt,info);
    DoblDobl_Complex_Linear_Solvers.lusolve(A,bigdim,ipvt,b);
    if verbose then
      put_line("The solution vector :"); put_line(b);
    end if;
    res.deg := deg;
    for k in 0..deg loop
      for i in 1..dim loop
        wrk(i) := b(k*dim+i);
      end loop;
      res.cff(k) := new DoblDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    return res;
  end Hermite_Interpolate;

  function Hermite_Laurent_Matrix
             ( m : DoblDobl_Complex_VecMats.VecMat )
             return DoblDobl_Complex_Matrices.Matrix is

    lmt : DoblDobl_Complex_Matrices.Link_to_Matrix := m(0);
    nbr : constant integer32 := lmt'last(1);
    nbc : constant integer32 := lmt'last(2);
    nrows : constant integer32 := nbr*(2*m'last+1);
    ncols : constant integer32 := nbc*(2*m'last+1);
    res : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    rowidx : constant integer32 := nbr*(m'last+1); 
    colidx : constant integer32 := nbc*(m'last+1); 

  begin
    for i in 1..nrows loop
      for j in 1..ncols loop
        res(i,j) := Create(integer32(0));
      end loop;
    end loop;
    for k in m'range loop -- fill up to rowidx and colidx
      lmt := m(k);
      for L in m'range loop
        for i in lmt'range(1) loop
          for j in lmt'range(2) loop
            res((k+L)*nbr+i,L*nbc+j) := lmt(i,j);
          end loop;
        end loop;
      end loop;
    end loop;
    for k in 0..m'last-1 loop
      lmt := m(k);
      for L in 0..m'last-1-k loop
        for i in lmt'range(1) loop
          for j in lmt'range(2) loop
            res(rowidx+(k+L)*nbr+i,colidx+L*nbc+j) := lmt(i,j);
          end loop;
        end loop;
      end loop;
    end loop;
    return res;
  end Hermite_Laurent_Matrix;

  function Hermite_Laurent_Vector
             ( v : DoblDobl_Complex_VecVecs.VecVec )
             return DoblDobl_Complex_Vectors.Vector is

    lv : DoblDobl_Complex_Vectors.Link_to_Vector := v(0);
    nvr : constant integer32 := lv'last;
    dim : constant integer32 := nvr*(2*v'last+1);
    res : DoblDobl_Complex_Vectors.Vector(1..dim);
    idx : constant integer32 := nvr*v'last;

  begin
    for k in 1..idx loop
      res(k) := Create(integer32(0));
    end loop;
    for k in v'range loop
      lv := v(k);
      for i in lv'range loop
        res(idx+k*nvr+i) := lv(i);
      end loop;
    end loop;
    return res;
  end Hermite_Laurent_Vector;

  function Hermite_Laurent_Interpolate
             ( mat : DoblDobl_Dense_Matrix_Series.Matrix;
               rhs : DoblDobl_Dense_Vector_Series.Vector;
               verbose : boolean := true )
             return DoblDobl_Dense_Vector_Series.Vector is

    res : DoblDobl_Dense_Vector_Series.Vector;
    deg : constant integer32 := mat.deg;
    nbr : constant integer32 := mat.cff(0)'last(1);
    nbc : constant integer32 := mat.cff(0)'last(2);
    nrows : constant integer32 := nbr*(2*deg+1);
    ncols : constant integer32 := nbc*(2*deg+1);
    A : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols)
      := Hermite_Laurent_Matrix(mat.cff(0..deg));
    b : DoblDobl_Complex_Vectors.Vector(1..nrows)
      := Hermite_Laurent_Vector(rhs.cff(0..deg));
    qraux : DoblDobl_Complex_Vectors.Vector(1..ncols)
          := (1..ncols => Create(integer32(0)));
    jpvt : Standard_Integer_Vectors.Vector(1..ncols) := (1..ncols => 0);
    sol : DoblDobl_Complex_Vectors.Vector(1..ncols);
    rsd,dum,dum2,dum3 : DoblDobl_Complex_Vectors.Vector(1..nrows);
    info,idx : integer32;
    wrk : DoblDobl_Complex_Vectors.Vector(1..nbc);

    use DoblDobl_Complex_QR_Least_Squares;

  begin
    if verbose then
      put_line("The coefficient matrix :"); put(A,3);
      put_line("The right hand side vector :"); put_line(b);
    end if;
    QRD(A,qraux,jpvt,false);
    QRLS(A,nrows,ncols,qraux,b,dum,dum2,sol,rsd,dum3,110,info);
    if verbose then
      put_line("The least squares solution :"); put_line(sol);
    end if;
    res.deg := deg;
    idx := nbc*deg; -- skip the negative degree entries
    for k in 0..res.deg loop
      for i in wrk'range loop
        wrk(i) := sol(idx+k*nbc+i);
      end loop;
      res.cff(k) := new DoblDobl_Complex_Vectors.Vector'(wrk);
    end loop;
    return res;
  end Hermite_Laurent_Interpolate;

end DoblDobl_Interpolating_Series;
