with text_io;                             use text_io;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with QuadDobl_Random_Numbers;
with QuadDobl_Random_Vectors;
with Standard_Integer_Vectors;
with Quad_Double_Vectors_io;            use Quad_Double_Vectors_io;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs_io;         use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Singular_Values;

package body QuadDobl_Interpolating_Series is

  function Eval ( v : QuadDobl_Dense_Vector_Series.Vector;
                  t : Complex_Number )
                return QuadDobl_Complex_Vectors.Vector is

    vec : QuadDobl_Complex_Vectors.Link_to_Vector := v.cff(0);
    res : QuadDobl_Complex_Vectors.Vector(vec'range) := vec.all;
    pwt : Complex_Number := Create(integer32(1));
 
    use QuadDobl_Complex_Vectors;

  begin
    for k in 1..v.deg loop
      pwt := pwt*t;
      res := res + pwt*v.cff(k).all;
    end loop;
    return res;
  end Eval;

  function Eval ( m : QuadDobl_Dense_Matrix_Series.Matrix;
                  t : Complex_Number )
                return QuadDobl_Complex_Matrices.Matrix is

    mat : QuadDobl_Complex_Matrices.Link_to_Matrix := m.cff(0);
    res : QuadDobl_Complex_Matrices.Matrix(mat'range(1),mat'range(2))
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

  function Full_Rank
             ( m : QuadDobl_Dense_Matrix_Series.Matrix;
               d : integer32; verbose : boolean := true ) return boolean is

    t : constant Complex_Number := QuadDobl_Random_Numbers.Random1;
    mat : QuadDobl_Complex_Matrices.Link_to_Matrix := m.cff(0);
    val : QuadDobl_Complex_Matrices.Matrix(mat'range(1),mat'range(2))
        := mat.all;
    pwt : Complex_Number := Create(integer32(1));
    n : constant integer32 := val'last(1);
    p : constant integer32 := val'last(2);
    mm : constant integer32 := QuadDobl_Complex_Singular_Values.Min0(n+1,p);
    s : QuadDobl_Complex_Vectors.Vector(1..mm);
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    u : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : QuadDobl_Complex_Matrices.Matrix(1..p,1..p);
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
    QuadDobl_Complex_Singular_Values.SVD(val,n,p,s,e,u,v,job,info);
    rnk := QuadDobl_Complex_Singular_Values.Rank(s,tol);
    if verbose then
      put_line("The singular values : "); put_line(s);
      put("The numerical rank : "); put(rnk,1); new_line;
    end if;
    return (rnk = n);
  end Full_Rank;

  function Full_Rank
             ( m : QuadDobl_Dense_Matrix_Series.Matrix;
               verbose : boolean := true ) return integer32 is
  begin
    for d in 0..m.deg loop
      if Full_Rank(m,d)
       then return d;
      end if;
    end loop;
    return -1;
  end Full_Rank;

  function Sample ( v : QuadDobl_Dense_Vector_Series.Vector;
                    t : QuadDobl_Complex_Vectors.Vector )
                  return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(t'range);
    lv : QuadDobl_Complex_Vectors.Link_to_Vector := v.cff(0);

  begin
    for i in t'range loop
      declare
        e : constant QuadDobl_Complex_Vectors.Vector(lv'range)
          := Eval(v,t(i));
      begin
        res(i) := new QuadDobl_Complex_Vectors.Vector'(e);
      end;
    end loop;
    return res;
  end Sample;

  function Sample ( m : QuadDobl_Dense_Matrix_Series.Matrix;
                    t : QuadDobl_Complex_Vectors.Vector )
                  return QuadDobl_Complex_VecMats.VecMat is

    res : QuadDobl_Complex_VecMats.VecMat(t'range);
    lm : QuadDobl_Complex_Matrices.Link_to_Matrix := m.cff(0);

  begin
    for i in t'range loop
      declare
        e : constant QuadDobl_Complex_Matrices.Matrix(lm'range(1),lm'range(2))
          := Eval(m,t(i));
      begin
        res(i) := new QuadDobl_Complex_Matrices.Matrix'(e);
      end;
    end loop;
    return res;
  end Sample;

  function Solve_Linear_Systems
             ( m : QuadDobl_Complex_VecMats.VecMat;
               v : QuadDobl_Complex_VecVecs.VecVec )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(m'range);
    dim : constant integer32 := v(0)'last;
    info : integer32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    sol : QuadDobl_Complex_Vectors.Vector(1..dim);
    wrk : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in v'range loop
      sol := v(i).all;
      wrk := m(i).all;
      QuadDobl_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
      QuadDobl_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,sol);
      res(i) := new QuadDobl_Complex_Vectors.Vector'(sol);
    end loop;
    return res;
  end Solve_Linear_Systems;

  function Residuals
             ( m : QuadDobl_Complex_VecMats.VecMat;
               v,x : QuadDobl_Complex_VecVecs.VecVec )
             return Quad_Double_Vectors.Vector is

    res : Quad_Double_Vectors.Vector(v'range);
    dim : constant integer32 := v(0)'last;
    lm : QuadDobl_Complex_Matrices.Link_to_Matrix;
    lv,lx : QuadDobl_Complex_Vectors.Link_to_Vector;
    wrk : QuadDobl_Complex_Vectors.Vector(1..dim);

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

  begin
    for i in res'range loop
      lm := m(i);
      lv := v(i);
      lx := x(i);
      wrk := lv.all - lm.all*lx.all;
      res(i) := QuadDobl_Complex_Vector_Norms.Max_Norm(wrk);
    end loop;
    return res;
  end Residuals;

  function Transpose ( x : QuadDobl_Complex_VecVecs.VecVec )
                     return QuadDobl_Complex_VecVecs.VecVec is

    deg : constant integer32 := x'last;
    dim : constant integer32 := x(0)'last;
    res : QuadDobl_Complex_VecVecs.VecVec(1..dim);

  begin
    for i in res'range loop
      res(i) := new QuadDobl_Complex_Vectors.Vector(1..deg+1);
    end loop;
    for i in x'range loop
      for j in x(i)'range loop
        res(j)(i+1) := x(i)(j);
      end loop;
    end loop;
    return res;
  end Transpose;

  function Vandermonde_Matrix
             ( t : QuadDobl_Complex_Vectors.Vector )
             return QuadDobl_Complex_Matrices.Matrix is

    dim : constant integer32 := t'last - t'first + 1;
    res : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
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
             ( v : QuadDobl_Complex_Matrices.Matrix;
               f : QuadDobl_Complex_VecVecs.VecVec )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(f'range);
    dim : constant integer32 := v'last(1);
    wrk : QuadDobl_Complex_Matrices.Matrix(v'range(1),v'range(2)) := v;
    sol : QuadDobl_Complex_Vectors.Vector(1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;

  begin
    QuadDobl_Complex_Linear_Solvers.lufac(wrk,dim,ipvt,info);
    for i in f'range loop
      sol := f(i).all;
      QuadDobl_Complex_Linear_Solvers.lusolve(wrk,dim,ipvt,sol);
      res(i) := new QuadDobl_Complex_Vectors.Vector'(sol);
    end loop;
    return res;
  end Solve_Interpolation_Systems;

  function Construct ( x : QuadDobl_Complex_VecVecs.VecVec )
                     return QuadDobl_Dense_Vector_Series.Vector is

    res : QuadDobl_Dense_Vector_Series.Vector;
    dim : constant integer32 := x'last;
    lx0 : constant QuadDobl_Complex_Vectors.Link_to_Vector := x(x'first); 
    deg : constant integer32 := lx0'last - lx0'first;

  begin
    res.deg := deg;
    for i in 0..res.deg loop
      res.cff(i) := new QuadDobl_Complex_Vectors.Vector(1..dim);
    end loop;
    for i in x'range loop
      for j in x(i)'range loop
        res.cff(j-1)(i) := x(i)(j);
      end loop;
    end loop;
    return res;
  end Construct;

  function Interpolate
             ( mat : QuadDobl_Dense_Matrix_Series.Matrix;
               rhs : QuadDobl_Dense_Vector_Series.Vector;
               verbose : boolean := true )
             return QuadDobl_Dense_Vector_Series.Vector is

    res : QuadDobl_Dense_Vector_Series.Vector;
    dim : constant integer32 := rhs.cff(0)'last;
    t : constant QuadDobl_Complex_Vectors.Vector(0..mat.deg)
      := QuadDobl_Random_Vectors.Random_Vector(0,mat.deg);
    m : QuadDobl_Complex_VecMats.VecMat(t'range) := Sample(mat,t);
    v : QuadDobl_Complex_VecVecs.VecVec(t'range) := Sample(rhs,t);
    x : QuadDobl_Complex_VecVecs.VecVec(t'range);
    r : Quad_Double_Vectors.Vector(t'range);
    xt : QuadDobl_Complex_VecVecs.VecVec(1..dim);
    vdm : QuadDobl_Complex_Matrices.Matrix(1..mat.deg+1,1..mat.deg+1);
    cff : QuadDobl_Complex_VecVecs.VecVec(1..dim);

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

end QuadDobl_Interpolating_Series;
