with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with DoblDobl_Complex_QR_Least_Squares;  use DoblDobl_Complex_QR_Least_Squares;
with DoblDobl_Complex_Singular_Values;   use DoblDobl_Complex_Singular_Values;
with DoblDobl_Interpolating_CSeries;
with DoblDobl_Echelon_Forms;             use DoblDobl_Echelon_Forms;

package body DoblDobl_Series_Matrix_Solvers is

  procedure Solve_Lead_by_lufac
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                a0lu : out DoblDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : DoblDobl_Complex_Vectors.Vector(1..dim) := b.cff(0).all;

  begin
    a0lu := lead.all;
    lufac(a0lu,dim,ipvt,info);
    if info = 0 then
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new DoblDobl_Complex_Vectors.Vector'(x0);
    end if;
  end Solve_Lead_by_lufac;

  procedure Solve_Lead_by_lufco
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                a0lu : out DoblDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : DoblDobl_Complex_Vectors.Vector(1..dim) := b.cff(0).all;
    one : constant double_double := create(1.0);

  begin
    a0lu := lead.all;
    lufco(a0lu,dim,ipvt,rcond);
    if one + rcond /= one then
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new DoblDobl_Complex_Vectors.Vector'(x0);
    end if;
  end Solve_Lead_by_lufco;

  procedure Solve_Lead_by_QRLS
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                a0qr : out DoblDobl_Complex_Matrices.Matrix;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);
    b0 : constant DoblDobl_Complex_Vectors.Vector(1..nrows) := b.cff(0).all;
    x0 : DoblDobl_Complex_Vectors.Vector(1..ncols);
    rsd,dum,dum2,dum3 : DoblDobl_Complex_Vectors.Vector(1..nrows);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    zero : constant double_double := create(0.0);

  begin
    a0qr := lead.all;
    qraux := (qraux'range => Create(zero));
    ipvt := (ipvt'range => 0);
    QRD(a0qr,qraux,ipvt,false);
    wrk := a0qr;
    QRLS(wrk,nrows,ncols,qraux,b0,dum2,dum3,x0,rsd,dum,110,info);
    x.cff(0) := new DoblDobl_Complex_Vectors.Vector'(x0);
  end Solve_Lead_by_QRLS;

  procedure Solve_Lead_by_SVD
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                S : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    n : constant integer32 := lead'last(1);
    p : constant integer32 := lead'last(2);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..n,1..p) := lead.all;
    e : DoblDobl_Complex_Vectors.Vector(1..p);
    job : constant integer32 := 11;
    b0 : constant DoblDobl_Complex_Vectors.Vector(1..n) := b.cff(0).all;
    x0 : DoblDobl_Complex_Vectors.Vector(1..p);

  begin
    SVD(wrk,n,p,S,e,U,V,job,info);
    rcond := Inverse_Condition_Number(S);
    x0 := Solve(U,V,S,b0);
    x.cff(0) := new DoblDobl_Complex_Vectors.Vector'(x0);
  end Solve_Lead_by_SVD;

  procedure Solve_Next_by_lusolve
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                a0lu : in DoblDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                idx : in integer32;
                x : in out DoblDobl_Complex_Vector_Series.Vector ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

    Aidx : DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    dim : constant integer32 := Aidx'last(1);
    wA : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim) := Aidx.all;
    wb : DoblDobl_Complex_Vectors.Vector(1..dim) := b.cff(idx).all;
    wx : DoblDobl_Complex_Vectors.Vector(1..dim) := x.cff(0).all;

  begin
    wb := wb - wA*wx;
    for k in 1..(idx-1) loop
      Aidx := A.cff(idx-k);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    lusolve(a0lu,dim,ipvt,wb);
    x.cff(idx) := new DoblDobl_Complex_Vectors.Vector'(wb);
  end Solve_Next_by_lusolve;

  procedure Solve_Next_by_QRLS
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                a0qr : in DoblDobl_Complex_Matrices.Matrix;
                qraux : in DoblDobl_Complex_Vectors.Vector;
                idx : in integer32; info : out integer32;
                x : in out DoblDobl_Complex_Vector_Series.Vector ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

    Aidx : DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    nrows : constant integer32 := Aidx'last(1);
    ncols : constant integer32 := Aidx'last(2);
    wA : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := Aidx.all;
    wb : DoblDobl_Complex_Vectors.Vector(1..nrows) := b.cff(idx).all;
    wx : DoblDobl_Complex_Vectors.Vector(1..ncols) := x.cff(0).all;
    rsd,dum,dum2,dum3 : DoblDobl_Complex_Vectors.Vector(1..nrows);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := a0qr;

  begin
    wb := wb - wA*wx;
    for k in 1..(idx-1) loop
      Aidx := A.cff(idx-k);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    QRLS(wrk,nrows,ncols,qraux,wb,dum2,dum3,wx,rsd,dum,110,info);
    x.cff(idx) := new DoblDobl_Complex_Vectors.Vector'(wx);
  end Solve_Next_by_QRLS;

  procedure Solve_Next_by_SVD
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                S : in DoblDobl_Complex_Vectors.Vector;
                U,V : in DoblDobl_Complex_Matrices.Matrix;
                idx : in integer32;
                x : in out DoblDobl_Complex_Vector_Series.Vector ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

    Aidx : DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    nrows : constant integer32 := Aidx'last(1);
    ncols : constant integer32 := Aidx'last(2);
    wA : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := Aidx.all;
    wb : DoblDobl_Complex_Vectors.Vector(1..nrows) := b.cff(idx).all;
    wx : DoblDobl_Complex_Vectors.Vector(1..ncols) := x.cff(0).all;

  begin
    wb := wb - wA*wx;
    for k in 1..(idx-1) loop
      Aidx := A.cff(idx-k);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    wx := Solve(U,V,S,wb);
    x.cff(idx) := new DoblDobl_Complex_Vectors.Vector'(wx);
  end Solve_Next_by_SVD;

  procedure Solve_by_lufac
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                info : out integer32;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
   
  begin
    Solve_Lead_by_lufac(A,b,lwrk,ipvt,info,x);
    if info = 0 then
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,k,x);
      end loop;
    end if;
  end Solve_by_lufac;

  procedure Solve_by_lufco
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                rcond : out double_double;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    one : constant double_double := Create(1.0);
   
  begin
    Solve_Lead_by_lufco(A,b,lwrk,ipvt,rcond,x);
    if one + rcond /= one then
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,k,x);
      end loop;
    end if;
  end Solve_by_lufco;

  procedure Solve_by_QRLS
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                info : out integer32;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    nrows : constant integer32 := A.cff(0)'last(1);
    ncols : constant integer32 := A.cff(0)'last(2);
    lwrk : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    qraux : DoblDobl_Complex_Vectors.Vector(1..ncols);
    ipvt : Standard_Integer_Vectors.Vector(1..ncols);

  begin
    Solve_Lead_by_QRLS(A,b,lwrk,qraux,ipvt,info,x);
    if info = 0 then
      for k in 1..b.deg loop
        Solve_Next_by_QRLS(A,b,lwrk,qraux,k,info,x);
      end loop;
    end if;
  end Solve_by_QRLS;

  procedure Solve_by_SVD
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                info : out integer32; rcond : out double_double;
                x : out DoblDobl_Complex_Vector_Series.Vector ) is

    nrows : constant integer32 := A.cff(0)'last(1);
    ncols : constant integer32 := A.cff(0)'last(2);
    mm : constant integer32
       := DoblDobl_Complex_Singular_Values.Min0(nrows+1,ncols);
    S : DoblDobl_Complex_Vectors.Vector(1..mm);
    U : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..nrows);
    V : DoblDobl_Complex_Matrices.Matrix(1..ncols,1..ncols);
    one : constant double_double := create(1.0);

  begin
    Solve_Lead_by_SVD(A,b,S,U,V,info,rcond,x);
    if one + rcond /= one then
      for k in 1..b.deg loop
        Solve_Next_by_SVD(A,b,S,U,V,k,x);
      end loop;
    end if;
  end Solve_by_SVD;

  procedure Echelon_Solve
              ( A : in DoblDobl_Complex_Matrix_Series.Matrix;
                b : in DoblDobl_Complex_Vector_Series.Vector;
                det : out Complex_Number;
                xp : out DoblDobl_Complex_Vector_Series.Vector;
                xn : out DoblDobl_Complex_Vector_Series.Vector ) is

    deg : constant integer32 := A.deg;
    nbr : constant integer32 := A.cff(0)'last(1);
    nbc : constant integer32 := A.cff(0)'last(2);
    nrows : constant integer32 := nbr*(2*deg+1);
    ncols : constant integer32 := nbc*(2*deg+1);
    hlm : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols)
        := DoblDobl_Interpolating_CSeries.Hermite_Laurent_Matrix(A.cff(0..deg));
    x : DoblDobl_Complex_Vectors.Vector(1..ncols);
    rhs : constant DoblDobl_Complex_Vectors.Vector(1..nrows)
        := DoblDobl_Interpolating_CSeries.Hermite_Laurent_Vector(b.cff(0..deg));
    U : DoblDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    row_ipvt : Standard_Integer_Vectors.Vector(1..nrows);
    col_ipvt,pivots : Standard_Integer_Vectors.Vector(1..ncols);
    wx : DoblDobl_Complex_Vectors.Vector(1..nbc);
    startidx : constant integer32 := nbc*deg; -- space for Laurent portion
    one : constant double_double := create(1.0);
    absdet : double_double;

  begin
    Lower_Triangular_Echelon_Form(nbc,hlm,U,row_ipvt,col_ipvt,pivots,false);
    Solve_with_Echelon_Form(hlm,rhs,x);
    Multiply_and_Permute(x,U,pivots);
    det := Determinant(hlm);
    absdet := AbsVal(det);
   -- xp.deg := deg;
    for i in 0..deg loop
      for j in 1..nbc loop
        wx(j) := x(startidx+i*nbc+j);
      end loop;
      xp.cff(i) := new DoblDobl_Complex_Vectors.Vector'(wx);
    end loop;
   -- if absdet + 1.0 /= 1.0 then
   --   xn.deg := -1; -- no Laurent portion
   -- else
    if absdet + one = one then
      for i in 1..deg loop
        for j in 1..nbc loop
          wx(j) := x(startidx-i*nbc+j);
        end loop;
        xn.cff(i) := new DoblDobl_Complex_Vectors.Vector'(wx);
      end loop;
    end if;
  end Echelon_Solve;

-- ON FLATTENED DATA STRUCTURES :

  procedure Solve_Lead_by_lufac
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    a0lu : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A(0);
    dim : constant integer32 := a0lu'last(1);
 
  begin
    lufac(a0lu.all,dim,ipvt,info);
    if info = 0
     then lusolve(a0lu.all,dim,ipvt,b(0).all);
    end if;
  end Solve_Lead_by_lufac;

  procedure Solve_Lead_by_lufco
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double ) is

    a0lu : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A(0);
    dim : constant integer32 := a0lu'last(1);
    one : constant double_double := create(1.0);
 
  begin
    lufco(a0lu.all,dim,ipvt,rcond);
    if one + rcond /= one
     then lusolve(a0lu.all,dim,ipvt,b(0).all);
    end if;
  end Solve_Lead_by_lufco;

  procedure Solve_Lead_by_QRLS
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x0 : in DoblDobl_Complex_Vectors.Link_to_Vector;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);
    zero : constant Complex_Number := create(integer32(0));

  begin
    qraux := (qraux'range => zero);
    ipvt := (ipvt'range => 0);
    QRD(lead.all,qraux,ipvt,false);
    w1 := b(0).all;
    QRLS(lead.all,nrows,ncols,qraux,w1,w2,w3,x0.all,w4,w5,110,info);
  end Solve_Lead_by_QRLS;

  procedure Solve_Lead_by_SVD
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x0 : in DoblDobl_Complex_Vectors.Link_to_Vector;
                S : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A(0);
    n : constant integer32 := lead'last(1);
    p : constant integer32 := lead'last(2);
    job : constant integer32 := 11;

  begin
    SVD(lead.all,n,p,S,ewrk.all,U,V,job,info,wrkv.all);
    rcond := Inverse_Condition_Number(S);
    x0.all := Solve(U,V,S,b(0).all);
  end Solve_Lead_by_SVD;

  procedure Matrix_Vector_Multiply
              ( A : in DoblDobl_Complex_Matrices.Link_to_Matrix;
                x,y : in DoblDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    for i in y'range loop
      y(i) := A(i,A'first(2))*x(x'first);
      for j in x'first+1..x'last loop
        y(i) := y(i) + A(i,j)*x(j);
      end loop;
    end loop;
  end Matrix_Vector_Multiply;

  procedure Subtract ( x,y : in DoblDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    for i in x'range loop
      x(i) := x(i) - y(i);
    end loop;
  end Subtract;

  procedure Solve_Next_by_lusolve
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                idx : in integer32;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

    dim : constant integer32 := wrk'last;

  begin
    Matrix_Vector_Multiply(A(idx),b(0),wrk); -- wrk = A(idx)*b(0)
    Subtract(b(idx),wrk);                    -- b(idx) := b(idx) - wrk
    for k in 1..(idx-1) loop
      Matrix_Vector_Multiply(A(idx-k),b(k),wrk); -- wrk = A(idx-k)*b(k)
      Subtract(b(idx),wrk);                      -- b(idx) := b(idx) - wrk
    end loop;
    lusolve(A(0).all,dim,ipvt,b(idx).all);
  end Solve_Next_by_lusolve;

  procedure Solve_Next_by_QRLS
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : in DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                idx : in integer32; info : out integer32;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    lead : constant DoblDobl_Complex_Matrices.Link_to_Matrix := A(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);

  begin
    Matrix_Vector_Multiply(A(idx),x(0),wrk);     -- wrk = A(idx)*x(0)
    Subtract(b(idx),wrk);                        -- b(idx) := b(idx) - wrk
    for k in 1..(idx-1) loop
      Matrix_Vector_Multiply(A(idx-k),x(k),wrk); -- wrk = A(idx-k)*x(k)
      Subtract(b(idx),wrk);                      -- b(idx) := b(idx) - wrk
    end loop;
    w1 := b(idx).all;
    QRLS(lead.all,nrows,ncols,qraux,w1,w2,w3,x(idx).all,w4,w5,110,info);
  end Solve_Next_by_QRLS;

  procedure Solve_Next_by_SVD
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : in DoblDobl_Complex_Vectors.Vector;
                U,V : in DoblDobl_Complex_Matrices.Matrix;
                idx : in integer32;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Matrix_Vector_Multiply(A(idx),x(0),wrk);     -- wrk = A(idx)*x(0)
    Subtract(b(idx),wrk);                        -- b(idx) := b(idx) - wrk
    for k in 1..(idx-1) loop
      Matrix_Vector_Multiply(A(idx-k),x(k),wrk); -- wrk = A(idx-k)*x(k)
      Subtract(b(idx),wrk);                      -- b(idx) := b(idx) - wrk
    end loop;
    x(idx).all := Solve(U,V,S,b(idx).all);
  end Solve_Next_by_SVD;

  procedure Solve_by_lufac
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0 then
      for k in 1..b'last loop
        Solve_Next_by_lusolve(A,b,ipvt,k,wrk);
      end loop;
    end if;
  end Solve_by_lufac;

  procedure Solve_by_lufco
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    one : constant double_double := create(1.0);

  begin
    Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if one + rcond /= one then
      for k in 1..b'last loop
        Solve_Next_by_lusolve(A,b,ipvt,k,wrk);
      end loop;
    end if;
  end Solve_by_lufco;

  procedure Solve_by_QRLS
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                qraux : out DoblDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in DoblDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Solve_Lead_by_QRLS(A,b,x(0),qraux,w1,w2,w3,w4,w5,ipvt,info);
    if info = 0 then
      for k in 1..b'last loop
        Solve_Next_by_QRLS(A,b,x,qraux,w1,w2,w3,w4,w5,k,info,wrk);
      end loop;
    end if;
  end Solve_by_QRLS;

  procedure Solve_by_SVD
              ( A : in DoblDobl_Complex_VecMats.VecMat;
                b : in DoblDobl_Complex_VecVecs.VecVec;
                x : in DoblDobl_Complex_VecVecs.VecVec;
                S : out DoblDobl_Complex_Vectors.Vector;
                U,V : out DoblDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out double_double;
                ewrk : in DoblDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    one : constant double_double := create(1.0);

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv);
    if one + rcond /= one then
      for k in 1..b'last loop
        Solve_Next_by_SVD(A,b,x,S,U,V,k,wrkv);
      end loop;
    end if;
  end Solve_by_SVD;

end DoblDobl_Series_Matrix_Solvers;
