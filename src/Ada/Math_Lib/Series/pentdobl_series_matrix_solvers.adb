with PentDobl_Complex_Linear_Solvers;    use PentDobl_Complex_Linear_Solvers;
with PentDobl_Complex_QR_Least_Squares;  use PentDobl_Complex_QR_Least_Squares;
with PentDobl_Complex_Singular_Values;   use PentDobl_Complex_Singular_Values;
with PentDobl_Interpolating_CSeries;
with PentDobl_Echelon_Forms;             use PentDobl_Echelon_Forms;

package body PentDobl_Series_Matrix_Solvers is

  procedure Solve_Lead_by_lufac
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                a0lu : out PentDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : PentDobl_Complex_Vectors.Vector(1..dim) := b.cff(0).all;

  begin
    a0lu := lead.all;
    lufac(a0lu,dim,ipvt,info);
    if info = 0 then
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new PentDobl_Complex_Vectors.Vector'(x0);
    end if;
  end Solve_Lead_by_lufac;

  procedure Solve_Lead_by_lufco
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                a0lu : out PentDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out penta_double;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : PentDobl_Complex_Vectors.Vector(1..dim) := b.cff(0).all;
    one : constant penta_double := create(1.0);

  begin
    a0lu := lead.all;
    lufco(a0lu,dim,ipvt,rcond);
    if one + rcond /= one then
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new PentDobl_Complex_Vectors.Vector'(x0);
    end if;
  end Solve_Lead_by_lufco;

  procedure Solve_Lead_by_QRLS
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                a0qr : out PentDobl_Complex_Matrices.Matrix;
                qraux : out PentDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);
    b0 : constant PentDobl_Complex_Vectors.Vector(1..nrows) := b.cff(0).all;
    x0 : PentDobl_Complex_Vectors.Vector(1..ncols);
    rsd,dum,dum2,dum3 : PentDobl_Complex_Vectors.Vector(1..nrows);
    wrk : PentDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    zero : constant penta_double := create(0.0);

  begin
    a0qr := lead.all;
    qraux := (qraux'range => Create(zero));
    ipvt := (ipvt'range => 0);
    QRD(a0qr,qraux,ipvt,false);
    wrk := a0qr;
    QRLS(wrk,nrows,ncols,qraux,b0,dum2,dum3,x0,rsd,dum,110,info);
    x.cff(0) := new PentDobl_Complex_Vectors.Vector'(x0);
  end Solve_Lead_by_QRLS;

  procedure Solve_Lead_by_SVD
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                S : out PentDobl_Complex_Vectors.Vector;
                U,V : out PentDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out penta_double;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    n : constant integer32 := lead'last(1);
    p : constant integer32 := lead'last(2);
    wrk : PentDobl_Complex_Matrices.Matrix(1..n,1..p) := lead.all;
    e : PentDobl_Complex_Vectors.Vector(1..p);
    job : constant integer32 := 11;
    b0 : constant PentDobl_Complex_Vectors.Vector(1..n) := b.cff(0).all;
    x0 : PentDobl_Complex_Vectors.Vector(1..p);

  begin
    SVD(wrk,n,p,S,e,U,V,job,info);
    rcond := Inverse_Condition_Number(S);
    x0 := Solve(U,V,S,b0);
    x.cff(0) := new PentDobl_Complex_Vectors.Vector'(x0);
  end Solve_Lead_by_SVD;

  procedure Solve_Next_by_lusolve
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                a0lu : in PentDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                idx : in integer32;
                x : in out PentDobl_Complex_Vector_Series.Vector ) is

    use PentDobl_Complex_Vectors;
    use PentDobl_Complex_Matrices;

    Aidx : PentDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    dim : constant integer32 := Aidx'last(1);
    wA : PentDobl_Complex_Matrices.Matrix(1..dim,1..dim) := Aidx.all;
    wb : PentDobl_Complex_Vectors.Vector(1..dim) := b.cff(idx).all;
    wx : PentDobl_Complex_Vectors.Vector(1..dim) := x.cff(0).all;

  begin
    wb := wb - wA*wx;
    for k in 1..(idx-1) loop
      Aidx := A.cff(idx-k);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    lusolve(a0lu,dim,ipvt,wb);
    x.cff(idx) := new PentDobl_Complex_Vectors.Vector'(wb);
  end Solve_Next_by_lusolve;

  procedure Solve_Next_by_QRLS
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                a0qr : in PentDobl_Complex_Matrices.Matrix;
                qraux : in PentDobl_Complex_Vectors.Vector;
                idx : in integer32; info : out integer32;
                x : in out PentDobl_Complex_Vector_Series.Vector ) is

    use PentDobl_Complex_Vectors;
    use PentDobl_Complex_Matrices;

    Aidx : PentDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    nrows : constant integer32 := Aidx'last(1);
    ncols : constant integer32 := Aidx'last(2);
    wA : PentDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := Aidx.all;
    wb : PentDobl_Complex_Vectors.Vector(1..nrows) := b.cff(idx).all;
    wx : PentDobl_Complex_Vectors.Vector(1..ncols) := x.cff(0).all;
    rsd,dum,dum2,dum3 : PentDobl_Complex_Vectors.Vector(1..nrows);
    wrk : PentDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := a0qr;

  begin
    wb := wb - wA*wx;
    for k in 1..(idx-1) loop
      Aidx := A.cff(idx-k);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    QRLS(wrk,nrows,ncols,qraux,wb,dum2,dum3,wx,rsd,dum,110,info);
    x.cff(idx) := new PentDobl_Complex_Vectors.Vector'(wx);
  end Solve_Next_by_QRLS;

  procedure Solve_Next_by_SVD
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                S : in PentDobl_Complex_Vectors.Vector;
                U,V : in PentDobl_Complex_Matrices.Matrix;
                idx : in integer32;
                x : in out PentDobl_Complex_Vector_Series.Vector ) is

    use PentDobl_Complex_Vectors;
    use PentDobl_Complex_Matrices;

    Aidx : PentDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    nrows : constant integer32 := Aidx'last(1);
    ncols : constant integer32 := Aidx'last(2);
    wA : PentDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := Aidx.all;
    wb : PentDobl_Complex_Vectors.Vector(1..nrows) := b.cff(idx).all;
    wx : PentDobl_Complex_Vectors.Vector(1..ncols) := x.cff(0).all;

  begin
    wb := wb - wA*wx;
    for k in 1..(idx-1) loop
      Aidx := A.cff(idx-k);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    wx := Solve(U,V,S,wb);
    x.cff(idx) := new PentDobl_Complex_Vectors.Vector'(wx);
  end Solve_Next_by_SVD;

  procedure Solve_by_lufac
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                info : out integer32;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : PentDobl_Complex_Matrices.Matrix(1..dim,1..dim);
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
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                rcond : out penta_double;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : PentDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    one : constant penta_double := Create(1.0);
   
  begin
    Solve_Lead_by_lufco(A,b,lwrk,ipvt,rcond,x);
    if one + rcond /= one then
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,k,x);
      end loop;
    end if;
  end Solve_by_lufco;

  procedure Solve_by_QRLS
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                info : out integer32;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    nrows : constant integer32 := A.cff(0)'last(1);
    ncols : constant integer32 := A.cff(0)'last(2);
    lwrk : PentDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    qraux : PentDobl_Complex_Vectors.Vector(1..ncols);
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
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                info : out integer32; rcond : out penta_double;
                x : out PentDobl_Complex_Vector_Series.Vector ) is

    nrows : constant integer32 := A.cff(0)'last(1);
    ncols : constant integer32 := A.cff(0)'last(2);
    mm : constant integer32
       := PentDobl_Complex_Singular_Values.Min0(nrows+1,ncols);
    S : PentDobl_Complex_Vectors.Vector(1..mm);
    U : PentDobl_Complex_Matrices.Matrix(1..nrows,1..nrows);
    V : PentDobl_Complex_Matrices.Matrix(1..ncols,1..ncols);
    one : constant penta_double := create(1.0);

  begin
    Solve_Lead_by_SVD(A,b,S,U,V,info,rcond,x);
    if one + rcond /= one then
      for k in 1..b.deg loop
        Solve_Next_by_SVD(A,b,S,U,V,k,x);
      end loop;
    end if;
  end Solve_by_SVD;

  procedure Echelon_Solve
              ( A : in PentDobl_Complex_Matrix_Series.Matrix;
                b : in PentDobl_Complex_Vector_Series.Vector;
                det : out Complex_Number;
                xp : out PentDobl_Complex_Vector_Series.Vector;
                xn : out PentDobl_Complex_Vector_Series.Vector ) is

    deg : constant integer32 := A.deg;
    nbr : constant integer32 := A.cff(0)'last(1);
    nbc : constant integer32 := A.cff(0)'last(2);
    nrows : constant integer32 := nbr*(2*deg+1);
    ncols : constant integer32 := nbc*(2*deg+1);
    hlm : PentDobl_Complex_Matrices.Matrix(1..nrows,1..ncols)
        := PentDobl_Interpolating_CSeries.Hermite_Laurent_Matrix(A.cff(0..deg));
    x : PentDobl_Complex_Vectors.Vector(1..ncols);
    rhs : constant PentDobl_Complex_Vectors.Vector(1..nrows)
        := PentDobl_Interpolating_CSeries.Hermite_Laurent_Vector(b.cff(0..deg));
    U : PentDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    row_ipvt : Standard_Integer_Vectors.Vector(1..nrows);
    col_ipvt,pivots : Standard_Integer_Vectors.Vector(1..ncols);
    wx : PentDobl_Complex_Vectors.Vector(1..nbc);
    startidx : constant integer32 := nbc*deg; -- space for Laurent portion
    one : constant penta_double := create(1.0);
    absdet : penta_double;

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
      xp.cff(i) := new PentDobl_Complex_Vectors.Vector'(wx);
    end loop;
   -- if absdet + 1.0 /= 1.0 then
   --   xn.deg := -1; -- no Laurent portion
   -- else
    if absdet + one = one then
      for i in 1..deg loop
        for j in 1..nbc loop
          wx(j) := x(startidx-i*nbc+j);
        end loop;
        xn.cff(i) := new PentDobl_Complex_Vectors.Vector'(wx);
      end loop;
    end if;
  end Echelon_Solve;

-- ON FLATTENED DATA STRUCTURES :

  procedure Solve_Lead_by_lufac
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    a0lu : constant PentDobl_Complex_Matrices.Link_to_Matrix := A(0);
    dim : constant integer32 := a0lu'last(1);
 
  begin
    lufac(a0lu.all,dim,ipvt,info);
    if info = 0
     then lusolve(a0lu.all,dim,ipvt,b(0).all);
    end if;
  end Solve_Lead_by_lufac;

  procedure Solve_Lead_by_lufco
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out penta_double ) is

    a0lu : constant PentDobl_Complex_Matrices.Link_to_Matrix := A(0);
    dim : constant integer32 := a0lu'last(1);
    one : constant penta_double := create(1.0);
 
  begin
    lufco(a0lu.all,dim,ipvt,rcond);
    if one + rcond /= one
     then lusolve(a0lu.all,dim,ipvt,b(0).all);
    end if;
  end Solve_Lead_by_lufco;

  procedure Solve_Lead_by_QRLS
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x0 : in PentDobl_Complex_Vectors.Link_to_Vector;
                qraux : out PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A(0);
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
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x0 : in PentDobl_Complex_Vectors.Link_to_Vector;
                S : out PentDobl_Complex_Vectors.Vector;
                U,V : out PentDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out penta_double;
                ewrk : in PentDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A(0);
    n : constant integer32 := lead'last(1);
    p : constant integer32 := lead'last(2);
    job : constant integer32 := 11;

  begin
    SVD(lead.all,n,p,S,ewrk.all,U,V,job,info,wrkv.all);
    rcond := Inverse_Condition_Number(S);
    x0.all := Solve(U,V,S,b(0).all);
  end Solve_Lead_by_SVD;

  procedure Matrix_Vector_Multiply
              ( A : in PentDobl_Complex_Matrices.Link_to_Matrix;
                x,y : in PentDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    for i in y'range loop
      y(i) := A(i,A'first(2))*x(x'first);
      for j in x'first+1..x'last loop
        y(i) := y(i) + A(i,j)*x(j);
      end loop;
    end loop;
  end Matrix_Vector_Multiply;

  procedure Subtract ( x,y : in PentDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    for i in x'range loop
      x(i) := x(i) - y(i);
    end loop;
  end Subtract;

  procedure Solve_Next_by_lusolve
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                idx : in integer32;
                wrk : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    use PentDobl_Complex_Vectors;
    use PentDobl_Complex_Matrices;

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
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : in PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_Vectors.Vector;
                idx : in integer32; info : out integer32;
                wrk : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    lead : constant PentDobl_Complex_Matrices.Link_to_Matrix := A(0);
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
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : in PentDobl_Complex_Vectors.Vector;
                U,V : in PentDobl_Complex_Matrices.Matrix;
                idx : in integer32;
                wrk : in PentDobl_Complex_Vectors.Link_to_Vector ) is
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
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in PentDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Solve_Lead_by_lufac(A,b,ipvt,info);
    if info = 0 then
      for k in 1..b'last loop
        Solve_Next_by_lusolve(A,b,ipvt,k,wrk);
      end loop;
    end if;
  end Solve_by_lufac;

  procedure Solve_by_lufco
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out penta_double;
                wrk : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    one : constant penta_double := create(1.0);

  begin
    Solve_Lead_by_lufco(A,b,ipvt,rcond);
    if one + rcond /= one then
      for k in 1..b'last loop
        Solve_Next_by_lusolve(A,b,ipvt,k,wrk);
      end loop;
    end if;
  end Solve_by_lufco;

  procedure Solve_by_QRLS
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                qraux : out PentDobl_Complex_Vectors.Vector;
                w1,w2,w3,w4,w5 : in out PentDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                wrk : in PentDobl_Complex_Vectors.Link_to_Vector ) is
  begin
    Solve_Lead_by_QRLS(A,b,x(0),qraux,w1,w2,w3,w4,w5,ipvt,info);
    if info = 0 then
      for k in 1..b'last loop
        Solve_Next_by_QRLS(A,b,x,qraux,w1,w2,w3,w4,w5,k,info,wrk);
      end loop;
    end if;
  end Solve_by_QRLS;

  procedure Solve_by_SVD
              ( A : in PentDobl_Complex_VecMats.VecMat;
                b : in PentDobl_Complex_VecVecs.VecVec;
                x : in PentDobl_Complex_VecVecs.VecVec;
                S : out PentDobl_Complex_Vectors.Vector;
                U,V : out PentDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out penta_double;
                ewrk : in PentDobl_Complex_Vectors.Link_to_Vector;
                wrkv : in PentDobl_Complex_Vectors.Link_to_Vector ) is

    one : constant penta_double := create(1.0);

  begin
    Solve_Lead_by_SVD(A,b,x(0),S,U,V,info,rcond,ewrk,wrkv);
    if one + rcond /= one then
      for k in 1..b'last loop
        Solve_Next_by_SVD(A,b,x,S,U,V,k,wrkv);
      end loop;
    end if;
  end Solve_by_SVD;

end PentDobl_Series_Matrix_Solvers;
