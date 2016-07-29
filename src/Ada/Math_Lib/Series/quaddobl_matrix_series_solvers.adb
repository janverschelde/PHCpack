with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_QR_Least_Squares;  use QuadDobl_Complex_QR_Least_Squares;
with QuadDobl_Complex_Singular_Values;   use QuadDobl_Complex_Singular_Values;

package body QuadDobl_Matrix_Series_Solvers is

  procedure Solve_Lead_by_lufac
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                a0lu : out QuadDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    lead : constant QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : QuadDobl_Complex_Vectors.Vector(1..dim) := b.cff(0).all;

  begin
    a0lu := lead.all;
    lufac(a0lu,dim,ipvt,info);
    if info /= 0 then
      x.deg := -1;
    else
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new QuadDobl_Complex_Vectors.Vector'(x0);
      x.deg := 0;
    end if;
  end Solve_Lead_by_lufac;

  procedure Solve_Lead_by_lufco
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                a0lu : out QuadDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    lead : constant QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : QuadDobl_Complex_Vectors.Vector(1..dim) := b.cff(0).all;
    one : constant quad_double := create(1.0);

  begin
    a0lu := lead.all;
    lufco(a0lu,dim,ipvt,rcond);
    if one + rcond = one then
      x.deg := -1;
    else
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new QuadDobl_Complex_Vectors.Vector'(x0);
      x.deg := 0;
    end if;
  end Solve_Lead_by_lufco;

  procedure Solve_Lead_by_QRLS
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                a0qr : out QuadDobl_Complex_Matrices.Matrix;
                qraux : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    lead : constant QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);
    b0 : constant QuadDobl_Complex_Vectors.Vector(1..nrows) := b.cff(0).all;
    x0 : QuadDobl_Complex_Vectors.Vector(1..ncols);
    rsd,dum,dum2,dum3 : QuadDobl_Complex_Vectors.Vector(1..nrows);
    zero : constant quad_double := create(0.0);
    wrk : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);

  begin
    a0qr := lead.all;
    qraux := (qraux'range => Create(zero));
    ipvt := (ipvt'range => 0);
    QRD(a0qr,qraux,ipvt,false);
    wrk := a0qr;
    QRLS(wrk,nrows,ncols,qraux,b0,dum2,dum3,x0,rsd,dum,110,info);
    x.cff(0) := new QuadDobl_Complex_Vectors.Vector'(x0);
    x.deg := 0;
  end Solve_Lead_by_QRLS;

  procedure Solve_Lead_by_SVD
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                S : out QuadDobl_Complex_Vectors.Vector;
                U,V : out QuadDobl_Complex_Matrices.Matrix;
                info : out integer32; rcond : out quad_double;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    lead : constant QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(0);
    n : constant integer32 := lead'last(1);
    p : constant integer32 := lead'last(2);
    wrk : QuadDobl_Complex_Matrices.Matrix(1..n,1..p) := lead.all;
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    job : constant integer32 := 11;
    b0 : constant QuadDobl_Complex_Vectors.Vector(1..n) := b.cff(0).all;
    x0 : QuadDobl_Complex_Vectors.Vector(1..p);

  begin
    SVD(wrk,n,p,S,e,U,V,job,info);
    rcond := Inverse_Condition_Number(S);
    x0 := Solve(U,V,S,b0);
    x.cff(0) := new QuadDobl_Complex_Vectors.Vector'(x0);
    x.deg := 0;
  end Solve_Lead_by_SVD;

  procedure Solve_Next_by_lusolve
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                a0lu : in QuadDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                x : in out QuadDobl_Dense_Vector_Series.Vector ) is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

    idx : integer32 := x.deg+1;
    Aidx : QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    dim : constant integer32 := Aidx'last(1);
    wA : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim) := Aidx.all;
    wb : QuadDobl_Complex_Vectors.Vector(1..dim) := b.cff(idx).all;
    wx : QuadDobl_Complex_Vectors.Vector(1..dim) := x.cff(0).all;

  begin
    wb := wb - wA*wx;
    for k in 1..x.deg loop
      idx := idx - 1;
      Aidx := A.cff(idx);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    lusolve(a0lu,dim,ipvt,wb);
    x.deg := x.deg + 1;
    x.cff(x.deg) := new QuadDobl_Complex_Vectors.Vector'(wb);
  end Solve_Next_by_lusolve;

  procedure Solve_Next_by_QRLS
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                a0qr : in QuadDobl_Complex_Matrices.Matrix;
                qraux : in QuadDobl_Complex_Vectors.Vector;
                info : out integer32;
                x : in out QuadDobl_Dense_Vector_Series.Vector ) is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

    idx : integer32 := x.deg+1;
    Aidx : QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    nrows : constant integer32 := Aidx'last(1);
    ncols : constant integer32 := Aidx'last(2);
    wA : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := Aidx.all;
    wb : QuadDobl_Complex_Vectors.Vector(1..nrows) := b.cff(idx).all;
    wx : QuadDobl_Complex_Vectors.Vector(1..ncols) := x.cff(0).all;
    rsd,dum,dum2,dum3 : QuadDobl_Complex_Vectors.Vector(1..nrows);
    wrk : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := a0qr;

  begin
    wb := wb - wA*wx;
    for k in 1..x.deg loop
      idx := idx - 1;
      Aidx := A.cff(idx);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    QRLS(wrk,nrows,ncols,qraux,wb,dum2,dum3,wx,rsd,dum,110,info);
    x.deg := x.deg + 1;
    x.cff(x.deg) := new QuadDobl_Complex_Vectors.Vector'(wx);
  end Solve_Next_by_QRLS;

  procedure Solve_Next_by_SVD
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                S : in QuadDobl_Complex_Vectors.Vector;
                U,V : in QuadDobl_Complex_Matrices.Matrix;
                x : in out QuadDobl_Dense_Vector_Series.Vector ) is

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

    idx : integer32 := x.deg+1;
    Aidx : QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    nrows : constant integer32 := Aidx'last(1);
    ncols : constant integer32 := Aidx'last(2);
    wA : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..ncols) := Aidx.all;
    wb : QuadDobl_Complex_Vectors.Vector(1..nrows) := b.cff(idx).all;
    wx : QuadDobl_Complex_Vectors.Vector(1..ncols) := x.cff(0).all;

  begin
    wb := wb - wA*wx;
    for k in 1..x.deg loop
      idx := idx - 1;
      Aidx := A.cff(idx);
      wA := Aidx.all;
      wx := x.cff(k).all;
      wb := wb - wA*wx;
    end loop;
    wx := Solve(U,V,S,wb);
    x.deg := x.deg + 1;
    x.cff(x.deg) := new QuadDobl_Complex_Vectors.Vector'(wx);
  end Solve_Next_by_SVD;

  procedure Solve_by_lufac
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                info : out integer32;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
   
  begin
    Solve_Lead_by_lufac(A,b,lwrk,ipvt,info,x);
    if info = 0 then
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,x);
      end loop;
    end if;
  end Solve_by_lufac;

  procedure Solve_by_lufco
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                rcond : out quad_double;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    one : constant quad_double := create(1.0);
   
  begin
    Solve_Lead_by_lufco(A,b,lwrk,ipvt,rcond,x);
    if one + rcond /= one then
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,x);
      end loop;
    end if;
  end Solve_by_lufco;

  procedure Solve_by_QRLS
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                info : out integer32;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    nrows : constant integer32 := A.cff(0)'last(1);
    ncols : constant integer32 := A.cff(0)'last(2);
    lwrk : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..ncols);
    qraux : QuadDobl_Complex_Vectors.Vector(1..ncols);
    ipvt : Standard_Integer_Vectors.Vector(1..ncols);

  begin
    Solve_Lead_by_QRLS(A,b,lwrk,qraux,ipvt,info,x);
    if info = 0 then
      for k in 1..b.deg loop
        Solve_Next_by_QRLS(A,b,lwrk,qraux,info,x);
      end loop;
    end if;
  end Solve_by_QRLS;

  procedure Solve_by_SVD
              ( A : in QuadDobl_Dense_Matrix_Series.Matrix;
                b : in QuadDobl_Dense_Vector_Series.Vector;
                info : out integer32; rcond : out quad_double;
                x : out QuadDobl_Dense_Vector_Series.Vector ) is

    nrows : constant integer32 := A.cff(0)'last(1);
    ncols : constant integer32 := A.cff(0)'last(2);
    mm : constant integer32
       := QuadDobl_Complex_Singular_Values.Min0(nrows+1,ncols);
    S : QuadDobl_Complex_Vectors.Vector(1..mm);
    U : QuadDobl_Complex_Matrices.Matrix(1..nrows,1..nrows);
    V : QuadDobl_Complex_Matrices.Matrix(1..ncols,1..ncols);
    one : constant quad_double := create(1.0);

  begin
    Solve_Lead_by_SVD(A,b,S,U,V,info,rcond,x);
    if one + rcond /= one then
      for k in 1..b.deg loop
        Solve_Next_by_SVD(A,b,S,U,V,x);
      end loop;
    end if;
  end Solve_by_SVD;

end QuadDobl_Matrix_Series_Solvers;
