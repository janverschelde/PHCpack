with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_QR_Least_Squares;  use Standard_Complex_QR_Least_Squares;

package body Standard_Matrix_Series_Solvers is

  procedure Solve_Lead_by_lufac
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0lu : out Standard_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector ) is

    lead : constant Standard_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : Standard_Complex_Vectors.Vector(1..dim) := b.cff(0).all;

  begin
    a0lu := lead.all;
    lufac(a0lu,dim,ipvt,info);
    if info /= 0 then
      x.deg := -1;
    else
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new Standard_Complex_Vectors.Vector'(x0);
      x.deg := 0;
    end if;
  end Solve_Lead_by_lufac;

  procedure Solve_Lead_by_lufco
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0lu : out Standard_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                x : out Standard_Dense_Vector_Series.Vector ) is

    lead : constant Standard_Complex_Matrices.Link_to_Matrix := A.cff(0);
    dim : constant integer32 := lead'last(1);
    x0 : Standard_Complex_Vectors.Vector(1..dim) := b.cff(0).all;

  begin
    a0lu := lead.all;
    lufco(a0lu,dim,ipvt,rcond);
    if 1.0 + rcond = 1.0 then
      x.deg := -1;
    else
      lusolve(a0lu,dim,ipvt,x0);
      x.cff(0) := new Standard_Complex_Vectors.Vector'(x0);
      x.deg := 0;
    end if;
  end Solve_Lead_by_lufco;

  procedure Solve_Lead_by_QRLS
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0qr : out Standard_Complex_Matrices.Matrix;
                qraux : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector ) is

    lead : constant Standard_Complex_Matrices.Link_to_Matrix := A.cff(0);
    nrows : constant integer32 := lead'last(1);
    ncols : constant integer32 := lead'last(2);
    b0 : Standard_Complex_Vectors.Vector(1..nrows) := b.cff(0).all;
    x0 : Standard_Complex_Vectors.Vector(1..ncols);
    rsd,dum,dum2,dum3 : Standard_Complex_Vectors.Vector(1..nrows);
    wrk : Standard_Complex_Matrices.Matrix(1..nrows,1..ncols);

  begin
    a0qr := lead.all;
    qraux := (qraux'range => Create(0.0));
    ipvt := (ipvt'range => 0);
    QRD(a0qr,qraux,ipvt,false);
    wrk := a0qr;
    QRLS(wrk,nrows,ncols,qraux,b0,dum2,dum3,x0,rsd,dum,110,info);
    x.cff(0) := new Standard_Complex_Vectors.Vector'(x0);
    x.deg := 0;
  end Solve_Lead_by_QRLS;

  procedure Solve_Next_by_lusolve
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0lu : in Standard_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                x : in out Standard_Dense_Vector_Series.Vector ) is

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

    idx : integer32 := x.deg+1;
    Aidx : Standard_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    dim : constant integer32 := Aidx'last(1);
    wA : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := Aidx.all;
    wb : Standard_Complex_Vectors.Vector(1..dim) := b.cff(idx).all;
    wx : Standard_Complex_Vectors.Vector(1..dim) := x.cff(0).all;

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
    x.cff(x.deg) := new Standard_Complex_Vectors.Vector'(wb);
  end Solve_Next_by_lusolve;

  procedure Solve_Next_by_QRLS
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                a0qr : in Standard_Complex_Matrices.Matrix;
                qraux : in Standard_Complex_Vectors.Vector;
                ipvt : in Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : in out Standard_Dense_Vector_Series.Vector ) is

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

    idx : integer32 := x.deg+1;
    Aidx : Standard_Complex_Matrices.Link_to_Matrix := A.cff(idx);
    nrows : constant integer32 := Aidx'last(1);
    ncols : constant integer32 := Aidx'last(2);
    wA : Standard_Complex_Matrices.Matrix(1..nrows,1..ncols) := Aidx.all;
    wb : Standard_Complex_Vectors.Vector(1..nrows) := b.cff(idx).all;
    wx : Standard_Complex_Vectors.Vector(1..ncols) := x.cff(0).all;
    rsd,dum,dum2,dum3 : Standard_Complex_Vectors.Vector(1..nrows);
    wrk : Standard_Complex_Matrices.Matrix(1..nrows,1..ncols) := a0qr;

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
    x.cff(x.deg) := new Standard_Complex_Vectors.Vector'(wx);
  end Solve_Next_by_QRLS;

  procedure Solve_by_lufac
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
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
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                rcond : out double_float;
                x : out Standard_Dense_Vector_Series.Vector ) is

    dim : constant integer32 := A.cff(0)'last;
    lwrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
   
  begin
    Solve_Lead_by_lufco(A,b,lwrk,ipvt,rcond,x);
    if 1.0 + rcond /= 1.0 then
      for k in 1..b.deg loop
        Solve_Next_by_lusolve(A,b,lwrk,ipvt,x);
      end loop;
    end if;
  end Solve_by_lufco;

  procedure Solve_by_QRLS
              ( A : in Standard_Dense_Matrix_Series.Matrix;
                b : in Standard_Dense_Vector_Series.Vector;
                info : out integer32;
                x : out Standard_Dense_Vector_Series.Vector ) is

    nrows : constant integer32 := A.cff(0)'last(1);
    ncols : constant integer32 := A.cff(0)'last(2);
    lwrk : Standard_Complex_Matrices.Matrix(1..nrows,1..ncols);
    qraux : Standard_Complex_Vectors.Vector(1..ncols);
    ipvt : Standard_Integer_Vectors.Vector(1..ncols);

  begin
    Solve_Lead_by_QRLS(A,b,lwrk,qraux,ipvt,info,x);
    if info = 0 then
      for k in 1..b.deg loop
        Solve_Next_by_QRLS(A,b,lwrk,qraux,ipvt,info,x);
      end loop;
    end if;
  end Solve_by_QRLS;

end Standard_Matrix_Series_Solvers;
