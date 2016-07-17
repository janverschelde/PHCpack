with Standard_Complex_Vectors;
with Standard_Complex_Linear_Solvers;     use Standard_Complex_Linear_Solvers;

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

end Standard_Matrix_Series_Solvers;
