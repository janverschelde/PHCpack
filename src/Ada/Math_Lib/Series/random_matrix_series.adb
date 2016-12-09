with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers; 
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Dense_Series;
with DoblDobl_Dense_Series;
with QuadDobl_Dense_Series;

package body Random_Matrix_Series is

  function Standard_Random_Matrix_Series
             ( deg,dim,lower,upper : integer32 )
             return Standard_Dense_Matrix_Series.Matrix is

    res : Standard_Dense_Matrix_Series.Matrix;

  begin
    if deg > Standard_Dense_Series.max_deg
     then res.deg := Standard_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
        rnd : integer32;
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            rnd := Standard_Random_Numbers.Random(lower,upper);
            mat(i,j) := Standard_Complex_Numbers.Create(double_float(rnd));
          end loop;
        end loop;
        res.cff(d) := new Standard_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Standard_Random_Matrix_Series;

  function DoblDobl_Random_Matrix_Series
             ( deg,dim,lower,upper : integer32 )
             return DoblDobl_Dense_Matrix_Series.Matrix is

    res : DoblDobl_Dense_Matrix_Series.Matrix;

  begin
    if deg > DoblDobl_Dense_Series.max_deg
     then res.deg := DoblDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
        rnd : integer32;
        ddr : double_double;
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            rnd := Standard_Random_Numbers.Random(lower,upper);
            ddr := Double_Double_Numbers.Create(rnd);
            mat(i,j) := DoblDobl_Complex_Numbers.Create(ddr);
          end loop;
        end loop;
        res.cff(d) := new DoblDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Matrix_Series;

  function QuadDobl_Random_Matrix_Series
             ( deg,dim,lower,upper : integer32 )
             return QuadDobl_Dense_Matrix_Series.Matrix is

    res : QuadDobl_Dense_Matrix_Series.Matrix;

  begin
    if deg > QuadDobl_Dense_Series.max_deg
     then res.deg := QuadDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
        rnd : integer32;
        qdr : quad_double;
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            rnd := Standard_Random_Numbers.Random(lower,upper);
            qdr := Quad_Double_Numbers.Create(rnd);
            mat(i,j) := QuadDobl_Complex_Numbers.Create(qdr);
          end loop;
        end loop;
        res.cff(d) := new QuadDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Matrix_Series;

  function Standard_Random_Matrix_Series
             ( deg,dim : integer32 )
             return Standard_Dense_Matrix_Series.Matrix is

    res : Standard_Dense_Matrix_Series.Matrix;

  begin
    if deg > Standard_Dense_Series.max_deg
     then res.deg := Standard_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := Standard_Random_Numbers.Random1;
          end loop;
        end loop;
        res.cff(d) := new Standard_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Standard_Random_Matrix_Series;

  function DoblDobl_Random_Matrix_Series
             ( deg,dim : integer32 )
             return DoblDobl_Dense_Matrix_Series.Matrix is

    res : DoblDobl_Dense_Matrix_Series.Matrix;

  begin
    if deg > DoblDobl_Dense_Series.max_deg
     then res.deg := DoblDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := DoblDobl_Random_Numbers.Random1;
          end loop;
        end loop;
        res.cff(d) := new DoblDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Matrix_Series;

  function QuadDobl_Random_Matrix_Series
             ( deg,dim : integer32 )
             return QuadDobl_Dense_Matrix_Series.Matrix is

    res : QuadDobl_Dense_Matrix_Series.Matrix;

  begin
    if deg > QuadDobl_Dense_Series.max_deg
     then res.deg := QuadDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for d in 0..res.deg loop
      declare
        mat : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := QuadDobl_Random_Numbers.Random1;
          end loop;
        end loop;
        res.cff(d) := new QuadDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Matrix_Series;

end Random_Matrix_Series;
