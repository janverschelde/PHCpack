with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;

package body Hessian_Convolution_Circuits is

  function Hessian ( c : Standard_Speelpenning_Convolutions.Circuit;
                     x : Standard_Complex_Vectors.Vector ) 
                   return Standard_Complex_Matrices.Matrix is

    dim : constant integer32 := c.dim;
    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);

    use Standard_Speelpenning_Convolutions;

  begin
    for i in 1..dim loop
      res(i,i) := Diff(c,x,i,i);
      for j in (i+1)..dim loop
        res(i,j) := Diff(c,x,i,j);
        res(j,i) := res(i,j);
      end loop;
    end loop;
    return res;
  end Hessian;

  function Hessian ( c : DoblDobl_Speelpenning_Convolutions.Circuit;
                     x : DoblDobl_Complex_Vectors.Vector ) 
                   return DoblDobl_Complex_Matrices.Matrix is

    dim : constant integer32 := c.dim;
    res : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);

    use DoblDobl_Speelpenning_Convolutions;

  begin
    for i in 1..dim loop
      res(i,i) := Diff(c,x,i,i);
      for j in (i+1)..dim loop
        res(i,j) := Diff(c,x,i,j);
        res(j,i) := res(i,j);
      end loop;
    end loop;
    return res;
  end Hessian;

  function Hessian ( c : QuadDobl_Speelpenning_Convolutions.Circuit;
                     x : QuadDobl_Complex_Vectors.Vector ) 
                   return QuadDobl_Complex_Matrices.Matrix is

    dim : constant integer32 := c.dim;
    res : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);

    use QuadDobl_Speelpenning_Convolutions;

  begin
    for i in 1..dim loop
      res(i,i) := Diff(c,x,i,i);
      for j in (i+1)..dim loop
        res(i,j) := Diff(c,x,i,j);
        res(j,i) := res(i,j);
      end loop;
    end loop;
    return res;
  end Hessian;

  function Hessian ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit;
                     x : Standard_Complex_Vectors.Vector ) 
                   return Standard_Complex_Matrices.Matrix is

    use Standard_Speelpenning_Convolutions;

  begin
    if c /= null then
      return Hessian(c.all,x);
    else
      declare
        dim : constant integer32 := x'last;
        res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            res(i,j) := Standard_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        return res;
      end;
    end if;
  end Hessian;

  function Hessian ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                     x : DoblDobl_Complex_Vectors.Vector ) 
                   return DoblDobl_Complex_Matrices.Matrix is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if c /= null then
      return Hessian(c.all,x);
    else
      declare
        dim : constant integer32 := x'last;
        res : DoblDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            res(i,j) := DoblDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        return res;
      end;
    end if;
  end Hessian;

  function Hessian ( c : QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                     x : QuadDobl_Complex_Vectors.Vector ) 
                   return QuadDobl_Complex_Matrices.Matrix is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    if c /= null then
      return Hessian(c.all,x);
    else
      declare
        dim : constant integer32 := x'last;
        res : QuadDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            res(i,j) := QuadDobl_Complex_Numbers.Create(integer(0));
          end loop;
        end loop;
        return res;
      end;
    end if;
  end Hessian;

  function Hessians ( c : Standard_Speelpenning_Convolutions.Circuits;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_VecMats.VecMat is

    res : Standard_Complex_VecMats.VecMat(c'range);

  begin
    for i in c'range loop
      res(i) := new Standard_Complex_Matrices.Matrix'(Hessian(c(i),x));
    end loop;
    return res;
  end Hessians;

  function Hessians ( c : DoblDobl_Speelpenning_Convolutions.Circuits;
                      x : DoblDobl_Complex_Vectors.Vector )
                    return DoblDobl_Complex_VecMats.VecMat is

    res : DoblDobl_Complex_VecMats.VecMat(c'range);

  begin
    for i in c'range loop
      res(i) := new DoblDobl_Complex_Matrices.Matrix'(Hessian(c(i),x));
    end loop;
    return res;
  end Hessians;

  function Hessians ( c : QuadDobl_Speelpenning_Convolutions.Circuits;
                      x : QuadDobl_Complex_Vectors.Vector )
                    return QuadDobl_Complex_VecMats.VecMat is

    res : QuadDobl_Complex_VecMats.VecMat(c'range);

  begin
    for i in c'range loop
      res(i) := new QuadDobl_Complex_Matrices.Matrix'(Hessian(c(i),x));
    end loop;
    return res;
  end Hessians;

  procedure Singular_Values
              ( c : in Standard_Speelpenning_Convolutions.Circuit;
                x : in Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    info : integer32;

  begin
    A := Hessian(c,x);
    Standard_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,0,info);
  end Singular_Values;

  procedure Singular_Values
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuit;
                x : in DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix;
                U : out DoblDobl_Complex_Matrices.Matrix;
                V : out DoblDobl_Complex_Matrices.Matrix;
                e : out DoblDobl_Complex_Vectors.Vector;
                s : out DoblDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    info : integer32;

  begin
    A := Hessian(c,x);
    DoblDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,0,info);
  end Singular_Values;

  procedure Singular_Values
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuit;
                x : in QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix;
                U : out QuadDobl_Complex_Matrices.Matrix;
                V : out QuadDobl_Complex_Matrices.Matrix;
                e : out QuadDobl_Complex_Vectors.Vector;
                s : out QuadDobl_Complex_Vectors.Vector ) is

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    job : constant integer32 := 0;
    info : integer32;

  begin
    A := Hessian(c,x);
    QuadDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,0,info);
  end Singular_Values;

  procedure Singular_Values
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                x : in Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix;
                U : out Standard_Complex_Matrices.Matrix;
                V : out Standard_Complex_Matrices.Matrix;
                e : out Standard_Complex_Vectors.Vector;
                s : out Standard_Complex_Vectors.Vector ) is

    use Standard_Speelpenning_Convolutions;

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    job : constant integer32 := 0;
    info : integer32;

  begin
    A := Hessian(c,x);
    if c /= null then
      Standard_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,0,info);
    else
      for i in s'range loop
        s(i) := Standard_Complex_Numbers.Create(integer(0));
      end loop;
    end if;
  end Singular_Values;

  procedure Singular_Values
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                x : in DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix;
                U : out DoblDobl_Complex_Matrices.Matrix;
                V : out DoblDobl_Complex_Matrices.Matrix;
                e : out DoblDobl_Complex_Vectors.Vector;
                s : out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Speelpenning_Convolutions;

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    job : constant integer32 := 0;
    info : integer32;

  begin
    A := Hessian(c,x);
    if c /= null then
      DoblDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,0,info);
    else
      for i in s'range loop
        s(i) := DoblDobl_Complex_Numbers.Create(integer(0));
      end loop;
    end if;
  end Singular_Values;

  procedure Singular_Values
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                x : in QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix;
                U : out QuadDobl_Complex_Matrices.Matrix;
                V : out QuadDobl_Complex_Matrices.Matrix;
                e : out QuadDobl_Complex_Vectors.Vector;
                s : out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Speelpenning_Convolutions;

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    job : constant integer32 := 0;
    info : integer32;

  begin
    A := Hessian(c,x);
    if c /= null then
      QuadDobl_Complex_Singular_Values.SVD(A,n,p,s,e,u,v,0,info);
    else
      for i in s'range loop
        s(i) := QuadDobl_Complex_Numbers.Create(integer(0));
      end loop;
    end if;
  end Singular_Values;

end Hessian_Convolution_Circuits;
