with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;

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

end Hessian_Convolution_Circuits;
