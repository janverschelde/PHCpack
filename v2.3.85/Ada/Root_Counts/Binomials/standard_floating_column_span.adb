with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_GramSchmidt;

package body Standard_Floating_Column_Span is

-- NOTE : the modified Gram-Schmidt method is applied to
--   decide whether a vector belongs to the span of vectors.

  function In_Span ( v : Standard_Integer_VecVecs.VecVec;
                     x : Standard_Integer_Vectors.Vector )
                   return boolean is

    res : boolean;
    fv : Standard_Floating_VecVecs.VecVec(v'range);
    fx : Standard_Floating_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      fx(i) := double_float(x(i));
    end loop;
    for j in v'range loop
      fv(j) := new Standard_Floating_Vectors.Vector(v(j)'range);
      for i in v(j)'range loop
        fv(j)(i) := double_float(v(j)(i));
      end loop;
    end loop;
    res := In_Span(fv,fx);
    Standard_Floating_VecVecs.Clear(fv);
    return res;
  end In_Span;

  function In_Span ( v : Standard_Floating_VecVecs.VecVec;
                     x : Standard_Floating_Vectors.Vector )
                   return boolean is
  begin
    return In_Span(v,x,1.0E-8);
  end In_Span;

  procedure Initialize ( v : in Standard_Floating_VecVecs.VecVec;
                         q,r : out Standard_Floating_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Copies v into q and initializes r to zero.

  -- REQUIRED : v'range = q'range = r'range.

  begin
    for i in v'range loop
      q(i) := new Standard_Floating_Vectors.Vector'(v(i).all);
      r(i) := new Standard_Floating_Vectors.Vector(r'range);
      for j in r'range loop
        r(i)(j) := 0.0;
      end loop;
    end loop;
  end Initialize;

  function Residual ( v : Standard_Floating_VecVecs.VecVec;
                      x,s : Standard_Floating_Vectors.Vector )
                    return double_float is

  -- DESCRIPTION :
  --   Returns the residual of the least squares solution s,
  --   as the 1-norm between the difference of the components
  --   of x and the linear combination of s with v.

    res : double_float := 0.0;
    sum : double_float;

  begin
    for i in x'range loop
      sum := 0.0;
      for j in s'range loop
        sum := sum + s(j)*v(j)(i);
      end loop;
      res := res + abs(sum - x(i));
    end loop;
    return res;
  end Residual;

  function In_Span ( v : Standard_Floating_VecVecs.VecVec;
                     x : Standard_Floating_Vectors.Vector;
                     tol : double_float ) return boolean is

    n : constant integer32 := x'last;
    m : constant integer32 := v'last;
    q,r : Standard_Floating_VecVecs.VecVec(v'range);
    qx,sol : Standard_Floating_Vectors.Vector(v'range);
    res : double_float;
   

  begin
    Initialize(v,q,r);
    Standard_Floating_GramSchmidt.QR(n,m,q,r);
    qx := Standard_Floating_GramSchmidt.Matrix_Projection(n,m,q,x);
    sol := Standard_Floating_GramSchmidt.Solve(m,r,qx);
    res := Residual(v,x,sol);
    Standard_Floating_VecVecs.Clear(q);
    Standard_Floating_VecVecs.Clear(r);
    return (res <= tol);
  end In_Span;

  function In_Span ( file : file_type;
                     v : Standard_Floating_VecVecs.VecVec;
                     x : Standard_Floating_Vectors.Vector;
                     tol : double_float ) return boolean is

    n : constant integer32 := x'last;
    m : constant integer32 := v'last;
    q,r : Standard_Floating_VecVecs.VecVec(v'range);
    qx,sol : Standard_Floating_Vectors.Vector(v'range);
    res,maxerr : double_float;
    fail : boolean;

  begin
    Initialize(v,q,r);
    Standard_Floating_GramSchmidt.QR(n,m,q,r);
    Standard_Floating_GramSchmidt.Test_Decomposition
      (n,m,v,q,r,tol,true,maxerr,fail);
    qx := Standard_Floating_GramSchmidt.Matrix_Projection(n,m,q,x);
    put(file,"the projected vector : "); put(file,qx); new_line(file);
    sol := Standard_Floating_GramSchmidt.Solve(m,r,qx);
    put(file,"the solution : "); put(file,sol); new_line(file);
    res := Residual(v,x,sol);
    put(file,"the residual : "); put(file,res); new_line(file);
    Standard_Floating_VecVecs.Clear(q);
    Standard_Floating_VecVecs.Clear(r);
    return (res <= tol);
  end In_Span;

end Standard_Floating_Column_Span;
