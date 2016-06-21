with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Dense_Series_io;
with Standard_Dense_Series_Matrices;
with Standard_Linear_Series_Solvers;    use Standard_Linear_Series_Solvers;
with Standard_Series_Polynomials;
with Series_and_Polynomials;
with Standard_Series_Poly_SysFun;
with Standard_Series_Jaco_Matrices;

package body Standard_Newton_Series is

  procedure LU_Newton_Step
              ( p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : Standard_Dense_Series_Vectors.Vector(x'range);
    px : Standard_Dense_Series_Vectors.Vector(p'range);
    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);
    jm : Standard_Dense_Series_Matrices.Matrix(p'range,x'range);
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := Standard_Series_Poly_SysFun.Eval(p,x);
    Standard_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := Standard_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    LUfac(jm,n,ipvt,info);
    if info = 0 then
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      Standard_Dense_Series_Vectors.Add(x,dx);
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                p : in Standard_Series_Poly_Systems.Poly_Sys;
                order : in integer32;
                x : in out Standard_Dense_Series_Vectors.Vector;
                info : out integer32 ) is

    dx : Standard_Dense_Series_Vectors.Vector(x'range);
    px : Standard_Dense_Series_Vectors.Vector(p'range);
    jp : Standard_Series_Jaco_Matrices.Jaco_Mat(p'range,x'range)
       := Standard_Series_Jaco_Matrices.Create(p);
    jm : Standard_Dense_Series_Matrices.Matrix(p'range,x'range);
    n : constant integer32 := p'last;
    ipvt : Standard_Integer_Vectors.Vector(1..n);

  begin
    px := Standard_Series_Poly_SysFun.Eval(p,x);
    put_line(file,"The evaluated series :");
    for i in px'range loop
      Standard_Dense_Series_io.put(file,px(i)); new_line(file);
    end loop;
    Standard_Dense_Series_Vectors.Min(px);
    Series_and_Polynomials.Set_Order(px,order);
    jm := Standard_Series_Jaco_Matrices.Eval(jp,x);
    Series_and_Polynomials.Set_Order(jm,order);
    LUfac(jm,n,ipvt,info);
    if info /= 0 then
      put(file,"LUfac info : "); put(file,info,1); new_line(file);
    else
      dx := px;
      LUsolve(jm,n,ipvt,dx);
      put_line(file,"The update to the series :");
      for i in dx'range loop
        Standard_Dense_Series_io.put(file,dx(i)); new_line(file);
      end loop;
      Standard_Dense_Series_Vectors.Add(x,dx);
    end if;
  end LU_Newton_Step;

end Standard_Newton_Series;
