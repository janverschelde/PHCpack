with integer_io;                        use integer_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Complex_Linear_Solvers;   use Standard_Complex_Linear_Solvers;
with Continuation_Parameters;
with Increment_and_Fix_Continuation;    use Increment_and_Fix_Continuation;
with Standard_Affine_Planes;            use Standard_Affine_Planes;
with Standard_Affine_Solutions;         use Standard_Affine_Solutions;

package body Intrinsic_Sampling_Machine is

  function Combine ( c : Vector; v : VecVec ) return Vector is

  -- DESCRIPTION :
  --   Returns the linear combination of the vectors in v,
  --   each multiplied with the corresponding coefficient in c.

    res : Vector(v(v'first)'range) := (v(v'first)'range => Create(0.0));

  begin
    for i in v'range loop
      for j in res'range loop
        res(j) := res(j) + c(i)*v(i)(j);
      end loop;
    end loop;
    return res;
  end Combine;

  function Eval ( p : Eval_Poly_Sys; v : VecVec; c : Vector )
                return Vector is

    x : constant Vector := Combine(c,v);

  begin
    return Eval(p,x);
  end Eval;

  function Eval ( p : Eval_Poly_Sys; b : Vector; v : VecVec; c : Vector )
                return Vector is

    x : constant Vector := Combine(c,b,v);

  begin
    return Eval(p,x);
  end Eval;

  function Chain_Rule ( eva : Matrix; v : VecVec ) return Matrix is

  -- DESCRIPTION :
  --   Applies the chain rule to the matrix obtained after evaluating
  --   the Jacobi matrix at some point, with respect to the coefficients
  --   in the linear combinations of the vectors in v.

    res : Matrix(eva'range(1),v'range);

  begin
    for i in res'range(1) loop
      for j in v'range loop
        res(i,j) := Create(0.0);
        for k in eva'range(2) loop
          res(i,j) := res(i,j) + eva(i,k)*v(j)(k);
        end loop;
      end loop;
    end loop;
    return res;
  end Chain_Rule;

  function Diff ( jm : Eval_Jaco_Mat; v : VecVec; c : Vector )
                return Matrix is

    x : constant Vector := Combine(c,v);
    eva : Matrix(jm'range(1),jm'range(2)) := Eval(jm,x);

  begin
    return Chain_Rule(eva,v);
  end Diff;

  function Diff ( jm : Eval_Jaco_Mat; b : Vector; v : VecVec; c : Vector )
                return Matrix is

    x : constant Vector := Combine(c,b,v);
    eva : Matrix(jm'range(1),jm'range(2)) := Eval(jm,x);

  begin
    return Chain_Rule(eva,v);
  end Diff;

  procedure Moving_Plane
              ( start_b,target_b : in Vector; start_v,target_v : in VecVec;
                t : in Complex_Number; b : out Vector; v : in out VecVec ) is

    one_min_t : constant Complex_Number := 1.0 - t;

  begin
    for i in b'range loop
      b(i) := one_min_t*start_b(i) + t*target_b(i);
    end loop;
    for i in v'range loop
      if v(i) = null
       then v(i) := new Standard_Complex_Vectors.Vector(b'range);
      end if;
      for j in v(i)'range loop
        v(i)(j) := one_min_t*start_v(i)(j) + t*target_v(i)(j);
      end loop;
    end loop;
  end Moving_Plane;

  procedure Affine_Newton_Refinement
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; c : in out Vector;
                epsxa,epsfa : in double_float; maxit : in natural ) is

    k : constant natural := c'length;
    y : Vector(p'range);
    jac : Matrix(jm'range(1),c'range);
    ipvt : Standard_Natural_Vectors.Vector(c'range);
    dc : Vector(c'range);
    nc,ny,rcond : double_float;

  begin
    y := Eval(p,b,v,c);
    for i in 1..maxit loop
      jac := Diff(jm,b,v,c);
      lufco(jac,k,ipvt,rcond);
      dc := -y;
      lusolve(jac,k,ipvt,dc);
      c := c + dc;
      y := Eval(p,b,v,c);
      nc := Max_Norm(c);
      ny := Max_Norm(y);
      put("  |dc| : "); put(nc,3);
      put("  |y| : "); put(ny,3); new_line;
      exit when (nc < epsxa) or (ny < epsfa);
    end loop;
    put_line("After Newton refinement : "); put_line(y);
  end Affine_Newton_Refinement;

  procedure Affine_Solution_Evaluation
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                sol,b : Vector; v : in VecVec ) is

    c : Vector(v'range) := Decompose(sol,b,v);
    y : Vector(p'range) := Eval(p,b,v,c);
    t : Vector(v'range) := (v'range => Create(0.0));

  begin
    put_line("The solution evaluated at computed basis : "); put_line(y);
    Affine_Newton_Refinement(p,jm,b,v,c,1.0E-16,1.0E-14,4);
    y := Eval(p,sol,v,t);
    put_line("The solution evaluated at its own basis : "); put_line(y);
  end Affine_Solution_Evaluation;

  procedure Silent_Affine_Sampler
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(p,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Sil_Cont is
      new Silent_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Sil_Cont(sols,false,Create(1.0));
  end Silent_Affine_Sampler;

  procedure Reporting_Affine_Sampler
              ( file : in file_type;
                p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                start_b,target_b : in Vector; start_v,target_v : in VecVec;
                sols : in out Solution_List ) is

    tb : Vector(start_b'range);
    tv : VecVec(start_v'range);

    function Eval ( c : Vector; t : Complex_Number ) return Vector is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Eval(p,tb,tv,c);
    end Eval;

    function Diff ( c : Vector; t : Complex_Number ) return Matrix is
    begin
      Moving_Plane(start_b,target_b,start_v,target_v,t,tb,tv);
      return Diff(jm,tb,tv,c);
    end Diff;

    function Diff ( c : Vector; t : Complex_Number ) return Vector is
    begin
      return c;
    end Diff;

    procedure Rep_Cont is
      new Reporting_Continue(Max_Norm,Eval,Diff,Diff);

  begin
    Rep_Cont(file,sols,false,Create(1.0));
  end Reporting_Affine_Sampler;

  procedure Silent_Newton_Refiner
              ( p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in out Solution_List;
                epsxa,epsfa : in double_float; maxit : in natural ) is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
      begin
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Silent_Newton_Refiner;

  procedure Reporting_Newton_Refiner
              ( file : in file_type;
                p : in Eval_Poly_Sys; jm : in Eval_Jaco_Mat;
                b : in Vector; v : in VecVec; sols : in out Solution_List;
                epsxa,epsfa : in double_float; maxit : in natural ) is
  begin
    null;
  end Reporting_Newton_Refiner;

end Intrinsic_Sampling_Machine;
