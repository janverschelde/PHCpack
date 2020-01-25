with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Exponent_Indices;
with Standard_Complex_Series;
with Standard_Complex_Series_io;         use Standard_Complex_Series_io;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_Vectors_io; use Standard_Complex_Series_Vectors_io;
with Standard_Complex_Series_Matrices;
with Standard_Random_Series_Vectors;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_io;         use DoblDobl_Complex_Series_io;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors_io; use DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_Random_Series_Vectors;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_io;         use QuadDobl_Complex_Series_io;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors_io; use QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_Complex_Series_Matrices;
with QuadDobl_Random_Series_Vectors;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Polynomials_io;    use Standard_CSeries_Polynomials_io;
with Standard_CSeries_Poly_Functions;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_Systems_io;
with Standard_CSeries_Poly_SysFun;
with Standard_CSeries_Jaco_Matrices;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Polynomials_io;    use DoblDobl_CSeries_Polynomials_io;
with DoblDobl_CSeries_Poly_Functions;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems_io;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Polynomials_io;    use QuadDobl_CSeries_Polynomials_io;
with QuadDobl_CSeries_Poly_Functions;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems_io;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Series_Coefficient_Vectors;         use Series_Coefficient_Vectors;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Complex_Series_and_Polynomials;     use Complex_Series_and_Polynomials;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;

procedure ts_speelcnv is

-- DESCRIPTION :
--   Tests the evaluation of the gradient of a polynomial in many variables,
--   in a power series of some fixed degree.

  procedure Standard_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in standard double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Speelpenning_Convolutions;

    xps : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Exponents(dim,nbr,pwr);
    idx : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Exponent_Index(xps);
    fac : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Factor_Index(xps);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Indices.Maxima(xps);
    polcff : constant Standard_Complex_Series_Vectors.Vector(1..nbr)
           := Standard_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    pol : constant Standard_CSeries_Polynomials.Poly
       -- := Standard_Polynomial(dim,deg,idx); -- all coefficients are one
       -- := Standard_Polynomial(dim,idx,polcff); -- all exponents are one
        := Standard_Polynomial(dim,xps,polcff,false);
    x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    y : Standard_Complex_Series.Link_to_Series;
    grad : Standard_Complex_Series_Vectors.Vector(1..dim);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    pcff : constant Standard_Complex_VecVecs.VecVec(1..nbr)
         := Standard_Series_Coefficients(polcff);
    forward : constant Standard_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : constant Standard_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : constant Standard_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    ygrad : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    work : constant Standard_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);
    acc : constant Standard_Complex_Vectors.Link_to_Vector
        := Allocate_Coefficients(deg);
    err,sumerr : double_float := 0.0;
    pwt : Link_to_VecVecVec := Create(xcff,mxe);

  begin
    put_line("Some random exponents :");
    Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :");
    Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :");
    Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := Standard_CSeries_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
   -- Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
    Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficient vector of the value of the product :");
    put_line(ygrad(ygrad'last));
    grad := Standard_Gradient(pol,x);
    err := Difference(y,ygrad(ygrad'last));
    put("The error :"); put(err,3); new_line;
    sumerr := err;
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :");
      put(grad(k)); new_line;
      put("The coefficient vector of derivative ");
      put(k,1); put_line(" :"); put_line(ygrad(k));
      err := Difference(grad(k),ygrad(k));
      put("The error :"); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("Sum of errors :"); put(sumerr,3); new_line;
    Clear(pwt);
  end Standard_Test;

  procedure DoblDobl_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Speelpenning_Convolutions;

    xps : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Exponents(dim,nbr,pwr);
    idx : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Exponent_Index(xps);
    fac : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Factor_Index(xps);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Indices.Maxima(xps);
    polcff : constant DoblDobl_Complex_Series_Vectors.Vector(1..nbr)
           := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    pol : constant DoblDobl_CSeries_Polynomials.Poly
       -- := DoblDobl_Polynomial(dim,deg,idx); -- all coefficients are one
       -- := DoblDobl_Polynomial(dim,idx,polcff); -- all exponents are one
        := DoblDobl_Polynomial(dim,xps,polcff,false);
    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    y : DoblDobl_Complex_Series.Link_to_Series;
    grad : DoblDobl_Complex_Series_Vectors.Vector(1..dim);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    pcff : constant DoblDobl_Complex_VecVecs.VecVec(1..nbr)
         := DoblDobl_Series_Coefficients(polcff);
    forward : constant DoblDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : constant DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : constant DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    ygrad : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    work : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);
    acc : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := Allocate_Coefficients(deg);
    err,sumerr : double_double;
    pwt : Link_to_VecVecVec := Create(xcff,mxe);

  begin
    put_line("Some random exponents :");
    Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :");
    Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :");
    Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := DoblDobl_CSeries_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
   -- Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
    Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficient vector of the value of the product :");
    put_line(ygrad(ygrad'last));
    grad := DoblDobl_Gradient(pol,x);
    err := Difference(y,ygrad(ygrad'last));
    put("The error : "); put(err,3); new_line;
    sumerr := err;
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :");
      put(grad(k)); new_line;
      put("The coefficient vector of derivative ");
      put(k,1); put_line(" :"); put_line(ygrad(k));
      err := Difference(grad(k),ygrad(k));
      put("The error : "); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("Sum of errors : "); put(sumerr,3); new_line;
    Clear(pwt);
  end DoblDobl_Test;

  procedure QuadDobl_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in quad double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Speelpenning_Convolutions;

    xps : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Random_Exponents(dim,nbr,pwr);
    idx : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Exponent_Index(xps);
    fac : constant Standard_Integer_VecVecs.VecVec(1..nbr)
        := Exponent_Indices.Factor_Index(xps);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Indices.Maxima(xps);
    polcff : constant QuadDobl_Complex_Series_Vectors.Vector(1..nbr)
           := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,nbr,deg);
    pol : constant QuadDobl_CSeries_Polynomials.Poly
       -- := QuadDobl_Polynomial(dim,deg,idx); -- all coefficients are one
       -- := QuadDobl_Polynomial(dim,idx,polcff); -- all exponents are one
        := QuadDobl_Polynomial(dim,xps,polcff,false);
    x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    y : QuadDobl_Complex_Series.Link_to_Series;
    grad : QuadDobl_Complex_Series_Vectors.Vector(1..dim);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    pcff : constant QuadDobl_Complex_VecVecs.VecVec(1..nbr)
         := QuadDobl_Series_Coefficients(polcff);
    forward : constant QuadDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    backward : constant QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    cross : constant QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    ygrad : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    work : constant QuadDobl_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);
    acc : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := Allocate_Coefficients(deg);
    err,sumerr : quad_double;
    pwt : Link_to_VecVecVec := Create(xcff,mxe);

  begin
    put_line("Some random exponents :");
    Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :");
    Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :");
    Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := QuadDobl_CSeries_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
   -- Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
    Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    put_line("The coefficient vector of the value of the product :");
    put_line(ygrad(ygrad'last));
    grad := QuadDobl_Gradient(pol,x);
    err := Difference(y,ygrad(ygrad'last));
    put("The error : "); put(err,3); new_line;
    sumerr := err;
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :");
      put(grad(k)); new_line;
      put("The coefficient vector of derivative ");
      put(k,1); put_line(" :"); put_line(ygrad(k));
      err := Difference(grad(k),ygrad(k));
      put("The error : "); put(err,3); new_line;
      sumerr := sumerr + err;
    end loop;
    put("Sum of errors : "); put(sumerr,3); new_line;
    Clear(pwt);
  end QuadDobl_Test;

  procedure Standard_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in standard double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Speelpenning_Convolutions;

    c : constant Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : Link_to_System := Create(c,dim,deg);
    p : Standard_CSeries_Poly_Systems.Poly_Sys(1..dim) := Standard_System(c);
    x : Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    px : Standard_Complex_Series_Vectors.Vector(p'range)
       := Standard_CSeries_Poly_SysFun.Eval(p,x);
    xcff : Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    jp : Standard_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
       := Standard_CSeries_Jaco_Matrices.Create(p);
    jm : Standard_Complex_Series_Matrices.Matrix(1..dim,1..dim)
       := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
    err : double_float;

  begin
    put_line("the polynomial system :");
    Standard_CSeries_Poly_Systems_io.put(p);
    Compute(s.pwt,s.mxe,xcff);
    EvalDiff(s,xcff);
    put_line("The evaluation values :"); put_line(px);
    put_line("The coefficient vectors of the evaluation :"); put_line(s.yv);
    err := Difference(px,s.yv);
    put("The error :"); put(err,3); new_line;
    for i in s.vm'range loop
      put("The matrix "); put(i,1); put_line(" :");
      put(s.vm(i).all);
    end loop;
    for i in 1..dim loop
      for j in 1..dim loop
        put("the series of the Jacobian at ");
        put(i,1); put(" and "); put(j,1); put_line(" :");
        put(jm(i,j));
      end loop;
    end loop;
    err := Difference(jm,s.vm);
    put("The error :"); put(err,3); new_line;
    Clear(s); -- Clear(c) is not needed, the Clear(s) does Clear(s.crc)
    Standard_CSeries_Poly_Systems.Clear(p);
    Standard_CSeries_Jaco_Matrices.Clear(jp);
    Standard_Complex_Series_Vectors.Clear(x);
    Standard_Complex_Series_Vectors.Clear(px);
    Standard_Complex_Series_Matrices.Clear(jm);
    Standard_Complex_VecVecs.Clear(xcff);
  end Standard_System_Test;

  procedure DoblDobl_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : Link_to_System := Create(c,dim,deg);
    p : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..dim) := DoblDobl_System(c);
    x : DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    px : DoblDobl_Complex_Series_Vectors.Vector(p'range)
       := DoblDobl_CSeries_Poly_SysFun.Eval(p,x);
    xcff : DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    jp : DoblDobl_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
       := DoblDobl_CSeries_Jaco_Matrices.Create(p);
    jm : DoblDobl_Complex_Series_Matrices.Matrix(1..dim,1..dim)
       := DoblDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    err : double_double;

  begin
    put_line("the polynomial system :");
    DoblDobl_CSeries_Poly_Systems_io.put(p);
    Compute(s.pwt,s.mxe,xcff);
    EvalDiff(s,xcff);
    put_line("The evaluation values :"); put_line(px);
    put_line("The coefficient vectors of the evaluation :"); put_line(s.yv);
    err := Difference(px,s.yv);
    put("The error : "); put(err,3); new_line;
    for i in s.vm'range loop
      put("The matrix "); put(i,1); put_line(" :");
      put(s.vm(i).all);
    end loop;
    for i in 1..dim loop
      for j in 1..dim loop
        put("the series of the Jacobian at ");
        put(i,1); put(" and "); put(j,1); put_line(" :");
        put(jm(i,j));
      end loop;
    end loop;
    err := Difference(jm,s.vm);
    put("The error : "); put(err,3); new_line;
    Clear(s);
    DoblDobl_CSeries_Poly_Systems.Clear(p);
    DoblDobl_CSeries_Jaco_Matrices.Clear(jp);
    DoblDobl_Complex_Series_Vectors.Clear(x);
    DoblDobl_Complex_Series_Vectors.Clear(px);
    DoblDobl_Complex_Series_Matrices.Clear(jm);
    DoblDobl_Complex_VecVecs.Clear(xcff);
  end DoblDobl_System_Test;

  procedure QuadDobl_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in double double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : Link_to_System := Create(c,dim,deg);
    p : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..dim) := QuadDobl_System(c);
    x : QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range)
       := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    xcff : QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(1..dim,1..dim)
       := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    err : quad_double;

  begin
    put_line("the polynomial system :");
    QuadDobl_CSeries_Poly_Systems_io.put(p);
    Compute(s.pwt,s.mxe,xcff);
    EvalDiff(s,xcff);
    put_line("The evaluation values :"); put_line(px);
    put_line("The coefficient vectors of the evaluation :"); put_line(s.yv);
    err := Difference(px,s.yv);
    put("The error : "); put(err,3); new_line;
    for i in s.vm'range loop
      put("The matrix "); put(i,1); put_line(" :");
      put(s.vm(i).all);
    end loop;
    for i in 1..dim loop
      for j in 1..dim loop
        put("the series of the Jacobian at ");
        put(i,1); put(" and "); put(j,1); put_line(" :");
        put(jm(i,j));
      end loop;
    end loop;
    err := Difference(jm,s.vm);
    put("The error : "); put(err,3); new_line;
    Clear(s);
    QuadDobl_CSeries_Poly_Systems.Clear(p);
    QuadDobl_CSeries_Jaco_Matrices.Clear(jp);
    QuadDobl_Complex_Series_Vectors.Clear(x);
    QuadDobl_Complex_Series_Vectors.Clear(px);
    QuadDobl_Complex_Series_Matrices.Clear(jm);
    QuadDobl_Complex_VecVecs.Clear(xcff);
  end QuadDobl_System_Test;

  procedure Standard_Input_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   makes convolution circuits, evaluates and differentiates
  --   at a vector of random series, in double precision.

    p : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    d : constant natural32 := natural32(deg);

    use Standard_Speelpenning_Convolutions;

  begin
    put_line("Reading a polynomial system ..."); get(p);
    declare
      s : constant Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
        := System_to_Series_System(p.all);
      dim : constant integer32 := p'last;
      c : constant Circuits(p'range) := Make_Convolution_Circuits(p.all,d);
      q : Link_to_System := Create(c,dim,deg);
      x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
        := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
      jp : constant Standard_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
         := Standard_CSeries_Jaco_Matrices.Create(s);
      jm : constant Standard_Complex_Series_Matrices.Matrix(1..dim,1..dim)
         := Standard_CSeries_Jaco_Matrices.Eval(jp,x);
      sx : constant Standard_Complex_Series_Vectors.Vector(p'range)
         := Standard_CSeries_Poly_SysFun.Eval(s,x);
      xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
           := Standard_Series_Coefficients(x);
      err : double_float;
    begin
      Compute(q.pwt,q.mxe,xcff);
      EvalDiff(q,xcff);
      err := Difference(sx,q.yv);
      put("The evaluation error : "); put(err,3); new_line;
      err := Difference(jm,q.vm);
      put("The differentiation error : "); put(err,3); new_line;
      Clear(q);
    end;
  end Standard_Input_Test;

  procedure DoblDobl_Input_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   makes convolution circuits, evaluates and differentiates
  --   at a vector of random series, in double double precision.

    p : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    d : constant natural32 := natural32(deg);

    use DoblDobl_Speelpenning_Convolutions;

  begin
    put_line("Reading a polynomial system ..."); get(p);
    declare
      s : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := System_to_Series_System(p.all);
      dim : constant integer32 := p'last;
      c : constant Circuits(p'range) := Make_Convolution_Circuits(p.all,d);
      q : Link_to_System := Create(c,dim,deg);
      x : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
        := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
      jp : constant DoblDobl_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
         := DoblDobl_CSeries_Jaco_Matrices.Create(s);
      jm : constant DoblDobl_Complex_Series_Matrices.Matrix(1..dim,1..dim)
         := DoblDobl_CSeries_Jaco_Matrices.Eval(jp,x);
      sx : constant DoblDobl_Complex_Series_Vectors.Vector(p'range)
         := DoblDobl_CSeries_Poly_SysFun.Eval(s,x);
      xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
           := DoblDobl_Series_Coefficients(x);
      err : double_double;
    begin
      Compute(q.pwt,q.mxe,xcff);
      EvalDiff(q,xcff);
      err := Difference(sx,q.yv);
      put("The evaluation error : "); put(err,3); new_line;
      err := Difference(jm,q.vm);
      put("The differentiation error : "); put(err,3); new_line;
      Clear(q);
    end;
  end DoblDobl_Input_Test;

  procedure QuadDobl_Input_Test ( deg : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   makes convolution circuits, evaluates and differentiates
  --   at a vector of random series, in quad double precision.

    p : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    d : constant natural32 := natural32(deg);

    use QuadDobl_Speelpenning_Convolutions;

  begin
    put_line("Reading a polynomial system ..."); get(p);
    declare
      s : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := System_to_Series_System(p.all);
      dim : constant integer32 := p'last;
      c : constant Circuits(p'range) := Make_Convolution_Circuits(p.all,d);
      q : Link_to_System := Create(c,dim,deg);
      x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
        := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
      jp : constant QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
         := QuadDobl_CSeries_Jaco_Matrices.Create(s);
      jm : constant QuadDobl_Complex_Series_Matrices.Matrix(1..dim,1..dim)
         := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
      sx : constant QuadDobl_Complex_Series_Vectors.Vector(p'range)
         := QuadDobl_CSeries_Poly_SysFun.Eval(s,x);
      xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
           := QuadDobl_Series_Coefficients(x);
      err : quad_double;
    begin
      Compute(q.pwt,q.mxe,xcff);
      EvalDiff(q,xcff);
      err := Difference(sx,q.yv);
      put("The evaluation error : "); put(err,3); new_line;
      err := Difference(jm,q.vm);
      put("The differentiation error : "); put(err,3); new_line;
      Clear(q);
    end;
  end QuadDobl_Input_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree, the dimension,
  --   the number of monomials, and the precision.  Then runs the tests.

    dim,deg,nbr,pwr : integer32 := 0;
    precision,random,answer : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    new_line;
    put("Random polynomials ? (y/n) "); Ask_Yes_Or_No(random);
    if random = 'n' then
      new_line;
      put("Give the degree of the series : "); get(deg);
      new_line;
      case precision is
        when '0' => Standard_Input_Test(deg);
        when '1' => DoblDobl_Input_Test(deg);
        when '2' => QuadDobl_Input_Test(deg);
        when others => null;
      end case;
    else
      new_line;
      put("Give the dimension : "); get(dim);
      put("Give the degree of the series : "); get(deg);
      put("Give the number of monomials : "); get(nbr);
      put("Give the highest power of each variable : "); get(pwr);
      new_line;
      put("Test system ? (y/n) "); Ask_Yes_or_No(answer);
      new_line;
      if answer = 'y' then
        case precision is
          when '0' => Standard_System_Test(dim,deg,nbr,pwr);
          when '1' => DoblDobl_System_Test(dim,deg,nbr,pwr);
          when '2' => QuadDobl_System_Test(dim,deg,nbr,pwr);
          when others => null;
        end case;
      else
        case precision is
          when '0' => Standard_Test(dim,deg,nbr,pwr);
          when '1' => DoblDobl_Test(dim,deg,nbr,pwr);
          when '2' => QuadDobl_Test(dim,deg,nbr,pwr);
          when others => null;
        end case;
      end if;
    end if;
  end Main;

begin
  Main;
end ts_speelcnv;
