with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Double_Double_Numbers_io;            use Double_Double_Numbers_io;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Quad_Double_Numbers_io;              use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Random_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with Exponent_Indices;
with Standard_Dense_Series;
with Standard_Dense_Series_Vectors;
with DoblDobl_Dense_Series;
with DoblDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series;
with QuadDobl_Dense_Series_Vectors;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Functions;
with Standard_Series_Poly_Systems;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Functions;
with DoblDobl_Series_Poly_Systems;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Functions;
with QuadDobl_Series_Poly_Systems;
with Standard_Random_Series;
with DoblDobl_Random_Series;
with QuadDobl_Random_Series;
with Series_and_Polynomials_io;           use Series_and_Polynomials_io;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Series_Polynomial_Gradients;         use Series_Polynomial_Gradients;

procedure ts_speelcnv is

-- DESCRIPTION :
--   Tests the evaluation of the gradient of a polynomial in many variables,
--   in a power series of some fixed degree.

  function Random_Exponents
             ( dim,nbr : integer32; pwr : integer32 := 1 ) 
             return Standard_Integer_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Generates as many random vectors of dimension dim
  --   as the value of nbr.  All vectors have a nonzero sum.
  --   The higher power equals pow.

    res : Standard_Integer_VecVecs.VecVec(1..nbr);
    nz : integer32;

  begin
    for i in 1..nbr loop
      loop
        declare
          xp : constant Standard_Integer_Vectors.Vector(1..dim)
             := Standard_Random_Vectors.Random_Vector(1,dim,0,pwr);
        begin
          nz := Standard_Integer_Vectors.Sum(xp);
          if nz > 0
           then res(i) := new Standard_Integer_Vectors.Vector'(xp);
          end if;
        end;
        exit when (nz > 0);
      end loop;
    end loop;
    return res;
  end Random_Exponents;

  function Difference ( s : Standard_Dense_Series.Series;
                        c : Standard_Complex_Vectors.Link_to_Vector )
                      return double_float is

  -- DESCRIPTION :
  --   Returns sum of the absolute value of the difference
  --   between the coefficients of s and the values in c.

    use Standard_Complex_Numbers;

    res : double_float := 0.0;
    avl : double_float;
    val : Complex_Number;

  begin
    for i in c'range loop
      val := s.cff(i) - c(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

  function Difference ( s : DoblDobl_Dense_Series.Series;
                        c : DoblDobl_Complex_Vectors.Link_to_Vector )
                      return double_double is

  -- DESCRIPTION :
  --   Returns sum of the absolute value of the difference
  --   between the coefficients of s and the values in c.

    use DoblDobl_Complex_Numbers;

    res : double_double := create(0.0);
    avl : double_double;
    val : Complex_Number;

  begin
    for i in c'range loop
      val := s.cff(i) - c(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

  function Difference ( s : QuadDobl_Dense_Series.Series;
                        c : QuadDobl_Complex_Vectors.Link_to_Vector )
                      return quad_double is

  -- DESCRIPTION :
  --   Returns sum of the absolute value of the difference
  --   between the coefficients of s and the values in c.

    use QuadDobl_Complex_Numbers;

    res : quad_double := create(0.0);
    avl : quad_double;
    val : Complex_Number;

  begin
    for i in c'range loop
      val := s.cff(i) - c(i);
      avl := AbsVal(val);
      res := res + avl;
    end loop;
    return res;
  end Difference;

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
    polcff : constant Standard_Dense_Series_Vectors.Vector(1..nbr)
           := Standard_Random_Series.Random_Series_Vector(1,nbr,deg);
    pol : constant Standard_Series_Polynomials.Poly
       -- := Standard_Polynomial(dim,deg,idx); -- all coefficients are one
       -- := Standard_Polynomial(dim,idx,polcff); -- all exponents are one
        := Standard_Polynomial(dim,xps,polcff,false);
    x : constant Standard_Dense_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series.Random_Series_Vector(1,dim,deg);
    y : Standard_Dense_Series.Series;
    grad : Standard_Dense_Series_Vectors.Vector(1..dim);
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
    y := Standard_Series_Poly_Functions.Eval(pol,x);
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
    polcff : constant DoblDobl_Dense_Series_Vectors.Vector(1..nbr)
           := DoblDobl_Random_Series.Random_Series_Vector(1,nbr,deg);
    pol : constant DoblDobl_Series_Polynomials.Poly
       -- := DoblDobl_Polynomial(dim,deg,idx); -- all coefficients are one
       -- := DoblDobl_Polynomial(dim,idx,polcff); -- all exponents are one
        := DoblDobl_Polynomial(dim,xps,polcff,false);
    x : constant DoblDobl_Dense_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    y : DoblDobl_Dense_Series.Series;
    grad : DoblDobl_Dense_Series_Vectors.Vector(1..dim);
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
    y := DoblDobl_Series_Poly_Functions.Eval(pol,x);
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
    polcff : constant QuadDobl_Dense_Series_Vectors.Vector(1..nbr)
           := QuadDobl_Random_Series.Random_Series_Vector(1,nbr,deg);
    pol : constant QuadDobl_Series_Polynomials.Poly
       -- := QuadDobl_Polynomial(dim,deg,idx); -- all coefficients are one
       -- := QuadDobl_Polynomial(dim,idx,polcff); -- all exponents are one
        := QuadDobl_Polynomial(dim,xps,polcff,false);
    x : constant QuadDobl_Dense_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series.Random_Series_Vector(1,dim,deg);
    y : QuadDobl_Dense_Series.Series;
    grad : QuadDobl_Dense_Series_Vectors.Vector(1..dim);
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
    y := QuadDobl_Series_Poly_Functions.Eval(pol,x);
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

  function Standard_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return Standard_Speelpenning_Convolutions.Convolution_Circuit is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in standard double precision and returns a convolution circuit.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Speelpenning_Convolutions;

    res : Convolution_Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant Standard_Dense_Series_Vectors.Vector(1..nbr)
           := Standard_Random_Series.Random_Series_Vector(1,nbr,deg);
    rancst : constant Standard_Dense_Series.Series
           := Standard_Random_Series.Random_Series(deg);
    cstcff : constant Standard_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := Standard_Series_Coefficients(polcff);
    res.cst := new Standard_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end Standard_Random_Convolution_Circuit;

  function DoblDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuit is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in double double precision and returns a convolution circuit.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Speelpenning_Convolutions;

    res : Convolution_Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant DoblDobl_Dense_Series_Vectors.Vector(1..nbr)
           := DoblDobl_Random_Series.Random_Series_Vector(1,nbr,deg);
    rancst : constant DoblDobl_Dense_Series.Series
           := DoblDobl_Random_Series.Random_Series(deg);
    cstcff : constant DoblDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := DoblDobl_Series_Coefficients(polcff);
    res.cst := new DoblDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end DoblDobl_Random_Convolution_Circuit;

  function QuadDobl_Random_Convolution_Circuit
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuit is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in quad double precision and returns a convolution circuit.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Speelpenning_Convolutions;

    res : Convolution_Circuit(nbr,dim,dim-1,dim-2);
    polcff : constant QuadDobl_Dense_Series_Vectors.Vector(1..nbr)
           := QuadDobl_Random_Series.Random_Series_Vector(1,nbr,deg);
    rancst : constant QuadDobl_Dense_Series.Series
           := QuadDobl_Random_Series.Random_Series(deg);
    cstcff : constant QuadDobl_Complex_Vectors.Vector(0..rancst.deg)
           := rancst.cff(0..rancst.deg);

  begin
    res.xps := Random_Exponents(dim,nbr,pwr);
    res.idx := Exponent_Indices.Exponent_Index(res.xps);
    res.fac := Exponent_Indices.Factor_Index(res.xps);
    res.cff := QuadDobl_Series_Coefficients(polcff);
    res.cst := new QuadDobl_Complex_Vectors.Vector'(cstcff);
    res.forward := Allocate_Coefficients(dim-1,deg);
    res.backward := Allocate_Coefficients(dim-2,deg);
    res.cross := Allocate_Coefficients(dim-2,deg);
    res.wrk := Allocate_Coefficients(deg);
    res.acc := Allocate_Coefficients(deg);
    return res;
  end QuadDobl_Random_Convolution_Circuit;

  function Standard_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return Standard_Speelpenning_Convolutions.Convolution_Circuits is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in standard double precision and returns as many convolution circuits
  --   as the value of dim.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use Standard_Speelpenning_Convolutions;

    res : Convolution_Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Convolution_Circuit(nbr,dim,dim-1,dim-2)
          := Standard_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Convolution_Circuit'(c);
      end;
    end loop;
    return res;
  end Standard_Random_Convolution_Circuits;

  function DoblDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return DoblDobl_Speelpenning_Convolutions.Convolution_Circuits is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in double double precision and returns as many convolution circuits
  --   as the value of dim.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use DoblDobl_Speelpenning_Convolutions;

    res : Convolution_Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Convolution_Circuit(nbr,dim,dim-1,dim-2)
          := DoblDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Convolution_Circuit'(c);
      end;
    end loop;
    return res;
  end DoblDobl_Random_Convolution_Circuits;

  function QuadDobl_Random_Convolution_Circuits
             ( dim,deg,nbr,pwr : in integer32 )
             return QuadDobl_Speelpenning_Convolutions.Convolution_Circuits is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and coefficients
  --   in quad double precision and returns as many convolution circuits
  --   as the value of dim.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

    use QuadDobl_Speelpenning_Convolutions;

    res : Convolution_Circuits(1..dim);

  begin
    for k in 1..dim loop
      declare
        c : constant Convolution_Circuit(nbr,dim,dim-1,dim-2)
          := QuadDobl_Random_Convolution_Circuit(dim,deg,nbr,pwr);
      begin
        res(k) := new Convolution_Circuit'(c);
      end;
    end loop;
    return res;
  end QuadDobl_Random_Convolution_Circuits;

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

    c : constant Convolution_Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    p : constant Standard_Series_Poly_Systems.Poly_Sys(1..dim)
      := Standard_System(c);

  begin
    put_line("the polynomial system :"); put(p,dim+1);
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

    c : constant Convolution_Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    p : constant DoblDobl_Series_Poly_Systems.Poly_Sys(1..dim)
      := DoblDobl_System(c);

  begin
    put_line("the polynomial system :"); put(p,dim+1);
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

    c : constant Convolution_Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    p : constant QuadDobl_Series_Poly_Systems.Poly_Sys(1..dim)
      := QuadDobl_System(c);

  begin
    put_line("the polynomial system :"); put(p,dim+1);
  end QuadDobl_System_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree, the dimension,
  --   the number of monomials, and the precision.  Then runs the tests.

    dim,deg,nbr,pwr : integer32 := 0;
    precision,answer : character;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the series : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    put("Give the highest power of each variable : "); get(pwr);
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
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
  end Main;

begin
  Main;
end ts_speelcnv;
