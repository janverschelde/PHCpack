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
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
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
with Standard_Vector_Splitters;
with Standard_Coefficient_Convolutions;

procedure ts_speelcnv is

-- DESCRIPTION :
--   Tests the evaluation of the gradient of a polynomial in many variables,
--   in a power series of some fixed degree.

  function Leading_Coefficients
             ( s : Standard_Complex_Series_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector of leading coefficients of the vector s.

    res : Standard_Complex_Vectors.Vector(s'range);

  begin
    for i in s'range loop
      res(i) := s(i).cff(0);
    end loop;
    return res;
  end Leading_Coefficients;

  function Leading_Coefficients
             ( s : DoblDobl_Complex_Series_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector of leading coefficients of the vector s.

    res : DoblDobl_Complex_Vectors.Vector(s'range);

  begin
    for i in s'range loop
      res(i) := s(i).cff(0);
    end loop;
    return res;
  end Leading_Coefficients;

  function Leading_Coefficients
             ( s : QuadDobl_Complex_Series_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector of leading coefficients of the vector s.

    res : QuadDobl_Complex_Vectors.Vector(s'range);

  begin
    for i in s'range loop
      res(i) := s(i).cff(0);
    end loop;
    return res;
  end Leading_Coefficients;

  function One_Coefficients
             ( nbr,deg : integer32 )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns a vector of range 1..nbr with vectors of range 0..deg
  --   to represent the constant one as the coefficients.

    res : Standard_Complex_VecVecs.VecVec(1..nbr);

  begin
    for k in 1..nbr loop
      declare
        cff : Standard_Complex_Vectors.Vector(0..deg)
            := (0..deg => Standard_Complex_Numbers.Create(0.0));
      begin
        cff(0) := Standard_Complex_Numbers.Create(1.0);
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end One_Coefficients;

  function Standard_Make_Polynomial
             ( dim,deg : integer32;
               idx : Standard_Integer_VecVecs.VecVec;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : Standard_Complex_Series_Vectors.Vector;
               expone,cffone : boolean )
             return Standard_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Wraps the construction of a polynomial with series coefficients,
  --   with flags for special cases if expone and cffone.

  -- ON ENTRY :
  --   dim     dimension of the exponent vectors, number of variables;
  --   deg     degree of the power series;
  --   idx     exponent indices;
  --   xps     exponent vectors;
  --   cff     coefficients of the monomials;
  --   expone  true if all exponents are equal to one;
  --   cffone  true if all coefficients are equal to one.

    res : Standard_CSeries_Polynomials.Poly;

  begin
    if expone then -- all exponents are equal to one
      if cffone then -- all coefficients are equal to one
        res := Standard_Polynomial(dim,deg,idx);
      else
        res := Standard_Polynomial(dim,idx,cff,true);
      end if;
    else  
      if cffone then -- all coefficients are equal to one
        res := Standard_Polynomial(dim,deg,xps,false);
      else
        res := Standard_Polynomial(dim,xps,cff,false);
      end if;
    end if;
    return res;
  end Standard_Make_Polynomial;

  procedure Standard_Test ( dim,deg,nbr,pwr : in integer32;
                            expone,cffone : in boolean ) is

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in standard double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.
  --   expone   true if all exponents are equal to one;
  --   cffone   true if all coefficients are equal to one.

    use Standard_Speelpenning_Convolutions;
    use Standard_Vector_Splitters;
    use Standard_Coefficient_Convolutions;

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
        := Standard_Make_Polynomial(dim,deg,idx,xps,polcff,expone,cffone);
    x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xpt : constant Standard_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
    y : Standard_Complex_Series.Link_to_Series;
    ypt : Standard_Complex_Numbers.Complex_Number;
    grad : Standard_Complex_Series_Vectors.Vector(1..dim);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Allocate_Floating_Coefficients(dim,deg);
    pcff : Standard_Complex_VecVecs.VecVec(1..nbr);
    rpcf : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(nbr,deg);
    ipcf : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(nbr,deg);
    forward : constant Standard_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    rfwd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-1,deg);
    ifwd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-1,deg);
    backward : constant Standard_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    rbck : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,deg);
    ibck : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,deg);
    cross : constant Standard_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    rcrs : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,deg);
    icrs : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,deg);
    ygrad : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    ygrad2 : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
           := Allocate_Coefficients(dim+1,deg);
    ryd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,deg);
    iyd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,deg);
    work : constant Standard_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);
    rwrk : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(deg);
    iwrk : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(deg);
    acc : constant Standard_Complex_Vectors.Link_to_Vector
        := Allocate_Coefficients(deg);
    racc : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(deg);
    iacc : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(deg);
    err,err2,sumerr,sumerr2 : double_float := 0.0;
    pwt : Standard_Speelpenning_Convolutions.Link_to_VecVecVec
        := Standard_Speelpenning_Convolutions.Create(xcff,mxe);
    rpwt : Standard_Coefficient_Convolutions.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    ipwt : Standard_Coefficient_Convolutions.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    crc : Standard_Speelpenning_Convolutions.Circuit(nbr,dim,dim+1,dim+2);

    use Standard_Complex_Vectors;

  begin
    if cffone 
     then pcff := One_Coefficients(nbr,deg);
     else pcff := Standard_Series_Coefficients(polcff);
    end if;
    crc.xps := xps; crc.idx := idx; crc.fac := fac; crc.cff := pcff;
    put_line("Some random exponents :"); Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :"); Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :"); Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := Standard_CSeries_Poly_Functions.Eval(pol,x);
    if expone and cffone then
      Speel(idx,xcff,forward,backward,cross,ygrad); -- all coefficients one
      Standard_Vector_Splitters.Complex_Parts(xcff,rx,ix);
      Speel(idx,rx.all,ix.all,rfwd.all,ifwd.all,rbck.all,ibck.all,
            rcrs.all,icrs.all,ryd.all,iyd.all);
      Standard_Vector_Splitters.Complex_Merge(ryd,iyd,ygrad2);
    elsif expone then
      Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
      Standard_Vector_Splitters.Complex_Parts(pcff,rpcf,ipcf);
      Standard_Vector_Splitters.Complex_Parts(xcff,rx,ix);
      Speel(idx,rpcf.all,ipcf.all,rx.all,ix.all,rfwd.all,ifwd.all,rbck.all,
            ibck.all,rcrs.all,icrs.all,ryd.all,iyd.all,rwrk,iwrk);
      Standard_Vector_Splitters.Complex_Merge(ryd,iyd,ygrad2);
    else
      Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
      Standard_Vector_Splitters.Complex_Parts(pcff,rpcf,ipcf);
      Standard_Vector_Splitters.Complex_Parts(xcff,rx,ix);
      Compute(rpwt,ipwt,mxe,rx,ix);
      Speel(xps,idx,fac,
            rpcf.all,ipcf.all,rx.all,ix.all,rfwd.all,ifwd.all,
            rbck.all,ibck.all,rcrs.all,icrs.all,ryd.all,iyd.all,
            rwrk,iwrk,racc,iacc,rpwt,ipwt);
      Standard_Vector_Splitters.Complex_Merge(ryd,iyd,ygrad2);
    end if;
    put_line("The value of the polynomial at the random series :");
    put(y); new_line;
    ypt := Eval(crc,xpt);
    put_line("The leading coefficient of the evaluated polynomial :");
    put(ypt); new_line;
    for i in ygrad'range loop
      ygrad(i)(0) := Standard_Complex_Numbers.Create(0.0);
    end loop;
    work(0) := Standard_Complex_Numbers.Create(0.0);
    acc(0) := Standard_Complex_Numbers.Create(0.0);
    if expone then
      Speel(idx,pcff,xpt,forward,backward,cross,ygrad,work); -- all powers 1
    else
      Speel(xps,idx,fac,pcff,xpt,forward,backward,cross,ygrad,work,acc,pwt);
    end if;
    put_line("The leading coefficients computed in reverse mode :");
    for i in ygrad'range loop
      put(ygrad(i)(0)); new_line;
    end loop;
    put_line("The coefficient vector of the value of the polynomial :");
    put_line(ygrad(ygrad'last));
    put_line("Recomputed coefficient vector of the value of the polynomial :");
    put_line(ygrad2(ygrad2'last));
    grad := Standard_Gradient(pol,x);
    err := Difference(y,ygrad(ygrad'last));
    put("The error :"); put(err,3); new_line;
    err2 := Difference(y,ygrad2(ygrad2'last));
    put("Recomputed error :"); put(err2,3); new_line;
    sumerr := err; sumerr2 := err2;
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :"); put(grad(k)); new_line;
      if ygrad(k) /= null then
        put("The coefficient vector of derivative ");
        put(k,1); put_line(" :"); put_line(ygrad(k));
        put("Recomputed coefficient vector of derivative ");
        put(k,1); put_line(" :"); put_line(ygrad2(k));
        err := Difference(grad(k),ygrad(k));
        put("The error :"); put(err,3); new_line;
        err2 := Difference(grad(k),ygrad2(k));
        put("Recomputed error :"); put(err2,3); new_line;
        sumerr := sumerr + err; sumerr2 := sumerr2 + err2;
      end if;
    end loop;
    put("Sum of errors :"); put(sumerr,3); new_line;
    put("Recomputed sum of errors :"); put(sumerr2,3); new_line;
    Standard_Speelpenning_Convolutions.Clear(pwt);
    Standard_Coefficient_Convolutions.Clear(rpwt);
    Standard_Coefficient_Convolutions.Clear(ipwt);
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
    xpt : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
    y : DoblDobl_Complex_Series.Link_to_Series;
    ypt : DoblDobl_Complex_Numbers.Complex_Number;
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
    crc : Circuit(nbr,dim,dim+1,dim+2);

  begin
    crc.xps := xps; crc.idx := idx; crc.fac := fac; crc.cff := pcff;
    put_line("Some random exponents :"); Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :"); Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :"); Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := DoblDobl_CSeries_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
   -- Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
    Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
    put_line("The value of the polynomial at the random series :");
    put(y); new_line;
    ypt := Eval(crc,xpt);
    put_line("The leading coefficient of the evaluated polynomial :");
    put(ypt); new_line;
    for i in ygrad'range loop
      ygrad(i)(0) := DoblDobl_Complex_Numbers.Create(integer(0));
    end loop;
    work(0) := DoblDobl_Complex_Numbers.Create(integer(0));
    acc(0) := DoblDobl_Complex_Numbers.Create(integer(0));
   -- Speel(idx,pcff,xpt,forward,backward,cross,ygrad,work); -- all powers 1
    Speel(xps,idx,fac,pcff,xpt,forward,backward,cross,ygrad,work,acc,pwt);
    put_line("The leading coefficients computed in reverse mode :");
    for i in ygrad'range loop
      put(ygrad(i)(0)); new_line;
    end loop;
    put_line("The coefficient vector of the value of the polynomial :");
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
    xpt : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
    y : QuadDobl_Complex_Series.Link_to_Series;
    ypt : QuadDobl_Complex_Numbers.Complex_Number;
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
    crc : Circuit(nbr,dim,dim+1,dim+2);

  begin
    crc.xps := xps; crc.idx := idx; crc.fac := fac; crc.cff := pcff;
    put_line("Some random exponents :"); Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :"); Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :"); Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := QuadDobl_CSeries_Poly_Functions.Eval(pol,x);
   -- Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
   -- Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
    Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
    put_line("The value of the product at the random series :");
    put(y); new_line;
    ypt := Eval(crc,xpt);
    put_line("The leading coefficient of the evaluated polynomial :");
    put(ypt); new_line;
    for i in ygrad'range loop
      ygrad(i)(0) := QuadDobl_Complex_Numbers.Create(integer(0));
    end loop;
    work(0) := QuadDobl_Complex_Numbers.Create(integer(0));
    acc(0) := QuadDobl_Complex_Numbers.Create(integer(0));
   -- Speel(idx,pcff,xpt,forward,backward,cross,ygrad,work); -- all powers 1
    Speel(xps,idx,fac,pcff,xpt,forward,backward,cross,ygrad,work,acc,pwt);
    put_line("The leading coefficients computed in reverse mode :");
    for i in ygrad'range loop
      put(ygrad(i)(0)); new_line;
    end loop;
    put_line("The coefficient vector of the value of the polynomial :");
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
    xpt : constant Standard_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
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
      put("The matrix "); put(i,1); put_line(" :"); put(s.vm(i).all);
    end loop;
    for i in 1..dim loop
      for j in 1..dim loop
        put("the series of the Jacobian at ");
        put(i,1); put(" and "); put(j,1); put_line(" :"); put(jm(i,j));
      end loop;
    end loop;
    err := Difference(jm,s.vm);
    put("The error :"); put(err,3); new_line;
    Compute(s.pwt,s.mxe,xpt);
    EvalDiff(s,xpt);
    put_line("The evaluation at a point : ");
    for k in s.yv'range loop
      put(s.yv(k)(0)); new_line;
      put(px(k).cff(0)); new_line;
    end loop;
    put_line("The Jacobian matrix at a point : "); put(s.vm(0).all);
    put_line("The leading coefficients of the Jacobian :");
    for i in 1..dim loop
      for j in 1..dim loop
        put(" "); put(jm(i,j).cff(0));
      end loop;
      new_line;
    end loop;
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
    xpt : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
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
      put("The matrix "); put(i,1); put_line(" :"); put(s.vm(i).all);
    end loop;
    for i in 1..dim loop
      for j in 1..dim loop
        put("the series of the Jacobian at ");
        put(i,1); put(" and "); put(j,1); put_line(" :"); put(jm(i,j));
      end loop;
    end loop;
    err := Difference(jm,s.vm);
    put("The error : "); put(err,3); new_line;
    Compute(s.pwt,s.mxe,xpt);
    EvalDiff(s,xpt);
    put_line("The evaluation at a point : ");
    for k in s.yv'range loop
      put(s.yv(k)(0)); new_line;
      put(px(k).cff(0)); new_line;
    end loop;
    put_line("The Jacobian matrix at a point : "); put(s.vm(0).all);
    put_line("The leading coefficients of the Jacobian :");
    for i in 1..dim loop
      for j in 1..dim loop
        put(" "); put(jm(i,j).cff(0));
      end loop;
      new_line;
    end loop;
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
    xpt : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
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
      put("The matrix "); put(i,1); put_line(" :"); put(s.vm(i).all);
    end loop;
    for i in 1..dim loop
      for j in 1..dim loop
        put("the series of the Jacobian at ");
        put(i,1); put(" and "); put(j,1); put_line(" :"); put(jm(i,j));
      end loop;
    end loop;
    err := Difference(jm,s.vm);
    put("The error : "); put(err,3); new_line;
    Compute(s.pwt,s.mxe,xpt);
    EvalDiff(s,xpt);
    put_line("The evaluation at a point : ");
    for k in s.yv'range loop
      put(s.yv(k)(0)); new_line;
      put(px(k).cff(0)); new_line;
    end loop;
    put_line("The Jacobian matrix at a point : "); put(s.vm(0).all);
    put_line("The leading coefficients of the Jacobian :");
    for i in 1..dim loop
      for j in 1..dim loop
        put(" "); put(jm(i,j).cff(0));
      end loop;
      new_line;
    end loop;
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
  --   Prompts the user for the precision, degree, the dimension,
  --   the number of monomials.  Then runs the tests.

    dim,deg,nbr,pwr : integer32 := 0;
    precision : constant character := Prompt_for_Precision;
    random,answer : character;
    expone,cffone : boolean;

  begin
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
      expone := (pwr = 1);
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
        put("All coefficients equal to one ? (y/n) ");
        Ask_Yes_or_No(answer); cffone := (answer = 'y');
        case precision is
          when '0' => Standard_Test(dim,deg,nbr,pwr,expone,cffone);
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
