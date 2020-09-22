with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Exponent_Indices;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_io;         use QuadDobl_Complex_Series_io;
with QuadDobl_Complex_Series_Vectors_io; use QuadDobl_Complex_Series_Vectors_io;
with QuadDobl_Complex_Series_Matrices;
with QuadDobl_Random_Series_Vectors;
with QuadDobl_CSeries_Polynomials_io;    use QuadDobl_CSeries_Polynomials_io;
with QuadDobl_CSeries_Poly_Functions;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems_io;
with QuadDobl_CSeries_Poly_SysFun;
with QuadDobl_CSeries_Jaco_Matrices;
with QuadDobl_Speelpenning_Convolutions;
with Series_Coefficient_Vectors;         use Series_Coefficient_Vectors;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Complex_Series_and_Polynomials;     use Complex_Series_and_Polynomials;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Standard_Vector_Splitters;
with QuadDobl_Vector_Splitters;
with Standard_Coefficient_Convolutions;
with QuadDobl_Coefficient_Convolutions;
with QuadDobl_Convolution_Splitters;

package body Test_QuadDobl_Coeff_Convolutions is

  function Leading_Coefficients
             ( s : QuadDobl_Complex_Series_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(s'range);

  begin
    for i in s'range loop
      res(i) := s(i).cff(0);
    end loop;
    return res;
  end Leading_Coefficients;

  function One_Coefficients
             ( nbr,deg : integer32 )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(1..nbr);

  begin
    for k in 1..nbr loop
      declare
        cff : QuadDobl_Complex_Vectors.Vector(0..deg)
            := (0..deg => QuadDobl_Complex_Numbers.Create(integer(0)));
      begin
        cff(0) := QuadDobl_Complex_Numbers.Create(integer(1.0));
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end One_Coefficients;

  function QuadDobl_Make_Polynomial
             ( dim,deg : integer32;
               idx : Standard_Integer_VecVecs.VecVec;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : QuadDobl_Complex_Series_Vectors.Vector;
               expone,cffone : boolean )
             return QuadDobl_CSeries_Polynomials.Poly is

    res : QuadDobl_CSeries_Polynomials.Poly;

  begin
    if expone then -- all exponents are equal to one
      if cffone then -- all coefficients are equal to one
        res := QuadDobl_Polynomial(dim,deg,idx);
      else
        res := QuadDobl_Polynomial(dim,idx,cff,true);
      end if;
    else  
      if cffone then -- all coefficients are equal to one
        res := QuadDobl_Polynomial(dim,deg,xps,false);
      else
        res := QuadDobl_Polynomial(dim,xps,cff,false);
      end if;
    end if;
    return res;
  end QuadDobl_Make_Polynomial;

  procedure QuadDobl_Test ( dim,deg,nbr,pwr : in integer32;
                            expone,cffone : in boolean ) is

    use Standard_Vector_Splitters;
    use QuadDobl_Speelpenning_Convolutions;
    use QuadDobl_Coefficient_Convolutions;

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
        := QuadDobl_Make_Polynomial(dim,deg,idx,xps,polcff,expone,cffone);
    x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xpt : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
    y : QuadDobl_Complex_Series.Link_to_Series;
    ypt : QuadDobl_Complex_Numbers.Complex_Number;
    grad : QuadDobl_Complex_Series_Vectors.Vector(1..dim);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    degdim : constant integer32 := 4*(deg+1)-1;
    xr : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Allocate_Floating_Coefficients(dim,degdim);
    xi : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Allocate_Floating_Coefficients(dim,degdim);
    pcff : QuadDobl_Complex_VecVecs.VecVec(1..nbr);
    rpcf : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(nbr,degdim);
    ipcf : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(nbr,degdim);
    forward : constant QuadDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    rfwd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-1,degdim);
    ifwd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-1,degdim);
    backward : constant QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    rbck : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,degdim);
    ibck : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,degdim);
    cross : constant QuadDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    rcrs : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,degdim);
    icrs : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim-2,degdim);
    ygrad : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    ygrad2 : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
           := Allocate_Coefficients(dim+1,deg);
    ryd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,degdim);
    iyd : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim+1,degdim);
    work : constant QuadDobl_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);
    rwrk : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(degdim);
    iwrk : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(degdim);
    acc : constant QuadDobl_Complex_Vectors.Link_to_Vector
        := Allocate_Coefficients(deg);
    racc : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(degdim);
    iacc : constant Standard_Floating_Vectors.Link_to_Vector
         := Allocate_Floating_Coefficients(degdim);
    err,sumerr : quad_double;
    pwt : Link_to_VecVecVec := Create(xcff,mxe);
    rpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(mxe,degdim);
    ipwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
         := Standard_Coefficient_Convolutions.Allocate(mxe,degdim);
    crc : QuadDobl_Speelpenning_Convolutions.Circuit(nbr,dim,dim+1,dim+2);
    ur,vr,wr : constant Standard_Floating_Vectors.Vector(0..3) 
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(ur);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(vr);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(wr);

  begin
    if cffone 
     then pcff := One_Coefficients(nbr,deg);
     else pcff := QuadDobl_Series_Coefficients(polcff);
    end if;
    crc.xps := xps; crc.idx := idx; crc.fac := fac; crc.cff := pcff;
    put_line("Some random exponents :"); Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :"); Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :"); Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := QuadDobl_CSeries_Poly_Functions.Eval(pol,x);
    if expone and cffone then
      Speel(idx,xcff,forward,backward,cross,ygrad); -- if all coefficients one
      QuadDobl_Vector_Splitters.Complex_Parts(xcff,xr,xi);
      Speel(idx,xr.all,xi.all,rfwd.all,ifwd.all,rbck.all,ibck.all,
            rcrs.all,icrs.all,ryd.all,iyd.all,u,v,w);
      QuadDobl_Vector_Splitters.Complex_Merge(ryd,iyd,ygrad2);
    elsif expone then
      Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
      QuadDobl_Vector_Splitters.Complex_Parts(pcff,rpcf,ipcf);
      QuadDobl_Vector_Splitters.Complex_Parts(xcff,xr,xi);
      Speel(idx,rpcf.all,ipcf.all,xr.all,xi.all,rfwd.all,ifwd.all,rbck.all,
            ibck.all,rcrs.all,icrs.all,ryd.all,iyd.all,rwrk,iwrk,u,v,w);
      QuadDobl_Vector_Splitters.Complex_Merge(ryd,iyd,ygrad2);
    else
      Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
      QuadDobl_Vector_Splitters.Complex_Parts(pcff,rpcf,ipcf);
      QuadDobl_Vector_Splitters.Complex_Parts(xcff,xr,xi);
      Compute(rpwt,ipwt,mxe,xr,xi,u,v,w);
      Speel(xps,idx,fac,rpcf.all,ipcf.all,xr.all,xi.all,rfwd.all,ifwd.all,
            rbck.all,ibck.all,rcrs.all,icrs.all,ryd.all,iyd.all,
            rwrk,iwrk,racc,iacc,rpwt,ipwt,u,v,w);
      QuadDobl_Vector_Splitters.Complex_Merge(ryd,iyd,ygrad2);
    end if;
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

  procedure QuadDobl_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

    use QuadDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : Link_to_System := Create(c,dim,deg);
    cs : constant QuadDobl_Coefficient_Convolutions.Link_to_System
       := QuadDobl_Convolution_Splitters.Split(s);
    p : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..dim) := QuadDobl_System(c);
    x : QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xpt : constant QuadDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
    px : QuadDobl_Complex_Series_Vectors.Vector(p'range)
       := QuadDobl_CSeries_Poly_SysFun.Eval(p,x);
    xcff : QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    degdim : constant integer32 := 4*(deg+1)-1;
    xr : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,degdim);
    xi : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,degdim);
    jp : QuadDobl_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
       := QuadDobl_CSeries_Jaco_Matrices.Create(p);
    jm : QuadDobl_Complex_Series_Matrices.Matrix(1..dim,1..dim)
       := QuadDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    err : quad_double;
    ans : character;
    ur,vr,wr : constant Standard_Floating_Vectors.Vector(0..3)
             := (0..3 => 0.0);
    u : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(ur);
    v : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(vr);
    w : constant Standard_Floating_Vectors.Link_to_Vector
      := new Standard_Floating_Vectors.Vector'(wr);

  begin
    put_line("the polynomial system :");
    QuadDobl_CSeries_Poly_Systems_io.put(p);
    Compute(s.pwt,s.mxe,xcff);
    EvalDiff(s,xcff);
    put_line("The evaluation values :"); put_line(px);
    put_line("The coefficient vectors of the evaluation :"); put_line(s.yv);
    err := Difference(px,s.yv);
    put("The error : "); put(err,3); new_line;
    QuadDobl_Vector_Splitters.Complex_Parts(xcff,xr,xi);
    QuadDobl_Coefficient_Convolutions.Compute
      (cs.rpwt,cs.ipwt,cs.mxe,xr,xi,u,v,w);
    QuadDobl_Coefficient_Convolutions.EvalDiff(cs,xr.all,xi.all,u,v,w);
    err := Difference(px,cs.yv);
    put("Recomputed error : "); put(err,3); new_line;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
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
    err := Difference(jm,cs.vm);
    put("Recomputed error : "); put(err,3); new_line;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
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

  procedure QuadDobl_Input_Test ( deg : in integer32 ) is

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

    dim,deg,nbr,pwr : integer32 := 0;
    random,answer : character;
    expone,cffone : boolean;

  begin
    new_line;
    put("Random polynomials ? (y/n) "); Ask_Yes_Or_No(random);
    if random = 'n' then
      new_line;
      put("Give the degree of the series : "); get(deg);
      new_line;
      QuadDobl_Input_Test(deg);
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
        QuadDobl_System_Test(dim,deg,nbr,pwr);
      else
        put("All coefficients equal to one ? (y/n) ");
        Ask_Yes_or_No(answer); cffone := (answer = 'y');
        QuadDobl_Test(dim,deg,nbr,pwr,expone,cffone);
      end if;
    end if;
  end Main;

end Test_QuadDobl_Coeff_Convolutions;
