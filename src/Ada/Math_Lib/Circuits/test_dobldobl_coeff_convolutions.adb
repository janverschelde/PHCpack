with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with Exponent_Indices;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_io;         use DoblDobl_Complex_Series_io;
with DoblDobl_Complex_Series_Vectors_io; use DoblDobl_Complex_Series_Vectors_io;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_Random_Series_Vectors;
with DoblDobl_CSeries_Polynomials_io;    use DoblDobl_CSeries_Polynomials_io;
with DoblDobl_CSeries_Poly_Functions;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems_io;
with DoblDobl_CSeries_Poly_SysFun;
with DoblDobl_CSeries_Jaco_Matrices;
with DoblDobl_Speelpenning_Convolutions;
with Series_Coefficient_Vectors;         use Series_Coefficient_Vectors;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Complex_Series_and_Polynomials;     use Complex_Series_and_Polynomials;
with System_Convolution_Circuits;        use System_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Standard_Vector_Splitters;
with DoblDobl_Vector_Splitters;
with Standard_Coefficient_Convolutions;
with DoblDobl_Coefficient_Convolutions;
with DoblDobl_Convolution_Splitters;

package body Test_DoblDobl_Coeff_Convolutions is

  function Leading_Coefficients
             ( s : DoblDobl_Complex_Series_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(s'range);

  begin
    for i in s'range loop
      res(i) := s(i).cff(0);
    end loop;
    return res;
  end Leading_Coefficients;

  function One_Coefficients
             ( nbr,deg : integer32 )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(1..nbr);

  begin
    for k in 1..nbr loop
      declare
        cff : DoblDobl_Complex_Vectors.Vector(0..deg)
            := (0..deg => DoblDobl_Complex_Numbers.Create(integer(0)));
      begin
        cff(0) := DoblDobl_Complex_Numbers.Create(integer(1.0));
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end One_Coefficients;

  function DoblDobl_Make_Polynomial
             ( dim,deg : integer32;
               idx : Standard_Integer_VecVecs.VecVec;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : DoblDobl_Complex_Series_Vectors.Vector;
               expone,cffone : boolean )
             return DoblDobl_CSeries_Polynomials.Poly is

    res : DoblDobl_CSeries_Polynomials.Poly;

  begin
    if expone then -- all exponents are equal to one
      if cffone then -- all coefficients are equal to one
        res := DoblDobl_Polynomial(dim,deg,idx);
      else
        res := DoblDobl_Polynomial(dim,idx,cff,true);
      end if;
    else  
      if cffone then -- all coefficients are equal to one
        res := DoblDobl_Polynomial(dim,deg,xps,false);
      else
        res := DoblDobl_Polynomial(dim,xps,cff,false);
      end if;
    end if;
    return res;
  end DoblDobl_Make_Polynomial;

  procedure DoblDobl_Test ( dim,deg,nbr,pwr : in integer32;
                            expone,cffone : in boolean ) is

    use Standard_Vector_Splitters;
    use DoblDobl_Speelpenning_Convolutions;
    use DoblDobl_Coefficient_Convolutions;

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
        := DoblDobl_Make_Polynomial(dim,deg,idx,xps,polcff,expone,cffone);
    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xpt : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
    y : DoblDobl_Complex_Series.Link_to_Series;
    ypt : DoblDobl_Complex_Numbers.Complex_Number;
    grad : DoblDobl_Complex_Series_Vectors.Vector(1..dim);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    rhx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim,deg);
    ihx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim,deg);
    rlx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim,deg);
    ilx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Allocate_Floating_Coefficients(dim,deg);
    pcff : DoblDobl_Complex_VecVecs.VecVec(1..nbr);
    rhpcf : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(nbr,deg);
    ihpcf : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(nbr,deg);
    rlpcf : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(nbr,deg);
    ilpcf : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(nbr,deg);
    forward : constant DoblDobl_Complex_VecVecs.VecVec(1..dim-1)
            := Allocate_Coefficients(dim-1,deg);
    rhfwd : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-1,deg);
    ihfwd : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-1,deg);
    rlfwd : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-1,deg);
    ilfwd : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-1,deg);
    backward : constant DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
             := Allocate_Coefficients(dim-2,deg);
    rhbck : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    ihbck : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    rlbck : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    ilbck : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    cross : constant DoblDobl_Complex_VecVecs.VecVec(1..dim-2)
          := Allocate_Coefficients(dim-2,deg);
    rhcrs : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    ihcrs : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    rlcrs : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    ilcrs : constant Standard_Floating_VecVecs.Link_to_VecVec
          := Allocate_Floating_Coefficients(dim-2,deg);
    ygrad : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
          := Allocate_Coefficients(dim+1,deg);
    ygrad2 : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
           := Allocate_Coefficients(dim+1,deg);
    rhyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    ihyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    rlyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    ilyd : constant Standard_Floating_VecVecs.Link_to_VecVec
         := Allocate_Floating_Coefficients(dim+1,deg);
    work : constant DoblDobl_Complex_Vectors.Link_to_Vector
         := Allocate_Coefficients(deg);
    rhwrk : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    ihwrk : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    rlwrk : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    ilwrk : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    acc : constant DoblDobl_Complex_Vectors.Link_to_Vector
        := Allocate_Coefficients(deg);
    rhacc : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    ihacc : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    rlacc : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    ilacc : constant Standard_Floating_Vectors.Link_to_Vector
          := Allocate_Floating_Coefficients(deg);
    err,err2,sumerr,sumerr2 : double_double;
    pwt : DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec
        := DoblDobl_Speelpenning_Convolutions.Create(xcff,mxe);
    rhpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    ihpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    rlpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    ilpwt : constant Standard_Floating_VecVecVecs.Link_to_VecVecVec
          := Standard_Coefficient_Convolutions.Allocate(mxe,deg);
    crc : DoblDobl_Speelpenning_Convolutions.Circuit(nbr,dim,dim+1,dim+2);

    use DoblDobl_Complex_Vectors;

  begin
    if cffone 
     then pcff := One_Coefficients(nbr,deg);
     else pcff := DoblDobl_Series_Coefficients(polcff);
    end if;
    crc.xps := xps; crc.idx := idx; crc.fac := fac; crc.cff := pcff;
    put_line("Some random exponents :"); Standard_Integer_VecVecs_io.put(xps);
    put_line("its exponent indices :"); Standard_Integer_VecVecs_io.put(idx);
    put_line("its factor indices :"); Standard_Integer_VecVecs_io.put(fac);
    put("its maxima :"); Standard_Integer_Vectors_io.put(mxe); new_line;
    put_line("the polynomial :"); put(pol); new_line;
    y := DoblDobl_CSeries_Poly_Functions.Eval(pol,x);
    if expone and cffone then
      Speel(idx,xcff,forward,backward,cross,ygrad); -- all coefficients one
      DoblDobl_Vector_Splitters.Complex_Parts(xcff,rhx,ihx,rlx,ilx);
      Speel(idx,rhx.all,ihx.all,rlx.all,ilx.all,rhfwd.all,ihfwd.all,
            rlfwd.all,ilfwd.all,rhbck.all,ihbck.all,rlbck.all,ilbck.all,
            rhcrs.all,ihcrs.all,rlcrs.all,ilcrs.all,
            rhyd.all,ihyd.all,rlyd.all,ilyd.all);
      DoblDobl_Vector_Splitters.Complex_Merge(rhyd,ihyd,rlyd,ilyd,ygrad2);
    elsif expone then
      Speel(idx,pcff,xcff,forward,backward,cross,ygrad,work); -- all powers 1
      DoblDobl_Vector_Splitters.Complex_Parts(pcff,rhpcf,ihpcf,rlpcf,ilpcf);
      DoblDobl_Vector_Splitters.Complex_Parts(xcff,rhx,ihx,rlx,ilx);
      Speel(idx,rhpcf.all,ihpcf.all,rlpcf.all,ilpcf.all,
                rhx.all,ihx.all,rlx.all,ilx.all,
            rhfwd.all,ihfwd.all,rlfwd.all,ilfwd.all,
            rhbck.all,ihbck.all,rlbck.all,ilbck.all,
            rhcrs.all,ihcrs.all,rlcrs.all,ilcrs.all,
            rhyd.all,ihyd.all,rlyd.all,ilyd.all,rhwrk,ihwrk,rlwrk,ilwrk);
      DoblDobl_Vector_Splitters.Complex_Merge(rhyd,ihyd,rlyd,ilyd,ygrad2);
    else
      Speel(xps,idx,fac,pcff,xcff,forward,backward,cross,ygrad,work,acc,pwt);
      DoblDobl_Vector_Splitters.Complex_Parts(pcff,rhpcf,ihpcf,rlpcf,ilpcf);
      DoblDobl_Vector_Splitters.Complex_Parts(xcff,rhx,ihx,rlx,ilx);
      Compute(rhpwt,ihpwt,rlpwt,ilpwt,mxe,rhx,ihx,rlx,ilx);
      Speel(xps,idx,fac,rhpcf.all,ihpcf.all,rlpcf.all,ilpcf.all,
            rhx.all,ihx.all,rlx.all,ilx.all,
            rhfwd.all,ihfwd.all,rlfwd.all,ilfwd.all,
            rhbck.all,ihbck.all,rlbck.all,ilbck.all,
            rhcrs.all,ihcrs.all,rlcrs.all,ilcrs.all,
            rhyd.all,ihyd.all,rlyd.all,ilyd.all,
            rhwrk,ihwrk,rlwrk,ilwrk,rhacc,ihacc,rlacc,ilacc,
            rhpwt,ihpwt,rlpwt,ilpwt);
      DoblDobl_Vector_Splitters.Complex_Merge(rhyd,ihyd,rlyd,ilyd,ygrad2);
    end if;
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
    grad := DoblDobl_Gradient(pol,x);
    err := Difference(y,ygrad(ygrad'last));
    put("The error : "); put(err,3); new_line;
    err2 := Difference(y,ygrad2(ygrad2'last));
    put("Recomputed error : "); put(err2,3); new_line;
    sumerr := err; sumerr2 := err2;
    for k in grad'range loop
      put("derivative "); put(k,1); put_line(" :");
      put(grad(k)); new_line;
      if ygrad(k) /= null then
        put("The coefficient vector of derivative ");
        put(k,1); put_line(" :"); put_line(ygrad(k));
        put("Recomputed coefficient vector of derivative ");
        put(k,1); put_line(" :"); put_line(ygrad2(k));
        err := Difference(grad(k),ygrad(k));
        put("The error : "); put(err,3); new_line;
        err2 := Difference(grad(k),ygrad2(k));
        put("Recomputed error : "); put(err2,3); new_line;
        sumerr := sumerr + err; sumerr2 := sumerr2 + err2;
      end if;
    end loop;
    put("Sum of errors : "); put(sumerr,3); new_line;
    put("Recomputed sum of errors : "); put(sumerr2,3); new_line;
    Clear(pwt);
  end DoblDobl_Test;

  procedure DoblDobl_System_Test ( dim,deg,nbr,pwr : in integer32 ) is

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : Link_to_System := Create(c,dim,deg);
    cs : constant DoblDobl_Coefficient_Convolutions.Link_to_System
       := DoblDobl_Convolution_Splitters.Split(s);
    p : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..dim) := DoblDobl_System(c);
    x : DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xpt : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := Leading_Coefficients(x);
    px : DoblDobl_Complex_Series_Vectors.Vector(p'range)
       := DoblDobl_CSeries_Poly_SysFun.Eval(p,x);
    xcff : DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    rhx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ihx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rlx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ilx : constant Standard_Floating_VecVecs.Link_to_VecVec
        := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    jp : DoblDobl_CSeries_Jaco_Matrices.Jaco_Mat(1..dim,1..dim)
       := DoblDobl_CSeries_Jaco_Matrices.Create(p);
    jm : DoblDobl_Complex_Series_Matrices.Matrix(1..dim,1..dim)
       := DoblDobl_CSeries_Jaco_Matrices.Eval(jp,x);
    err : double_double;
    ans : character;

  begin
    put_line("the polynomial system :");
    DoblDobl_CSeries_Poly_Systems_io.put(p);
    Compute(s.pwt,s.mxe,xcff);
    EvalDiff(s,xcff);
    put_line("The evaluation values :"); put_line(px);
    put_line("The coefficient vectors of the evaluation :"); put_line(s.yv);
    err := Difference(px,s.yv);
    put("The error : "); put(err,3); new_line;
    DoblDobl_Vector_Splitters.Complex_Parts(xcff,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.Compute
      (cs.rhpwt,cs.ihpwt,cs.rlpwt,cs.ilpwt,cs.mxe,rhx,ihx,rlx,ilx);
    DoblDobl_Coefficient_Convolutions.EvalDiff
      (cs,rhx.all,ihx.all,rlx.all,ilx.all);
    err := Difference(px,cs.yv);
    put("Recomputed error : "); put(err,3); new_line;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
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
    put("Recomputed error :"); put(err,3); new_line;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    DoblDobl_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,xpt);
    DoblDobl_Speelpenning_Convolutions.EvalDiff(s,xpt);
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

  procedure DoblDobl_Input_Test ( deg : in integer32 ) is

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
      DoblDobl_Input_Test(deg);
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
        DoblDobl_System_Test(dim,deg,nbr,pwr);
      else
        put("All coefficients equal to one ? (y/n) ");
        Ask_Yes_or_No(answer); cffone := (answer = 'y');
        DoblDobl_Test(dim,deg,nbr,pwr,expone,cffone);
      end if;
    end if;
  end Main;

end Test_DoblDobl_Coeff_Convolutions;
