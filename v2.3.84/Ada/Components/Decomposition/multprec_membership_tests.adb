with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Standard_Integer_Vectors;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Multprec_Complex_Vector_Tools;     use Multprec_Complex_Vector_Tools;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Functions;   use Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_SysFun;      use Multprec_Complex_Poly_SysFun;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Multprec_Linear_Projections;       use Multprec_Linear_Projections;
with Multprec_Central_Projections;      use Multprec_Central_Projections;

package body Multprec_Membership_Tests is

  function In_Subspace
               ( file : in file_type;
                 p : Multprec_Complex_Poly_Systems.Poly_Sys;
                 ind : integer32; x : Standard_Complex_Vectors.Vector;
                 tol : double_float ) return boolean is

    mpx : Multprec_Complex_Vectors.Vector(x'range) := Create(x);
    res : constant boolean := In_Subspace(file,p,ind,mpx,tol);

  begin
    Multprec_Complex_Vectors.Clear(mpx);
    return res;
  end In_Subspace;

  function In_Subspace
               ( file : in file_type;
                 p : Multprec_Complex_Poly_Systems.Poly_Sys;
                 ind : integer32; x : Multprec_Complex_Vectors.Vector;
                 tol : double_float ) return boolean is

    eva : Multprec_Complex_Vectors.Vector(p'range) := Eval(p,x);
    val : Floating_Number := Max_Norm(eva);
    fltval : constant double_float := Trunc(val);
    res : constant boolean := (fltval <= tol);

  begin
    put(file,"Residual in span of component "); put(file,ind,1);
    put(file," : "); put(file,fltval,3);
    if res
     then put_line(file,"  success");
     else put_line(file,"  failure");
    end if;
    Multprec_Complex_Vectors.Clear(eva);
    Clear(val);
    return res;
  end In_Subspace;

  function On_Component
             ( file : file_type;
               p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
               projp : Multprec_Complex_Vectors.Vector;
               tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Auxiliary routine for the on component test for solutions.
  --   The projp vector is the resulting vector after projection.

    res : boolean := false;
    eva : Complex_Number;
    fltacc : Floating_Number;
    tmp : double_float;

  begin
    for i in 1..(ind-1) loop
      if p(i) /= Null_Poly then
        eva := Eval(p(i),projp);
        fltacc := AbsVal(eva);
        tmp := Trunc(fltacc);
        Clear(fltacc);
        put(file,"Residual in membership test at component ");
        put(file,i,1); put(file," : ");
        put(file,tmp,3,2,3);
        if tmp < tol then
          res := true;
          put_line(file,"  success");
        else
          put_line(file,"  failure");
        end if;
      end if;
      exit when res;
    end loop;
    return res;
  end On_Component;

  function On_Component
             ( file : file_type;
               p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
               sols : Standard_Complex_Solutions.Solution_List;
               hyp : Multprec_Complex_VecVecs.VecVec;
               level : integer32; size : natural32; tol : double_float )
             return boolean is

    res : boolean;
    use Standard_Complex_Solutions;

  begin
    if ind = 1 then
      res := false;
    else
      declare
        point : constant Standard_Complex_Vectors.Vector
              := Retrieve(sols,natural32(ind)).v;
        projp : Multprec_Complex_Vectors.Vector(1..level+1)
              := Evaluate(hyp,point,level+1);
      begin
        res := On_Component(file,p,ind,projp,tol);
        Multprec_Complex_Vectors.Clear(projp);
      end;
    end if;
    return res;
  end On_Component;

  function On_Component
             ( file : file_type;
                p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
                sols : Multprec_Complex_Solutions.Solution_List;
                hyp : Multprec_Complex_VecVecs.VecVec;
                level : integer32; size : natural32; tol : double_float )
              return boolean is

    res : boolean;
    use Multprec_Complex_Solutions;

  begin
    if ind = 1 then
      res := false;
    else
      declare
        point : constant Multprec_Complex_Vectors.Vector
              := Retrieve(sols,natural32(ind)).v;
        extpt : Multprec_Complex_Vectors.Vector(point'range);
        projp : Multprec_Complex_Vectors.Vector(1..level+1);
      begin
        Multprec_Complex_Vectors.Copy(point,extpt);
        Set_Size(extpt,2*size);
        projp := Evaluate(hyp,extpt,level+1);
        res := On_Component(file,p,ind,projp,tol);
        Multprec_Complex_Vectors.Clear(projp);
        Multprec_Complex_Vectors.Clear(extpt);
      end;
    end if;
    return res;
  end On_Component;

  function Evaluate ( p : Multprec_Complex_Polynomials.Poly;
                      x : Standard_Complex_Vectors.Vector;
                      level : integer32;
                      hyp : Multprec_Complex_VecVecs.VecVec;
                      pivots : Standard_Integer_Vectors.Link_to_Vector;
                      base : Multprec_Complex_VecVecs.Link_to_VecVec;
                      basecard : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns p(r(x)), with r the suitable projection.

    use Standard_Integer_Vectors;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_VecVecs;

    res : Complex_Number;
    projx : Multprec_Complex_Vectors.Vector(1..level+1);

  begin
    if pivots = null then
      projx := Evaluate(hyp,x,level+1);
      res := Eval(p,projx);
    else
      declare
        dim : constant integer32 := pivots'last+level;
        restx : constant Standard_Complex_Vectors.Vector
              := Remove_Variables(x,level,dim,pivots.all);
        mprestx : Multprec_Complex_Vectors.Vector(restx'range);
      begin
        if base = null then
          projx := Evaluate(hyp,restx,level+1);
          res := Eval(p,projx);
        else
          mprestx := Create(restx);
          projx := Intersect(hyp(1..basecard),base(1..basecard),
                             mprestx,level+1);
          res := Eval(p,projx);
          Multprec_Complex_Vectors.Clear(mprestx);
        end if;
      end;
     end if;
     Multprec_Complex_Vectors.Clear(projx);
     return res;
  end Evaluate;

  function Evaluate ( p : Multprec_Complex_Polynomials.Poly;
                      x : Multprec_Complex_Vectors.Vector;
                      level : integer32;
                      hyp : Multprec_Complex_VecVecs.VecVec;
                      pivots : Standard_Integer_Vectors.Link_to_Vector;
                      base : Multprec_Complex_VecVecs.Link_to_VecVec;
                      basecard : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns p(r(x)), with r the suitable projection.

    use Standard_Integer_Vectors;
    use Multprec_Complex_Vectors;
    use Multprec_Complex_VecVecs;

    res : Complex_Number;
    projx : Multprec_Complex_Vectors.Vector(1..level+1);

  begin
    if pivots = null then
      projx := Evaluate(hyp,x,level+1);
      res := Eval(p,projx);
    else
      declare
        dim : constant integer32 := pivots'last+level;
        restx : Multprec_Complex_Vectors.Vector(1..dim)
              := Remove_Variables(x,level,dim,pivots.all);
      begin
        if base = null then
          projx := Evaluate(hyp,restx,level+1);
          res := Eval(p,projx);
        else
          projx := Intersect(hyp(1..basecard),base(1..basecard),
                             restx,level+1);
          res := Eval(p,projx);
        end if;
        Multprec_Complex_Vectors.Clear(restx);
      end;
    end if;
    Multprec_Complex_Vectors.Clear(projx);
    return res;
  end Evaluate;

  function On_Component
               ( file : file_type;
                 p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
                 subspaces : Multprec_Complex_Poly_Systems.Array_of_Poly_Sys;
                 pivots : Standard_Integer_VecVecs.VecVec;
                 basepts : Multprec_Complex_VecVecs.Array_of_VecVecs;
                 basecard : Standard_Natural_Vectors.Vector;
                 sols : Standard_Complex_Solutions.Solution_List;
                 hyp : Multprec_Complex_VecVecs.VecVec;
                 level : integer32; size : natural32; tol : double_float )
               return boolean is

    res,sub : boolean;
    use Standard_Integer_Vectors;
    use Standard_Complex_Solutions;

  begin
    if ind = 1 then
      res := false;
    else
      declare
        point : constant Standard_Complex_Vectors.Vector
              := Retrieve(sols,natural32(ind)).v;
        eva : Complex_Number;
        fltacc : Floating_Number;
        tmp : double_float;
      begin
        res := false;
        for i in 1..(ind-1) loop
          if p(i) /= Null_Poly then
            if pivots(i) = null
             then sub := true;
             else sub := In_Subspace(file,subspaces(i).all,i,point,tol);
            end if;
            if sub then
              eva := Evaluate(p(i),point,level,hyp,pivots(i),
                              basepts(i),integer32(basecard(i)));
              fltacc := AbsVal(eva);
              tmp := Trunc(fltacc);
              Clear(fltacc);
              put(file,"Residual in membership test at component ");
              put(file,i,1); put(file," : ");
              put(file,tmp,3,2,3);
              if tmp < tol then
                res := true;
                put_line(file,"  success");
              else
                put_line(file,"  failure");
              end if;
            end if;
          end if;
          exit when res;
        end loop;
      end;
    end if;
    return res;
  end On_Component;

  function On_Component
             ( file : file_type;
               p : Multprec_Complex_Poly_Systems.Poly_Sys; ind : integer32;
               subspaces : Multprec_Complex_Poly_Systems.Array_of_Poly_Sys;
               pivots : Standard_Integer_VecVecs.VecVec;
               basepts : Multprec_Complex_VecVecs.Array_of_VecVecs;
               basecard : Standard_Natural_Vectors.Vector;
               sols : Multprec_Complex_Solutions.Solution_List;
               hyp : Multprec_Complex_VecVecs.VecVec;
               level : in integer32; size : natural32; tol : double_float )
             return boolean is

    res,sub : boolean;
    use Standard_Integer_Vectors;
    use Multprec_Complex_Solutions;

  begin
    if ind = 1 then
      res := false;
    else
      declare
        point : constant Multprec_Complex_Vectors.Vector
              := Retrieve(sols,natural32(ind)).v;
        extpt : Multprec_Complex_Vectors.Vector(point'range);
        eva : Complex_Number;
        fltacc : Floating_Number;
        tmp : double_float;
      begin
        res := false;
        Multprec_Complex_Vectors.Copy(point,extpt);
        Set_Size(extpt,2*size);
        for i in 1..(ind-1) loop
          if p(i) /= Null_Poly then
            if pivots(i) = null
             then sub := true;
             else sub := In_Subspace(file,subspaces(i).all,i,extpt,tol);
            end if;
            if sub then
              eva := Evaluate(p(i),extpt,level,hyp,pivots(i),
                              basepts(i),integer32(basecard(i)));
              fltacc := AbsVal(eva);
              tmp := Trunc(fltacc);
              Clear(fltacc);
              put(file,"Residual in membership test at component ");
              put(file,i,1); put(file," : ");
              put(file,tmp,3,2,3);
              if tmp < tol then
                res := true;
                put_line(file,"  success");
              else
                put_line(file,"  failure");
              end if;
            end if;
          end if;
          exit when res;
        end loop;
      end;
    end if;
    return res;
  end On_Component;

end Multprec_Membership_Tests;
