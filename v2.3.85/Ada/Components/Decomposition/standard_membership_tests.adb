with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Planes_and_Polynomials;             use Planes_and_Polynomials;
with Standard_Linear_Projections;        use Standard_Linear_Projections;
with Standard_Central_Projections;       use Standard_Central_Projections;

package body Standard_Membership_Tests is

  function In_Subspace
               ( file : in file_type;
                 p : Standard_Complex_Poly_Systems.Poly_Sys;
                 ind : integer32; x : Standard_Complex_Vectors.Vector;
                 tol : double_float ) return boolean is

    eva : constant Standard_Complex_Vectors.Vector(p'range) := Eval(p,x);
    val : constant double_float := Max_Norm(eva);
    res : constant boolean := (val <= tol);

  begin
    put(file,"Residual in span of component "); put(file,ind,1);
    put(file," : "); put(file,val,3);
    if res
     then put_line(file,"  success");
     else put_line(file,"  failure");
    end if;
    return res;
  end In_Subspace;

  function On_Component
             ( file : file_type;
               p : Standard_Complex_Poly_Systems.Poly_Sys; ind : integer32;
               projp : Standard_Complex_Vectors.Vector;
               tol : double_float ) return boolean is

  -- DESCRIPTION :
  --   Auxiliary routine for the on component test for solutions.
  --   The projp vector is the resulting vector after projection.

    res : boolean := false;
    eva : Complex_Number;
    tmp : double_float;

  begin
    for i in 1..(ind-1) loop
      if p(i) /= Null_Poly then
        eva := Eval(p(i),projp);
        tmp := AbsVal(eva);
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
             ( file : file_type; intpols : Poly_Sys; ind : integer32;
               sols : Solution_List; hyp : Standard_Complex_VecVecs.VecVec;
               level : integer32; tol : double_float ) return boolean is

    use Standard_Complex_Vectors;

  begin
    if ind = 1 then
      return false;
    else
      declare
        point : constant Vector := Retrieve(sols,natural32(ind)).v;
        projp : constant Vector := Evaluate(hyp,point,level+1);
      begin
        return On_Component(file,intpols,ind,projp,tol);
      end;
    end if;
  end On_Component;

  function On_Component
             ( file : file_type; intpols : Poly_Sys; ind : integer32;
               pivots : Standard_Integer_VecVecs.VecVec;
               sols : Solution_List; hyp : Standard_Complex_VecVecs.VecVec;
               level : integer32; tol : double_float ) return boolean is

    use Standard_Integer_Vectors;

  begin
    if ind = 1 then
      return false;
    else
      declare
        point : constant Standard_Complex_Vectors.Vector
              := Retrieve(sols,natural32(ind)).v;
        projp : Standard_Complex_Vectors.Vector(1..level+1);
        eva : Complex_Number;
        tmp : double_float;
      begin
        for i in 1..(ind-1) loop
          if intpols(i) /= Null_Poly then
            if pivots(i) = null then
              projp := Evaluate(hyp,point,level+1);
              eva := Eval(intpols(i),projp);
            else
              declare
                dim : constant integer32 := pivots(i)'last+level; 
                pivpt : constant Standard_Complex_Vectors.Vector 
                      := Remove_Variables(point,level,dim,pivots(i).all);
                prp : constant Standard_Complex_Vectors.Vector
                    := Evaluate(hyp,pivpt,level+1);
              begin
                eva := Eval(intpols(i),prp);
              end;
            end if;
            tmp := AbsVal(eva);
            put(file,"Residual in membership test at component ");
            put(file,i,1); put(file," : ");
            put(file,tmp,3,2,3);
            if tmp < tol then
              put_line(file,"  success");
              return true;
            else 
              put_line(file,"  failure");
            end if;
          end if;
        end loop;
        return false;
      end;
    end if;
  end On_Component;

  function Evaluate ( p : Standard_Complex_Polynomials.Poly;
                      x : Standard_Complex_Vectors.Vector;
                      level : integer32;
                      hyp : Standard_Complex_VecVecs.VecVec;
                      pivots : Standard_Integer_Vectors.Link_to_Vector;
                      base : Standard_Complex_VecVecs.Link_to_VecVec;
                      basecard : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns p(r(x)), with r the suitable projection.

    use Standard_Integer_Vectors;
    use Standard_Complex_Vectors;
    use Standard_Complex_VecVecs;

    res : Complex_Number;
    projx : Standard_Complex_Vectors.Vector(1..level+1);

  begin
    if pivots = null then
      projx := Evaluate(hyp,x,level+1);
      res := Eval(p,projx);
    else
      declare
        dim : constant integer32 := pivots'last+level;
        restx : constant Standard_Complex_Vectors.Vector(1..dim)
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
      end;
    end if;
    return res;
  end Evaluate;

  function On_Component
               ( file : file_type;
                 p : Standard_Complex_Poly_Systems.Poly_Sys; ind : integer32;
                 subspaces : Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                 pivots : Standard_Integer_VecVecs.VecVec;
                 basepts : Standard_Complex_VecVecs.Array_of_VecVecs;
                 basecard : Standard_Natural_Vectors.Vector;
                 sols : Standard_Complex_Solutions.Solution_List;
                 hyp : Standard_Complex_VecVecs.VecVec;
                 level : integer32; tol : double_float )
               return boolean is

    res,sub : boolean;
    use Standard_Integer_Vectors;

  begin
    if ind = 1 then
      res := false;
    else
      declare
        point : constant Standard_Complex_Vectors.Vector
              := Retrieve(sols,natural32(ind)).v;
        eva : Complex_Number;
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
              tmp := AbsVal(eva);
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

end Standard_Membership_Tests;
