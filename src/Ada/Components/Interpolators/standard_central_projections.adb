package body Standard_Central_Projections is

-- PROJECTION OF ONE POINT :

  function Intersect ( hyp,base,point : Vector; evabas : Complex_Number;
                       dim : integer32 ) return Vector is

    res : Vector(1..dim);
    evapnt,t : Complex_Number;

  begin
    evapnt := hyp(point'range)*point;
    t := (hyp(0) + evabas)/(evapnt - evabas);
    for i in 1..dim loop
      res(i) := base(i) - point(i);
    end loop;
    Mul(res,t);
    for i in 1..dim loop
      res(i) := res(i) + base(i);
    end loop;
    return res;
  end Intersect;

  function Intersect ( hyp,base,point : Vector; dim : integer32 )
                     return Vector is

    evabas : Complex_Number := hyp(base'range)*base;

  begin
    return Intersect(hyp,base,point,evabas,dim);
  end Intersect;

  procedure Intersect_Base_Points ( hyp : in VecVec; base : in out VecVec ) is

    dim : constant integer32 := base(base'first)'last;
    projpt : Vector(1..dim);

  begin
    for i in base'first..base'last-1 loop    -- use base(i) to
      for j in i+1..base'last loop           -- project base(j) onto hyp(i)
        projpt := Intersect(hyp(i).all,base(i).all,base(j).all,dim);
        Clear(base(j));
        base(j) := new Vector'(projpt);
      end loop;
    end loop;
  end Intersect_Base_Points;

  function Intersect ( hyp,base : VecVec; point : Vector; dim : integer32 )
                     return Vector is

    res : Vector(1..dim);
    wrk : Vector(point'range) := point;

  begin
    for i in base'range loop
      exit when ((hyp(i) = null) or (base(i) = null));
      wrk := Intersect(hyp(i).all,base(i).all,wrk,wrk'last);
    end loop;
    res := wrk(1..dim);
    return res;
  end Intersect;

-- PROJECTION OF A SEQUENCE OF POINTS :

  function Intersect ( hyp,base : Vector; evabas : Complex_Number;
                       points : VecVec; dim : integer32 ) return VecVec is

    res : VecVec(points'range);

  begin
    for i in res'range loop
      res(i) := new Vector'(Intersect(hyp,base,points(i).all,evabas,dim));
    end loop;
    return res;
  end Intersect;

  function Intersect ( hyp,base : Vector; points : VecVec; dim : integer32 )
                     return VecVec is

    evabas : Complex_Number := hyp(base'range)*base;

  begin
    return Intersect(hyp,base,evabas,points,dim);
  end Intersect;

end Standard_Central_Projections;
