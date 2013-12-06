package body Standard_Linear_Projections is

  function Evaluate ( hyp,point : Vector ) return Complex_Number is
 
    res : Complex_Number := hyp(0);

  begin
    for i in point'range loop
      res := res + hyp(i)*point(i);
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : VecVec; point : Vector; n : integer32 )
                    return Vector is

    res : Vector(hyps'first..n);

  begin
    for i in res'range loop
      res(i) := Evaluate(hyps(i).all,point);
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps,pts : VecVec; n : integer32 ) return VecVec is

    res : VecVec(pts'range);

  begin
    for i in pts'range loop
      if pts(i) /= null
       then res(i) := new Vector'(Evaluate(hyps,pts(i).all,n));
      end if;
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : VecVec; arvv : Array_of_VecVecs; n : integer32 )
                    return Array_of_VecVecs is

    res : Array_of_VecVecs(arvv'range);

  begin
    for i in arvv'range loop
      res(i) := new VecVec'(Evaluate(hyps,arvv(i).all,n));
     -- res(i) := new VecVec(arvv(i)'range);
     -- for j in arvv(i)'range loop
     --   res(i)(j) := new Vector'(Evaluate(hyps,arvv(i)(j).all,n));
     -- end loop;
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : VecVec; sols : Solution_List; n : integer32 )
                    return VecVec is
  
    res : VecVec(1..integer32(Length_Of(sols)));
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := new Vector'(Evaluate(hyps,ls.v,n));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Evaluate;

end Standard_Linear_Projections;
