with Multprec_Complex_Vector_Tools;    use Multprec_Complex_Vector_Tools;

package body Multprec_Linear_Projections is

  function Evaluate ( hyp : Multprec_Complex_Vectors.Vector;
                      point : Standard_Complex_Vectors.Vector )
                    return Complex_Number is
 
    mpt : Multprec_Complex_Vectors.Vector(point'range) := Create(point);
    res : constant Complex_Number := Evaluate(hyp,mpt);
    
  begin
    Multprec_Complex_Vectors.Clear(mpt);
    return res;
  end Evaluate;

  function Evaluate ( hyp,point : Multprec_Complex_Vectors.Vector )
                    return Complex_Number is
 
    res,acc : Complex_Number;

  begin
    Copy(hyp(0),res);
    for i in point'range loop
      acc := hyp(i)*point(i);
      Add(res,acc);
      Clear(acc);
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      point : Standard_Complex_Vectors.Vector; n : integer32 )
                    return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(hyps'first..n);

  begin
    for i in res'range loop
      res(i) := Evaluate(hyps(i).all,point);
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      point : Multprec_Complex_Vectors.Vector; n : integer32 )
                    return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(hyps'first..n);

  begin
    for i in res'range loop
      res(i) := Evaluate(hyps(i).all,point);
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      points : Standard_Complex_VecVecs.VecVec;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(points'range);
 
  begin
    for i in points'range loop
      res(i) := new Multprec_Complex_Vectors.Vector'
                      (Evaluate(hyps,points(i).all,n));
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps,points : Multprec_Complex_VecVecs.VecVec;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(points'range);

  begin
    for i in points'range loop
      res(i) := new Multprec_Complex_Vectors.Vector'
                      (Evaluate(hyps,points(i).all,n));
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec; 
                      arvv : Standard_Complex_VecVecs.Array_of_VecVecs;
                      n : integer32 )
                    return Multprec_Complex_VecVecs.Array_of_VecVecs is

    res : Multprec_Complex_VecVecs.Array_of_VecVecs(arvv'range);

  begin
    for i in arvv'range loop
      res(i) := new Multprec_Complex_VecVecs.VecVec(arvv(i)'range);
      for j in arvv(i)'range loop
        res(i)(j) := new Multprec_Complex_Vectors.Vector'
                           (Evaluate(hyps,arvv(i)(j).all,n));
      end loop;
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec; 
                      arvv : Multprec_Complex_VecVecs.Array_of_VecVecs;
                      n : integer32 )
                    return Multprec_Complex_VecVecs.Array_of_VecVecs is

    res : Multprec_Complex_VecVecs.Array_of_VecVecs(arvv'range);

  begin
    for i in arvv'range loop
      res(i) := new Multprec_Complex_VecVecs.VecVec(arvv(i)'range);
      for j in arvv(i)'range loop
        res(i)(j) := new Multprec_Complex_Vectors.Vector'
                           (Evaluate(hyps,arvv(i)(j).all,n));
      end loop;
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      sols : Standard_Complex_Solutions.Solution_List;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec is
  
    use Standard_Complex_Solutions;
    res : Multprec_Complex_VecVecs.VecVec(1..integer32(Length_Of(sols)));
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := new Multprec_Complex_Vectors.Vector'(Evaluate(hyps,ls.v,n));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Evaluate;

  function Evaluate ( hyps : Multprec_Complex_VecVecs.VecVec;
                      sols : Multprec_Complex_Solutions.Solution_List;
                      n : integer32 ) return Multprec_Complex_VecVecs.VecVec is
  
    use Multprec_Complex_Solutions;
    res : Multprec_Complex_VecVecs.VecVec(1..integer32(Length_Of(sols)));
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := new Multprec_Complex_Vectors.Vector'(Evaluate(hyps,ls.v,n));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Evaluate;

end Multprec_Linear_Projections;
