with Standard_Point_Coordinates;        use Standard_Point_Coordinates;

package body Standard_Affine_Solutions is

  function Rewrite ( sol : Solution; n : integer32; b : Vector; v,w : VecVec )
                   return Solution is

    res : Solution(sol.n) := sol;
    x : constant Vector(1..n) := Affine_Expand(sol.v,b,v);

  begin
    res.v := Project(x,b,w);
    return res;
  end Rewrite;

  function Rewrite ( sols : Solution_List; n : integer32;
                     b : Vector; v,w : VecVec ) return Solution_List is

    res,res_last : Solution_List;
    ls : Link_to_Solution;
    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Rewrite(ls.all,n,b,v,w));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Rewrite;

  function Rewrite ( sols : Array_of_Solution_Lists; n : integer32;
                     b : Vector; v,w : Array_of_VecVecs ) 
                   return Array_of_Solution_Lists is

    res : Array_of_Solution_Lists(sols'range);

  begin
    for i in sols'range loop
      res(i) := Rewrite(sols(i),n,b,v(i).all,w(i).all);
    end loop;
    return res;
  end Rewrite;

end Standard_Affine_Solutions;
