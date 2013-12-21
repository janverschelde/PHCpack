with unchecked_deallocation;
--with Standard_Affine_Planes;           use Standard_Affine_Planes;
with Standard_Plane_Operations;        use Standard_Plane_Operations;
with Standard_Affine_Solutions;        use Standard_Affine_Solutions;

package body Multihomogeneous_Solutions is

-- AUXILIARY OPERATIONS :

  function Remove_Empty_Lists ( sols : Array_of_Solution_Lists )
                              return Array_of_Solution_Lists is

    res : Array_of_Solution_Lists(sols'range);
    cnt : natural := res'first-1;

  begin
    for i in sols'range loop
      if not Is_Null(sols(i))
       then cnt := cnt + 1;
            res(cnt) := sols(i);
      end if;
    end loop;
    return res(res'first..cnt);
  end Remove_Empty_Lists;

  function Remove_Empty_VecVecs ( w : Array_of_VecVecs )
                                return Array_of_VecVecs is

    res : Array_of_VecVecs(w'range);
    cnt : natural := res'first-1;

  begin
    for i in w'range loop
      if w(i) /= null
       then cnt := cnt + 1;
            res(cnt) := w(i);
      end if;
    end loop;
    return res(res'first..cnt);
  end Remove_Empty_VecVecs;

-- CREATORS :

  function Create ( n : natural; b : Vector; w : Array_of_VecVecs;
                    s : Array_of_Solution_Lists )
                  return Multihomogeneous_Solution is

    new_w : constant Array_of_VecVecs := Remove_Empty_VecVecs(w);
    aux_w : Array_of_VecVecs(new_w'range) := Truncate(new_w,n);
    ort_w : Array_of_VecVecs(aux_w'range) := Orthogonalize(aux_w);
    new_s : constant Array_of_Solution_Lists := Remove_Empty_Lists(s);
    rwt_s : Array_of_Solution_Lists(new_s'range)
          := Rewrite(new_s,n,b,aux_w,ort_w);
    res : Multihomogeneous_Solution(n,new_w'length);

  begin
    res.b := b;
    res.k := ort_w(ort_w'first)'length;
    res.w := ort_w;
    res.s := rwt_s;
    return res;
  end Create;

  function Create ( n : natural; b : Vector; w : Array_of_VecVecs;
                    s : Array_of_Solution_Lists )
                  return Link_to_Multihomogeneous_Solution is
  begin
    return new Multihomogeneous_Solution'(Create(n,b,w,s));
  end Create;

-- DESTRUCTORS :

  procedure Clear ( mhs : in out Multihomogeneous_Solution ) is
  begin
    for i in mhs.w'range loop
      Deep_Clear(mhs.w(i));
      Clear(mhs.s(i));
    end loop;
  end Clear;

  procedure Clear ( mhs : in out Link_to_Multihomogeneous_Solution ) is

    procedure free is 
      new unchecked_deallocation(Multihomogeneous_Solution,
                                 Link_to_Multihomogeneous_Solution);
  begin
    if mhs /= null
     then Clear(mhs.all);
          free(mhs);
    end if;
  end Clear;

end Multihomogeneous_Solutions;
