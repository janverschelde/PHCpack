with Standard_Floating_Numbers;          use Standard_Floating_Numbers;

package body DEMiCs_Output_Convertors is

  function Apply_Lifting
             ( sup : Lists_of_Integer_Vectors.List;
               lif : Standard_Floating_Vectors.Vector )
             return Lists_of_Floating_Vectors.List is

    res,res_last : Lists_of_Floating_Vectors.List;
    tmp : Lists_of_Integer_Vectors.List := sup;
    ilv : Standard_Integer_Vectors.Link_to_Vector;
    idx : integer32 := 0;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      idx := idx + 1;
      ilv := Lists_of_Integer_Vectors.Head_Of(tmp);
      declare
        flv : Standard_Floating_Vectors.Vector(ilv'first..ilv'last+1);
      begin
        for k in ilv'range loop
          flv(k) := double_float(ilv(k));
        end loop;
        flv(flv'last) := lif(idx);
        Lists_of_Floating_Vectors.Append(res,res_last,flv);
      end;
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Apply_Lifting;

  function Apply_Lifting
             ( sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lif : Standard_Floating_VecVecs.VecVec )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range);

  begin
    for k in sup'range loop
      res(k) := Apply_Lifting(sup(k),lif(k).all);
    end loop;
    return res;
  end Apply_Lifting;

  function Make_Mixed_Cell
             ( dim : integer32;
               lbl : Standard_Integer_Vectors.Vector;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
             return Mixed_Cell is

    res : Mixed_Cell;
    idxfirst,idxsecond : integer32;
    tmp,last : Lists_of_Floating_Vectors.List;
    idx,done : integer32;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(1..dim+1 => 0.0);
    res.pts := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..dim);
    res.sub := null;
    for k in 1..dim loop
      last := res.pts(k);
      idxfirst := lbl(2*(k-1)+1);
      idxsecond := lbl(2*(k-1)+2);
      idx := 0;
      done := 2;
      tmp := lifsup(k);
      while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
        idx := idx + 1;
        if idx = idxfirst then
          lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
          Lists_of_Floating_Vectors.Append(res.pts(k),last,lpt.all);
          done := done - 1;
        elsif idx = idxsecond then
          lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
          Lists_of_Floating_Vectors.Append(res.pts(k),last,lpt.all);
          done := done - 1;
        end if;
        exit when (done = 0);
        tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Make_Mixed_Cell;

  function Make_Mixed_Cells
             ( dim : integer32;
               labels : Lists_of_Integer_Vectors.List;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
             return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Lists_of_Integer_Vectors.List := labels;
    lbl : Standard_Integer_Vectors.Link_to_Vector;
    mic : Mixed_Cell;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      lbl := Lists_of_Integer_Vectors.Head_Of(tmp);
      mic := Make_Mixed_Cell(dim,lbl.all,lifsup);
      Append(res,res_last,mic);
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Make_Mixed_Cells;

end DEMiCs_Output_Convertors;
