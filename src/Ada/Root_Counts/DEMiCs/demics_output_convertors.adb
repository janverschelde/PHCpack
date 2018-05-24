with text_io;                            use text_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;

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

  function Apply_Lifting
              ( mix : Standard_Integer_Vectors.Vector;
                sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : Standard_Floating_VecVecs.VecVec )
              return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    idx : integer32 := sup'first;

  begin
    for k in mix'range loop
      res(k) := Apply_Lifting(sup(idx),lif(k).all);
      idx := idx + mix(k);
    end loop;
    return res;
  end Apply_Lifting;

  function Minimum ( lifpts : Lists_of_Floating_Vectors.List;
                     normal : Standard_Floating_Vectors.Vector )
                   return double_float is

    use Standard_Floating_Vectors;

    lpt : Standard_Floating_Vectors.Link_to_Vector
        := Lists_of_Floating_Vectors.Head_Of(lifpts);
    res : double_float := lpt.all*normal;
    ipr : double_float;
    tmp : Lists_of_Floating_Vectors.List
        := Lists_of_Floating_Vectors.Tail_Of(lifpts);

  begin
    while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
      lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
      ipr := lpt.all*normal;
      if ipr < res
       then res := ipr;
      end if;
      tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Minimum;

  function Arguments_of_Minima
             ( lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               normal : Standard_Floating_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    use Standard_Floating_Vectors;

    dim : constant integer32 := lifsup'last;
    res : Standard_Integer_Vectors.Vector(1..2*dim);
    min,ipr : double_float;
    tmp : Lists_of_Floating_Vectors.List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
    idxlpt,idxres : integer32;
 
  begin
    idxres := 0;
    for k in lifsup'range loop
      min := Minimum(lifsup(k),normal);
      idxlpt := 0;
      tmp := lifsup(k);
      while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
        idxlpt := idxlpt + 1;
        lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
        ipr := lpt.all*normal;
        if abs(ipr - min) < 1.0E-8 then
          idxres := idxres + 1;
          res(idxres) := idxlpt;
        end if;
        tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Arguments_of_Minima;

  function Arguments_of_Minima
             ( mix : Standard_Integer_Vectors.Vector;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               normal : Standard_Floating_Vectors.Vector )
             return Standard_Integer_Vectors.Vector is

    use Standard_Floating_Vectors;

    dim : constant integer32
        := Standard_Integer_Vectors.Sum(mix) + mix'last; 
    res : Standard_Integer_Vectors.Vector(1..dim);
    min,ipr : double_float;
    tmp : Lists_of_Floating_Vectors.List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
    idxlpt,idxres : integer32;
 
  begin
    idxres := 0;
    for k in lifsup'range loop
      min := Minimum(lifsup(k),normal);
      idxlpt := 0;
      tmp := lifsup(k);
      while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
        idxlpt := idxlpt + 1;
        lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
        ipr := lpt.all*normal;
        if abs(ipr - min) < 1.0E-8 then
          idxres := idxres + 1;
          res(idxres) := idxlpt;
        end if;
        tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
      end loop;
    end loop;
    return res;
  end Arguments_of_Minima;

  function Sort_Labels
              ( labels : Standard_Integer_Vectors.Vector )
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(labels'range) := labels;
    dim : constant integer32 := labels'last/2;
    idx,tmp : integer32;

  begin
    for k in 1..dim loop
      idx := 2*(k-1) + 1;
      if res(idx) > res(idx+1) then
        tmp := res(idx);
        res(idx) := res(idx+1);
        res(idx+1) := tmp;
      end if;
    end loop;
    return res;
  end Sort_Labels;

  procedure Sort ( x : in out Standard_Integer_Vectors.Vector;
                   xfirst,xlast : in integer32 ) is

    min,idx : integer32;

  begin
    for i in xfirst..xlast-1 loop 
      idx := i; min := x(idx);
      for j in i+1..xlast loop         -- invariant: min = x(idx)
        if x(j) < min 
         then idx := j; min := x(idx); -- new minimum found
        end if;
      end loop;
      if idx /= i then  -- swap min = x(idx) to x(i)
        x(idx) := x(i); -- overwrite x(idx) first, backup is min 
        x(i) := min;    -- assign min to x(i)
      end if;
    end loop;
  end Sort;

  function Sort_Labels
              ( mix,labels : Standard_Integer_Vectors.Vector )
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(labels'range) := labels;
    offset : integer32 := res'first;

  begin
    for k in mix'range loop
      Sort(res,offset,offset+mix(k));
      offset := offset + mix(k) + 1;
    end loop;
    return res;
  end Sort_Labels;

  procedure Test_Inner_Normal
              ( lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                normal : in Standard_Floating_Vectors.Vector;
                labels : in Standard_Integer_Vectors.Vector;
                fail : out boolean ) is

    argmin : constant Standard_Integer_Vectors.Vector
           := Arguments_of_Minima(lifsup,normal);
    sorted : constant Standard_Integer_Vectors.Vector
           := Sort_Labels(labels);

  begin
    put("The labels of demics : "); put(sorted); new_line;
    put("The computed labels  : "); put(argmin);
    fail := false;
    for k in sorted'range loop
      if sorted(k) /= argmin(k)
       then fail := true;
      end if;
      exit when fail;
    end loop;
    if fail
     then put_line("  wrong!?");
     else put_line("  okay.");
    end if;
  end Test_Inner_Normal;

  procedure Test_Inner_Normal
              ( mix : in Standard_Integer_Vectors.Vector;
                lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                normal : in Standard_Floating_Vectors.Vector;
                labels : in Standard_Integer_Vectors.Vector;
                fail : out boolean ) is

    argmin : constant Standard_Integer_Vectors.Vector
           := Arguments_of_Minima(mix,lifsup,normal);
    sorted : constant Standard_Integer_Vectors.Vector
           := Sort_Labels(mix,labels);

  begin
    put("The labels of demics : "); put(labels); new_line;
    put("The sorted labels    : "); put(sorted); new_line;
    put("The computed labels  : "); put(argmin);
    fail := false;
    for k in sorted'range loop
      if sorted(k) /= argmin(k)
       then fail := true;
      end if;
      exit when fail;
    end loop;
    if fail
     then put_line("  wrong!?");
     else put_line("  okay.");
    end if;
  end Test_Inner_Normal;

  function Make_Mixed_Cell
              ( dim : integer32;
                lbl : Standard_Integer_Vectors.Vector;
                lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                verbose : boolean := true )
              return Mixed_Cell is

    res : Mixed_Cell;
    idxfirst,idxsecond : integer32;
    tmp,last : Lists_of_Floating_Vectors.List;
    idx,done : integer32;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
    first,second : Standard_Floating_Vectors.Vector(1..dim+1);
    mat : Standard_Floating_Matrices.Matrix(1..dim,1..dim);
    rhs : Standard_Floating_Vectors.Vector(1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    fail : boolean;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(1..dim+1 => 0.0);
    res.pts := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..dim);
    res.sub := null;
    for row in 1..dim loop
      last := res.pts(row);
      idxfirst := lbl(2*(row-1)+1);
      idxsecond := lbl(2*(row-1)+2);
      idx := 0;
      done := 2;
      tmp := lifsup(row);
      while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
        idx := idx + 1;
        if idx = idxfirst then
          lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
          Lists_of_Floating_Vectors.Append(res.pts(row),last,lpt.all);
          for i in lpt'range loop
            first(i) := lpt(i);
          end loop;
          done := done - 1;
        elsif idx = idxsecond then
          lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
          Lists_of_Floating_Vectors.Append(res.pts(row),last,lpt.all);
          for i in lpt'range loop
            second(i) := lpt(i);
          end loop;
          done := done - 1;
        end if;
        if done = 0 then
          for j in 1..dim loop
            mat(row,j) := first(j) - second(j);
          end loop;
          rhs(row) := second(dim+1) - first(dim+1);
        end if;
        exit when (done = 0);
        tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
      end loop;
    end loop;
    lufac(mat,dim,ipvt,info);
    lusolve(mat,dim,ipvt,rhs);
    res.nor(1..dim) := rhs;
    res.nor(dim+1) := 1.0;
    if verbose
     then Test_Inner_Normal(lifsup,res.nor.all,lbl,fail);
    end if;
    return res;
  end Make_Mixed_Cell;

  function Is_In ( labels : Standard_Integer_Vectors.Vector;
                   lbfirst,lblast,idx : integer32 ) return boolean is
  begin
    for k in lbfirst..lblast loop
      if labels(k) = idx
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Make_Mixed_Cell
              ( dim : integer32;
                mix,lbl : Standard_Integer_Vectors.Vector;
                lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                verbose : boolean := true )
              return Mixed_Cell is

    res : Mixed_Cell;
    tmp,last : Lists_of_Floating_Vectors.List;
    lpt : Standard_Floating_Vectors.Link_to_Vector;
    first : Standard_Floating_Vectors.Vector(1..dim+1);
    mat : Standard_Floating_Matrices.Matrix(1..dim,1..dim);
    rhs : Standard_Floating_Vectors.Vector(1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info,pntidx,lblidx,rowidx : integer32;
    fail,atfirstpt : boolean;

  begin
    res.nor := new Standard_Floating_Vectors.Vector'(1..dim+1 => 0.0);
    res.pts := new Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    res.sub := null;
    lblidx := 1;
    rowidx := 0;
    for k in mix'range loop
      tmp := lifsup(k);
      atfirstpt := true;
      pntidx := 0;
      while not Lists_of_Floating_Vectors.Is_Null(tmp) loop
        pntidx := pntidx + 1;
        if Is_In(lbl,lblidx,lblidx+mix(k),pntidx) then
          lpt := Lists_of_Floating_Vectors.Head_Of(tmp);
          if atfirstpt then
            for i in lpt'range loop
              first(i) := lpt(i);
            end loop;
            atfirstpt := false;
          else
            rowidx := rowidx + 1;   
            for j in 1..dim loop
              mat(rowidx,j) := first(j) - lpt(j);
            end loop;
            rhs(rowidx) := lpt(dim+1) - first(dim+1);
          end if;
          Lists_of_Floating_Vectors.Append(res.pts(k),last,lpt.all);
        end if;
        tmp := Lists_of_Floating_Vectors.Tail_Of(tmp);
      end loop;
      lblidx := lblidx + mix(k) + 1;
    end loop;
    lufac(mat,dim,ipvt,info);
    lusolve(mat,dim,ipvt,rhs);
    res.nor(1..dim) := rhs;
    res.nor(dim+1) := 1.0;
    if verbose
     then Test_Inner_Normal(mix,lifsup,res.nor.all,lbl,fail);
    end if;
    return res;
  end Make_Mixed_Cell;

  procedure Make_Mixed_Cell
              ( mic : in out Mixed_Cell; dim : in integer32;
                mix,lbl : in Standard_Integer_Vectors.Vector;
                lifsup : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mat : in out Standard_Floating_Matrices.Matrix;
                rhs : in out Standard_Floating_Vectors.Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                work : in out Standard_Floating_Vectors.Vector;
                verbose : in boolean := true ) is

    lifptr,micptr : Lists_of_Floating_Vectors.List;
    liflpt,miclpt : Standard_Floating_Vectors.Link_to_Vector;
    info,pntidx,lblidx,rowidx : integer32;
    fail,atfirstpt : boolean;

  begin
    lblidx := 1;
    rowidx := 0;
    for k in mix'range loop
      lifptr := lifsup(k);
      micptr := mic.pts(k);
      atfirstpt := true;
      pntidx := 0;
      while not Lists_of_Floating_Vectors.Is_Null(lifptr) loop
        pntidx := pntidx + 1;
        if Is_In(lbl,lblidx,lblidx+mix(k),pntidx) then
          liflpt := Lists_of_Floating_Vectors.Head_Of(lifptr);
          if atfirstpt then
            for i in liflpt'range loop
              work(i) := liflpt(i);
            end loop;
            atfirstpt := false;
          else
            rowidx := rowidx + 1;   
            for j in 1..dim loop
              mat(rowidx,j) := work(j) - liflpt(j);
            end loop;
            rhs(rowidx) := liflpt(dim+1) - work(dim+1);
          end if;
          miclpt := Lists_of_Floating_Vectors.Head_Of(micptr);
          for i in liflpt'range loop
            miclpt(i) := liflpt(i); -- copy coordinates of lifted point
          end loop;
          micptr := Lists_of_Floating_Vectors.Tail_Of(micptr);
        end if;
        lifptr := Lists_of_Floating_Vectors.Tail_Of(lifptr);
      end loop;
      lblidx := lblidx + mix(k) + 1;
    end loop;
    lufac(mat,dim,ipvt,info);
    lusolve(mat,dim,ipvt,rhs);
    mic.nor(1..dim) := rhs;
    mic.nor(dim+1) := 1.0;
    if verbose
     then Test_Inner_Normal(mix,lifsup,mic.nor.all,lbl,fail);
    end if;
  end Make_Mixed_Cell;

  function Make_Mixed_Cells
             ( dim : integer32;
               labels : Lists_of_Integer_Vectors.List;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               verbose : boolean := true )
             return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Lists_of_Integer_Vectors.List := labels;
    lbl : Standard_Integer_Vectors.Link_to_Vector;
    mic : Mixed_Cell;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      lbl := Lists_of_Integer_Vectors.Head_Of(tmp);
      mic := Make_Mixed_Cell(dim,lbl.all,lifsup,verbose);
      Append(res,res_last,mic);
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Make_Mixed_Cells;

  function Make_Mixed_Cells
             ( dim : integer32;
               mix : Standard_Integer_Vectors.Vector;
               labels : Lists_of_Integer_Vectors.List;
               lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               verbose : boolean := true )
             return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Lists_of_Integer_Vectors.List := labels;
    lbl : Standard_Integer_Vectors.Link_to_Vector;
    mic : Mixed_Cell;

  begin
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      lbl := Lists_of_Integer_Vectors.Head_Of(tmp);
      mic := Make_Mixed_Cell(dim,mix,lbl.all,lifsup,verbose);
      Append(res,res_last,mic);
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Make_Mixed_Cells;

end DEMiCs_Output_Convertors;
