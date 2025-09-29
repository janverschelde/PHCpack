with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Lists_of_Floating_Vectors;
with Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;
with DEMiCs_Simplex;
with DEMiCs_Output_Cells;

package body DEMiCs_Translated_Setup is

  procedure Make_Supports 
              ( data : in out DEMiCs_Input_Data.class_dataSet.dataSet;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    mixidx,supidx,size : integer32;
    tmp : Lists_of_Integer_Vectors.List;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated_Setup.make_supports ...");
    end if;
    data.termSumNum := 0;
    data.termStart := new Standard_Integer_Vectors.Vector(0..data.supN);
    data.termStart(0) := 0;
    for k in 0..data.supN-1 loop
      data.termSumNum := data.termSumNum + data.termSet(k);
      data.termStart(k+1) := data.termSumNum;
    end loop;
    if vrblvl > 0
     then put("termStart : "); put(data.termStart); new_line;
    end if;
    data.termStart(data.supN) := data.termSumNum;
    size := data.termSumNum*data.dim;
    data.support := new Standard_Floating_Vectors.Vector(0..size-1);
    supidx := 0; -- index in data.support
    mixidx := 1; -- index in sup
    for i in mix'range loop
      tmp := sup(mixidx);
      for j in 1..data.termSet(i-1) loop
        lpt := Lists_of_Integer_Vectors.Head_Of(tmp);
        for k in lpt'range loop
          data.support(supidx) := double_float(lpt(k));
          supidx := supidx + 1;
        end loop;
        tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
      end loop;
      mixidx := mixidx + mix(i);
    end loop;
  end Make_Supports;

  procedure Make_Data
              ( res : out DEMiCs_Input_Data.class_dataSet.dataSet;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                storemix : in boolean; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Returns the data object for input to DEMiCs.

    mix,prm : Standard_Integer_Vectors.Link_to_Vector;
    nbr : natural32;
    idx : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in DEMiCs_Translated_Setup.make_data ...");
    end if;
    Mixed_Volume_Computation.Compute_Mixture(sup,mix,prm);
    if storemix
     then DEMiCs_Output_Cells.Store_Dimension_and_Mixture(sup'last,mix);
    end if;
    if vrblvl > 0 then
      put("number of different supports : "); put(mix'last,1); new_line;
      put("type of mixture : "); put(mix);
      put(", permutation : "); put(prm); new_line;
    end if;
    res.dim := sup'last;
    res.supN := mix'last;
    res.supType := new Standard_Integer_Vectors.Vector(0..mix'last-1);
    res.termSet := new Standard_Integer_Vectors.Vector(0..mix'last-1);
    idx := 1;
    for i in mix'range loop
      res.supType(i-1) := mix(i);
      nbr := Lists_of_Integer_Vectors.Length_Of(sup(idx));
      res.termSet(i-1) := integer32(nbr);
      idx := idx + mix(i);
    end loop;
    res.termMax := res.termSet(0);
    res.typeMax := res.supType(0);
    for i in 1..res.supType'last loop
      if res.termSet(i) > res.termMax
       then res.termMax := res.termSet(i);
      end if;
      if res.supType(i) > res.typeMax
       then res.typeMax := res.supType(i);
      end if;
    end loop;
    Make_Supports(res,sup,mix,vrblvl-1);
    if not storemix
     then Standard_Integer_Vectors.Clear(mix);
    end if;
    Standard_Integer_Vectors.Clear(prm);
  end Make_Data;

  function Make_Data ( p : Poly_sys; storemix : boolean;
                       vrblvl : integer32 := 0 ) 
                     return DEMiCs_Input_Data.class_dataSet.dataSet is

    res : DEMiCs_Input_Data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);

  begin
    if vrblvl > 0 then
      put_line("-> in DEMiCs_Translated_Setup.make_data for polynomials ...");
      put_line("the support sets : "); put(sup);
    end if;
    Make_Data(res,sup,storemix,vrblvl-1);
    Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    return res;
  end Make_Data;

  function Make_Data ( p : Laur_sys; storemix : boolean;
                       vrblvl : integer32 := 0 ) 
                     return DEMiCs_Input_Data.class_dataSet.dataSet is

    res : DEMiCs_Input_Data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);

  begin
    if vrblvl > 0 then
      put("-> in DEMiCs_Translated_Setup.make_data for Laurent polynomials");
      put_line(" ...");
      put_line("the support sets : "); put(sup);
    end if;
    Make_Data(res,sup,storemix,vrblvl-1);
    Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    return res;
  end Make_Data;

  function Extract_Support
             ( dim : integer32;
               sup : Standard_Floating_Vectors.Link_to_Vector;
               vrblvl : integer32 := 0 )
             return Lists_of_Integer_Vectors.List is

    res : Lists_of_Integer_Vectors.List;
    res_last : Lists_of_Integer_Vectors.List;
    idx : integer32 := sup'first;
    pnt : Standard_Integer_Vectors.Vector(1..dim);

  begin
    if vrblvl > 0 then
      put_line("-> in DEMiCs_Translated_Setup.extract_support ...");
      put_line("the support vector : ");
      put(sup); new_line;
    end if;
    loop
      for i in 1..dim loop
        pnt(i) := integer32(sup(idx));
        idx := idx + 1;
      end loop;
      Lists_of_Integer_Vectors.Append(res,res_last,pnt);
      exit when (idx >= sup'last);
    end loop;
    return res;
  end Extract_Support;

  function Extract_Supports
             ( ptr2MVC : DEMiCs_MVC.class_mvc.Link_to_mvc;
               vrblvl : integer32 := 0 )
             return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    sup : constant Standard_Floating_VecVecs.Link_to_VecVec
        := ptr2MVC.the_Simplex.oriSupp; -- original supports start at 0
    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..sup'last+1);
    dim : constant integer32 := ptr2MVC.dim;

  begin
    if vrblvl > 0 then
      put("-> in DEMiCs_Translated.extract_supports, dim : ");
      put(dim,1); put_line(" ...");
      put_line("The original support sets : ");
      DEMiCs_Simplex.class_simplex.info_oriSup(ptr2MVC.the_Simplex);
    end if;
    for i in res'range loop
      res(i) := Extract_Support(dim,sup(i-1),vrblvl-1);
    end loop;
    return res;
  end Extract_Supports;

  function Apply_Lifting
             ( sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists;
               lft : Standard_Floating_Vectors.Link_to_Vector;
               vrblvl : integer32 := 0 )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

  -- DESCRIPTION :
  --   Extends the supports sup with the lifting in lft.

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range);
    res_last : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range);
    idx : integer32 := lft'first;
    ipt : Standard_Integer_Vectors.Link_to_Vector;
    fpt : Standard_Floating_Vectors.Link_to_Vector;
    p2sup : Lists_of_Integer_Vectors.List;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.apply_lifting ...");
    end if;
    for i in sup'range loop
      p2sup := sup(i);
      while not Lists_of_Integer_Vectors.Is_Null(p2sup) loop
        ipt := Lists_of_Integer_Vectors.Head_Of(p2sup);
        fpt := new Standard_Floating_Vectors.Vector(ipt'first..ipt'last+1);
        for j in ipt'range loop
          fpt(j) := double_float(ipt(j));
        end loop;
        fpt(fpt'last) := lft(idx);
        idx := idx + 1;
        Lists_of_Floating_Vectors.Append(res(i),res_last(i),fpt);
        p2sup := Lists_of_Integer_Vectors.Tail_Of(p2sup);
      end loop;
    end loop;
    return res;
  end Apply_Lifting;

end DEMiCs_Translated_Setup;
