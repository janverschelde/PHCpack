with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;
with Mixed_Volume_Computation;
with demics_input_data;
with demics_mvc;

package body DEMiCs_Translated is

  procedure Make_Supports 
              ( data : in out demics_input_data.class_dataSet.dataSet;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    mixidx,supidx,size : integer32;
    tmp : Lists_of_Integer_Vectors.List;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.make_supports ...");
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

  function Make_Data ( p : Poly_sys; vrblvl : integer32 := 0 ) 
             return demics_input_data.class_dataSet.dataSet is

  -- DESCRIPTION :
  --   Returns the data object for input to DEMiCs.

    res : demics_input_data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range)
        := Supports_of_Polynomial_Systems.Create(p);
    mix,prm : Standard_Integer_Vectors.Link_to_Vector;
    nbr : natural32;
    idx : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in DEMiCs_Translated.make_data ...");
      put_line("the support sets : "); put(sup);
    end if;
    Mixed_Volume_Computation.Compute_Mixture(sup,mix,prm);
    if vrblvl > 0 then
      put("number of different supports : "); put(mix'last,1); new_line;
      put("type of mixture : "); put(mix);
      put(", permutation : "); put(prm); new_line;
    end if;
    res.dim := p'last;
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
    Standard_Integer_Vectors.Clear(mix);
    Standard_Integer_Vectors.Clear(prm);
    Arrays_of_Integer_Vector_Lists.Deep_Clear(sup);
    return res;
  end Make_Data;

  procedure Compute_Mixed_Volume
              ( data : in demics_input_data.class_dataSet.dataSet;
                mixvol : out integer32;
                vrblvl : in integer32 := 0 ) is

    use demics_mvc;

    ptr2MVC : constant class_mvc.Link_to_mvc
            := new class_mvc.mvc'(class_mvc.new_mvc);
    seed : constant integer32 := Standard_Random_Numbers.Get_Seed;

  begin
    if vrblvl > 0
     then put("the seed : "); put(seed,1); new_line;
    end if;
    if vrblvl > 0
     then class_mvc.allocateAndIni(ptr2MVC,data,seed,1,vrblvl);
     else class_mvc.allocateAndIni(ptr2MVC,data,seed,0);
    end if;
    class_mvc.Enum(ptr2MVC,vrblvl);
    mixvol := integer32(ptr2MVC.the_Simplex.mixedvol);
  end Compute_Mixed_Volume;

  function Mixed_Volume
             ( p : Poly_Sys; vrblvl : integer32 := 0 ) return integer32 is

    res : integer32;
    data : demics_input_data.class_dataSet.dataSet;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.mixed_volume ...");
    end if;
    data := Make_Data(p,vrblvl-1);
    if vrblvl > 0 then
      put_line("the preamble of the DEMiCs input data : ");
      demics_input_data.class_dataSet.info_preamble(data);
      put_line("the supports of the DEMiCs input data : ");
      demics_input_data.class_dataSet.info_supports(data);
    end if;
    Compute_Mixed_Volume(data,res,vrblvl-1);
    return res;
  exception
    when others => return -1;
  end Mixed_Volume;

end DEMiCs_Translated;
