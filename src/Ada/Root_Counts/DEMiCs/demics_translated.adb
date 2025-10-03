with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_VecVecs_io;      use Standard_Floating_VecVecs_io;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Supports_of_Polynomial_Systems;
with Floating_Lifting_Functions;
with Floating_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;
with DEMiCs_Input_Data;
with DEMiCs_MVC;
with DEMiCs_Output_Cells;
with DEMiCs_Output_Convertors;
with DEMiCs_Translated_Setup;

package body DEMiCs_Translated is

  ptr2MVC : DEMiCs_MVC.class_mvc.Link_to_mvc;

  procedure Compute_Mixed_Volume
              ( data : in demics_input_data.class_dataSet.dataSet;
                mixvol : out integer32; seednbr : in integer32 := 0;
                uselif : in Standard_Floating_Vectors.Link_to_Vector := null;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Makes an mvc object and dynamically enumerates the mixed cells.

    use DEMiCs_MVC;
    seed : integer32;

  begin
    if vrblvl > 0 then
      put("-> in DEMiCs_Translated.compute_mixed_volume, seednbr : ");
      put(seednbr,1); put_line(" ...");
    end if;
    if seednbr = 0
     then seed := Standard_Random_Numbers.Get_Seed;
     else seed := seednbr;
    end if;
    ptr2MVC := new class_mvc.mvc'(class_mvc.new_mvc);
    if vrblvl > 0
     then class_mvc.allocateAndIni(ptr2MVC,data,seed,1,uselif,vrblvl);
     else class_mvc.allocateAndIni(ptr2MVC,data,seed,0,uselif);
    end if;
    class_mvc.Enum(ptr2MVC,vrblvl);
    mixvol := integer32(ptr2MVC.the_Simplex.mixedvol);
    DEMiCs_Output_Cells.mixed_volume := mixvol;
  end Compute_Mixed_Volume;

  procedure Compute_Mixed_Labels
              ( data : in demics_input_data.class_dataSet.dataSet;
                mixvol : out integer32; seednbr : in integer32 := 0;
                uselif : in Standard_Floating_Vectors.Link_to_Vector := null;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Stores the labels to the mixed cells.

    use DEMiCs_MVC;
    seed : integer32;

  begin
    if vrblvl > 0 then
      put("-> in DEMiCs_Translated.compute_mixed_labels, seednbr : ");
      put(seednbr,1); put_line(" ...");
    end if;
    if seednbr = 0
     then seed := Standard_Random_Numbers.Get_Seed;
     else seed := seednbr;
    end if;
    ptr2MVC := new class_mvc.mvc'(class_mvc.new_mvc);
    if vrblvl > 0
     then class_mvc.allocateAndIni(ptr2MVC,data,seed,2,uselif,vrblvl);
     else class_mvc.allocateAndIni(ptr2MVC,data,seed,2,uselif);
    end if;
    class_mvc.Enum(ptr2MVC,vrblvl);
    mixvol := integer32(ptr2MVC.the_Simplex.mixedvol);
    DEMiCs_Output_Cells.mixed_volume := mixvol;
  end Compute_Mixed_Labels;

-- EXPORTED OPERATIONS :

  function Mixed_Volume
             ( p : Laur_Sys; seednbr : integer32 := 0;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32 is

    res : integer32 := -1;
    data : DEMiCs_Input_Data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    mix,prm : Standard_Integer_Vectors.Link_to_Vector;
    uselif : Standard_Floating_Vectors.Link_to_Vector := null;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.mixed_volume ...");
    end if;
    if not userlifting then
      data := DEMiCs_Translated_Setup.Make_Data(p,false,vrblvl-1);
    else
      sup := Supports_of_Polynomial_Systems.Create(p);
      Mixed_Volume_Computation.Compute_Mixture(sup,mix,prm);
      DEMiCs_Translated_Setup.Make_Data(data,sup,mix,vrblvl-1);
      DEMiCs_Output_Cells.Store_Dimension_and_Mixture(sup'last,mix);
      uselif := DEMiCs_Translated_Setup.User_Lifting(mix,sup);
    end if;
    if vrblvl > 0 then
      put_line("the preamble of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_preamble(data);
      put_line("the supports of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_supports(data);
    end if;
    Compute_Mixed_Volume(data,res,seednbr,uselif,vrblvl-1);
    return res;
  end Mixed_Volume;

  function Mixed_Volume
             ( p : Poly_Sys; seednbr : integer32 := 0;
               stablemv : boolean := false;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32 is

    res : integer32 := -1;
    data : DEMiCs_Input_Data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    mix,prm : Standard_Integer_Vectors.Link_to_Vector;
    uselif : Standard_Floating_Vectors.Link_to_Vector := null;
    stlb : double_float := 0.0;
    dim : constant integer32 := sup'last;
    nbrad : integer32;
    added : Standard_Integer_Vectors.Vector(sup'range);

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.mixed_volume ...");
    end if;
    if (not stablemv) and (not userlifting) then
      data := DEMiCs_Translated_Setup.Make_Data(p,false,vrblvl-1);
    else
      sup := Supports_of_Polynomial_Systems.Create(p);
      Mixed_Volume_Computation.Compute_Mixture(sup,mix,prm);
      DEMiCs_Output_Cells.Store_Dimension_and_Mixture(sup'last,mix);
      if userlifting then
        uselif := DEMiCs_Translated_Setup.User_Lifting(mix,sup);
      elsif stablemv then
        stlb := Floating_Lifting_Functions.Lifting_Bound(p);
        DEMiCs_Translated_Setup.Add_Artificial_Origins(dim,sup,nbrad,added);
        uselif := DEMiCs_Translated_Setup.Random_Lifting(mix,sup,stlb,added);
        DEMiCs_Output_Cells.stable := true;
        DEMiCs_Output_Cells.stlb := stlb;
      end if;
      DEMiCs_Translated_Setup.Make_Data(data,sup,mix,vrblvl-1);
    end if;
    if vrblvl > 0 then
      put_line("the preamble of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_preamble(data);
      put_line("the supports of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_supports(data);
    end if;
    Compute_Mixed_Volume(data,res,seednbr,uselif,vrblvl-1);
    return res;
  end Mixed_Volume;

  function Mixed_Labels
             ( p : Poly_Sys; monitor : boolean := true;
               seednbr : integer32 := 0;
               stablemv : boolean := false;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32 is

    res : integer32 := -1;
    data : DEMiCs_Input_Data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    mix,prm : Standard_Integer_Vectors.Link_to_Vector;
    uselif : Standard_Floating_Vectors.Link_to_Vector := null;
    stlb : double_float := 0.0;
    dim : constant integer32 := sup'last;
    nbrad : integer32;
    added : Standard_Integer_Vectors.Vector(sup'range);

  begin
    if vrblvl > 0 then
      put("-> in DEMiCs_Translated.mixed_labels, seednbr : ");
      put(seednbr,1); put_line(" ...");
    end if;
    if (not stablemv) and (not userlifting) then
      data := DEMiCs_Translated_Setup.Make_Data(p,true,vrblvl-1);
    else
      sup := Supports_of_Polynomial_Systems.Create(p);
      Mixed_Volume_Computation.Compute_Mixture(sup,mix,prm);
      DEMiCs_Output_Cells.Store_Dimension_and_Mixture(sup'last,mix);
      if userlifting then
        uselif := DEMiCs_Translated_Setup.User_Lifting(mix,sup);
      elsif stablemv then
        stlb := Floating_Lifting_Functions.Lifting_Bound(p);
        DEMiCs_Translated_Setup.Add_Artificial_Origins(dim,sup,nbrad,added);
        uselif := DEMiCs_Translated_Setup.Random_Lifting(mix,sup,stlb,added);
        DEMiCs_Output_Cells.stable := true;
        DEMiCs_Output_Cells.stlb := stlb;
      end if;
      DEMiCs_Translated_Setup.Make_Data(data,sup,mix,vrblvl-1);
    end if;
    if vrblvl > 0 then
      put_line("the preamble of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_preamble(data);
      put_line("the supports of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_supports(data);
    end if;
    DEMiCs_Output_Cells.monitor := monitor;
    Compute_Mixed_Labels(data,res,seednbr,uselif,vrblvl-1);
    return res;
  end Mixed_Labels;

  function Mixed_Labels
             ( p : Laur_Sys; monitor : boolean := true;
               seednbr : integer32 := 0;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32 is

    res : integer32 := -1;
    data : DEMiCs_Input_Data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    mix,prm : Standard_Integer_Vectors.Link_to_Vector;
    uselif : Standard_Floating_Vectors.Link_to_Vector := null;

  begin
    if vrblvl > 0 then
      put("-> in DEMiCs_Translated.mixed_labels, seednbr : ");
      put(seednbr,1); put_line(" ...");
    end if;
    if not userlifting then
      data := DEMiCs_Translated_Setup.Make_Data(p,true,vrblvl-1);
    else
      sup := Supports_of_Polynomial_Systems.Create(p);
      Mixed_Volume_Computation.Compute_Mixture(sup,mix,prm);
      DEMiCs_Translated_Setup.Make_Data(data,sup,mix,vrblvl-1);
      DEMiCs_Output_Cells.Store_Dimension_and_Mixture(sup'last,mix);
      uselif := DEMiCs_Translated_Setup.User_Lifting(mix,sup);
    end if;
    if vrblvl > 0 then
      put_line("the preamble of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_preamble(data);
      put_line("the supports of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_supports(data);
    end if;
    DEMiCs_Output_Cells.monitor := monitor;
    Compute_Mixed_Labels(data,res,seednbr,uselif,vrblvl-1);
    return res;
  end Mixed_Labels;

  function Mixed_Cells ( vrblvl : integer32 := 0 ) return Mixed_Subdivision is

    res : Mixed_Subdivision;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.mixed_cells ...");
    end if;
    declare
      dim : constant integer32 := ptr2MVC.dim;
      mix : constant Standard_Integer_Vectors.Link_to_Vector
          := DEMiCs_Output_Cells.Get_Mixture;
      sup : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists
          := DEMiCs_Translated_Setup.Extract_Supports(ptr2MVC,vrblvl-1);
      lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(sup'range);
      lft : constant Standard_Floating_Vectors.Link_to_Vector
          := ptr2MVC.the_Simplex.lifting;
      labels : constant Lists_of_Integer_Vectors.List
             := DEMiCs_Output_Cells.Retrieve_Cell_Indices;

    begin
      if vrblvl > 0
       then put_line("the extracted supports : "); put(sup);
      end if;
      lifsup := DEMiCs_Translated_Setup.Apply_Lifting(mix,sup,lft,vrblvl-1);
      if vrblvl > 0 then
        put_line("the lifted supports : ");
        Floating_Mixed_Subdivisions_io.put(lifsup);
        res := DEMiCs_Output_Convertors.Make_Mixed_Cells
                 (dim,mix.all,labels,lifsup);
      else
        res := DEMiCs_Output_Convertors.Make_Mixed_Cells
                 (dim,mix.all,labels,lifsup,false);
      end if;
    end;
    return res;
  end Mixed_Cells;

-- ORIGINAL INTERFACE :

  procedure Extract_Supports 
               ( p : in Poly_Sys;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 verbose : in boolean := true ) is

    perm : Standard_Integer_Vectors.Link_to_Vector;

  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
    end if;
  end Extract_Supports;

  procedure Extract_Supports 
               ( p : in Poly_Sys;
                 mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                 supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                 verbose : in boolean := true ) is
  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
      put("The permutation : "); put(perm.all); new_line;
    end if;
  end Extract_Supports;

  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is

    perm : Standard_Integer_Vectors.Link_to_Vector;

  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
    end if;
  end Extract_Supports;

  procedure Extract_Supports 
              ( p : in Laur_Sys;
                mix,perm : out Standard_Integer_Vectors.Link_to_Vector;
                supports : out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                verbose : in boolean := true ) is
  begin
    supports := Supports_of_Polynomial_Systems.Create(p);
    Mixed_Volume_Computation.Compute_Mixture(supports,mix,perm);
    if verbose then
      put_line("The supports : "); put(supports);
      put("The mixture type : "); put(mix.all); new_line;
      put("The permutation : "); put(perm.all); new_line;
    end if;
  end Extract_Supports;

  function Random_Lifting
             ( mix : Standard_Integer_Vectors.Link_to_Vector;
               sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
             return Standard_Floating_VecVecs.Link_to_VecVec is

    res : Standard_Floating_VecVecs.Link_to_VecVec;
    resrep : Standard_Floating_VecVecs.VecVec(mix'range);
    idx : integer32 := 1;
    len : integer32;

  begin
    for i in resrep'range loop
      len := integer32(Lists_of_Integer_Vectors.Length_Of(sup(idx)));
      declare
        vals : Standard_Floating_Vectors.Vector(1..len);
      begin
        for j in 1..len loop
          vals(j) := Standard_Random_Numbers.Random;
        end loop;
        resrep(i) := new Standard_Floating_Vectors.Vector'(vals);
      end;
      idx := idx + mix(i);
    end loop;
    res := new Standard_Floating_VecVecs.VecVec'(resrep);
    return res;
  end Random_Lifting;

  function Size ( v : Standard_Floating_VecVecs.Link_to_VecVec )
                return integer32 is

  -- DESCRIPTION :
  --   Returns the total number of values in v.

    res : integer32 := 0;

  begin
    for i in v'range loop
      res := res + v(i)'last;
    end loop;
    return res;
  end Size;

  function Flatten ( v : Standard_Floating_VecVecs.Link_to_VecVec )
                   return Standard_Floating_Vectors.Link_to_Vector is

    res : Standard_Floating_Vectors.Link_to_Vector;
    res_rep : Standard_Floating_Vectors.Vector(1..Size(v));
    idx : integer32 := 0;

  begin
    for i in v'range loop
      for j in v(i)'range loop
        idx := idx + 1;
        res_rep(idx) := v(i)(j);
      end loop;
    end loop;
    res := new Standard_Floating_Vectors.Vector'(res_rep);
    return res;
  end Flatten;

  procedure Call_DEMiCs
              ( mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                stlb : in double_float := 0.0;
                lft : in Standard_Floating_Vectors.Link_to_Vector := null;
                vrblvl : in integer32 := 0 ) is

    dim : constant integer32 := sup'last;
    data : DEMiCs_Input_Data.class_dataSet.dataSet;
    mv : integer32 := -1;
    nbrad : integer32;
    added : Standard_Integer_Vectors.Vector(mix'range);
    rndlif : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.call_DEMiCs ...");
    end if;
    if stlb /= 0.0 then
      DEMiCs_Translated_Setup.Add_Artificial_Origins(dim,sup,nbrad,added);
      rndlif := DEMiCs_Translated_Setup.Random_Lifting(mix,sup,stlb,added);
      DEMiCs_Output_Cells.stable := true;
      DEMiCs_Output_Cells.stlb := stlb;
    end if;
    DEMiCs_Translated_Setup.Make_Data(data,sup,mix,vrblvl-1);
    DEMiCs_Output_Cells.Store_Dimension_and_Mixture(dim,mix);
    if vrblvl > 0 then
      put_line("the preamble of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_preamble(data);
      put_line("the supports of the DEMiCs input data : ");
      DEMiCs_Input_Data.class_dataSet.info_supports(data);
    end if;
    if vrblvl > 0 -- be aware that it may have been already true
     then DEMiCs_Output_Cells.monitor := true;
    end if;
    if stlb = 0.0
     then Compute_Mixed_Labels(data,mv,0,lft,vrblvl-1);
     else Compute_Mixed_Labels(data,mv,0,rndlif,vrblvl-1);
    end if;
    if vrblvl > 0 then
      if stlb = 0.0
       then put("the mixed volume : "); put(mv,1); new_line;
       else put("the stable mixed volume : "); put(mv,1); new_line;
      end if;
    end if;
  end Call_DEMiCs;

  procedure Show_Output is

    lifting : constant Standard_Floating_VecVecs.Link_to_VecVec
            := DEMiCs_Output_Cells.Lifting_Values;
    cells : constant Lists_of_Integer_Vectors.List
          := DEMiCs_Output_Cells.Retrieve_Cell_Indices;
    tmp : Lists_of_Integer_Vectors.List := cells;
    lpt : Standard_Integer_Vectors.Link_to_Vector;
    mv : constant integer32 := DEMiCs_Output_Cells.mixed_volume;
    cnt : integer32 := 0;

  begin
    put_line("The lifting values :");
    put(lifting.all);
    put_line("The labels to the mixed cells :");
    while not Lists_of_Integer_Vectors.Is_Null(tmp) loop
      lpt := Lists_of_Integer_Vectors.Head_Of(tmp);
      cnt := cnt + 1;
      put("cell "); put(cnt,1); put(" : ");
      put(lpt); new_line;
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    put("The mixed volume : "); put(mv,1); new_line;
  end Show_Output;

  procedure Process_Output
              ( dim : in integer32;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                sup : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lif : out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mcc : out Mixed_Subdivision; vrblvl : in integer32 := 0 ) is

    lft : constant Standard_Floating_Vectors.Link_to_Vector
        := ptr2MVC.the_Simplex.lifting;
    labels : constant Lists_of_Integer_Vectors.List
           := DEMiCs_Output_Cells.Retrieve_Cell_Indices;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.process_output ...");
    end if;
    lif := DEMiCs_Translated_Setup.Apply_Lifting(mix,sup,lft,vrblvl-1);
    if vrblvl > 0 then
      put_line("the lifted supports : ");
      Floating_Mixed_Subdivisions_io.put(lif);
      mcc := DEMiCs_Output_Convertors.Make_Mixed_Cells
               (dim,mix.all,labels,lif,true);
    else
      mcc := DEMiCs_Output_Convertors.Make_Mixed_Cells
               (dim,mix.all,labels,lif,false);
    end if;
  end Process_Output;

  procedure Clear ( vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.clear ...");
    end if;
    DEMiCs_MVC.class_mvc.delete_mvc(ptr2MVC);
  end Clear;

end DEMiCs_Translated;
