with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
with Arrays_of_Floating_Vector_Lists;
with Supports_of_Polynomial_Systems;
with Floating_Mixed_Subdivisions_io;
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
                stlb : in double_float := 0.0;
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
  end Compute_Mixed_Volume;

  procedure Compute_Mixed_Labels
              ( data : in demics_input_data.class_dataSet.dataSet;
                mixvol : out integer32; seednbr : in integer32 := 0;
                stlb : in double_float := 0.0;
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
  end Compute_Mixed_Labels;

-- EXPORTED OPERATIONS :

  function Mixed_Volume
             ( p : Laur_Sys; seednbr : integer32 := 0;
               userlifting : boolean := false;
               vrblvl : integer32 := 0 ) return integer32 is

    res : integer32 := -1;
    data : DEMiCs_Input_Data.class_dataSet.dataSet;
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(p'range);
    mix : Standard_Integer_Vectors.Link_to_Vector;
    uselif : Standard_Floating_Vectors.Link_to_Vector := null;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.mixed_volume ...");
    end if;
    if not userlifting then
      data := DEMiCs_Translated_Setup.Make_Data(p,false,vrblvl-1);
    else
      sup := Supports_of_Polynomial_Systems.Create(p);
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
    Compute_Mixed_Volume(data,res,seednbr,0.0,uselif,vrblvl-1);
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
    mix : Standard_Integer_Vectors.Link_to_Vector;
    uselif : Standard_Floating_Vectors.Link_to_Vector := null;
    stlb : double_float := 0.0;

  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.mixed_volume ...");
    end if;
    if (not stablemv) and (not userlifting) then
      data := DEMiCs_Translated_Setup.Make_Data(p,false,vrblvl-1);
    else
      sup := Supports_of_Polynomial_Systems.Create(p);
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
    Compute_Mixed_Volume(data,res,seednbr,stlb,uselif,vrblvl-1);
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
    mix : Standard_Integer_Vectors.Link_to_Vector;
    uselif : Standard_Floating_Vectors.Link_to_Vector := null;
    stlb : double_float := 0.0;

  begin
    if vrblvl > 0 then
      put("-> in DEMiCs_Translated.mixed_labels, seednbr : ");
      put(seednbr,1); put_line(" ...");
    end if;
    if (not stablemv) and (not userlifting) then
      data := DEMiCs_Translated_Setup.Make_Data(p,true,vrblvl-1);
    else
      sup := Supports_of_Polynomial_Systems.Create(p);
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
    Compute_Mixed_Labels(data,res,seednbr,stlb,uselif,vrblvl-1);
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
    mix : Standard_Integer_Vectors.Link_to_Vector;
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
    Compute_Mixed_Labels(data,res,seednbr,0.0,uselif,vrblvl-1);
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
      lifsup := DEMiCs_Translated_Setup.Apply_Lifting(sup,lft,vrblvl-1);
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

  procedure Clear ( vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in DEMiCs_Translated.clear ...");
    end if;
    DEMiCs_MVC.class_mvc.delete_mvc(ptr2MVC);
  end Clear;

end DEMiCs_Translated;
