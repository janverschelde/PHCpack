with Ada.text_io;                       use Ada.text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Communications_with_User;
with demics_input_data;
with demics_input_main;
with demics_mvc;

procedure demics_main is

  procedure Compute_Mixed_Volume
              ( data : in demics_input_data.class_dataSet.dataSet;
                vrblvl : in integer32 := 0 ) is

    use demics_mvc;

    ptr2MVC : constant class_mvc.Link_to_mvc
            := new class_mvc.mvc'(class_mvc.new_mvc);
    seed : constant integer32 := Standard_Random_Numbers.Get_Seed;

  begin
    if vrblvl > 0
     then put("the seed : "); put(seed,1); new_line;
    end if;
    class_mvc.allocateAndIni(ptr2MVC,data,seed,1,vrblvl);
    class_mvc.Enum(ptr2MVC,vrblvl);
  end Compute_Mixed_Volume;
 
  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the input file name and
  --   launches the mixed volume computation.

    use demics_input_data.class_dataSet;

    ans : character;
    data : dataSet;
    fail : boolean;
    vrblvl : integer32 := 0;

  begin
    put("Intermediate output wanted ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    new_line;
    if ans /= 'y' then
      demics_input_main.read_data_from_file(data,fail,0);
    else
      demics_input_main.read_data_from_file(data,fail,1);
      demics_input_data.class_dataSet.info_preamble(data);
      demics_input_data.class_dataSet.info_supports(data);
      vrblvl := 99;
    end if;
    if not fail
     then Compute_Mixed_Volume(data,vrblvl);
    end if;
  end Main;

begin
  Main;
end demics_main;
