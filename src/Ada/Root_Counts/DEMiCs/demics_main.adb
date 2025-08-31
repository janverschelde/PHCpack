with text_io;                           use text_io;
with Communications_with_User;
with demics_input_data;
with demics_input_main;
with demics_mvc;

procedure demics_main is

  procedure Compute_Mixed_Volume
              ( data : in demics_input_data.class_dataSet.dataSet ) is

    use demics_mvc;

    ptr2MVC : constant class_mvc.Link_to_mvc := new class_mvc.mvc;

  begin
    class_mvc.allocateAndIni(ptr2MVC,data,1,1,99);
  end Compute_Mixed_Volume;
 
  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the input file name and
  --   launches the mixed volume computation.

    use demics_input_data.class_dataSet;

    ans : character;
    data : dataSet;
    fail : boolean;

  begin
    put("Intermediate output wanted ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y'
     then demics_input_main.read_data_from_file(data,fail,1);
     else demics_input_main.read_data_from_file(data,fail,0);
    end if;
    if not fail
     then Compute_Mixed_Volume(data);
    end if;
  end Main;

begin
  Main;
end demics_main;
