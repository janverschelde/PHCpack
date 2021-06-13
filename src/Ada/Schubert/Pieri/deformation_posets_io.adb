with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_VecMats;           use Standard_Complex_VecMats;

package body Deformation_Posets_io is

  procedure put_size ( poset : in Array_of_Array_of_VecMats ) is
  begin
    put_size(Standard_Output,poset);
  end put_size;

  procedure put_size 
              ( file : in file_type; poset : in Array_of_Array_of_VecMats ) is

    np : natural32;
    lavm : Link_to_VecMat;

  begin
    if poset'last < 10
     then np := 1;
     else np := 2;
    end if;
    for i in poset'range loop
      put(file,"n = "); put(file,i,np); put(file," : ");
      if poset(i) /= null then
        for j in poset(i)'range loop
          lavm := poset(i)(j);
          if lavm = null
           then put(file," 0");
           else put(file," "); put(file,integer32(lavm'length),1);
          end if;
        end loop;
      end if;
      new_line(file);
    end loop;
  end put_size;

end Deformation_Posets_io;
