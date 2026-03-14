with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Standard_Random_Vectors;
with DCMPLX_VecVecs_Container;

package body Test_DCMPLX_VecVecs_Container is

  procedure Prompt_Dimensions is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    DCMPLX_VecVecs_Container.Initialize(dim);
    put("-> the size : "); put(DCMPLX_VecVecs_Container.size,1); new_line;
    declare
      sizes : Standard_Integer_Vectors.Vector(1..dim) := (1..dim => 0);
    begin
      for i in 1..dim loop
        put("Give size of array "); put(i,1); put(" : ");
        get(sizes(i));
      end loop;
      DCMPLX_VecVecs_Container.Initialize(sizes);
    end;
    for i in 1..dim loop
      put("-> size of array "); put(i,1); put(" : ");
      put(DCMPLX_VecVecs_Container.size(i),1); new_line;
    end loop;  
  end Prompt_Dimensions;

  procedure Add_Random_Vectors ( dim : in integer32 ) is

    data : Standard_Complex_Vectors.Vector(1..dim);

  begin
    for i in 1..DCMPLX_VecVecs_Container.size loop
      for j in 1..DCMPLX_VecVecs_Container.size(i) loop
        data := Standard_Random_Vectors.Random_Vector(1,dim);
        put("-> adding "); put(j,1); put(" to component "); put(i,1); 
        put_line(" : "); put_line(data);
        DCMPLX_VecVecs_Container.Store_Copy(i,j,data);
      end loop;
    end loop;
  end Add_Random_Vectors;

  procedure Get_Vectors is

    data : Standard_Complex_Vectors.Link_to_Vector;

  begin
    for i in 1..DCMPLX_VecVecs_Container.size loop
      for j in 1..DCMPLX_VecVecs_Container.size(i) loop
        data := DCMPLX_VecVecs_Container.Get(i,j);
        put("-> vector "); put(j,1); put(" of component "); put(i,1); 
        put_line(" : "); put_line(data);
      end loop;
    end loop;
  end Get_Vectors;

  procedure Main is

    size : integer32 := 0;

  begin
    new_line;
    put_line("Testing the double vectors of vectors container ...");
    new_line;
    Prompt_Dimensions;
    new_line;
    put("Give the size of the vectors : "); get(size);
    Add_Random_Vectors(size);
    new_line;
    put_line("Retrieving the vectors ...");
    Get_Vectors;
    DCMPLX_VecVecs_Container.Clear;
  end Main;

end Test_DCMPLX_VecVecs_Container;
