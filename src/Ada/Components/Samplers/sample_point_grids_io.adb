with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Sample_Point_Lists_io;             use Sample_Point_Lists_io;

package body Sample_Point_Grids_io is

  procedure get ( grid,grid_last : in out Standard_Sample_Grid ) is
  begin
    get(Standard_Input,grid,grid_last);
  end get;

  procedure get ( file : in file_type;
                  grid,grid_last : in out Standard_Sample_Grid ) is

    nb : natural32 := 0;

  begin
    get(file,nb);
    for i in 1..nb loop
      declare
        samples,samples_last : Standard_Sample_List;
      begin
        get(file,samples,samples_last);
        Append(grid,grid_last,samples);
      end;
    end loop;
  end get;

  procedure get ( grid,grid_last : in out Multprec_Sample_Grid ) is
  begin
    get(Standard_Input,grid,grid_last);
  end get;

  procedure get ( file : in file_type;
                  grid,grid_last : in out Multprec_Sample_Grid ) is

    nb : natural32 := 0;

  begin
    get(file,nb);
    for i in 1..nb loop
      declare
        samples,samples_last : Multprec_Sample_List;
      begin
        get(file,samples,samples_last);
        Append(grid,grid_last,samples);
      end;
    end loop;
  end get;

  procedure put ( grid : in Standard_Sample_Grid ) is
  begin
    put(Standard_Output,grid);
  end put;

  procedure put ( file : in file_type; grid : in Standard_Sample_Grid ) is

    tmp : Standard_Sample_Grid := grid;

  begin
    put(file,Length_Of(grid),1); new_line(file);
    while not Is_Null(tmp) loop
      declare
        samples : constant Standard_Sample_List := Head_Of(tmp);
        l,n,k : natural32;
        s : Standard_Sample; 
      begin
        if not Is_Null(samples) then
          l := Length_Of(samples);
          s := Head_Of(samples);
          n := natural32(Number_of_Variables(s));
          k := natural32(Number_of_Slices(s));
          put(file,l,n,k,samples);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( grid : in Multprec_Sample_Grid ) is
  begin
    put(Standard_Output,grid);
  end put;

  procedure put ( file : in file_type; grid : in Multprec_Sample_Grid ) is

    tmp : Multprec_Sample_Grid := grid;

  begin
    put(file,Length_Of(grid),1); new_line(file);
    while not Is_Null(tmp) loop
      declare
        samples : constant Multprec_Sample_List := Head_Of(tmp);
        l,n,k : natural32;
        s : Multprec_Sample;
      begin
        if not Is_Null(samples) then
          l := Length_Of(samples);
          s := Head_Of(samples);
          n := natural32(Number_of_Variables(s));
          k := natural32(Number_of_Slices(s));
          put(file,l,n,k,samples);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end put;

end Sample_Point_Grids_io;
