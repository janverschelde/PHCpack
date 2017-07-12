with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with DoblDobl_Sample_Points;            use DoblDobl_Sample_Points;

package body DoblDobl_Sample_Grids is

-- THE STATE IS POLYNOMIAL OR LAURENT :

  procedure Set_Polynomial_Type ( laurent : in boolean ) is
  begin
    DoblDobl_Sample_Points.Set_Polynomial_Type(laurent);
    DoblDobl_Sample_Lists.Set_Polynomial_Type(laurent);
  end Set_Polynomial_Type;

-- CREATORS AS TYPE CONVERTORS :

  function Create ( grid : DoblDobl_Sample_Grid )
                  return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(1..integer32(Length_Of(grid)));
    tmp : DoblDobl_Sample_Grid := grid;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Create;

  procedure Create ( samples : in Array_of_DoblDobl_Sample_Lists;
                     grid,grid_last : in out DoblDobl_Sample_Grid ) is
  begin
    for i in samples'range loop
      Append(grid,grid_last,samples(i));
    end loop;
  end Create;

-- SAMPLERS :

  procedure Sample ( samples : in DoblDobl_Sample_List; nb : in natural32;
                     grid,grid_last : in out DoblDobl_Sample_Grid ) is

    tmp : DoblDobl_Sample_List := samples;

  begin
    while not Is_Null(tmp) loop
      declare
        spt : constant DoblDobl_Sample := Head_Of(tmp);
        sps,sps_last : DoblDobl_Sample_List;
      begin
        Sample(spt,nb,sps,sps_last);
        Append(grid,grid_last,sps);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Sample;

  procedure Sample ( file : in file_type; full_output : in boolean;
                     samples : in DoblDobl_Sample_List; nb : in natural32;
                     grid,grid_last : in out DoblDobl_Sample_Grid ) is

  -- DESCRIPTION :
  --   Creates a grid of sample points, generating from every point in
  --   the list nb new sample points.

    tmp : DoblDobl_Sample_List := samples;

  begin
    while not Is_Null(tmp) loop
      declare
        spt : constant DoblDobl_Sample := Head_Of(tmp);
        sps,sps_last : DoblDobl_Sample_List;
      begin
        Sample(file,full_output,spt,nb,sps,sps_last);
        Append(grid,grid_last,sps);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Sample;

-- DESTRUCTORS :

  procedure Shallow_Clear ( g : in out DoblDobl_Sample_Grid ) is
  begin
    Lists_of_DoblDobl_Sample_Lists.Clear
      (Lists_of_DoblDobl_Sample_Lists.List(g));
  end Shallow_Clear;

  procedure Deep_Clear ( g : in out DoblDobl_Sample_Grid ) is

    tmp : DoblDobl_Sample_Grid := g;

  begin
    while not Is_Null(tmp) loop
      declare
        samples : DoblDobl_Sample_List := Head_Of(tmp);
      begin
        Deep_Clear(samples);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(g);
  end Deep_Clear;

end DoblDobl_Sample_Grids;
