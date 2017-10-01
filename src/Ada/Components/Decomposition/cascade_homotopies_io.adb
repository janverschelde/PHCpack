with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;

package body Cascade_Homotopies_io is

  procedure Write_Super_Witness_Points
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Super_Witness_Points;

  procedure Write_Super_Witness_Points
              ( file : in file_type;
                sols : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Super_Witness_Points;

  procedure Write_Super_Witness_Points
              ( file : in file_type;
                sols : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

  begin
    if not Is_Null(sols) then
      new_line(file);
      put_line(file,"THE SOLUTIONS with zz = 0 :");
      put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    end if;
  end Write_Super_Witness_Points;

  function Append_ck ( name : string; k : natural32 ) return string is

    nbk : constant string := Convert(integer32(k));

  begin
     return name & "_sw" & nbk;
  end Append_ck;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32 ) is

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
    create(file,out_file,filename);
    put_line(file,p);
    Write_Super_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32 ) is

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
    create(file,out_file,filename);
    put_line(file,p);
    Write_Super_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
    create(file,out_file,filename);
    put_line(file,p);
    Write_Super_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
    create(file,out_file,filename);
    put_line(file,p);
    Write_Super_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
    create(file,out_file,filename);
    put_line(file,p);
    Write_Super_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

  procedure Write_Witness_Superset
              ( name : in string;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32 ) is

    filename : constant string := Append_ck(name,k);
    file : file_type;

  begin
    create(file,out_file,filename);
    put_line(file,p);
    Write_Super_Witness_Points(file,sols);
    close(file);
  end Write_Witness_Superset;

end Cascade_Homotopies_io;
