with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;

package body Crude_Path_Trackers is

  procedure Standard_Track_Paths ( verbose : in boolean := false ) is

    use Standard_Complex_Solutions;
 
    sols,tmp : Solution_List;
    len : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    crash : boolean;
    cnt : natural32 := 0;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    Standard_Solutions_Container.Clear;
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        PHCpack_Operations.Silent_Path_Tracker
          (ls,len,nbstep,nbfail,nbiter,nbsyst,crash);
        if verbose then
          cnt := cnt + 1;
          put("Solution "); put(cnt,1); put_line(" :");
          put_vector(ls.v);
        end if;
        Standard_Solutions_Container.Append(ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    if verbose then
      put("Number of solutions in the solution container : ");
      put(Standard_Solutions_Container.Length,1);
      new_line;
    end if;
  end Standard_Track_Paths;

  procedure DoblDobl_Track_Paths ( verbose : in boolean := false ) is

    use DoblDobl_Complex_Solutions;
 
    sols,tmp : Solution_List;
    len : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    crash : boolean;
    cnt : natural32 := 0;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    DoblDobl_Solutions_Container.Clear;
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        PHCpack_Operations.Silent_Path_Tracker
          (ls,len,nbstep,nbfail,nbiter,nbsyst,crash);
        if verbose then
          cnt := cnt + 1;
          put("Solution "); put(cnt,1); put_line(" :");
          put_vector(ls.v);
        end if;
        DoblDobl_Solutions_Container.Append(ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    if verbose then
      put("Number of solutions in the solution container : ");
      put(DoblDobl_Solutions_Container.Length,1);
      new_line;
    end if;
  end DoblDobl_Track_Paths;

  procedure QuadDobl_Track_Paths ( verbose : in boolean := false ) is

    use QuadDobl_Complex_Solutions;
 
    sols,tmp : Solution_List;
    len : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    crash : boolean;
    cnt : natural32 := 0;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    QuadDobl_Solutions_Container.Clear;
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
        ls : constant Link_to_Solution := Head_Of(tmp);
      begin
        PHCpack_Operations.Silent_Path_Tracker
          (ls,len,nbstep,nbfail,nbiter,nbsyst,crash);
        if verbose then
          cnt := cnt + 1;
          put("Solution "); put(cnt,1); put_line(" :");
          put_vector(ls.v);
        end if;
        QuadDobl_Solutions_Container.Append(ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    if verbose then
      put("Number of solutions in the solution container : ");
      put(QuadDobl_Solutions_Container.Length,1);
      new_line;
    end if;
  end QuadDobl_Track_Paths;

end Crude_Path_Trackers;
