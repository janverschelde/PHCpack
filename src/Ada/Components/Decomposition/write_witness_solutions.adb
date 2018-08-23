with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Path_Counts_Table;
with Standard_Witness_Solutions;
with DoblDobl_Witness_Solutions;
with QuadDobl_Witness_Solutions;

package body Write_Witness_Solutions is

  procedure Standard_Write ( topdim,lowdim : in natural32 ) is

    use Standard_Complex_Solutions;

    sols : Solution_List;
    len : natural32;

  begin
    for k in lowdim..topdim loop
      sols := Standard_Witness_Solutions.Load_Witness_Points(k);
      len := Length_Of(sols);
      put("Number of solutions at dimension "); put(k,1);
      put(" : "); put(len,1); new_line;
    end loop;
  end Standard_Write;

  procedure DoblDobl_Write ( topdim,lowdim : in natural32 ) is

    use DoblDobl_Complex_Solutions;

    sols : Solution_List;
    len : natural32;

  begin
    for k in lowdim..topdim loop
      sols := DoblDobl_Witness_Solutions.Load_Witness_Points(k);
      len := Length_Of(sols);
      put("Number of solutions at dimension "); put(k,1);
      put(" : "); put(len,1); new_line;
    end loop;
  end DoblDobl_Write;

  procedure QuadDobl_Write ( topdim,lowdim : in natural32 ) is

    use QuadDobl_Complex_Solutions;

    sols : Solution_List;
    len : natural32;

  begin
    for k in lowdim..topdim loop
      sols := QuadDobl_Witness_Solutions.Load_Witness_Points(k);
      len := Length_Of(sols);
      put("Number of solutions at dimension "); put(k,1);
      put(" : "); put(len,1); new_line;
    end loop;
  end QuadDobl_Write;

  procedure Write_Counts
              ( filter,factor : in boolean;
                pc,fc : in Standard_Natural_VecVecs.Link_to_VecVec;
                idx : in Standard_Natural_VecVecs.Link_to_Array_of_VecVecs ) is
  begin
    Path_Counts_Table.Write_Path_Counts(standard_output,pc.all);
    if filter
     then Path_Counts_Table.Write_Filter_Counts(standard_output,fc.all);
    end if;
    if factor
     then Path_Counts_Table.Write_Decomposition(standard_output,idx.all);
    end if;
  end Write_Counts;

end Write_Witness_Solutions;
