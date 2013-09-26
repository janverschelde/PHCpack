with File_Scanning;                      use File_Scanning;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;

package body Irreducible_Decompositions_io is

  procedure get ( dc : in out Standard_Irreducible_Decomposition ) is
  begin
    get(Standard_Input,dc);
  end get;

  procedure get ( file : in file_type;
                  dc : in out Standard_Irreducible_Decomposition ) is

    nb : natural32 := 0;
    k : integer32 := 0;
    found : boolean;

  begin
    get(file,k);
    dc := Create(k);
    for i in 0..k loop
      declare
        lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
        sols : Standard_Complex_Solutions.Solution_List;
        len,n : natural32 := 0;
      begin
        get(file,lp);
        Add_Embedding(dc,i,lp);
        Scan_and_Skip(file,"SOLUTIONS",found);
        exit when not found;
        get(file,len);
        if len = 0 then
          skip_line(file);
        else
          get(file,n); skip_line(file);
          get(file,len,n,sols);
          declare
            cpsols : Standard_Complex_Solutions.Solution_List;
          begin
            Standard_Complex_Solutions.Copy(sols,cpsols);
            Add_Generic_Points(dc,i,cpsols);
          end;
        end if;
        if i < k then
          nb := Symbol_Table.Number;
          if nb < natural32(lp'last)+1
           then Symbol_Table.Enlarge(natural32(lp'last)+1-nb);
          end if;
        end if;
      end;
    end loop;
  end get;

  procedure put ( dc : in Standard_Irreducible_Decomposition ) is
  begin
    put(Standard_Output,dc);
  end put;

  procedure put ( file : in file_type;
                  dc : in Standard_Irreducible_Decomposition ) is

    k : constant integer32 := Top_Dimension(dc);

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

  begin
    put(file,k,1); new_line(file);
    for i in 0..k loop
      declare
        lp : constant Link_to_Poly_Sys := Embedding(dc,i);
        sols : constant Solution_List := Generic_Points(dc,i);
      begin
        if lp = null then
          put(file,"0");
        else
          put(file,lp'last,1); new_line(file);
          put(file,lp.all);
        end if;
        new_line(file);
        put_line(file,"THE SOLUTIONS :");
        if Is_Null(sols)
         then put_line(file,"0 0");
         else put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
        end if;
      end;
      new_line(file);
    end loop;
  end put;

end Irreducible_Decompositions_io;
