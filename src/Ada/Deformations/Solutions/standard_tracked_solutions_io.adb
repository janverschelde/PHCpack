with File_Scanning;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;

package body Standard_Tracked_Solutions_io is

  procedure get ( file : in file_type; 
                  lp,lq : out Link_to_Poly_Sys;
                  psols,qsols : out Solution_List;
                  verbose : in boolean := false ) is

    found : boolean;

  begin
    get(file,lp);
    if verbose then
      new_line;
      put_line("The target system :");
      put(lp.all);
    end if;
    File_Scanning.Scan_and_Skip(file,"START SYSTEM",found);
    if found then
      get(file,lq);
      if verbose then
        new_line;
        put_line("The start system :");
        put(lq.all);
      end if;
      File_Scanning.Scan_and_Skip(file,"START SOLUTIONS",found);
      if found then
        get(file,qsols);
        if verbose then
          new_line;
          put("Read ");
          put(Standard_Complex_Solutions.Length_Of(qsols),1);
          put_line(" start solutions.");
        end if;
        File_Scanning.Scan_and_Skip(file,"SOLUTIONS",found);
        if found then
          get(file,psols);
          if verbose then
            new_line;
            put("Read ");
            put(Standard_Complex_Solutions.Length_Of(psols),1);
            put_line(" solutions.");
          end if;
        end if;
      end if;
    end if;
  end get;

  procedure get ( file : in file_type; 
                  lp,lq : out Link_to_Laur_Sys;
                  psols,qsols : out Solution_List;
                  verbose : in boolean := false ) is

    found : boolean;

  begin
    get(file,lp);
    if verbose then
      new_line;
      put_line("The target system :");
      put(lp.all);
    end if;
    File_Scanning.Scan_and_Skip(file,"START SYSTEM",found);
    if found then
      get(file,lq);
      if verbose then
        new_line;
        put_line("The start system :");
        put(lq.all);
      end if;
      File_Scanning.Scan_and_Skip(file,"START SOLUTIONS",found);
      if found then
        get(file,qsols);
        if verbose then
          new_line;
          put("Read ");
          put(Standard_Complex_Solutions.Length_Of(qsols),1);
          put_line(" start solutions.");
        end if;
        File_Scanning.Scan_and_Skip(file,"SOLUTIONS",found);
        if found then
          get(file,psols);
          if verbose then
            new_line;
            put("Read ");
            put(Standard_Complex_Solutions.Length_Of(psols),1);
            put_line(" solutions.");
          end if;
        end if;
      end if;
    end if;
  end get;

end Standard_Tracked_Solutions_io;
