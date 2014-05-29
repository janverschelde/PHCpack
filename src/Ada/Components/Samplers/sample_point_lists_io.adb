with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Sample_Points;                     use Sample_Points;
with Sample_Points_io;                  use Sample_Points_io;

package body Sample_Point_Lists_io is

  procedure get ( samples,samples_last : in out Standard_Sample_List ) is
  begin
    get(Standard_Input,samples,samples_last);
  end get;

  procedure get ( len,n,k : in natural32;
                  samples,samples_last : in out Standard_Sample_List ) is
  begin
    get(Standard_Input,len,n,k,samples,samples_last);
  end get;

  procedure get ( file : in file_type;
                  samples,samples_last : in out Standard_Sample_List ) is

    len,n,k : natural32 := 0;

  begin
    get(file,len); get(file,n); get(file,k);
    get(file,len,n,k,samples,samples_last);
  end get;

  procedure get ( file : in file_type; len,n,k : in natural32;
                  samples,samples_last : in out Standard_Sample_List ) is

    c : character;

  begin
    for i in 1..len loop
      get(file,c); skip_line(file);            -- skip the "Sample x :" line
      declare
        s : Standard_Sample;
      begin
        get(file,integer32(n),integer32(k),s);
        Append(samples,samples_last,s);
      end;
    end loop;
  end get;

  procedure get ( samples,samples_last : in out Multprec_Sample_List ) is
  begin
    get(Standard_Input,samples,samples_last);
  end get;

  procedure get ( len,n,k : in natural32;
                  samples,samples_last : in out Multprec_Sample_List ) is
  begin
    get(Standard_Input,len,n,k,samples,samples_last);
  end get;

  procedure get ( file : in file_type;
                  samples,samples_last : in out Multprec_Sample_List ) is

    len,n,k : natural32 := 0;

  begin
    get(file,len); get(file,n); get(file,k);
    get(file,len,n,k,samples,samples_last);
  end get;

  procedure get ( file : in file_type; len,n,k : in natural32;
                  samples,samples_last : in out Multprec_Sample_List ) is

    c : character;

  begin
    for i in 1..len loop
      get(file,c); skip_line(file);           -- skip the "Sample x :" line
      declare
        s : Multprec_Sample;
      begin
        get(file,integer32(n),integer32(k),s);
        Append(samples,samples_last,s);
      end;
    end loop;
  end get;

  procedure put ( samples : in Standard_Sample_List ) is
  begin
    put(Standard_Output,samples);
  end put;

  procedure put ( len,n,k : in natural32;
                  samples : in Standard_Sample_List ) is
  begin
    put(Standard_Output,len,n,k,samples);
  end put;

  procedure put ( file : in file_type; samples : in Standard_Sample_List ) is

    tmp : Standard_Sample_List := samples;

  begin
    for i in 1..Length_Of(samples) loop
      put(file,"Sample "); put(file,i,1); put_line(file," :");
      put(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( file : in file_type; len,n,k : in natural32;
                  samples : in Standard_Sample_List ) is
  begin
    put(file,len,1); put(file," ");
    put(file,n,1); put(file," ");
    put(file,k,1); new_line(file);
    put(file,samples);
  end put;

  procedure put ( samples : in Multprec_Sample_List ) is
  begin
    put(Standard_Output,samples);
  end put;

  procedure put ( len,n,k : in natural32;
                  samples : in Multprec_Sample_List ) is
  begin
    put(Standard_Output,len,n,k,samples);
  end put;

  procedure put ( file : in file_type; samples : in Multprec_Sample_List ) is

    tmp : Multprec_Sample_List := samples;

  begin
    for i in 1..Length_Of(samples) loop
      put(file,"Sample " ); put(file,i,1); put_line(file," : ");
      put(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( file : in file_type; len,n,k : in natural32;
                  samples : in Multprec_Sample_List ) is
  begin
    put(file,len,1); put(file," ");
    put(file,n,1); put(file," ");
    put(file,k,1); new_line(file);
    put(file,samples);
  end put;

  procedure Write_Summaries 
                ( file : in file_type; samples : in Standard_Sample_List ) is

    tmp : Standard_Sample_List := samples;

  begin
    for i in 1..Length_Of(samples) loop
      put(file,i,4); put(file," : ");
      declare
        spt : constant Standard_Sample := Head_Of(tmp);
      begin
        put(file," err : ");
        put(file,Sample_Point(spt).err,3);
        put(file,"   rco : ");
        put(file,Sample_Point(spt).rco,3);
        put(file,"   res : ");
        put(file,Sample_Point(spt).res,3);
      end;
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Summaries;

end Sample_Point_Lists_io;
