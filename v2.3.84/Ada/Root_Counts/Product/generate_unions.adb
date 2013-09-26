procedure Generate_Unions ( k,first,last : in natural32 ) is

  acc : boolean_array(first..last) := (first..last => false);
  cont : boolean := true;

  procedure Generate ( k,start,last : in natural32;
                       acc : in out boolean_array ) is
  begin
    if k = 0 then
      process(acc,cont);
    elsif k > last - start + 1 then
      return;
    else
      for i in start..last loop
        acc(i) := true;
        generate(k-1,i+1,last,acc);
        exit when not cont;
        acc(i) := false;
      end loop;
    end if;
  end Generate;

begin
  Generate(k,first,last,acc);
end Generate_Unions;
