with String_Splitters;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Characters_and_Numbers;
with Partitions_of_Sets_Of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with Partitions_of_Sets_Strings;
with Set_Structure,Set_Structure_io;
with Set_Structure_Strings;

package body Root_Counters_Output is

  procedure Write_Root_Counts
               ( file : in file_type; no_mv : in boolean;
                 d : in natural64; mp_d : in Natural_Number;
                 m : in natural32; bz,bs : in natural64;
                 mv,smv : in natural32; z : in Partition ) is
  begin
    new_line(file);
    put_line(file,"ROOT COUNTS :");
    new_line(file);
    put(file,"total degree : ");
    if d = 0
     then put(file,mp_d,1); new_line(file);
     else put(file,d,1); new_line(file);
    end if;
    if m >= 1 then
      put(file,m,1); put(file,"-homogeneous Bezout number : ");
      put(file,bz,1); new_line(file);
      put(file,"  with partition : "); put(file,z(1..m)); new_line(file);
    end if;
    put(file,"general linear-product Bezout number : ");
    put(file,bs,1); new_line(file);
    if bs > 0 then
      put_line(file,"  based on the set structure :");
      Set_Structure_io.put(file);
    end if;
    if not no_mv then
      put(file,"mixed volume : "); put(file,mv,1); new_line(file);
      put(file,"stable mixed volume : "); put(file,smv,1); new_line(file);
    end if;
  end Write_Root_Counts;

  function Root_Counts_to_String
               ( no_mv : in boolean;
                 d : in natural64; mp_d : in Natural_Number;
                 m : in natural32; bz,bs : in natural64;
                 mv,smv : in natural32; z : in Partition ) return string is

    use String_Splitters;

    res : Link_to_String;
    nlc : character := ASCII.LF;
    snl : constant string(1..1) := (1..1 => nlc);

  begin
    Append(res,"total degree : ");
    Append(res,Characters_and_Numbers.nConvert(d));
    Append(res,snl);
    if m >= 1 then
      Append(res,Characters_and_Numbers.nConvert(m));
      Append(res,"-homogeneous Bezout number : ");
      Append(res,Characters_and_Numbers.nConvert(bz));
      Append(res,snl);
      Append(res,"  with with partition : ");
      Append(res,Partitions_of_Sets_Strings.to_string(z(1..m)));
      Append(res,snl);
    end if;
    Append(res,"general linear-product Bezout number : ");
    Append(res,Characters_and_Numbers.nConvert(bs));
    Append(res,snl);
    if bs > 0 then
      Append(res,"  based on the set structure :");
      Append(res,snl);
      for i in 1..Set_Structure.Dimension loop
        Append(res,"     ");
        Append(res,Set_Structure_Strings.to_String(i));
        Append(res,snl);
      end loop;
    end if;
    Append(res,"mixed volume : ");
    Append(res,Characters_and_Numbers.nConvert(mv));
    Append(res,snl);
    Append(res,"stable mixed volume : ");
    Append(res,Characters_and_Numbers.nConvert(smv));
    return res.all;
  end Root_Counts_to_String;

end Root_Counters_Output;
