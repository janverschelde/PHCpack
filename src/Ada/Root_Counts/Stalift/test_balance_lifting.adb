with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Floating_Integer_Convertors;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;
with Lifted_Configurations;              use Lifted_Configurations;

package body Test_Balance_Lifting is

  procedure Read_Subdivision
             ( n : out natural32;
               mix : out Standard_Integer_Vectors.Link_to_Vector;
               mixsub : out Mixed_Subdivision ) is

    insubft : file_type;
    m : natural32 := 0;

  begin
    new_line;
    put_line("Reading the name of the file where subdivision is.");
    Read_Name_and_Open_File(insubft);
    get(insubft,n,m,mix,mixsub);
    Close(insubft);
  exception
    when DATA_ERROR
       => put_line("Data not in correct format.  Will ignore it...");
          Close(insubft);
  end Read_Subdivision;

  function Inner_Product
             ( point : Standard_Floating_Vectors.Link_to_Vector;
               lifvals,normal : Standard_Floating_Vectors.Link_to_Vector )
             return double_float is

   res : double_float := 0.0;
   ind : constant integer32 := integer32(point(point'last));

  begin
    for i in point'first..point'last-1 loop
      res := res + point(i)*normal(i);
    end loop;
    res := res + lifvals(ind)*normal(normal'last);
    return res;
  end Inner_Product;

  procedure Extremal_Inner_Products
             ( points : in Lists_of_Floating_Vectors.List;
               lifvals,normal : in Standard_Floating_Vectors.Link_to_Vector;
               min,max : out double_float ) is

   use Lists_of_Floating_Vectors;

   tmp : List;
   lpt : Standard_Floating_Vectors.Link_to_Vector;
   ip : double_float;

  begin
    lpt := Head_Of(points);
    max := Inner_Product(lpt,lifvals,normal);
    min := max;
    tmp := Tail_Of(points);
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      ip := Inner_Product(lpt,lifvals,normal);
      if ip > max then
        max := ip;
      elsif ip < min then
        min := ip;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Extremal_Inner_Products;

  function Maximal_Inner_Product
             ( points : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               lifvals,normal : Standard_Floating_Vectors.Link_to_Vector )
             return double_float is

    res,min,max,ip_min,ip_max : double_float;

  begin
    Extremal_Inner_Products(points(points'first),lifvals,normal,min,max);
    for i in points'first+1..points'last loop
      Extremal_Inner_Products(points(i),lifvals,normal,ip_min,ip_max);
      if ip_max > max
       then max := ip_max;
      end if;
      if ip_min < min
       then min := ip_min;
      end if;
    end loop;
    res := max - min;
    return res;
  end Maximal_Inner_Product;

  procedure Balance_Lifting_Values
              ( file : in file_type;
                n : in integer32;
                mix : in Standard_Integer_Vectors.Vector;
                mcc : in out Mixed_Subdivision ) is

    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(mix'range);
    lif : Standard_Floating_Vectors.Link_to_Vector;
    nbp : natural32;
    tmp : Mixed_Subdivision;
    mic : Mixed_Cell;
    max,ipr : double_float;
    ind : integer32 := 1;

  begin
    Collect_Points_and_Lifting(n,mix,mcc,pts,nbp,lif);
    sup := Floating_Integer_Convertors.Convert(pts);
    put("Number of lifted points : "); put(nbp,1); put_line(".");
    put(file,"There are "); put(file,nbp,1);
    put_line(file," lifted points :"); put(file,sup);
    put_line(file,"The "); put(file,nbp,1);
    put_line(file," lifting values :"); put_line(file,lif);
    mic := Head_Of(mcc);
    max := Maximal_Inner_Product(pts,lif,mic.nor);
    put("maximal inner product at cell "); put(ind,1);
    put(" : "); put(max); new_line;
    tmp := Tail_Of(mcc);
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      ipr := Maximal_Inner_Product(pts,lif,mic.nor);
      ind := ind + 1;
      put("maximal inner product at cell "); put(ind,1);
      put(" : "); put(ipr); new_line;
      if ipr > max
       then max := ipr;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put("Maximal inner product over all cells : ");
    put(max); new_line;
  end Balance_Lifting_Values;

  procedure Main is

    n : integer32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    mixsub : Mixed_Subdivision;
    outfile : file_type;

  begin
    Read_Subdivision(natural32(n),mix,mixsub);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    Balance_Lifting_Values(outfile,n,mix.all,mixsub);
    new_line(outfile);
    Close(outfile);
  end Main;

end Test_Balance_Lifting;
