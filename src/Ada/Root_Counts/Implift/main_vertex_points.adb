with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Integer32_Vertices;        use Standard_Integer32_Vertices;

package body Main_Vertex_Points is

  procedure Vertex_Points ( file : in file_type; L : in out List ) is

    timer : Timing_Widget;
    res : List;

  begin
    put_line(file,"****  THE SUPPORT  ****");
    new_line(file); put(file,L); new_line(file);
    tstart(timer);
    res := Vertex_Points(L);
    tstop(timer);
    new_line(file);
    put_line(file,"****  THE VERTEX POINTS  ****");
    new_line(file); put(file,res); new_line(file);
    put_line(file,"****  REDUCTION OF POINTS  ****");
    new_line(file);
    put(file,"The number of points in the support   : ");
    put(file,Length_Of(L),1); new_line(file);
    put(file,"The number of remaining vertex points : ");
    put(file,Length_Of(res),1); new_line(file);
    new_line(file);
    print_times(file,timer,"computing the vertex points");
    new_line(file);
    Clear(L);
    L := res;
  end Vertex_Points;

  procedure Vertex_Points 
               ( file : in file_type; L : in out Array_of_Lists ) is

    timer : Timing_Widget;
    res : Array_of_Lists(l'range);

  begin
    tstart(timer);
    for i in res'range loop
      res(i) := Vertex_Points(L(i));
    end loop;
    tstop(timer);
    new_line(file);
    put_line(file,"****  THE VERTEX POINTS  ****");
    new_line(file);
    for i in res'range loop
      put(file,res(i)); new_line(file);
    end loop;
    put_line(file,"****  REDUCTION OF POINTS  ****");
    new_line(file);
    for i in L'range loop
      put(file,"The #points versus #vertices for support ");
      put(file,i,1); put(file," : ");
      put(file,Length_Of(L(i)),1); put(file,"  ->  ");
      put(file,Length_Of(res(i)),1); new_line(file);
    end loop;
    new_line(file);
    print_times(file,timer,"computing the vertex points");
    new_line(file);
    for i in l'range loop
      Clear(L(i));
      L(i) := res(i);
    end loop;
  end Vertex_Points;

  procedure Vertex_Points 
               ( file : in file_type; mix : in Link_to_Vector;
                 L : in out Array_of_Lists ) is

    timer : Timing_Widget;
    res : Array_of_Lists(L'range);
    cnt : integer32;

  begin
    tstart(timer);
    cnt := 1;
    for i in mix'range loop
      res(cnt) := Vertex_Points(L(cnt));
      cnt := cnt + mix(i);
    end loop;
    tstop(timer);
    new_line(file);
    put_line(file,"****  THE VERTEX POINTS  ****");
    new_line(file);
    cnt := 1;
    for i in mix'range loop
      put(file,res(cnt)); new_line(file);
      cnt := cnt + mix(i);
    end loop;
    put_line(file,"****  REDUCTION OF POINTS  ****");
    new_line(file);
    cnt := 1;
    for i in mix'range loop
      put(file,"The #points versus #vertices for support ");
      put(file,cnt,1); put(file," : ");
      put(file,Length_Of(L(cnt)),1); put(file,"  ->  ");
      put(file,Length_Of(res(cnt)),1); new_line(file);
      cnt := cnt + mix(i);
    end loop;
    new_line(file);
    print_times(file,timer,"computing the vertex points");
    new_line(file);
    cnt := 1;
    for i in mix'range loop
      for j in 0..(mix(i)-1) loop
        Clear(L(cnt+j));
        L(cnt+j) := res(cnt);
      end loop;
      cnt := cnt + mix(i);
    end loop;
  end Vertex_Points;

  procedure Vertex_Points
                ( file : in file_type; p : in out Poly_Sys ) is

    L : Array_of_Lists(p'range) := Create(p);
    res : Poly_Sys(p'range);

  begin
    Vertex_Points(file,l);
    res := Select_Terms(p,l);
    Deep_Clear(l);
    Clear(p); p := res;
  end Vertex_Points;

  procedure Vertex_Points
                ( file : in file_type; mix : in Link_to_Vector;
                  p : in out Poly_Sys ) is

    L : Array_of_Lists(p'range) := Create(p);
    res : Poly_Sys(p'range);

  begin
    Vertex_Points(file,mix,L);
    res := Select_Terms(p,l);
    Deep_Clear(L);
    Clear(p); p := res;
  end Vertex_Points;

  procedure Main is

    lp : Link_to_Poly_Sys;
    file : file_type;

  begin
    new_line;
    put_line("Extracting vertex points of support sets ...");
    new_line;
    get(lp);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Vertex_Points(file,lp.all);
  end Main;

end Main_Vertex_Points;
