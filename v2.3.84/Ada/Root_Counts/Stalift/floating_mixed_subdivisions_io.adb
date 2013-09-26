with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors_io;
with Lists_of_Floating_Vectors_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;

package body Floating_Mixed_Subdivisions_io is

-- I/O of coordinate representations
-- INPUT ROUTINES :

  procedure get ( n,m : in natural32; mic : out Mixed_Cell ) is
  begin
    get(Standard_Input,n,m,mic);
  end get;

  procedure get ( file : in file_type;
	          n,m : in natural32; mic : out Mixed_Cell ) is

    adl : Array_of_Lists(1..integer32(m));
    L : natural32 := 0;

  begin
    Standard_Floating_Vectors_io.get(file,n+1,mic.nor);
    for k in 1..integer32(m) loop
      get(file,L);
      Lists_of_Floating_Vectors_io.get(file,n+1,L,adl(k));
    end loop;
    mic.pts := new Array_of_Lists'(adl);
    get(file,L);
    if l /= 0 then
      declare
        nn,mm : natural32 := 0;
        mix : Standard_Integer_Vectors.Link_to_Vector;
        sub : Mixed_Subdivision;
      begin
        get(file,nn,mm,mix,sub);
        if not Is_Null(sub)
         then mic.sub := new Mixed_Subdivision'(sub);
        end if;
      end;
    end if;
  end get;

  procedure get ( n,m : out natural32;
                  mixed_type : out Standard_Integer_Vectors.Link_to_Vector;
                  mixsub : out Mixed_Subdivision ) is
  begin
    get(Standard_Input,n,m,mixed_type,mixsub);
  end get;

  procedure get ( file : in file_type; n,m : out natural32;
		  mixed_type : out Standard_Integer_Vectors.Link_to_Vector;
                  mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    L,nn,mm : natural32 := 0;

  begin
    get(file,nn); n := nn;
    get(file,mm); m := mm;
    Standard_Integer_Vectors_io.get(file,mm,mixed_type);
    get(file,L);
    for k in 1..L loop
      declare
	mic : Mixed_Cell;
      begin
	get(file,nn,mm,mic);
	Append(res,res_last,mic);
      end;
    end loop;
    mixsub := res;
  end get;

-- OUTPUT ROUTINES :

  procedure put ( lifvec : in Standard_Floating_Vectors.Vector ) is
  begin
    put(Standard_Output,lifvec);
  end put;

  procedure put ( file : in file_type;
                  lifvec : in Standard_Floating_Vectors.Vector ) is
  begin
    for i in lifvec'first..lifvec'last-1 loop
      text_io.put(file,' '); put(file,lifvec(i),1,0,0);
    end loop;
    text_io.put(file,' '); put(file,lifvec(lifvec'last));
  end put;

  procedure put ( lifsup : in List ) is
  begin
    put(Standard_Output,lifsup);
  end put;

  procedure put ( file : in file_type; lifsup : in List ) is

    tmp : List := lifsup;

  begin
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp).all); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( lifsup : in Array_of_Lists ) is
  begin
    put(Standard_Output,lifsup);
  end put;

  procedure put ( file : in file_type; lifsup : in Array_of_Lists ) is
  begin
    for i in lifsup'range loop
      put(file,lifsup(i)); new_line(file);
    end loop;
  end put;

  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mic : in Mixed_Cell ) is
  begin
    put(Standard_Output,n,mix,mic);
  end put;

  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mic : in out Mixed_Cell; mv : out natural32;
                  multprec_hermite : in boolean := false ) is
  begin
    put(Standard_Output,n,mix,mic,mv,multprec_hermite);
  end put;

  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  mic : in Mixed_Cell ) is
  begin
    for i in mic.nor'range loop
      put(file,mic.nor(i)); new_line(file);
    end loop;
    for k in mic.pts'range loop
      put(file,Length_Of(mic.pts(k)),1); new_line(file);
      put(file,mic.pts(k));
    end loop;
    if mic.sub = null then
      put(file,natural32(0),1); new_line(file);
    else
      put(file,natural32(1),1); new_line(file);
      put(file,n,mix,mic.sub.all);
    end if;
  end put;
  
  procedure put ( file : in file_type;
                  n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mic : in out Mixed_Cell; mv : out natural32;
                  multprec_hermite : in boolean := false ) is
  begin
    text_io.put_line(file," normal to cell : ");
    for i in mic.nor'range loop
      put(file,mic.nor(i)); new_line(file);
    end loop;
    text_io.put_line(file," the points in the cell : ");
    for k in mic.pts'range loop
      text_io.put(file,"  component "); put(file,k,1); 
      text_io.put(file," with ");
      put(file,Length_Of(mic.pts(k)),1); text_io.put_line(file," points :");
      put(file,mic.pts(k));
    end loop;
    Mixed_Volume(integer32(n),mix,mic,mv,multprec_hermite);
    if mic.sub /= null
     then text_io.put_line(file," with refinement : ");
          put(file,n,mix,mic.sub.all,mv);
    end if;
  end put;

  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
		  mixsub : in Mixed_Subdivision ) is
  begin
    put(Standard_Output,n,mix,mixsub);
  end put;

  procedure put ( n : in natural32; mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural32;
                  multprec_hermite : in boolean := false ) is
  begin
    put(Standard_Output,n,mix,mixsub,mv,multprec_hermite);
  end put;

  procedure put ( file : in file_type; n : in natural32;
		  mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;

  begin
    put(file,n,1); new_line(file);
    put(file,mix'last,1); new_line(file);
    Standard_Integer_Vectors_io.put(file,mix); new_line(file);
    put(file,Length_Of(mixsub),1); new_line(file);
    while not Is_Null(tmp) loop
      put(file,n,mix,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural32;
                  multprec_hermite : in boolean := false ) is

    tmp : Mixed_Subdivision := mixsub;
    cnt,res : natural32 := 0;
    vol : natural32;

  begin
    text_io.put(file,"Dimension without lifting : ");
    put(file,n,1); new_line(file);
    text_io.put(file,"Number of different supports : ");
    put(file,mix'last,1); new_line(file);
    text_io.put(file,"Type of mixture : ");
    Standard_Integer_Vectors_io.put(file,mix); new_line(file);
    put_line(file,"The cells in the subdivision :");
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      text_io.put(file,"Cell "); put(file,cnt,1); text_io.put_line(file," :");
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        put(file,n,mix,mic,vol,multprec_hermite);
        Set_Head(tmp,mic);
      exception
        when others => vol := 0; Set_Head(tmp,mic);
                       put_line("patched exception in volume computation...");
      end;
      text_io.put(file,"==> Volume : "); put(file,vol,1); put_line(file,".");
      res := res + vol;
      tmp := Tail_Of(tmp);
    end loop;
    mv := res;
  end put;

  procedure put ( file : in file_type; n : in natural32; b : in double_float;
                  mix : in Standard_Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; 
                  mv,smv,tmv : out natural32;
                  multprec_hermite : in boolean := false ) is

    tmp : Mixed_Subdivision := mixsub;
    cnt : natural32 := 0;
    vol : natural32;

  begin
    mv := 0; smv := 0; tmv := 0;
    text_io.put(file,"Dimension without lifting : ");
    put(file,n,1); new_line(file);
    text_io.put(file,"Number of different supports : ");
    put(file,mix'last,1); new_line(file);
    text_io.put(file,"Type of mixture : ");
    Standard_Integer_Vectors_io.put(file,mix); new_line(file);
    put_line(file,"The cells in the subdivision :");
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      text_io.put(file,"Cell "); put(file,cnt,1); text_io.put_line(file," :");
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        put(file,n,mix,mic,vol,multprec_hermite);
        Set_Head(tmp,mic);
        text_io.put(file,"==> Volume : "); put(file,vol,1);
        if Is_Original(mic,b) then
          mv := mv + vol;
          smv := smv + vol;
          put_line(file," original.");
        elsif Is_Stable(mic.nor.all,b,mic.pts.all) then
          smv := smv + vol;
          put_line(file," stable.");
        else
          put_line(file," superfluous.");
        end if;
        tmv := tmv + vol;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end put;

-- INCREMENTAL READ and WRITE for coordinate representations :

  procedure Read_Dimensions
               ( file : in file_type; n,r,m : out natural32;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 fail : out boolean ) is
  begin
    n := 0; get(file,n);
    r := 0; get(file,r);
    mix := new Standard_Integer_Vectors.Vector(1..integer32(r));
    Standard_Integer_Vectors_io.get(file,mix.all);
    m := 0; get(file,m);
    fail := false;
  exception
    when others =>
      put_line("An exception occurred when reading the dimensions.");
      fail := true;
  end Read_Dimensions;

  procedure Write_Dimensions
               ( file : in file_type; n,r,m : in natural32;
                 mix : in Standard_Integer_Vectors.Vector ) is
  begin
    put(file,n,1); new_line(file);
    put(file,r,1); new_line(file);
    Standard_Integer_Vectors_io.put(file,mix); new_line(file);
    put(file,m,1); new_line(file);
  end Write_Dimensions;

  procedure Read_Next
               ( file : in file_type; n,r : in natural32;
                 mic : out Mixed_Cell; fail : out boolean ) is
  begin
    get(file,n,r,mic);
    fail := false;
  exception
    when others =>
      put_line("An exception occurred while reading a mixed cell.");
      fail := true;
  end Read_Next;

-- I/O of labeled representations :

  procedure get ( file : in file_type; n,r : in natural32;
                  mlb : out Mixed_Labels ) is

    m : natural32 := 0;

  begin
    Standard_Floating_Vectors_io.get(file,n+1,mlb.nor);
    mlb.lab := new Standard_Integer_VecVecs.VecVec(1..integer32(r));
    for i in 1..integer32(r) loop
      get(file,m);
      Standard_Integer_Vectors_io.get(file,m,mlb.lab(i));
    end loop;
    get(file,m);
    if m /= 0 then
      declare
        nn,rr : natural32;
        mix : Standard_Integer_Vectors.Link_to_Vector;
        sub : Mixed_Sublabeling;
      begin
        get(file,nn,rr,mix,sub);
        mlb.sub := new Mixed_Sublabeling'(sub);
      end;
    end if;
  end get;

  procedure get ( file : in file_type; n,r : out natural32;
                  mix : out Standard_Integer_Vectors.Link_to_Vector;
                  sub : out Mixed_Sublabeling ) is

    m : natural32 := 0;
    mlb : Mixed_Labels;
    last : Lists_of_Mixed_Labels.List;

  begin
    n := 0; r := 0;
    get(file,n);
    get(file,r);
    Standard_Integer_Vectors_io.get(file,r,mix);
    sub.pts := new Standard_Floating_VecVecs.Array_of_VecVecs(1..integer32(r));
    declare
      point : Standard_Floating_Vectors.Vector(1..integer32(n)+1);
    begin
      for i in 1..integer32(r) loop
        get(file,m);
        sub.pts(i) := new Standard_Floating_VecVecs.VecVec(1..integer32(m));
        for j in 1..integer32(m) loop
          Standard_Floating_Vectors_io.get(file,point);
          sub.pts(i)(j) := new Standard_Floating_Vectors.Vector'(point);
        end loop;
      end loop;
    end;
    get(file,m);
    for i in 1..m loop
      get(file,n,r,mlb);
      Lists_of_Mixed_Labels.Append(sub.cells,last,mlb);
    end loop;
  end get;

  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  mlb : in Mixed_Labels ) is
  begin
    Standard_Floating_Vectors_io.put_line(file,mlb.nor.all);
    for i in mlb.lab'range loop
      put(file,mlb.lab(i)'last,1); new_line(file);
      Standard_Integer_Vectors_io.put(file,mlb.lab(i).all);
      new_line(file);
    end loop;
    if mlb.sub = null
     then put_line(file,"0");
     else put_line(file,"1");
          put(file,n,mix,mlb.sub.all);
    end if;
  end put;

  procedure put ( file : in file_type; n : in natural32;
                  mix : in Standard_Integer_Vectors.Vector;
                  sub : in Mixed_Sublabeling ) is

    use Standard_Floating_VecVecs;
    use Lists_of_Mixed_Labels;
    tmp : Lists_of_Mixed_Labels.List := sub.cells;
    mlb : Mixed_Labels;

  begin
    put(file,n,1); new_line(file);
    put(file,mix'last,1); new_line(file);
    Standard_Integer_Vectors_io.put(file,mix); new_line(file);
    if sub.pts /= null then
      for i in sub.pts'range loop
        put(file,sub.pts(i)'last,1); new_line(file);
        for j in sub.pts(i)'range loop
          Floating_Mixed_Subdivisions_io.put(file,sub.pts(i)(j).all);
          new_line(file);
        end loop;
      end loop;
    end if;
    put(file,Length_Of(sub.cells),1); new_line(file);
    while not Is_Null(tmp) loop
      mlb := Head_Of(tmp);
      put(file,n,mix,mlb);
      tmp := Tail_Of(tmp);
    end loop;
  end put;

-- INCREMENTAL READ and WRITE for labeled representations :

  procedure Read_Dimensions
               ( file : in file_type; n,r : out natural32;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 fail : out boolean ) is
  begin
    n := 0; get(file,n);
    r := 0; get(file,r);
    mix := new Standard_Integer_Vectors.Vector(1..integer32(r));
    Standard_Integer_Vectors_io.get(file,mix.all);
    fail := false;
  exception
    when others =>
      put_line("An exception occurred when reading the dimensions.");
      fail := true;
  end Read_Dimensions;

  procedure Write_Dimensions
               ( file : in file_type; n,r : in natural32;
                 mix : in Standard_Integer_Vectors.Vector ) is
  begin
    put(file,n,1); new_line(file);
    put(file,r,1); new_line(file);
    Standard_Integer_Vectors_io.put(file,mix); new_line(file);
  end Write_Dimensions;

  procedure Read_Lifted_Supports
               ( file : in file_type; n,r : in natural32;
                 ls : out Standard_Floating_VecVecs.Array_of_VecVecs;
                 fail : out boolean ) is

    m : natural32 := 0;
    point : Standard_Floating_Vectors.Vector(1..integer32(n)+1);

  begin
    for i in ls'range loop
      get(file,m);
      ls(i) := new Standard_Floating_VecVecs.VecVec(1..integer32(m));
      for j in 1..integer32(m) loop
        Standard_Floating_Vectors_io.get(file,point);
        ls(i)(j) := new Standard_Floating_Vectors.Vector'(point);
      end loop;
    end loop;
    fail := false;
  exception 
    when others =>
      put_line("An exception occurred when reading the lifted point sets.");
      fail := true;
  end Read_Lifted_Supports;

  procedure Write_Lifted_Supports
               ( file : in file_type;
                 ls : in Standard_Floating_VecVecs.Array_of_VecVecs ) is
  begin
    for i in ls'range loop
      put(file,ls(i)'last,1); new_line(file);
      for j in ls(i)'range loop
        Floating_Mixed_Subdivisions_io.put(file,ls(i)(j).all);
        new_line(file);
      end loop;
    end loop;
  end Write_Lifted_Supports;

  procedure Read_Next
               ( file : in file_type; n,r,k : in natural32;
                 mlb : out Mixed_Labels; fail : out boolean ) is
  begin
    get(file,n,r,mlb);
    fail := false;
  exception
    when others => 
      put("An exception occurred when reading cell "); put(k,1);
      put_line(".");
      fail := true;
  end Read_Next;

end Floating_Mixed_Subdivisions_io;
