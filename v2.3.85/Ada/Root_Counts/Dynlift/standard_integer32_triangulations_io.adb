with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs_io;        use Standard_Integer_VecVecs_io;
with Standard_Integer32_Simplices;       use Standard_Integer32_Simplices;
with Standard_Integer32_Simplices_io;    use Standard_Integer32_Simplices_io;

package body Standard_Integer32_Triangulations_io is

  procedure get ( t : in out Triangulation ) is
  begin
    get(Standard_Input,t);
  end get;

  procedure get ( n,m : in natural32; t : in out Triangulation ) is
  begin
    get(Standard_Input,n,m,t);
  end get;

  procedure get ( file : in file_type; t : in out Triangulation ) is

    n,m : natural32 := 0;
 
  begin
    get(file,n); get(file,m);
    get(file,n,m,t);
  end get;

  procedure get ( file : in file_type; n,m : in natural32;
                  t : in out Triangulation ) is
  begin
    for k in 1..m loop
      declare
        s : Simplex;
      begin
        get(file,n,s);  
        Construct(s,t);
      end;
    end loop;
    Connect(t);
  end get;

  procedure put ( n : in natural32; t : in Triangulation ) is
  begin
    put(Standard_Output,n,t);
  end put;

  procedure put ( n : in natural32; t : in Triangulation;
                  v : out natural32 ) is
  begin
    put(Standard_Output,n,t,v);
  end put;

  procedure put ( file : in file_type;
                  n : in natural32; t : in Triangulation ) is

    tmp : Triangulation := t;

  begin
    put(file,n,1); new_line(file);
    put(file,Length_Of(t),1); new_line(file);
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp));                    -- write the cell
      put(file,natural32(0),1); new_line(file);  -- write refinement of cell
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  n : in natural32; t : in Triangulation;
                  v : out natural32 ) is

    tmp : Triangulation := t;
    res,cnt : natural32 := 0;
    s : Simplex;
    vol : natural32;

  begin
    put(file,"Dimension without lifting : "); put(file,n,1); new_line(file);
    put(file,"Number of simplices : "); put(file,Length_Of(t),1);
    new_line(file);
    put_line(file,"The simplices in the triangulation :");
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      put(file,"Simplex "); put(file,cnt,1); put_line(file," :");
      s := Head_Of(tmp);
      put(file," with normal : "); put(file,Normal(s)); new_line(file);
      put_line(file," spanned by the points :"); put(file,Vertices(s));
      vol := Volume(s);
      put(file," ==> volume : "); put(file,vol,1); put_line(file,".");
      res := res + vol;
      tmp := Tail_Of(tmp);
    end loop;
    v := res;
  end put;

  function Position ( t : Triangulation; s : Simplex ) return natural32 is

  -- DESCRIPTION :
  --   Returns the position number of the simplex in the triangulation.
  --   Counting starts from one.  If the simplex s does not occur in t,
  --   then Length(t)+1 will be returned.

    res : natural32 := 1;
    tmp : Triangulation := t;
    s1 : Simplex;

  begin
    while not Is_Null(tmp) loop
      s1 := Head_Of(tmp);
      if Equal(s1,s)
       then return res;
       else tmp := Tail_Of(tmp);
            res := res + 1;
      end if;
    end loop;
    return res;
  end Position;

  function Connectivity ( t : Triangulation; s : Simplex ) return Vector is

  -- DESCRIPTION :
  --   Returns the connectivity vector of the given simplex w.r.t. the
  --   given triangulation.  This connectivity vector, name it cv, is
  --   defined as follows:
  --    cv(i) = 0 if Neighbor(s,i) = Null_Simplex
  --    cv(i) = k if Neighbor(s,i) /= Null_Simplex
  --         and Position(t,Neighbor(s,i)) = k.

    res : Vector(1..integer32(Dimension(s)));
    nei : Simplex;

  begin
    for i in res'range loop
      nei := Neighbor(s,i);
      if nei = Null_Simplex
       then res(i) := 0;
       else res(i) := integer32(Position(t,nei));
      end if;
    end loop;
    return res;
  end Connectivity;

  procedure put ( n : natural32; t : in Triangulation;
                  convecs : in out List; v : out natural32 ) is
  begin
    put(Standard_Output,n,t,convecs,v);
  end put;

  procedure put ( file : in file_type; n : natural32; t : in Triangulation;
                  convecs : in out List; v : out natural32 ) is

    tmp : Triangulation := t;
    s : Simplex;
    res,cnt : natural32 := 0;
    last : List;
    vol : natural32;

  begin
    put(file,"Dimension without lifting : "); put(file,n,1); new_line(file);
    put(file,"Number of simplices : "); put(file,Length_Of(t),1);
    new_line(file);
    put_line(file,"The simplices in the triangulation :");
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      put(file,"Simplex "); put(file,cnt,1); put_line(file," :");
      s := Head_Of(tmp);
      put(file," with normal : "); put(file,Normal(s)); new_line(file);
      put_line(file," spanned by the points :"); put(file,Vertices(s));
      declare
        cv : constant Vector := Connectivity(t,s);
        lcv : Link_to_Vector := new Vector'(cv);
      begin
        put(file," connectivity vector : "); put(file,cv); new_line(file);
        Append(convecs,last,lcv);
      end;
      vol := Volume(s);
      put(file," ==> volume : "); put(file,vol,1); put_line(file,".");
      res := res + vol;
      tmp := Tail_Of(tmp);
    end loop;
    v := res;
  end put;

end Standard_Integer32_Triangulations_io;
