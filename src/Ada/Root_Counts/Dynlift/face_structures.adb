package body Face_Structures is

-- CONSTRUCTOR :

  procedure Compute_New_Faces
                 ( fs : in out Face_Structure; k,n : in natural;
                   point : in out Link_to_Vector; newfs : out Faces ) is

    procedure Append ( fc : in VecVec ) is

      f : Face;

    begin
      f.points := Shallow_Create(fc);
      Append(res,res_last,f);
    end Append;
    procedure EnumLis is new Enumerate_Faces_in_List(Append);
    procedure EnumTri is new Enumerate_Faces_in_Triangulation(Append);

  begin
   -- COMPUTE THE NEW FACES AND UPDATE fs :
    if Is_Null(fs.t)
     then
       if Length_Of(fs.l) >= n
        then
          fs.t := Initial_Triangulation(n,fs.l,point);
          if Is_Null(fs.t)
           then EnumLis(fs.l,point,k);
           else EnumTri(fs.t,point,k);
          end if;
        else
          EnumLis(fs.l,point,k);
       end if;
     else
       declare
         newt : Triangulation;
       begin
         point(point'last) := Next_Lifting(fs.t,point.all);
         Update(fs.t,point,newt);
         Enumtri(newt,point,k);
       end;
    end if;
    Append(fs.l,fs.last,point);
    newfs := res;
  end Compute_New_Faces;

-- FLATTENING :

  procedure Flatten ( f : in out Face ) is
  begin
    if f.normal /= null
     then f.normal.all := (f.normal'range => 0);
          f.normal(f.normal'last) := 1;
    end if;
    Flatten(f.points);
  end Flatten;

  procedure Flatten ( fs : in out Faces ) is

    tmp : Faces := fs;
    f : Face;

  begin
    while not Is_Null(tmp) loop
      f := Head_Of(tmp);
      Flatten(f);
      Set_Head(tmp,f);
      tmp := Tail_Of(tmp);
    end loop;
  end Flatten;

  procedure Flatten ( fs : in out Face_Structure ) is
  begin
    Flatten(fs.l); Flatten(fs.t); Flatten(fs.f);
  end Flatten;

  procedure Flatten ( fs : in out Array_of_Face_Structures ) is
  begin
    for i in fs'range loop
      Flatten(fs(i));
    end loop;
  end Flatten;

-- SELECTORS :

  function Is_In ( fs : Face_Structure; point : vector ) return boolean is
  begin
    if Is_Null(fs.t)
     then return Is_In(fs.l,point);
     else return Is_In(fs.t,point);
    end if;
  end Is_In;

-- DESTRUCTORS :

  procedure Deep_Clear ( fs : in out Face_Structure ) is
  begin
    Deep_Clear(fs.l); Deep_Clear(fs.f); Deep_Clear(fs.t);
  end Deep_Clear;

  procedure Shallow_Clear ( fs : in out Face_Structure ) is
  begin
    Shallow_Clear(fs.l); Shallow_Clear(fs.f); Shallow_Clear(fs.t);
  end Shallow_Clear;

  procedure Deep_Clear ( fs : in out Array_of_Face_Structures ) is
  begin
    for i in fs'range loop
      Deep_Clear(fs(i));
    end loop;
  end Deep_Clear;

  procedure Shallow_Clear ( fs : in out Array_of_Face_Structures ) is
  begin
    for i in fs'range loop
      Shallow_Clear(fs(i));
    end loop;
  end Shallow_Clear;

end Face_Structures;
