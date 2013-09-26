with Integer_Face_Enumerators;           use Integer_Face_Enumerators;

package body Enumerate_Faces_of_Polytope is

-- UTILITIES :

  function Is_In ( i : integer32; v : Vector ) return boolean is
  begin
    for j in v'range loop
      if v(j) = i
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

-- TARGET ROUTINES :

  procedure Enumerate_Faces_in_List
              ( L : in List; point : in Link_to_Vector; k : in integer32 ) is

    verpts : constant VecVec := Shallow_Create(L);
    face : VecVec(0..k);

    procedure Enumerate_Vertices ( ell,start : in integer32 ) is

    -- DESCRIPTION :
    --   Already l vectors have been choosen out of verpts and
    --   are in face(0..l).  In indices(0..l) the position of the
    --   ith face in the vector of vertices is given by indices(i).
    --   The parameter start indices from which vertex on, the
    --   enumeration has to start.
    --   In this way, the faces are enumerated in lexicographic order.

    begin
      if ell >= k then
        Process(face);
      elsif (verpts'last-start+1) >= k-ell then
        for j in start..verpts'last loop
          face(ell+1) := verpts(j);
          Enumerate_Vertices(ell+1,j+1);
        end loop;
      end if;
    end Enumerate_Vertices;

  begin
    face(0) := point;
    Enumerate_Vertices(0,verpts'first);
  end Enumerate_Faces_in_List;

  procedure Enumerate_Lower_Faces_in_List
              ( L : in List; point : in Link_to_Vector; k : in integer32 ) is

    verpts : VecVec(1..integer32(Length_Of(L))+1);
    face : VecVec(0..k);

  begin
    verpts(1) := point;
    verpts(2..verpts'last) := Shallow_Create(L);
    face(0) := point;
    if k = 1 then
      declare
        procedure Low_Edge ( i,j : in integer32; cont : out boolean ) is
        begin
          if i = 1 then
            face(1) := verpts(j);
            Process(face); cont := true;
          else
            cont := false;
          end if;
        end Low_Edge;
        procedure Enum_Low_Edges is new Enumerate_Lower_Edges(Low_Edge);
      begin
        Enum_Low_Edges(verpts);
      end;
    else
      declare
        procedure Low_Face ( f : in Vector; cont : out boolean ) is
        begin
          if f(f'first) = 1 then
            for i in 1..f'last loop
              face(i) := verpts(f(i));
            end loop;
            Process(face); cont := true;
          else
            cont := false;
          end if;
        end Low_Face;
        procedure Enum_Low_Faces is new Enumerate_Lower_Faces(Low_Face);
       begin
         Enum_Low_Faces(k,verpts);
       end;
    end if;
  end Enumerate_Lower_Faces_in_List;

  procedure Enumerate_Faces_in_Simplex
              ( s : in Simplex; point : in Link_to_Vector;
                k : in integer32 ) is

    verpts : constant VecVec := Vertices(s);
    face : VecVec(0..k);
    indices : Vector(0..k);

    procedure Enumerate_Vertices ( ell,start : in integer32 ) is

    -- DESCRIPTION :
    --   Already l vectors have been choosen out of verpts and
    --   are in face(0..l).  In indices(0..l) the position of the
    --   ith face in the vector of vertices is given by indices(i).
    --   The parameter start indices from which vertex on, the
    --   enumeration has to start.
    --   In this way, the faces are enumerated in lexicographic order.

    begin
      if ell >= k then
        Process(face);
      elsif (verpts'last-start+1) >= k-ell then
        for j in start..verpts'last loop
          if not Is_In(j,indices(0..ell)) then
            indices(ell+1) := j;
            face(ell+1) := verpts(j);
            Enumerate_Vertices(ell+1,j+1);
          end if;
        end loop;
      end if;
    end Enumerate_Vertices;

  begin
   -- SEARCH FOR THE INDEX OF THE POINT :
    indices(0) := verpts'first - 1;
    for i in verpts'range loop
      if verpts(i).all = point.all
       then indices(0) := i;
      end if;
      exit when indices(0) >= verpts'first;
    end loop;
   -- ENUMERATE ALL OTHER k CHOICES OF POINTS :
    if indices(0) >= verpts'first
     then 
       face(0) := point;
       Enumerate_Vertices(0,verpts'first);
    end if;
  end Enumerate_Faces_in_Simplex;

 -- procedure Enumerate_Lower_Faces_in_Simplex
 --               ( s : in Simplex; point : in Link_to_Vector;
--                  k : in natural ) is

 --   verpts : constant VecVec := Vertices(s);
 --   face : VecVec(0..k);
 --   indices : Vector(0..k);

 -- begin
 --   null;
 -- end Enumerate_Lower_Faces_in_Simplex;

  procedure Enumerate_Faces_in_Triangulation
                ( t : in Triangulation;
                  point : in Link_to_Vector; k : in integer32 ) is

    verpts : constant VecVec := Vertices(t,point.all);
    face : VecVec(0..k);
    indices : Vector(0..k);

    procedure Enumerate_Vertices ( ell,start : in integer32 ) is

    -- DESCRIPTION :
    --   Already l vectors have been choosen out of verpts and
    --   are in face(0..l).  In indices(0..l) the position of the
    --   ith face in the vector of vertices is given by indices(i).
    --   The parameter start indices from which vertex on, the
    --   enumeration has to start.
    --   In this way, the faces are enumerated in lexicographic order.

    begin
      if ell >= k then
        Process(face);
      elsif (verpts'last-start+1) >= k-ell then
        for j in start..verpts'last loop
          if not Is_In(j,indices(0..ell)) then
            indices(ell+1) := j;
            face(ell+1) := verpts(j);
            Enumerate_Vertices(ell+1,j+1);
          end if;
        end loop;
      end if;
    end Enumerate_Vertices;

  begin
   -- SEARCH FOR THE INDEX OF THE POINT :
    indices(0) := verpts'first - 1;
    for i in verpts'range loop
      if verpts(i).all = point.all
       then indices(0) := i;
      end if;
      exit when indices(0) >= verpts'first;
    end loop;
   -- ENUMERATE ALL OTHER k CHOICES OF POINTS :
    if indices(0) >= verpts'first
     then
       face(0) := point;
       Enumerate_Vertices(0,verpts'first);
    end if;
  end Enumerate_Faces_in_Triangulation;

--  procedure Enumerate_Lower_Faces_in_Triangulation
--                ( t : in Triangulation;
--                  point : in Link_to_Vector; k : in natural ) is
--
--    verpts : constant VecVec := Vertices(t,point.all);
--    face : VecVec(0..k);
--    indices : Vector(0..k);
--
--  begin
--    null;
--  end Enumerate_Lower_Faces_in_Triangulation;

end Enumerate_Faces_of_Polytope;
