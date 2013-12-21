with unchecked_deallocation;

package body Sets_of_Unknowns is

-- REPRESENTATION OF A SET :

  type Set_Rep is array ( natural32 range <> ) of boolean;

  procedure free is new unchecked_deallocation(Set_Rep,Set);

-- CREATORS :

  function Create ( n : natural32 ) return Set is

    s : constant Set := new Set_Rep'(1..n => false);

  begin
    return s;
  end Create;

  function Create ( s : Set ) return Set is

    s1 : Set;

  begin
    if s = null
     then s1 := s;
     else s1 := new Set_Rep'(s.all);
    end if;
    return s1;
  end Create;

  function Universe ( n : natural32 ) return Set is

    s : Set := Create(n);

  begin
    for i in 1..n loop
      Add(s,i);
    end loop;
    return s;
  end Universe;
    
-- CONSTRUCTORS :

  procedure Add ( s : in out Set; i : in natural32 ) is
  begin
    s(i) := true;
  end Add;

  procedure Union ( s1 : in out Set; s2 : in Set ) is
  begin
    for i in 1..Dimension(s2) loop
      if Is_In(s2,i)
       then Add(s1,i);
      end if;
    end loop;
  end Union;

  function Union ( s1,s2 : Set ) return Set is

    s : Set := Create(s1);

  begin
    Union(s,s2);
    return s;
  end Union;

  procedure Remove ( s : in out Set; i : in natural32 ) is
  begin
    s(i) := false;
  end Remove;

  procedure Difference ( s1 : in out Set; s2 : in Set ) is
  begin
    for i in 1..Dimension(s2) loop
      if Is_In(s2,i)
       then Remove(s1,i);
      end if;
    end loop;
  end Difference;

  function Difference ( s1,s2 : Set ) return Set is

    s : Set := Create(s1);
 
  begin
    Difference(s,s2);
    return s;
  end Difference;

  procedure Intersection ( s1 : in out Set; s2 : in Set ) is
  begin
    for i in 1..Dimension(s1) loop
      if Is_In(s1,i) and then not Is_In(s2,i)
       then Remove(s1,i);
      end if;
    end loop;
  end Intersection;

  function Intersection ( s1,s2 : Set ) return Set is

    s : Set := Create(s1);

  begin
    Intersection(s,s2);
    return s;
  end Intersection;

-- SELECTORS :

  function Dimension ( s : Set ) return natural32 is
  begin
    if s = null
     then return 0;
     else return s'last;
    end if;
  end Dimension;

  function Extent_Of ( s : Set ) return natural32 is

    cnt : natural32 := 0;

  begin
    for i in 1..Dimension(s) loop
      if Is_In(s,i)
       then cnt := cnt + 1;
      end if;
    end loop;
    return cnt;
  end Extent_Of;

  function Is_In ( s : Set; i : natural32 ) return boolean is
  begin
    return s(i);
  end Is_In;

  function Is_Subset ( s1,s2 : Set ) return boolean is
  begin
    for i in 1..Dimension(s1) loop
      if Is_In(s1,i) and then not Is_In(s2,i)
       then return false;
      end if;
    end loop;
    return true;
  end Is_Subset;

  function Is_Equal ( s1,s2 : Set ) return boolean is
  begin
    return (Is_Subset(s1,s2) and then Is_Subset(s2,s1));
  end Is_Equal;

  procedure Generate_Subsets ( s : in Set; k : in natural32 ) is

    n : constant natural32 := Dimension(s);
    sub : Set := Create(n);
    cont : boolean;

    procedure Generate ( level,start : in natural32 ) is
    begin
      if level = 0 then
        Process(sub,cont);
      else
        for i in start..n-level+1 loop
          if Is_In(s,i) then
            Add(sub,i);
            Generate(level-1,i+1);
            Remove(sub,i);
          end if;
          exit when not cont;
        end loop;
      end if;
    end Generate;

  begin
    Generate(k,1);
    Clear(sub);
  end Generate_Subsets;
   
  procedure Generate_All_Subsets ( s : in Set ) is

    n : constant natural32 := Dimension(s);
    sub : Set := Create(n);
    cont : boolean;

    procedure Generate ( level,start : in natural32 ) is
    begin
      if level > 0 then
        for i in start..n loop
          if Is_In(s,i) then
            Add(sub,i);
            Process(sub,cont);
            if cont then
              Generate(level-1,i+1);
              Remove(sub,i);
            end if;
          end if;
          exit when not cont;
        end loop;
      end if;
    end Generate;

  begin
    Generate(n,1);
    Clear(sub);
  end Generate_All_Subsets;

-- DESTRUCTOR :

  procedure Clear ( s : in out Set ) is
  begin
    free(s);
  end Clear;

end Sets_of_Unknowns;
