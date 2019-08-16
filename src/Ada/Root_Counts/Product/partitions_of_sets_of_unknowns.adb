with unchecked_deallocation;

package body Partitions_of_Sets_of_Unknowns is

-- CREATORS :

  procedure Create ( p : in out Partition; n : in natural32 ) is
  begin
    for i in p'range loop
      p(i) := Create(n);
    end loop;
  end Create;

  function Create ( p : Partition ) return Partition is

    res : Partition(p'range);

  begin
    for i in p'range loop
      res(i) := Create(p(i));
    end loop;
    return res;
  end Create;

-- CONSTRUCTOR :

  procedure Generate_Partitions ( s : in Set ) is

  -- NOTE :
  --   The algorithm below is a rather unelegant construction.
  --   The VADS compiler for IBM RS/6000 had problems with the nested
  --   generics, so the generation of all subsets is repeated here in full.

    n : constant natural32 := Dimension(s);
    continue : boolean := true;
    p : Partition(1..n);
    cnt : natural32 := 0;

    procedure Generate ( v : in Set; cont : out boolean );

    -- DESCRIPTION :
    --   Generation of all partitions makes use of a double recursive process.

    procedure Empty_Subsets ( w : in Set; cont : out boolean ) is
  
      rest : Set := Difference(w,p(cnt));

    begin
      if Extent_of(rest) = 0
       then Process(p(1..cnt),cont);
       else Generate(rest,cont);
      end if;
      Clear(rest);
    end Empty_Subsets;

    procedure All_Subsets ( w : in Set; cont : out boolean ) is

      sb : Set := Create(n);

      procedure Create_Partition ( sub : in Set; cont : out boolean ) is

        rest : Set;
        back : Set := Create(p(cnt));   -- back up copy needed to restore

      begin
        Union(p(cnt),sub);
        rest := Difference(w,p(cnt));
        if Extent_Of(rest) = 0
         then Process(p(1..cnt),cont);
         else Generate(rest,cont);
        end if;
        Clear(p(cnt)); p(cnt) := Create(back);
        Clear(rest); Clear(back);
      end Create_Partition;

      procedure Generate_Subset ( level,start : in natural32 ) is
      begin
        if level > 0 then
          for i in start..n loop
            if Is_In(w,i) then
              Add(sb,i);
              Create_Partition(sb,continue);
              if continue then
                Generate_Subset(level-1,i+1);
                Remove(sb,i);
              end if;
            end if;
            exit when not continue;
          end loop;
          cont := continue;
        end if;
      end Generate_Subset;

    begin
      Generate_Subset(n,1);
      Clear(sb);
    end All_Subsets;

    procedure Generate ( v : in Set; cont : out boolean ) is
    begin
      for i in 1..n loop
        if Is_In(v,i) then
          cnt := cnt + 1;
          p(cnt) := Create(n); Add(p(cnt),i); 
          Empty_Subsets(v,continue);
          if continue then
            declare
              w : Set := Create(v);
            begin
              Remove(w,i);
              All_Subsets(w,cont);
              Clear(w);
            end;
          end if;
          Clear(p(cnt)); cnt := cnt - 1;
          cont := continue;
        end if;
        exit when Is_In(v,i);
      end loop;
    end Generate;

  begin
    Generate(s,continue);
  end Generate_Partitions;

-- SELECTOR :

  function Number_of_Partitions ( n : natural32 ) return natural32 is

    sum : natural32;

    function comb ( n,i : natural32 ) return natural32 is

      n1,n2 : natural32 := 1;

    begin
      if (i = 0) or (i = n) then
        return 1;
      else
        for k in 1..i loop
          n1 := n1 * (n - k + 1);
          n2 := n2 * k;
        end loop;
        return (n1/n2);
      end if;
    end comb;

  begin
    if (n = 0) or (n = 1) then
      return 1;
    else
      sum := 0;
      for k in 0..(n-1) loop
        sum := sum + comb(n-1,k) * Number_Of_Partitions(n-1-k);
      end loop;
      return sum;
    end if;
  exception
    when CONSTRAINT_ERROR => return 0;
  end Number_of_Partitions;

-- DESTRUCTOR :

  procedure Clear ( p : in out Partition ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Partition ) is

    procedure free is
      new unchecked_deallocation(Partition,Link_to_Partition);

  begin
    if p /= null then
      Clear(p.all);
      free(p);
    end if;
  end Clear;

end Partitions_of_Sets_of_Unknowns;
