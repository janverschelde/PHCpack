with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Lists_of_Integer_Vectors;
with Lists_of_Floating_Vectors;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Permutations,Permute_Operations;    use Permutations,Permute_Operations;

package body Symmetric_Lifting_Functions is

  procedure Classify_Orbits
              ( supports : in Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                mix : in Standard_Integer_Vectors.Vector;
                v,w : in List_of_Permutations; norb : out natural32;
                orbits : out Arrays_of_Integer_Vector_Lists.Array_of_Lists ) is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    use Arrays_of_Integer_Vector_Lists;

    res,done,res_last,done_last : Array_of_Lists(mix'range);
    k,inmixk,orbit,cnt : integer32;
    n : constant integer32 := supports'last - supports'first + 1;
    tmp : List;
    lv : Link_to_Vector;

    function Lift ( lv : Link_to_Vector ) return Link_to_Vector is

      res : Link_to_Vector;

    begin
      res := new Vector(lv'first..lv'last+1);
      res(lv'range) := lv.all;
      res(res'last) := orbit;
      return res;
    end Lift;

    function Search_and_Lift ( lv : Link_to_Vector; l : List )
                             return Link_to_Vector is
    begin
      if Is_In(l,lv)
       then return Lift(lv);
       else return lv;
      end if;
    end Search_and_Lift;

    procedure Update ( k : in integer32; lv,liftlv : in Link_to_Vector ) is
    begin
      if not Is_In(done(k),lv) then
        Append(done(k),done_last(k),lv.all);
        Append(res(k),res_last(k),liftlv.all);
      end if;
    end Update;

  begin
    orbit := 0;
    k := supports'first;
    for i in mix'range loop
      tmp := supports(k);
      inmixk := Compute_Index(k,mix);
      while not Is_Null(tmp) loop
        lv := Head_Of(tmp);
        if not Is_In(done(inmixk),lv) then
          orbit := orbit + 1; -- new orbit
          declare
            tmpv,tmpw : List_of_Permutations;
            liftlv : Link_to_Vector := Lift(lv);
          begin
            Update(inmixk,lv,liftlv); Clear(liftlv);
            tmpv := v; tmpw := w;
            while not Is_Null(tmpv) loop   -- construct the orbit
              declare
                plv : Link_to_Vector := new Vector(lv'range);
                index : constant integer32 := Head_Of(tmpw)(k);
                inmix : integer32 := Compute_Index(index,mix);
              begin
                plv.all := Permutation(Head_Of(tmpv).all)*lv.all;
                liftlv := Search_and_Lift(plv,supports(index));
                if liftlv'last = n+1 then
                  inmix := Compute_Index(index,mix);
                  Update(inmix,plv,liftlv); Clear(liftlv);
                end if;
                Clear(plv);
              end;
              tmpv := Tail_Of(tmpv);
              tmpw := Tail_Of(tmpw);
            end loop;
          end;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      k := k + mix(i);
    end loop;
    Deep_Clear(done);
    cnt := 1;
    for i in res'range loop
      for j in 1..mix(i) loop
        orbits(cnt) := res(i);
        cnt := cnt + 1;
      end loop;
    end loop;
    norb := natural32(orbit);
  end Classify_Orbits;

  procedure Float_Lift_Orbits
              ( orbits : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                lifting : in Standard_Floating_Vectors.Vector ) is

    use Standard_Floating_Vectors;
    use Lists_of_Floating_Vectors;
    tmp : List;

  begin
    for k in orbits'range loop
      tmp := orbits(k);
      while not Is_Null(tmp) loop
        declare
          lv : Link_to_Vector := Head_Of(tmp);
        begin
          lv(lv'last) := lifting(integer32(lv(lv'last)));
          Set_Head(tmp,lv);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Float_Lift_Orbits;

  procedure Integer_Lift_Orbits
              ( orbits : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                lifting : in Standard_Integer_Vectors.Vector ) is

    use Standard_Integer_Vectors;
    use Lists_of_Integer_Vectors;
    tmp : List;

  begin
    for k in orbits'range loop
      tmp := orbits(k);
      while not Is_Null(tmp) loop
        declare
          lv : Link_to_Vector := Head_Of(tmp);
        begin
          lv(lv'last) := lifting(lv(lv'last));
          Set_Head(tmp,lv);
        end;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Integer_Lift_Orbits;

  procedure Float_Random_Lift_Orbits
              ( orbits : in out Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                norb : in natural32; lower,upper : in double_float ) is

    use Standard_Floating_Vectors;
    rv : Vector(1..integer32(norb));

  begin
    for k in rv'range loop
      rv(k) := (Random + 1.0)*(upper - lower)/2.0;
    end loop;
    Float_Lift_Orbits(orbits,rv);
  end Float_Random_Lift_Orbits;

  procedure Integer_Random_Lift_Orbits
              ( orbits : in out Arrays_of_Integer_Vector_Lists.Array_of_Lists;
                norb : in natural32; lower,upper : in integer32 ) is

    use Standard_Integer_Vectors;
    rv : Vector(1..integer32(norb));

  begin
    for k in rv'range loop
      rv(k) := Random(lower,upper);
    end loop;
    Integer_Lift_Orbits(orbits,rv);
  end Integer_Random_Lift_Orbits;

end Symmetric_Lifting_Functions;
