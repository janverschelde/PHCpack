with unchecked_deallocation;

package body Generic_NesVecs is

-- DEEP COPY :

  procedure Copy ( v : in NesVec; w : in out NesVec ) is
  begin
    if v.n = 1
     then Copy(v.v,w.v);
     else Copy(v.w,w.w);
    end if;
  end Copy;
 
  procedure Copy ( v : in Link_to_NesVec; w : in out Link_to_NesVec ) is
  begin
    if v /= null then
      declare
        w_rep : NesVec(v.n,v.a,v.b);
      begin
        Copy(v.all,w_rep);
        w := new NesVec'(w_rep);
      end;
    end if;
  end Copy;

  procedure Copy ( v : in Array_of_NesVecs; w : in out Array_of_NesVecs ) is
  begin
    for i in v'range loop
      Copy(v(i),w(i));
    end loop;
  end Copy;

-- DESTRUCTORS :

  procedure Clear ( v : in out NesVec ) is
  begin
    if v.n > 1
     then Clear(v.w);
    end if;
  end Clear;

  procedure Clear ( v : in out Link_to_NesVec ) is

    procedure free is new unchecked_deallocation(NesVec,Link_to_NesVec);

  begin
    if v /= null
     then Clear(v.all); free(v);
    end if;
  end Clear;

  procedure Clear ( v : in out Array_of_NesVecs ) is
  begin
    for i in v'range loop
      Clear(v(i));
    end loop;
  end Clear;

end Generic_NesVecs;
