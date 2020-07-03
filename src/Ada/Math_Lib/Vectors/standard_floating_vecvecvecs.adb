with unchecked_deallocation;

package body Standard_Floating_VecVecVecs is

  procedure Copy ( v_from : in Link_to_VecVecVec;
                   v_to : out Link_to_VecVecVec ) is

    use Standard_Floating_VecVecs;

  begin
    Clear(v_to);
    if v_from /= null then
      declare
        pwt : VecVecVec(v_from'range);
      begin
        for k in v_from'range loop
          if v_from(k) /= null then
            declare
              vv : Standard_Floating_VecVecs.VecVec(v_from(k)'range);
            begin
              Standard_Floating_VecVecs.Copy(v_from(k).all,vv);
              pwt(k) := new Standard_Floating_VecVecs.VecVec'(vv);
            end;
          end if;
        end loop;
        v_to := new VecVecVec'(pwt);
      end;
    end if;
  end Copy;


  procedure Clear ( v : in out VecVecVec ) is
  begin
    for k in v'range loop
      Standard_Floating_VecVecs.Deep_Clear(v(k));
    end loop;
  end Clear;

  procedure Clear ( v : in out Link_to_VecVecVec ) is

    procedure free is new unchecked_deallocation(VecVecVec,Link_to_VecVecVec);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

  procedure Clear ( v : in out VecVecVec_Array ) is
  begin
    for k in v'range loop
      Clear(v(k));
    end loop;
  end Clear;

end Standard_Floating_VecVecVecs;
