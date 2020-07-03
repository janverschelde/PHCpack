with Standard_Floating_Vectors;

package body Standard_Coefficient_Storage is

  procedure Allocate_and_Store
              ( c : in Link_to_Circuit;
                rcf : out Standard_Floating_VecVecs.Link_to_VecVec;
                icf : out Standard_Floating_VecVecs.Link_to_VecVec ) is

    rpc,ipc : Standard_Floating_VecVecs.VecVec(0..c.nbr);
    lnk : Standard_Floating_Vectors.Link_to_Vector;

    use Standard_Floating_Vectors;

  begin
    if c.rct /= null
     then rpc(0) := new Standard_Floating_Vectors.Vector'(c.rct.all);
    end if;
    if c.ict /= null
     then ipc(0) := new Standard_Floating_Vectors.Vector'(c.ict.all);
    end if;
    for k in 1..c.nbr loop
      if c.rcf(k) /= null then
        lnk := c.rcf(k);
        declare
          vck : Standard_Floating_Vectors.Vector(lnk'range);
        begin
          for i in lnk'range loop
            vck(i) := lnk(i);
          end loop;
          rpc(k) := new Standard_Floating_Vectors.Vector'(vck);
        end;
      end if;
      if c.icf(k) /= null then
        lnk := c.icf(k);
        declare
          vck : Standard_Floating_Vectors.Vector(lnk'range);
        begin
          for i in lnk'range loop
            vck(i) := lnk(i);
          end loop;
          ipc(k) := new Standard_Floating_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    rcf := new Standard_Floating_VecVecs.VecVec'(rpc);
    icf := new Standard_Floating_VecVecs.VecVec'(ipc);
  end Allocate_and_Store;

  procedure Allocate_and_Store
              ( c : in Circuits;
                rcf : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                icf : out Standard_Floating_VecVecVecs.Link_to_VecVecVec ) is

    rcp,icp : Standard_Floating_VecVecVecs.VecVecVec(c'range);

  begin
    for k in c'range loop
      if c(k) /= null
       then Allocate_and_Store(c(k),rcp(k),icp(k));
      end if;
    end loop;
    rcf := new Standard_Floating_VecVecVecs.VecVecVec'(rcp);
    icf := new Standard_Floating_VecVecVecs.VecVecVec'(icp);
  end Allocate_and_Store;

  procedure Restore
              ( rcf : in Standard_Floating_VecVecs.Link_to_VecVec;
                icf : in Standard_Floating_VecVecs.Link_to_VecVec;
                c : in Link_to_Circuit ) is

    crclnk,storelnk : Standard_Floating_Vectors.Link_to_Vector;

    use Standard_Floating_Vectors;

  begin
    if rcf(0) /= null then
      storelnk := rcf(0);
      crclnk := c.rct;
      for k in storelnk'range loop
        crclnk(k) := storelnk(k);
      end loop;
    end if;
    if icf(0) /= null then
      storelnk := icf(0);
      crclnk := c.ict;
      for k in storelnk'range loop
        crclnk(k) := storelnk(k);
      end loop;
    end if;
    for k in 1..c.nbr loop
      if rcf(k) /= null then
        storelnk := rcf(k);
        crclnk := c.rcf(k);
        for i in storelnk'range loop
          crclnk(i) := storelnk(i);
        end loop;
      end if;
      if icf(k) /= null then
        storelnk := icf(k);
        crclnk := c.icf(k);
        for i in storelnk'range loop
          crclnk(i) := storelnk(i);
        end loop;
      end if;
    end loop;
  end Restore;

  procedure Restore
              ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                c : in Circuits ) is
  begin
    for k in c'range loop
      if c(k) /= null
       then Restore(rcf(k),icf(k),c(k));
      end if;
    end loop;
  end Restore;

end Standard_Coefficient_Storage;
