with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with DoblDobl_Vector_Splitters;

package body DoblDobl_Convolution_Splitters is

  function Split ( c : DoblDobl_Speelpenning_Convolutions.Circuit )
                 return DoblDobl_Coefficient_Convolutions.Circuit is

    res : DoblDobl_Coefficient_Convolutions.Circuit(c.nbr,c.dim,c.dim1,c.dim2);

    use DoblDobl_Complex_Vectors;

  begin
    Standard_Integer_VecVecs.Copy(c.xps,res.xps);
    Standard_Integer_VecVecs.Copy(c.idx,res.idx);
    Standard_Integer_VecVecs.Copy(c.fac,res.fac);
    DoblDobl_Vector_Splitters.Split_Complex
      (c.cff,res.rhcf,res.ihcf,res.rlcf,res.ilcf);
    if c.cst /= null then
      DoblDobl_Vector_Splitters.Split_Complex
        (c.cst,res.rhct,res.ihct,res.rlct,res.ilct);
    end if;
    DoblDobl_Vector_Splitters.Split_Complex
      (c.forward,res.rhfwd,res.ihfwd,res.rlfwd,res.ilfwd);
    DoblDobl_Vector_Splitters.Split_Complex
      (c.backward,res.rhbck,res.ihbck,res.rlbck,res.ilbck);
    DoblDobl_Vector_Splitters.Split_Complex
      (c.cross,res.rhcrs,res.ihcrs,res.rlcrs,res.ilcrs);
    DoblDobl_Vector_Splitters.Split_Complex
      (c.wrk,res.rhwrk,res.ihwrk,res.rlwrk,res.ilwrk);
    DoblDobl_Vector_Splitters.Split_Complex
      (c.acc,res.rhacc,res.ihacc,res.rlacc,res.ilacc);
    return res;
  end Split;

  function Split ( c : DoblDobl_Speelpenning_Convolutions.Link_to_Circuit )
                 return DoblDobl_Coefficient_Convolutions.Link_to_Circuit is

    res : DoblDobl_Coefficient_Convolutions.Link_to_Circuit;

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if c /= null then
      res := new DoblDobl_Coefficient_Convolutions.Circuit'(Split(c.all));
    end if;
    return res;
  end Split;

  function Split ( c : DoblDobl_Speelpenning_Convolutions.Circuits )
                 return DoblDobl_Coefficient_Convolutions.Circuits is

    res : DoblDobl_Coefficient_Convolutions.Circuits(c'range);

  begin
    for k in c'range loop
      res(k) := Split(c(k));
    end loop;
    return res;
  end Split;

  procedure Split
              ( p : in DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                rhp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ihp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rlp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ilp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec
              ) is

    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

  begin
    if p /= null then
      declare
        rhpwt : Standard_Floating_VecVecVecs.VecVecVec(p'range);
        ihpwt : Standard_Floating_VecVecVecs.VecVecVec(p'range);
        rlpwt : Standard_Floating_VecVecVecs.VecVecVec(p'range);
        ilpwt : Standard_Floating_VecVecVecs.VecVecVec(p'range);
      begin
        for k in p'range loop
          if p(k) /= null then
            DoblDobl_Vector_Splitters.Split_Complex
              (p(k),rhpwt(k),ihpwt(k),rlpwt(k),ilpwt(k));
          end if;
        end loop;
        rhp := new Standard_Floating_VecVecVecs.VecVecVec'(rhpwt);
        ihp := new Standard_Floating_VecVecVecs.VecVecVec'(ihpwt);
        rlp := new Standard_Floating_VecVecVecs.VecVecVec'(rlpwt);
        ilp := new Standard_Floating_VecVecVecs.VecVecVec'(ilpwt);
      end;
    end if;
  end Split;

  function Split ( s : DoblDobl_Speelpenning_Convolutions.System )
                 return DoblDobl_Coefficient_Convolutions.System is

    res : DoblDobl_Coefficient_Convolutions.System
            (s.neq,s.neq1,s.dim,s.dim1,s.deg);

  begin
    res.crc := Split(s.crc);
    Standard_Integer_Vectors.Copy(s.mxe,res.mxe);
    Split(s.pwt,res.rhpwt,res.ihpwt,res.rlpwt,res.ilpwt);
    DoblDobl_Vector_Splitters.Split_Complex
      (s.yd,res.rhyd,res.ihyd,res.rlyd,res.ilyd);
    DoblDobl_Complex_VecVecs.Copy(s.vy,res.vy);
    DoblDobl_Complex_VecVecs.Copy(s.yv,res.yv);
    DoblDobl_Complex_VecMats.Copy(s.vm,res.vm);
    return res;
  end Split;

  function Split ( s : DoblDobl_Speelpenning_Convolutions.Link_to_System )
                 return DoblDobl_Coefficient_Convolutions.Link_to_System is

    res : DoblDobl_Coefficient_Convolutions.Link_to_System;

    use DoblDobl_Speelpenning_Convolutions;

  begin
    if s /= null then
      res := new DoblDobl_Coefficient_Convolutions.System'(Split(s.all));
    end if;
    return res;
  end Split;

end DoblDobl_Convolution_Splitters;
