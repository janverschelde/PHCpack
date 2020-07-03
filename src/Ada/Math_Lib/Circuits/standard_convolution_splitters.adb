with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with Standard_Vector_Splitters;

package body Standard_Convolution_Splitters is

  function Split ( c : Standard_Speelpenning_Convolutions.Circuit )
                 return Standard_Coefficient_Convolutions.Circuit is

    res : Standard_Coefficient_Convolutions.Circuit(c.nbr,c.dim,c.dim1,c.dim2);

    use Standard_Complex_Vectors;

  begin
    Standard_Integer_VecVecs.Copy(c.xps,res.xps);
    Standard_Integer_VecVecs.Copy(c.idx,res.idx);
    Standard_Integer_VecVecs.Copy(c.fac,res.fac);
    Standard_Vector_Splitters.Split_Complex(c.cff,res.rcf,res.icf);
    if c.cst /= null
     then Standard_Vector_Splitters.Split_Complex(c.cst,res.rct,res.ict);
    end if;
    Standard_Vector_Splitters.Split_Complex(c.forward,res.rfwd,res.ifwd);
    Standard_Vector_Splitters.Split_Complex(c.backward,res.rbck,res.ibck);
    Standard_Vector_Splitters.Split_Complex(c.cross,res.rcrs,res.icrs);
    Standard_Vector_Splitters.Split_Complex(c.wrk,res.rwrk,res.iwrk);
    Standard_Vector_Splitters.Split_Complex(c.acc,res.racc,res.iacc);
    return res;
  end Split;

  function Split ( c : Standard_Speelpenning_Convolutions.Link_to_Circuit )
                 return Standard_Coefficient_Convolutions.Link_to_Circuit is

    res : Standard_Coefficient_Convolutions.Link_to_Circuit;

    use Standard_Speelpenning_Convolutions;

  begin
    if c /= null then
      res := new Standard_Coefficient_Convolutions.Circuit'(Split(c.all));
    end if;
    return res;
  end Split;

  function Split ( c : Standard_Speelpenning_Convolutions.Circuits )
                 return Standard_Coefficient_Convolutions.Circuits is

    res : Standard_Coefficient_Convolutions.Circuits(c'range);

  begin
    for k in c'range loop
      res(k) := Split(c(k));
    end loop;
    return res;
  end Split;

  procedure Split
              ( p : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                rp : out Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                ip : out Standard_Floating_VecVecVecs.Link_to_VecVecVec ) is

    use Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

  begin
    if p /= null then
      declare
        rpwt,ipwt : Standard_Floating_VecVecVecs.VecVecVec(p'range);
      begin
        for k in p'range loop
          if p(k) /= null
           then Standard_Vector_Splitters.Split_Complex(p(k),rpwt(k),ipwt(k));
          end if;
        end loop;
        rp := new Standard_Floating_VecVecVecs.VecVecVec'(rpwt);
        ip := new Standard_Floating_VecVecVecs.VecVecVec'(ipwt);
      end;
    end if;
  end Split;

  function Split ( s : Standard_Speelpenning_Convolutions.System )
                 return Standard_Coefficient_Convolutions.System is

    res : Standard_Coefficient_Convolutions.System
            (s.neq,s.neq1,s.dim,s.dim1,s.deg);

  begin
    res.crc := Split(s.crc);
    Standard_Integer_Vectors.Copy(s.mxe,res.mxe);
    Split(s.pwt,res.rpwt,res.ipwt);
    Standard_Vector_Splitters.Split_Complex(s.yd,res.ryd,res.iyd);
    Standard_Complex_VecVecs.Copy(s.vy,res.vy);
    Standard_Complex_VecVecs.Copy(s.yv,res.yv);
    Standard_Complex_VecMats.Copy(s.vm,res.vm);
    return res;
  end Split;

  function Split ( s : Standard_Speelpenning_Convolutions.Link_to_System )
                 return Standard_Coefficient_Convolutions.Link_to_System is

    res : Standard_Coefficient_Convolutions.Link_to_System;

    use Standard_Speelpenning_Convolutions;

  begin
    if s /= null then
      res := new Standard_Coefficient_Convolutions.System'(Split(s.all));
    end if;
    return res;
  end Split;

end Standard_Convolution_Splitters;
