with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with DoblDobl_Complex_Linear_Solvers;    use DoblDobl_Complex_Linear_Solvers;
with QuadDobl_Complex_Linear_Solvers;    use QuadDobl_Complex_Linear_Solvers;
with Standard_Mixed_Residuals;
with DoblDobl_Mixed_Residuals;
with QuadDobl_Mixed_Residuals;

package body Corrector_Convolutions is

-- MANAGEMENT OF LEADING COEFFICIENTS :

  procedure Allocate_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                lead : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    cff : Standard_Complex_VecVecs.VecVec(c'range);

    use Standard_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null then
        declare
          vck : Standard_Complex_Vectors.Vector(0..c(k).nbr);
        begin
          vck := (vck'range => Standard_Complex_Numbers.Create(0.0));
          cff(k) := new Standard_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    lead := new Standard_Complex_VecVecs.VecVec'(cff);
  end Allocate_Leading_Coefficients;

  procedure Allocate_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                lead : out DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    cff : DoblDobl_Complex_VecVecs.VecVec(c'range);

    use DoblDobl_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null then
        declare
          vck : DoblDobl_Complex_Vectors.Vector(0..c(k).nbr);
        begin
          vck := (vck'range => DoblDobl_Complex_Numbers.Create(integer(0)));
          cff(k) := new DoblDobl_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    lead := new DoblDobl_Complex_VecVecs.VecVec'(cff);
  end Allocate_Leading_Coefficients;

  procedure Allocate_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                lead : out QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    cff : QuadDobl_Complex_VecVecs.VecVec(c'range);

    use QuadDobl_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null then
        declare
          vck : QuadDobl_Complex_Vectors.Vector(0..c(k).nbr);
        begin
          vck := (vck'range => QuadDobl_Complex_Numbers.Create(integer(0)));
          cff(k) := new QuadDobl_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    lead := new QuadDobl_Complex_VecVecs.VecVec'(cff);
  end Allocate_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in Standard_Complex_Vectors.Link_to_Vector ) is

    use Standard_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst = null
     then lead(0) := Standard_Complex_Numbers.Create(0.0);
     else lead(0) := c.cst(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lead(k) := lnk(0);
    end loop;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in DoblDobl_Complex_Vectors.Link_to_Vector ) is

    use DoblDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst = null
     then lead(0) := DoblDobl_Complex_Numbers.Create(integer(0));
     else lead(0) := c.cst(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lead(k) := lnk(0);
    end loop;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                lead : in QuadDobl_Complex_Vectors.Link_to_Vector ) is

    use QuadDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst = null
     then lead(0) := QuadDobl_Complex_Numbers.Create(integer(0));
     else lead(0) := c.cst(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lead(k) := lnk(0);
    end loop;
  end Store_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in Standard_Complex_Vectors.Link_to_Vector;
                c : in Standard_Speelpenning_Convolutions.Link_to_Circuit ) is

    use Standard_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst /= null
     then c.cst(0) := lead(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lnk(0) := lead(k);
    end loop;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in DoblDobl_Complex_Vectors.Link_to_Vector;
                c : in  DoblDobl_Speelpenning_Convolutions.LInk_to_Circuit ) is

    use DoblDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst /= null
     then c.cst(0) := lead(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lnk(0) := lead(k);
    end loop;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in QuadDobl_Complex_Vectors.Link_to_Vector;
                c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit ) is

    use QuadDobl_Complex_Vectors;
    lnk : Link_to_Vector;

  begin
    if c.cst /= null
     then c.cst(0) := lead(0);
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      lnk(0) := lead(k);
    end loop;
  end Restore_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                lead : in Standard_Complex_VecVecs.Link_to_VecVec ) is

    use Standard_Complex_Vectors,Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant Standard_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Store_Leading_Coefficients(c(k),vck);
          end;
        end if;
      end loop;
    end if;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                lead : in DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    use DoblDobl_Complex_Vectors,DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant DoblDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Store_Leading_Coefficients(c(k),vck);
          end;
        end if;
      end loop;
    end if;
  end Store_Leading_Coefficients;

  procedure Store_Leading_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                lead : in QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    use QuadDobl_Complex_Vectors,QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant QuadDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Store_Leading_Coefficients(c(k),vck);
          end;
        end if;
      end loop;
    end if;
  end Store_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in Standard_Complex_VecVecs.Link_to_VecVec;
                c : in Standard_Speelpenning_Convolutions.Circuits ) is

    use Standard_Complex_Vectors,Standard_Complex_VecVecs;
    use Standard_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant Standard_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Restore_Leading_Coefficients(vck,c(k));
          end;
        end if;
      end loop;
    end if;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits ) is

    use DoblDobl_Complex_Vectors,DoblDobl_Complex_VecVecs;
    use DoblDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant DoblDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Restore_Leading_Coefficients(vck,c(k));
          end;
        end if;
      end loop;
    end if;
  end Restore_Leading_Coefficients;

  procedure Restore_Leading_Coefficients
              ( lead : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits ) is

    use QuadDobl_Complex_Vectors,QuadDobl_Complex_VecVecs;
    use QuadDobl_Speelpenning_Convolutions;

  begin
    if lead /= null then
      for k in c'range loop
        if c(k) /= null and lead(k) /= null then
          declare
            vck : constant QuadDobl_Complex_Vectors.Link_to_Vector := lead(k);
          begin
            Restore_Leading_Coefficients(vck,c(k));
          end;
        end if;
      end loop;
    end if;
  end Restore_Leading_Coefficients;

-- MANAGEMENT OF ALL COEFFICIENTS :

  procedure Allocate_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec ) is

    res : Standard_Complex_VecVecs.VecVec(0..c.nbr);

    use Standard_Complex_Vectors;

  begin
    if c.cst /= null
     then res(0) := new Standard_Complex_Vectors.Vector'(c.cst.all);
    end if;
    for k in c.cff'range loop
      if c.cff(k) /= null then
        declare
          vck : Standard_Complex_Vectors.Vector(c.cff(k)'range);
        begin
          vck := (vck'range => Standard_Complex_Numbers.Create(0.0));
          res(k) := new Standard_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    cff := new Standard_Complex_VecVecs.VecVec'(res);
  end Allocate_Coefficients;

  procedure Allocate_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : out DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    res : DoblDobl_Complex_VecVecs.VecVec(0..c.nbr);

    use DoblDobl_Complex_Vectors;

  begin
    if c.cst /= null
     then res(0) := new DoblDobl_Complex_Vectors.Vector'(c.cst.all);
    end if;
    for k in c.cff'range loop
      if c.cff(k) /= null then
        declare
          vck : DoblDobl_Complex_Vectors.Vector(c.cff(k)'range);
        begin
          vck := (vck'range => DoblDobl_Complex_Numbers.Create(integer(0)));
          res(k) := new DoblDobl_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    cff := new DoblDobl_Complex_VecVecs.VecVec'(res);
  end Allocate_Coefficients;

  procedure Allocate_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : out QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    res : QuadDobl_Complex_VecVecs.VecVec(0..c.nbr);

    use QuadDobl_Complex_Vectors;

  begin
    if c.cst /= null
     then res(0) := new QuadDobl_Complex_Vectors.Vector'(c.cst.all);
    end if;
    for k in c.cff'range loop
      if c.cff(k) /= null then
        declare
          vck : QuadDobl_Complex_Vectors.Vector(c.cff(k)'range);
        begin
          vck := (vck'range => QuadDobl_Complex_Numbers.Create(integer(0)));
          res(k) := new QuadDobl_Complex_Vectors.Vector'(vck);
        end;
      end if;
    end loop;
    cff := new QuadDobl_Complex_VecVecs.VecVec'(res);
  end Allocate_Coefficients;

  procedure Allocate_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                cff : out Standard_Speelpenning_Convolutions.
                          Link_to_VecVecVec ) is

    use Standard_Speelpenning_Convolutions;

    res : VecVecVec(c'range);

  begin
    for k in c'range loop
      if c(k) /= null
       then Allocate_Coefficients(c(k),res(k));
      end if;
    end loop;
    cff := new Standard_Speelpenning_Convolutions.VecVecVec'(res);
  end Allocate_Coefficients;

  procedure Allocate_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                cff : out DoblDobl_Speelpenning_Convolutions.
                          Link_to_VecVecVec ) is

    use DoblDobl_Speelpenning_Convolutions;

    res : VecVecVec(c'range);

  begin
    for k in c'range loop
      if c(k) /= null
       then Allocate_Coefficients(c(k),res(k));
      end if;
    end loop;
    cff := new DoblDobl_Speelpenning_Convolutions.VecVecVec'(res);
  end Allocate_Coefficients;

  procedure Allocate_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                cff : out QuadDobl_Speelpenning_Convolutions.
                          Link_to_VecVecVec ) is

    use QuadDobl_Speelpenning_Convolutions;

    res : VecVecVec(c'range);

  begin
    for k in c'range loop
      if c(k) /= null
       then Allocate_Coefficients(c(k),res(k));
      end if;
    end loop;
    cff := new QuadDobl_Speelpenning_Convolutions.VecVecVec'(res);
  end Allocate_Coefficients;

  procedure Store_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                cff : in Standard_Complex_VecVecs.Link_to_VecVec ) is

    use Standard_Complex_Vectors;
    lnk,cfflnk : Link_to_Vector;

  begin
    if c.cst /= null then
      lnk := cff(0);
      for i in c.cst'range loop
        lnk(i) := c.cst(i);
      end loop;
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      cfflnk := cff(k);
      for i in lnk'range loop
        cfflnk(i) := lnk(i);
      end loop;
    end loop;
  end Store_Coefficients;

  procedure Store_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : in DoblDobl_Complex_VecVecs.Link_to_VecVec ) is

    use DoblDobl_Complex_Vectors;
    lnk,cfflnk : Link_to_Vector;

  begin
    if c.cst /= null then
      lnk := cff(0);
      for i in c.cst'range loop
        lnk(i) := c.cst(i);
      end loop;
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      cfflnk := cff(k);
      for i in lnk'range loop
        cfflnk(i) := lnk(i);
      end loop;
    end loop;
  end Store_Coefficients;

  procedure Store_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                cff : in QuadDobl_Complex_VecVecs.Link_to_VecVec ) is

    use QuadDobl_Complex_Vectors;
    lnk,cfflnk : Link_to_Vector;

  begin
    if c.cst /= null then
      lnk := cff(0);
      for i in c.cst'range loop
        lnk(i) := c.cst(i);
      end loop;
    end if;
    for k in c.cff'range loop
      lnk := c.cff(k);
      cfflnk := cff(k);
      for i in lnk'range loop
        cfflnk(i) := lnk(i);
      end loop;
    end loop;
  end Store_Coefficients;

  procedure Restore_Coefficients
              ( cff : in Standard_Complex_VecVecs.Link_to_VecVec;
                c : in Standard_Speelpenning_Convolutions.
                       Link_to_Circuit ) is

    use Standard_Complex_Vectors;
    lnk,cfflnk : Link_to_Vector;

  begin
    if cff(0) /= null then
      if c.cst /= null then
        lnk := cff(0);
        for i in c.cst'range loop
          c.cst(i) := lnk(i);
        end loop;
      end if;
    end if;
    for k in c.cff'range loop -- note : c.cff'first = 1, cff'first = 0
      if cff(k) /= null then
        lnk := c.cff(k);
        cfflnk := cff(k);
        for i in lnk'range loop
          lnk(i) := cfflnk(i);
        end loop;
      end if;
    end loop;
  end Restore_Coefficients;

  procedure Restore_Coefficients
              ( cff : in DoblDobl_Complex_VecVecs.Link_to_VecVec;
                c : in DoblDobl_Speelpenning_Convolutions.
                       Link_to_Circuit ) is

    use DoblDobl_Complex_Vectors;
    lnk,cfflnk : Link_to_Vector;

  begin
    if cff(0) /= null then
      if c.cst /= null then
        lnk := cff(0);
        for i in c.cst'range loop
          c.cst(i) := lnk(i);
        end loop;
      end if;
    end if;
    for k in c.cff'range loop -- note : c.cff'first = 1, cff'first = 0
      if cff(k) /= null then
        lnk := c.cff(k);
        cfflnk := cff(k);
        for i in lnk'range loop
          lnk(i) := cfflnk(i);
        end loop;
      end if;
    end loop;
  end Restore_Coefficients;

  procedure Restore_Coefficients
              ( cff : in QuadDobl_Complex_VecVecs.Link_to_VecVec;
                c : in QuadDobl_Speelpenning_Convolutions.
                       Link_to_Circuit ) is

    use QuadDobl_Complex_Vectors;
    lnk,cfflnk : Link_to_Vector;

  begin
    if cff(0) /= null then
      if c.cst /= null then
        lnk := cff(0);
        for i in c.cst'range loop
          c.cst(i) := lnk(i);
        end loop;
      end if;
    end if;
    for k in c.cff'range loop -- note : cff'first = 1, cff'first = 0
      if cff(k) /= null then
        lnk := c.cff(k);
        cfflnk := cff(k);
        for i in lnk'range loop
          lnk(i) := cfflnk(i);
        end loop;
      end if;
    end loop;
  end Restore_Coefficients;

  procedure Store_Coefficients
              ( c : in Standard_Speelpenning_Convolutions.Circuits;
                cff : in Standard_Speelpenning_Convolutions.
                         Link_to_VecVecVec ) is

    use Standard_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null
       then Store_Coefficients(c(k),cff(k));
      end if;
    end loop;
  end Store_Coefficients;

  procedure Store_Coefficients
              ( c : in DoblDobl_Speelpenning_Convolutions.Circuits;
                cff : in DoblDobl_Speelpenning_Convolutions.
                         Link_to_VecVecVec ) is

    use DoblDobl_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null
       then Store_Coefficients(c(k),cff(k));
      end if;
    end loop;
  end Store_Coefficients;

  procedure Store_Coefficients
              ( c : in QuadDobl_Speelpenning_Convolutions.Circuits;
                cff : in QuadDobl_Speelpenning_Convolutions.
                         Link_to_VecVecVec ) is

    use QuadDobl_Speelpenning_Convolutions;

  begin
    for k in c'range loop
      if c(k) /= null
       then Store_Coefficients(c(k),cff(k));
      end if;
    end loop;
  end Store_Coefficients;

  procedure Restore_Coefficients
              ( cff : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                c : in Standard_Speelpenning_Convolutions.Circuits ) is

    use Standard_Complex_VecVecs,Standard_Speelpenning_Convolutions;

  begin
    if cff /= null then
      for k in cff'range loop
        if cff(k) /= null and c(k) /= null
         then Restore_Coefficients(cff(k),c(k));
        end if;
      end loop;
    end if;
  end Restore_Coefficients;

  procedure Restore_Coefficients
              ( cff : in DoblDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                c : in DoblDobl_Speelpenning_Convolutions.Circuits ) is

    use DoblDobl_Complex_VecVecs,DoblDobl_Speelpenning_Convolutions;

  begin
    if cff /= null then
      for k in cff'range loop
        if cff(k) /= null and c(k) /= null
         then Restore_Coefficients(cff(k),c(k));
        end if;
      end loop;
    end if;
  end Restore_Coefficients;

  procedure Restore_Coefficients
              ( cff : in QuadDobl_Speelpenning_Convolutions.Link_to_VecVecVec;
                c : in QuadDobl_Speelpenning_Convolutions.Circuits ) is

    use QuadDobl_Complex_VecVecs,QuadDobl_Speelpenning_Convolutions;

  begin
    if cff /= null then
      for k in cff'range loop
        if cff(k) /= null and c(k) /= null
         then Restore_Coefficients(cff(k),c(k));
        end if;
      end loop;
    end if;
  end Restore_Coefficients;

-- EVALUATION OF STEP SIZE :

  function Step_Coefficient
              ( c : Standard_Complex_Vectors.Vector; t : double_float )
              return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for k in reverse 0..c'last-1 loop
      res := res*t + c(k);
    end loop;
    return res;
  end Step_Coefficient;

  function Step_Coefficient
              ( c : DoblDobl_Complex_Vectors.Vector; t : double_double )
              return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for k in reverse 0..c'last-1 loop
      res := res*t + c(k);
    end loop;
    return res;
  end Step_Coefficient;

  function Step_Coefficient
              ( c : QuadDobl_Complex_Vectors.Vector; t : quad_double )
              return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;

    res : Complex_Number := c(c'last);

  begin
    for k in reverse 0..c'last-1 loop
      res := res*t + c(k);
    end loop;
    return res;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( c : in Standard_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_float ) is

    use Standard_Complex_Vectors;

  begin
    if c.cst /= null
     then c.cst(0) := Step_Coefficient(c.cst.all,t);
    end if;
    for k in c.cff'range loop
      c.cff(k)(0) := Step_Coefficient(c.cff(k).all,t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( c : in DoblDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in double_double ) is

    use DoblDobl_Complex_Vectors;

  begin
    if c.cst /= null
     then c.cst(0) := Step_Coefficient(c.cst.all,t);
    end if;
    for k in c.cff'range loop
      c.cff(k)(0) := Step_Coefficient(c.cff(k).all,t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( c : in QuadDobl_Speelpenning_Convolutions.Link_to_Circuit;
                t : in quad_double ) is

    use QuadDobl_Complex_Vectors;

  begin
    if c.cst /= null
     then c.cst(0) := Step_Coefficient(c.cst.all,t);
    end if;
    for k in c.cff'range loop
      c.cff(k)(0) := Step_Coefficient(c.cff(k).all,t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                t : in double_float ) is
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                t : in double_double ) is
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

  procedure Step_Coefficient
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                t : in quad_double ) is
  begin
    for k in hom.crc'range loop
      Step_Coefficient(hom.crc(k),t);
    end loop;
  end Step_Coefficient;

-- RUNNING ONE STEP OF NEWTON'S METHOD :

  procedure LU_Newton_Step
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    use Standard_Complex_Numbers;

  begin
    Standard_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    Standard_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float ) is

    use Standard_Complex_Numbers;

  begin
    Standard_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    Standard_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    lufco(hom.vm(0).all,hom.dim,ipvt,rcond);
    if 1.0 + rcond /= 1.0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true ) is

    use Standard_Complex_Numbers;

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    Standard_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    Standard_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                sol : in out Standard_Complex_Vectors.Vector;
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float; verbose : in boolean := true ) is

    use Standard_Complex_Numbers;

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    Standard_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    Standard_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufco(hom.vm(0).all,hom.dim,ipvt,rcond);
    if 1.0 + rcond /= 1.0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    use DoblDobl_Complex_Numbers;

  begin
    DoblDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    DOblDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double ) is

    use DoblDobl_Complex_Numbers;

    one : constant double_double := create(1.0);

  begin
    DoblDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    DOblDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    lufco(hom.vm(0).all,hom.dim,ipvt,rcond);
    if one + rcond /= one then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true ) is

    use DoblDobl_Complex_Numbers;

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    DoblDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    DOblDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out DoblDobl_Complex_Vectors.Vector;
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double; verbose : in boolean := true ) is

    use DoblDobl_Complex_Numbers;

    one : constant double_double := create(1.0);

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    DoblDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    DOblDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufco(hom.vm(0).all,hom.dim,ipvt,rcond);
    if one + rcond /= one then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 ) is

    use QuadDobl_Complex_Numbers;

  begin
    QuadDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double ) is

    use QuadDobl_Complex_Numbers;

    one : constant quad_double := create(1.0);

  begin
    QuadDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    lufco(hom.vm(0).all,hom.dim,ipvt,rcond);
    if one + rcond /= one then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; verbose : in boolean := true ) is

    use QuadDobl_Complex_Numbers;

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    QuadDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufac(hom.vm(0).all,hom.dim,ipvt,info);
    if info = 0 then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

  procedure LU_Newton_Step
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double; verbose : in boolean := true ) is

    use QuadDobl_Complex_Numbers;

    one : constant quad_double := create(1.0);

  begin
    if verbose
     then put_line(file,"The solution on input : "); put_line(file,sol);
    end if;
    QuadDobl_Speelpenning_Convolutions.Compute(hom.pwt,hom.mxe,sol);
    QuadDobl_Speelpenning_Convolutions.EvalDiff(hom,sol);
    for k in dx'range loop 
      dx(k) := -hom.yv(k)(0);
    end loop;
    if verbose
     then put_line(file,"The function value :"); put_line(file,dx);
    end if;
    lufco(hom.vm(0).all,hom.dim,ipvt,rcond);
    if one + rcond /= one then
      lusolve(hom.vm(0).all,hom.dim,ipvt,dx);
      if verbose
       then put_line(file,"The update : "); put_line(file,dx);
      end if;
      for k in dx'range loop
        sol(k) := sol(k) + dx(k);
      end loop;
      if verbose
       then put_line(file,"The updated solution : "); put_line(file,sol);
      end if;
    end if;
  end LU_Newton_Step;

-- RUNNING MANY STEPS OF NEWTON'S METHOD :

  procedure LU_Newton_Steps
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0 ) is

    use Standard_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(hom,psv.sol,dx,ipvt,info);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := Standard_Complex_Vector_Norms.Max_Norm(dx);
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (maxdx = 0.0) or (mixres = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float; fail : out boolean;
                extra : in integer32 := 0 ) is

    use Standard_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true;
    for k in 1..maxit+extra loop
      LU_Newton_Step(hom,psv.sol,dx,ipvt,rcond);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := Standard_Complex_Vector_Norms.Max_Norm(dx);
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (maxdx = 0.0) or (mixres = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
    nbrit := maxit;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true ) is

    use Standard_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,hom,psv.sol,dx,ipvt,info,verbose);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := Standard_Complex_Vector_Norms.Max_Norm(dx);
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put(file,"after step "); put(file,k,1);
        put(file,", maxdx :"); put(file,maxdx,3);
        put(file,", mixres :"); put(file,mixres,3); new_line(file);
      end if;
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (maxdx = 0.0) or (mixres = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                psv : in out Standard_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_float; 
                dx : out Standard_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true ) is

    use Standard_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,hom,psv.sol,dx,ipvt,rcond,verbose);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := Standard_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := Standard_Complex_Vector_Norms.Max_Norm(dx);
      mixres := Standard_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put(file,"after step "); put(file,k,1);
        put(file,", maxdx :"); put(file,maxdx,3);
        put(file,", mixres :"); put(file,mixres,3); new_line(file);
      end if;
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (maxdx = 0.0) or (mixres = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0 ) is

    use DoblDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(hom,psv.sol,dx,ipvt,info);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := DoblDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := DoblDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := DoblDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hi_part(maxdx) = 0.0)
                              or (hi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double; fail : out boolean;
                extra : in integer32 := 0 ) is

    use DoblDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(hom,psv.sol,dx,ipvt,rcond);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := DoblDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := DoblDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := DoblDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hi_part(maxdx) = 0.0)
                              or (hi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true ) is

    use DoblDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,hom,psv.sol,dx,ipvt,info,verbose);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := DoblDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := DoblDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := DoblDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put(file,"after step "); put(file,k,1);
        put(file,", maxdx : "); put(file,maxdx,3);
        put(file,", mixres : "); put(file,mixres,3); new_line(file);
      end if;
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hi_part(maxdx) = 0.0)
                              or (hi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out DoblDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out double_double; 
                dx : out DoblDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_double; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true ) is

    use DoblDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,hom,psv.sol,dx,ipvt,rcond,verbose);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := DoblDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := DoblDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := DoblDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put(file,"after step "); put(file,k,1);
        put(file,", maxdx : "); put(file,maxdx,3);
        put(file,", mixres : "); put(file,mixres,3); new_line(file);
      end if;
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hi_part(maxdx) = 0.0)
                              or (hi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0 ) is

    use QuadDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(hom,psv.sol,dx,ipvt,info);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := QuadDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := QuadDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := QuadDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hihi_part(maxdx) = 0.0)
                              or (hihi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double; fail : out boolean;
                extra : in integer32 := 0 ) is

    use QuadDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(hom,psv.sol,dx,ipvt,rcond);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := QuadDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := QuadDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := QuadDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hihi_part(maxdx) = 0.0)
                              or (hihi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true ) is

    use QuadDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,hom,psv.sol,dx,ipvt,info,verbose);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := QuadDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := QuadDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := QuadDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put(file,"after step "); put(file,k,1);
        put(file,", maxdx : "); put(file,maxdx,3);
        put(file,", mixres : "); put(file,mixres,3); new_line(file);
      end if;
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hihi_part(maxdx) = 0.0)
                              or (hihi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

  procedure LU_Newton_Steps
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                psv : in out QuadDobl_Predictor_Convolutions.Predictor_Vectors;
                maxit : in integer32; nbrit : out integer32;
                tol : in double_float; maxdx,mixres : out quad_double; 
                dx : out QuadDobl_Complex_Vectors.Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out quad_double; fail : out boolean;
                extra : in integer32 := 0;
                verbose : in boolean := true ) is

    use QuadDobl_Speelpenning_Convolutions;

    cntextra : integer32 := 0;

  begin
    fail := true; nbrit := maxit;
    for k in 1..maxit+extra loop
      LU_Newton_Step(file,hom,psv.sol,dx,ipvt,rcond,verbose);
      psv.res := Eval(hom.crc,psv.sol);
      psv.radsol := QuadDobl_Mixed_Residuals.AbsVal(psv.sol);
      psv.radres := Eval(abh.crc,psv.radsol);
      maxdx := QuadDobl_Complex_Vector_Norms.Max_Norm(dx);
      mixres := QuadDobl_Mixed_Residuals.Mixed_Residual(psv.res,psv.radres);
      if verbose then
        put(file,"after step "); put(file,k,1);
        put(file,", maxdx : "); put(file,maxdx,3);
        put(file,", mixres : "); put(file,mixres,3); new_line(file);
      end if;
      if maxdx <= tol and mixres <= tol then -- convergence
        if (cntextra = extra) or (hihi_part(maxdx) = 0.0)
                              or (hihi_part(mixres) = 0.0)
         then nbrit := k; fail := false; exit;
        end if;
        cntextra := cntextra + 1; -- do an extra Newton step
      end if;
    end loop;
  end LU_Newton_Steps;

end Corrector_Convolutions;
