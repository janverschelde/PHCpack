with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Integer_Vectors_IO;
with Standard_Integer_VecVecs_IO;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_LaurSys_Container;
with Double_VecVecs_Container;
with DCMPLX_VecVecs_Container;
with Real_Powered_Homotopy_IO;

package body Double_Puiseux_Structures is

  procedure Indexing_Series 
              ( powers : out Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                coeffs : out Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                vrblvl : in integer32 := 0) is

    p : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
      := Standard_LaurSys_Container.Retrieve;
    pwr : constant Standard_Floating_VecVecs.Link_to_Array_of_VecVecs
        := Double_VecVecs_Container.Get(vrblvl-1);
    cff : constant Standard_Complex_VecVecs.Link_to_Array_of_VecVecs
        := DCMPLX_VecVecs_Container.Get(vrblvl-1);

    use Standard_Complex_Laur_Systems;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_structures.Indexing_Series ...");
    end if;
    declare
      afv : Standard_Floating_VecVecs.Array_of_VecVecs(pwr'range);
      acv : Standard_Complex_VecVecs.Array_of_VecVecs(cff'range);
    begin
      for i in afv'range loop
        declare
          afvi : Standard_Floating_VecVecs.VecVec(pwr(i)'range);
          acvi : Standard_Complex_VecVecs.VecVec(cff(i)'range);
        begin
          for j in afvi'range loop -- j-th power series coefficient
            declare
              afvij : constant Standard_Floating_Vectors.Link_to_Vector
                    := pwr(i)(j);
              acvij : constant Standard_Complex_Vectors.Link_to_Vector
                    := cff(i)(j);
              dim : constant integer32 := afvij'last;
              npw : Standard_Floating_Vectors.Vector(1..dim-1);
              ncf : Standard_Complex_Vectors.Vector(0..dim-1);
            begin
              for k in npw'range loop
                npw(k) := afvij(k+1);
              end loop;
              afvi(j) := new Standard_Floating_Vectors.Vector'(npw);
              for k in ncf'range loop
                ncf(k) := acvij(k+1);
              end loop;
              acvi(j) := new Standard_Complex_Vectors.Vector'(ncf);
            end;
          end loop;
          afv(i) := new Standard_Floating_VecVecs.VecVec'(afvi);
          acv(i) := new Standard_Complex_VecVecs.VecVec'(acvi);
        end;
      end loop;
      powers := new Standard_Floating_VecVecs.Array_of_VecVecs'(afv);
      coeffs := new Standard_Complex_VecVecs.Array_of_VecVecs'(acv);
    end;
    if p = null then
      put_line("The Laurent system is not defined!");
    elsif vrblvl > 0 then
      put_line("The Laurent system :"); put(p.all);
    end if;
    if vrblvl > 0 then
      Real_Powered_Homotopy_IO.put_line
        (p'last,p'last,pwr(pwr'first)'last,p.all,coeffs.all,powers.all,
         vrblvl=>vrblvl-1);
    end if;
  end Indexing_Series;

  procedure Show_Data ( vrblvl : in integer32 := 0) is

    p : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
      := Standard_LaurSys_Container.Retrieve;
    pwr : constant Standard_Floating_VecVecs.Link_to_Array_of_VecVecs
        := Double_VecVecs_Container.Get(vrblvl-1);
    cff : constant Standard_Complex_VecVecs.Link_to_Array_of_VecVecs
        := DCMPLX_VecVecs_Container.Get(vrblvl-1);

    use Standard_Complex_Laur_Systems;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_structures.Show_Data ...");
    end if;
    if p = null
     then put_line("The Laurent system is not defined!");
     else put_line("The Laurent system :"); put(p.all);
    end if;
    Real_Powered_Homotopy_IO.put_line
      (p'last,p'last,pwr(pwr'first)'last,p.all,cff.all,pwr.all,
       vrblvl=>vrblvl-1);
  end Show_Data;

  function Is_Zero ( v : Standard_Integer_Vectors.Link_to_Vector )
                   return boolean is
  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Is_Variable ( deg : Standard_Integer_Vectors.Link_to_Vector;
                         var : integer32 ) return boolean is
  begin
    if deg(var) /= 1 then
      return false;
    else
      for i in deg'first..var-1 loop
        if deg(i) /= 0
         then return false;
        end if;
      end loop;
      for i in var+1..deg'last loop
        if deg(i) /= 0
         then return false;
        end if;
      end loop;
      return true;
    end if;
  end Is_Variable;

  function Is_Diagonal_Binomial_Homotopy
              ( p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                vrblvl : integer32 := 0 ) return boolean is

    cnt,idx : integer32;

    procedure Visit_Term ( t : in Standard_Complex_Laurentials.Term;
                           continue : out boolean ) is

      deg : constant Standard_Integer_Vectors.Link_to_Vector
          := Standard_Integer_Vectors.Link_to_Vector(t.dg);

    begin
      if Is_Zero(deg) then
        cnt := cnt + 1;
      elsif Is_Variable(deg,idx) then
        cnt := cnt + 1;
      end if;
      continue := true;
    end Visit_Term;

    procedure Visit_Terms is
      new Standard_Complex_Laurentials.Visiting_Iterator(Visit_Term);

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_structures.");
      put_line("Is_Diagonal_Binomial_Homotopy ...");
    end if;
    for i in p'range loop
      idx := i; cnt := 0;
      Visit_Terms(p(i));
      if cnt < 2 then
        if vrblvl > 0
         then put_line("Returning false.");
        end if;
        return false;
      end if;
    end loop;
    if vrblvl > 0
     then put_line("Returning true.");
    end if;
    return true;
  end Is_Diagonal_Binomial_Homotopy;

  function Extract_Leading_Coefficients
             ( cffs : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
               hdg : Standard_Integer_VecVecs.Array_of_VecVecs;
               vrblvl : integer32 := 0 )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(cffs'range);

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_structures.");
      put_line("Extract_Leading_Coefficients ...");
    end if;
    for i in cffs'range loop
      declare
        moncff : constant Standard_Complex_VecVecs.Link_to_VecVec := cffs(i);
        rescff : Standard_Complex_Vectors.Vector(moncff'range);
        ideg : constant Standard_Integer_VecVecs.Link_to_VecVec := hdg(i);
        sup : Standard_Integer_Vectors.Link_to_Vector;
      begin
        for j in moncff'range loop
          sup := ideg(j);
          if not Is_Zero(sup)
           then rescff(j) := moncff(j)(1);
           else rescff(j) := moncff(j)(0);
          end if;
        end loop;
        res(i) := new Standard_Complex_Vectors.Vector'(rescff);
      end;
    end loop;
    return res;
  end Extract_Leading_Coefficients;

  function Extract_Leading_Powers
             ( pwrs : Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
               hdg : Standard_Integer_VecVecs.Array_of_VecVecs;
               vrblvl : integer32 := 0 )
             return Standard_Floating_VecVecs.VecVec is

    res : Standard_Floating_VecVecs.VecVec(pwrs'range);

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_structures.");
      put_line("Extract_Leading_Powers ...");
    end if;
    for i in pwrs'range loop
      declare
        monpwr : constant Standard_Floating_VecVecs.Link_to_VecVec := pwrs(i);
        respwr : Standard_Floating_Vectors.Vector(monpwr'range);
        ideg : constant Standard_Integer_VecVecs.Link_to_VecVec := hdg(i);
        sup : Standard_Integer_Vectors.Link_to_Vector;
      begin
        for j in monpwr'range loop
          sup := ideg(j);
          if not Is_Zero(sup)
           then respwr(j) := monpwr(j)(1);
           else respwr(j) := 0.0;
          end if;
        end loop;
        res(i) := new Standard_Floating_Vectors.Vector'(respwr);
      end;
    end loop;
    return res;
  end Extract_Leading_Powers;

  procedure Normalize_Binomial_Homotopy 
              ( hdg : in out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : in out Standard_Complex_VecVecs.VecVec;
                hct : in out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_structures.");
      put_line("Normalize_Binomial_Homotopy ...");
      put_line("before normalizing ...");
      for i in hdg'range loop
        put("support of polynomial "); put(i,1); put_line(" :");
        Standard_Integer_VecVecs_IO.put(hdg(i));
      end loop;
      for i in hcf'range loop
        put("leading terms of polynomial "); put(i,1); put_line(" :");
        for j in hcf(i)'range loop
          put(hcf(i)(j)); put(" t^"); put(hct(i)(j)); new_line;
        end loop;
      end loop;
    end if;
    for i in hdg'range loop -- fixing the variable term
      declare
        idg : constant Standard_Integer_VecVecs.Link_to_VecVec := hdg(i);
        sup : Standard_Integer_Vectors.Link_to_Vector;
        idx : integer32 := 0;
      begin
        for j in idg'range loop
          sup := idg(j);
          if Is_Variable(sup,i)
           then idx := j;
          end if;
          exit when idx > 0;
        end loop;
        if idx /= idg'first then -- swapping is needed
          idg(idx) := idg(idg'first);
          idg(idg'first) := sup; 
          hcf(i)(idx) := hcf(i)(hcf(i)'first);
          hcf(i)(hcf(i)'first) := Create(1.0); -- fixed coefficient
          hct(i)(idx) := hct(i)(hct(i)'first);
          hct(i)(hct(i)'first) := 0.0; -- fixed power
        else -- fixing must happen regardless of swapping
          hcf(i)(hcf(i)'first) := Create(1.0); -- fixed coefficient
          hct(i)(hct(i)'first) := 0.0; -- fixed power
        end if;
      end;
    end loop;
    if vrblvl > 0 then
      put_line("after normalizing ...");
      for i in hdg'range loop
        put("support of polynomial "); put(i,1); put_line(" :");
        Standard_Integer_VecVecs_IO.put(hdg(i));
      end loop;
      for i in hcf'range loop
        put("leading terms of polynomial "); put(i,1); put_line(" :");
        for j in hcf(i)'range loop
          put(hcf(i)(j)); put(" t^"); put(hct(i)(j)); new_line;
        end loop;
      end loop;
    end if;
  end Normalize_Binomial_Homotopy;

  function Truncation_Index
              ( pwrs : Standard_Floating_VecVecs.Link_to_VecVec )
              return integer32 is

    pwr : Standard_Floating_Vectors.Link_to_Vector := pwrs(pwrs'first);
    res : integer32 := pwr'last;

  begin
    for i in pwrs'first+1..pwrs'last loop
      pwr := pwrs(i);
      if pwr'last > res
       then res := pwr'last;
      end if;
    end loop;
    return res;
  end Truncation_Index;

  function Is_Zero ( c : Complex_Number ) return boolean is
  begin
    return (REAL_PART(c) = 0.0 and IMAG_PART(c) = 0.0);
  end Is_Zero;

  procedure Product_Monomials
              ( p : in Standard_Complex_Laurentials.Poly;
                cffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                pwrs : in Standard_Floating_VecVecs.Link_to_VecVec;
                hdg : out Standard_Integer_VecVecs.Link_to_VecVec;
                hcf : out Standard_Complex_Vectors.Link_to_Vector;
                hct : out Standard_Floating_Vectors.Link_to_Vector;
                vrblvl : in integer32 := 0 ) is

    nbr : constant integer32
        := integer32(Standard_Complex_Laurentials.Number_of_Terms(p));
    deg : constant integer32 := Truncation_Index(pwrs);
    nbterms : constant integer32 := nbr*deg;
    vdg : Standard_Integer_VecVecs.VecVec(1..nbterms);
    vcf : Standard_Complex_Vectors.Vector(1..nbterms)
        := (1..nbterms => Standard_Complex_Numbers.create(0.0));
    vct : Standard_Floating_Vectors.Vector(1..nbterms)
        := (1..nbterms => 0.0);
    idxtrm : integer32 := 0; -- index in hdg
    idxcff : integer32 := cffs'first-1; -- index in cffs and pwrs

    procedure Visit_Term ( t : in Standard_Complex_Laurentials.Term;
                           continue : out boolean ) is

      deg : constant Standard_Integer_Vectors.Link_to_Vector
          := Standard_Integer_Vectors.Link_to_Vector(t.dg);
      cff : Standard_Complex_Vectors.Link_to_Vector;
      pwr : Standard_Floating_Vectors.Link_to_Vector;

    begin
      if vrblvl > 0 then
        put("idxcff : "); put(idxcff,1); 
        put(", idxtrm : "); put(idxtrm,1); new_line;
      end if;
      idxcff := idxcff + 1;
      cff := cffs(idxcff);
      pwr := pwrs(idxcff);
      if vrblvl > 0 then
        put("coefficient series "); put(idxcff,1); put_line(" :");
        put(cff(0)); new_line;
        for i in pwr'range loop
          put(cff(i)); put("  ");
          put(pwr(i)); new_line;
        end loop;
      end if;
      if not Is_Zero(cff(0)) then -- separate treatment for constant
        idxtrm := idxtrm + 1;
        vdg(idxtrm) := new Standard_Integer_Vectors.Vector'(deg.all);
        vcf(idxtrm) := cff(0);
        vct(idxtrm) := 0.0;
        if vrblvl > 0 then
          put("at idxterm "); put(idxtrm,1); put_line(" :");
          put(vcf(idxtrm)); put("  ");
          put(vct(idxtrm)); put("  ");
          Standard_Integer_Vectors_IO.put(vdg(idxtrm)); new_line;
        end if;
      end if;
      for i in pwr'range loop
        if not Is_Zero(cff(i)) then
          idxtrm := idxtrm + 1;
          vdg(idxtrm) := new Standard_Integer_Vectors.Vector'(deg.all);
          vcf(idxtrm) := cff(i);
          vct(idxtrm) := pwr(i);
          if vrblvl > 0 then
            put("at idxterm "); put(idxtrm,1); put_line(" :");
            put(vcf(idxtrm)); put("  ");
            put(vct(idxtrm)); put("  ");
            Standard_Integer_Vectors_IO.put(vdg(idxtrm)); new_line;
          end if;
        end if;
      end loop;
      continue := true;
    end Visit_Term;

    procedure Visit_Terms is
      new Standard_Complex_Laurentials.Visiting_Iterator(Visit_Term);

  begin
    if vrblvl > 0 then
      put_line("-> in double_puiseux_structures.Product_Monomials ...");
      put("number of monomials : "); put(nbr,1);
      put(", truncation index : "); put(deg,1); new_line;
      put("number of terms : "); put(nbterms,1); new_line;
    end if;
    Visit_Terms(p);
    if vrblvl > 0
     then put("idxtrm : "); put(idxtrm,1); new_line;
    end if;
    hdg := new Standard_Integer_VecVecs.VecVec'(vdg(1..idxtrm));
    hcf := new Standard_Complex_Vectors.Vector'(vcf(1..idxtrm));
    hct := new Standard_Floating_Vectors.Vector'(vct(1..idxtrm));
  end Product_Monomials;

  procedure Normalize_Product_Homotopy
              ( p : in Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                cffs : in Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
                pwrs : in Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                hct : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_structures.");
      put_line("Product_Coefficients_Powers ...");
    end if;
    for i in p'range loop
      Product_Monomials(p(i),cffs(i),pwrs(i),hdg(i),hcf(i),hct(i),vrblvl);
    end loop;
  end Normalize_Product_Homotopy;

end Double_Puiseux_Structures;
