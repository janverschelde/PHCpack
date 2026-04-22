with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_IO;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_IO;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Double_Newton_Puiseux;
with Assignments_in_Ada_and_C;
with Standard_LaurSys_Container;
with Double_VecVecs_Container;
with DCMPLX_VecVecs_Container;
with Real_Powered_Homotopy;
with Real_Powered_Homotopy_IO;

package body Double_Puiseux_Interface is

  powers : Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
  coeffs : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;

  procedure Indexing_Series ( vrblvl : in integer32 := 0) is

  -- DESCRIPTION :
  --   Copies the arrays of vectors of vectors into powers
  --   and coefficients, adjusting the start of the indices
  --   of the coefficient vectors and shifting the powers.

    p : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
      := Standard_LaurSys_Container.Retrieve;
    pwr : constant Standard_Floating_VecVecs.Link_to_Array_of_VecVecs
        := Double_VecVecs_Container.Get(vrblvl-1);
    cff : constant Standard_Complex_VecVecs.Link_to_Array_of_VecVecs
        := DCMPLX_VecVecs_Container.Get(vrblvl-1);

    use Standard_Complex_Laur_Systems;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Copy_Series ...");
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

  -- DESCRIPTION :
  --   Retrieves the Laurent polynomial system and the corresponding
  --   coefficients of the series coefficients and the arrays of vectors.
  --   Writes the series system to screen as a test.

    p : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
      := Standard_LaurSys_Container.Retrieve;
    pwr : constant Standard_Floating_VecVecs.Link_to_Array_of_VecVecs
        := Double_VecVecs_Container.Get(vrblvl-1);
    cff : constant Standard_Complex_VecVecs.Link_to_Array_of_VecVecs
        := DCMPLX_VecVecs_Container.Get(vrblvl-1);

    use Standard_Complex_Laur_Systems;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Show_Data ...");
    end if;
    if p = null
     then put_line("The Laurent system is not defined!");
     else put_line("The Laurent system :"); put(p.all);
    end if;
    Real_Powered_Homotopy_IO.put_line
      (p'last,p'last,pwr(pwr'first)'last,p.all,cff.all,pwr.all,
       vrblvl=>vrblvl-1);
  end Show_Data;

  function Linear_Solver
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Linear_Solver ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      nbr : constant integer32 := integer32(v_a(v_a'first));
    begin
      if vrblvl > 0 then
        put("  nbr : "); put(nbr,1); put_line(" ...");
        Show_Data(vrblvl);
      end if;
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_puiseux_interface.");
        put_line("Linear_Solver.");
      end if;
      return -1;
  end Linear_Solver;

  function Extract_Constant_Coefficients
             ( cffs : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
               vrblvl : integer32 := 0 )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   In a binomial Laurent homotopy, in diagonal test format,
  --   the constant coefficients of the solution series are exposed
  --   in the constant coefficients of the homotopy, given by cffs.
  --   Returns the constants of the solution series.

  -- REQUIRED :
  --   The coefficients have been properly indexed at zero.

    res : Standard_Complex_Vectors.Vector(cffs'range);

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_interface.");
      put_line("Extract_Constant_Coefficients ...");
    end if;
    for i in cffs'range loop
      declare
        moncff : constant Standard_Complex_VecVecs.Link_to_VecVec := cffs(i);
        cst : Complex_Number;
        done : boolean := false;
      begin
        for j in moncff'range loop
          cst := moncff(j)(0);
         -- constant should not be zero and not be equal to one
          if REAL_PART(cst) /= 0.0 or IMAG_PART(cst) /= 0.0 then
            if REAL_PART(cst) /= 1.0 and IMAG_PART(cst) /= 0.0
             then res(i) := -cst; done := true;
            end if;
          end if;
          exit when done;
        end loop;
      end;
    end loop;
    return res;
  end Extract_Constant_Coefficients;

  function Extract_Solution_Constants
             ( dim : integer32; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   For Newton's method on power series to start properly,
  --   the constant coefficients of the solution series need
  --   to be provided in the c parameter of the interface.
  --   Return a vector of dimension dim, extracted from the
  --   2*dim double values of c.

    res : Standard_Complex_Vectors.Vector(1..dim)
        := (1..dim => Standard_Complex_Numbers.create(0.0));
    ddm : constant Interfaces.C.size_t := Interfaces.C.size_t(2*dim);
    use Interfaces.C;
    sol : C_Double_Array(0..ddm-1)
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(ddm));
    idx : Interfaces.C.size_t := 0;
    cre,cim : double_float;

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_interface.");
      put_line("Extract_Constant_Coefficients ...");
    end if;
    for i in res'range loop
      cre := double_float(sol(idx)); idx := idx + 1;
      cim := double_float(sol(idx)); idx := idx + 1;
      res(i) := Standard_Complex_Numbers.create(cre,cim);
    end loop;
    if vrblvl > 0 then
      put_line("The extracted constants :");
      put_line(res);
    end if;
    return res;
  end Extract_Solution_Constants;

  function Is_Zero ( v : Standard_Integer_Vectors.Link_to_Vector )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if all elements in v are zero.

  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function Is_Variable ( v : Standard_Integer_Vectors.Link_to_Vector )
                       return boolean is

  -- DESCRIPTION :
  --   Returns true if all elements in v are zero,
  --   except for exactly one component which should be one.

    cnt : integer32 := 0;

  begin
    for i in v'range loop
      if v(i) /= 0 then
        if v(i) /= 1
         then return false;
         else cnt := cnt + 1;
        end if;
        if cnt > 1
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Is_Variable;

  function Extract_Leading_Coefficients
             ( cffs : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
               hdg : Standard_Integer_VecVecs.Array_of_VecVecs;
               vrblvl : integer32 := 0 )
             return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Given properly indexed coefficients of a Laurent homotopy,
  --   extracts the coefficients at index 1.
  --   Although this is not correct when the monomial is a variable,
  --   it gets fixed when normalizing the binomial homotopy.
  --   If the corresponding degrees in hdg are zero,
  --   that is: we have a constant, then the power should be zero.

    res : Standard_Complex_VecVecs.VecVec(cffs'range);

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_interface.");
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

  -- DESCRIPTION :
  --   Given properly indexed powers of a Laurent homotopy,
  --   extracts the powers at index 1.
  --   Although this is not correct when the monomial is a variable,
  --   it gets fixed when normalizing the binomial homotopy.
  --   If the corresponding degrees in hdg are zero,
  --   that is: we have a constant, then the power should be zero.

    res : Standard_Floating_VecVecs.VecVec(pwrs'range);

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_interface.");
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

  function Is_Variable ( deg : Standard_Integer_Vectors.Link_to_Vector;
                         var : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the degrees in deg correspond to the variable
  --   with index given in var.

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

  function Is_Binomial_Homotopy
              ( hdg : Standard_Integer_VecVecs.Array_of_VecVecs;
                vrblvl : integer32 := 0 ) return boolean is

  -- DESCRIPTION :
  --   A binomial contains in its k-th equation the monomials x(k)
  --   and a constant.  This function runs through the supports
  --   and checks if this condition is satified.

  begin
    for k in hdg'range loop
      declare
        kdg : constant Standard_Integer_VecVecs.Link_to_VecVec := hdg(k);
        cnt : integer32 := 0; -- counts zero and variable degrees
      begin
        for j in kdg'range loop
          declare
            deg : constant Standard_Integer_Vectors.Link_to_Vector := kdg(j);
          begin
            if Is_Zero(deg) then
              cnt := cnt + 1;
            elsif Is_Variable(deg) then
              cnt := cnt + 1;
            end if;
          end;
        end loop;
        if cnt < 2
         then return false;
        end if;
      end;
    end loop;
    return true;
  end Is_Binomial_Homotopy;

  procedure Normalize_Binomial_Homotopy 
              ( hdg : in out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : in out Standard_Complex_VecVecs.VecVec;
                hct : in out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   In a binomial homotopy, the first term is a variable,
  --   with index equal to the index of the polynomial,
  --   and the second term is the constant term, the solution series.

  begin
    if vrblvl > 0 then
      put("-> in double_puiseux_interface.");
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

  procedure Run_Newton_Steps
              ( nbr : in integer32; c : C_dblarrs.Pointer;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Extracts the real powered Laurent homotopy data
  --   and then runs as many Newton steps as the value of nbr.
  --   Assigns the coefficients of the computed series to c.

    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    tol : constant double_float := 1.0E-12;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Run_Newton_Steps ...");
    end if;
    p := Standard_LaurSys_Container.Retrieve;
    declare
      hdg : Standard_Integer_VecVecs.Array_of_VecVecs(p'range)
          := Real_Powered_Homotopy.Supports(p.all,vrblvl-1);
      hcf : Standard_Complex_VecVecs.VecVec(p'range)
          := Extract_Leading_Coefficients(coeffs,hdg,vrblvl-1);
      hct : Standard_Floating_VecVecs.VecVec(p'range)
          := Extract_Leading_Powers(powers,hdg,vrblvl-1);
      zc0 : constant Standard_Complex_Vectors.Vector(p'range)
         -- := Extract_Constant_Coefficients(coeffs,vrblvl-1);
          := Extract_Solution_Constants(p'last,c,vrblvl-1);
      zc1,zc2,zc3,zc4 : Standard_Complex_Vectors.Vector(p'range);
      pw1,pw2,pw3,pw4 : Standard_Floating_Vectors.Vector(p'range);
      isb : constant boolean := Is_Binomial_Homotopy(hdg,vrblvl);
    begin
      if isb
       then Normalize_Binomial_Homotopy(hdg,hcf,hct,vrblvl);
      end if;
      if vrblvl > 0 then
        put_line("The constants of the solution series :");
        put_line(zc0);
      end if;
      Double_Newton_Puiseux.Diagonal_Newton_Steps
        (hcf,hct,hdg,zc0,nbr,zc1,zc2,zc3,zc4,pw1,pw2,pw3,pw4,tol,vrblvl+1);
      if vrblvl > 0 then
        put_line("constant terms of solution series :");
        for j in zc0'range loop
          put(zc0(j)); new_line;
        end loop;
        put_line("first terms of solution series :");
        for j in zc1'range loop
          put(zc1(j)); put(" t^"); put(pw1(j)); new_line;
        end loop;
        if nbr >= 2 then
          put_line("second terms of solution :");
          for j in zc2'range loop
            put(zc2(j)); put(" t^"); put(pw2(j)); new_line;
          end loop;
          if nbr >= 3 then
            put_line("third terms of solution :");
            for j in zc3'range loop
              put(zc3(j)); put(" t^"); put(pw3(j)); new_line;
            end loop;
            if nbr >= 4 then
              put_line("fourth terms of solution :");
              for j in zc4'range loop
                put(zc4(j)); put(" t^"); put(pw4(j)); new_line;
              end loop;
            end if;        
          end if;
        end if;
      end if;
      declare
        dim : constant integer32 := zc0'last;
        outsize : constant integer32 := 2*dim + 3*dim*nbr;
        result : Standard_Floating_Vectors.Vector(1..outsize);
        idx : integer32 := 1;
      begin
        if vrblvl > 0 then
          put("Assignment to vector of size "); put(outsize,1);
          put_line(" ...");
        end if;
        for i in 1..dim loop
          result(idx) := REAL_PART(zc0(i)); idx := idx + 1;
          result(idx) := IMAG_PART(zc0(i)); idx := idx + 1;
        end loop;
        for i in 1..dim loop
          result(idx) := REAL_PART(zc1(i)); idx := idx + 1;
          result(idx) := IMAG_PART(zc1(i)); idx := idx + 1;
          result(idx) := pw1(i); idx := idx + 1;
        end loop;
        if nbr >= 2 then
          for i in 1..dim loop
            result(idx) := REAL_PART(zc2(i)); idx := idx + 1;
            result(idx) := IMAG_PART(zc2(i)); idx := idx + 1;
            result(idx) := pw2(i); idx := idx + 1;
          end loop;
          if nbr >= 3 then
            for i in 1..dim loop
              result(idx) := REAL_PART(zc3(i)); idx := idx + 1;
              result(idx) := IMAG_PART(zc3(i)); idx := idx + 1;
              result(idx) := pw3(i); idx := idx + 1;
            end loop;
            if nbr >= 4 then
              for i in 1..dim loop
                result(idx) := REAL_PART(zc4(i)); idx := idx + 1;
                result(idx) := IMAG_PART(zc4(i)); idx := idx + 1;
                result(idx) := pw4(i); idx := idx + 1;
              end loop;
            end if;   
          end if;
        end if;
        if vrblvl > 0 then
          put_line("the resulting numbers :");
          Standard_Floating_Vectors_IO.put_line(result);
        end if;
        Assignments_in_Ada_and_C.Assign(result,c);
      end;
    end; 
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_puiseux_interface.");
        put_line("Run_Newton_Steps."); raise;
      end if;
  end Run_Newton_Steps;

  function Newton_Steps
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Newton_Steps ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      nbr : constant integer32 := integer32(v_a(v_a'first));
    begin
      if vrblvl > 0 then
        put("  nbr : "); put(nbr,1); put_line(" ...");
      end if;
     -- Show_Data(vrblvl);
      Indexing_Series(vrblvl);
      Run_Newton_Steps(nbr,c,vrblvl);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_puiseux_interface.");
        put_line("Newton_Steps.");
      end if;
      return -1;
  end Newton_Steps;

end Double_Puiseux_Interface;
