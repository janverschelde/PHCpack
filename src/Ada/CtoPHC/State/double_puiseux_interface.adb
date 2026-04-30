with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_IO;      use Standard_Floating_Numbers_IO;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_IO;       use Standard_Complex_Numbers_IO;
with Standard_Integer_Vectors_IO;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_IO;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;
with Double_Newton_Puiseux;
with Assignments_in_Ada_and_C;
with Standard_LaurSys_Container;
with Real_Powered_Homotopy;
with Double_Puiseux_Structures;

package body Double_Puiseux_Interface is

  powers : Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
  coeffs : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;

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
        Double_Puiseux_Structures.Show_Data(vrblvl);
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
    sol : constant C_Double_Array(0..ddm-1)
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

  procedure Run_Newton_Steps
              ( nbr : in integer32; c : C_dblarrs.Pointer;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Extracts the real powered Laurent homotopy data
  --   and then runs as many Newton steps as the value of nbr.
  --   Assigns the coefficients of the computed series to c.

    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    tol : constant double_float := 1.0E-12;
    isbinhom : boolean;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Run_Newton_Steps ...");
    end if;
    p := Standard_LaurSys_Container.Retrieve;
    isbinhom
      := Double_Puiseux_Structures.Is_Diagonal_Binomial_Homotopy(p,vrblvl);
    declare
      hdg : Standard_Integer_VecVecs.Array_of_VecVecs(p'range);
      hcf : Standard_Complex_VecVecs.VecVec(p'range);
      hct : Standard_Floating_VecVecs.VecVec(p'range);
      zc0 : constant Standard_Complex_Vectors.Vector(p'range)
          := Extract_Solution_Constants(p'last,c,vrblvl-1);
      zc1,zc2,zc3,zc4 : Standard_Complex_Vectors.Vector(p'range);
      pw1,pw2,pw3,pw4 : Standard_Floating_Vectors.Vector(p'range);
      use Double_Puiseux_Structures;
    begin
      if isbinhom then
        hdg := Real_Powered_Homotopy.Supports(p.all,vrblvl-1);
        hcf := Extract_Leading_Coefficients(coeffs,hdg,vrblvl-1);
        hct := Extract_Leading_Powers(powers,hdg,vrblvl-1);
        Normalize_Binomial_Homotopy(hdg,hcf,hct,vrblvl);
      else
        Product_Coefficients_Powers(p,coeffs,powers,hdg,hcf,hct,vrblvl);
      end if;
      if vrblvl > 0 then
        put_line("coefficients, powers, supports :");
        for i in hcf'range loop
          put("polynomial "); put(i,1); put_line(" :");
          for j in hcf(i)'range loop
            put(hcf(i)(j)); put("  ");
            put(hct(i)(j)); put("  ");
            Standard_Integer_Vectors_IO.put(hdg(i)(j)); new_line;
          end loop;
        end loop;
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
      Double_Puiseux_Structures.Indexing_Series(powers,coeffs,vrblvl);
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
