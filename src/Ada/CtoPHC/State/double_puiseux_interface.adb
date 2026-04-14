with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Integer_VecVecs;
with Standard_Integer_VecVecs_IO;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
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

  procedure Run_Newton_Steps
              ( nbr : in integer32; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Extracts the real powered Laurent homotopy data
  --   and then runs as many Newton steps as the value of nbr.

    p : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    if vrblvl > 0
     then put_line("-> in double_puiseux_interface.Run_Newton_Steps ...");
    end if;
    p := Standard_LaurSys_Container.Retrieve;
    declare
      hdg : Standard_Integer_VecVecs.Array_of_VecVecs(p'range)
          := Real_Powered_Homotopy.Supports(p.all);
    begin
      if vrblvl > 0 then
        for i in hdg'range loop
          put("support of polynomial "); put(i,1); put_line(" :");
          Standard_Integer_VecVecs_IO.put(hdg(i));
        end loop;
      end if;
    end; 
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
       -- Show_Data(vrblvl);
        Indexing_Series(vrblvl);
        Run_Newton_Steps(nbr,vrblvl);
      end if;
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
