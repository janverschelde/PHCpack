with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Floating_Vectors;
with Multprec_Complex_Vectors;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Multprec_Random_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;
with Standard_Durand_Kerner;
with DoblDobl_Durand_Kerner;
with QuadDobl_Durand_Kerner;
with Multprec_Durand_Kerner;
with Black_Box_Univariate_Solvers;       use Black_Box_Univariate_Solvers;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Multprec_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Multprec_Solutions_Container;

package body Univariate_Solvers_Interface is

  function Standard_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Random_Vectors;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;
    use Standard_Durand_Kerner;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    max: constant natural32 := natural32(v_a(v_a'first));
    v_c : constant C_Double_Array := C_dblarrs.Value(c);
    eps : constant double_float := double_float(v_c(v_c'first));
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in univariate_solvers_interface.");
      put_line("Standard_Univariate_Solver ...");
    end if;
    if lp = null then
      return -1;
    else
      declare
        p : constant Poly := lp(lp'first);
        deg : constant integer32 := Degree(p);
        cff : constant Standard_Complex_Vectors.Vector
            := Coefficient_Vector(natural32(deg),p);
        dcp : constant Standard_Complex_Vectors.Vector := Derivative(cff);
        rts : Standard_Complex_Vectors.Vector(1..deg) := Random_Vector(1,deg);
        rsi : Standard_Complex_Vectors.Vector(1..deg) := rts;
        err,rco,res : Standard_Floating_Vectors.Vector(rts'range);
        nb : natural32 := 0;
        fail : boolean;
      begin
        if deg = 1 then
          sols := Create_Solution_List(-cff(0));
        elsif deg > 1 then
          Silent_Durand_Kerner(cff,rts,rsi,max,eps,nb,fail);
          Newton(cff,dcp,rts,err,rco,res);
          sols := Create_Solution_List(rts,err,rco,res);
        end if;
        Assign(integer32(nb),b);
        Standard_Solutions_Container.Initialize(sols);
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in univariate_solvers_interface.");
        put_line("Standard_Univariate_Solver.");
      end if;
      return 272;
  end Standard_Univariate_Solver;

  function DoblDobl_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Random_Vectors;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Durand_Kerner;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    max: constant natural32 := natural32(v_a(v_a'first));
    v_c : constant C_Double_Array := C_dblarrs.Value(c);
    eps : constant double_float := double_float(v_c(v_c'first));
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in univariate_solvers_interface.");
      put_line("DoblDobl_Univariate_Solver ...");
    end if;
    if lp = null then
      return -1;
    else
      declare
        p : constant Poly := lp(lp'first);
        deg : constant integer32 := Degree(p);
        cff : constant DoblDobl_Complex_Vectors.Vector
            := Coefficient_Vector(natural32(deg),p);
        dcp : constant DoblDobl_Complex_Vectors.Vector := Derivative(cff);
        rts : DoblDobl_Complex_Vectors.Vector(1..deg) := Random_Vector(1,deg);
        rsi : DoblDobl_Complex_Vectors.Vector(1..deg) := rts;
        err,rco,res : Standard_Floating_Vectors.Vector(rts'range);
        nb : natural32 := 0;
        fail : boolean;
      begin
        if deg = 1 then
          sols := Create_Solution_List(-cff(0));
        elsif deg > 1 then
          Silent_Durand_Kerner(cff,rts,rsi,max,eps,nb,fail);
          Newton(cff,dcp,rts,err,rco,res);
          sols := Create_Solution_List(rts,err,rco,res);
        end if;
        Assign(integer32(nb),b);
        DoblDobl_Solutions_Container.Initialize(sols);
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in univariate_solvers_interface.");
        put_line("DoblDobl_Univariate_Solver.");
      end if;
      return 273;
  end DoblDobl_Univariate_Solver;

  function QuadDobl_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Random_Vectors;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Durand_Kerner;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    max: constant natural32 := natural32(v_a(v_a'first));
    v_c : constant C_Double_Array := C_dblarrs.Value(c);
    eps : constant double_float := double_float(v_c(v_c'first));
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in univariate_solvers_interface.");
      put_line("QuadDobl_Univariate_Solver ...");
    end if;
    if lp = null then
      return -1;
    else
      declare
        p : constant Poly := lp(lp'first);
        deg : constant integer32 := Degree(p);
        cff : constant QuadDobl_Complex_Vectors.Vector
            := Coefficient_Vector(natural32(deg),p);
        dcp : constant QuadDobl_Complex_Vectors.Vector := Derivative(cff);
        rts : QuadDobl_Complex_Vectors.Vector(1..deg) := Random_Vector(1,deg);
        rsi : QuadDobl_Complex_Vectors.Vector(1..deg) := rts;
        err,rco,res : Standard_Floating_Vectors.Vector(rts'range);
        nb : natural32 := 0;
        fail : boolean;
      begin
        if deg = 1 then
          sols := Create_Solution_List(-cff(0));
        elsif deg > 1 then
          Silent_Durand_Kerner(cff,rts,rsi,max,eps,nb,fail);
          Newton(cff,dcp,rts,err,rco,res);
          sols := Create_Solution_List(rts,err,rco,res);
        end if;
        Assign(integer32(nb),b);
        QuadDobl_Solutions_Container.Initialize(sols);
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in univariate_solvers_interface.");
        put_line("QuadDobl_Univariate_Solver.");
      end if;
      return 274;
  end QuadDobl_Univariate_Solver;

  function Multprec_Univariate_Solver
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;
    use Multprec_Floating_Numbers;
    use Multprec_Complex_Numbers;
    use Multprec_Complex_Vectors;
    use Multprec_Random_Vectors;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Solutions;
    use Multprec_Durand_Kerner;

    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    deci : constant natural32 := natural32(v_a(v_a'first));
    max : constant natural32 := natural32(v_a(v_a'first+1));
    size : constant natural32 := Decimal_to_Size(deci);
    v_c : constant C_Double_Array := C_dblarrs.Value(c);
    d_eps : constant double_float := double_float(v_c(v_c'first));
    sols : Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in univariate_solvers_interface.");
      put_line("Multprec_Univariate_Solver ...");
    end if;
    if lp = null then
      return -1;
    else
      declare
        p : constant Poly := lp(lp'first);
        deg : constant integer32 := Degree(p);
        cff : constant Multprec_Complex_Vectors.Vector
            := Coefficient_Vector(natural32(deg),p);
        dcp : constant Multprec_Complex_Vectors.Vector := Derivative(cff);
        rts : Multprec_Complex_Vectors.Vector(1..deg)
            := Random_Vector(1,deg,size);
        rsi : Multprec_Complex_Vectors.Vector(1..deg);
        eps : Floating_Number := create(d_eps);
        err,rco,res : Multprec_Floating_Vectors.Vector(rts'range);
        nb : natural32 := 0;
        fail : boolean;
      begin
        if deg = 1 then
          sols := Create_Solution_List(-cff(0));
        elsif deg > 1 then
          Multprec_Complex_Vectors.Copy(rts,rsi);
          Silent_Durand_Kerner(cff,rts,rsi,max,eps,nb,fail);
          Newton(cff,dcp,rts,err,rco,res);
          sols := Create_Solution_List(rts,err,rco,res);
          Multprec_Complex_Vectors.Clear(rsi);
          Multprec_Floating_Vectors.Clear(err);
          Multprec_Floating_Vectors.Clear(rco);
          Multprec_Floating_Vectors.Clear(res);
        end if;
        Assign(integer32(nb),b);
        Multprec_Solutions_Container.Initialize(sols);
        Multprec_Floating_Numbers.Clear(eps);
      end;
      return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in univariate_solvers_interface.");
        put_line("Multprec_Univariate_Solver.");
      end if;
      return 275;
  end Multprec_Univariate_Solver;

end Univariate_Solvers_Interface;
