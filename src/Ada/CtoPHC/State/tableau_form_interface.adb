with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;

package body Tableau_Form_Interface is

-- AUXILIARY FUNCTIONS :

  procedure Extract_Dimensions
              ( a : in C_intarrs.Pointer;
                neq,nvr : out integer32 ) is

  -- DESCRIPTION :
  --   Extracts the dimensions out of the two numbers of a.

  -- ON ENTRY :
  --   a        in a[0] is the number of equations,
  --            in a[1] is the number of variables.

  -- ON RETURN :
  --   neq      the number of equations,
  --   nvr      the number of variables.

    v : constant C_Integer_Array
      := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));

    use Interfaces.C;

  begin
    neq := integer32(v(v'first));
    nvr := integer32(v(v'first+1));
  end Extract_Dimensions;

  --procedure Extract_Dimensions
  --            ( a : in C_intarrs.Pointer;
  --              neq,nvr,nbt : out integer32; verbose : out boolean ) is

  -- DESCRIPTION :
  --   Extracts the dimensions out of the four numbers of a.

  -- ON ENTRY :
  --   a        in a[0] is the number of equations,
  --            in a[1] is the number of variables,
  --            in a[2] is the total number of terms,
  --            in a[3] is the value for the verbose flag.

  -- ON RETURN :
  --   neq      the number of equations;
  --   nvr      the number of variables;
  --   nbt      total number of terms;
  --   verbose  is the verbose flag.

  --  v : constant C_Integer_Array
  --    := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(4));

  --  use Interfaces.C;

  --begin
  --  neq := integer32(v(v'first));
  --  nvr := integer32(v(v'first+1));
  --  nbt := integer32(v(v'first+2));
  --  verbose := (integer32(v(v'first+3)) = 1);
  --end Extract_Dimensions;

  procedure Write_Tableau
              ( neq,nvr : in integer32;
                nbterms : in Standard_Integer_Vectors.Vector;
                coefficients : in Standard_Complex_Vectors.Vector;
                exponents : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Writes the tableau form of a polynomial system.

  -- ON ENTRY :
  --   neq      the number of equations;
  --   nvr      the number of variables;
  --   nbterms  nbterms(k) is the number of terms in the k-th polynomial;
  --   exponents are the exponents of the monomials;
  --   coefficients are the coefficients of the terms.

    cffidx : integer32 := coefficients'first-1;
    expidx : integer32 := exponents'first-1;

  begin
    put(neq,1);
    if neq /= nvr
     then put("  "); put(nvr,1);
    end if;
    new_line;
    for i in nbterms'range loop
      put(nbterms(i),1); new_line;
      for j in 1..nbterms(i) loop
        cffidx := cffidx + 1;
        put(coefficients(cffidx));
        for k in 1..nvr loop
          expidx := expidx + 1;
          put(" "); put(exponents(expidx),1);
        end loop;
        new_line;
      end loop;
    end loop;
  end Write_Tableau;

  function Make_System
             ( neq,nvr : integer32;
               nbterms : Standard_Integer_Vectors.Vector;
               coefficients : Standard_Complex_Vectors.Vector;
               exponents : Standard_Integer_Vectors.Vector;
               verbose : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the polynomial system that corresponds
  --   to the give tableau form.

  -- ON ENTRY :
  --   neq      the number of equations;
  --   nvr      the number of variables;
  --   nbterms  nbterms(k) is the number of terms in the k-th polynomial;
  --   exponents are the exponents of the monomials;
  --   coefficients are the coefficients of the terms.

    res : Standard_Complex_Poly_Systems.Poly_Sys(1..neq);
    trm : Standard_Complex_Polynomials.Term;
    cffidx : integer32 := coefficients'first-1;
    expidx : integer32 := exponents'first-1;

  begin
    if verbose then
      put("-> the number of equations : "); put(neq,1); new_line;
      put("-> the number of variables : "); put(nvr,1); new_line;
    end if;
    trm.dg := new Standard_Natural_Vectors.Vector(1..nvr);
    for i in nbterms'range loop
      res(i) := Standard_Complex_Polynomials.Null_Poly;
      if verbose then
        put("-> the number of terms of polynomial ");
        put(i,1); put(" : "); put(nbterms(i),1); new_line;
      end if;
      for j in 1..nbterms(i) loop
        cffidx := cffidx + 1;
        trm.cf := coefficients(cffidx);
        if verbose
         then put(coefficients(cffidx));
        end if;
        for k in 1..nvr loop
          expidx := expidx + 1;
          trm.dg(k) := natural32(exponents(expidx));
          if verbose
           then put(" "); put(exponents(expidx),1);
          end if;
        end loop;
        Standard_Complex_Polynomials.Add(res(i),trm);
        if verbose
         then new_line;
        end if;
      end loop;
    end loop;
    return res;
  end Make_System;

  procedure Extract_Tableau
              ( b : in C_intarrs.Pointer;
                c : in C_dblarrs.Pointer;
                neq,nvr,nbt : in integer32;
                nbterms : in Standard_Integer_Vectors.Vector;
                verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Extracts the coefficients and exponents from b and c.

    nbrexp : constant integer32 := nvr*nbt;
    nbrcff : constant integer32 := 2*nbt;
    expdata : Standard_Integer_Vectors.Vector(1..nbrexp);
    cffdata : Standard_Complex_Vectors.Vector(1..nbrcff);
    p : Standard_Complex_Poly_Systems.Poly_Sys(1..neq);

  begin
    Assign(natural32(nbrexp),b,expdata);
    Assign(natural32(nbrcff),c,cffdata);
    if verbose
     then Write_Tableau(neq,nvr,nbterms,cffdata,expdata);
    end if;
    p := Make_System(neq,nvr,nbterms,cffdata,expdata,verbose);
    if verbose then
      put_line("The polynomial system made from the tableau form :");
      put(p);
    end if;
    Standard_PolySys_Container.Initialize(p);
  end Extract_Tableau;

-- EXPORTED FUNCTIONS :

  function Tableau_Form_Store
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    neq,nvr,nbt : integer32;
    verbose : boolean;

  begin
    if vrblvl > 0 then
      put_line("-> in tableau_form_interface.Tableau_Form_Store ...");
    end if;
    Extract_Dimensions(a,neq,nvr);
    declare
      dim : constant integer32 := 3+neq+1;
      dimdata : Standard_Integer_Vectors.Vector(1..dim);
      nbterms : Standard_Integer_Vectors.Vector(1..neq);
    begin
      Assign(natural32(dim),a,dimdata);
      verbose := (dimdata(dim) > 0);
      if verbose then
        put("-> the number of equations : "); put(neq,1); new_line;
        put("-> the number of variables : "); put(nvr,1); new_line;
      end if;
      for i in 1..neq loop
        nbterms(i) := dimdata(i+3);
      end loop;      
      nbt := Standard_Integer_Vectors.Sum(nbterms);
      if verbose then
        put("-> total number of terms : "); put(nbt,1); new_line;
        for i in 1..neq loop
          put("-> number of terms in polynomial ");
          put(i,1); put(" : "); put(nbterms(i),1); new_line;
        end loop;
      end if;
      Extract_Tableau(b,c,neq,nvr,nbt,nbterms,verbose);
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in tableau_form_interface.");
        put_line("Tableau_Form_Store.");
      end if;
      return 889;
  end Tableau_Form_Store;

  function Tableau_Form_Dimensions
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    neq,nvr,nbt : integer32 := 0;
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    res : Standard_Integer_Vectors.Vector(1..3);

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

  begin
    if vrblvl > 0 then
      put_line("-> in tableau_form_interface.Tableau_Form_Dimensions ...");
    end if;
    if lp /= null then
      neq := lp'last;
      nvr := integer32(Number_of_Unknowns(lp(lp'first)));
      for k in 1..neq loop
        nbt := nbt + integer32(Standard_PolySys_Container.Number_of_Terms(k));
      end loop;
    end if;
    res(1) := neq;
    res(2) := nvr;
    res(3) := nbt;
    Assign(res,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in tableau_form_interface.");
        put_line("Tableau_Form_Dimension.");
      end if;
      return 890;
  end Tableau_Form_Dimensions;

end Tableau_Form_Interface;
