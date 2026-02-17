with Ada.Text_IO;                       use Ada.Text_IO;
with String_Splitters;
with Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Complex_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_IO;      use Standard_Complex_Matrices_IO;
with Standard_Numerical_Rank;
with Symbol_Table;
with Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;
with Standard_Random_Laurentials;       use Standard_Random_Laurentials;
with Standard_Complex_Laur_Systems;
with Real_Powered_Series_IO;
with Random_Laurent_Homotopy;
with Test_Real_Powered_Series;
with Real_Powered_Homotopy;
with Real_Powered_Homotopy_io;

package body Test_Real_Powered_Homotopy is

  procedure Random_Series_Vector
              ( size : in integer32;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwt : out Standard_Floating_VecVecs.VecVec ) is
  begin
    for i in cff'range loop
      declare
        icf : Standard_Complex_Vectors.Vector(0..size);
        ipw : Standard_Floating_Vectors.Vector(1..size);
      begin
        Test_Real_Powered_Series.Random_Series(size,icf,ipw);
        cff(i) := new Standard_Complex_Vectors.Vector'(icf);
        pwt(i) := new Standard_Floating_Vectors.Vector'(ipw);
      end;
    end loop;
  end Random_Series_Vector;

  procedure Test_Random_Polynomial ( nbr,nvr,size : in integer32 ) is

    q : Standard_Complex_Laurentials.Poly;
    cff : Standard_Complex_VecVecs.VecVec(1..nbr);
    pwt : Standard_Floating_VecVecs.VecVec(1..nbr);
    ans : character;
    file : file_type;
    name : String_Splitters.Link_to_String;
    nsyms : integer32;
  
  begin
    new_line;
    put("-> generating "); put(nbr,1); put(" monomials in ");
    put(nvr,1); put_line(" variables ...");
    q := Random_Laurent_Polynomial(natural32(nvr),natural32(nbr),-9,+9);
    Standard_Complex_Laurentials_io.put(q); new_line;
    Random_Series_Vector(size,cff,pwt);
    for i in cff'range loop
      put("random series "); put(i,1); put_line(" :");
      Real_Powered_Series_IO.put_line(cff(i).all,pwt(i).all);
    end loop;
    new_line;
    put("new line format ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then Real_Powered_Homotopy_IO.put_line(q,cff,pwt);
     else Real_Powered_Homotopy_IO.put(q,cff,pwt);
    end if;
    new_line;
    put_line("Reading file name for output ...");
    Communications_with_User.Read_Name_and_Create_File(file,name);
    if ans = 'y'
     then Real_Powered_Homotopy_IO.put_line(file,q,cff,pwt);
     else Real_Powered_Homotopy_IO.put(file,q,cff,pwt);
    end if;
    Close(file);
    new_line;
    put_line("Closed file.  Reopening again for reading ...");
    new_line;
    Communications_with_User.Open_Input_File(file,name.all);
    nsyms := integer32(Symbol_Table.number);
    put("Number of symbols : "); put(nsyms,1); new_line;
    if nsyms = 0 
     then Symbol_Table.Init(natural32(nvr));
    end if;
    declare
      q2 : Standard_Complex_Laurentials.Poly;
      c2 : Standard_Complex_VecVecs.VecVec(1..nbr);
      p2 : Standard_Floating_VecVecs.VecVec(1..nbr);
    begin
      Real_Powered_Homotopy_IO.get(file,nvr,size,q2,c2,p2,vrblvl=>2);
      for i in 1..nbr loop
        put("-> power series "); put(i,1); put_line(" :");
        Real_Powered_Series_IO.put_line(cff(i).all,pwt(i).all);
      end loop;
      put_line("Laurent polynomial read from file :");
      Standard_Complex_Laurentials_io.put(q2); new_line;
      new_line;
      put_line("-> the real powered Laurent homotopy polynomial :");
      if ans = 'y'
       then Real_Powered_Homotopy_IO.put_line(q,cff,pwt);
       else Real_Powered_Homotopy_IO.put(q,cff,pwt);
      end if;
    end;
  end Test_Random_Polynomial;

  procedure Test_String_Polynomial ( nbr,nvr,size : in integer32 ) is

    q : Standard_Complex_Laurentials.Poly;
    cff : Standard_Complex_VecVecs.VecVec(1..nbr);
    pwt : Standard_Floating_VecVecs.VecVec(1..nbr);
    nsyms : constant integer32 := integer32(Symbol_Table.number);
  
  begin
    new_line;
    put("-> generating "); put(nbr,1); put(" monomials in ");
    put(nvr,1); put_line(" variables ...");
    q := Random_Laurent_Polynomial(natural32(nvr),natural32(nbr),-9,+9);
    Standard_Complex_Laurentials_io.put(q); new_line;
    Random_Series_Vector(size,cff,pwt);
    for i in cff'range loop
      put("random series "); put(i,1); put_line(" :");
      Real_Powered_Series_IO.put_line(cff(i).all,pwt(i).all);
    end loop;
    put("Number of symbols : "); put(nsyms,1); new_line;
    if nsyms = 0 then
      declare
        asb : constant Symbol_Table.Array_of_Symbols(1..nvr)
            := Symbol_Table.Standard_Symbols(nvr);
      begin
        Symbol_Table.Init(asb);
      end;
    end if;
    declare
      s : constant string
        := Real_Powered_Homotopy_IO.to_string(q,cff,pwt,vrblvl=>1);
    begin
      new_line;
      put_line("The string representation of the polynomial :");
      put_line(s);
      new_line;
      put_line("Parsing the string for the data ...");
      declare
        m : constant integer32
          := Real_Powered_Homotopy_IO.number_of_terms(s,2);
      begin
        put("number of terms in string representation : ");
        put(m,1);
        if m /= nbr then
          put(" /= "); put(nbr,1); put_line(", bug!");
        else
          put_line(", okay.");
          declare
            q2 : Standard_Complex_Laurentials.Poly;
            c2 : Standard_Complex_VecVecs.VecVec(1..nbr);
            p2 : Standard_Floating_VecVecs.VecVec(1..nbr);
          begin
            Real_Powered_Homotopy_IO.parse_string(s,nvr,q2,c2,p2,vrblvl=>2);
            for i in cff'range loop
              put("parsed series "); put(i,1); put_line(" :");
              Real_Powered_Series_IO.put_line(c2(i).all,p2(i).all);
            end loop;
            put_line("parsed Laurent polynomial :");
            Standard_Complex_Laurentials_io.put(q2); new_line;
          end;
        end if;
      end;
    end;
  end Test_String_Polynomial;

  procedure Random_System
             ( dim,low,upp,size : in integer32;
               mbn : in Standard_Integer_Vectors.Vector;
               intpow : in boolean;
               eqs : out Standard_Complex_Laur_Systems.Laur_Sys;
               cff : out Standard_Complex_VecVecs.Array_of_VecVecs;
               pwt : out Standard_Floating_VecVecs.Array_of_VecVecs ) is

  -- DESCRIPTION :
  --   Generates a random Laurent system of dimension dim.

  -- ON ENTRY :
  --   dim     number of equations and variables;
  --   low     lower bound on the exponents of the monomials;
  --   upp     upper bound on the exponents of the monomials;
  --   size    size of the power series coefficients;
  --   mbn     mbn(i) equals the number of monomials in polynomial i;
  --   intpow  flag for integer powers of the power series.

  -- ON RETURN :
  --   eqs     random Laurent polynomials;
  --   cff     coefficients of each power series;
  --   pwt     powers of the each power series.

    deg : Standard_Integer_VecVecs.Array_of_VecVecs(1..dim);

  begin
    Random_Laurent_Homotopy.Random_Laurent_System
      (dim,dim,low,upp,size,mbn,deg,cff,pwt,intpow);
    for i in deg'range loop
      eqs(i) := Standard_Complex_Laurentials.Null_Poly;
      for j in deg(i)'range loop
        declare
          trm : Standard_Complex_Laurentials.Term;
        begin
          trm.cf := Standard_Complex_Numbers.Create(1.0);
          trm.dg := Standard_Complex_Laurentials.Degrees(deg(i)(j));
          Standard_Complex_Laurentials.Add(eqs(i),trm);
        end;
      end loop;
    end loop;
  end Random_System;

  procedure Test_Random_System ( dim,nbr,size : in integer32 ) is

    q : Standard_Complex_Laur_Systems.Laur_Sys(1..dim);
    c : Standard_Complex_VecVecs.Array_of_VecVecs(1..dim);
    p : Standard_Floating_VecVecs.Array_of_VecVecs(1..dim);
    m : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => nbr);
    ans : character;
    integer_powers : boolean;
    file : file_type;
    name : String_Splitters.Link_to_String;

  begin
    put("Integer powers of the series ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    integer_powers := (ans = 'y');
    Random_System(dim,-9,9,size,m,integer_powers,q,c,p);
    new_line;
    put("new line format ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    if ans = 'y'
     then Real_Powered_Homotopy_IO.put_line(dim,dim,size,q,c,p);
     else Real_Powered_Homotopy_IO.put(dim,dim,size,q,c,p);
    end if;
    new_line;
    put_line("Reading file name for output ...");
    Communications_with_User.Read_Name_and_Create_File(file,name);
    if ans = 'y'
     then Real_Powered_Homotopy_IO.put_line(file,dim,dim,size,q,c,p);
     else Real_Powered_Homotopy_IO.put(file,dim,dim,size,q,c,p);
    end if;
    Close(file);
    new_line;
    put_line("Closed file.  Reopening again for reading ...");
    new_line;
    Communications_with_User.Open_Input_File(file,name.all);
    declare
      lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
      lc : Standard_Complex_VecVecs.Link_to_Array_of_VecVecs;
      lp : Standard_Floating_VecVecs.Link_to_Array_of_VecVecs;
    begin
      Real_Powered_Homotopy_IO.get(file,lq,lc,lp,vrblvl=>2);
      new_line;
      put_line("The system read from file :");
      Real_Powered_Homotopy_IO.put_line(dim,dim,size,lq.all,lc.all,lp.all);
    end;
  end Test_Random_System;

  procedure Test_Regularity ( dim,nbr,size : in integer32 ) is

    q : Standard_Complex_Laur_Systems.Laur_Sys(1..dim);
    c : Standard_Complex_VecVecs.Array_of_VecVecs(1..dim);
    p : Standard_Floating_VecVecs.Array_of_VecVecs(1..dim);
    m : constant Standard_Integer_Vectors.Vector(1..dim) := (1..dim => nbr);
    ans : character;
    integer_powers : boolean;
    A : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    b : Standard_Complex_Vectors.Vector(1..dim);
    tol : constant double_float := 1.0E-14;
    rnk : integer32;

  begin
    put("Integer powers of the series ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    integer_powers := (ans = 'y');
    Random_System(dim,0,1,size,m,integer_powers,q,c,p);
    new_line;
    Real_Powered_Homotopy_IO.put_line(dim,dim,size,q,c,p);
    Real_Powered_Homotopy.Get_Constant_Coefficients(q,c,A,b,1);
    put_line("The coefficient matrix :"); put(A);
    put_line("The right hand side vector :"); put_line(b);
    rnk := integer32(Standard_Numerical_Rank.Numerical_Rank(A,tol));
    put("-> the numerical rank : "); put(rnk,1); new_line;
    if rnk = dim
     then put_line("The initial solution is regular.");
     else put_line("The initial solution is singular.");
    end if;
  end Test_Regularity;

  procedure Main is

    nbr,nvr,size : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the number of variables : "); get(nvr);
    put("Give the number of monomials : "); get(nbr);
    put("Give the size of each series : "); get(size);
    new_line;
    put_line("MENU for testing real powered homotopies :");
    put_line("  1. test string output/input");
    put_line("  2. test output/input of one Laurent polynomial");
    put_line("  3. test output/input of a Laurent system");
    put_line("  4. test regularity of a random system");
    put("Type 1, 2, 3, or 4 to select a test : ");
    Communications_with_User.Ask_Alternative(ans,"1234");
    new_line;
    case ans is
      when '1' => Test_String_Polynomial(nbr,nvr,size);
      when '2' => Test_Random_Polynomial(nbr,nvr,size);
      when '3' => Test_Random_System(nvr,nbr,size);
      when '4' => Test_Regularity(nvr,nbr,size);
      when others => null;
    end case;
  end Main;

end Test_Real_Powered_Homotopy;
