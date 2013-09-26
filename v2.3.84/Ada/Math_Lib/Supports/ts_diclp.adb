with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Dictionaries,Linear_Programming;    use Dictionaries,Linear_Programming;

procedure ts_diclp is

-- DESCRIPTION :
--   This procedure can be used for solving the primal problem:
--
--      max <c,x>
--          <a,x> <= 0
--   
--   where x = (1,x1,x2,..,xn)
--
--   and for solving the dual problem :
--
--      min <c,x>
--          <a,x> >= 0
--
--   where x = (1,x1,x2,..,xn).

 
  n,m : integer32;
  ans : character;
  eps : constant double_float := 10.0**(-12);

-- OUTPUT ROUTINES :

  procedure Write ( file : in file_type; dic : in Matrix;
                    in_bas,out_bas : in Standard_Integer_Vectors.Vector;
                    fo,af,ex : in natural32 ) is

  -- DESCRIPTION :
  --   This procedure writes the dictionary on standard output or on file.

  -- ON ENTRY :
  --   dic          matrix for the dictionary;
  --   in_bas       unknowns in the basis;
  --   out_bas      unknowns out the basis;
  --   fo           the number of decimal literals before the dot;
  --   af           the number of decimal literals after the dot;
  --   ex           the number of decimal literals in the exponent.

  begin
    new_line(file);
    put_line(file,"*****  THE DICTIONARY  *****");
    new_line(file);
    put_line(file,"The elements in the basis : ");
    for i in in_bas'range loop
      put(file," "); put(file,in_bas(i),1);
    end loop;
    new_line(file);
    put_line(file,"The dictionary : ");
    for i in dic'range(1) loop
      for j in dic'range(2) loop
        put(file,dic(i,j), fore => fo, aft => af, exp => ex );
      end loop;
      new_line(file);
    end loop;
  end Write;

  procedure Report ( dic : in Matrix;
                     in_bas,out_bas : in Standard_Integer_Vectors.Vector ) is
  begin
    Write(Standard_Output,dic,in_bas,out_bas,3,3,3);
  end Report;

  procedure plp is new Generic_Primal_Simplex(Report);
  procedure dlp is new Generic_Dual_Simplex(Report);

begin

-- INPUT OF THE DATA :

  put("Give the number of unknowns : "); get(n);
  put("Give the number of constraints : "); get(m);

  declare
    c : Standard_Floating_Vectors.vector(0..n);
    a : matrix(1..m,0..n);
    dic : matrix(0..m,0..n);
    inbas : Standard_Integer_Vectors.Vector(1..m) := (1..m => 0);
    outbas : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0);
    primsol : Standard_Floating_Vectors.Vector(1..n);
    dualsol : Standard_Floating_Vectors.Vector(1..m);
    nit : natural32 := 0;
    ok : boolean := false;
    feasibound : boolean;
  begin

    put_line("Give the coefficients of the target :");
    for i in c'range loop
      get(c(i));
    end loop;
    put_line("Give the coefficients of the constraints :");
    for i in a'range(1) loop
      for j in a'range(2) loop
	get(a(i,j));
      end loop;
    end loop;
    loop
      put("Primal of Dual simplex problem ? (p/d) "); get(ans);
      if ans = 'p'
       then Dictionaries.Primal_Init(c,a,dic,inbas,outbas); ok := true;
       elsif ans = 'd'
  	   then Dictionaries.Dual_Init(c,a,dic,inbas,outbas); ok := true;
  	   else put_line("try again ...");  ok := false;
      end if;
      exit when ok;
    end loop;

    if ans = 'p'
     then plp(dic,eps,inbas,outbas,nit,feasibound);
     else dlp(dic,eps,inbas,outbas,nit,feasibound);
    end if;

    primsol := Dictionaries.Primal_Solution(dic,inbas,outbas);
    new_line;
    put_line("THE PRIMAL SOLUTION :");
    for i in primsol'range loop
      put("  x"); put(i,1); put(" : ");
      put(primsol(i),3,3,3); new_line;
    end loop;

    if ans = 'd'
     then dualsol := -Dictionaries.Dual_Solution(dic,inbas,outbas);
     else dualsol :=  Dictionaries.Dual_Solution(dic,inbas,outbas);
    end if;
    new_line;
    put_line("THE DUAL SOLUTION :");
    for i in dualsol'range loop
      put("  x"); put(i,1); put(" : ");
      put(dualsol(i),3,3,3); new_line;
    end loop;

    new_line;
    put("OPTIMUM : "); put(Optimum(dic)); new_line;
    put("NUMBER OF ITERATIONS : "); put(nit,1); new_line;
    if ans = 'd'
     then if feasibound
           then put_line("THE PROBLEM IS FEASIBLE");
           else put_line("THE PROBLEM IS NOT FEASIBLE");
          end if;
     else if feasibound
           then put_line("THE PROBLEM IS UNBOUNDED");
           else put_line("THE PROBLEM IS NOT UNBOUNDED");
          end if;
    end if;

  end;

end ts_diclp;
