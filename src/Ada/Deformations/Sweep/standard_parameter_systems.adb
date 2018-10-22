with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Symbol_Table,Symbol_Table_io;      use Symbol_Table;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;

package body Standard_Parameter_Systems is

  procedure Sort ( v : in out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   The labels of the parameters must be sorted in increasing order,
  --   otherwise, the substitute will not work.

    ind : integer32;
    min : integer32;

  begin
    for i in v'first..v'last-1 loop
      min := v(i);
      ind := i;
      for j in v'first+1..v'last loop
        if v(j) < min then
          ind := j;
          min := v(j);
        end if;
      end loop;
      if ind /= i then
        v(ind) := v(i);
        v(i) := min;
      end if;
    end loop;
  end Sort;

  procedure Read_Solution_Parameters
              ( infile : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 ) is

    found : boolean;

  begin
    nb_equ := integer32(p'last);
    nb_unk := integer32(Standard_Complex_Polynomials.Number_of_Unknowns(p(1)));
    nb_par := nb_unk - nb_equ;
    Scan_and_Skip(infile,"THE SOLUTIONS",found);
    if found
     then get(infile,sols);
     else Read(sols);
    end if;
  end Read_Solution_Parameters;

  procedure Read_Solution_Parameters
              ( infile : in file_type; outfile : out file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 ) is
  begin
    Read_Solution_Parameters(infile,p,sols,nb_equ,nb_unk,nb_par);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    put(outfile,natural32(nb_equ),natural32(nb_unk),p);
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS : ");
    put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Read_Solution_Parameters;

  procedure Read_Parameter_Homotopy
              ( lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the file name for a polynomial system.");
    Read_Name_and_Open_File(infile);
    get(infile,lp);
    Read_Solution_Parameters(infile,lp.all,sols,nb_equ,nb_unk,nb_par);
  end Read_Parameter_Homotopy;

  procedure Read_Parameter_Homotopy
              ( outfile : out file_type;
                lp : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                sols : out Solution_List;
                nb_equ,nb_unk,nb_par : out integer32 ) is

    infile : file_type;

  begin
    new_line;
    put_line("Reading the file name for a polynomial system.");
    Read_Name_and_Open_File(infile);
    get(infile,lp);
    Read_Solution_Parameters
      (infile,outfile,lp.all,sols,nb_equ,nb_unk,nb_par);
  end Read_Parameter_Homotopy;

  function Define_Parameters ( nb_equ,nb_unk,nb_par : integer32 )
                             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..nb_par);
    ind : integer32 := 0;
    ans : character;
    valid : boolean;

  begin
    new_line;
    put("  + number of unknowns   : "); put(nb_unk,1); new_line;
    put("  - number of equations  : "); put(nb_equ,1); new_line;
    put("  = number of parameters : "); put(nb_par,1); new_line;
    new_line;
    put("The "); put(nb_unk,1); put(" unknowns :");
    Symbol_Table_io.Write; new_line;
    put("We have to define ");
    if nb_par = 1
     then put_line("one parameter.");
     else put(nb_par,1); put_line(" parameters.");
    end if;
    for i in 1..nb_par loop
      declare
        sb : Symbol;
      begin
        loop
          put("Give the symbol for parameter ");
          put(i,1); put(" : ");
          Symbol_Table_io.get(sb);
          ind := integer32(Symbol_Table.Get(sb));
          if (ind < 1) or (ind > nb_unk) then
            valid := false;
            put_line("  Invalid symbol.  Please try again...");
          else
            put("> Is "); put('"');
            Symbol_Table_io.put(Symbol_Table.Get(natural32(ind)));
            put('"'); put(" going to be a parameter ? (y/n) ");
            Ask_Yes_or_No(ans);
            valid := (ans = 'y');
	  end if;
          exit when valid;
        end loop;
      end;
      res(i) := ind;
    end loop;
    Sort(res);
    return res;
  end Define_Parameters;

  function Complement ( n : integer32; v : Standard_Integer_Vectors.Vector )
                      return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(v'first..n-v'last+v'first-1);
    ind : integer32 := res'first-1;
    ind_v : integer32 := v'first;

  begin
    for i in 1..n loop
      if ((ind_v > v'last) or else (i < v(ind_v))) then
        ind := ind+1;
        res(ind) := i;
      else
        ind_v := ind_v + 1;
      end if;
    end loop;
    return res;
  end Complement;

  function Substitute
	      ( t : Standard_Complex_Polynomials.Term;
                pars : Standard_Integer_Vectors.Vector;
                vals : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Term is

    use Standard_Complex_Polynomials;

    res : Term;
    n : constant integer32 := t.dg'last - pars'last;
    ind_pars : integer32 := pars'first;
    ind_vars : integer32 := t.dg'first-1;

  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..n);
    for i in t.dg'range loop
      if ((ind_pars > pars'last)
          or else (i < pars(ind_pars))) then 
        ind_vars := ind_vars + 1;
        res.dg(ind_vars) := t.dg(i);
      else -- i = pars(ind_pars)
        for j in 1..t.dg(pars(ind_pars)) loop
          res.cf := res.cf*vals(ind_pars);
        end loop;
        ind_pars := ind_pars + 1;
      end if;
    end loop;
    return res;
  end Substitute;

  function Substitute
              ( p : Standard_Complex_Polynomials.Poly;
                pars : Standard_Integer_Vectors.Vector;
                vals : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Substitute_Term ( t : in Term; continue : out boolean ) is

      st : constant Term := Substitute(t,pars,vals);

    begin
      Add(res,st);
      continue := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator(Substitute_Term);

  begin
    Substitute_Terms(p);
    return res;
  end Substitute;

  function Substitute
               ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                 pars : Standard_Integer_Vectors.Vector;
                 vals : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Substitute(p(i),pars,vals);
    end loop;
    return res;
  end Substitute;

end Standard_Parameter_Systems; 
