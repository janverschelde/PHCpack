with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Matrices;
with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Symbolic_Minor_Equations;           use Symbolic_Minor_Equations;

procedure ts_canocurv is

-- DESCRIPTION :
--   Test for localization patterns of q-curves.

  procedure Write_Separator ( p,k : in natural32 ) is

  -- DESCRIPTION :
  --   Writes a separating bar between blocks in a stretched representation
  --   of a q-curve that produces p-planes.
  --   The elements in the columns occupy each k spaces.

  begin
    for j in 1..p loop
      for i in 1..k loop
        put("-");
      end loop;
    end loop;
    put_line("-");
  end Write_Separator;

  procedure Write_Pattern ( m,p,q : in natural32; top,bottom : in Bracket;
                            locpat : in Standard_Natural_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the pattern for top and bottom pivots of a q-curve that
  --   produces p-planes into (m+p)-dimensional space.

  begin
    put("The pattern of the "); put(q,1); put("-curve into G(");
    put(m,1); put(","); put(m+p,1); put(") for (");
    put(top); put(","); put(bottom); put_line(") :");
    Write_Separator(p,2);
    for i in locpat'range(1) loop
      for j in locpat'range(2) loop
        if locpat(i,j) = 0
         then put(" 0");
         else put(" *");
        end if;
      end loop;
      new_line;
      if natural32(i) mod (m+p) = 0
       then Write_Separator(p,2);
      end if;
    end loop;
  end Write_Pattern;

  procedure Write_Variables ( m,p,q : in natural32; top,bottom : in Bracket;
                              locpat : in Standard_Natural_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the variables in the pattern for top and bottom pivots of a
  --   q-curve that produces p-planes into (m+p)-dimensional space.

    row,deg : natural32 := 0;

  begin
    put("Variable pattern of the "); put(q,1); put("-curve into G(");
    put(m,1); put(","); put(m+p,1); put(") for (");
    put(top); put(","); put(bottom); put_line(") :");
    Write_Separator(p,6);
    for i in locpat'range(1) loop
      row := row + 1;
      for j in locpat'range(2) loop
        put(" ");
        if locpat(i,j) = 0 then
          for k in 1..4 loop
            put(" ");
          end loop;
          put("0");
        else
          put("x"); put(row,1); put(j,1); put("s"); put(deg,1);
        end if;
      end loop;
      new_line;
      if natural32(i) mod (m+p) = 0 then
        Write_Separator(p,6);
        row := 0; deg := deg+1;
      end if;
    end loop;
  end Write_Variables;

  procedure Write_Localization_Pattern
              ( m,p,q : in natural32; top,bottom : in Bracket ) is

  -- DESCRIPTION :
  --   Writes the variables in the pattern for top and bottom pivots of a
  --   q-curve that produces p-planes into (m+p)-dimensional space.

    rws : constant natural32 := (m+p)*(q+1);
    row,deg : natural32 := 0;

  begin
    put("Variable pattern of the "); put(q,1); put("-curve into G(");
    put(m,1); put(","); put(m+p,1); put(") for (");
    put(top); put(","); put(bottom); put_line(") :");
    Write_Separator(p,6);
    for i in 1..rws loop
      row := row + 1;
      for j in 1..integer32(p) loop
        put(" ");
        if (i < top(j) or i > bottom(j)) then
          for k in 1..4 loop
            put(" ");
          end loop;
          put("0");
        else
          put("x"); put(row,1); put(j,1); put("s"); put(deg,1);
        end if;
      end loop;
      new_line;
      if i mod (m+p) = 0 then
        Write_Separator(p,6);
        row := 0; deg := deg+1;
      end if;
    end loop;
  end Write_Localization_Pattern;

  function Number_of_Variables ( top,bottom : Bracket ) return natural32 is

  -- DESCRIPTION :
  --   Counts the number of x_ij-variables needed to represent the pattern
  --   prescribed with top and bottom pivots.
  --   Note that the actual dimension of the corresponding space is p less.

    cnt : natural32 := 0;

  begin
    for j in top'range loop 
      cnt := cnt + (bottom(j) - top(j) + 1);
    end loop;
    return cnt;
  end Number_of_Variables;

  procedure Set_up_Symbol_Table
              ( m,p,q : natural32; top,bottom : in Bracket ) is

  -- DESCRIPTION :
  --   Fills the symbol table with those symbols needed to represent a
  --   q-curve into the Grassmannian of p-planes into (m+p)-space,
  --   represented in the pattern prescribed by top and bottom pivots.
  --   The "s" and the "t" variables are the first ones that occur.

    cnt : constant natural32 := 2 + Number_of_Variables(top,bottom);
    rws : constant natural32 := (m+p)*(q+1);
    row,deg : natural32;
    sb : Symbol;

  begin
    Symbol_Table.Init(cnt);           -- initialization with #variables
    sb := (sb'range => ' ');
    sb(1) := 's';
    Symbol_Table.Add(sb);             -- adding "s"
    sb(1) := 't';
    Symbol_Table.Add(sb);             -- adding "t"
    sb(1) := 'x';
    sb(4) := 's';
    for j in 1..p loop                -- adding the rest columnwise
      row := 0; deg := 0;
      for i in 1..rws loop
        row := row + 1;
        if i >= top(integer32(j)) and i <= bottom(integer32(j)) then
          sb(2) := Convert_Hexadecimal(row);
          sb(3) := Convert_Hexadecimal(j);
          sb(5) := Convert_Hexadecimal(deg);
          Symbol_Table.Add(sb);
        end if;
        if i mod (m+p) = 0
         then row := 0; deg := deg+1;
        end if;
      end loop;
    end loop;
  end Set_up_Symbol_Table;

  procedure Write_Symbol_Table is

  -- DESCRIPTION :
  --   Writes the content of the symbol table.

  begin
    put_line("The symbols currently in the symbol table : ");
    for i in 1..Symbol_Table.Number loop
      Symbol_Table_io.put(Symbol_Table.Get(i)); new_line;
    end loop;
  end Write_Symbol_Table;

  function Polynomial_Pattern ( m,p,q : natural32; top,bottom : Bracket )
                              return Standard_Complex_Poly_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the representation of a q-curve that produces p-planes in 
  --   (m+p)-dimensional space, in the localization pattern prescribed by
  --   top and bottom pivots.  The columns of the matrix on return contain
  --   the polynomial functions that generate the p-plane.

    res : Standard_Complex_Poly_Matrices.Matrix
            (1..integer32(m+p),1..integer32(p));
    rws : constant natural32 := (m+p)*(q+1);
    n : constant natural32 := 2 + Number_of_Variables(top,bottom);
    s_deg,t_deg : natural32;
    row,ind : integer32;
    t : Term;

  begin
    for i in res'range(1) loop                    -- initialization
      for j in res'range(2) loop
        res(i,j) := Null_Poly;
      end loop;
    end loop;
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := Create(1.0);
    ind := 2;                                     -- ind counts #variables
    for j in 1..integer32(p) loop                 -- assign columnwise
      t_deg := (bottom(j)-1)/(m+p);               -- degree in t to homogenize
      row := 0; s_deg := 0;
      for i in 1..rws loop
        row := row + 1;
        if i >= top(j) and i <= bottom(j) then
          ind := ind+1;
          t.dg(1) := s_deg; t.dg(2) := t_deg;
          t.dg(ind) := 1;
          Add(res(row,j),t);
          t.dg(1) := 0; t.dg(2) := 0;
          t.dg(ind) := 0;
        end if;
        if i mod (m+p) = 0
         then row := 0; s_deg := s_deg+1; t_deg := t_deg-1;
        end if;
      end loop;
    end loop;
    Clear(t);
    return res;
  end Polynomial_Pattern;

  procedure Main is

    m,p,q : natural32 := 0;
    ans : character;

  begin
    put("Give m : "); get(m);
    put("Give p : "); get(p);
    put("Give q : "); get(q);
    put("  m = "); put(m,1);
    put("  p = "); put(p,1); 
    put("  q = "); put(q,1); new_line;
    loop
      declare
        top,bottom : Bracket(1..integer32(p));
      begin
        put("Give top pivots : "); get(top);
        put("Give bottom pivots : "); get(bottom);
        put("  top pivots : "); put(top);
        put("  bottom pivots : "); put(bottom); new_line;
        Write_Localization_Pattern(m,p,q,top,bottom);
        Set_up_Symbol_Table(m,p,q,top,bottom);
        Write_Symbol_Table;
        put_line("The representation as a matrix of polynomials : ");
        put(Polynomial_Pattern(m,p,q,top,bottom));
      end;
      put("Do you want patterns for other pivots ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Main;

begin
  new_line;
  put_line("Test on canonical form of a curve into the Grassmannian");
  new_line;
  Main;
end ts_canocurv;
