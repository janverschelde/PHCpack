with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Standard_Integer64_Linear_Solvers;  use Standard_Integer64_Linear_Solvers;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;        use Degrees_in_Sets_of_Unknowns;

package body m_Homogeneous_Bezout_Numbers is

-- UTILITIY :

  function Create ( p : Poly_Sys ) return Set is

  -- DESCRIPTION :
  --   Returns the set of the unknowns of the polynomial system p.

    s : Set := Create(p'length);

  begin
    for i in p'range loop
      Add(s,natural32(i));
    end loop;
    return s;
  end Create;

-- TARGET ROUTINES :

  function Cardinalities ( z : Partition ) return Vector is

    res : Vector(integer32(z'first)..integer32(z'last));

  begin
    for i in z'range loop
      res(integer32(i)) := integer32(Extent_Of(z(i)));
    end loop;
    return res;
  end Cardinalities;

  function Total_Degree ( p : Poly_Sys ) return natural32 is

    d : natural32 := 1;

  begin
    for i in p'range loop
      d := d*natural32(Degree(p(i)));
    end loop;
    return d;
  end Total_Degree;

  function Total_Degree ( p : Poly_Sys ) return natural64 is

    d : natural64 := 1;

  begin
    for i in p'range loop
      d := d*natural64(Degree(p(i)));
    end loop;
    return d;
  end Total_Degree;

  function Total_Degree ( p : Poly_Sys ) return Natural_Number is

    d : Natural_Number := Create(natural32(1));

  begin
    for i in p'range loop
      Mul(d,natural32(Degree(p(i))));
    end loop;
    return d;
  end Total_Degree;

  function Bezout_Number ( p : Poly_Sys; z : Partition ) return natural64 is

    k : constant Vector := Cardinalities(z);
    d : constant Matrix := Degree_Table(p,z);
    r : constant integer64 := Per(d,k); -- permanent of degree table

  begin
    return natural64(r);
  end Bezout_Number;

  function Bezout_Number ( p : Poly_Sys; z : Partition; max : natural64 )
                         return natural64 is

    k : constant Vector := Cardinalities(z);
    d : constant Matrix := Degree_Table(p,z);
    r : constant integer64 := Per(d,k,integer64(max));

  begin
    return natural64(r);
  end Bezout_Number;

  function Bezout_Number ( p : Poly_Sys ) return natural64 is

    s : Set := Create(p);
    res : natural64 := Total_Degree(p);

    procedure Evaluate ( z : in Partition; cont : out boolean ) is

      b : constant natural64 := Bezout_Number(p,z,res);

    begin
      if b < res
       then res := b;
      end if;
      cont := true;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  function Bezout_Number ( max : natural64; p : Poly_Sys ) return natural64 is

    s : Set := Create(p);
    res : natural64 := Total_Degree(p);
    cnt : natural64 := 0;

    procedure Evaluate ( z : in Partition; cont : out boolean ) is

      b : constant natural64 := Bezout_Number(p,z,res);

    begin
      if b < res
       then res := b;
      end if;
      cnt := cnt + 1;
      if cnt < max
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  function Bezout_Number ( p : Poly_Sys; min : natural64 ) return natural64 is

    s : Set := Create(p);
    res : natural64 := Total_Degree(p);

    procedure Evaluate ( z : in Partition; cont : out boolean ) is

      b : constant natural64 := Bezout_Number(p,z,res);

    begin
      if b < res
       then res := b;
      end if;
      if res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  function Bezout_Number ( max : natural64; p : Poly_Sys; min : natural64 )
                         return natural64 is

    s : Set := Create(p);
    res : natural64 := Total_Degree(p);
    cnt : natural64 := 0;   

    procedure Evaluate ( z : in Partition; cont : out boolean ) is

      b : constant natural64 := Bezout_Number(p,z,res);

    begin
      if b < res
       then res := b;
      end if;
      cnt := cnt + 1;
      if cnt < max and then res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  procedure Bezout_Number 
               ( p : in Poly_Sys; b : out natural64;
                 m : out natural32; z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural64 := Total_Degree(p);
    res : natural64 := tdg;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is

      nb : constant natural64 := Bezout_Number(p,nz,res);

    begin
      if nb < res then
        res := nb;
        m := nz'length; Clear(z);
        z(nz'range) := Create(nz);
      end if;
      cont := true;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;
 
  procedure Bezout_Number 
               ( max : in natural64; p : in Poly_Sys; b : out natural64;
                 m : out natural32; z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural64 := Total_Degree(p);
    res : natural64 := tdg;
    cnt : natural64 := 0;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is

      nb : constant natural64 := Bezout_Number(p,nz,res);

    begin
      if nb < res then
        res := nb;
        m := nz'length; Clear(z);
        z(nz'range) := Create(nz);
      end if;
      cnt := cnt + 1;
      if cnt < max
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;

  procedure Bezout_Number
               ( p : in Poly_Sys; min : in natural64; b : out natural64;
                 m : out natural32; z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural64 := Total_Degree(p);
    res : natural64 := tdg;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is

      nb : constant natural64 := Bezout_Number(p,nz,res);

    begin
      if nb < res then
        res := nb;
        m := nz'length; Clear(z);
        z(nz'range) := Create(nz);
      end if;
      if res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;

  procedure Bezout_Number
              ( max : in natural64; p : in Poly_Sys; min : in natural64;
                b : out natural64; m : out natural32; z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural64 := Total_Degree(p);
    res : natural64 := tdg;
    cnt : natural64 := 0;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is

      nb : constant natural64 := Bezout_Number(p,nz,res);

    begin
      if nb < res then
        res := nb;
        m := nz'length; Clear(z);
        z(nz'range) := Create(nz);
      end if;
      cnt := cnt + 1;
      if cnt < max and then res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;

  function Evaluate ( z : partition; m : natural32; p : Poly_Sys )
                    return natural64 is

    n : constant natural32 := natural32(p'length);
    d : constant Matrix := Degree_Table(p,z);

    function Bezout_number 
               ( n,m : natural32; z: partition; d : Matrix ) return natural64 is

    -- DESCRIPTION : the Bezout number is computed

      type boolean_array is array ( natural32 range <> ) of boolean;
      Ii,Iacc : boolean_array(1..n) := (1..n => false);
      b,b_mult : natural64;

      procedure column ( j,start,number : in natural32;
                         Ii,Iacc : in out boolean_array );

      -- DESCRIPTION : the computation of a term coming from a column
      --               of the degree table;

      procedure row ( j : in natural32; Ii : in out boolean_array;
                      Iacc : in boolean_array );

      -- DESCRIPTION : the computation of a row in the degree table

      procedure column ( j,start,number : in natural32;
                         Ii,Iacc : in out boolean_array ) is
      begin
        if number > (n - start + 1) then
          return;
        elsif number = 0 then
          row(j,Ii,Iacc);
        else
          for i in start..n loop
            if not Ii(i) then
              Ii(i) := true; Iacc(i) := true;
              column(j,i+1,number-1,Ii,Iacc);
              Ii(i) := false; Iacc(i) := false;
            end if;
          end loop;
        end if;
      end column;

      procedure row ( j : in natural32; Ii : in out boolean_array;
                      Iacc : in boolean_array ) is

        temp : natural64 := 1;
        Iacc1 : boolean_array(1..n) := (1..n => false);

      begin
        for k in 1..n loop
          if Iacc(k)
           then temp := temp*natural64(d(integer32(k),integer32(j)));
          end if;
        end loop;
        if (j /= m) and (temp /= 0) then
          b_mult := b_mult*temp;
          column(j+1,1,Extent_Of(z(j+1)),Ii,Iacc1);
          b_mult := b_mult/temp;
        elsif j = m then
          temp := temp*b_mult;
          b := b + temp;
        end if;
      end row;

    begin
      b := 0; b_mult := 1;
      column(1,1,Extent_Of(z(1)),Ii,Iacc);
      return b;
    end Bezout_number;

  begin
    return Bezout_number(n,m,z,d);
  end Evaluate;

  procedure PB ( p : in Poly_Sys; b : out natural64;
                 m : out natural32; z : in out Partition ) is

    n : constant natural32 := natural32(p'length);
    wb,b_min : natural64;
    wz,z_min : partition(1..n);
    wm,m_min : natural32 := 0;

    procedure pcopy ( p1 : in partition; p2 : in out partition;
                      p1n : in natural32; p2n : in out natural32 ) is

    -- DESCRIPTION : the partition p1 is copied to p2

    begin
      Clear(p2); p2 := Create(p1);
      p2n := p1n;
    end pcopy;

  begin
    b_min := Total_Degree(p);
    for i in 1..n loop
      for k in 1..wm loop
        Add(wz(k),i); 
        wb := Evaluate(wz,wm,p);
        if (k = 1) or else (wb < b_min) then
          pcopy(wz,z_min,wm,m_min);
          b_min := wb;
        end if;
        Remove(wz(k),i);
      end loop;
      wm := wm + 1;
      wz(wm) := Create(n);
      Add(wz(wm),i); 
      wb := Evaluate(wz,wm,p);
      if wb < b_min then
        pcopy(wz,z_min,wm,m_min);
        b_min := wb;
      end if;
      pcopy(z_min,wz,m_min,wm);
    end loop;
    b := b_min;
    m := m_min;
    pcopy(z_min,z,m_min,wm);
    Clear(wz);
  end PB;

  procedure Patch ( p : in Poly_Sys; z : in out Partition;
                    nz : in natural32; bz : in out natural64 ) is

    dim : natural32 := natural32(p'last);

  begin
    if nz = 1 then
      if Extent_Of(z(1)) < dim then
        Clear(z(1));
        z(1) := Universe(dim);
        bz := Total_Degree(p);
      end if;
    end if;
  end Patch;

end m_Homogeneous_Bezout_Numbers;
