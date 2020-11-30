with Numbers_io;                        use Numbers_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Complex_VecMats_io;       use Standard_Complex_VecMats_io;
with Osculating_Planes;                 use Osculating_Planes;
with Complex_Osculating_Planes;         use Complex_Osculating_Planes;

package body Make_Input_Planes is

  function Finite ( dim : Bracket; fin_sum : natural32 ) return boolean is

    sum : natural32 := 0;

  begin
    for i in dim'range loop
      sum := sum + dim(i);
    end loop;
    if sum = fin_sum
     then return true;
     else return false;
    end if;
  end Finite;

  function Read_Codimensions ( m,p,q : natural32 ) return Bracket is

    mpq : constant natural32 := m*p + q*(m+p);
    codim : Bracket(1..integer32(mpq));
    n : natural32 := 0;

  begin
    loop
      put("Give number of intersection conditions : "); get(n);
      put("Give "); put(n,1); put(" co-dimensions (sum = ");
      put(mpq,1); put(") : ");
      for i in 1..integer32(n) loop
        codim(i) := 0; get(codim(i));
      end loop;
      for i in 1..integer32(n)-1 loop
        put(codim(i),1); put(" + ");
      end loop;
      put(codim(integer32(n)),1);
      if Finite(codim(1..integer32(n)),mpq) then
        put(" = "); put(mpq,1); put_line("  Finite #sols.");
        exit;
      else 
        put(" /= "); put(mpq,1);
        put_line("  Please try again.");
      end if;
    end loop;
    return codim(1..integer32(n));
  end Read_Codimensions;

  function Convert ( mat : Standard_Floating_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Converts a real matrix into a complex one.

    res : Standard_Complex_Matrices.Matrix(mat'range(1),mat'range(2));

  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        res(i,j) := Create(mat(i,j));
      end loop;
    end loop;
    return res;
  end Convert;

  function Random_Complex_Planes ( m,p : natural32 ) return VecMat is

    res : VecMat(1..integer32(m*p));
    n : constant natural32 := m+p;
    use Standard_Complex_Matrices;

  begin
    for i in res'range loop
      res(i) := new Matrix'(Random_Orthogonal_Matrix(n,m));
    end loop;
    return res;
  end Random_Complex_Planes;

  function Random_Complex_Planes
             ( m,p : natural32; k : Bracket ) return VecMat is

    res : VecMat(k'range);
    n : constant natural32 := m+p;
    use Standard_Complex_Matrices;

  begin
    for i in res'range loop
      res(i) := new Matrix'(Random_Orthogonal_Matrix(n,m+1-k(i)));
    end loop;
    return res;
  end Random_Complex_Planes;

  function Random_Complex_Planes ( m,p,q : natural32 ) return VecMat is

    res : VecMat(1..integer32(m*p+q*(m+p)));
    n : constant natural32 := m+p;
    use Standard_Complex_Matrices;

  begin
    for i in res'range loop
      res(i) := new Matrix'(Random_Orthogonal_Matrix(n,m));
    end loop;
    return res;
  end Random_Complex_Planes;      

  function Random_Real_Planes ( m,p : natural32 ) return VecMat is

    res : VecMat(1..integer32(m*p));
    n : constant natural32 := m+p;
    fltmat : Standard_Floating_Matrices.Matrix
               (1..integer32(n),1..integer32(m));
    cmpmat : Standard_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range loop
      fltmat := Random_Orthogonal_Matrix(n,m);
      cmpmat := Convert(fltmat);
      res(i) := new Standard_Complex_Matrices.Matrix'(cmpmat);
    end loop;
    return res;
  end Random_Real_Planes;

  function Random_Real_Planes ( m,p : natural32; k : Bracket ) return VecMat is

    res : VecMat(k'range);
    n : constant natural32 := m+p;

  begin
    for i in res'range loop
      declare
        fltmat : Standard_Floating_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-k(i)));
        cmpmat : Standard_Complex_Matrices.Matrix
                   (1..integer32(n),1..integer32(m+1-k(i)));
      begin
        fltmat := Random_Orthogonal_Matrix(n,m+1-k(i));
        cmpmat := Convert(fltmat);
        res(i) := new Standard_Complex_Matrices.Matrix'(cmpmat);
      end;
    end loop;
    return res;
  end Random_Real_Planes;

  function Random_Real_Planes ( m,p,q : natural32 ) return VecMat is

    res : VecMat(1..integer32(m*p+q*(m+p)));
    n : constant natural32 := m+p;
    fltmat : Standard_Floating_Matrices.Matrix(1..integer32(n),1..integer32(m));
    cmpmat : Standard_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range loop
      fltmat := Random_Orthogonal_Matrix(n,m);
      cmpmat := Convert(fltmat);
      res(i) := new Standard_Complex_Matrices.Matrix'(cmpmat);
    end loop;
    return res;
  end Random_Real_Planes;      

  function Equidistant_Interpolation_Points ( n : natural32 ) return Vector is

    res : Vector(1..integer32(n));
    s : double_float := Random;
    inc : constant double_float := 2.0/double_float(n);

  begin
    for i in res'range loop
      res(i) := s;
      s := s+inc;
      if s >= 1.0
       then s := s - 2.0;
      end if;
    end loop;
    return res;
  end Equidistant_Interpolation_Points;

  function Read_Interpolation_Points ( n : natural32 ) return Vector is

    res : Vector(1..integer32(n));

  begin
    new_line;
    put("Give "); put(n,1); put_line(" distinct real interpolation points : ");
    for i in res'range loop
      Read_Double_Float(res(i));
    end loop;       
    return res;
  end Read_Interpolation_Points;

  function Osculating_Input_Planes ( m,p : natural32 ) return VecMat is

    s : constant Vector := Equidistant_Interpolation_Points(m*p);

  begin
    return Osculating_Input_Planes(m,p,s);
  end Osculating_Input_Planes;

  function Complex_Osculating_Input_Planes ( m,p : natural32 ) return VecMat is

    s : constant Vector := Equidistant_Interpolation_Points(m*p);

  begin
    return Complex_Osculating_Input_Planes(m,p,s);
  end Complex_Osculating_Input_Planes;

  function Osculating_Input_Planes
             ( m,p : natural32; s : Vector ) return VecMat is

    res : VecMat(1..integer32(m*p));
    n : constant natural32 := m+p;
    realmat : Standard_Floating_Matrices.Matrix
                (1..integer32(n),1..integer32(m));

  begin
    for i in res'range loop
      realmat := Orthogonal_Basis(n,m,s(i));
      res(i) := new Standard_Complex_Matrices.Matrix'(Convert(realmat));
    end loop;
    return res;
  end Osculating_Input_Planes;

  function Complex_Osculating_Input_Planes
             ( m,p : natural32; s : Vector ) return VecMat is

    res : VecMat(1..integer32(m*p));
    n : constant natural32 := m+p;
    mat : Standard_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range loop
      mat := Standard_Basis(n,m,Create(0.0,s(i)));
      res(i) := new Standard_Complex_Matrices.Matrix'(mat);
    end loop;
    return res;
  end Complex_Osculating_Input_Planes;

  function Osculating_Input_Planes
             ( m,p : natural32; k : bracket ) return VecMat is

    s : constant Vector
      := Equidistant_Interpolation_Points(natural32(k'length));

  begin
    return Osculating_Input_Planes(m,p,k,s);
  end Osculating_Input_Planes;

  function Osculating_Input_Planes
             ( m,p : natural32; k : bracket; s : Vector ) return VecMat is

    res : VecMat(k'range);
    n : constant natural32 := m+p;

  begin
    for i in res'range loop
      declare
        realmat : Standard_Floating_Matrices.Matrix
                    (1..integer32(n),1..integer32(m+1-k(i)));
      begin
        realmat := Orthogonal_Basis(n,m+1-k(i),s(i));
        res(i) := new Standard_Complex_Matrices.Matrix'(Convert(realmat));
      end;
    end loop;
    return res;
  end Osculating_Input_Planes;

  function Osculating_Input_Planes ( m,p,q : natural32 ) return VecMat is

    dim : constant natural32 := m*p + q*(m+p);
    s : constant Vector := Equidistant_Interpolation_Points(dim);

  begin
    return Osculating_Input_Planes(m,p,q,s);
  end Osculating_Input_Planes;

  function Osculating_Input_Planes
             ( m,p,q : natural32; s : Vector ) return VecMat is

    res : VecMat(1..integer32(m*p+q*(m+p)));
    n : constant natural32 := m+p;
    realmat : Standard_Floating_Matrices.Matrix
                (1..integer32(n),1..integer32(m));

  begin
    for i in res'range loop
      realmat := Orthogonal_Basis(n,m,s(i));
      res(i) := new Standard_Complex_Matrices.Matrix'(Convert(realmat));
    end loop;
    return res;
  end Osculating_Input_Planes;

  function Read_Input_Planes ( m,p : natural32 ) return VecMat is

    res : VecMat(1..integer32(m*p));
    planesfile : file_type;
    n : constant natural32 := m+p;
    realmat : Standard_Floating_Matrices.Matrix
                (1..integer32(n),1..integer32(m));

  begin
    new_line;
    put_line("Reading the name of the file with the input planes.");
    Read_Name_and_Open_File(planesfile);
    for i in res'range loop
      get(planesfile,realmat);
     -- realmat := Orthogonalize(realmat);
      res(i) := new Standard_Complex_Matrices.Matrix'(Convert(realmat));
    end loop;
    Close(planesfile);
    return res;
  end Read_Input_Planes;

  function Read_Complex_Input_Planes ( m,p : natural32 ) return VecMat is

    res : VecMat(1..integer32(m*p));
    planesfile : file_type;
    n : constant natural32 := m+p;
    compmat : Standard_Complex_Matrices.Matrix
                (1..integer32(n),1..integer32(m));

  begin
    new_line;
    put_line("Reading the name of the file with the input planes.");
    Read_Name_and_Open_File(planesfile);
    for i in res'range loop
      get(planesfile,compmat);
     -- compmat := Orthogonalize(compmat);
      res(i) := new Standard_Complex_Matrices.Matrix'(compmat);
    end loop;
    Close(planesfile);
    return res;
  end Read_Complex_Input_Planes;

  function Read_Input_Planes 
             ( m,p : natural32; k : in Bracket ) return VecMat is

    res : VecMat(k'range);
    planesfile : file_type;
    n : constant natural32 := m+p;

  begin
    new_line;
    put_line("Reading the name of the file with the input planes.");
    Read_Name_and_Open_File(planesfile);
    for i in res'range loop
      declare
        realmat : Standard_Floating_Matrices.Matrix
                    (1..integer32(n),1..integer32(m+1-k(i)));
      begin
        get(planesfile,realmat);
        realmat := Orthogonalize(realmat);
        res(i) := new Standard_Complex_Matrices.Matrix'(Convert(realmat));
      end;
    end loop;
    Close(planesfile);
    return res;
  end Read_Input_Planes;

  function Read_Input_Planes ( m,p,q : natural32 ) return VecMat is

    res : VecMat(1..integer32(m*p+q*(m+p)));
    planesfile : file_type;
    n : constant natural32 := m+p;
    realmat : Standard_Floating_Matrices.Matrix
                (1..integer32(n),1..integer32(m));

  begin
    new_line;
    put_line("Reading the name of the file with the real input planes.");
    Read_Name_and_Open_File(planesfile);
    for i in res'range loop
      get(planesfile,realmat);
      realmat := Orthogonalize(realmat);
      res(i) := new Standard_Complex_Matrices.Matrix'(Convert(realmat));
    end loop;
    Close(planesfile);
    return res;
  end Read_Input_Planes;

-- MAIN INTERACTIVE DRIVERS :

  function Select_Input_Choice ( hypersurface : boolean ) return character is

  -- DESCRIPTION :
  --   Displays the menu to obtain a choice for the kind of input.
  --   If hypersurface, then the character on return is in 0..6,
  --   otherwise the character on return is in 0..4.

    res : character;

  begin
    new_line;
    put_line("MENU for generating real input planes : ");
    put_line("  0. Run only Pieri homotopies, no cheater homotopy.");
    put_line("  1. Random real planes at equidistant interpolation points.");
    put_line("  2. Generate input planes osculating a "
                               & "real rational normal curve.");
    put_line("  3. Interactively give real s-values for the "
                               & "osculating input planes.");
    put_line("  4. Give the name of the file with real input planes.");
    if hypersurface then
      put_line("  5.                           with complex input planes.");
      put_line("  6. Test variant of Shapiro^2 conjecture.");
      put("Type 0, 1, 2, 3, 4, 5, or 6 to select : ");
      Ask_Alternative(res,"0123456");
    else
      put("Type 0, 1, 2, 3, or 4 to select : ");
      Ask_Alternative(res,"01234");
    end if;
    return res;
  end Select_Input_Choice;

  procedure Main ( file : file_type; m,p : in natural32; planes : out VecMat;
                   nocheater : out boolean ) is

    input_choice : constant character := Select_Input_Choice(true);
    dim : constant natural32 := m*p;
    svals : Vector(1..integer32(dim));

  begin
    nocheater := false;
    case input_choice is
      when '0' => nocheater := true;
      when '1' => planes := Random_Real_Planes(m,p);
      when '2' => planes := Osculating_Input_Planes(m,p);
      when '3' => svals := Read_Interpolation_Points(dim);
                  put_line(file,"The interpolation points : ");
                  put_line(file,svals);
                  planes := Osculating_Input_Planes(m,p,svals);
      when '4' => planes := Read_Input_Planes(m,p);
      when '5' => planes := Read_Complex_Input_Planes(m,p);
      when '6' => svals := Read_Interpolation_Points(dim);
                  put_line(file,"The interpolation points : ");
                  put_line(file,svals);
                  planes := Complex_Osculating_Input_Planes(m,p,svals);
      when others => null;
    end case;
    if not nocheater then
      put_line(file,"THE INPUT PLANES (eventually after othogonalization) : ");
      put(file,planes);
    end if;
  end Main;

  procedure Main ( m,p : in natural32; k : in Bracket; planes : out VecMat;
                   nocheater : out boolean ) is

    input_choice : constant character := Select_Input_Choice(false);
    svals : Vector(k'range);

  begin
    nocheater := false;
    case input_choice is
      when '0' => nocheater := true;
      when '1' => planes := Random_Real_Planes(m,p,k);
      when '2' => planes := Osculating_Input_Planes(m,p,k);
      when '3' => svals := Read_Interpolation_Points(k'length);
                  planes := Osculating_Input_Planes(m,p,k,svals);
      when '4' => planes := Read_Input_Planes(m,p,k);
      when others => null;
    end case;
  end Main;

  procedure Main ( file : in file_type;
                   m,p : in natural32; k : in Bracket; planes : out VecMat;
                   nocheater : out boolean ) is
  begin
    Main(m,p,k,planes,nocheater);
    if not nocheater then
      put_line(file,"THE INPUT PLANES (eventually after othogonalization) : ");
      put(file,planes);
    end if;
  end Main;

  procedure Main ( file : in file_type; m,p,q : in natural32;
                   s : out Vector; planes : out VecMat;
                   nocheater : out boolean ) is

    input_choice : constant character := Select_Input_Choice(false);
    dim : constant natural32 := m*p + q*(m+p);
    svals : Vector(1..integer32(dim));

  begin
    nocheater := false;
    case input_choice is
      when '0' => nocheater := true;
      when '1' => svals := Equidistant_Interpolation_Points(dim);
                  s := svals;
                  planes := Random_Real_Planes(m,p,q);
      when '2' => svals := Equidistant_Interpolation_Points(dim);
                  s := svals;
                  planes := Osculating_Input_Planes(m,p,q,svals);
      when '3' => svals := Read_Interpolation_Points(dim);
                  s := svals;
                  planes := Osculating_Input_Planes(m,p,q,svals);
      when '4' => svals := Read_Interpolation_Points(dim);
                  s := svals;
                  planes := Read_Input_Planes(m,p,q);
      when others => null;
    end case;
    if not nocheater then
      put_line(file,"THE INTERPOLATION POINTS : ");
      put_line(file,svals);
      put_line(file,"THE TARGET PLANES : ");
      put(file,planes);
    end if;
  end Main;

end Make_Input_Planes;
