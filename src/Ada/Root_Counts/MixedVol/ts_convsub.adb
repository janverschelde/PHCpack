with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Floating_Mixed_Subdivisions;        use Floating_Mixed_Subdivisions;
with Floating_Mixed_Subdivisions_io;     use Floating_Mixed_Subdivisions_io;

procedure ts_convsub is

-- DESCRIPTION :
--   Converts the output file of Xing Li his program into
--   a floating-point mixed subdivision.

  procedure Read_Lifted_Supports
              ( file : in file_type; n : in integer32;
                lifsup : out Standard_Floating_VecVecs.Array_of_VecVecs ) is

  -- DESCRIPTION :
  --   Reads the lifted supports from file.

    k : integer32 := 0;

  begin
    for i in lifsup'range loop
      get(file,k);
      declare
        pts : Standard_Floating_VecVecs.VecVec(1..k);
        lpt : Standard_Floating_Vectors.Vector(1..n+1);
      begin
        for j in 1..k loop
          get(file,lpt);
          pts(j) := new Standard_Floating_Vectors.Vector'(lpt);
        end loop;
        lifsup(i) := new Standard_Floating_VecVecs.VecVec'(pts);
      end;
    end loop;
  end Read_Lifted_Supports;

  function Minimum ( lifsup : in Standard_Floating_VecVecs.VecVec;
                     v : Standard_Floating_Vectors.Vector ) 
                   return double_float is

  -- DESCRIPTION :
  --   Returns the minimal inner product the points in lifsup make
  --   with the given vector v.

    use Standard_Floating_Vectors;

    min : double_float := lifsup(1).all*v;
    ipr : double_float;

  begin
    for i in 2..lifsup'last loop
      ipr := lifsup(i).all*v;
      if ipr < min
       then min := ipr;
      end if;
    end loop;
    return min;
  end Minimum;

  procedure Test_Inner_Normal
               ( n : in integer32; v : in Standard_Floating_Vectors.Vector;
                 lifsup : in Standard_Floating_VecVecs.Array_of_VecVecs;
                 pts : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks whether the minimum really occurs at the points in pts.

    min,ipr : double_float;
    cnt : integer32 := 0;
    use Standard_Floating_Vectors;

  begin
    for i in 1..n loop
      min := Minimum(lifsup(i).all,v);
      put("min : "); put(min,3);
      cnt := cnt+1;
      put(" at "); put(pts(cnt)+1,1); put(" : ");
      ipr := lifsup(i)(pts(cnt)+1).all*v;
      put(ipr,3);
      cnt := cnt+1;
      put(" at "); put(pts(cnt)+1,1); put(" : ");
      ipr := lifsup(i)(pts(cnt)+1).all*v;
      put(ipr,3);
      new_line;
    end loop;
  end Test_Inner_Normal;

  procedure Create_Cell
               ( file : in file_type; n : in integer32;
                 lifsup : in Standard_Floating_VecVecs.Array_of_VecVecs;
                 mic : out Mixed_Cell ) is

  -- DESCRIPTION :
  --   Reads one line from file and returns the corresponding mixed cell.
  --   This involves also the computation of the inner normal.

    pt1,pt2 : integer32 := 0;
    mat : Matrix(1..n,1..n);
    rhs : Standard_Floating_Vectors.Vector(1..n);
    ipvt : Standard_Integer_Vectors.Vector(1..n);
    info : integer32;
    cnt : integer32 := 0;
    pts : Standard_Integer_Vectors.Vector(1..2*n);

  begin
    mic.pts := new Array_of_Lists(1..n);
    for i in 1..n loop
      get(file,pt1); cnt := cnt+1; pts(cnt) := pt1;
      get(file,pt2); cnt := cnt+1; pts(cnt) := pt2;
      declare
        last : List;
        point1 : constant Standard_Floating_Vectors.Vector
               := lifsup(i)(pt1+1).all;
        point2 : constant Standard_Floating_Vectors.Vector
               := lifsup(i)(pt2+1).all;
      begin
        Append(mic.pts(i),last,point1);
        Append(mic.pts(i),last,point2);
        for j in 1..n loop
          mat(i,j) := point1(j) - point2(j);
        end loop;
        rhs(i) := point2(n+1) - point1(n+1);
      end;
    end loop;
    skip_line(file);
    lufac(mat,n,ipvt,info);
    lusolve(mat,n,ipvt,rhs);
    mic.nor := new Standard_Floating_Vectors.Vector(1..n+1);
    mic.nor(1..n) := rhs;
    mic.nor(n+1) := 1.0;
    Test_Inner_Normal(n,mic.nor.all,lifsup,pts);
  end Create_Cell;

  procedure Create_Subdivision
               ( file : in file_type; n,m : in integer32;
                 lifsup : in Standard_Floating_VecVecs.Array_of_VecVecs;
                 sub : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Reads m lines from file and creates the mixed subdivision.

    mic : Mixed_Cell;
    sub_last : Mixed_Subdivision;

  begin
    for i in 1..m loop
      Create_Cell(file,n,lifsup,mic);
      Append(sub,sub_last,mic);
    end loop;
  end Create_Subdivision;

  procedure Read_Subdivision
               ( file : in file_type; n : in integer32;
                 mix : out Standard_Integer_Vectors.Link_to_Vector;
                 sub : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Reads the data from file and constructs a mixed subdivision.
  --   The input dimension n has alredy been read.

    lifsup : Standard_Floating_VecVecs.Array_of_VecVecs(1..n);
    m : integer32 := 0;
    mv : natural32;
    ans : character;

  begin
    put("The dimension is "); put(n,1); put_line(".");
    mix := new Standard_Integer_Vectors.Vector'(1..n => 1);
    put("The type of mixture : "); put(mix.all); put_line(".");
    Read_Lifted_Supports(file,n,lifsup);
    put_line("The lifted supports : ");
    for i in lifsup'range loop
      put(lifsup(i)'last,1); new_line;
      for j in lifsup(i)'range loop
        for k in 1..n loop
          put(" ");
          put(integer32(lifsup(i)(j)(k)),1);
        end loop;
        put(lifsup(i)(j)(n+1)); new_line;
      end loop;
    end loop;
    get(file,m);
    put("The number of cells is "); put(m,1); put_line(".");
    Create_Subdivision(file,n,m,lifsup,sub);
    put("Do you want to continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put(natural32(n),mix.all,sub,mv);
      put("The mixed volume is "); put(mv,1); put_line(".");
    end if;
  end Read_Subdivision;

  procedure Main is

    infile,outfile : file_type;
    n : integer32 := 0;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    sub : Mixed_Subdivision;

  begin
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(infile);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
   -- new_line;
   -- put_line("See the output file for results...");
   -- new_line;
    get(infile,n);
    Read_Subdivision(infile,n,mix,sub);
    put(outfile,natural32(n),mix.all,sub);
  end Main;

begin
  new_line;
  put_line("Converting output file of Xing Li's program"
            & " into mixed subdivision.");
  new_line;
  Main;
end ts_convsub;
