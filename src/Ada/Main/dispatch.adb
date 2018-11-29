with text_io;                            use text_io;
with Unix_Command_Line;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Random_Numbers; -- to handle zero seed
with Greeting_Banners;
with mainscal,mainred;        -- scaling and reduction
with mainred2,mainred4;       -- reduction for double double and quad doubles
with mainroco,bablroco;       -- general root counting
with mainsmvc,babldmvc;       -- mixed-volume computation
with maintrack;               -- path tracking
with mainpoco,bablpoco;       -- polynomial continuation
with bablpoco2,bablpoco4;     -- double double and quad double continuation
with mainadep;                -- path tracking with algorithmic differentiation
with mainphc,bablphc;         -- main phc driver + blackbox
with bablphc2,bablphc4;       -- blackbox in double double and quad double
with compsolve;               -- numerical irreducible decomposition as -B
with compsolve2,compsolve4;   -- dobldobl and quaddobl versions, -B2 and -B4
with mainvali,bablvali;       -- verification tool
with bablvali2,bablvali4;     -- verification with dd and qd
with mainenum,bablenum;       -- numerical Schubert calculus
with mainfeed;                -- realization of pole placing feedback
with maindict;                -- converts solutions into Python dictionary
with mainzip;                 -- putting solutions in Maple format
with mainsam;                 -- computing samples given a witness set
with mainfac,maindeco;        -- factorization and decomposition
with mainsolve;               -- equation-by-equation solver
with mainwit;                 -- witness set intersection
with maingood;                -- to test if a system is good
with mainsymb;                -- to get the symbol table contents
with mainhyp;                 -- witness set for hypersurface
with mainhyp2,mainhyp4;       -- double double and quad double versions
with mainseries;              -- Newton's method for power series solutions
-- NOTE (added for pieri_solver.ali) :
with Interfaces.C;
with Complex_Polynomial_Matrices;        use Complex_Polynomial_Matrices;
with Complex_Polynomial_Matrices_io;     use Complex_Polynomial_Matrices_io;
with Verify_Solution_Maps;
with C_to_Ada_Arrays;                    use C_to_Ada_Arrays;
with Write_Seed_Number;

with Standard_Natural_Numbers_io;
 use Standard_Natural_Numbers_io;

procedure Dispatch is

-- BANNERS WITH INFORMATION TO START DIALOGUE WITH USER :

  welcome : constant string := Greeting_Banners.welcome;
  compban : constant string := Greeting_Banners.compban;
  enumban : constant string := Greeting_Banners.enumban;
  facban  : constant string := Greeting_Banners.facban;
  goodban : constant string := Greeting_Banners.goodban;
  symbban : constant string := Greeting_Banners.symbban;
  hypban  : constant string := Greeting_Banners.hypban;
  mvcban  : constant string := Greeting_Banners.mvcban;
  pocoban : constant string := Greeting_Banners.pocoban;
  reduban : constant string := Greeting_Banners.reduban;
  rocoban : constant string := Greeting_Banners.rocoban;
  samban  : constant string := Greeting_Banners.samban;
  scalban : constant string := Greeting_Banners.scalban;
  slvban  : constant string := Greeting_Banners.slvban;
  adepban : constant string := Greeting_Banners.adepban;
  trackban : constant string := Greeting_Banners.trackban;
  seriesban : constant string := Greeting_Banners.seriesban;
  veriban : constant string := Greeting_Banners.veriban;
  witban  : constant string := Greeting_Banners.witban;

-- AVAILABLE OPTIONS :

  options : constant string := "bB0asdpqmrvekcxyzftwlgo-huj";
  -- 0 : zero seed for repeatable runs
  -- a : solve => equation-by-equation solver
  -- b : batch or black box processing
  -- c : comp => numerical irreducible decomposition
  -- d : redu => reduction w.r.t. the total degree
  -- e : enum => numerical Schubert calculus
  -- f : fac  => factor pure dimensional solution set into irreducibles
  -- g : good => check if the input is a system in the valid format
  -- h : help => write information about a certain option
  -- j : adep => algorithmic differentiation path trackers
  -- k : feed => dynamic output feedback to control linear systems
  -- l : hyp  => witness set for hypersurface cutting with random line
  -- m : mvc  => mixed-volume computation
  -- o : symbOls -> get symbols using as variables in a system
  -- p : poco => polynomial continuation
  -- q : track => tracking solution paths
  -- r : roco => root counting methods
  -- s : scal => scaling of a polynomial system
  -- t : task => use multitasking
  -- u : series => Newton's method for power series solutions
  -- v : vali => validation of solutions
  -- w : wit  => intersection of witness sets using diagonal homotopies
  -- x : dix  => Python dictionary output for solutions
  -- y : sam  => sampling points from an algebraic set
  -- z : zip  => strip output for solutions into Maple format

  argc : constant integer32
       := integer32(Unix_Command_Line.Number_of_Arguments);

  fixed_seed : constant natural32 := 68717;

-- UTILITIES FOR PROCESSING THE ARGUMENTS AND OPTIONS :

  function Read_Argument ( k : in integer32 ) return string is

  -- DESCRIPTION :
  --   Reads the kth argument from the command line.
  --   An argument is a string not proceeded by a `-' character.
  --   The empty string is returned when there is no argument.

    null_string : constant string := "";
    cnt : integer32 := 0;

  begin
    if argc >= 1 then
      for i in 1..argc loop
        declare
          s : constant string := Unix_Command_Line.Argument(integer(i));
        begin
          if s(1) /= '-' then
            cnt := cnt + 1;
            if k = cnt
             then return s;
            end if;
          end if;
        end;
      end loop;
    end if;
    return null_string;
  end Read_Argument;

  function Position ( c : character; s : string ) return integer32 is

  -- DESCRIPTION :
  --   If the the string contains the character c, then its position
  --   in the string will be returned.  Otherwise s'first-1 will be returned.

  begin
    for i in s'range loop
      if s(i) = c
       then return integer32(i);
      end if;
    end loop;
    return integer32(s'first-1);
  end Position;

  procedure Read_Next_Option
              ( pos : in out integer32; legal : in string;
                option : out character ) is

  -- DESCRIPTION :
  --   Reads the next option from the command line arguments.

  -- ON ENTRY :
  --   pos      position in the command line of the last option
  --            that has been read;
  --   legal    string which contains all legal options.

  -- ON RETURN :
  --   pos      the position in the command line of the last option read;
  --   option   is blank when no legal option could be read, otherwise it
  --            contains the next legal option.

    res : character := ' ';
    start : constant integer32 := pos+1;

  begin
    if argc >= 1 then
      for i in start..argc loop
        declare
          s : constant string := Unix_Command_Line.Argument(integer(i));
        begin
          if s(1) = '-' then
            pos := Position(s(2),legal);
            if pos >= integer32(legal'first) then
              res := legal(integer(pos));
            else
              put("The option '"); put(s);
              put_line("' is not recognised.  Will ignore it...");
            end if;
          end if;
        end;
        pos := i;
        exit when (res /= ' ');
      end loop;
    end if;
    option := res;
  end Read_Next_Option;

  function Number_of_Tasks ( i : integer32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the number of tasks of the i-th argument on the command line.

  -- REQUIRED :
  --   The i-th argument is a string starting with "-t"
  --   and followed by a natural number.

    s : constant string := Unix_Command_Line.Argument(integer(i));
    res : constant natural32 := Convert(s(3..s'last));

  begin
    return res;
  end Number_of_Tasks;

  function Number_of_Tasks return natural32 is

  -- DESCRIPTION :
  --   Returns the number of tasks of the argument -t on the command line.
  --   If there is no argument -t, then 0 is returned.

    res : natural32 := 0;

  begin
    for i in 1..Unix_Command_Line.Number_of_Arguments loop
      declare
        s : constant string := Unix_Command_Line.Argument(i);
      begin
        if s(2) = 't'
         then res := Convert(s(3..s'last)); exit;
        end if;
      end;
    end loop;
    return res;
  end Number_of_Tasks;

  function Scan_Precision ( opt : character ) return natural32 is

  -- DESCRIPTION :
  --   Returns the precision of the option defined by the charactor opt.
  --   1 : the -opt is followed by a space (or nothing);
  --   2 : double double precision, as we have -opt2 at the command line;
  --   4 : quad double precision is given as -opt4 at the command line.

    res : natural32 := 1;

  begin
    for i in 1..Unix_Command_Line.Number_of_Arguments loop
      declare
        s : constant string := Unix_Command_Line.Argument(i);
      begin
        if s(2) = opt then
          if s'last > 2 
           then res := Convert(s(3..s'last)); exit;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Scan_Precision;

  function Scan_Series_Precision ( opt : character ) return natural32 is

  -- DESCRIPTION :
  --   Returns the precision of the option defined by the charactor opt,
  --   for the power series solutions with Newton's method.
  --   0 : the -opt is followed by a space (or nothing);
  --   1 : standard double precision, as we have -opt1 at the command line;
  --   2 : double double precision, as we have -opt2 at the command line;
  --   4 : quad double precision is given as -opt4 at the command line.

    res : natural32 := 0;

  begin
    for i in 1..Unix_Command_Line.Number_of_Arguments loop
      declare
        s : constant string := Unix_Command_Line.Argument(i);
      begin
        if s(2) = opt then
          if s'last > 2 
           then res := Convert(s(3..s'last)); exit;
          end if;
        end if;
      end;
    end loop;
    return res;
  end Scan_Series_Precision;

  function Find_Seed return natural32 is

  -- DESCRIPTION :
  --   Reads the digits after the '-0' and returns the seed.
  --   If there is nothing after the '-0' then 0 is returned.

    res : natural32 := 0;

  begin
    for i in 1..Unix_Command_Line.Number_of_Arguments loop
      declare
        s : constant string := Unix_Command_Line.Argument(i);
      begin
        if s(2) = '0'
         then res := Convert(s(3..s'last)); exit;
        end if;
      end;
    end loop;
    return res;
  end Find_Seed;

  procedure Find_and_Set_Seed is

  -- DESCRIPTION :
  --   In case the '-0' option is on, the value of the seed
  --   will be extracted from the command line, and if nonzero
  --   it will be used to set the seed.
  --   If there is no value at the command line after '-0',
  --   then a fixed constant value is used as seed.

    seed_found : constant natural32 := Find_Seed;

  begin
    if seed_found = 0
     then Standard_Random_Numbers.Set_Seed(fixed_seed);
     else Standard_Random_Numbers.Set_Seed(seed_found);
    end if;
  end Find_and_Set_Seed;

-- DISPATCHING ACCORDING TO OPTIONS :

  procedure Black_Box_Dispatcher
              ( option2,option3 : in character;
                file1,file2,file3 : in string ) is

  -- DESCRIPTION :
  --   When the first option is 'b', then this routine handles the
  --   second option and calls the appropriate main driver.

    bbprc : constant natural32 := Scan_Precision('b');
    valiprc : constant natural32 := Scan_Precision('v');
    contprc : constant natural32 := Scan_Precision('p');
    redprc : constant natural32 := Scan_Precision('d');

  begin
    case option2 is
      when 'h' | '-' => Greeting_Banners.help4blackbox;
      when 's'    => mainscal(file1,file2);
      when 'd'    => 
        if redprc = 2 then
          mainred2(file1,file2);
        elsif redprc = 4 then
          mainred4(file1,file2);
        else
          mainred(file1,file2);
        end if;
      when 'r'    =>
        case option3 is
          when 't' => bablroco(Number_of_Tasks,file1,file2);
          when others => bablroco(0,file1,file2);
        end case;
      when 'm'    =>
        case option3 is
          when 't' => babldmvc(Number_of_Tasks,file1,file2);
          when others => babldmvc(0,file1,file2);
        end case;
      when 'p' =>
        if bbprc = 2 or contprc = 2 then
          bablpoco2(file1,file2,file3);
        elsif bbprc = 4 or contprc = 4 then
          bablpoco4(file1,file2,file3);
        else
          bablpoco(file1,file2,file3);
        end if;
      when 'v' => 
        if bbprc = 2 or valiprc = 2 then
          bablvali2(file1,file2);
        elsif bbprc = 4 or valiprc = 4 then
          bablvali4(file1,file2);
        else
          bablvali(file1,file2);
        end if;
      when 'e'    => bablenum(file1,file2);
      when 't'    => 
        case option3 is
          when 'm' => babldmvc(Number_of_Tasks,file1,file2);
          when others =>
            case bbprc is
              when 2 => bablphc2(Number_of_Tasks,file1,file2);
              when 4 => bablphc4(Number_of_Tasks,file1,file2);
              when others => bablphc(Number_of_Tasks,file1,file2);
            end case;
        end case;
      when others =>
        case bbprc is
          when 2 => bablphc2(0,file1,file2);
          when 4 => bablphc4(0,file1,file2);
          when others => bablphc(0,file1,file2);
        end case;
    end case;
  end Black_Box_Dispatcher;

  procedure Component_Solver
              ( option2 : in character; file1,file2 : in string ) is

  -- DESCRIPTION :
  --   This procedure handles the option '-B' to compute the numerical
  --   irreducible decomposition in a blackbox solver.
  --   Double double and quad double precisions are applied
  --   respectively for -B2 and -B4.
  --   The second option is either '-t' for the number of tasks,
  --   or '-h' for help.

    compprc : constant natural32 := Scan_Precision('B');

  begin
    if option2 = 't' then
      case compprc is
        when 2 => compsolve2(Number_of_Tasks,file1,file2);
        when 4 => compsolve4(Number_of_Tasks,file1,file2);
        when others => compsolve(Number_of_Tasks,file1,file2);
      end case;
    elsif option2 = 'h' then
      Greeting_Banners.help4compsolve;
    else 
      case compprc is
        when 2 => compsolve2(0,file1,file2);
        when 4 => compsolve4(0,file1,file2);
        when others => compsolve(0,file1,file2);
      end case;
    end if;
  end Component_Solver;

  procedure Scaling_Dispatcher ( infile,outfile : in string ) is

  -- DESCRIPTION :
  --   This routine is called when the first option is 's'.

  begin
    put_line(welcome); put_line(scalban);
    mainscal(infile,outfile);
  end Scaling_Dispatcher;

  procedure Reduction_Dispatcher ( infile,outfile : in string ) is

  -- DESCRIPTION :
  --   This dispatcher is called when the first option is 'd'.

    redprc : constant natural32 := Scan_Precision('d');

  begin
    put_line(welcome); put_line(reduban);
    if redprc = 2 then
      mainred2(infile,outfile);
    elsif redprc = 4 then
      mainred4(infile,outfile);
    else
      mainred(infile,outfile);
    end if;
  end Reduction_Dispatcher;

  procedure Root_Count_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Calls the root couting facilities in PHCpack,
  --   when the first option was 'r'.

  begin
    case o2 is
      when 'h' | '-' => Greeting_Banners.help4rootcounts;
      when 'b' => bablroco(0,infile,outfile);
      when 't' =>
        declare
          nt : constant natural32 := Number_of_Tasks;
          ns : constant string := Convert(integer32(nt));
        begin
          put_line(welcome); put_line(rocoban & ", with " & ns & " tasks.");
          mainroco(nt,infile,outfile);
        end;
      when others =>
        put_line(welcome); put_line(rocoban & ", no multitasking.");
        mainroco(0,infile,outfile);
    end case;
  end Root_Count_Dispatcher;

  procedure Mixed_Volume_Dispatcher
              ( o2,o3 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Calls the mixed-volume computation menu, as the first option o1 = 'm'.
  --   Recognized options for o2 and o3 are
  --   't' : for multitasking polyhedral path trackers; or
  --   'b' : for blackbox construction of a random coefficient system.
  --   The options 't' and 'b' may be combined.

  begin
    case o2 is
      when 'h' | '-' => Greeting_Banners.help4mixvol;
      when 'b' =>
        case o3 is
          when 't' =>
            declare
             -- nt : constant natural := Number_of_Tasks(3);
              nt : constant natural32 := Number_of_Tasks;
            begin
              babldmvc(nt,infile,outfile);
            end;
          when others => babldmvc(0,infile,outfile);
        end case;
      when 't' =>
        declare
         -- nt : constant natural := Number_of_Tasks(2);
          nt : constant natural32 := Number_of_Tasks;
          ns : constant string := Convert(integer32(nt));
        begin
          if o3 = 'b' then
            babldmvc(nt,infile,outfile);
          else
            put_line(welcome); put_line(mvcban & ", with " & ns & " tasks.");
            mainsmvc(nt,infile,outfile);
          end if;
        end;
      when others =>
       put_line(welcome); put_line(mvcban & ", no multitasking.");
       mainsmvc(0,infile,outfile);
    end case;
  end Mixed_Volume_Dispatcher;

  procedure Continuation_Dispatcher 
              ( o2 : in character; file1,file2,file3 : in string ) is

  -- DESCRIPTION :
  --   Invokes the path trackers in PHCpack.

    contprc : constant natural32 := Scan_Precision('p');
    bbprc : constant natural32 := Scan_Precision('b');

  begin
    case o2 is
      when 'h' | '-' => Greeting_Banners.help4continuation;
      when 'b' =>
        if contprc = 2 or bbprc = 2 then
          bablpoco2(file1,file2,file3);
        elsif contprc = 4 or bbprc = 4 then
          bablpoco4(file1,file2,file3);
        else
          bablpoco(file1,file2,file3);
        end if;
      when 't'
        => declare
            -- nt : constant natural := Number_of_Tasks(2);
             nt : constant natural32 := Number_of_Tasks;
             ns : constant string := Convert(integer32(nt));
           begin
             put_line(welcome);
             put_line(pocoban & ", with " & ns & " tasks.");
             mainpoco(nt,file1,file2,contprc);
           end;
      when others
        => put_line(welcome); put_line(pocoban & ", no multitasking.");
           mainpoco(0,file1,file2,contprc);
    end case;
  end Continuation_Dispatcher;

  procedure Verification_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Validates the computed results.

    veriprc : constant natural32 := Scan_Precision('v');
    bbprc : constant natural32 := Scan_Precision('b');

  begin
    case o2 is
      when 'b' =>
        if veriprc = 2 or bbprc = 2 then
          bablvali2(infile,outfile);
        elsif veriprc = 4 or bbprc = 4 then
          bablvali4(infile,outfile);
        else
          bablvali(infile,outfile);
        end if;
      when 'h' | '-' => Greeting_Banners.help4verification;
      when others => put_line(welcome); put_line(veriban);
                     mainvali(infile,outfile);
    end case;
  end Verification_Dispatcher;

  procedure Series_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Applies Newton's method to compute power series solutions.

    prc : constant natural32 := Scan_Series_Precision('u');

  begin
    case o2 is
      when 'h' | '-' =>
        Greeting_Banners.help4series;
      when others =>
        put_line(welcome);
        case prc is
          when 1 =>
            put_line(seriesban & ", in double precision.");
            mainseries('1',infile,outfile);
          when 2 =>
            put_line(seriesban & ", with double doubles.");
            mainseries('2',infile,outfile);
          when 4 =>
            put_line(seriesban & ", with quad doubles.");
            mainseries('4',infile,outfile);
          when others =>
            put_line(seriesban & ".");
            mainseries('0',infile,outfile);
        end case;
    end case;
  end Series_Dispatcher;

  procedure Enumeration_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Invokes the numerical Schubert calculus.

  begin
    case o2 is
      when 'h' | '-' => Greeting_Banners.help4enumeration;
      when 'b'    => bablenum(infile,outfile);
      when others => put_line(welcome); put_line(enumban);
                     mainenum;
    end case;
  end Enumeration_Dispatcher;

  procedure Decomposition_Dispatcher 
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Witness point generation facilities.

  begin
    case o2 is
      when 'h' | '-' => Greeting_Banners.help4components;
      when 't' =>
        declare
          nt : constant natural32 := Number_of_Tasks;
          ns : constant string := Convert(integer32(nt));
        begin
          put_line(welcome);
          put_line(compban & " with " & ns & " tasks");
          maindeco(nt,infile,outfile);
        end;
      when others =>
        put_line(welcome); put_line(compban);
        maindeco(0,infile,outfile);
    end case;
  end Decomposition_Dispatcher;

  procedure Factorization_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Filtering and factorization capabilities, when o1 = 'f'.

    nt : natural32;

  begin
    put_line(welcome);
    case o2 is
      when 't' =>
        nt := Number_of_Tasks;
        declare
          ns : constant string := Convert(integer32(nt));
        begin
          put_line(facban & ", with " & ns & " tasks.");
        end;
      when others =>
        nt := 0;
        put_line(facban & ", no multitasking.");
    end case;
    mainfac(nt,infile,outfile);
  end Factorization_Dispatcher;

  procedure Tasking_Dispatcher
              ( o2,o3 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   For multitasking, when the first option o1 = 't'.

   -- nt : constant natural := Number_of_Tasks(1);
    nt : constant natural32 := Number_of_Tasks;
    ns : constant string := Convert(integer32(nt));
    bbprc : constant natural32 := Scan_Precision('b');
    contprc : constant natural32 := Scan_Precision('p');

  begin
    case o2 is
      when 'h' | '-' => Greeting_Banners.help4tasking;
      when 'm' =>
        case o3 is
          when 'b' => babldmvc(nt,infile,outfile);
          when others
            => put_line(welcome); put_line(mvcban & ", with " & ns & " tasks");
               mainsmvc(nt,infile,outfile);
        end case;
      when 'p' =>
        put_line(welcome); put_line(pocoban & ", with " & ns & " tasks");
        mainpoco(nt,infile,outfile,contprc);
      when 'b' =>
        case o3 is
          when 'm' => babldmvc(nt,infile,outfile);
          when others =>
            case bbprc is
              when 2 => bablphc2(nt,infile,outfile);
              when 4 => bablphc4(nt,infile,outfile);
              when others => bablphc(nt,infile,outfile);
            end case;
        end case;
      when 'c' =>
        put_line(welcome);
        put_line(compban & " with " & ns & " tasks");
        maindeco(nt,infile,outfile);
      when 'f' => 
        put_line(welcome); put_line(facban & ", with " & ns & " tasks.");
        mainfac(nt,infile,outfile);
      when 'r' =>
        put_line(welcome); put_line(rocoban & ", with " & ns & " tasks.");
        mainroco(nt,infile,outfile);
      when others =>
        mainphc(nt,infile,outfile);
    end case;
  end Tasking_Dispatcher;

  procedure Witness_Set_Intersection_Dispatcher
              ( witset_one,witset_two,logfile : in string ) is

  -- DESCRIPTION :
  --   Intersection of two witness sets using diagonal homotopies.

  begin
    if witset_one = "" or witset_two = "" or logfile = ""
     then put_line(welcome); put_line(witban);
    end if;
    mainwit(witset_one,witset_two,logfile);
  end Witness_Set_Intersection_Dispatcher;

  procedure Witness_Set_for_Hypersurface_Dispatcher
              ( polyfile,logfile : in string ) is

  -- DESCRIPTION :
  --   Scans the command line arguments for the precision
  --   and then makes a witness set for a hypersurface.

    precision : constant natural32 := Scan_Precision('l');

  begin
    if polyfile = "" or logfile = ""
     then put_line(welcome); put_line(hypban);
    end if;
    case precision is
      when 1 => mainhyp(polyfile,logfile);
      when 2 => mainhyp2(polyfile,logfile);
      when 4 => mainhyp4(polyfile,logfile);
      when others => null;
    end case;
  end Witness_Set_for_Hypersurface_Dispatcher;

  procedure Test_if_System_is_Good ( infile,outfile : in string ) is

  -- DESCRIPTION :
  --   The -g option provides a syntactical analysis of the system
  --   given on the input file.

  begin
    if infile = "" or outfile = ""  -- still interactive input needed
     then put_line(welcome); put_line(goodban);
    end if;
    maingood(infile,outfile);
  exception
    when others => raise; -- put_line("exception caught ..."); raise;
  end Test_if_System_is_Good;

  procedure General_Help ( opt : in character ) is

  -- DESCRIPTION :
  --   Writes general help about an option to screen.

  begin
    case opt is
      when '0' => Greeting_Banners.help4setseed;
      when 'a' => Greeting_Banners.help4eqnbyeqn;
      when 'b' => Greeting_Banners.help4blackbox;
      when 'B' => Greeting_Banners.help4compsolve;
      when 'c' => Greeting_Banners.help4components;
      when 'd' => Greeting_Banners.help4reduction;
      when 'e' => Greeting_Banners.help4enumeration;
      when 'f' => Greeting_Banners.help4factor;
      when 'g' => Greeting_Banners.help4goodformat;
      when 'h' | '-' => Greeting_Banners.help4help;
      when 'j' => Greeting_Banners.help4adepath;
      when 'k' => Greeting_Banners.help4feedback;
      when 'l' => Greeting_Banners.help4hypersurface;
      when 'm' => Greeting_Banners.help4mixvol;
      when 'o' => Greeting_Banners.help4symbols;
      when 'p' => Greeting_Banners.help4continuation;
      when 'q' => Greeting_Banners.help4jumpstart;
      when 'r' => Greeting_Banners.help4rootcounts;
      when 's' => Greeting_Banners.help4scaling;
      when 't' => Greeting_Banners.help4tasking;
      when 'u' => Greeting_Banners.help4series;
      when 'v' => Greeting_Banners.help4verification;
      when 'w' => Greeting_Banners.help4witsetinsect;
      when 'x' => Greeting_Banners.help4pythondict;
      when 'y' => Greeting_Banners.help4sampler;
      when 'z' => Greeting_Banners.help4mapleform;
      when others => Greeting_Banners.show_help;
    end case;
  end General_Help;

  procedure General_Dispatcher
              ( o1,o2,o3 : in character; f1,f2,f3 : in string ) is

  -- DESCRIPTION :
  --   This is the general dispatcher for the first option o1.

  begin
    case o1 is
      when 'a'    => if o2 = 'h' or o2 = '-' then
                       Greeting_Banners.help4eqnbyeqn;
                     else
                       put_line(welcome); put_line(slvban);
                       mainsolve(f1,f2);
                     end if;
      when 'b'    => Black_Box_Dispatcher(o2,o3,f1,f2,f3);
      when 'B'    => Component_Solver(o2,f1,f2);
      when 'c'    => Decomposition_Dispatcher(o2,f1,f2);
      when 'd'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4reduction;
                      else Reduction_Dispatcher(f1,f2);
                     end if;
      when 'e'    => Enumeration_Dispatcher(o2,f1,f2);
      when 'f'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4factor;
                      else Factorization_Dispatcher(o2,f1,f2);
                     end if;
      when 'g'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4goodformat;
                      else Test_if_System_is_Good(f1,f2);
                     end if;
      when 'h'    => General_Help(o2);
      when 'j'    => if o2 = 'h' or o2 = '-' then
                       Greeting_Banners.help4adepath;
                     else
                       put_line(welcome); put_line(adepban);
                       mainadep(f1,f2,f3);
                     end if;
      when 'k'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4feedback;
                      else mainfeed(f1,f2);
                     end if;
      when 'l'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4hypersurface;
                      else Witness_Set_for_Hypersurface_Dispatcher(f1,f2);
                     end if;
      when 'm'    => Mixed_Volume_Dispatcher(o2,o3,f1,f2);
      when 'o'    => if o2 = 'h' or o2 = '-' then
                       Greeting_Banners.help4symbols;
                     else
                       put_line(welcome); put_line(symbban);
                       mainsymb(f1,f2);
                     end if;
      when 'p'    => Continuation_Dispatcher(o2,f1,f2,f3);
      when 'q'    => if o2 = 'h' or o2 = '-' then
                       Greeting_Banners.help4jumpstart;
                     else
                       put_line(welcome); put_line(trackban);
                       maintrack(f1,f2,f3);
                     end if;
      when 'r'    => Root_Count_Dispatcher(o2,f1,f2);
      when 's'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4scaling;
                      else Scaling_Dispatcher(f1,f2);
                     end if;
      when 't'    => Tasking_Dispatcher(o2,o3,f1,f2);
      when 'u'    => Series_Dispatcher(o2,f1,f2);
      when 'v'    => Verification_Dispatcher(o2,f1,f2);
      when 'w'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4witsetinsect;
                      else Witness_Set_Intersection_Dispatcher(f1,f2,f3);
                     end if;
      when 'x'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4pythondict;
                      else maindict(f1,f2);
                     end if;
      when 'y'    => if o2 = 'h' or o2 = '-' then
                       Greeting_Banners.help4sampler;
                     else
                       put_line(welcome); put_line(samban);
                       mainsam(f1,f2);
                     end if;
      when 'z'    => if o2 = 'h' or o2 = '-'
                      then Greeting_Banners.help4mapleform;
                      else mainzip(f1,f2);
                     end if;
      when others => put_line(welcome); mainphc(0,f1,f2);
    end case;
  exception
    when others => raise; -- put_line("exception in dispatch..."); raise;
  end General_Dispatcher;

  procedure Help_Version_License is

  -- DESCRIPTION :
  --   Deals with the options --help, --version, and --license.
  --   Writes the version string to file if file1 /= ""
  --   or to screen if file = "".

    arg : constant string := Unix_Command_Line.Argument(1);

  begin
   -- put_line("arg = " & arg);
    if arg = "--help" then
     -- put("argc = "); put(natural32(argc),1); new_line;
      if argc = 1 then
        Greeting_Banners.show_help;
      else
        declare
          opt2 : constant string := Unix_Command_Line.Argument(2);
        begin
         -- put("The second option : "); put(opt2); new_line;
          General_Help(opt2(2));
        end;
      end if;
    elsif arg = "--version" then
      if argc = 1 then
        put_line(Greeting_Banners.Version);
      else
        declare
          name : constant string := Unix_Command_Line.Argument(2);
          file : file_type;
        begin
         -- put_line("name = " & name);
          create(file,out_file,name);
          put_line(file,Greeting_Banners.Version);
          close(file);
        end;
      end if;
    elsif arg = "--license" then
      put_line("PHCpack is free and open source software.");
      put_line("You can redistribute the code and/or modify it under");
      put_line("the GNU General Pulic License as published by");
      put_line("the Free Software Foundation.");
    elsif arg = "--cite" then
      Greeting_Banners.How_to_Cite;
    else
      put_line(arg & " is not recognized.");
    end if;
  end Help_Version_License;

  procedure Main is

  -- DESCRIPTION :
  --   Parses the arguments on the command line, separating the arguments
  --   into options (at most 4) and file names (at most 3).

    posi : integer32 := 0;
    option1,option2,option3,option4 : character;
    file1 : constant string := Read_Argument(1);
    file2 : constant string := Read_Argument(2);
    file3 : constant string := Read_Argument(3);

  begin
    Read_Next_Option(posi,options,option1);
    if option1 = '-'
     then Help_Version_License; return;
    end if;
    if option1 /= ' ' then
      Read_Next_Option(posi,options,option2);
      if option2 /= ' ' then
        Read_Next_Option(posi,options,option3);        
        if option3 /= ' '
         then Read_Next_Option(posi,options,option4);        
         else option4 := ' ';
        end if;
      else
        option3 := ' ';
        option4 := ' ';
      end if;
    else
      option2 := ' ';
      option3 := ' ';
      option4 := ' ';
    end if;
    if (file1 /= "") and then (file1 = file2) then
      new_line;
      put_line("The first two file names are the same.");
      put_line("Will ignore the second file name...");
      if option1 = '0' then
        if option2 = 'h' or option2 = '-' then
          Greeting_Banners.help4setseed;
        else
          Find_and_Set_Seed;
          General_Dispatcher(option2,option3,option4,file1,"",file3);
        end if;
      elsif option2 = '0' then
        if option1 = 'h' or option1 = '-' then
          Greeting_Banners.help4setseed;
        else
          Find_and_Set_Seed;
          General_Dispatcher(option1,option3,option4,file1,"",file3);
        end if;
      elsif option3 = '0' then
        Find_and_Set_Seed;
        General_Dispatcher(option1,option2,option4,file1,"",file3);
      elsif option4 = '0' then
        Find_and_Set_Seed;
        General_Dispatcher(option1,option2,option3,file1,"",file3);
      else
        General_Dispatcher(option1,option2,option3,file1,"",file3);
      end if;
    else
      if option1 = '0' then
        if option2 = 'h' or option2 = '-' then
          Greeting_Banners.help4setseed;
        else
          Find_and_Set_Seed;
          General_Dispatcher(option2,option3,option4,file1,file2,file3);
        end if;
      elsif option2 = '0' then
        if option1 = 'h' or option1 = '-' then
          Greeting_banners.help4setseed;
        else
          Find_and_Set_Seed;
          General_Dispatcher(option1,option3,option4,file1,file2,file3);
        end if;
      elsif option3 = '0' then
        Find_and_Set_Seed;
        General_Dispatcher(option1,option2,option4,file1,file2,file3);
      elsif option4 = '0' then
        Find_and_Set_Seed;
        General_Dispatcher(option1,option2,option3,file1,file2,file3);
      else
        General_Dispatcher(option1,option2,option3,file1,file2,file3);
      end if;
    end if;
  exception
    when others => raise; -- put_line("exception in main ..."); raise;
  end Main;

begin
  Main;
exception
  when others =>
    put_line("Oops, an unhandled exception occurred.");
    Write_Seed_Number(standard_output);
    put_line("Use the seed number to reproduce the error.");
    raise; -- put_line("remain silent ..."); return;
end Dispatch;
