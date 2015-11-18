with Standard_Natural_Numbers_io;
 use Standard_Natural_Numbers_io;

with text_io;                            use text_io;
with Unix_Command_Line;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Random_Numbers; -- to handle zero seed
with Greeting_Banners;
with mainscal,mainred;        -- scaling and reduction
with mainroco,bablroco;       -- general root counting
with mainsmvc,babldmvc;       -- mixed-volume computation
with maintrack;               -- path tracking
with mainpoco,bablpoco;       -- polynomial continuation
with bablpoco2,bablpoco4;     -- double double and quad double continuation
with mainphc,bablphc;         -- main phc driver + blackbox
with bablphc2,bablphc4;       -- blackbox in double double and quad double
with mainvali,bablvali;       -- validation tool
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
-- NOTE (added for pieri_solver.ali) :
with Interfaces.C;
with Complex_Polynomial_Matrices;        use Complex_Polynomial_Matrices;
with Complex_Polynomial_Matrices_io;     use Complex_Polynomial_Matrices_io;
with Verify_Solution_Maps;
with C_to_Ada_Arrays;                    use C_to_Ada_Arrays;

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
  trackban : constant string := Greeting_Banners.trackban;
  valiban : constant string := Greeting_Banners.valiban;
  witban  : constant string := Greeting_Banners.witban;

-- AVAILABLE OPTIONS :

  options : constant string := "0asdpqmrvbekcxyzftwlgo";
  -- 0 : zero seed for repeatable runs
  -- a : solve => equation-by-equation solver
  -- b : batch or black box processing
  -- c : comp => numerical irreducible decomposition
  -- d : redu => reduction w.r.t. the total degree
  -- e : enum => numerical Schubert calculus
  -- f : fac  => factor pure dimensional solution set into irreducibles
  -- g : good => check if the input is a system in the valid format
  -- k : feba => dynamic output feedback to control linear systems
  -- l : hyp  => witness set for hypersurface cutting with random line
  -- m : mvc  => mixed-volume computation
  -- o : symbOls -> get symbols using as variables in a system
  -- p : poco => polynomial continuation
  -- q : track => tracking solution paths
  -- r : roco => root counting methods
  -- s : scal => scaling of a polynomial system
  -- t : task => use multitasking
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
    contprc : constant natural32 := Scan_Precision('p');

  begin
   -- put("The blackbox precision : "); put(bbprc,1); new_line;
    case option2 is
      when 's'    => mainscal(file1,file2);
      when 'd'    => mainred(file1,file2);
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
      when 'v'    => bablvali(file1,file2);
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

  begin
    put_line(welcome); put_line(reduban);
    mainred(infile,outfile);
  end Reduction_Dispatcher;

  procedure Root_Count_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Calls the root couting facilities in PHCpack,
  --   when the first option was 'r'.

  begin
    case o2 is
      when 'b'    => bablroco(0,infile,outfile);
     -- when 't'    => bablroco(Number_of_Tasks(2),infile,outfile);
      when 't'    => bablroco(Number_of_Tasks,infile,outfile);
      when others => put_line(welcome); put_line(rocoban);
                     mainroco(infile,outfile);
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
            put_line(welcome); put_line(mvcban & ", with " & ns & " tasks");
            mainsmvc(nt,infile,outfile);
          end if;
        end;
      when others =>
       put_line(welcome); put_line(mvcban & ", no multitasking");
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
      when 'b' =>
        put("The value of contprc : "); put(contprc,1); new_line;
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
             put_line(welcome); put_line(pocoban & ", with " & ns & " tasks");
             mainpoco(nt,file1,file2,contprc);
           end;
      when others
        => put_line(welcome); put_line(pocoban & ", no multitasking");
           mainpoco(0,file1,file2,contprc);
    end case;
  end Continuation_Dispatcher;

  procedure Validation_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Validates the computed results.

  begin
    case o2 is
      when 'b'    => bablvali(infile,outfile);
      when others => put_line(welcome); put_line(valiban);
                     mainvali(infile,outfile);
    end case;
  end Validation_Dispatcher;

  procedure Enumeration_Dispatcher
              ( o2 : in character; infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Invokes the numerical Schubert calculus.

  begin
    case o2 is
      when 'b'    => bablenum(infile,outfile);
      when others => put_line(welcome); put_line(enumban);
                     mainenum;
    end case;
  end Enumeration_Dispatcher;

  procedure Decomposition_Dispatcher ( infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Witness point generation facilities.

  begin
    put_line(welcome); put_line(compban);
    maindeco(infile,outfile);
  end Decomposition_Dispatcher;

  procedure Factorization_Dispatcher ( infile,outfile : in string ) is

  -- DESCRIPTION :
  --   Filtering and factorization capabilities.

  begin
    put_line(welcome); put_line(facban);
    mainfac(infile,outfile);
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
      when others => mainphc(nt,infile,outfile);
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

  procedure General_Dispatcher
              ( o1,o2,o3 : in character; f1,f2,f3 : in string ) is

  -- DESCRIPTION :
  --   This is the general dispatcher for the first option o1.

  begin
    case o1 is
      when 'a'    => put_line(welcome); put_line(slvban);
                     mainsolve(f1,f2);
      when 'b'    => Black_Box_Dispatcher(o2,o3,f1,f2,f3);
      when 'c'    => Decomposition_Dispatcher(f1,f2);
      when 'd'    => Reduction_Dispatcher(f1,f2);
      when 'e'    => Enumeration_Dispatcher(o2,f1,f2);
      when 'g'    => Test_if_System_is_Good(f1,f2);
      when 'm'    => Mixed_Volume_Dispatcher(o2,o3,f1,f2);
      when 'k'    => mainfeed(f1,f2);
      when 'l'    => Witness_Set_for_Hypersurface_Dispatcher(f1,f2);
      when 's'    => Scaling_Dispatcher(f1,f2);
      when 'f'    => Factorization_Dispatcher(f1,f2);
      when 'o'    => put_line(welcome); put_line(symbban);
                     mainsymb(f1,f2);
      when 'p'    => Continuation_Dispatcher(o2,f1,f2,f3);
      when 'q'    => put_line(welcome); put_line(trackban);
                     maintrack(f1,f2,f3);
      when 't'    => Tasking_Dispatcher(o2,o3,f1,f2);
      when 'v'    => Validation_Dispatcher(o2,f1,f2);
      when 'r'    => Root_Count_Dispatcher(o2,f1,f2);
      when 'w'    => Witness_Set_Intersection_Dispatcher(f1,f2,f3);
      when 'x'    => maindict(f1,f2);
      when 'y'    => put_line(welcome); put_line(samban); mainsam(f1,f2);
      when 'z'    => mainzip(f1,f2);
      when others => put_line(welcome); mainphc(0,f1,f2);
    end case;
  exception
    when others => raise; -- put_line("exception in dispatch..."); raise;
  end General_Dispatcher;

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
        Find_and_Set_Seed;
        General_Dispatcher(option2,option3,option4,file1,"",file3);
      elsif option2 = '0' then
        Find_and_Set_Seed;
        General_Dispatcher(option1,option3,option4,file1,"",file3);
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
        Find_and_Set_Seed;
        General_Dispatcher(option2,option3,option4,file1,file2,file3);
      elsif option2 = '0' then
        Find_and_Set_Seed;
        General_Dispatcher(option1,option3,option4,file1,file2,file3);
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
  when others => raise; -- put_line("remain silent ..."); return;
end Dispatch;
