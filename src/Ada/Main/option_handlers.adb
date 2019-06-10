with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Greeting_Banners;
with Actions_and_Options;
with mainfac;
with maindeco;
with mainenum,bablenum;

package body Option_Handlers is

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

-- THE HANDLERS :

  procedure General_Help ( opt : in character ) is
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

  procedure Enumeration_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string ) is

    hpos : constant integer32 := Actions_and_Options.Position(opts,'h');
    bpos : constant integer32 := Actions_and_Options.Position(opts,'b');

  begin
    if hpos >= integer32(opts'first) then
      Greeting_Banners.help4enumeration;
    elsif bpos >= integer32(opts'first) then
      bablenum(infile,outfile);
    else
      put_line(welcome); put_line(enumban);
      mainenum;
    end if;
  end Enumeration_Handler;

  procedure Decomposition_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string ) is

    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    pos : constant integer32 := Actions_and_Options.Position(opts,'h');

  begin
    if pos >= integer32(opts'first) then
      Greeting_Banners.help4components;
    else
      if nt > 0 then
        declare
          ns : constant string := Convert(integer32(nt));
        begin
          put_line(welcome);
          put_line(compban & " with " & ns & " tasks");
          maindeco(nt,infile,outfile);
        end;
      else
        put_line(welcome); put_line(compban);
        maindeco(0,infile,outfile);
      end if;
    end if;
  end Decomposition_Handler;

  procedure Factorization_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string ) is

    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    pos : constant integer32 := Actions_and_Options.Position(opts,'h');

  begin
    if pos >= integer32(opts'first) then
      Greeting_Banners.help4components;
    else
      put_line(welcome);
      if nt = 0 then
        put_line(facban & ", no multitasking.");
      else
        declare
          ns : constant string := Convert(integer32(nt));
        begin
          put_line(facban & ", with " & ns & " tasks.");
        end;
      end if;
      mainfac(nt,infile,outfile);
    end if;
  end Factorization_Handler;

end Option_Handlers;
