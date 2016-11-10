with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Characters_and_Numbers;           use Characters_and_Numbers;
with Brackets_io;                      use Brackets_io;

package body Bracket_Monomials_io is

  procedure Skip_Spaces ( file : in file_type; ch : out character ) is

  -- DESCRIPTION :
  --   As long as the file is not at the end and ch is a space,
  --   characters are read from file.  Returns with the first nonspace
  --   character or when the end of file is reached.

  begin
    if not end_of_file(file) then
      get(file,ch);
      while ((ch = ' ') and not end_of_file(file)) loop
        get(file,ch);
      end loop;
    end if;
  end Skip_Spaces;

  procedure Read_Number ( file : in file_type; n : out natural32;
                          ch : in out character ) is

  -- DESCRIPTION :
  --   Given in ch the first character of a number,
  --   continues reading characters till the first character that
  --   is not a single digit is encountered. 

  -- ON RETURN :
  --   n         number read from file;
  --   ch        first character that is not a single digit.

  begin
    n := Convert(ch);
    loop
      get(file,ch);
      exit when (Convert(ch) = 10);
      n := 10*n + Convert(ch);
    end loop;
  end Read_Number;

  procedure get ( file : in file_type; k,nb : in natural32;
                  ch : in character; b : out Link_to_Bracket ) is

  -- DESCRIPTION :
  --   The start of k-th number in the bracket b starts with the
  --   single digit number nb.

    new_ch : character;
    new_nb : natural32;

  begin
    if ch = ']' then
      b := new Bracket(1..integer32(k));
    elsif ch = ' ' then
      new_ch := ch;
      Skip_Spaces(file,new_ch);
      if new_ch = ']' then
        b := new Bracket(1..integer32(k));
      else
        Read_Number(file,new_nb,new_ch);
        get(file,k+1,new_nb,new_ch,b);
      end if;
    end if;
    b(integer32(k)) := nb;
  end get;

-- TARGET PROCEDURES :

  procedure get ( b : out Link_to_Bracket ) is
  begin
    get(standard_input,b);
  end get;

  procedure get ( file : in file_type; b : out Link_to_Bracket ) is

    ch : character;
    nb : natural32;

  begin
    Skip_Spaces(file,ch);
    if ch = '[' then
      Skip_Spaces(file,ch);
      if ch /= ']' then
        Read_Number(file,nb,ch);
        get(file,1,nb,ch,b);
      end if;
    end if;
  end get;

  procedure get ( bm : out Bracket_Monomial ) is
  begin
    get(standard_input,bm);
  end get;

  procedure get ( file : in file_type; bm : out Bracket_Monomial ) is

    b : Link_to_Bracket;
    ch : character;
    e : natural32;

  begin
    get(file,b);
    if b /= null then
     -- bm := Create(b.all);
      loop
        get(file,ch);
        if ch = ' '
         then Skip_Spaces(file,ch);
        end if;
        if (ch = ';')
         then Append(bm,b.all);
        end if;
        exit when (ch = ';');
        if ch = '^' then
          get(file,ch);
          Read_Number(file,e,ch);
          for i in 1..e loop
            Append(bm,b.all);
          end loop;
        else
          Append(bm,b.all);
        end if;
        exit when (ch = ';');
        Clear(b);
        get(file,b);
        exit when (b = null);
       -- Append(bm,b.all);
      end loop;
    end if;
  end get;

  procedure put ( bm : in Bracket_Monomial ) is
  begin
    put(Standard_Output,bm);
  end put;

  procedure put ( file : in file_type; bm : in Bracket_Monomial ) is

    procedure Write_Bracket ( b : in Bracket; cont : out boolean ) is
    begin
      put(file,b);
      cont := true;
    end Write_Bracket;
    procedure Write_Brackets is new Enumerate_Brackets(Write_Bracket);

  begin
    Write_Brackets(bm);
  end put;

end Bracket_Monomials_io;
