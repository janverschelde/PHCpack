with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;

package body Test_Deca_Doubles is

  procedure write ( x : deca_double ) is
  begin
    put("  thumb right  : "); put(thumb_right(x),2,17,3); new_line;
    put("  index right  : "); put(index_right(x),2,17,3); new_line;
    put("  middle right : "); put(middle_right(x),2,17,3); new_line;
    put("  ring right   : "); put(ring_right(x),2,17,3); new_line;
    put("  pink right   : "); put(pink_right(x),2,17,3); new_line;
    put("  thumb left   : "); put(thumb_left(x),2,17,3); new_line;
    put("  index left   : "); put(index_left(x),2,17,3); new_line;
    put("  middle left  : "); put(middle_left(x),2,17,3); new_line;
    put("  ring left    : "); put(ring_left(x),2,17,3); new_line;
    put("  pink left    : "); put(pink_left(x),2,17,3); new_line;
  end Write;

  function random return deca_double is

    res : deca_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..10 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end random;

  procedure Test_Addition_and_Subtraction is

    x : constant deca_double := random;
    y : constant deca_double := random;
    z : constant deca_double := x + y;
    v : constant deca_double := z - y;

  begin
   put_line("All parts of a random deca double x :"); Write(x);
   put_line("All parts of a random deca double y :"); Write(y);
   put_line("All parts of x + y :"); Write(z);
   put_line("All parts of (x + y) - y :"); Write(v);
  end Test_Addition_and_Subtraction;

  procedure Test_Multiplication_and_Division is

    x : constant deca_double := random;
    y : constant deca_double := random;
    z : constant deca_double := x*y;
    v : constant deca_double := z/y;

  begin
   put_line("All parts of a random deca double x :"); Write(x);
   put_line("All parts of a random deca double y :"); Write(y);
   put_line("All parts of x * y :"); Write(z);
   put_line("All parts of (x * y) / y :"); Write(v);
  end Test_Multiplication_and_Division;

  procedure Test_Read is

  --   >>> from sympy import evalf, sqrt
  --   >>> s2 = sqrt(2).evalf(160)
  --   >>> s2
  --   1.414213562373095048801688724209698078569671875376948073176679737
  --   99073247846210703885038753432764157273501384623091229702492483605
  --   5850737212644121497099935831413
    
    sqrt2 : constant string
      := "1.414213562373095048801688724209698078569671875376948073176679737"
       & "99073247846210703885038753432764157273501384623091229702492483605"
       & "5850737212644121497099935831413";
    x,r : deca_double;
    two : constant deca_double := create(2.0);
    fail : boolean;

  begin
    read(sqrt2,x,fail);
    new_line;
    if fail then
      put_line("The read procedure reports failure!");
    else
      put_line("All parts of the sqrt(2) read :"); Write(x);
      r := x*x - 2.0;
      put_line("All parts of x*x - 2.0 : "); Write(r);
      r := x*x - two;
      put_line("All parts of x*x - two : "); Write(r);
    end if;
  end Test_Read;

  procedure Test_Write is

    dfx : double_float := 1.0E-16;
    dax : deca_double := create(dfx);
    ans : character;

  begin
    loop
      new_line;
      put("A double float : "); put(dfx); new_line;
      put_line("All words in the corresponding deca double :");
      write(dax);
      put("Written as a deca double : "); put(dax); new_line;
      put("Give your own double float ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a double float : "); get(dfx);
      dax := create(dfx);
    end loop;
  end Test_Write;

  procedure Log10log2exp1_doubles is

  -- >>> from sympy import log, exp, evalf
  -- >>> log(10).evalf(160)
  -- 2.3025850929940456840179914546843642076011014886287729760333
  -- 279009675726096773524802359972050895982983419677840422862486
  -- 33409525465082806756666287369098781689483
  -- >>> log(2).evalf(160)
  -- 0.6931471805599453094172321214581765680755001343602552541206
  -- 800094933936219696947156058633269964186875420014810205706857
  -- 336855202357581305570326707516350759619307
  -- >>> exp(1).evalf(160)
  -- 2.7182818284590452353602874713526624977572470936999595749669
  -- 676277240766303535475945713821785251664274274663919320030599
  -- 21817413596629043572900334295260595630738

    log10 : constant string
          := "2.3025850929940456840179914546843642076011014886287729760333"
           & "279009675726096773524802359972050895982983419677840422862486"
           & "33409525465082806756666287369098781689483";
    log2 : constant string
         := "0.6931471805599453094172321214581765680755001343602552541206"
          & "800094933936219696947156058633269964186875420014810205706857"
          & "336855202357581305570326707516350759619307";
    exp1 : constant string
         := "2.7182818284590452353602874713526624977572470936999595749669"
          & "676277240766303535475945713821785251664274274663919320030599"
          & "21817413596629043572900334295260595630738";

    x : deca_double;
    fail : boolean;

  begin
    read(log10,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the log(10) read :"); Write(x);
    end if;
    read(log2,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the log(2) read :"); Write(x);
    end if;
    read(exp1,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the exp(1) read :"); Write(x);
    end if;
  end Log10log2exp1_doubles;

  procedure inverse_factorials is

  -- >>> from sympy import factorial, evalf
  -- >>> f = [1/factorial(k) for k in range(3,15+3)]
  -- >>> for x in f: print(x.evalf(160))

    f0 : constant string
       := "0.1666666666666666666666666666666666666666666666666666666666"
        & "666666666666666666666666666666666666666666666666666666666666"
        & "666666666666666666666666666666666666666667";
    f1 : constant string
       := "0.0416666666666666666666666666666666666666666666666666666666"
        & "666666666666666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666667";
    f2 : constant string
       := "0.00833333333333333333333333333333333333333333333333333333333"
        & "3333333333333333333333333333333333333333333333333333333333333"
        & "333333333333333333333333333333333333333333";
    f3 : constant string
       := "0.00138888888888888888888888888888888888888888888888888888888"
        & "8888888888888888888888888888888888888888888888888888888888888"
        & "888888888888888888888888888888888888888889";
    f4 : constant string
       := "0.00019841269841269841269841269841269841269841269841269841269"
        & "8412698412698412698412698412698412698412698412698412698412698"
        & "4126984126984126984126984126984126984126984";
    f5 : constant string
       := "0.00002480158730158730158730158730158730158730158730158730158"
        & "7301587301587301587301587301587301587301587301587301587301587"
        & "30158730158730158730158730158730158730158730";
    f6 : constant string
       := "0.00000275573192239858906525573192239858906525573192239858906"
        & "5255731922398589065255731922398589065255731922398589065255731"
        & "922398589065255731922398589065255731922398589";
    f7 : constant string
       := "0.00000027557319223985890652557319223985890652557319223985890"
        & "6525573192239858906525573192239858906525573192239858906525573"
        & "1922398589065255731922398589065255731922398589";
    f8 : constant string
       := "0.00000002505210838544171877505210838544171877505210838544171"
        & "8775052108385441718775052108385441718775052108385441718775052"
        & "10838544171877505210838544171877505210838544172";
    f9 : constant string
       := "0.00000000208767569878680989792100903212014323125434236545347"
        & "6564587675698786809897921009032120143231254342365453476564587"
        & "675698786809897921009032120143231254342365453477";
    f10 : constant string
       := "0.00000000016059043836821614599392377170154947932725710503488"
        & "2812660590438368216145993923771701549479327257105034882812660"
        & "5904383682161459939237717015494793272571050348828";
    f11 : constant string
       := "0.00000000001147074559772972471385169797868210566623265035963"
        & "4486618613602740586867570994555121539248523375507502491629475"
        & "75645988344401042813741226439639138051836464534877";
    f12 : constant string
       := "0.00000000000076471637318198164759011319857880704441551002397"
        & "5632441240906849372457838066303674769283234891700500166108631"
        & "7170973255629340285424941509597594253678909763565848";
    f13 : constant string
       := "0.00000000000004779477332387385297438207491117544027596937649"
        & "8477027577556678085778614879143979673080202180731281260381789"
        & "48231858284768337678390588443498496408549318602228655";
    f14 : constant string
       := "0.00000000000000281145725434552076319894558301032001623349273"
        & "5204531033973922240339918522302587039592953069454781250610693"
        & "498959916638099022163759169672646174357970187413075679";

    x : deca_double;
    fail : boolean;
    ans : character;

  begin
    read(f0,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(0) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f1,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(1) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f2,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(2) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f3,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(3) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f4,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(4) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f5,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(5) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f6,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(6) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f7,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(7) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f8,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(8) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f9,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(9) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f10,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(10) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f11,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(11) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f12,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(12) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f13,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(13) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f14,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(14) read :"); Write(x);
    end if;
  end inverse_factorials;

  procedure Test_io is

    x,y : deca_double;
    ans : character;
 
  begin
    new_line;
    loop
      put("Give x : "); get(x);
      put(" --> x : "); put(x); new_line;
      put("Give pair x y : "); get(x,y); 
      put(" --> x : "); put(x); new_line;
      put(" --> y : "); put(y); new_line;
      put("More tests ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_io;

  procedure Test_sqrt2 is

    n,x,y,z,e,a : deca_double;
    max_steps : constant natural32 := 9;
    sqrt2 : constant string
      := "1.414213562373095048801688724209698078569671875376948073176679737"
       & "99073247846210703885038753432764157273501384623091229702492483605"
       & "5850737212644121497099935831413";
    fail : boolean;

  begin
    n := Create(2.0); Copy(n,x);
    new_line;
    put_line("running Newton's method for sqrt(2) ...");
    read(sqrt2,y,fail);
    if fail
     then put_line("reading value for sqrt2 from string failed!");
    end if;
    put(" sqrt2: "); put(y); new_line;
    put("step 0: "); put(x); new_line;
    for i in 1..max_steps loop
      z := x*x;
      z := z+n;
      z := z/x;
      z := 0.5*z;
      put("step "); put(i,1); put(": "); put(z); new_line;
      copy(z,x);
      e := x - y;
      a := abs(e);
      put("  error : "); put(a,3); new_line;
    end loop;
  end Test_sqrt2;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing deca double arithmetic ...");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication and division");
    put_line("  3. test reading from a string");
    put_line("  4. test writing a double as deca double");
    put_line("  5. write 10 leading doubles of log(10), log(2), exp(1)");
    put_line("  6. write 10 leading doubles for inverse factorials");
    put_line("  7. input and output");
    put_line("  8. Newton's method for sqrt(2)");
    put("Type 1, 2, 3, 4, 5, 6, 7, or 8 to select a test : ");
    Ask_Alternative(ans,"12345678");
    case ans is
      when '1' => Test_Addition_and_Subtraction;
      when '2' => Test_Multiplication_and_Division;
      when '3' => Test_Read;
      when '4' => Test_Write;
      when '5' => Log10log2exp1_doubles;
      when '6' => inverse_factorials;
      when '7' => Test_io;
      when '8' => Test_sqrt2;
      when others => null;
    end case;
  end Main;

end Test_Deca_Doubles;
