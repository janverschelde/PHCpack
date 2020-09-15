with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Deca_Double_Numbers_io;             use Deca_Double_Numbers_io;
with Deca_Double_Constants;

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
    -- >>> from sympy import factorial, evalf
    -- >>> f = [1/factorial(k) for k in range(17,17+24)]
    -- >>> for x in f: print(x.evalf(160))
    -- the first number is f14, as a check ...
    f15 : constant string
       := "0.00000000000000015619206968586226462216364350057333423519404"
        & "0844696168554106791129995473461254835532941837191932291700594"
        & "0832755509243388345646532872040358985754427881896153155";
    f16 : constant string
       := "0.00000000000000000822063524662432971695598123687228074922073"
        & "8991826114134426673217368182813750254501733780904838541668452"
        & "320172397417070464977087015116001889398707515167874490291";
    f17 : constant string
       := "0.00000000000000000041103176233121648584779906184361403746103"
        & "6949591305706721333660868409140687512725086689045241927083422"
        & "6160086198708535232488543507558000944699353757583937245145";
    f18 : constant string
       := "0.00000000000000000001957294106339126123084757437350543035528"
        & "7473790062176510539698136590911461310129766032811678187003972"
        & "50552421999385016777375496908360952830809216075039970116736";
    f19 : constant string
       := "0.00000000000000000000088967913924505732867488974425024683433"
        & "1248808639189841388168097117768702786824080274218712644863816"
        & "9320692827269931894442615895038004331049132800341090773257891";
    f20 : constant string
       := "0.00000000000000000000003868170170630684037716911931522812323"
        & "1793426462573471364702960744250813164644525229313857071515818"
        & "12748127316204318214975050389146958404803970782756995988372996";
    f21 : constant string
       := "0.00000000000000000000000161173757109611834904871330480117180"
        & "1324726102607227973529290031010450548526855217888077377979825"
        & "755311719715085132589572937662144566002001654492815414995155415";
    f22 : constant string
       := "0.00000000000000000000000006446950284384473396194853219204687"
        & "2052989044104289118941171601240418021941074208715523095119193"
        & "03021246878860340530358291750648578264008006617971261659980621660";
    f23 : constant string
       := "0.00000000000000000000000000247959626322479746007494354584795"
        & "6617422655542472658420814292355400693151579777258289349812276"
        & "655008171876484746357830112211787914716926156391527408330761777561";
    f24 : constant string
       := "0.000000000000000000000000000091836898637955461484257168364739133"
        & "97861687194343179336349230945928493153999175030701295601024648178"
        & "414357350912436407823006621906358985764413064475299117694672";
    f25 : constant string
       := "0.000000000000000000000000000003279889237069837910152041727312111"
        & "92780774542655113547726758248068874755499970536810760557179451720"
        & "6576556196754441574222502364966556780630147523026892542033811";
    f26 : constant string
       := "0.000000000000000000000000000000113099628864477169315587645769383"
        & "16992440501470865984404370974071340508810343811614164157144119024"
        & "85026398688536014335938793918953985096769016387250652600701314";
    f27 : constant string
       := "0.000000000000000000000000000000003769987628815905643852921525646"
        & "10566414683382362199480145699135711350293678127053805471904803967"
        & "4950087995628453381119795979729846616989230054624168842002337714";
    f28 : constant string
       := "0.0000000000000000000000000000000001216125041553517949629974685692"
        & "292149724785104394191871437739147455968689284280818727328725174088"
        & "693576772783372058425740638622531166770719372459409303871721843";
    f29 : constant string
       := "0.0000000000000000000000000000000000038003907548547435925936708927"
        & "884129678899534512318495982429348357999021540133775585229022661690"
        & "27167427414948037682580439495695409896158498038935654074599130760";
    f30 : constant string
       := "0.0000000000000000000000000000000000001151633562077195028058688149"
        & "329822111481804076130863514619071162363606713337387138946334020051"
        & "220353765883317587176539527119907699968532878193616864871090645685";
    f31 : constant string
       := "0.00000000000000000000000000000000000000338715753552116184723143"
        & "5733323006210240600223914304454761974006951784450992315114548041"
        & "2354447657463702450517269898221385879638234368614064518143084443"
        & "84252015";
    f32 : constant string
       := "0.00000000000000000000000000000000000000009677592958631890992089"
        & "8163809228748864017149254694412993199257341479555742637574701372"
        & "6067269933070391498586207711377753882275378124817544700518373841"
        & "2526434328";
    f33 : constant string
       := "0.00000000000000000000000000000000000000000268822026628663638669"
        & "1615661367465246222698590408178138699979370596654326184377075038"
        & "1279646387029733097182950214204937607840982725689376241681065940"
        & "03479565091";
    f34 : constant string
       := "0.00000000000000000000000000000000000000000007265460179153071315"
        & "3827450307228790438451313254275084829729172178287954761739920946"
        & "9764314767217019813437377032816349665076783316910523682207596376"
        & "7576971797543";
    f35 : constant string
       := "0.00000000000000000000000000000000000000000000191196320504028192"
        & "5100722376506020801011876664586186442887609794165472493729997919"
        & "6572745125453079468774667816653061833291494297813434833742305167"
        & "80941308367775";
    f36 : constant string
       := "0.00000000000000000000000000000000000000000000004902469756513543"
        & "3976941599397590276949022478579132985715066917799114679326410203"
        & "0681352439114181524840376097862899021366448571738806021378007824"
        & "8156259765045576";
    f37 : constant string
       := "0.00000000000000000000000000000000000000000000000122561743912838"
        & "5849423539984939756923725561964478324642876672944977866983160255"
        & "0767033810977854538121009402446572475534161214293470150534450195"
        & "62039064941261394";

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
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f15,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(15) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f16,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(16) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f17,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(17) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f18,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(18) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f19,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(19) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f20,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(20) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f21,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(21) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f22,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(22) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f23,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(23) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f24,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(24) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f25,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(25) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f26,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(26) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f27,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(27) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f28,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(28) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f29,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(29) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f30,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(30) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f31,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(31) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f32,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(32) read :"); Write(x);
    end if;
    read(f33,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(33) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f34,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(34) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f35,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(35) read :"); Write(x);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then return; end if;
    read(f36,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(36) read :"); Write(x);
    end if;
    read(f37,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else put_line("All parts of the i_fac(37) read :"); Write(x);
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

  procedure Test_da_eps is

    one : constant deca_double := create(1.0);
    eps : constant double_float := Deca_Double_Constants.da_eps;
    one_plus_da_eps : constant deca_double := one + eps;
    inc : constant double_float := (eps/2.0);
    one_plus_da_eps_half : constant deca_double := one + inc;

  begin
    new_line;
    put("    da_eps   :"); put(eps); new_line;
    put_line("1 + da_eps   : "); put(one_plus_da_eps,159); new_line;
    put_line("1 + da_eps/2 : "); put(one_plus_da_eps_half,159); new_line;
  end Test_da_eps;

  procedure Write_Pi is

  -- To evaluate the constants with 160 decimal places, sympy is used.
  -- >>> from sympy import pi
  -- >>> from sympy import evalf
  -- >>> pi.evalf(160)
  -- 3.141592653589793238462643383279502884197169399375105820
  -- 97494459230781640628620899862803482534211706798214808651
  -- 3282306647093844609550582231725359408128481117450
    pi : constant string
       := "3.141592653589793238462643383279502884197169399375105820"
        & "97494459230781640628620899862803482534211706798214808651"
        & "3282306647093844609550582231725359408128481117450";
  -- >>> twopi = 2*pi
  -- >>> twopi.evalf(160)
  -- 6.283185307179586476925286766559005768394338798750211641
  -- 94988918461563281257241799725606965068423413596429617302
  -- 6564613294187689219101164463450718816256962234901
    twopi : constant string
          := "6.283185307179586476925286766559005768394338798750211641"
           & "94988918461563281257241799725606965068423413596429617302"
           & "6564613294187689219101164463450718816256962234901";
  -- >>> pi2 = pi/2
  -- >>> pi2.evalf(160)
  -- 1.570796326794896619231321691639751442098584699687552910
  -- 48747229615390820314310449931401741267105853399107404325
  -- 6641153323546922304775291115862679704064240558725
    pi2 : constant string
        := "1.570796326794896619231321691639751442098584699687552910"
         & "48747229615390820314310449931401741267105853399107404325"
         & "6641153323546922304775291115862679704064240558725";
  -- >>> pi4 = pi/4
  -- >>> pi4.evalf(160)
  -- 0.785398163397448309615660845819875721049292349843776455
  -- 24373614807695410157155224965700870633552926699553702162
  -- 83205766617734611523876455579313398520321202793626
    pi4 : constant string
        := "0.785398163397448309615660845819875721049292349843776455"
         & "24373614807695410157155224965700870633552926699553702162"
         & "83205766617734611523876455579313398520321202793626";
  -- >>> threepi4 = 3*pi/4
  -- >>> threepi4.evalf(160)
  -- 2.356194490192344928846982537459627163147877049531329365
  -- 73120844423086230471465674897102611900658780098661106488
  -- 4961729985320383457162936673794019556096360838088
    threepi4 : constant string
             := "2.356194490192344928846982537459627163147877049531329365"
              & "73120844423086230471465674897102611900658780098661106488"
              & "4961729985320383457162936673794019556096360838088";
  -- >>> pi1024 = pi/1024
  -- >>> pi1024.evalf(160)
  -- 0.003067961575771282459436175178983889535348798241577251
  -- 77829584432842560195926387597522269025912316119920131649
  -- 0735627252585052582626514240460669296297000469841260
    pi1024 : constant string
           := "0.003067961575771282459436175178983889535348798241577251"
            & "77829584432842560195926387597522269025912316119920131649"
            & "0735627252585052582626514240460669296297000469841260";

    x : deca_double;
    fail : boolean;

  begin
    new_line;
    read(pi,x,fail);
    if fail
     then put_line("Reading pi from a string failed!");
     else put_line("The deca double expansion of pi :"); Write(x);
    end if;
    read(twopi,x,fail);
    if fail
     then put_line("Reading twopi from a string failed!");
     else put_line("The deca double expansion of 2*pi :"); Write(x);
    end if;
    read(pi2,x,fail);
    if fail
     then put_line("Reading pi2 from a string failed!");
     else put_line("The deca double expansion of pi/2 :"); Write(x);
    end if;
    read(pi4,x,fail);
    if fail
     then put_line("Reading pi4 from a string failed!");
     else put_line("The deca double expansion of pi/4 :"); Write(x);
    end if;
    read(threepi4,x,fail);
    if fail
     then put_line("Reading threepi4 from a string failed!");
     else put_line("The deca double expansion of 3*pi/4 :"); Write(x);
    end if;
    read(pi1024,x,fail);
    if fail
     then put_line("Reading pi1024 from a string failed!");
     else put_line("The deca double expansion of pi/1024 :"); Write(x);
    end if;
  end Write_Pi;

  procedure Write_Create ( s : in string; x : in deca_double ) is

  -- DESCRIPTION :
  --   Writes the statement of the constant definition of x,
  --   as the variable s.

  begin
    put("  "); put(s); put_line(" : constant deca_double");
    put("           := create(");
    put(thumb_right(x),2,17,3); put(",");
    put(index_right(x),2,17,3); put_line(",");
    put("                     ");
    put(middle_right(x),2,17,3); put(",");
    put(ring_right(x),2,17,3); put_line(",");
    put("                     ");
    put(pink_right(x),2,17,3); put(",");
    put(thumb_left(x),2,17,3); put_line(",");
    put("                     ");
    put(index_left(x),2,17,3); put(",");
    put(middle_left(x),2,17,3); put_line(",");
    put("                     ");
    put(ring_left(x),2,17,3); put(",");
    put(pink_left(x),2,17,3); put_line(");");
  end Write_Create;

  procedure Sine_Cosine_Table is

  -- A Python script with the evalf(160) of sympy was applied
  -- to get the constant strings as below.

    sint000 : constant string
            := "0.0030679567629659762701453654909198425189446102134519"
             & "953971468958989757975630206511538160929633594748814215"
             & "83059462621697743351395786515363562739875353550803226116";
    sint001 : constant string
            := "0.0061358846491544753596402345903725809170578863173913"
             & "293567324603009566244544704848058004531862493196429611"
             & "51889379576099507772535469616056036418297879841971709135";
    sint002 : constant string
            := "0.0092037547820598193151023784151914288481834267293798"
             & "136069270434437469160436908613450519342656077893272758"
             & "13511492090875736703323751740271437464687576766353433532";
    sint003 : constant string
            := "0.0122715382857199260794082619510032121403723195917692"
             & "500382416732467478717678330299039555616638847753459942"
             & "4828715724288165194023287218311959930019595126720510144";
    sint004 : constant string
            := "0.0153392062849881010441518676024626213017450937386568"
             & "930147905504326236057284121964514572183109886420042821"
             & "9679840418609631226321628713164460400993085612797686264";
    sint005 : constant string
            := "0.0184067299058048209273663130148401265504987492083374"
             & "102969438498348877742857037604055581045570636510261339"
             & "6751442374538643342337834081750889653747276371767553117";
    sint006 : constant string
            := "0.0214740802754695074183748977980622605209359681178804"
             & "897079228466408716203524793993911834218045420257543395"
             & "3791980095674525649626648191208084577969516787728236814";
    sint007 : constant string
            := "0.0245412285229122880317345294592829250654661192394514"
             & "775767567737884929481373947340741505810955139589365788"
             & "4873323319293879341156639761420670724898789962377250913";
    sint008 : constant string
            := "0.0276081457789657416123548717439784658641636872579351"
             & "767458936847395238360869290068361719678590186040785858"
             & "7171093092179215921425777860104174757883287857318732066";
    sint009 : constant string
            := "0.0306748031766366259340210275652237129977158501934989"
             & "490493518369635002304189534882726798972483577776422987"
             & "5527026519908152294758226765895144375487544462701996629";
    sint010 : constant string
            := "0.0337411718513775848337161124064239499621908244133321"
             & "147176627523000057304347553116142168312633657961743319"
             & "4418828651835682132661349363119754987221306578405922652";
    sint011 : constant string
            := "0.0368072229413588323243326909279513010552317719584573"
             & "991270406521191349411640973866059376600351670385798327"
             & "8415529700479250605049586002432611067418082700464844422";
    sint012 : constant string
            := "0.0398729275877398111285787376798802131956085404986545"
             & "607168397113812318940862671905262026720712197533124666"
             & "5840877552931935442635344861730854930936718048393901559";
    sint013 : constant string
            := "0.0429382569349408230771245402817839553939332746003164"
             & "055238510235477711595141985694950520660025025156667624"
             & "2980453548446766572505880373673637786611712568927220555";
    sint014 : constant string
            := "0.0460031821309146288143017879102414734403645628801823"
             & "322697973836392193145630803181413947004672597150693000"
             & "7342447406965629362757539908653614075295334966440279421";
    sint015 : constant string
            := "0.0490676743274180142549549769426826583147453630257529"
             & "202101225532691659395639607568222207671118007237540668"
             & "6097746558346131875748097793857856904513601493805535142";
    sint016 : constant string
            := "0.0521317046802833212363582164233367485583059324482465"
             & "249282482988320236995506600617153979447938540006359257"
             & "7962881628467177910520957287381358117604735916964148795";
    sint017 : constant string
            := "0.0551952443496899398094475256977309118980986639932421"
             & "058784004356548895024883146963677035872057828874079553"
             & "2011618958409260073379466333277814206900968809916805032";
    sint018 : constant string
            := "0.0582582645004357596139797819343505758564073174794686"
             & "451223567110178390567524446152906342167667826943804277"
             & "3691145482287203007632639578531915003109649077208426897";
    sint019 : constant string
            := "0.0613207363022085777826145929172350071939802796769731"
             & "652631930233073550398959628364542397455942871312109220"
             & "0445476683861131982124358568737585832404545847710164956";
    sint020 : constant string
            := "0.0643826309298574608193245368265668401635526079771437"
             & "676241010823076329196026827690263403204294406233089670"
             & "5353758036931119596147125911081500052247388902067988677";
    sint021 : constant string
            := "0.0674439195636640578979724218744933647690956421490120"
             & "689450947962656047478096050136701068308041516702435356"
             & "3420436813185021772927868604600653450707572096109499367";
    sint022 : constant string
            := "0.0705045733896138630273514705529420755992491548764074"
             & "945726869361035649510365154767643657508608334975274536"
             & "3957766740141632122219177860022977706003857646883012937";
    sint023 : constant string
            := "0.0735645635996674235294656215752343218132992656778872"
             & "442091692001525149465249903079306352217560656216324609"
             & "0544308763448679450301879516441627447204482134771117772";
    sint024 : constant string
            := "0.0766238613920314922783324628237898766848511890571246"
             & "492549418814500150275364564006640774834474325215384700"
             & "2703065971448445738948378704926946280492350776172041544";
    sint025 : constant string
            := "0.0796824379714301211471206559958854013639452967811396"
             & "339998965642878222905090589215211582134818509289718966"
             & "1140924858209229398320318386263863555840117042369125680";
    sint026 : constant string
            := "0.0827402645493756931119870831858619740379930597113823"
             & "915502305770040283986563622765088795014629951914770441"
             & "6082594013513036876942819169914938027719954538576340577";
    sint027 : constant string
            := "0.0857973123444398904615563321468461719405315153026428"
             & "988435681912372625497215351651236721322280041695357898"
             & "8548530304480226083586774234618773753794888040187493964";
    sint028 : constant string
            := "0.0888535525825245965615865350033739514293329062219376"
             & "145362942707169527566098113804037732924765816896850521"
             & "2707002325458784148551462945509817781217979910241254278";
    sint029 : constant string
            := "0.0919089564971327286249909790776994870092419232930348"
             & "413509974936530677796700136083525674292855313317630259"
             & "5936745699805288897745406088034883219049260270805659297";
    sint030 : constant string
            := "0.0949634953296389989380343123604925943134232383578768"
             & "199767976229285364047420420757855103861804731497390995"
             & "1391455949768672410229728628752980152275348724165865118";
    sint031 : constant string
            := "0.0980171403295606019941955638886418458611366731675005"
             & "672572649798093873027890875368071110771463185595540742"
             & "0652644410215687524684302442393412859968714808558841196";
    sint032 : constant string
            := "0.1010698627548278249878875845507100762213559985945653"
             & "472145443054814579275842923554740620153545693653091531"
             & "829221536761599796765098187057452103546021610404272109";
    sint033 : constant string
            := "0.1041216338720545791209438800601790207959749274903091"
             & "947185613660056851894964907800660776969635807967793521"
             & "043619393205843427394709461672155472470226970844019631";
    sint034 : constant string
            := "0.1071724249568088491755291482281967503779193796583499"
             & "779700652020405492612745401432526400745163454984421670"
             & "320580743551347148463376499902506803525031691646537783";
    sint035 : constant string
            := "0.1102222072938830588078991402156777252744746231149868"
             & "762164343487878463726789387838462304925486091225683518"
             & "891808825321976681351136528770483671615701102679472633";
    sint036 : constant string
            := "0.1132709521775643490182287329082850660797838690096898"
             & "863359229809609257813975054426989624128029689741772784"
             & "126723856065215742063210089346277329992855128753860970";
    sint037 : constant string
            := "0.1163186309119047672525443194705125409238566764922737"
             & "103284823008114364699561563219959396264027906063467476"
             & "428403041620588019002394443189820501050181923500380510";
    sint038 : constant string
            := "0.1193652148109913645936377898047947467110184824121012"
             & "990740624445550955097871079956241295827888490384056980"
             & "312542388026903365873077616922457649905955944757878296";
    sint039 : constant string
            := "0.1224106751992161984987044741509457875752236090851072"
             & "485773187024063510194872256293642941583891435347409939"
             & "546081853122649616772937850341507948866397664753406951";
    sint040 : constant string
            := "0.1254549834115462385423364532675945494325282092800867"
             & "244102468434687594777940698904771232954950995269288980"
             & "458298673555649579429096528808638400295735879299821245";
    sint041 : constant string
            := "0.1284981107937931726244155891727575476975072849139067"
             & "467049686590181812159401129116307688479310509669045348"
             & "150715432392326481143583722247380779798238448201362483";
    sint042 : constant string
            := "0.1315400287028831111033874926922301510074010766035916"
             & "890220473051757135560484005290470858672195362241148884"
             & "302214104376328787765629454684894489745662233772976326";
    sint043 : constant string
            := "0.1345807085071261863163584092539792556314246096030667"
             & "794184456873039562949256817967758427591489535350168066"
             & "922671824341747333827186486653293605070056698959002399";
    sint044 : constant string
            := "0.1376201215864860449484416634310973343662606497934429"
             & "611714989898675814618394878177612574414178926679401689"
             & "650428985376505635603267135356629586877378461315574959";
    sint045 : constant string
            := "0.1406582393328492307147888464071430911197878898999696"
             & "127538723881418578855043578625789877994288147213066043"
             & "478974379031297686872655070320394462509779357058667328";
    sint046 : constant string
            := "0.1436950331502944548197733493230505826203393800429713"
             & "977001593704205048177084676510985739254973016784891207"
             & "091301743130332593932601433318232007503781481102542885";
    sint047 : constant string
            := "0.1467304744553617516588501296467178197062153165293920"
             & "668973938768923171676867758277059618817711050948125690"
             & "558476379337616896667576365627798992375193448833326562";
    sint048 : constant string
            := "0.1497645346773215172296957373427385287768305173155205"
             & "944444038558652026598191276150124939773224839797345312"
             & "570326512200707512553366125599860711241909187453516799";
    sint049 : constant string
            := "0.1527971852584434277203366125438313132217351043555721"
             & "869605997356764876132128246479495405330992447897472145"
             & "326702096385561918941783925779088822673452488739252197";
    sint050 : constant string
            := "0.1558283976542652357431014862462223493747859198377059"
             & "500956243530417661606847322898498362634696292005485400"
             & "660968208270900407390008362078124824964360543303344109";
    sint051 : constant string
            := "0.1588581433338614416843853596530813017855264604464552"
             & "249536981739616921261330337502389330825502347697979805"
             & "001620277873683828850510474928183192216361237075115526";
    sint052 : constant string
            := "0.1618863937801118376413879953333824608490014042822816"
             & "066257217443298522299444388617851003822303615474930434"
             & "841464796142832371074545838486274371166911628026771518";
    sint053 : constant string
            := "0.1649131204899699214181891132844101245584977531171477"
             & "283189379971883900593364668195125454277509079100427221"
             & "712548185852670384363234664410584761879606107650195324";
    sint054 : constant string
            := "0.1679382949747311780547455359966578549999671775156040"
             & "891137248480506205841504574665679184078199332049201103"
             & "520028610731185279563706752889587328963301759504319001";
    sint055 : constant string
            := "0.1709618887603012263636423572082635319663290591445910"
             & "012596894557734102929076287581584670239143184040811875"
             & "435497336675697952852373347815029074446204950059668691";
    sint056 : constant string
            := "0.1739838733874638279507008074667317502779946746329410"
             & "658045649990593665397244083898604295139382845848944327"
             & "780482637304591538337412952036844765709172661513075638";
    sint057 : constant string
            := "0.1770042204121487561968398439291657307197507865185694"
             & "054959712981727525072009635977053491894129344306025520"
             & "576628876718194873732523220888932215432740164096241686";
    sint058 : constant string
            := "0.1800229014056995226799065898456067867505057226925682"
             & "852574957995762482466538285297798910671393036846386862"
             & "703291663011632115912097538903578505788092772248273898";
    sint059 : constant string
            := "0.1830398879551409585165325784769200138775632254052333"
             & "280565068804504059490554431319073868536585392761499461"
             & "027105354717875439522251555395721962205454182320526933";
    sint060 : constant string
            := "0.1860551516634466481054383041691617772795158180940969"
             & "962866574421467013568062246697956552897317269569331482"
             & "367632973662863946295951616450809587911975473674657387";
    sint061 : constant string
            := "0.1890686641498062127549978370874602018495823335434643"
             & "580564464916344443836388470574133991029539869581228698"
             & "172865196875361556939155160452506329351211238237667910";
    sint062 : constant string
            := "0.1920803970498924416792882046279485355947208168729749"
             & "228413443541811512053491245791265876135572673897645255"
             & "843973394487309229842331353089630746875740030204880000";
    sint063 : constant string
            := "0.1950903220161282678482848684770222409276916177519548"
             & "077545020894947633187859245802253253092340903817309920"
             & "701055366117589657998348476073888265914344275520716095";
    sint064 : constant string
            := "0.1980984107179535861793249181510733849221204375709235"
             & "023812777319344216088053886274663426277947312680136000"
             & "982196113444302643852171352778723529102710797534972445";
    sint065 : constant string
            := "0.2011046348420919115584435458820671774655273168648817"
             & "352323864815781581188420293117431478902484563306608924"
             & "162163486463145323060741320960191768733194022693343927";
    sint066 : constant string
            := "0.2041089660928168741816969499474189576619453422964934"
             & "729950150565701709863573527480651505149501379071287229"
             & "939319257473383169911550172746049339150215578256240664";
    sint067 : constant string
            := "0.2071113761922185497081160197866674167674139345372761"
             & "775960795601783456612437330322511583414223902359857791"
             & "028401115584741293058922090368444614494912337377916361";
    sint068 : constant string
            := "0.2101118368804696217174899720901250451970695512599436"
             & "652108556451898259504376085813325649770923123670942402"
             & "296479455089230639808555125946624689424916878332924930";
    sint069 : constant string
            := "0.2131103199160913739677575178515008423979589821936905"
             & "244000263807712401308857601312582098440534995357671946"
             & "476457314255418587673490619590764316750261943951436587";
    sint070 : constant string
            := "0.2161067970762195099483851312908295838456946552854169"
             & "793529831931866315990199882245978569335662345437291519"
             & "313656282395853273777053172277024153243943778586547984";
    sint071 : constant string
            := "0.2191012401568697972277375474973577988483607967055921"
             & "085026246329098112435823550704068816745768548328961855"
             & "395039949531274554839545365041418372412606732403429514";
    sint072 : constant string
            := "0.2220936209732035340940947213139774485664983603248351"
             & "303370813517104400706786900738912195230628208055526427"
             & "635552612383947562023820420103411198034583547578192139";
    sint073 : constant string
            := "0.2250839113597928359916421198633534633451336324953449"
             & "929123582830508862734414983307301963032022948069632054"
             & "243872926583990779258320011255789689040176801332574006";
    sint074 : constant string
            := "0.2280720831708857392544573794575372441626994680061817"
             & "717327915962597737135789155123583246830012698831926187"
             & "690978109363291607343971403156220065609660520016326734";
    sint075 : constant string
            := "0.2310581082806711196432360184727066652304106304332499"
             & "034754787627016066610915480229880729291823228214413711"
             & "947833422784072545470803196652838243969480501494933870";
    sint076 : constant string
            := "0.2340419585835434231912420449226212800691004582593853"
             & "317524045774735367698614878370152011238056976716759899"
             & "296341588527702813485902374565522109962745009689941965";
    sint077 : constant string
            := "0.2370236059943672068677359145212644641039880051746337"
             & "859577920504773580535296815838446941170393114248363487"
             & "287172272781497999602026265786719359455676804526398600";
    sint078 : constant string
            := "0.2400030224487414865689223653588895716013486111028527"
             & "542367150976672180677630304581666346936817822285389694"
             & "055052717095204075958405714163449471019194094142221257";
    sint079 : constant string
            := "0.2429801799032638899482741620774711183209907832838321"
             & "260423208417946981663958870487838239414776082706567113"
             & "809212419797439857436275167386038607246447923207968505";
    sint080 : constant string
            := "0.2459550503357946115999247085519682118632786652695591"
             & "815243556320971753901406174334260709052943087037003758"
             & "868440900660697435006974522990692189208064907565233739";
    sint081 : constant string
            := "0.2489276057457201681106828162729877050769317717225214"
             & "148620243492995773755201416183215765857822061422872445"
             & "175780301496531831751370206372035276482653783339658757";
    sint082 : constant string
            := "0.2518978181542169504981066283742714258673799187698320"
             & "950278639112246802345119190755183708115090393417316191"
             & "419001853581146727945127780510161810642567920964262384";
    sint083 : constant string
            := "0.2548656596045145715539807788247035490779248168185416"
             & "691847029015242716092767782953110555934130876439041182"
             & "937402705182787951330575453135979999689649679507543582";
    sint084 : constant string
            := "0.2578311021621590056144712947590814441497126759390575"
             & "943305003334694421679329753588356942119043161244808915"
             & "630639807111173031622253543230624493013092610465256488";
    sint085 : constant string
            := "0.2607941179152755182801865090847883426487919206749109"
             & "574696319510928797532549056279570119960704052461438751"
             & "257556024085926528646163364782215669166121996931393470";
    sint086 : constant string
            := "0.2637546789748313836113493219832235967119934070276339"
             & "532296885820048967988787995909934505185063017054696760"
             & "207946376998182446890249667969547375557195878896455639";
    sint087 : constant string
            := "0.2667127574748983863252865151164363940421169883561562"
             & "081989024345080212831231418753245630552591635424722296"
             & "970341971436103120233338649903679175939660494374539714";
    sint088 : constant string
            := "0.2696683255729151065254644624226725396346525991952706"
             & "398329285603879918784591839717586333477069003868284519"
             & "288108507036379790919840414837440409844927657957219584";
    sint089 : constant string
            := "0.2726213554499489844933474772922102408041755580770619"
             & "026594726427921093737282625702489119160101349605901676"
             & "341992112685083878024714705085895399065584807044705397";
    sint090 : constant string
            := "0.2755718193109581630764251683907659018326603860153325"
             & "642733575371349810524607630141097230611496034073543127"
             & "405812274472377540226685218718471742735953664373445332";
    sint091 : constant string
            := "0.2785196893850531052078485259570194793626584607718518"
             & "665295610233074975007968133888519732729596325275491261"
             & "150862709086548130437666236208319032222344214285437035";
    sint092 : constant string
            := "0.2814649379257579840952310073400376703951576887475711"
             & "932180251708887761310156844315512315350430052172224214"
             & "066170092027372439983515516284687156891265323447200594";
    sint093 : constant string
            := "0.2844075372112718436183106149398162274615807409638922"
             & "373945413262728690190259537764158871113645460130107533"
             & "083835695899356704063070873893120753615774534285314141";
    sint094 : constant string
            := "0.2873474595447295264773318414299190824647101685673589"
             & "622734061895346010495925554047417553950305864075616078"
             & "182247317656111540716627721413082305462594448566429522";
    sint095 : constant string
            := "0.2902846772544623676361923758173952746914762783241511"
             & "114206671131253928982994574468532436606893438468421981"
             & "569266985974726632161936158765975844925431147293592719";
    sint096 : constant string
            := "0.2932191626942586506066085989561645990749759199300713"
             & "653651177739816498514982610098314366431463802979988291"
             & "065963619930004728395109048801084193351040707420102401";
    sint097 : constant string
            := "0.2961508882436238241217861277826588797865038519096636"
             & "779958151129536903429743522448388946435785091030929190"
             & "235355439125946217441183235181498849881735706117656099";
    sint098 : constant string
            := "0.2990798263080404767503369727760827896234515469302238"
             & "094621249196809632193386632939258427223608021178241407"
             & "958947660643493862579352952416641628024837104604457714";
    sint099 : constant string
            := "0.3020059493192280670034632317324243912848169556015980"
             & "256384425145434035973129048750553850354021934377521859"
             & "926303677858550070242817916888960925671089608320816501";
    sint100 : constant string
            := "0.3049292297354024064907286334365223463192519271026275"
             & "918567468129123724603911315960759021201024090157123366"
             & "254791230512527199479340994613222001782536372830403246";
    sint101 : constant string
            := "0.3078496400415348936820636455585202015094184912737627"
             & "553150503322770208324630417182556808257602291804930959"
             & "198088831064734309553913282524033525352458482240890320";
    sint102 : constant string
            := "0.3107671527496114958359972502117632942196132072684204"
             & "428620846281659509548355038269353084251540200658468351"
             & "705514520539825219448442043732174700328541197150642184";
    sint103 : constant string
            := "0.3136817403988914766564788459941003099933775094565467"
             & "851932847388033398137344679026583246149717961566235954"
             & "198732386085078444263815068892657515312999186282843991";
    sint104 : constant string
            := "0.3165933755561658672430470346829517504096969166492718"
             & "883162986204146192975248430663290601076797398273531595"
             & "785673673232792825909766963497280251525934332039184841";
    sint105 : constant string
            := "0.3195020308160156779015182715397565772982056106134798"
             & "799855121117012355370919647825868253578023825912881893"
             & "441156290508573647875448029678151471909056012464306250";
    sint106 : constant string
            := "0.3224076788010698483848074776591158480870902495070670"
             & "106805740766909232272165044162817919786326175238291365"
             & "805299618719579486595649462130206202821498970180897828";
    sint107 : constant string
            := "0.3253102921622629341359547080141967703657466600660294"
             & "553285869535761121159776102343671524645071413736559779"
             & "131711234838468962730339977932466888263798560197803753";
    sint108 : constant string
            := "0.3282098435790925261079168166295056612837898419418725"
             & "944234413214229647379607651698191151528544065564816165"
             & "380066829323599834582312058457523599322425628129700122";
    sint109 : constant string
            := "0.3311063057598764017371907372666501496458861183171805"
             & "551051946784530741512178795861156507461644458143553377"
             & "268515519956600112161840882832345317335804895979773647";
    sint110 : constant string
            := "0.3339996514420094046508654805349790667399323599263156"
             & "946221339881600014597038541584093163038424464117165422"
             & "345402942207689664860133289652187762264869057936681671";
    sint111 : constant string
            := "0.3368898533922200506892532126191475704777667796712222"
             & "840515309421523592933598104109594588644731514875566646"
             & "455342412120020059137247381125522626800604920655000914";
    sint112 : constant string
            := "0.3397768844068268578288258028174093267566821180251100"
             & "390776254190116119713334976331110001903719930996126744"
             & "385025482773285054720339441091322021015242110437372129";
    sint113 : constant string
            := "0.3426607173119943975927819825612868796115591641438863"
             & "646620820902565103339096773237999505132004361568708638"
             & "964487908550150105906557533294004777672240633927104908";
    sint114 : constant string
            := "0.3455413249639890655391917230787393996813489848402525"
             & "087688677092509107784472764210007539829478194774765170"
             & "063260130437397825102845735880846274702067792661127485";
    sint115 : constant string
            := "0.3484186802494345684193085876948116364531209085280093"
             & "224679374582799681831732945144230449272152019415834710"
             & "443118508175920770082552196458983517741342187953065798";
    sint116 : constant string
            := "0.3512927560855671256013076230482732902746930770598977"
             & "253119215370091420174209160787857292082854299676870584"
             & "126632072870259970950803752173976736304610116919933890";
    sint117 : constant string
            := "0.3541635254204903823573957961389289436032296670206147"
             & "338379015563202835449557421950661038042795152235634498"
             & "249296458801321798764353513761013463355676224755876316";
    sint118 : constant string
            := "0.3570309612334300326149540357940367878148326224386568"
             & "290621142200414124655942644050783002348336093305425407"
             & "695875424224171316443481546952504823107767903467141355";
    sint119 : constant string
            := "0.3598950365349881487751045723267564202023174211290258"
             & "497630100778756401389522703121629916592490241555351151"
             & "531079351839697550218630851277756334889268487566649029";
    sint120 : constant string
            := "0.3627557243673972162048544621153241297874011048202520"
             & "155940190663014948554028365186983075440044111335827799"
             & "881434392552984491206801634131858371646854795524844938";
    sint121 : constant string
            := "0.3656129978047738700117459086069781700085140664294972"
             & "528903290581509153243479155962255861874122016102017092"
             & "598724863162982430436751006030357307971349089084457711";
    sint122 : constant string
            := "0.3684668299533723317127462216816982941984166163223511"
             & "021535623999869634686700797144574672689157817468955998"
             & "006689115009732886045956264516413146410978908226324906";
    sint123 : constant string
            := "0.3713171939518375434119349670219232661775725301769655"
             & "376587020886674745864278377441012715035639620971063021"
             & "496578278136708216168846279182766718786129485213981366";
    sint124 : constant string
            := "0.3741640629714579971043930195383235683167971603106387"
             & "167716978305976400246386965264533249658849709630537819"
             & "958693450669225015128400126843381152690200790484195367";
    sint125 : constant string
            := "0.3770074102164182567265678231998572323015378836479092"
             & "748594578595133254809758399906775242029870683486807553"
             & "121426412740061805531281941009270258496256196834755035";
    sint126 : constant string
            := "0.3798472089240511705762811467990666753860053150308724"
             & "931364148298657109440662952902056236223596557586600602"
             & "724894311144766027744437971613620457127339644146722897";
    sint127 : constant string
            := "0.3826834323650897717284599840303988667613445624856270"
             & "414338006356275460339600896922370137853422835471484242"
             & "886614935559007560102009675979208442091777288702111640";
    sint128 : constant string
            := "0.3855160538439188640756079493391946817280674679258861"
             & "679626726611129511902898642133121046856319143705318248"
             & "486199997255852469189307950333625124015310223484442734";
    sint129 : constant string
            := "0.3883450466988262916249935406705281014962833494394689"
             & "315956273329869069584940237001847685357330624351347671"
             & "191899549031927966660847610501348041435309082406644024";
    sint130 : constant string
            := "0.3911703843022538886875129486588618994434759784047975"
             & "256572286445776044968494016656469734591967348379205194"
             & "332085953878914898607410940176466723796537697174430766";
    sint131 : constant string
            := "0.3939920400610481085961886608903134244856522968110642"
             & "025862597880650858688180689520062438121892504346348068"
             & "717173583751488773992105397816400287258478083786057334";
    sint132 : constant string
            := "0.3968099874167103285952909113684558637493841960279746"
             & "464374813170959603329888314450244585686125729661948025"
             & "477930883642129567432813096298359270732102739940165520";
    sint133 : constant string
            := "0.3996241998456468285441170307420208607299235371692063"
             & "270556043397293889007199834800706515987380019972789747"
             & "779575879154006348694283918979140116946187140502906882";
    sint134 : constant string
            := "0.4024346508594184410825339335196296727887695860787983"
             & "330348426745285289846508941450120316989952994853837895"
             & "327016564286535127514480518951059237440075902245842106";
    sint135 : constant string
            := "0.4052413140049898709084813055050524665119477541050890"
             & "138058289893394311023351135448467297414812047645030267"
             & "707081423608221032123397377611618879087991399983894792";
    sint136 : constant string
            := "0.4080441628649786808207474989303193628209072906898636"
             & "067407340343582383381235642012728683126785912433359263"
             & "776235704395158977205530209034970854882544978857114091";
    sint137 : constant string
            := "0.4108431710579039421834666749289471041367041399440958"
             & "550780249699787653543310983875337898182315867885892238"
             & "589329360610218067491164625110717907138488936579325326";
    sint138 : constant string
            := "0.4136383122384345474719443235553664460012225907222602"
             & "846707167034147035999858249255474217540530440558651808"
             & "123523654794896620643981070336838076496120544035731250";
    sint139 : constant string
            := "0.4164295600976371825625989107894802468556074172857931"
             & "058987646494032866402575286675034298735441291868556500"
             & "209651269189808501225573440790067038653625822746904457";
    sint140 : constant string
            := "0.4192168883632239564330100199299511912473338408393713"
             & "432982946760067996254990238688798935469245422415359713"
             & "724822873518408583942948711200645465367464904064252948";
    sint141 : constant string
            := "0.4220002707997996859412879413320511405711334230154070"
             & "671722377080329788737291963477463296638801290977153227"
             & "400580753386771554776748238186498548336771731859797165";
    sint142 : constant string
            := "0.4247796812091088333572261892346645556946332007941068"
             & "757557751999991434480866013883200286352102348576607470"
             & "781492876080640501489970740378210211680063767872820249";
    sint143 : constant string
            := "0.4275550934302820943209668568887985343045786293424586"
             & "393648472046518438249098519379985739495457378037624943"
             & "453661015350997872508849588593069594248458435132139985";
    sint144 : constant string
            := "0.4303264813400826339081990305829805911285178100623346"
             & "043952892612779095579648710350617740521672820193384191"
             & "217738666278541222034291852637642052253186963310268374";
    sint145 : constant string
            := "0.4330938188531519684842226384895773616555246753722248"
             & "816024451496419591000693454963887390401225672107089892"
             & "364853991300388950510119087763080461828777372742320119";
    sint146 : constant string
            := "0.4358570799222554910325440803550755497283539552018526"
             & "669107410880240005586214739967150233125821154632754489"
             & "099702500383748546190017179556733962714043986734532534";
    sint147 : constant string
            := "0.4386162385385276376470257375461343598431014545442370"
             & "392629572239749453814826748644613989040384611131123410"
             & "196502588788254702626969625774843348408799712335325982";
    sint148 : constant string
            := "0.4413712687317166928799889682561164890337329042919533"
             & "093469963042270510760356721922852655307727922950686614"
             & "048011744896885023391417213295622280684976455446675822";
    sint149 : constant string
            := "0.4441221445704292316420694179834663691061023974011178"
             & "866841277090821315679521870047100800232390493591703751"
             & "676799014459898569203824455253477294847969891989621442";
    sint150 : constant string
            := "0.4468688401623741953530443887189265744266000023182254"
             & "650989703709025050177480703385323936744124241248871376"
             & "990986306008179682665447256732513069256941034845626974";
    sint151 : constant string
            := "0.4496113296546066000462945794242270758831870483778572"
             & "623655778812241904757362305784403432955948240634275609"
             & "799906262292838825302304420748030347222668022335182491";
    sint152 : constant string
            := "0.4523495872337708741330267029477754348290775509171694"
             & "613782517666006691883227550741759973290738878598811441"
             & "759203948643191538521683515034197150465402751083870537";
    sint153 : constant string
            := "0.4550835871263438235358692678967193602180057693617736"
             & "026727590793653638793805343698251691131176235850745160"
             & "203728291917576025159159642845407264129570349101640041";
    sint154 : constant string
            := "0.4578133035988772219049611553600256307008096136593152"
             & "933162216188529724108517564474185625293410506884199943"
             & "104202337441852908170850654256723535845547910965158853";
    sint155 : constant string
            := "0.4605387109582400236331814867414884243654576306456521"
             & "936829498960857226391413424386167986485136657837965191"
             & "702727966372051804085781579581756604719443496506876041";
    sint156 : constant string
            := "0.4632597835518601973907196370998114210798243014169286"
             & "606058491378684376946398398267872155880210157891011550"
             & "275174764651373480655709321564286006858887245408322753";
    sint157 : constant string
            := "0.4659764957679661779027560648877787034268950367328910"
             & "379702059564331413700435780236467855527735494041939504"
             & "615996446027026045789171115650755344424379050581648000";
    sint158 : constant string
            := "0.4686888220358279336976178702147423143397742323839151"
             & "459274730391300826704635523173813987940312770929397456"
             & "597377974871993462959538553651352664889663714755280343";
    sint159 : constant string
            := "0.4713967368259976485563876259052543776574603189324806"
             & "214016140310088352216651617529075631968879121177453582"
             & "548889551861368392179583452027771212751033145034183893";
    sint160 : constant string
            := "0.4741002146505500143985800146693645451235605768098899"
             & "666996257525821400268772851507418240100459495701337335"
             & "600033881275919380282470901269082362388249957030095149";
    sint161 : constant string
            := "0.4767992300633221333421581174135755721060696321798402"
             & "275257402669866139120279932334723702748725821837837707"
             & "479032569544807786828871882053164199959542031437157737";
    sint162 : constant string
            := "0.4794937576601530266798397976816903187312904532233566"
             & "297927930487177784233100846998650310625757221028380438"
             & "953111194389195607738922082295955649841715559325906886";
    sint163 : constant string
            := "0.4821837720791227485173444807974086112395892088444544"
             & "430747966272186590561952851698670976509806014910801589"
             & "907839273718328208480332619258661654673289417017315951";
    sint164 : constant string
            := "0.4848692480007911018229516986611111799490454472883253"
             & "494358003366183788613769070560940843463262055170566439"
             & "430812178022139872842741919720303066134013664158024974";
    sint165 : constant string
            := "0.4875501601484359546414850273076499075501098248882845"
             & "855719163327501235851994996111443061515185390876195029"
             & "731963854170254297528986071127060846021974998782760713";
    sint166 : constant string
            := "0.4902264832882911542295984490366087200183623387998826"
             & "273274744020450973640760878845770222081060328612474265"
             & "634926207172254336911739581012679931201281449304598154";
    sint167 : constant string
            := "0.4928981922297840368730266887588092682396873065483635"
             & "811114318744863825867427922036946116984077256090878411"
             & "750930004452956458840869839113436170005989757227050318";
    sint168 : constant string
            := "0.4955652618257725311502666695414869190434385500595884"
             & "044377154166428885268609119064802692431388651958958193"
             & "092953383819754425364895687921601249585092512194940568";
    sint169 : constant string
            := "0.4982276669727818524109838693598428202796319260431717"
             & "415486557398915976030519668069428635918823632385358216"
             & "565213621356169709805662784362369345688036144825622401";
    sint170 : constant string
            := "0.5008853826112407862412850037568218781730189009464244"
             & "529972233464239431259981189002467037347644511193338840"
             & "379929673124018057643141833273147882471650258478758053";
    sint171 : constant string
            := "0.5035383837257175586918670712604959834510947782836330"
             & "000052646722559809208656301981901541887020363397998796"
             & "998482700707057613025220886495807215450775645703418785";
    sint172 : constant string
            := "0.5061866453451552910489423435964958612549020365952933"
             & "947357664755626033671560566358284122682689221931802143"
             & "845120211865269589137355867202539744646628245946889972";
    sint173 : constant string
            := "0.5088301425431070369317493243516528643705928238912837"
             & "603929570591425500949957261576747316699101598933636570"
             & "325570106318731814468263146526485898860684194788133223";
    sint174 : constant string
            := "0.5114688504379703995043910009878180819645876215508495"
             & "736727383302681374994910764658317456346090814265081485"
             & "191413425074293606622319930602026446273847947500874296";
    sint175 : constant string
            := "0.5141027441932217265936938389688157726080491204162178"
             & "042879364539860379272742201183345765668289017940589653"
             & "702775218553661014349955692424614950660832991151061171";
    sint176 : constant string
            := "0.5167317990176498815087538760497663859531390303097521"
             & "104882035920141093319497283077194072311739602127232125"
             & "496375998029131489227838634555602176734112170566302994";
    sint177 : constant string
            := "0.5193559901655895873618299320920460363301387807596576"
             & "524517362128678295448363841925459578449457646464884168"
             & "655477615655466299466985684836967422975133929250163254";
    sint178 : constant string
            := "0.5219752929371543426942583175191106190747666027324104"
             & "558844314629004700707423659574541922168676975124514908"
             & "693267691421782355471703421827198712895225292940516145";
    sint179 : constant string
            := "0.5245896826784689062150984639335456725563206448119373"
             & "588726363875959195963160091506546822732900757611815021"
             & "897390506302939880801036215065561381218662489214706402";
    sint180 : constant string
            := "0.5271991347819013484642745754946705568645495625871803"
             & "817371024059366745883698908225064781859809257946743177"
             & "657178578892368195005972026898866614274425069706096914";
    sint181 : constant string
            := "0.5298036246862946682160546712352691063333056111012520"
             & "340517860983261173999865508431625968047726255931363212"
             & "118239425588202029961273572056067525659638490642360363";
    sint182 : constant string
            := "0.5324031278771979714428052182081180752112809608216542"
             & "138573406891345035927216559117170727126154960716978085"
             & "956499154509764892942151171976222857819121186196694174";
    sint183 : constant string
            := "0.5349976198870972106630769046370179155602656921719002"
             & "684556738510206831412598337502398630470672918217663346"
             & "744546653075097480586334339642430449559792901419082865";
    sint184 : constant string
            := "0.5375870762956454825022149323489381597209519052688873"
             & "708295028039411929637002733699816284603309066279852551"
             & "675360607949876364572790988846475571376965425482212859";
    sint185 : constant string
            := "0.5401714727298928812978454797368391243804773726575731"
             & "328559542076742614251127228693834535014490957262471044"
             & "685523369104702249940763652938638851656514379030760394";
    sint186 : constant string
            := "0.5427507848645159065867686612074833318665104452394303"
             & "870138832970441219322200386415236444326872372056700616"
             & "651039546977208778873417161825165819232663294367409789";
    sint187 : constant string
            := "0.5453249884220464223139873471738261971863535360472649"
             & "722480657047027389676435677651676327386251232354422477"
             & "892582883818991882337412619549802139889527022938112208";
    sint188 : constant string
            := "0.5478940591731001656088205706344303158431835263620271"
             & "505432991519430723608866511100780725112451994454175155"
             & "309410869742519474891911459633158184249777712159135720";
    sint189 : constant string
            := "0.5504579729366048029772898925285391908407957993772653"
             & "460290228409481027739716884170783773450786742276397881"
             & "710743384297383834289429705741080609359850486732460113";
    sint190 : constant string
            := "0.5530167055800275317642269884598102702774565804435897"
             & "678269512669499655247677167205902773101266730008728096"
             & "714408388038337237059491298958630518071912351585840828";
    sint191 : constant string
            := "0.5555702330196022247428308139485328743749371907548040"
             & "459241535282029492475774800683831248803512639343555379"
             & "458249282081463006489743886788062649174003545193643593";
    sint192 : constant string
            := "0.5581185312205561156937029638155638457395303538738867"
             & "342363454889817711309235287972510008872208141190248685"
             & "098635056244014106178472174392967235839413989520450995";
    sint193 : constant string
            := "0.5606615761973360238397102231455323464061827104244866"
             & "484907053386673352487699305395072211804866093502257934"
             & "635032183827706624038579211404216282607364765265567096";
    sint194 : constant string
            := "0.5631993440138341150073637718570158237455984152604557"
             & "338053465946957573353776284664373062059806169926968300"
             & "909939060488541433361950079519774928295983456150116970";
    sint195 : constant string
            := "0.5657318107836131973897650113692660904743014432468325"
             & "846677441006973264595406666136283809542103063533494605"
             & "064479803995331102355444455393368970326780370817631851";
    sint196 : constant string
            := "0.5682589526701315497905484891559202039443500356891325"
             & "755674634753283773618214707633616938651378335280668918"
             & "869772038592305594758396485776304332307982795874142096";
    sint197 : constant string
            := "0.5707807458869672802326528638849623650065504461592043"
             & "109527031644842099082840634183914531040168660657726878"
             & "317569441340161492007543507695867138504065662072735513";
    sint198 : constant string
            := "0.5732971666980422128201712389421399336043776624972067"
             & "178580572285003471764163548099061656682637589078026644"
             & "957261023230214086707748360486449313672393775610266725";
    sint199 : constant string
            := "0.5758081914178453007459724538157308417760084553140966"
             & "022029889262442367885227670337478439943268291258751515"
             & "043273020343810366353419781830506440676659671518980652";
    sint200 : constant string
            := "0.5783137964116555633422450192905771729677952159265018"
             & "988362914963433729100905933673967337558978004133322020"
             & "887956465265115479203698979690676709262350418450910828";
    sint201 : constant string
            := "0.5808139580957645450755952716785138763990615817533023"
             & "272667800751035973385574586595142774992220216564099392"
             & "657515904338935539591361775870307427383175907605427131";
    sint202 : constant string
            := "0.5833086529376982943928309612343084453657295517112574"
             & "064951557779390422051524371047513406771389330365941422"
             & "055182955930529572558471403432675094486967404858145505";
    sint203 : constant string
            := "0.5857978574564388603280808381186622857671368129040376"
             & "588390429865154945995197178191115892068783134378389897"
             & "809831485312898895908113774202896765076021187448111199";
    sint204 : constant string
            := "0.5882815482226453047864398132348808771714497340642036"
             & "716783410379977402712991684314747423022365860932521524"
             & "553563599245086591297387111149268132473028145909790636";
    sint205 : constant string
            := "0.5907597018588742284238879082605695714561338505893823"
             & "218551275930807683842486667611678670213622261462234644"
             & "164803768142010447659695958077037241339037850653884092";
    sint206 : constant string
            := "0.5932322950397998080478094263125274426679727983917349"
             & "942720027926289085759972243915559112345549981837045587"
             & "381074719924616012942706336918588360541129214599293559";
    sint207 : constant string
            := "0.5956993044924333434670365288299698895119263384375047"
             & "868977947345129738936287256664315349110825853769003100"
             & "949115109402206769593983646369984745112052270470367907";
    sint208 : constant string
            := "0.5981607069963423117249586521624986237340592960663008"
             & "127624009663732205155437014182821138648550182198041386"
             & "214761290188460137129959174105597437081328772865829962";
    sint209 : constant string
            := "0.6006164793838689266538758955455595386439545513218809"
             & "953821578919554504540907427737751843393314327215945265"
             & "858934856117472240085206193594067951699686182455088950";
    sint210 : constant string
            := "0.6030665985403482016934306169951154427009111934547324"
             & "024945770349033566458640528168599835260756119262962517"
             & "658281010384013079282803380722831364569429293098700972";
    sint211 : constant string
            := "0.6055110414043255139206269413298796581825983175372228"
             & "999106201976322202011281189863287846441187291491144768"
             & "902529327548649328942384523869387978363799601326420206";
    sint212 : constant string
            := "0.6079497849677736672436426710264256112166043249865144"
             & "038810334510795077714632028975622717157298095742720056"
             & "800440952301363401509647756614414322049635281090647914";
    sint213 : constant string
            := "0.6103828062763094527163521517406882689736974200427536"
             & "464451341914606160555843296894896186928517017990988353"
             & "693387566641824121223213100613749723080635504048118257";
    sint214 : constant string
            := "0.6128100824294097039352119357182669348358081983250578"
             & "624528845973159804129812242222699436589255619098191725"
             & "929940859932847990147660344848133892125946602452859017";
    sint215 : constant string
            := "0.6152315905806268454849135634139842776594300077644242"
             & "274544609298006215692903192269110107483029334541586100"
             & "695585596470941202999304458162281938483573217527500617";
    sint216 : constant string
            := "0.6176473079378039324039794017162291302857021824689483"
             & "307698284270294303693774715731033102023654274543958319"
             & "838812736584804598084988454276662950315835092285946962";
    sint217 : constant string
            := "0.6200572117632891786462681913114239630541487502709073"
             & "619040197912932014748016069110565428674319983132222097"
             & "438528919372298218860694966927203210985132345273255633";
    sint218 : constant string
            := "0.6224612793741499725191667208364860040790529154497645"
             & "824751449366331572160239533105771395867911164168771250"
             & "394018974719867934245588820067335320752002583734646554";
    sint219 : constant string
            := "0.6248594881423863770840728162810527125578830486090799"
             & "052802523184906598474771427390155904186015826647617026"
             & "029911673514394811774103683191113849024921186484219358";
    sint220 : constant string
            := "0.6272518154951441135096225651662877879063712406198640"
             & "212778253309966804512228971185442291205552774555657881"
             & "279188778552260090017495531390796834412686163033932350";
    sint221 : constant string
            := "0.6296382389149270253729813407145820848289792322209340"
             & "930839029102992388619721426669229583297251365555110948"
             & "908165368339932360188977238321854255769101683375579667";
    sint222 : constant string
            := "0.6320187359398090219094037057276920411962095885686819"
             & "279506708617217019092134143227160512764879393447329236"
             & "693861329000262629370333807299618185558429748937116184";
    sint223 : constant string
            := "0.6343932841636454982151716132254933706756870948417216"
             & "064338247186966867261235483941717674035365368534697200"
             & "773200211345837390831366418354602759052326453092637223";
    sint224 : constant string
            := "0.6367618612362842304139434349020790264832412203675084"
             & "467069993495150241594271522426303585879571135169367191"
             & "294864066920773244889594834911217035565775619559417372";
    sint225 : constant string
            := "0.6391244448637757438014881927921738391132953281173675"
             & "848149560865914244967504258663418668696043187461478702"
             & "621533231863461467143009762610058800425620851866674245";
    sint226 : constant string
            := "0.6414810128085831519887398976942527869152617376142429"
             & "124131192750126848703077404111091027113814630904504147"
             & "014763687655094282760651716129085470366369040415366465";
    sint227 : constant string
            := "0.6438315428897914650680860631769553884419570576009011"
             & "868107220178840894113968623573184042765517183682450448"
             & "352411396657497562775897490129876494644964000816913186";
    sint228 : constant string
            := "0.6461760129833163648328022195365852888821737660677396"
             & "609185129137581387590635623783601916288482155253193726"
             & "184521766962974259176327270335371302172596500582465131";
    sint229 : constant string
            := "0.6485144010221124450845605508348932123722212627783756"
             & "061610926447170148487295631545935086705193037931445000"
             & "007908909205541130526176559198014173621508160763533858";
    sint230 : constant string
            := "0.6508466849963809150689755729126486795480597434835441"
             & "866516678657682140513317415725453050144101835682623652"
             & "040969979332109209195256277297234732547788926936227542";
    sint231 : constant string
            := "0.6531728429537767640842030136563054150768600237141605"
             & "475133036582798174067452799172122116779246648060015260"
             & "864579856859481800978000664178278943480117592094318091";
    sint232 : constant string
            := "0.6554928529996153853126797012293059676952069517977622"
             & "662051420064331176687319945239297910626840576399522864"
             & "201321381812401887460807318626824450371896659145720646";
    sint233 : constant string
            := "0.6578066932970786569311822637300037965084167495920388"
             & "765198309430201853225379130424029761764228519113297327"
             & "399596575881098084770108719084505256897724373201442275";
    sint234 : constant string
            := "0.6601143420674204785594907468958086011814235036758911"
             & "913415855637633579784223784688232368289228594303332696"
             & "685255437464966910094400861663964018278153058263345643";
    sint235 : constant string
            := "0.6624157775901717611130698169566881852974425440511929"
             & "825277633930492309617675878204032686382703335051241103"
             & "811789726622294901122683449046583346924057913681758343";
    sint236 : constant string
            := "0.6647109782033448681303249852974479320666974744887784"
             & "859826918528432748115336806696568784741934041869177852"
             & "012879568406387837201809358295485207271013615879988541";
    sint237 : constant string
            := "0.6669999223036375066501542217927266025220265822372297"
             & "847236967136134419351497348575039360969155879917812548"
             & "695728478079612080433103207948244115441275757087370138";
    sint238 : constant string
            := "0.6692825883466360657206963659359292635496952150734586"
             & "104944256345583445023027408017425936581339185043079886"
             & "948720985820335257903354418306733839520813458486559794";
    sint239 : constant string
            := "0.6715589548470184006253768504274218032287506321997944"
             & "988832130731094649118790014151672904931784041606634999"
             & "442618036977416712557804483405115265096832534678659043";
    sint240 : constant string
            := "0.6738290003787560609175683717822752457223626629661524"
             & "063954096971439681915157772764276076471916838191306576"
             & "818726804259491470871459745591261248843121908795844543";
    sint241 : constant string
            := "0.6760927035753159603604192276581525580927731452766324"
             & "608588630133987883403023735467098779767622733026503399"
             & "383041087631843176711454740342176588656697748384653027";
    sint242 : constant string
            := "0.6783500431298614868736550417149607364422033657653815"
             & "716528878847599683946010960859016307172864331593347954"
             & "459439370447413089591940172289229497118424322988906250";
    sint243 : constant string
            := "0.6806009977954530505944304644563996185348546094786518"
             & "119251422058873009963578370745862507508405292100544647"
             & "750739483472645850368808805793462255083137411871602171";
    sint244 : constant string
            := "0.6828455463852480681645961230581124832961157832365145"
             & "886408415108394666365972970657536333661388673231852993"
             & "856399362785345121021406942123697647586624007859535276";
    sint245 : constant string
            := "0.6850836677727003813620525448786688535463171180026453"
             & "736937486985984569486549936815874603508114767262435249"
             & "031549368866246419210743354969037001314237583346614302";
    sint246 : constant string
            := "0.6873153408917591081991869482317427466776116545951892"
             & "451570873902027449124012540074630637094460104829544116"
             & "298456885928283651033453289770738663215570271324788604";
    sint247 : constant string
            := "0.6895405447370669246167306299574847028455368442791232"
             & "258618203924978160457288778645751098290599853791044920"
             & "716402292097131514817222691129798460371270071349750121";
    sint248 : constant string
            := "0.6917592583641577749067341320888287838269067705299095"
             & "911102144465317252155279438740777303584453751094186114"
             & "454525661629646869871671607742040001411735433013327733";
    sint249 : constant string
            := "0.6939714608896540090037343890191158870916147727358502"
             & "188425787927708104317601256810398247547895826050658317"
             & "959047665610448728670515686964997842841618985434147996";
    sint250 : constant string
            := "0.6961771314914629447885825914304004643727746960111956"
             & "668567824901156531946713995050976716898204124082610358"
             & "261379267682119479805636575235936989655810247788728924";
    sint251 : constant string
            := "0.6983762494089728535548135030617225398363276188521383"
             & "524865554105532064124710882403167866541746350997289509"
             & "601103265248176422468895195282044082833735633424012053";
    sint252 : constant string
            := "0.7005687939432483667928663802436673373586637428573520"
             & "206839033491748725997381000686901569091625534354810081"
             & "059523545103522952518474088354656031641208472320580851";
    sint253 : constant string
            := "0.7027547444572253024529144208959899415321728118683075"
             & "951906864810262181316229355224222327794546693269444181"
             & "399226491945322533783088388505865779654358104480808122";
    sint254 : constant string
            := "0.7049340803759049088525237581118148771392938248127094"
             & "351776507713929956213381924396801086866010559808095496"
             & "744752437725362683727308516154922272269289357978123894";
    sint255 : constant string
            := "0.7071067811865475244008443621048490392848359376884740"
             & "365883398689953662392310535194251937671638207863675069"
             & "231154561485124624180279253686063220607485499679157066";
    cost000 : constant string
            := "0.9999952938095761715115801257001198995529876336221876"
             & "541107369867438233011131433523910096238281304075117980"
             & "405044947912524999914655753762623999310202055873148097";
    cost001 : constant string
            := "0.9999811752826011426569904377285677161739172509443350"
             & "919401576950808210694765900655445062260270025721684823"
             & "688320219869651080170042578869600381875222785540644988";
    cost002 : constant string
            := "0.9999576445519638663331209195316303717368608499531550"
             & "935450877445651128811655924872174346546444639073767111"
             & "893453169630997624575282367272537767373291883662394376";
    cost003 : constant string
            := "0.9999247018391445409216464911963832243506064688022178"
             & "302440836730014837213868841334793501741589368376047625"
             & "856712191820356825849276874109509686173466167182467053";
    cost004 : constant string
            := "0.9998823474542125256330496265059156608292573548366226"
             & "134361084596641443980967303008001314853916560931653349"
             & "438935413501465519468406015071132852475665060715419680";
    cost005 : constant string
            := "0.9998305817958234220157222749226655145858957850267552"
             & "231789803669895796702500070420021220880128273906036424"
             & "979212831386200654845653655334494162240918881848823504";
    cost006 : constant string
            := "0.9997694053512153216576170363329930430512377077464978"
             & "564183800381509957301124499121640493646798982950476495"
             & "622091690884331858233452567964065700794539002409759437";
    cost007 : constant string
            := "0.9996988186962042201157656496661721968500610812577296"
             & "246441237879263816687547147750419016055276630732636855"
             & "255741913959553359622019704652270908328616194184995415";
    cost008 : constant string
            := "0.9996188224951785971168306373539819130463131387659250"
             & "077173328177194347814975904941807779361241973360657195"
             & "271241146753220649392753923961106698931614077660145574";
    cost009 : constant string
            := "0.9995294175010931630797033221567409693597944190919770"
             & "211776509087587144040740069696910600031388104767199499"
             & "128371718432551783033492858286378624638537230492895164";
    cost010 : constant string
            := "0.9994306045554617720190083272854685356212209690526418"
             & "374968986454828960440596227698215137450856146726631144"
             & "136593173506955823653212055921772406492571394239591050";
    cost011 : constant string
            := "0.9993223845883495008962210111399103525911368291690488"
             & "414700783462084728541753948606086943047741244546247752"
             & "525150425402250406721323980636886260726142064423734165";
    cost012 : constant string
            := "0.9992047586183638954929500005057034450442347068512060"
             & "424350343548518605611270845924963148470537207840039198"
             & "472117990071370362009854753445537643697419082171574745";
    cost013 : constant string
            := "0.9990777277526453828887819968641261429898648213239274"
             & "790163487757353577643585245814293732138210092255607248"
             & "773365678627962252735812865240556147732684878415527760";
    cost014 : constant string
            := "0.9989412931868568506339302657244654427194216613726419"
             & "750560313946908062727733108435230619910243468558963144"
             & "016761286257850713734910312328796958730775475544884985";
    cost015 : constant string
            := "0.9987954562051723927147716047591006944432036147046117"
             & "943428170736856099059521574038311892442663941710411642"
             & "753012920518949380253048468143861730196515371856672868";
    cost016 : constant string
            := "0.9986402181802652224181990491800163114016797180866067"
             & "141858045441497622818959219527227366797691180266339608"
             & "635264481060322820131985631270884937699434764617462383";
    cost017 : constant string
            := "0.9984755805732947522085590384263713301081889047769385"
             & "719273966753791254822045208650361164603629810427519787"
             & "474875849291463751931448660581845465112182671794624607";
    cost018 : constant string
            := "0.9983015449338928407387821630291092908158931696678409"
             & "510585852725396188756408788652597022456553249756536259"
             & "139950325342203317882486682318273607928146503559777877";
    cost019 : constant string
            := "0.9981181129001492071251558606836062321140632113778141"
             & "992456004556118395173713637433388786880504218344134758"
             & "266407544204991274602293829914522218134807354213256771";
    cost020 : constant string
            := "0.9979252861985960126230254623090194905748672013469828"
             & "081331772000341233114162448742677642034971763837781191"
             & "341880924709057161639297482352084193960696642959641130";
    cost021 : constant string
            := "0.9977230666441916098485467284287755605433413968690056"
             & "141311132658119507200350651581780478796198078540327066"
             & "797917348435717121668413859680351119580935073430527356";
    cost022 : constant string
            := "0.9975114561403034596994483898081441976323805321324117"
             & "935165395169731353813446088763183372462016758985638536"
             & "590586494222438991225955821495722159970377387108365192";
    cost023 : constant string
            := "0.9972904566786902161355971401825678211716886791662213"
             & "732075785240505426808841642150637983837742537777780594"
             & "907648121931286592117752244348771536463150986839597914";
    cost024 : constant string
            := "0.9970600703394829789879899493683807050365624334870277"
             & "165261250437243719240791124805542147053632904292330599"
             & "730055985331381188345307585010791733780719731254139983";
    cost025 : constant string
            := "0.9968202992911657149726293983440371335985977721848728"
             & "384982186920648330428731597999991292167879077734762468"
             & "168290967488196749554676822760362949288937862482848974";
    cost026 : constant string
            := "0.9965711457905548470935669103181862989689003194724270"
             & "469534710706320546528330397481108548428793895736022311"
             & "566723065236243782258042843656218626337656784677728796";
    cost027 : constant string
            := "0.9963126121827780126272261896697316108763234867490231"
             & "742314941322048183052286663034324352828581139530266115"
             & "339499361137423672006257261462429436424245357605214766";
    cost028 : constant string
            := "0.9960447009012519898879448102795669989215824446049551"
             & "206694524513484560542172664028865954502540119396265430"
             & "745718330229660482999919993521378803476796341766622444";
    cost029 : constant string
            := "0.9957674144676597939824956425161862154313075283397895"
             & "108131787511431958335863805779519167633221522049342421"
             & "162073319386192616725265370054096093868393771355195868";
    cost030 : constant string
            := "0.9954807554919269417691716003477196872654293223779239"
             & "383944372761458421190178183284205502344645098952402737"
             & "971598364364472531522370316540935872522019737441457536";
    cost031 : constant string
            := "0.9951847266721968862448369531094799215754748687298570"
             & "618336129657848901668945865379725290842696483902877244"
             & "931182981983461417865099877607189105687761328478748617";
    cost032 : constant string
            := "0.9948793307948056205911661067562028172051981133855192"
             & "006176197000846235633873621879262700855760064836581561"
             & "971022653583756017657107177477441949365146222698659706";
    cost033 : constant string
            := "0.9945645707342554521191062433890623879591963997686689"
             & "929497134733049764975365965519646053125923923732021389"
             & "239686013965974093242507455858116483971953831918448216";
    cost034 : constant string
            := "0.9942404494531879463584134419068988176869989775529323"
             & "156918600807582181213225653677072673736343168314268617"
             & "930921506414506004743990880827454260285572680900146734";
    cost035 : constant string
            := "0.9939069700023560415469228132477998214355953506226764"
             & "763674651868930403405085985075889181078039373099271825"
             & "577816977486659233044089260044331143156447598151510259";
    cost036 : constant string
            := "0.9935641355205953337820216973419474904243639784598839"
             & "099021476841949370539161558293602154930940078513184895"
             & "026119051905690968865168891844911739474676477637613462";
    cost037 : constant string
            := "0.9932119492347945331046010120927827784469943659755387"
             & "491148048822176047274351681217436733955647559801150203"
             & "577147384792449415431093005934236920592093259958116389";
    cost038 : constant string
            := "0.9928504144598650907935633439675960681742676607802046"
             & "356346969945330992396221153599285142673680410430885893"
             & "914298543467143410996471326024756859037727193284370105";
    cost039 : constant string
            := "0.9924795345987099981567672516611178200108206546341595"
             & "464907092780513666698265584048621348108183070855256895"
             & "261610517163417141668727170182189760817570592799036146";
    cost040 : constant string
            := "0.9920993131421917571120854453716869282950854130243223"
             & "527003934927291294853563690001933717867323780960189967"
             & "829459768651360246884263104785409447201970977929123197";
    cost041 : constant string
            := "0.9917097536690995228600499310995717595784264084693254"
             & "023626144153838723120274582293706301438463523444900482"
             & "292754535282253009148337231907881817928311956007473742";
    cost042 : constant string
            := "0.9913108598461154189573497986674833550034585194173984"
             & "391230798458510145478268668717638783054219086403442548"
             & "241938425462697343938183241390617203660937062296338470";
    cost043 : constant string
            := "0.9909026354277800251082370105274337521976237991149927"
             & "256992081100643996553792994128615689906479927039920767"
             & "795591021692784317564445054405244349341666236043564363";
    cost044 : constant string
            := "0.9904850842564570379986822425364353778167451152296616"
             & "438562948984225063638873064009260991606615290055476237"
             & "967462416992014053166584334804817442806539683980918594";
    cost045 : constant string
            := "0.9900582102622971055059064644647793938596186386186491"
             & "597796296552743739230798941892049947506658785750631376"
             & "064947051335215042077090071685325072867424775949377155";
    cost046 : constant string
            := "0.9896220174632008346236944537822198667249656540092215"
             & "756965962402622914860215558207185344718543851355413399"
             & "927592972581569434239086407081160444424448774148610749";
    cost047 : constant string
            := "0.9891765099647809734516737380162430639836895333369074"
             & "010191503644849872398622341144003243385877090478889318"
             & "178176746209443665298299387621976957746370975252269810";
    cost048 : constant string
            := "0.9887216919603237676045164854897910426445959783492790"
             & "177046789209879760118680682887001814430722574416927381"
             & "675121349900082925023049699309086508928992069550881301";
    cost049 : constant string
            := "0.9882575677307494914047925383512565371364538532078629"
             & "899363640808687435572386103632870868405863409365007178"
             & "603879751614794979760897337986868276575329976737663714";
    cost050 : constant string
            := "0.9877841416445721542309690323667278618121977647853197"
             & "604087723990264825241998838387775404939975137203179032"
             & "224245572518557642721076285050165480504853636931483705";
    cost051 : constant string
            := "0.9873014181578583823998158018450177283203725606333990"
             & "723464856282547202161680413296923198214066593567839494"
             & "814302672049182982737111733992398650406808427384122445";
    cost052 : constant string
            := "0.9868094018141854769702359522345500231768165633873553"
             & "340954664046433480249332223007001772523804215731333174"
             & "559781516602153882728947906451694319286080093092235962";
    cost053 : constant string
            := "0.9863080972445986478632975243258948530479327400784137"
             & "071185702432957470187041762660276955538002550988534046"
             & "454887452544769644761909465704276975486168985082021563";
    cost054 : constant string
            := "0.9857975091675674247009949996071128769398239356131899"
             & "318496370830312513571190685468020293868455907788213635"
             & "083950973458550760897373551459099508883020888274273455";
    cost055 : constant string
            := "0.9852776423889412447740184331785477871601291558128148"
             & "744442534259035894601318364005524044913961710172645813"
             & "729490339621374249967889805539357479462643228312957307";
    cost056 : constant string
            := "0.9847485018019042185565531758024046506443781277586363"
             & "605051367748947393685626232759770893447373403575237732"
             & "519398557636531278280699897034853488153506974183776303";
    cost057 : constant string
            := "0.9842100923869290731938743872398332594066332257994171"
             & "627485610797752359154338893194184379871857008731683279"
             & "821884143514920955761479073290483050621126550115800997";
    cost058 : constant string
            := "0.9836624192117302743962377761507951158350349472112216"
             & "953632389302090396314694800686040230697660489981458701"
             & "520749269608038139948505755658680759337483085648224596";
    cost059 : constant string
            := "0.9831054874312163271803011546739453386706019222933975"
             & "526577941595694195922260171330812506503594499562047512"
             & "196040219313048888188959025917036728127809965083569162";
    cost060 : constant string
            := "0.9825393022874412559070403955777269171829204613654459"
             & "166324238584350164968892806025018566800114381842599915"
             & "053934868553996362921116463132372765456762443105170212";
    cost061 : constant string
            := "0.9819638691095552640728481538315649202651768988504381"
             & "096017890537567277169108198098255310590228840645490771"
             & "307225266195515116675838583246726078658549006913032376";
    cost062 : constant string
            := "0.9813791933137545743182241898789480320709097648168140"
             & "978905287205202548896135627687603661635922259410372869"
             & "709676264878779898126423173742012660437686024689189921";
    cost063 : constant string
            := "0.9807852804032304491261822361342390369739337308933360"
             & "950029160885453065135496050639150649858533007632598948"
             & "662798775784681310960848381701091485451909052981223580";
    cost064 : constant string
            := "0.9801821359681173926902100086453528464396877109657589"
             & "618763359728242709151596636564777583396340507656320361"
             & "087734413494625869959685527414592148435028826369869078";
    cost065 : constant string
            := "0.9795697656854405344393261098798955052132344937016664"
             & "832965963308279433648744956009086340504990027179548115"
             & "018951496677702645785702384834566450609502297588178844";
    cost066 : constant string
            := "0.9789481753190621947154801236603956128550203763183727"
             & "240459934244801542948952559810222201621349275438624680"
             & "701677759948708870253317777641747325987002414531431298";
    cost067 : constant string
            := "0.9783173707196276331062400968954989486662500360631230"
             & "816495786015811018698451759744903891812532253407320754"
             & "949995683594864480028801618367471278574333021002890242";
    cost068 : constant string
            := "0.9776773578245099799434047624729313055872572249685078"
             & "999434977059277796082462615007685729231323363141915540"
             & "440286638582988895306281316673512424827272043548366775";
    cost069 : constant string
            := "0.9770281426577543514858662110857144252619956039677714"
             & "046346396481316448798722415793590892244462790233410184"
             & "391327185338324466681844988843168595862489460074535355";
    cost070 : constant string
            := "0.9763697313300211493127321944898351365000802823869514"
             & "605372914800699461585511590539771473127256441116325056"
             & "991560270462048421802762207460402366263261098120954054";
    cost071 : constant string
            := "0.9757021300385285444603957664195279716440122657920431"
             & "654131848601497388340819490924121720223405848031097546"
             & "663046407981482597590735321123648043591111645161960194";
    cost072 : constant string
            := "0.9750253450669941468449134678423469945313388750890490"
             & "217564694191829259000462786791014357315903556766970311"
             & "520154341694557283627586910004799155262365946912099857";
    cost073 : constant string
            := "0.9743393827855758605187216681943645931425572615457181"
             & "085075592055635059417511288773973236613467944281869388"
             & "415754556232372036308846129323885963079968196994851207";
    cost074 : constant string
            := "0.9736442496508119253183839115181956093635632173765468"
             & "164373125622288663310232607292745308186902523029637332"
             & "944977404210854582273598151621885367904820822564290260";
    cost075 : constant string
            := "0.9729399522055601454677201139037965700243897536956235"
             & "170889072514855566027279337974215662072989261852619521"
             & "956867667444228330752269100458433589402847680578701129";
    cost076 : constant string
            := "0.9722264970789363057083211442241431612174246543477987"
             & "177757107756439266440148511065719788243213119213086448"
             & "139107788660366100462035231761667231962259603724351801";
    cost077 : constant string
            := "0.9715038909862517755370996218349531511232991382862402"
             & "580222206632399626060785245148587469449986959431076234"
             & "847968790375495922685552146179777627776215383449521848";
    cost078 : constant string
            := "0.9707721407289503021381696106902808737283328953573998"
             & "065283373776964952818598324852026948909859916561538882"
             & "401569383999909828682004588409285042030238420485724950";
    cost079 : constant string
            := "0.9700312531945439926039842072861002514568659622480741"
             & "009834974506799112392593276800890798848932125239392165"
             & "702750461151864748558447393744389813988141225790573092";
    cost080 : constant string
            := "0.9692812353565484860482907381059832428042708423117998"
             & "036815894241546156347813079270456268886928875822272657"
             & "880043213836783185074221875151196906016147217782655103";
    cost081 : constant string
            := "0.9685220942744173162210883289834431285217750218321895"
             & "003012659367825618090665948528965074262923873720341546"
             & "235384941242411192457498850496205425817698807814102011";
    cost082 : constant string
            := "0.9677538370934754652433919122446032943536820495674942"
             & "035253925875806217946428974688566568351055063505094747"
             & "857299358448159645456991294374610744466910880157938679";
    cost083 : constant string
            := "0.9669764710448521090872202259367862730479839136104195"
             & "989303513041730219611083739525922661999501009065522955"
             & "892289033260477852355688934187463186728482105825275671";
    cost084 : constant string
            := "0.9661900034454125554338329612223419786343840872174320"
             & "922991761673887630240952669599693548990044195387761079"
             & "957248905702662899363490193309108552632882405303560185";
    cost085 : constant string
            := "0.9653944416976893745508438575169030701867730131426343"
             & "205542520119051642586856569821790410914804860161691422"
             & "407432868923619611421339771447521420719177767520330600";
    cost086 : constant string
            := "0.9645897932898127238364321586277055315770498729752100"
             & "975463058254447312938729556496907601426695453131598804"
             & "803610392359981067067043997717501944742253053948140043";
    cost087 : constant string
            := "0.9637760657954398666864643555078351536630838488266327"
             & "043089160414011547947386650792172514930553151595866123"
             & "289104277313575261318315272060684430627529457062679572";
    cost088 : constant string
            := "0.9629532668736838863479214808508748735204565905897460"
             & "682687117760262068253087088095449390223470855219494725"
             & "702861942010088949852782967616529514917874167439551236";
    cost089 : constant string
            := "0.9621214042690415954296043162301533683259328575220033"
             & "671705138807821048031823343814048250043454930134788895"
             & "545768029996974527354458910297296282875199311774848673";
    cost090 : constant string
            := "0.9612804858113206417486596525191235450943555695075193"
             & "414309501472414000188115250119301129932441617631827040"
             & "776129881443195618116971562133632822954049915467895033";
    cost091 : constant string
            := "0.9604305194155658111990351376552656335496268220003038"
             & "114349196996602195860983666547482838877672058391302116"
             & "193311886523342575990980422927759357493161554506715386";
    cost092 : constant string
            := "0.9595715130819845283355281812303626134485965186116516"
             & "972353870707142388981665457541897260610856850318891067"
             & "501494265122392465594008443460070852214300889065392844";
    cost093 : constant string
            := "0.9587034748958715553746457917668690903737111090162149"
             & "144316197753884929659828284853532259509025070689421279"
             & "574207141621812708310802615855474367464078260050306846";
    cost094 : constant string
            := "0.9578264130275328903210370287966757208083863485438522"
             & "062087553172321154698789902443960686854432584119800485"
             & "449246420489158184398767090607894457241012704733297446";
    cost095 : constant string
            := "0.9569403357322088649357978869802699694828492056300372"
             & "613012071998841601453681608249571934246962711819922698"
             & "493296146407927856166538403476343873898787038483921127";
    cost096 : constant string
            := "0.9560452513499964432704798225393117130005746467262919"
             & "480204953574988727611552191831339080422060870827081434"
             & "477470634912163597183392771004520010605784706744426814";
    cost097 : constant string
            := "0.9551411683057707214981577123356394246483804498203158"
             & "529646521887671458409193643991259473977487539849435622"
             & "826524164591301473623416405255290739035812991979941060";
    cost098 : constant string
            := "0.9542280951091056297804307321904861426773098032184780"
             & "401644871465492339712759831325969298147172342623947736"
             & "803642481726033045773183251931480787730797727260345751";
    cost099 : constant string
            := "0.9533060403541938369167403827397938349040969010171632"
             & "867949529339188348896363636206990880373294669276509252"
             & "572913893836171225334338271505790460258061567090335585";
    cost100 : constant string
            := "0.9523750127197658585298936075710087775910962246224436"
             & "027367604950942596168858603174563383596614003172550013"
             & "590984907991364131254913524556327349749967706363191239";
    cost101 : constant string
            := "0.9514350209690083695491755689557821707500772607476043"
             & "079392068176881272993806457738474895188569845765008128"
             & "961530878853247035986429508819934679121008756186987956";
    cost102 : constant string
            := "0.9504860739494817217599261006205415490837005595036579"
             & "598658877148852072280006495223904425141991699141856344"
             & "011234229157694408761497404703246004046906927891700599";
    cost103 : constant string
            := "0.9495281805930366671959360741893450282522241538324108"
             & "524439709653973932083707234524215787997564023190217594"
             & "060381506442301498403297053621440694324110218388017346";
    cost104 : constant string
            := "0.9485613499157302881584948257653042498390321546836932"
             & "216198241182734181066441995003027368498367635506693068"
             & "170640346930117303752068893748922119846083671779405012";
    cost105 : constant string
            := "0.9475855910177411346533873212314649157949189830936963"
             & "643808307761880154723549967661047076255098680343782408"
             & "329870219563520003003446625318713874104284489098441144";
    cost106 : constant string
            := "0.9466009130832835700445998229621097795147633836406999"
             & "140748968915553641624483601833887564938312934209814167"
             & "196816630770799048435887144599691115596487590470997106";
    cost107 : constant string
            := "0.9456073253805213257309453865238450864520238660965341"
             & "455629718085724350291286140933010689868381632885393497"
             & "071967001543043999345728765389608974707499477461498314";
    cost108 : constant string
            := "0.9446048372614802656592654934661592984649847828579160"
             & "291202279954325714120987907086831346001983686749156146"
             & "103124427974532948537295428686786222221556389598157775";
    cost109 : constant string
            := "0.9435934581619603614953014453784386693433755609167296"
             & "540524324871497281597136463378068067271732514632983988"
             & "873848100704043524874787401245106704428865783954892414";
    cost110 : constant string
            := "0.9425731976014468792807587350218082231326560438123738"
             & "601186867913336940993537302886286541050626022328707271"
             & "042244841108527402025371096107592677010749286745416226";
    cost111 : constant string
            := "0.9415440651830207784125094025995023571855897958251828"
             & "675468258789699971203630496978846001890548036847188371"
             & "045666972199780856472480360553974495519123755747961701";
    cost112 : constant string
            := "0.9405060705932683237872913092520213689163721527655615"
             & "618922449578789587957300079664333610328080693847265879"
             & "187455280033723136566760870170901780571025833550123693";
    cost113 : constant string
            := "0.9394592236021899119626692458704222131583983521196722"
             & "128837739195491790912870070251174422800036033831138594"
             & "089180324772370327292565646115271794692284301075110569";
    cost114 : constant string
            := "0.9384035340631081121924207736047628846643663901209036"
             & "726806226122421818740213601685058494060320800211099946"
             & "492216710648648637025282833745645162809297389101614150";
    cost115 : constant string
            := "0.9373390119125749232018995933723808790288308132995546"
             & "963325657458411778002447696447032083161920987636106024"
             & "800686367006129621211582352120543261802601092901372093";
    cost116 : constant string
            := "0.9362656671702782465763109956857752575590818942129451"
             & "192645921308548966316671964632560149229495942432466206"
             & "336423295286475685661589549239608080638624193149618098";
    cost117 : constant string
            := "0.9351835099389475776422074797424988136309455357656205"
             & "372815936361740075933392808785627155405552851652442630"
             & "934265715323015815117931469537084349952485935705006433";
    cost118 : constant string
            := "0.9340925504042589147298778825471074110553903629225686"
             & "445110875327636676881405172620186886075039348755511766"
             & "167491278739005701322617707214672493910674629504458570";
    cost119 : constant string
            := "0.9329927988347388877116602555433024982950155205122950"
             & "488923147712766684426703203108827359068660545814659569"
             & "432560336407190523005268397275774895773985253415942843";
    cost120 : constant string
            := "0.9318842655816681067185571985770264225996390618719338"
             & "789353596792493529343904461903246532570809745433789965"
             & "979393590490867259922579310006773901396822629959595753";
    cost121 : constant string
            := "0.9307669610789837319448723398218081380297515646350503"
             & "346297894971826961242543993839102265110455496911467003"
             & "639029353153841378388064160058413847744218621672080784";
    cost122 : constant string
            := "0.9296408958431812654579180664894332423899170154979732"
             & "937462342707783857564402827043132235152837168917667843"
             & "601635326851408785943104927022830404754987214434651217";
    cost123 : constant string
            := "0.9285060804732155659371673957159456406234109867828877"
             & "479861615510769100582911511828005233458856097702717418"
             & "630343658668283731906080996935255721239666335111882702";
    cost124 : constant string
            := "0.9273625256504010872745369590302413069028138258771417"
             & "244020837490472943330669356551880528131282370273941772"
             & "838448474171009009393309223139221586885216848574681891";
    cost125 : constant string
            := "0.9262102421383113419747933884371432823411412506173913"
             & "153462233814000053246707029808628138042140535213919362"
             & "929011240216327826317891709572055238713450528008607986";
    cost126 : constant string
            := "0.9250492407826775903023718686184477409811358598791183"
             & "891915584516614361127689994280999165891352126804051199"
             & "332561263768040491751050419262693000841448494785540228";
    cost127 : constant string
            := "0.9238795325112867561281831893967882868224166258636424"
             & "861150977312805350075011023587148399348503445960979630"
             & "257822478830308691775799042014275332219995578278983938";
    cost128 : constant string
            := "0.9227011283338785704372642268248986908391320246055309"
             & "606847655752376964394813901218803654534356020026270710"
             & "758250956768690826558758442194204204385066230942094895";
    cost129 : constant string
            := "0.9215140393420419434653963315480622263672465115990260"
             & "347766029679984936339535886574875348797982385112818102"
             & "009493514376737352339371821916229709747935375923293331";
    cost130 : constant string
            := "0.9203182767091105664400765410164427383299927934101930"
             & "909504779609513481251059397152247817000497521287604495"
             & "910517751570283375685078534248854817400476223919263791";
    cost131 : constant string
            := "0.9191138516900577439084777893585916788760254170155826"
             & "476013631182379061388510154227887534912418462755963565"
             & "031760114255395288822124178644437523686130697308191512";
    cost132 : constant string
            := "0.9179007756213904576422762970161218427527792101607703"
             & "185971662225004315031684208421819211416714707249659731"
             & "243077828120808897256767198429723187607886354060559830";
    cost133 : constant string
            := "0.9166790599210426631164570134177923250274471032269164"
             & "055587242622222926817895136578874369584703071366161133"
             & "635270531032986293594408835338321442614595982894837225";
    cost134 : constant string
            := "0.9154487160882678195664312919622163783863060352317183"
             & "605341850852170136655263558567909094958677925637054898"
             & "529874722935142519247384985389511124572545787395947212";
    cost135 : constant string
            := "0.9142097557035306546350148293935774010446911156821770"
             & "013565662410591304975518946726379074291674786932682938"
             & "611970114452390696598450028374233207993449586046161934";
    cost136 : constant string
            := "0.9129621904283981646280182333945881639288239722808657"
             & "867041241923202400392828472730809165261553569508411384"
             & "188554634255113090726752087766639352935909359587023054";
    cost137 : constant string
            := "0.9117060320054298514043973250755404989166873134114448"
             & "679139587726688467450341916215372487880953275350734449"
             & "029108069701741902307232499220158583679788263185527926";
    cost138 : constant string
            := "0.9104412922580671969340953692880071699802655492372571"
             & "006413260539468868768253516712935180550071088981642490"
             & "061003236003411218499297743386388560440866363379201044";
    cost139 : constant string
            := "0.9091679830905223765638847877078063304794860514149418"
             & "013493945486915413433040565043848415036564172092068602"
             & "274258642336248010015988971781204759789307377243686338";
    cost140 : constant string
            := "0.9078861164876662120386814798769818177571351230944045"
             & "028629650232143051085886448449027401107299512114804836"
             & "802691744754718201181112208358112779749186864098740127";
    cost141 : constant string
            := "0.9065957045149153653329605884237134126506369603843768"
             & "155054695993610692483759009919065695462661868140788511"
             & "698218036608866630200951702020083868998448130230247058";
    cost142 : constant string
            := "0.9052967593181187743540483291399726543292442596242562"
             & "214947388471444646149536217812978623617410026166949741"
             & "810274832805908461140320228406725219037503536989274690";
    cost143 : constant string
            := "0.9039892931234433315862002972305370487101320250506080"
             & "496646735759958654405663192393572459493188167709686373"
             & "412858922707859253343703234618218445090546064837391161";
    cost144 : constant string
            := "0.9026733182372588067515023906888894471608558820282872"
             & "818208427247102924789114202515037631176859530826863072"
             & "234237108576692885412297515227043581242041317696618685";
    cost145 : constant string
            := "0.9013488470460220145707460933354787941262707344023840"
             & "430179336531474003903532558241555053161255234923638620"
             & "946759720839918452658395888714256390824053396741922549";
    cost146 : constant string
            := "0.9000158920161602287145352665970513170516412731865543"
             & "879793180050977045286955066038426810044169324302716509"
             & "044047758313084815333473895846061407344841519427170781";
    cost147 : constant string
            := "0.8986744656939538430419767437334850971686375082110164"
             & "316056962924551981783325052487723219296180128524264508"
             & "501372651635826358481686922415569591622209021735482263";
    cost148 : constant string
            := "0.8973245807054182812313918361485724100227964935611997"
             & "250374640428573507988124307865275701413632768739634443"
             & "684821500488585011281460344587293677179773958640241516";
    cost149 : constant string
            := "0.8959662497561851559145602819685074912350183654394935"
             & "322797390060548191757512708320552930334719315116383282"
             & "954855960719239546316518094141468289325540929722612004";
    cost150 : constant string
            := "0.8945994856313826784330721256493119981033722667560281"
             & "534398173861636405663426352853061555401913127539125011"
             & "229424759933894216042650654370393721585120774479602027";
    cost151 : constant string
            := "0.8932243011955153203424164474933979780006255889988727"
             & "896079334615180005880405975123811164719357263991195187"
             & "766929392076457777412648507331833725128916872764364733";
    cost152 : constant string
            := "0.8918407093923427277964786972263580580505297928777078"
             & "756047695710567971745844367413729639262820804910986579"
             & "633662184349547246249321321275924712903683690004237603";
    cost153 : constant string
            := "0.8904487232447578899521505599180370203448336066920747"
             & "170908443080644308771945984098359575822232874916336758"
             & "053355711883600666491396597653030615421281990077705412";
    cost154 : constant string
            := "0.8890483558546645625407777293374767964898556757438253"
             & "030803481144320686554397899980082427075664283992961044"
             & "178603066188631928990097906713501803184482077373306870";
    cost155 : constant string
            := "0.8876396204028539477601816172209067692593291468142088"
             & "200144762682216257218234955710341178697470831918157328"
             & "295087262605310090189786060451396876998788811085926119";
    cost156 : constant string
            := "0.8862225301488806316479908209186334139362462577180094"
             & "050642547035792985516040336626859369658422637591303295"
             & "486232801297491900231883467693215529293305297239569905";
    cost157 : constant string
            := "0.8847970984309377801040070405857408897590863819469723"
             & "775926315124450646871598995194801771557952983844030665"
             & "668740156902220224770101743965130196734257678842975234";
    cost158 : constant string
            := "0.8833633386657315947363080147110533129482287388003985"
             & "211877201370905411934243708024170510672647602797983291"
             & "728415926548255064364204197109382693187094829444270035";
    cost159 : constant string
            := "0.8819212643483550297127568636603883495084426206747279"
             & "806325386167120666470450034970776205819055333733858182"
             & "337103195704203212736540957959522145572052095447098687";
    cost160 : constant string
            := "0.8804708890521607708065429294693790235531168140817586"
             & "696008925956573889244321730609730506078150494142579117"
             & "597149350001483360848768538676238457331918091610516743";
    cost161 : constant string
            := "0.8790122264286334778313237108883723441412213430684602"
             & "518469326246002622438112612767452207407805984759724359"
             & "400005161495820232508892177440625370557720355402292896";
    cost162 : constant string
            := "0.8775452902072612916684707502927493687195095184952529"
             & "155629858048283418734779969913947608768631341774315837"
             & "977144313256104496284145404622478506504866151533630085";
    cost163 : constant string
            := "0.8760700941954066070958442682679904961152886604629681"
             & "197407004228278540027820264873629675303398073130325121"
             & "817520452017768736484533940854332242980139346330648742";
    cost164 : constant string
            := "0.8745866522781761126344318973080041997486952406045097"
             & "017221737198284701982278791467158548106604733500243706"
             & "860395435591659391466244834954015991879172346620612129";
    cost165 : constant string
            := "0.8730949784182900986360859730828656214301581468010053"
             & "820263561488419767988141208542808846027449901986621568"
             & "637702861575167992616473675334434038370663675271947537";
    cost166 : constant string
            := "0.8715950866559510348424814352011506891190164630450746"
             & "728682633457394707186954452361118727204777353454175325"
             & "625503650235158682768758781400007892028624623180447275";
    cost167 : constant string
            := "0.8700869911087114186522924044838488439108277895298254"
             & "871093938228539361812818416287793461042575249436405558"
             & "234224224766132274228313542272339368380284157170666654";
    cost168 : constant string
            := "0.8685707059713408953404498757617203040593764213300616"
             & "888403051165700583584971260315630190900402787686606141"
             & "138178610181801520822263318866364454320785757298665628";
    cost169 : constant string
            := "0.8670462455156926514801956294958485611965851201621297"
             & "660798725181487863996483140811564586056150781468485088"
             & "891407574136942268928652325669606270709248056079205756";
    cost170 : constant string
            := "0.8655136240905690828254883576021393973799296090407015"
             & "002253905303982477455750015399792114601557558816311728"
             & "976300416185433815982839302068394072690010084502369510";
    cost171 : constant string
            := "0.8639728561215867379181470543455525270742850950242317"
             & "811555335678079167517202899177077972415950589790376901"
             & "322859688742463731012624713627807876094611736940382185";
    cost172 : constant string
            := "0.8624239561110405386909338777812499507348754626281270"
             & "153948166998718123396320322869767530182078194690535814"
             & "356515798814791806043933767343433983971726827247913537";
    cost173 : constant string
            := "0.8608669386377672793445838767919507564613366474625723"
             & "387199471186396489662014555571438079539062719344711524"
             & "565345691568293408078626674836309963052773137117358184";
    cost174 : constant string
            := "0.8593018183570084047835821392505594300776504877935064"
             & "551046923612187469093103698470644702895615236871627258"
             & "739117232923418895798382057773699086770035917475591656";
    cost175 : constant string
            := "0.8577286100002720699022699842847701370424907994337340"
             & "186047185423495190873294644473503467443827372978075241"
             & "933570661095405205680771980831058388087056657416312480";
    cost176 : constant string
            := "0.8561473283751944810196307322098535510157021598510809"
             & "975165737872013594207560187965013852508157240295969495"
             & "589085045055356485352241037344562891506790312370270095";
    cost177 : constant string
            := "0.8545579883654005207678622757156412315482088645418588"
             & "110512777579331631426671302978508560211359238233785652"
             & "148547697285290049303292081949746679019154961550923501";
    cost178 : constant string
            := "0.8529606049303636577465880817472955612247429334215427"
             & "919944020427523011457280123176945055587586258150967346"
             & "884391913724395546699436321107760990816403452140254557";
    cost179 : constant string
            := "0.8513551931052651422612903117258743556752017362384884"
             & "190667743706881798048993733629191335831105877273749353"
             & "298792744066798179376960853386224848393073801759240725";
    cost180 : constant string
            := "0.8497417680008524894712683949492727307512283411144292"
             & "886330832639998779886779551721857643521885904581879050"
             & "984940178687228307032853146962157971311277727302227335";
    cost181 : constant string
            := "0.8481203448032972512791335629489720926244788385512470"
             & "982499046933124708185509683309949120128627571915130218"
             & "316031682816771366873017819651088760664528186826719758";
    cost182 : constant string
            := "0.8464909387740520783005444881226811347222124933243732"
             & "803616230977976007739805772335103119046083987670200933"
             & "555266817880758564577811191835561922760351072653603753";
    cost183 : constant string
            := "0.8448535652497070732595712051049570977197859813891086"
             & "266626210143385371786920747277456118093123217801672284"
             & "127815307980250456681241186968318046987479170149687666";
    cost184 : constant string
            := "0.8432082396418454371617438651833876812949431645086799"
             & "259749346880129364719994550891835869905317673165238820"
             & "679271058907453505840236527843443102797466841136755189";
    cost185 : constant string
            := "0.8415549774368984096034995198422270646026529096920897"
             & "822933798480571048569463580896837136732489074750928916"
             & "766346523473312179072962381478955311593518841372103857";
    cost186 : constant string
            := "0.8398937941959995045833839865636283006705890551261649"
             & "383618330564158183111843403396582187935747266575296388"
             & "100943757636812646643192896759091248597200636005404555";
    cost187 : constant string
            := "0.8382247055548380431869968558042860626362709704203287"
             & "236075646985850546882622337805329101311213571178587164"
             & "549514256673277018441313948628286024228673699612373821";
    cost188 : constant string
            := "0.8365477272235119845242857901788464603635954918007066"
             & "490321657778572992257517778545426331240553942111453337"
             & "458115429376977130624425004009409809485689035921739375";
    cost189 : constant string
            := "0.8348628749863800563044013830288509569496169938396501"
             & "691698247621759545705834412655429270889159855589619476"
             & "910317250613553453935567219791610215498667567275981229";
    cost190 : constant string
            := "0.8331701647019131864399159215673905735273240019246449"
             & "123600064517642821648358427136153090160455324920141559"
             & "552077003278017431320396872868721318555891637550434894";
    cost191 : constant string
            := "0.8314696123025452370787883776179057567385608119872499"
             & "634461245902276379201446423439817749190079806503360029"
             & "402384434497288691970476296799430585600237731382115302";
    cost192 : constant string
            := "0.8297612337945230424690237646868315884264976428624741"
             & "896421886786974252694594756932508255874420987909567013"
             & "516541602032822213007055640646341451236946411785804663";
    cost193 : constant string
            := "0.8280450452577557520675275919240992974059373827702435"
             & "013402251673004249376370670847831166593190374881776133"
             & "698714113092717035508578000835834315107376592811548869";
    cost194 : constant string
            := "0.8263210628456634803111954517573735930748846670815290"
             & "526874387859085748608413637484304305595770407856863984"
             & "003157901665804313354600797205983924603387421868486277";
    cost195 : constant string
            := "0.8245893027850252644748037370848904463286468248009647"
             & "017983878967148285245862882494503860782365080786744153"
             & "591441870277785237871272095706568443595470149557667723";
    cost196 : constant string
            := "0.8228497813758263320467800344530491588480961622299734"
             & "493920894871901082223292255327703007795221417197598951"
             & "416093810991185342867924671869863499926815324217659678";
    cost197 : constant string
            := "0.8211025149911046790604308203298999031279085421450958"
             & "038286485691644319044054916388437027171445749242223983"
             & "561181943393282397082563412587060825660418888681602100";
    cost198 : constant string
            := "0.8193475200767969608246896372425308158363671409975222"
             & "751515759479898122243335921307009387027201996769810715"
             & "805870996139514132541894343134371980348305916089332126";
    cost199 : constant string
            := "0.8175848131515836965049208841306338094710425175669140"
             & "941589457011734569854789100989156894233606203180184873"
             & "664026944389114208199220536843189171199820539095831802";
    cost200 : constant string
            := "0.8158144108067337890107726598636859187152307591231952"
             & "886209737757258913335176731035680630728185876238524870"
             & "970146697720946384311332405963558723141532360692998610";
    cost201 : constant string
            := "0.8140363297059483616545166896872007386115265776746003"
             & "774044464669170594417949778338506169691073909138906816"
             & "176203157764999551897059034008881773435794550300479938";
    cost202 : constant string
            := "0.8122505865852039130497441807454391551551357211260330"
             & "985546644175991769320472814227448099419441185352713506"
             & "799518619654331777976069908300935610688145295127187010";
    cost203 : constant string
            := "0.8104571982525947917267034342445271670501183629908300"
             & "867998224525663956923539724574114552269252066104086877"
             & "539263619116755142330548458119971504775022574541741728";
    cost204 : constant string
            := "0.8086561815881749919469681278720381651798162762814637"
             & "093736461951585846591523845537267628945522920492806004"
             & "259693272600790253879857786242279611896846717024632206";
    cost205 : constant string
            := "0.8068475535437992722065143125084730510426323611548135"
             & "326799453325201134607217077195700453616968496817170194"
             & "097655426935719583202649925697054916621689643417554499";
    cost206 : constant string
            := "0.8050313311429635979226592819157852444787123174644970"
             & "190924831606995920885418544645815127031817692053367970"
             & "235907538532236848955656190584681607442292605516876043";
    cost207 : constant string
            := "0.8032075314806449098066765129631419238795694271704608"
             & "349765046543457740533480153272939564002841188636971250"
             & "153578526552530010885661343062055398088116519105861847";
    cost208 : constant string
            := "0.8013761717231402194302477771779672800629108667097545"
             & "508088530614950692727725548326704591481024484530711160"
             & "868288867797317350303472555809945172044613475560915258";
    cost209 : constant string
            := "0.7995372691079050335002462322515077643195102380441946"
             & "207778170017861806312346974341810108662847197345623768"
             & "915238715133112538694287904190134057758494284153573435";
    cost210 : constant string
            := "0.7976908409433911083626627549768352735344843743754048"
             & "329459382966025495650914084700581212228151074172464570"
             & "943941616037990599135563457366014255013219791442050942";
    cost211 : constant string
            := "0.7958369046088835362627919154816773610504779970095616"
             & "551123623411140461660007086782718225784605044266041172"
             & "990870198131525633514907456420893344881277651020550193";
    cost212 : constant string
            := "0.7939754775543371648950837572017825000808051746992113"
             & "525926753850058151174287248220758884883389652453035952"
             & "544896981820160492333033273433833218163143909037936160";
    cost213 : constant string
            := "0.7921065773002123517823428786210334566942030608291357"
             & "291633653650859464030368960855346219609943892685144273"
             & "844792870410830375068049655087346761471843072158088417";
    cost214 : constant string
            := "0.7902302214373100550302171523164021706369885772062786"
             & "870431947462777460459496932499438892559012534068801756"
             & "679515681621336886574275409276943229739710840863706545";
    cost215 : constant string
            := "0.7883464276266062620091647053596892826564931371496486"
             & "506948917380682596142739840278122197078706868444118194"
             & "652198692285559603098321321285698681291877082575229653";
    cost216 : constant string
            := "0.7864552135990857575223194638513642098426265082159146"
             & "577896986477583694630211488069540998437191961291223603"
             & "950270116029744001396408473888470479973323149981649742";
    cost217 : constant string
            := "0.7845565971555752330238925746397839911870155934056389"
             & "736686405662061245531157961897734539388161559437799075"
             & "059399312373791264536599711222251327760093927155326267";
    cost218 : constant string
            := "0.7826505961665757384589493005947525340922521134473689"
             & "623010314942817276545016964798441487213708274377967957"
             & "082220100633201521973210983610987892650128690633330499";
    cost219 : constant string
            := "0.7807372285720944783015884837795336966664392172732182"
             & "218575541571912949186280102646373014695612349156112113"
             & "482611362136388953620037599936843188763244034234030135";
    cost220 : constant string
            := "0.7788165123814759533747243252609644040752852140782341"
             & "266186630414670788747709445718121059506246577092029584"
             & "415733180816147319944829799435477522539652947557146459";
    cost221 : constant string
            := "0.7768884656732324500408279830138537082772169948111833"
             & "586492039017904288027304081652093260699741500773302338"
             & "513387247644368291355576552858279287047437866127695877";
    cost222 : constant string
            := "0.7749531065948738783591292824542866253702789416274744"
             & "528423192510499442521528468115690457394934866699693413"
             & "736913174207723046734513046686499067651910390339631629";
    cost223 : constant string
            := "0.7730104533627369608109066097584698009710412929008096"
             & "093564028966879506053059873020379548140068698393742614"
             & "841345656141239599015933191604613844325890155820694712";
    cost224 : constant string
            := "0.7710605242618137732006057586124016581017888104106236"
             & "957344745001546892119844197616344638790252564130343893"
             & "163048287054135342585468538861934390286780304564310775";
    cost225 : constant string
            := "0.7691033376455796393466260688578576671915689470561573"
             & "767814460706652638146767956577733764554233516231180327"
             & "946495696356196159565122469911571754039372493912027922";
    cost226 : constant string
            := "0.7671389119358203811816945732593211757841089136593752"
             & "564979783986311402782719432800690089805888915770964238"
             & "084871983571454629579077555039203396915710348578217942";
    cost227 : constant string
            := "0.7651672656224589258888159990649059180489680014138812"
             & "013888589144682411874807370646026874414570169458384590"
             & "161001362772241935555694149962299107816061549568773434";
    cost228 : constant string
            := "0.7631884172633812717048382970658654583073195549428169"
             & "546356156115837775787597651957490307989885858463266061"
             & "093738270172833041986570791411085455240872898105979922";
    cost229 : constant string
            := "0.7612023854842618140297098355118759758054588851057932"
             & "668222299793647276679764876301976803922542436070157902"
             & "264777075658697396242589931328292097923392271575120218";
    cost230 : constant string
            := "0.7592091889783880334855254426954076862902140668166461"
             & "023657046837885468643354082282363798126974557708643235"
             & "337239088695687528975677039675099168546663350791445102";
    cost231 : constant string
            := "0.7572088465064845475754640536057844730404337157316168"
             & "500555117689736318985200425229373510021761200010792622"
             & "571275917205513778900219117993025457748532455374704743";
    cost232 : constant string
            := "0.7552013768965365275987107562452439785001995161895532"
             & "408754923372853750076790503787509113637768265447350714"
             & "715893963324234637813109439851455629690502119616586675";
    cost233 : constant string
            := "0.7531867990436124824834304856149119103704206458977201"
             & "688199814238436253236772063638125057986905394760682303"
             & "663068979816634917333308120069765761277335181710956996";
    cost234 : constant string
            := "0.7511651319096864112058194217842731511004673418053341"
             & "235480639161350598565816620550059271121534413610077120"
             & "229376113427063999134073311940308893910124263054723680";
    cost235 : constant string
            := "0.7491363945234593254692032567668843148057429825117396"
             & "427265396568131702532271355663097737255290303954302284"
             & "845481425704745725791008995315191669647628894885130592";
    cost236 : constant string
            := "0.7471006059801801443230788471989288460126118161006954"
             & "948669300443836435085533034929445145991521141590539775"
             & "544396063187312614990823528550362087255118856478383848";
    cost237 : constant string
            := "0.7450577854414659624079073102652695805364073315923834"
             & "408792947418477155037262988420534798000506329271671063"
             & "501378323616602418026806745887622256925022462522881451";
    cost238 : constant string
            := "0.7430079521351216935173622932982347968199697288722255"
             & "960808796217890160970908919259333209551781889205644885"
             & "132582932926955210449387882453646229767194942704097784";
    cost239 : constant string
            := "0.7409511253549590911756168974951627297289553093090900"
             & "457364120438466825355639597576967103438140682232065883"
             & "642126995332937542226601359824256320509819474040054556";
    cost240 : constant string
            := "0.7388873244606151479331165079192798134314610826232796"
             & "435862715351515555688712044074906819612060735618750814"
             & "918894902791475724825261302909442893734288696137938294";
    cost241 : constant string
            := "0.7368165688773698750901325201727469468678844583869454"
             & "714583259896983881579200727962109789520953452696376811"
             & "865263815458028203372031667581338514096314513984534943";
    cost242 : constant string
            := "0.7347388780959634645632236038195365703151066276956077"
             & "451749980610634161053426603278410790659650783109111981"
             & "021484960037020256529977576329814136485864733516315789";
    cost243 : constant string
            := "0.7326542716724128346155466488994934632929541703531315"
             & "202018267179836348103639635123735795806023642575919801"
             & "044736355985436647018667509777873647669258392420223803";
    cost244 : constant string
            := "0.7305627692278275611777588499757241467643701736445517"
             & "067231339860665244616691304998033233471852716385818866"
             & "236733863542060798962517124566870653937734388212348170";
    cost245 : constant string
            := "0.7284643904482251964920354375100570986874932060934284"
             & "710721012965549394854709948367875391766589461991270418"
             & "515452596884368954019335324379030382596624790224992085";
    cost246 : constant string
            := "0.7263591550843459768174943145333992297606624731249970"
             & "521865590593986899205018907678919699763277902271624006"
             & "232812226985254659142758536174503036125514162673896429";
    cost247 : constant string
            := "0.7242470829514669209410692432905531674830930048004368"
             & "801650713787740884271114355639629308132763842398050531"
             & "783730636670194118500994201534533574744390214268391804";
    cost248 : constant string
            := "0.7221281939292153212436071976676250998279170867769836"
             & "496209438781382698583633964719256070148798860019224375"
             & "686390535053150143538824036674811823412636035894879789";
    cost249 : constant string
            := "0.7200025079613816290766829987842242033652881679560536"
             & "078689347656127839871312872838425296149569909680334951"
             & "399313675775948716762484993235014432810392588381783942";
    cost250 : constant string
            := "0.7178700450557317362113253293369297123442181563095726"
             & "148884065016740071136026826328360396420438539133177318"
             & "139657416671738505589733005852735839329063103906027998";
    cost251 : constant string
            := "0.7157308252838186541255326234552024122145960838025807"
             & "150765644749748984228732427893691391169957748675480690"
             & "950464559258239573343440357467530335143757047514471953";
    cost252 : constant string
            := "0.7135848687807935929031250994722950161934721469966256"
             & "124426691340233749164458502113745317734878817331891444"
             & "960841829336511702840732662368863020819875540290923209";
    cost253 : constant string
            := "0.7114321957452164415221302897745546729364570497435620"
             & "758261057567115582716221411872946223730266309968851886"
             & "195702205050546319843337468277673990219056645176479905";
    cost254 : constant string
            := "0.7092728264388656513165337715826299609990577295530415"
             & "148739423231633240870745750155623682388672793876320376"
             & "776854071233970627638237615713485708291508646592014686";
    cost255 : constant string
            := "0.7071067811865475244008443621048490392848359376884740"
             & "365883398689953662392310535194251937671638207863675069"
             & "231154561485124624180279253686063220607485499679157066";

    x : deca_double;
    fail : boolean;

  begin
    read(sint000,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t0  ",x);
    end if;
    read(sint001,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t1  ",x);
    end if;
    read(sint002,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t2  ",x);
    end if;
    read(sint003,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t3  ",x);
    end if;
    read(sint004,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t4  ",x);
    end if;
    read(sint005,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t5  ",x);
    end if;
    read(sint006,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t6  ",x);
    end if;
    read(sint007,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t7  ",x);
    end if;
    read(sint008,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t8  ",x);
    end if;
    read(sint009,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t9  ",x);
    end if;
    read(sint010,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t10 ",x);
    end if;
    read(sint011,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t11 ",x);
    end if;
    read(sint012,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t12 ",x);
    end if;
    read(sint013,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t13 ",x);
    end if;
    read(sint014,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t14 ",x);
    end if;
    read(sint015,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t15 ",x);
    end if;
    read(sint016,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t16 ",x);
    end if;
    read(sint017,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t17 ",x);
    end if;
    read(sint018,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t18 ",x);
    end if;
    read(sint019,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t19 ",x);
    end if;
    read(sint020,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t20 ",x);
    end if;
    read(sint021,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t21 ",x);
    end if;
    read(sint022,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t22 ",x);
    end if;
    read(sint023,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t23 ",x);
    end if;
    read(sint024,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t24 ",x);
    end if;
    read(sint025,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t25 ",x);
    end if;
    read(sint026,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t26 ",x);
    end if;
    read(sint027,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t27 ",x);
    end if;
    read(sint028,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t28 ",x);
    end if;
    read(sint029,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t29 ",x);
    end if;
    read(sint030,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t30 ",x);
    end if;
    read(sint031,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t31 ",x);
    end if;
    read(sint032,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t32 ",x);
    end if;
    read(sint033,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t33 ",x);
    end if;
    read(sint034,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t34 ",x);
    end if;
    read(sint035,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t35 ",x);
    end if;
    read(sint036,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t36 ",x);
    end if;
    read(sint037,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t37 ",x);
    end if;
    read(sint038,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t38 ",x);
    end if;
    read(sint039,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t39 ",x);
    end if;
    read(sint040,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t40 ",x);
    end if;
    read(sint041,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t41 ",x);
    end if;
    read(sint042,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t42 ",x);
    end if;
    read(sint043,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t43 ",x);
    end if;
    read(sint044,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t44 ",x);
    end if;
    read(sint045,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t45 ",x);
    end if;
    read(sint046,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t46 ",x);
    end if;
    read(sint047,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t47 ",x);
    end if;
    read(sint048,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t48 ",x);
    end if;
    read(sint049,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t49 ",x);
    end if;
    read(sint050,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t50 ",x);
    end if;
    read(sint051,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t51 ",x);
    end if;
    read(sint052,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t52 ",x);
    end if;
    read(sint053,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t53 ",x);
    end if;
    read(sint054,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t54 ",x);
    end if;
    read(sint055,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t55 ",x);
    end if;
    read(sint056,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t56 ",x);
    end if;
    read(sint057,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t57 ",x);
    end if;
    read(sint058,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t58 ",x);
    end if;
    read(sint059,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t59 ",x);
    end if;
    read(sint060,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t60 ",x);
    end if;
    read(sint061,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t61 ",x);
    end if;
    read(sint062,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else  Write_Create("sin_t62 ",x);
    end if;
    read(sint063,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t63 ",x);
    end if;
    read(sint064,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t64 ",x);
    end if;
    read(sint065,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t65 ",x);
    end if;
    read(sint066,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t66 ",x);
    end if;
    read(sint067,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t67 ",x);
    end if;
    read(sint068,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t68 ",x);
    end if;
    read(sint069,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t69 ",x);
    end if;
    read(sint070,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t70 ",x);
    end if;
    read(sint071,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t71 ",x);
    end if;
    read(sint072,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t72 ",x);
    end if;
    read(sint073,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t73 ",x);
    end if;
    read(sint074,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t74 ",x);
    end if;
    read(sint075,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t75 ",x);
    end if;
    read(sint076,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t76 ",x);
    end if;
    read(sint077,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t77 ",x);
    end if;
    read(sint078,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t78 ",x);
    end if;
    read(sint079,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t79 ",x);
    end if;
    read(sint080,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t80 ",x);
    end if;
    read(sint081,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t81 ",x);
    end if;
    read(sint082,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t82 ",x);
    end if;
    read(sint083,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t83 ",x);
    end if;
    read(sint084,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t84 ",x);
    end if;
    read(sint085,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t85 ",x);
    end if;
    read(sint086,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t86 ",x);
    end if;
    read(sint087,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t87 ",x);
    end if;
    read(sint088,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t88 ",x);
    end if;
    read(sint089,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t89 ",x);
    end if;
    read(sint090,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t90 ",x);
    end if;
    read(sint091,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t91 ",x);
    end if;
    read(sint092,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t92 ",x);
    end if;
    read(sint093,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t93 ",x);
    end if;
    read(sint094,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t94 ",x);
    end if;
    read(sint095,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t95 ",x);
    end if;
    read(sint096,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t96 ",x);
    end if;
    read(sint097,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t97 ",x);
    end if;
    read(sint098,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t98 ",x);
    end if;
    read(sint099,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t99 ",x);
    end if;
    read(sint100,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t100",x);
    end if;
    read(sint101,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t101",x);
    end if;
    read(sint102,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t102",x);
    end if;
    read(sint103,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t103",x);
    end if;
    read(sint104,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t104",x);
    end if;
    read(sint105,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t105",x);
    end if;
    read(sint106,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t106",x);
    end if;
    read(sint107,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t107",x);
    end if;
    read(sint108,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t108",x);
    end if;
    read(sint109,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t109",x);
    end if;
    read(sint110,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t110",x);
    end if;
    read(sint111,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t111",x);
    end if;
    read(sint112,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t112",x);
    end if;
    read(sint113,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t113",x);
    end if;
    read(sint114,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t114",x);
    end if;
    read(sint115,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t115",x);
    end if;
    read(sint116,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t116",x);
    end if;
    read(sint117,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t117",x);
    end if;
    read(sint118,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t118",x);
    end if;
    read(sint119,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t119",x);
    end if;
    read(sint120,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t120",x);
    end if;
    read(sint121,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t121",x);
    end if;
    read(sint122,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t122",x);
    end if;
    read(sint123,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t123",x);
    end if;
    read(sint124,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t124",x);
    end if;
    read(sint125,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t125",x);
    end if;
    read(sint126,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else  Write_Create("sin_t126",x);
    end if;
    read(sint127,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t127",x);
    end if;
    read(sint128,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t128",x);
    end if;
    read(sint129,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t129",x);
    end if;
    read(sint130,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t130",x);
    end if;
    read(sint131,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t131",x);
    end if;
    read(sint132,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t132",x);
    end if;
    read(sint133,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t133",x);
    end if;
    read(sint134,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t134",x);
    end if;
    read(sint135,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t135",x);
    end if;
    read(sint136,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t136",x);
    end if;
    read(sint137,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t137",x);
    end if;
    read(sint138,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t138",x);
    end if;
    read(sint139,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t139",x);
    end if;
    read(sint140,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t140",x);
    end if;
    read(sint141,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t141",x);
    end if;
    read(sint142,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t142",x);
    end if;
    read(sint143,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t143",x);
    end if;
    read(sint144,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t144",x);
    end if;
    read(sint145,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t145",x);
    end if;
    read(sint146,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t146",x);
    end if;
    read(sint147,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t147",x);
    end if;
    read(sint148,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t148",x);
    end if;
    read(sint149,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t149",x);
    end if;
    read(sint150,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t150",x);
    end if;
    read(sint151,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t151",x);
    end if;
    read(sint152,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t152",x);
    end if;
    read(sint153,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t153",x);
    end if;
    read(sint154,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t154",x);
    end if;
    read(sint155,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t155",x);
    end if;
    read(sint156,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t156",x);
    end if;
    read(sint157,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t157",x);
    end if;
    read(sint158,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t158",x);
    end if;
    read(sint159,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t159",x);
    end if;
    read(sint160,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t160",x);
    end if;
    read(sint161,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t161",x);
    end if;
    read(sint162,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t162",x);
    end if;
    read(sint163,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t163",x);
    end if;
    read(sint164,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t164",x);
    end if;
    read(sint165,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t165",x);
    end if;
    read(sint166,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t166",x);
    end if;
    read(sint167,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t167",x);
    end if;
    read(sint168,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t168",x);
    end if;
    read(sint169,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t169",x);
    end if;
    read(sint170,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t170",x);
    end if;
    read(sint171,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t171",x);
    end if;
    read(sint172,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t172",x);
    end if;
    read(sint173,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t173",x);
    end if;
    read(sint174,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t174",x);
    end if;
    read(sint175,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t175",x);
    end if;
    read(sint176,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t176",x);
    end if;
    read(sint177,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t177",x);
    end if;
    read(sint178,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t178",x);
    end if;
    read(sint179,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t179",x);
    end if;
    read(sint180,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t180",x);
    end if;
    read(sint181,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t181",x);
    end if;
    read(sint182,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t182",x);
    end if;
    read(sint183,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t183",x);
    end if;
    read(sint184,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t184",x);
    end if;
    read(sint185,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t185",x);
    end if;
    read(sint186,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t186",x);
    end if;
    read(sint187,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t187",x);
    end if;
    read(sint188,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t188",x);
    end if;
    read(sint189,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t189",x);
    end if;
    read(sint190,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else  Write_Create("sin_t190",x);
    end if;
    read(sint191,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t191",x);
    end if;
    read(sint192,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t192",x);
    end if;
    read(sint193,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t193",x);
    end if;
    read(sint194,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t194",x);
    end if;
    read(sint195,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t195",x);
    end if;
    read(sint196,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t196",x);
    end if;
    read(sint197,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t197",x);
    end if;
    read(sint198,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t198",x);
    end if;
    read(sint199,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t199",x);
    end if;
    read(sint200,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t200",x);
    end if;
    read(sint201,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t201",x);
    end if;
    read(sint202,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t202",x);
    end if;
    read(sint203,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t203",x);
    end if;
    read(sint204,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t204",x);
    end if;
    read(sint205,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t205",x);
    end if;
    read(sint206,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t206",x);
    end if;
    read(sint207,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t207",x);
    end if;
    read(sint208,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t208",x);
    end if;
    read(sint209,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t209",x);
    end if;
    read(sint210,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t210",x);
    end if;
    read(sint211,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t211",x);
    end if;
    read(sint212,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t212",x);
    end if;
    read(sint213,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t213",x);
    end if;
    read(sint214,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t214",x);
    end if;
    read(sint215,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t215",x);
    end if;
    read(sint216,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t216",x);
    end if;
    read(sint217,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t217",x);
    end if;
    read(sint218,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t218",x);
    end if;
    read(sint219,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t219",x);
    end if;
    read(sint220,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t220",x);
    end if;
    read(sint221,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t221",x);
    end if;
    read(sint222,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t222",x);
    end if;
    read(sint223,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t223",x);
    end if;
    read(sint224,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t224",x);
    end if;
    read(sint225,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t225",x);
    end if;
    read(sint226,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t226",x);
    end if;
    read(sint227,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t227",x);
    end if;
    read(sint228,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t228",x);
    end if;
    read(sint229,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t229",x);
    end if;
    read(sint230,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t230",x);
    end if;
    read(sint231,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t231",x);
    end if;
    read(sint232,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t232",x);
    end if;
    read(sint233,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t233",x);
    end if;
    read(sint234,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t234",x);
    end if;
    read(sint235,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t235",x);
    end if;
    read(sint236,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t236",x);
    end if;
    read(sint237,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t237",x);
    end if;
    read(sint238,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t238",x);
    end if;
    read(sint239,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t239",x);
    end if;
    read(sint240,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t240",x);
    end if;
    read(sint241,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t241",x);
    end if;
    read(sint242,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t242",x);
    end if;
    read(sint243,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t243",x);
    end if;
    read(sint244,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t244",x);
    end if;
    read(sint245,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t245",x);
    end if;
    read(sint246,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t246",x);
    end if;
    read(sint247,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t247",x);
    end if;
    read(sint248,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t248",x);
    end if;
    read(sint249,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t249",x);
    end if;
    read(sint250,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t250",x);
    end if;
    read(sint251,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t251",x);
    end if;
    read(sint252,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t252",x);
    end if;
    read(sint253,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t253",x);
    end if;
    read(sint254,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t254",x);
    end if;
    read(sint255,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("sin_t255",x);
    end if;
    read(cost000,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t0  ",x);
    end if;
    read(cost001,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t1  ",x);
    end if;
    read(cost002,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t2  ",x);
    end if;
    read(cost003,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t3  ",x);
    end if;
    read(cost004,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t4  ",x);
    end if;
    read(cost005,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t5  ",x);
    end if;
    read(cost006,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t6  ",x);
    end if;
    read(cost007,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t7  ",x);
    end if;
    read(cost008,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t8  ",x);
    end if;
    read(cost009,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t9  ",x);
    end if;
    read(cost010,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t10 ",x);
    end if;
    read(cost011,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t11 ",x);
    end if;
    read(cost012,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t12 ",x);
    end if;
    read(cost013,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t13 ",x);
    end if;
    read(cost014,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t14 ",x);
    end if;
    read(cost015,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t15 ",x);
    end if;
    read(cost016,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t16 ",x);
    end if;
    read(cost017,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t17 ",x);
    end if;
    read(cost018,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t18 ",x);
    end if;
    read(cost019,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t19 ",x);
    end if;
    read(cost020,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t20 ",x);
    end if;
    read(cost021,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t21 ",x);
    end if;
    read(cost022,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t22 ",x);
    end if;
    read(cost023,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t23 ",x);
    end if;
    read(cost024,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t24 ",x);
    end if;
    read(cost025,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t25 ",x);
    end if;
    read(cost026,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t26 ",x);
    end if;
    read(cost027,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t27 ",x);
    end if;
    read(cost028,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t28 ",x);
    end if;
    read(cost029,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t29 ",x);
    end if;
    read(cost030,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t30 ",x);
    end if;
    read(cost031,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t31 ",x);
    end if;
    read(cost032,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t32 ",x);
    end if;
    read(cost033,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t33 ",x);
    end if;
    read(cost034,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t34 ",x);
    end if;
    read(cost035,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t35 ",x);
    end if;
    read(cost036,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t36 ",x);
    end if;
    read(cost037,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t37 ",x);
    end if;
    read(cost038,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t38 ",x);
    end if;
    read(cost039,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t39 ",x);
    end if;
    read(cost040,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t40 ",x);
    end if;
    read(cost041,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t41 ",x);
    end if;
    read(cost042,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t42 ",x);
    end if;
    read(cost043,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t43 ",x);
    end if;
    read(cost044,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t44 ",x);
    end if;
    read(cost045,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t45 ",x);
    end if;
    read(cost046,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t46 ",x);
    end if;
    read(cost047,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t47 ",x);
    end if;
    read(cost048,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t48 ",x);
    end if;
    read(cost049,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t49 ",x);
    end if;
    read(cost050,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t50 ",x);
    end if;
    read(cost051,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t51 ",x);
    end if;
    read(cost052,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t52 ",x);
    end if;
    read(cost053,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t53 ",x);
    end if;
    read(cost054,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t54 ",x);
    end if;
    read(cost055,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t55 ",x);
    end if;
    read(cost056,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t56 ",x);
    end if;
    read(cost057,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t57 ",x);
    end if;
    read(cost058,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t58 ",x);
    end if;
    read(cost059,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t59 ",x);
    end if;
    read(cost060,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t60 ",x);
    end if;
    read(cost061,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t61 ",x);
    end if;
    read(cost062,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else  Write_Create("cos_t62 ",x);
    end if;
    read(cost063,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t63 ",x);
    end if;
    read(cost064,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t64 ",x);
    end if;
    read(cost065,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t65 ",x);
    end if;
    read(cost066,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t66 ",x);
    end if;
    read(cost067,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t67 ",x);
    end if;
    read(cost068,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t68 ",x);
    end if;
    read(cost069,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t69 ",x);
    end if;
    read(cost070,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t70 ",x);
    end if;
    read(cost071,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t71 ",x);
    end if;
    read(cost072,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t72 ",x);
    end if;
    read(cost073,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t73 ",x);
    end if;
    read(cost074,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t74 ",x);
    end if;
    read(cost075,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t75 ",x);
    end if;
    read(cost076,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t76 ",x);
    end if;
    read(cost077,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t77 ",x);
    end if;
    read(cost078,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t78 ",x);
    end if;
    read(cost079,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t79 ",x);
    end if;
    read(cost080,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t80 ",x);
    end if;
    read(cost081,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t81 ",x);
    end if;
    read(cost082,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t82 ",x);
    end if;
    read(cost083,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t83 ",x);
    end if;
    read(cost084,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t84 ",x);
    end if;
    read(cost085,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t85 ",x);
    end if;
    read(cost086,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t86 ",x);
    end if;
    read(cost087,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t87 ",x);
    end if;
    read(cost088,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t88 ",x);
    end if;
    read(cost089,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t89 ",x);
    end if;
    read(cost090,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t90 ",x);
    end if;
    read(cost091,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t91 ",x);
    end if;
    read(cost092,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t92 ",x);
    end if;
    read(cost093,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t93 ",x);
    end if;
    read(cost094,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t94 ",x);
    end if;
    read(cost095,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t95 ",x);
    end if;
    read(cost096,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t96 ",x);
    end if;
    read(cost097,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t97 ",x);
    end if;
    read(cost098,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t98 ",x);
    end if;
    read(cost099,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t99 ",x);
    end if;
    read(cost100,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t100",x);
    end if;
    read(cost101,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t101",x);
    end if;
    read(cost102,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t102",x);
    end if;
    read(cost103,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t103",x);
    end if;
    read(cost104,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t104",x);
    end if;
    read(cost105,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t105",x);
    end if;
    read(cost106,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t106",x);
    end if;
    read(cost107,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t107",x);
    end if;
    read(cost108,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t108",x);
    end if;
    read(cost109,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t109",x);
    end if;
    read(cost110,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t110",x);
    end if;
    read(cost111,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t111",x);
    end if;
    read(cost112,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t112",x);
    end if;
    read(cost113,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t113",x);
    end if;
    read(cost114,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t114",x);
    end if;
    read(cost115,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t115",x);
    end if;
    read(cost116,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t116",x);
    end if;
    read(cost117,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t117",x);
    end if;
    read(cost118,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t118",x);
    end if;
    read(cost119,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t119",x);
    end if;
    read(cost120,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t120",x);
    end if;
    read(cost121,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t121",x);
    end if;
    read(cost122,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t122",x);
    end if;
    read(cost123,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t123",x);
    end if;
    read(cost124,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t124",x);
    end if;
    read(cost125,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t125",x);
    end if;
    read(cost126,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else  Write_Create("cos_t126",x);
    end if;
    read(cost127,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t127",x);
    end if;
    read(cost128,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t128",x);
    end if;
    read(cost129,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t129",x);
    end if;
    read(cost130,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t130",x);
    end if;
    read(cost131,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t131",x);
    end if;
    read(cost132,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t132",x);
    end if;
    read(cost133,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t133",x);
    end if;
    read(cost134,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t134",x);
    end if;
    read(cost135,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t135",x);
    end if;
    read(cost136,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t136",x);
    end if;
    read(cost137,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t137",x);
    end if;
    read(cost138,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t138",x);
    end if;
    read(cost139,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t139",x);
    end if;
    read(cost140,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t140",x);
    end if;
    read(cost141,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t141",x);
    end if;
    read(cost142,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t142",x);
    end if;
    read(cost143,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t143",x);
    end if;
    read(cost144,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t144",x);
    end if;
    read(cost145,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t145",x);
    end if;
    read(cost146,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t146",x);
    end if;
    read(cost147,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t147",x);
    end if;
    read(cost148,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t148",x);
    end if;
    read(cost149,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t149",x);
    end if;
    read(cost150,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t150",x);
    end if;
    read(cost151,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t151",x);
    end if;
    read(cost152,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t152",x);
    end if;
    read(cost153,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t153",x);
    end if;
    read(cost154,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t154",x);
    end if;
    read(cost155,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t155",x);
    end if;
    read(cost156,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t156",x);
    end if;
    read(cost157,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t157",x);
    end if;
    read(cost158,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t158",x);
    end if;
    read(cost159,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t159",x);
    end if;
    read(cost160,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t160",x);
    end if;
    read(cost161,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t161",x);
    end if;
    read(cost162,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t162",x);
    end if;
    read(cost163,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t163",x);
    end if;
    read(cost164,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t164",x);
    end if;
    read(cost165,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t165",x);
    end if;
    read(cost166,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t166",x);
    end if;
    read(cost167,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t167",x);
    end if;
    read(cost168,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t168",x);
    end if;
    read(cost169,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t169",x);
    end if;
    read(cost170,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t170",x);
    end if;
    read(cost171,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t171",x);
    end if;
    read(cost172,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t172",x);
    end if;
    read(cost173,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t173",x);
    end if;
    read(cost174,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t174",x);
    end if;
    read(cost175,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t175",x);
    end if;
    read(cost176,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t176",x);
    end if;
    read(cost177,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t177",x);
    end if;
    read(cost178,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t178",x);
    end if;
    read(cost179,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t179",x);
    end if;
    read(cost180,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t180",x);
    end if;
    read(cost181,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t181",x);
    end if;
    read(cost182,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t182",x);
    end if;
    read(cost183,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t183",x);
    end if;
    read(cost184,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t184",x);
    end if;
    read(cost185,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t185",x);
    end if;
    read(cost186,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t186",x);
    end if;
    read(cost187,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t187",x);
    end if;
    read(cost188,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t188",x);
    end if;
    read(cost189,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t189",x);
    end if;
    read(cost190,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else  Write_Create("cos_t190",x);
    end if;
    read(cost191,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t191",x);
    end if;
    read(cost192,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t192",x);
    end if;
    read(cost193,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t193",x);
    end if;
    read(cost194,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t194",x);
    end if;
    read(cost195,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t195",x);
    end if;
    read(cost196,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t196",x);
    end if;
    read(cost197,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t197",x);
    end if;
    read(cost198,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t198",x);
    end if;
    read(cost199,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t199",x);
    end if;
    read(cost200,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t200",x);
    end if;
    read(cost201,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t201",x);
    end if;
    read(cost202,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t202",x);
    end if;
    read(cost203,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t203",x);
    end if;
    read(cost204,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t204",x);
    end if;
    read(cost205,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t205",x);
    end if;
    read(cost206,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t206",x);
    end if;
    read(cost207,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t207",x);
    end if;
    read(cost208,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t208",x);
    end if;
    read(cost209,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t209",x);
    end if;
    read(cost210,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t210",x);
    end if;
    read(cost211,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t211",x);
    end if;
    read(cost212,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t212",x);
    end if;
    read(cost213,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t213",x);
    end if;
    read(cost214,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t214",x);
    end if;
    read(cost215,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t215",x);
    end if;
    read(cost216,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t216",x);
    end if;
    read(cost217,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t217",x);
    end if;
    read(cost218,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t218",x);
    end if;
    read(cost219,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t219",x);
    end if;
    read(cost220,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t220",x);
    end if;
    read(cost221,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t221",x);
    end if;
    read(cost222,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t222",x);
    end if;
    read(cost223,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t223",x);
    end if;
    read(cost224,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t224",x);
    end if;
    read(cost225,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t225",x);
    end if;
    read(cost226,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t226",x);
    end if;
    read(cost227,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t227",x);
    end if;
    read(cost228,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t228",x);
    end if;
    read(cost229,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t229",x);
    end if;
    read(cost230,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t230",x);
    end if;
    read(cost231,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t231",x);
    end if;
    read(cost232,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t232",x);
    end if;
    read(cost233,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t233",x);
    end if;
    read(cost234,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t234",x);
    end if;
    read(cost235,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t235",x);
    end if;
    read(cost236,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t236",x);
    end if;
    read(cost237,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t237",x);
    end if;
    read(cost238,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t238",x);
    end if;
    read(cost239,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t239",x);
    end if;
    read(cost240,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t240",x);
    end if;
    read(cost241,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t241",x);
    end if;
    read(cost242,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t242",x);
    end if;
    read(cost243,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t243",x);
    end if;
    read(cost244,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t244",x);
    end if;
    read(cost245,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t245",x);
    end if;
    read(cost246,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t246",x);
    end if;
    read(cost247,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t247",x);
    end if;
    read(cost248,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t248",x);
    end if;
    read(cost249,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t249",x);
    end if;
    read(cost250,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t250",x);
    end if;
    read(cost251,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t251",x);
    end if;
    read(cost252,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t252",x);
    end if;
    read(cost253,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t253",x);
    end if;
    read(cost254,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t254",x);
    end if;
    read(cost255,x,fail);
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("cos_t255",x);
    end if;
  end Sine_Cosine_Table;

  procedure Log_exp_of_Pi is

    da_pi : constant deca_double
          := create( 3.14159265358979312E+00, 1.22464679914735321E-16,
                    -2.99476980971833967E-33, 1.11245422086336528E-49,
                     5.67223197964031574E-66, 1.74498621613524860E-83,
                     6.02937273224953984E-100, 1.91012354687998999E-116,
                     3.04397816534429330E-133, -4.71430003123113544E-150);
    exp_of_pi,log_of_exp_of_pi : deca_double;
    log_of_pi,exp_of_log_of_pi : deca_double;
    ans : character;
    x,expx,logx,logexpx,explogx : deca_double;

  begin
    new_line;
    put_line("testing log(exp(pi)) and exp(log(pi)) ...");
    exp_of_pi := exp(da_pi);
    put("     exp(pi) : "); put(exp_of_pi); new_line;
    log_of_exp_of_pi := log(exp_of_pi);
    put("log(exp(pi)) : "); put(log_of_exp_of_pi); new_line;
    put("        pi   : "); put(da_pi); new_line;
    put("  difference : "); put(log_of_exp_of_pi - da_pi,3); new_line;
    log_of_pi := log(da_pi);
    put("     log(pi) : "); put(log_of_pi); new_line;
    exp_of_log_of_pi := exp(log_of_pi);
    put("exp(log(pi)) : "); put(exp_of_log_of_pi); new_line;
    put("        pi   : "); put(da_pi); new_line;
    put("  difference : "); put(exp_of_log_of_pi - da_pi,3); new_line;
    loop
      put("Test a random number ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      x := abs(random);
      expx := exp(x); logexpx := log(expx);
      put("log(exp(x)) : "); put(logexpx); new_line;
      put("        x   : "); put(x); new_line;
      put(" difference : "); put(logexpx - x,3); new_line;
      logx := log(x); explogx := exp(logx);
      put("exp(log(x)) : "); put(explogx); new_line;
      put("        x   : "); put(x); new_line;
      put(" difference : "); put(explogx - x,3); new_line;
    end loop;
  end Log_exp_of_Pi;

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
    put_line("  9. test the value of da_eps");
    put_line("  A. write 10 leading doubles for pi and multiples");
    put_line("  B. write 10 leading doubles for sine and cosine table");
    put_line("  C. test log(exp(pi)) = pi = exp(log(pi))");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, or C to select a test : ");
    Ask_Alternative(ans,"123456789ABC");
    case ans is
      when '1' => Test_Addition_and_Subtraction;
      when '2' => Test_Multiplication_and_Division;
      when '3' => Test_Read;
      when '4' => Test_Write;
      when '5' => Log10log2exp1_doubles;
      when '6' => inverse_factorials;
      when '7' => Test_io;
      when '8' => Test_sqrt2;
      when '9' => Test_da_eps;
      when 'A' => Write_Pi;
      when 'B' => Sine_Cosine_Table;
      when 'C' => Log_exp_of_Pi;
      when others => null;
    end case;
  end Main;

end Test_Deca_Doubles;
