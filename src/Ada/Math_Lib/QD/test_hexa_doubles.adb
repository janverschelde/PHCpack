with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Hexa_Double_Numbers_io;             use Hexa_Double_Numbers_io;
with Hexa_Double_Constants;

package body Test_Hexa_Doubles is

  function random return hexa_double is

    res : hexa_double;
    first : constant double_float := Standard_Random_Numbers.Random; 
    second : double_float := Standard_Random_Numbers.Random; 
    eps : constant double_float := 2.0**(-52);
    multiplier : double_float := eps;

  begin
    res := create(first);
    res := res + eps*second;
    for k in 3..16 loop
      multiplier := eps*multiplier;
      second := Standard_Random_Numbers.Random;
      res := res + multiplier*second;
    end loop;
    return res;
  end random;

  procedure Write ( x : in hexa_double ) is
  begin
    put("  hihihihi : "); put(hihihihi_part(x),2,17,3); new_line;
    put("  lohihihi : "); put(lohihihi_part(x),2,17,3); new_line;
    put("  hilohihi : "); put(hilohihi_part(x),2,17,3); new_line;
    put("  lolohihi : "); put(lolohihi_part(x),2,17,3); new_line;
    put("  hihilohi : "); put(hihilohi_part(x),2,17,3); new_line;
    put("  lohilohi : "); put(lohilohi_part(x),2,17,3); new_line;
    put("  hilolohi : "); put(hilolohi_part(x),2,17,3); new_line;
    put("  lololohi : "); put(lololohi_part(x),2,17,3); new_line;
    put("  hihihilo : "); put(hihihilo_part(x),2,17,3); new_line;
    put("  lohihilo : "); put(lohihilo_part(x),2,17,3); new_line;
    put("  hilohilo : "); put(hilohilo_part(x),2,17,3); new_line;
    put("  lolohilo : "); put(lolohilo_part(x),2,17,3); new_line;
    put("  hihilolo : "); put(hihilolo_part(x),2,17,3); new_line;
    put("  lohilolo : "); put(lohilolo_part(x),2,17,3); new_line;
    put("  hilololo : "); put(hilololo_part(x),2,17,3); new_line;
    put("  lolololo : "); put(lolololo_part(x),2,17,3); new_line;
  end Write;

  procedure Write_Create ( s : in string; x : in hexa_double ) is

  -- DESCRIPTION :
  --   Writes the statement of the constant definition of x,
  --   as the variable s.

  begin
    put("  "); put(s); put_line(" : constant hexa_double");
    put("           := create(");
    put(hihihihi_part(x),2,17,3); put(",");
    put(lohihihi_part(x),2,17,3); put_line(",");
    put("                     ");
    put(hilohihi_part(x),2,17,3); put(",");
    put(lolohihi_part(x),2,17,3); put_line(",");
    put("                     ");
    put(hihilohi_part(x),2,17,3); put(",");
    put(lohilohi_part(x),2,17,3); put_line(",");
    put("                     ");
    put(hilolohi_part(x),2,17,3); put(",");
    put(lololohi_part(x),2,17,3); put_line(",");
    put("                     ");
    put(hihihilo_part(x),2,17,3); put(",");
    put(lohihilo_part(x),2,17,3); put_line(",");
    put("                     ");
    put(hilohilo_part(x),2,17,3); put(",");
    put(lolohilo_part(x),2,17,3); put_line(",");
    put("                     ");
    put(hihilolo_part(x),2,17,3); put(",");
    put(lohilolo_part(x),2,17,3); put_line(",");
    put("                     ");
    put(hilololo_part(x),2,17,3); put(",");
    put(lolololo_part(x),2,17,3); put_line(");");
  end Write_Create;

  procedure Test_Add_and_Subtract is

    x : constant hexa_double := random;
    y : constant hexa_double := random;
    z : constant hexa_double := x + y;
    v : constant hexa_double := z - y;
    w : constant double_float := Standard_Random_Numbers.Random; 
    r : constant hexa_double := w + x;
    s : constant hexa_double := r - w;

  begin
   put_line("All parts of a random hexa double x :"); Write(x);
   put_line("All parts of a random hexa double y :"); Write(y);
   put_line("All parts of x + y :"); Write(z);
   put_line("All parts of (x + y) - y :"); Write(v);
   put("A random double w : "); put(w,2,17,3); new_line;
   put_line("All parts of w + x :"); Write(r);
   put_line("All parts of (w + x) - w :"); Write(s);
  end Test_Add_and_Subtract;

  procedure Test_Multiplication_and_Division is

    x : constant hexa_double := random;
    y : constant hexa_double := random;
    z : constant hexa_double := x*y;
    v : constant hexa_double := z/y;
    w : constant double_float := Standard_Random_Numbers.Random; 
    r : constant hexa_double := w*x;
    s : constant hexa_double := r/w;

  begin
   put_line("All parts of a random hexa double x :"); Write(x);
   put_line("All parts of a random hexa double y :"); Write(y);
   put_line("All parts of x * y :"); Write(z);
   put_line("All parts of (x * y) / y :"); Write(v);
   put("A random double w : "); put(w,2,17,3); new_line;
   put_line("All parts of w * x :"); Write(r);
   put_line("All parts of (w * x) / w :"); Write(s);
  end Test_Multiplication_and_Division;

  procedure Test_Read is

  --   >>> from sympy import evalf, sqrt
  --   >>> s2 = sqrt(2).evalf(256)
  --   >>> s2
  --   1.414213562373095048801688724209698078569671875376948073176679737
  --   99073247846210703885038753432764157273501384623091229702492483605
  --   58507372126441214970999358314132226659275055927557999505011527820
  --   60571470109559971605970274534596862014728517418640889198609552
    
    sqrt2 : constant string
      := "1.414213562373095048801688724209698078569671875376948073176679737"
       & "99073247846210703885038753432764157273501384623091229702492483605"
       & "58507372126441214970999358314132226659275055927557999505011527820"
       & "60571470109559971605970274534596862014728517418640889198609552";

    x,r : hexa_double;
    two : constant hexa_double := create(2.0);
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

  procedure Log10log2exp1_doubles is

  -- >>> from sympy import log, evalf
  -- >>> x = log(10)
  -- >>> x.evalf(256)
  -- 2.3025850929940456840179914546843642076011014886287729760333279009
  -- 675726096773524802359972050895982983419677840422862486334095254650
  -- 828067566662873690987816894829072083255546808437998948262331985283
  -- 93505308965377732628846163366222287698219886746543667474404
  -- >>> from sympy import exp, evalf
  -- >>> x = exp(1)
  -- >>> x.evalf(256)
  -- 2.7182818284590452353602874713526624977572470936999595749669676277
  -- 240766303535475945713821785251664274274663919320030599218174135966
  -- 290435729003342952605956307381323286279434907632338298807531952510
  -- 19011573834187930702154089149934884167509244761460668082265
  -- >>> from sympy import log, evalf
  -- >>> x = log(2)
  -- >>> x.evalf(256)
  -- 0.6931471805599453094172321214581765680755001343602552541206800094
  -- 933936219696947156058633269964186875420014810205706857336855202357
  -- 581305570326707516350759619307275708283714351903070386238916734711
  -- 233501153644979552391204751726815749320651555247341395258830

    log10 : constant string
      := "2.3025850929940456840179914546843642076011014886287729760333279009"
       & "675726096773524802359972050895982983419677840422862486334095254650"
       & "828067566662873690987816894829072083255546808437998948262331985283"
       & "93505308965377732628846163366222287698219886746543667474404";
    log2 : constant string
      := "0.6931471805599453094172321214581765680755001343602552541206800094"
       & "933936219696947156058633269964186875420014810205706857336855202357"
       & "581305570326707516350759619307275708283714351903070386238916734711"
       & "233501153644979552391204751726815749320651555247341395258830";
    exp1 : constant string
      := "2.7182818284590452353602874713526624977572470936999595749669676277"
       & "240766303535475945713821785251664274274663919320030599218174135966"
       & "290435729003342952605956307381323286279434907632338298807531952510"
       & "19011573834187930702154089149934884167509244761460668082265";

    x : hexa_double;
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
  -- >>> f = [1/factorial(k) for k in range(3,41)]
  -- >>> for x in f: print(x.evalf(256))

    f0 : constant string
       := "0.16666666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666666666666666666"
        & "66666666666666666666666666666666666667";
    f1 : constant string
       := "0.04166666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666666666666666666"
        & "6666666666666666666666666666666666666666666666666666666"
        & "666666666666666666666666666666666666667";
    f2 : constant string
       := "0.00833333333333333333333333333333333333333333333333333"
        & "3333333333333333333333333333333333333333333333333333333"
        & "3333333333333333333333333333333333333333333333333333333"
        & "3333333333333333333333333333333333333333333333333333333"
        & "3333333333333333333333333333333333333333";
    f3 : constant string
       := "0.00138888888888888888888888888888888888888888888888888"
        & "8888888888888888888888888888888888888888888888888888888"
        & "8888888888888888888888888888888888888888888888888888888"
        & "8888888888888888888888888888888888888888888888888888888"
        & "8888888888888888888888888888888888888889";
    f4 : constant string
       := "0.00019841269841269841269841269841269841269841269841269"
        & "8412698412698412698412698412698412698412698412698412698"
        & "4126984126984126984126984126984126984126984126984126984"
        & "1269841269841269841269841269841269841269841269841269841"
        & "26984126984126984126984126984126984126984";
    f5 : constant string
       := "0.00002480158730158730158730158730158730158730158730158"
        & "7301587301587301587301587301587301587301587301587301587"
        & "3015873015873015873015873015873015873015873015873015873"
        & "0158730158730158730158730158730158730158730158730158730"
        & "158730158730158730158730158730158730158730";
    f6 : constant string
       := "0.00000275573192239858906525573192239858906525573192239"
        & "8589065255731922398589065255731922398589065255731922398"
        & "5890652557319223985890652557319223985890652557319223985"
        & "8906525573192239858906525573192239858906525573192239858"
        & "9065255731922398589065255731922398589065256";
    f7 : constant string
       := "0.00000027557319223985890652557319223985890652557319223"
        & "9858906525573192239858906525573192239858906525573192239"
        & "8589065255731922398589065255731922398589065255731922398"
        & "5890652557319223985890652557319223985890652557319223985"
        & "89065255731922398589065255731922398589065256";
    f8 : constant string
       := "0.00000002505210838544171877505210838544171877505210838"
        & "5441718775052108385441718775052108385441718775052108385"
        & "4417187750521083854417187750521083854417187750521083854"
        & "4171877505210838544171877505210838544171877505210838544"
        & "171877505210838544171877505210838544171877505";
    f9 : constant string
       := "0.00000000208767569878680989792100903212014323125434236"
        & "5453476564587675698786809897921009032120143231254342365"
        & "4534765645876756987868098979210090321201432312543423654"
        & "5347656458767569878680989792100903212014323125434236545"
        & "3476564587675698786809897921009032120143231254";
    f10 : constant string
       := "0.00000000016059043836821614599392377170154947932725710"
        & "5034882812660590438368216145993923771701549479327257105"
        & "0348828126605904383682161459939237717015494793272571050"
        & "3488281266059043836821614599392377170154947932725710503"
        & "48828126605904383682161459939237717015494793273";
    f11 : constant string
        := "0.00000000001147074559772972471385169797868210566623265"
         & "0359634486618613602740586867570994555121539248523375507"
         & "5024916294757564598834440104281374122643963913805183646"
         & "4534877233289931702630115328528026940725353423766122178"
         & "820591519004217416915829614242312655011067709480";
    f12 : constant string
        := "0.00000000000076471637318198164759011319857880704441551"
         & "0023975632441240906849372457838066303674769283234891700"
         & "5001661086317170973255629340285424941509597594253678909"
         & "7635658482219328780175341021901868462715023561584408145"
         & "25470610126694782779438864094948751033407118063203";
    f13 : constant string
        := "0.00000000000004779477332387385297438207491117544027596"
         & "9376498477027577556678085778614879143979673080202180731"
         & "2812603817894823185828476833767839058844349849640854931"
         & "8602228655138708048760958813868866778919688972599025509"
         & "078419131329184239237149290059342969395879448789502";
    f14 : constant string
        := "0.00000000000000281145725434552076319894558301032001623"
         & "3492735204531033973922240339918522302587039592953069454"
         & "7812506106934989599166380990221637591696726461743579701"
         & "8741307567949335767574174047874639222289393468976413265"
         & "2399070077252461317198323111799613511409340852229119";
    f15 : constant string
        := "0.00000000000000015619206968586226462216364350057333423"
         & "5194040844696168554106791129995473461254835532941837191"
         & "9322917005940832755509243388345646532872040358985754427"
         & "8818961531552740875976343002659702179016077414943134070"
         & "29110594487362478509554623950999785284116300473460621";
    f16 : constant string
        := "0.00000000000000000822063524662432971695598123687228074"
         & "9220738991826114134426673217368182813750254501733780904"
         & "8385416684523201723974170704649770870151160018893987075"
         & "1516787449029091625051386473824194851527161969207533372"
         & "1205845234144013044787129599742104133074296318281371692";
    f17 : constant string
        := "0.00000000000000000041103176233121648584779906184361403"
         & "7461036949591305706721333660868409140687512725086689045"
         & "2419270834226160086198708535232488543507558000944699353"
         & "7575839372451454581252569323691209742576358098460376668"
         & "60602922617072006522393564799871052066537148159140685846";
    f18 : constant string
        := "0.00000000000000000001957294106339126123084757437350543"
         & "0355287473790062176510539698136590911461310129766032811"
         & "6781870039725055242199938501677737549690836095283080921"
         & "6075039970116735932440598539223390940122683718974303650"
         & "886001391722415241201139792761843358126922451504352707546";
    f19 : constant string
        := "0.00000000000000000000088967913924505732867488974425024"
         & "6834331248808639189841388168097117768702786824080274218"
         & "7126448638169320692827269931894442615895038004331049132"
         & "8003410907732578906020027206328335951823758350862468347"
         & "76754551780556432914550635421644742536940556597747057761571";
    f20 : constant string
        := "0.00000000000000000000003868170170630684037716911931522"
         & "8123231793426462573471364702960744250813164644525229313"
         & "8570715158181274812731620431821497505038914695840480397"
         & "0782756995988372995913914226362101563122772102211411667"
         & "294241109469807144745456798009410757624756763738150894678944";
    f21 : constant string
        := "0.00000000000000000000000161173757109611834904871330480"
         & "1171801324726102607227973529290031010450548526855217888"
         & "0773779798257553117197150851325895729376621445660020016"
         & "5449281541499515541496413092765087565130115504258808819"
         & "4705933795612419643643940332503921149010315318224229539449560";
    f22 : constant string
        := "0.00000000000000000000000006446950284384473396194853219"
         & "2046872052989044104289118941171601240418021941074208715"
         & "5230951191930302124687886034053035829175064857826400800"
         & "6617971261659980621659856523710603502605204620170352352"
         & "778823735182449678574575761330015684596041261272896918157798240";
    f23 : constant string
        := "0.00000000000000000000000000247959626322479746007494354"
         & "5847956617422655542472658420814292355400693151579777258"
         & "2893498122766550081718764847463578301122117879147169261"
         & "5639152740833076177756148327835023211638661716160398167"
         & "4145701436608634491759452215896159878690785100489575737752999323";
    f24 : constant string
        := "0.00000000000000000000000000009183689863795546148425716"
         & "8364739133978616871943431793363492309459284931539991750"
         & "3070129560102464817841435735091243640782300662190635898"
         & "5764413064475299117694672160290186044875505989487422154"
         & "3486877830985504981176276007996154069581140188907021323"
         & "62048145641";
    f25 : constant string
        := "0.00000000000000000000000000000327988923706983791015204"
         & "1727312111927807745426551135477267582480688747554999705"
         & "3681076055717945172065765561967544415742225023649665567"
         & "8063014752302689254203381148581792358745553785338836505"
         & "5124531351106625177899152714571291216770755006746679332"
         & "986445766300";
    f26 : constant string
        := "0.00000000000000000000000000000011309962886447716931558"
         & "7645769383169924405014708659844043709740713405088103438"
         & "1161416415714411902485026398688536014335938793918953985"
         & "0967690163872506526007013143054544564094674268459959879"
         & "5004294184520918109582729403950734179888646724370575149"
         & "4133257160793";
    f27 : constant string
        := "0.00000000000000000000000000000000376998762881590564385"
         & "2921525646105664146833823621994801456991357113502936781"
         & "2705380547190480396749500879956284533811197959797298466"
         & "1698923005462416884200233771435151485469822475615331995"
         & "9833476472817363936986090980131691139329621557479019171"
         & "647110857202644";
    f28 : constant string
        := "0.00000000000000000000000000000000012161250415535179496"
         & "2997468569229214972478510439419187143773914745596868928"
         & "4280818727328725174088693576772783372058425740638622531"
         & "1667707193724594093038717218433391983402252337923075225"
         & "6768821821703785933451164225165538423849342630886419973"
         & "2789390599097627";
    f29 : constant string
        := "0.00000000000000000000000000000000000380039075485474359"
         & "2593670892788412967889953451231849598242934835799902154"
         & "0133775585229022661690271674274149480376825804394956954"
         & "0989615849803893565407459913076043499481320385560096100"
         & "8024025681928243310420348882036423075745291957215200624"
         & "164966845622180085";
    f30 : constant string
        := "0.00000000000000000000000000000000000011516335620771950"
         & "2805868814932982211148180407613086351461907116236360671"
         & "3337387138946334020051220353765883317587176539527119907"
         & "6999685328781936168648710906456849803014585466229093821"
         & "2364364414603886160921828754001103729568039150218642443"
         & "1565141468370357601";
    f31 : constant string
        := "0.00000000000000000000000000000000000000338715753552116"
         & "1847231435733323006210240600223914304454761974006951784"
         & "4509923151145480412354447657463702450517269898221385879"
         & "6382343686140645181430844438425201464794546631359679230"
         & "0363657776900114298850642022176503050869648210300548307"
         & "151662180789324581181";
    f32 : constant string
        := "0.00000000000000000000000000000000000000009677592958631"
         & "8909920898163809228748864017149254694412993199257341479"
         & "5557426375747013726067269933070391498586207711377753882"
         & "2753781248175447005183738412526434327565558475181705120"
         & "8581818793625717551395732629205042944310561377437158523"
         & "06147606230826641660517";
    f33 : constant string
        := "0.00000000000000000000000000000000000000000268822026628"
         & "6636386691615661367465246222698590408178138699979370596"
         & "6543261843770750381279646387029733097182950214204937607"
         & "8409827256893762416810659400347956509099043290977269586"
         & "6905050522045158820872103684144584526230848927151032181"
         & "196152112841896289350143";
    f34 : constant string
        := "0.00000000000000000000000000000000000000000007265460179"
         & "1530713153827450307228790438451313254275084829729172178"
         & "2879547617399209469764314767217019813437377032816349665"
         & "0767833169105236822075963767576971797543217386242628907"
         & "7483920284379598887050597396868772554762995916950027896"
         & "78908519223897016998243631";
    f35 : constant string
        := "0.00000000000000000000000000000000000000000000191196320"
         & "5040281925100722376506020801011876664586186442887609794"
         & "1654724937299979196572745125453079468774667816653061833"
         & "2914942978134348337423051678094130836777453089111648129"
         & "1512734744325778918080278878864967698809552524130263892"
         & "020765399795762372894274640";
    f36 : constant string
        := "0.00000000000000000000000000000000000000000000004902469"
         & "7565135433976941599397590276949022478579132985715066917"
         & "7991146793264102030681352439114181524840376097862899021"
         & "3664485717388060213780078248156259765045575720233632003"
         & "3115711147290404587643084073817050453815629551900775997"
         & "23130167691784006084344293948";
    f37 : constant string
        := "0.00000000000000000000000000000000000000000000000122561"
         & "7439128385849423539984939756923725561964478324642876672"
         & "9449778669831602550767033810977854538121009402446572475"
         & "5341612142934701505344501956203906494126139393005840800"
         & "0827892778682260114691077101845426261345390738797519399"
         & "930782541922946001521086073487";

    x : hexa_double;
    fail : boolean;

  begin
    read(f0,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac0 ",x);
    end if;
    read(f1,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac1 ",x);
    end if;
    read(f2,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac2 ",x);
    end if;
    read(f3,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac3 ",x);
    end if;
    read(f4,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac4 ", x);
    end if;
    read(f5,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac5 ", x);
    end if;
    read(f6,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac6 ", x);
    end if;
    read(f7,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac7 ", x);
    end if;
    read(f8,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac8 ", x);
    end if;
    read(f9,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac9 ", x);
    end if;
    read(f10,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac10", x);
    end if;
    read(f11,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac11", x);
    end if;
    read(f12,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac12", x);
    end if;
    read(f13,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac13", x);
    end if;
    read(f14,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac14", x);
    end if;
    read(f15,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac15", x);
    end if;
    read(f16,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac16", x);
    end if;
    read(f17,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac17", x);
    end if;
    read(f18,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac18", x);
    end if;
    read(f19,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac19", x);
    end if;
    read(f20,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac20", x);
    end if;
    read(f21,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac21", x);
    end if;
    read(f22,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac22", x);
    end if;
    read(f23,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac23", x);
    end if;
    read(f24,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac24", x);
    end if;
    read(f25,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac25", x);
    end if;
    read(f26,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac26", x);
    end if;
    read(f27,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac27", x);
    end if;
    read(f28,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac28", x);
    end if;
    read(f29,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac29", x);
    end if;
    read(f30,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac30", x);
    end if;
    read(f31,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac31", x);
    end if;
    read(f32,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac32", x);
    end if;
    read(f33,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac33", x);
    end if;
    read(f34,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac34", x);
    end if;
    read(f35,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac35", x);
    end if;
    read(f36,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac36", x);
    end if;
    read(f37,x,fail);
    new_line;
    if fail
     then put_line("The read procedure reports failure!");
     else Write_Create("i_fac37", x);
    end if;
  end inverse_factorials;

  procedure Test_io is

    x,y : hexa_double;
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

    n,x,y,z,e,a : hexa_double;
    max_steps : constant natural32 := 10;
    sqrt2 : constant string
      := "1.414213562373095048801688724209698078569671875376948073176679737"
       & "99073247846210703885038753432764157273501384623091229702492483605"
       & "58507372126441214970999358314132226659275055927557999505011527820"
       & "60571470109559971605970274534596862014728517418640889198609552";
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

  procedure Test_hd_eps is

    one : constant hexa_double := create(1.0);
    eps : constant double_float := Hexa_Double_Constants.hd_eps;
    one_plus_hd_eps : constant hexa_double := one + eps;
    inc : constant double_float := (eps/2.0);
    one_plus_hd_eps_half : constant hexa_double := one + inc;

  begin
    new_line;
    put("    hd_eps   :"); put(eps); new_line;
    put_line("1 + hd_eps   : "); put(one_plus_hd_eps,255); new_line;
    put_line("1 + hd_eps/2 : "); put(one_plus_hd_eps_half,255); new_line;
  end Test_hd_eps;

  procedure Log_exp_of_Pi is

    hd_pi : constant hexa_double
          := create( 3.14159265358979312E+00, 1.22464679914735321E-16,
                    -2.99476980971833967E-33, 1.11245422086336528E-49,
                     5.67223197964031574E-66, 1.74498621613524860E-83,
                     6.02937273224953984E-100, 1.91012354687998999E-116,
                     3.04397816534429330E-133, -4.71430003094702878E-150,
                     1.00154916943553891E-166, 6.21040447841589472E-183,
                    -1.71681323916116030E-199, -5.49577495896909932E-216,
                    -1.74900899480240872E-232, -3.59150954550134434E-249);

    exp_of_pi,log_of_exp_of_pi : hexa_double;
    log_of_pi,exp_of_log_of_pi : hexa_double;
    ans : character;
    x,expx,logx,logexpx,explogx : hexa_double;

  begin
    new_line;
    put_line("testing log(exp(pi)) and exp(log(pi)) ...");
    exp_of_pi := exp(hd_pi);
    put("     exp(pi) : "); put(exp_of_pi); new_line;
    log_of_exp_of_pi := log(exp_of_pi);
    put("log(exp(pi)) : "); put(log_of_exp_of_pi); new_line;
    put("        pi   : "); put(hd_pi); new_line;
    put("  difference : "); put(log_of_exp_of_pi - hd_pi,3); new_line;
    log_of_pi := log(hd_pi);
    put("     log(pi) : "); put(log_of_pi); new_line;
    exp_of_log_of_pi := exp(log_of_pi);
    put("exp(log(pi)) : "); put(exp_of_log_of_pi); new_line;
    put("        pi   : "); put(hd_pi); new_line;
    put("  difference : "); put(exp_of_log_of_pi - hd_pi,3); new_line;
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

  procedure Write_Pi is

  -- To evaluate the constants with 160 decimal places, sympy is used.
  -- >>> from sympy import pi
  -- >>> from sympy import evalf
  -- >>> pi.evalf(256)
  -- 3.141592653589793238462643383279502884197169399375105820
  -- 97494459230781640628620899862803482534211706798214808651
  -- 32823066470938446095505822317253594081284811174502841027
  -- 01938521105559644622948954930381964428810975665933446128
  -- 475648233786783165271201909145649
    pi : constant string
       := "3.141592653589793238462643383279502884197169399375105820"
        & "97494459230781640628620899862803482534211706798214808651"
        & "32823066470938446095505822317253594081284811174502841027"
        & "01938521105559644622948954930381964428810975665933446128"
        & "475648233786783165271201909145649";
  -- >>> twopi = 2*pi
  -- >>> twopi.evalf(256)
  -- 6.283185307179586476925286766559005768394338798750211641
  -- 94988918461563281257241799725606965068423413596429617302
  -- 65646132941876892191011644634507188162569622349005682054
  -- 03877042211119289245897909860763928857621951331866892256
  -- 951296467573566330542403818291297
    twopi : constant string
          := "6.283185307179586476925286766559005768394338798750211641"
           & "94988918461563281257241799725606965068423413596429617302"
           & "65646132941876892191011644634507188162569622349005682054"
           & "03877042211119289245897909860763928857621951331866892256"
           & "951296467573566330542403818291297";
  -- >>> pi2 = pi/2
  -- >>> pi2.evalf(256)
  -- 1.570796326794896619231321691639751442098584699687552910
  -- 48747229615390820314310449931401741267105853399107404325
  -- 66411533235469223047752911158626797040642405587251420513
  -- 50969260552779822311474477465190982214405487832966723064
  -- 237824116893391582635600954572824
    pi2 : constant string
        := "1.570796326794896619231321691639751442098584699687552910"
         & "48747229615390820314310449931401741267105853399107404325"
         & "66411533235469223047752911158626797040642405587251420513"
         & "50969260552779822311474477465190982214405487832966723064"
         & "237824116893391582635600954572824";
  -- >>> pi4 = pi/4
  -- >>> pi4.evalf(256)
  -- 0.785398163397448309615660845819875721049292349843776455
  -- 24373614807695410157155224965700870633552926699553702162
  -- 83205766617734611523876455579313398520321202793625710256
  -- 75484630276389911155737238732595491107202743916483361532
  -- 1189120584466957913178004772864121
    pi4 : constant string
        := "0.785398163397448309615660845819875721049292349843776455"
         & "24373614807695410157155224965700870633552926699553702162"
         & "83205766617734611523876455579313398520321202793625710256"
         & "75484630276389911155737238732595491107202743916483361532"
         & "1189120584466957913178004772864121";
  -- >>> threepi4 = 3*pi/4
  -- >>> threepi4.evalf(256)
  -- 2.356194490192344928846982537459627163147877049531329365
  -- 73120844423086230471465674897102611900658780098661106488
  -- 49617299853203834571629366737940195560963608380877130770
  -- 26453890829169733467211716197786473321608231749450084596
  -- 356736175340087373953401431859236
    threepi4 : constant string
             := "2.356194490192344928846982537459627163147877049531329365"
              & "73120844423086230471465674897102611900658780098661106488"
              & "49617299853203834571629366737940195560963608380877130770"
              & "26453890829169733467211716197786473321608231749450084596"
              & "356736175340087373953401431859236";
  -- >>> pi1024 = pi/1024
  -- >>> pi1024.evalf(256)
  -- 0.003067961575771282459436175178983889535348798241577251
  -- 77829584432842560195926387597522269025912316119920131649
  -- 07356272525850525826265142404606692962970004698412600430
  -- 69044861837017148090452098588799201137137510718423763130
  -- 984839500228307405434835158114400047
    pi1024 : constant string
           := "0.003067961575771282459436175178983889535348798241577251"
            & "77829584432842560195926387597522269025912316119920131649"
            & "07356272525850525826265142404606692962970004698412600430"
            & "69044861837017148090452098588799201137137510718423763130"
            & "984839500228307405434835158114400047";

    x : hexa_double;
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

  procedure Sine_Cosine_Table is

  -- A Python script with the evalf(160) of sympy was applied
  -- to get the constant strings as below.

    sint000 : constant string
            := "0.00306795676296597627014536549091984251894461021345"
             & "1995397146895898975797563020651153816092963359474881"
             & "4215830594626216977433513957865153635627398753535508"
             & "0322611590445868485306599512716862439319443122460819"
             & "9801290665767118647436056186075033334804304548797573";
    sint001 : constant string
            := "0.00613588464915447535964023459037258091705788631739"
             & "1329356732460300956624454470484805800453186249319642"
             & "9611518893795760995077725354696160560364182978798419"
             & "7170913458379002530241158166175165004060371009282128"
             & "0573717942734091225007562451633302925160537935313896";
    sint002 : constant string
            := "0.00920375478205981931510237841519142884818342672937"
             & "9813606927043443746916043690861345051934265607789327"
             & "2758135114920908757367033237517402714374646875767663"
             & "5343353242221616312144103715137717862631539429830708"
             & "0618885021089357144706731083897788982488430008440714";
    sint003 : constant string
            := "0.0122715382857199260794082619510032121403723195917"
             & "692500382416732467478717678330299039555616638847753"
             & "459942482871572428816519402328721831195993001959512"
             & "672051014416878185831425659567890131142951801219352"
             & "9982280535578289142648581602718225716596996534323632428";
    sint004 : constant string
            := "0.0153392062849881010441518676024626213017450937386"
             & "568930147905504326236057284121964514572183109886420"
             & "042821967984041860963122632162871316446040099308561"
             & "279768626396311608970618436586513591309765787976481"
             & "3074157912927856187446088067921240344836289829680249151";
    sint005 : constant string
            := "0.0184067299058048209273663130148401265504987492083"
             & "374102969438498348877742857037604055581045570636510"
             & "261339675144237453864334233783408175088965374727637"
             & "176755311741479238835352443883026546163977474723851"
             & "1681092153159332405211252411021119927469856942059895614";
    sint006 : constant string
            := "0.0214740802754695074183748977980622605209359681178"
             & "804897079228466408716203524793993911834218045420257"
             & "543395379198009567452564962664819120808457796951678"
             & "772823681365274215102138236979413419832905098499024"
             & "1719574966240060074390379416739218676618602659348243060";
    sint007 : constant string
            := "0.0245412285229122880317345294592829250654661192394"
             & "514775767567737884929481373947340741505810955139589"
             & "365788487332331929387934115663976142067072489878996"
             & "237725091251284090293371335173172901079602678960474"
             & "8354282010595009404887009511244390105868368333663908763";
    sint008 : constant string
            := "0.0276081457789657416123548717439784658641636872579"
             & "351767458936847395238360869290068361719678590186040"
             & "785858717109309217921592142577786010417475788328785"
             & "731873206629137147924485705869715824662316714617523"
             & "4688749259804587427465728562981265957322361227464447105";
    sint009 : constant string
            := "0.0306748031766366259340210275652237129977158501934"
             & "989490493518369635002304189534882726798972483577776"
             & "422987552702651990815229475822676589514437548754446"
             & "270199662938329878026483203334312586655586332703820"
             & "9321553448010535164572644644985042700095912074573007685";
    sint010 : constant string
            := "0.0337411718513775848337161124064239499621908244133"
             & "321147176627523000057304347553116142168312633657961"
             & "743319441882865183568213266134936311975498722130657"
             & "840592265171367438567019353637305634738629329202529"
             & "6662652417337763068978514620475523016952260799839358951";
    sint011 : constant string
            := "0.0368072229413588323243326909279513010552317719584"
             & "573991270406521191349411640973866059376600351670385"
             & "798327841552970047925060504958600243261106741808270"
             & "046484442163837801971719830648690645640166226946454"
             & "2451130578282972010854273247468703731035981849550518730";
    sint012 : constant string
            := "0.0398729275877398111285787376798802131956085404986"
             & "545607168397113812318940862671905262026720712197533"
             & "124666584087755293193544263534486173085493093671804"
             & "839390155920623951979059740797695437906239831689839"
             & "5430047468252147408027413980594086538104570558598826514";
    sint013 : constant string
            := "0.0429382569349408230771245402817839553939332746003"
             & "164055238510235477711595141985694950520660025025156"
             & "667624298045354844676657250588037367363778661171256"
             & "892722055475232702810452101744532283597430544899704"
             & "3618276863038126963925618173235568362353014568663897644";
    sint014 : constant string
            := "0.0460031821309146288143017879102414734403645628801"
             & "823322697973836392193145630803181413947004672597150"
             & "693000734244740696562936275753990865361407529533496"
             & "644027942100000827169673244235320781772279965497199"
             & "8687550792436000588850907304095419642498648430639406400";
    sint015 : constant string
            := "0.0490676743274180142549549769426826583147453630257"
             & "529202101225532691659395639607568222207671118007237"
             & "540668609774655834613187574809779385785690451360149"
             & "380553514179867719318954888179493718900818971600222"
             & "8623542241319750939917966848266897030702994819069823630";
    sint016 : constant string
            := "0.0521317046802833212363582164233367485583059324482"
             & "465249282482988320236995506600617153979447938540006"
             & "359257796288162846717791052095728738135811760473591"
             & "696414879537947303039086954168762147453759182724702"
             & "3402085789402809407225495605550743388707545197716635035";
    sint017 : constant string
            := "0.0551952443496899398094475256977309118980986639932"
             & "421058784004356548895024883146963677035872057828874"
             & "079553201161895840926007337946633327781420690096880"
             & "991680503248073570021880877098573093405359505653857"
             & "3224462138417844469731692323391734852435362674221734806";
    sint018 : constant string
            := "0.0582582645004357596139797819343505758564073174794"
             & "686451223567110178390567524446152906342167667826943"
             & "804277369114548228720300763263957853191500310964907"
             & "720842689675332472233564398028205706958894205283551"
             & "7694426643131961195841351247263500205725793290879291110";
    sint019 : constant string
            := "0.0613207363022085777826145929172350071939802796769"
             & "731652631930233073550398959628364542397455942871312"
             & "109220044547668386113198212435856873758583240454584"
             & "771016495581247931244038593828827252917330996342522"
             & "2679445761990224552633549034954831198443142418749833400";
    sint020 : constant string
            := "0.0643826309298574608193245368265668401635526079771"
             & "437676241010823076329196026827690263403204294406233"
             & "089670535375803693111959614712591108150005224738890"
             & "206798867727948528315780127480063777353426921743501"
             & "5078737406723360884096510690574597506260579369211709750";
    sint021 : constant string
            := "0.0674439195636640578979724218744933647690956421490"
             & "120689450947962656047478096050136701068308041516702"
             & "435356342043681318502177292786860460065345070757209"
             & "610949936711139079431736169585234270834935031487883"
             & "4752768056050406264620934478720332170104165192837483092";
    sint022 : constant string
            := "0.0705045733896138630273514705529420755992491548764"
             & "074945726869361035649510365154767643657508608334975"
             & "274536395776674014163212221917786002297770600385764"
             & "688301293746652358473163226213626864605749895159464"
             & "2845693271544392397435913542775871013170384564658317606";
    sint023 : constant string
            := "0.0735645635996674235294656215752343218132992656778"
             & "872442091692001525149465249903079306352217560656216"
             & "324609054430876344867945030187951644162744720448213"
             & "477111777219469329366409740312702320827652139380915"
             & "3628542243981086758525772363352685627108031163857891800";
    sint024 : constant string
            := "0.0766238613920314922783324628237898766848511890571"
             & "246492549418814500150275364564006640774834474325215"
             & "384700270306597144844573894837870492694628049235077"
             & "617204154430662641556682501212470627535298816953422"
             & "8617967916615766513229918075175320919207840530933106715";
    sint025 : constant string
            := "0.0796824379714301211471206559958854013639452967811"
             & "396339998965642878222905090589215211582134818509289"
             & "718966114092485820922939832031838626386355584011704"
             & "236912568047700514643385029003483882006425453044430"
             & "8648320643960833803693447416929472431235803485660023154";
    sint026 : constant string
            := "0.0827402645493756931119870831858619740379930597113"
             & "823915502305770040283986563622765088795014629951914"
             & "770441608259401351303687694281916991493802771995453"
             & "857634057678981313535216878106130737244710298260371"
             & "2198928074120547104613473217161915806696621573259276388";
    sint027 : constant string
            := "0.0857973123444398904615563321468461719405315153026"
             & "428988435681912372625497215351651236721322280041695"
             & "357898854853030448022608358677423461877375379488804"
             & "018749396370728194842009283062309820465584110238854"
             & "6823715831039211974039826054842959344438883912115125972";
    sint028 : constant string
            := "0.0888535525825245965615865350033739514293329062219"
             & "376145362942707169527566098113804037732924765816896"
             & "850521270700232545878414855146294550981778121797991"
             & "024125427811273016779262335444583803575123100171307"
             & "9822737776806968925496072523186674280801778043777321794";
    sint029 : constant string
            := "0.0919089564971327286249909790776994870092419232930"
             & "348413509974936530677796700136083525674292855313317"
             & "630259593674569980528889774540608803488321904926027"
             & "080565929731058677988099566467868506113278143631088"
             & "2864862862874499269909547282807115618483424481860774126";
    sint030 : constant string
            := "0.0949634953296389989380343123604925943134232383578"
             & "768199767976229285364047420420757855103861804731497"
             & "390995139145594976867241022972862875298015227534872"
             & "416586511827963849816496448902883696329234907844196"
             & "2554754630639650259621845548329323919603248610292393502";
    sint031 : constant string
            := "0.0980171403295606019941955638886418458611366731675"
             & "005672572649798093873027890875368071110771463185595"
             & "540742065264441021568752468430244239341285996871480"
             & "855884119585703243445125198387172878236055925742268"
             & "9670227492748768318525076200382595233254235685667725190";
    sint032 : constant string
            := "0.1010698627548278249878875845507100762213559985945"
             & "653472145443054814579275842923554740620153545693653"
             & "091531829221536761599796765098187057452103546021610"
             & "404272109336424468592186738155298739938467065159281"
             & "572531318917833503548997772466150745253269350958319098";
    sint033 : constant string
            := "0.1041216338720545791209438800601790207959749274903"
             & "091947185613660056851894964907800660776969635807967"
             & "793521043619393205843427394709461672155472470226970"
             & "844019631320489043047699587581479654803957464144350"
             & "669571288752133952790657286633128439635728934024415955";
    sint034 : constant string
            := "0.1071724249568088491755291482281967503779193796583"
             & "499779700652020405492612745401432526400745163454984"
             & "421670320580743551347148463376499902506803525031691"
             & "646537782532576257495344124364369528458126951325619"
             & "162419985224385260038189344018789257142454884104150879";
    sint035 : constant string
            := "0.1102222072938830588078991402156777252744746231149"
             & "868762164343487878463726789387838462304925486091225"
             & "683518891808825321976681351136528770483671615701102"
             & "679472632710649835293937287438348307532705717320882"
             & "955710697100171962351703623782582749391145241585500571";
    sint036 : constant string
            := "0.1132709521775643490182287329082850660797838690096"
             & "898863359229809609257813975054426989624128029689741"
             & "772784126723856065215742063210089346277329992855128"
             & "753860969794354388534626425217940600765177709784267"
             & "085464898698575472808861153609264262984212648383588543";
    sint037 : constant string
            := "0.1163186309119047672525443194705125409238566764922"
             & "737103284823008114364699561563219959396264027906063"
             & "467476428403041620588019002394443189820501050181923"
             & "500380509693199150768950723961444478213668053416220"
             & "720397297871454789337099384719149127074617325704620226";
    sint038 : constant string
            := "0.1193652148109913645936377898047947467110184824121"
             & "012990740624445550955097871079956241295827888490384"
             & "056980312542388026903365873077616922457649905955944"
             & "757878296399199798918176198984656267800034876909715"
             & "913659191738745124207102320863187257267525654653343752";
    sint039 : constant string
            := "0.1224106751992161984987044741509457875752236090851"
             & "072485773187024063510194872256293642941583891435347"
             & "409939546081853122649616772937850341507948866397664"
             & "753406951047556362671238322436962799348839486529642"
             & "348015668517234046237169354524835298908599115995491652";
    sint040 : constant string
            := "0.1254549834115462385423364532675945494325282092800"
             & "867244102468434687594777940698904771232954950995269"
             & "288980458298673555649579429096528808638400295735879"
             & "299821244874561252176574503554383609754638073350062"
             & "815145482840202379122567763443470647945932846098012470";
    sint041 : constant string
            := "0.1284981107937931726244155891727575476975072849139"
             & "067467049686590181812159401129116307688479310509669"
             & "045348150715432392326481143583722247380779798238448"
             & "201362483019964964829039730822790532227940101889563"
             & "822124182581246447329872805273776991227802205623606524";
    sint042 : constant string
            := "0.1315400287028831111033874926922301510074010766035"
             & "916890220473051757135560484005290470858672195362241"
             & "148884302214104376328787765629454684894489745662233"
             & "772976325824571755920416493309221702035581837788575"
             & "779436733704254446208495920417787079326644790342964455";
    sint043 : constant string
            := "0.1345807085071261863163584092539792556314246096030"
             & "667794184456873039562949256817967758427591489535350"
             & "168066922671824341747333827186486653293605070056698"
             & "959002398601621324445413238111431934149052490755822"
             & "002144470173810982123970045055231035442448640070769435";
    sint044 : constant string
            := "0.1376201215864860449484416634310973343662606497934"
             & "429611714989898675814618394878177612574414178926679"
             & "401689650428985376505635603267135356629586877378461"
             & "315574959420748765046964170632247463024633930392008"
             & "138141846253929312456211456256717205270980542397526644";
    sint045 : constant string
            := "0.1406582393328492307147888464071430911197878898999"
             & "696127538723881418578855043578625789877994288147213"
             & "066043478974379031297686872655070320394462509779357"
             & "058667327591398212874691374585182123512089274480069"
             & "654291338824921570826090151954376983606783895569430064";
    sint046 : constant string
            := "0.1436950331502944548197733493230505826203393800429"
             & "713977001593704205048177084676510985739254973016784"
             & "891207091301743130332593932601433318232007503781481"
             & "102542884967815009705879657819801257794855949732234"
             & "780240849668598473587196278917834417016304211112973987";
    sint047 : constant string
            := "0.1467304744553617516588501296467178197062153165293"
             & "920668973938768923171676867758277059618817711050948"
             & "125690558476379337616896667576365627798992375193448"
             & "833326561631393463597419204096728033984501101937783"
             & "614323535729485119716127600641621327886953017830762523";
    sint048 : constant string
            := "0.1497645346773215172296957373427385287768305173155"
             & "205944444038558652026598191276150124939773224839797"
             & "345312570326512200707512553366125599860711241909187"
             & "453516799341266358874683733212830759271865840900341"
             & "136269877363862446061410481987660301622448646734957815";
    sint049 : constant string
            := "0.1527971852584434277203366125438313132217351043555"
             & "721869605997356764876132128246479495405330992447897"
             & "472145326702096385561918941783925779088822673452488"
             & "739252196716577357670366164825274514274857991824966"
             & "636899151356257697010867245179045321923076702055563765";
    sint050 : constant string
            := "0.1558283976542652357431014862462223493747859198377"
             & "059500956243530417661606847322898498362634696292005"
             & "485400660968208270900407390008362078124824964360543"
             & "303344109289002991395586506808451659453163240952638"
             & "288261487490112205024867793377242268570847925455518167";
    sint051 : constant string
            := "0.1588581433338614416843853596530813017855264604464"
             & "552249536981739616921261330337502389330825502347697"
             & "979805001620277873683828850510474928183192216361237"
             & "075115525862737834644173067745803004463514109016380"
             & "553759324216540653168408424462708618539678111308095288";
    sint052 : constant string
            := "0.1618863937801118376413879953333824608490014042822"
             & "816066257217443298522299444388617851003822303615474"
             & "930434841464796142832371074545838486274371166911628"
             & "026771518347780174983632872004108593410599549821608"
             & "866782100492586084816434288594742360653864259622191673";
    sint053 : constant string
            := "0.1649131204899699214181891132844101245584977531171"
             & "477283189379971883900593364668195125454277509079100"
             & "427221712548185852670384363234664410584761879606107"
             & "650195324179545177461408255006277064443638940859714"
             & "746833172945912359012809362482617034810914509764908717";
    sint054 : constant string
            := "0.1679382949747311780547455359966578549999671775156"
             & "040891137248480506205841504574665679184078199332049"
             & "201103520028610731185279563706752889587328963301759"
             & "504319000652561991169576322091783882271347230088386"
             & "055827291446312266363795354878758495346558958876606724";
    sint055 : constant string
            := "0.1709618887603012263636423572082635319663290591445"
             & "910012596894557734102929076287581584670239143184040"
             & "811875435497336675697952852373347815029074446204950"
             & "059668690728955902116739590027623967956353026267919"
             & "498878869306408619679829474981903547134184038961093511";
    sint056 : constant string
            := "0.1739838733874638279507008074667317502779946746329"
             & "410658045649990593665397244083898604295139382845848"
             & "944327780482637304591538337412952036844765709172661"
             & "513075638264612846984370135953854587192006565975729"
             & "841702069190697145594978414654696102213530880116822207";
    sint057 : constant string
            := "0.1770042204121487561968398439291657307197507865185"
             & "694054959712981727525072009635977053491894129344306"
             & "025520576628876718194873732523220888932215432740164"
             & "096241685701804784388585075603069845836759642905804"
             & "200954132951522250274580286290061652600102050997113270";
    sint058 : constant string
            := "0.1800229014056995226799065898456067867505057226925"
             & "682852574957995762482466538285297798910671393036846"
             & "386862703291663011632115912097538903578505788092772"
             & "248273897533910230192994132080892486456413391748500"
             & "164511255885877904061849050005890592863701810079172420";
    sint059 : constant string
            := "0.1830398879551409585165325784769200138775632254052"
             & "333280565068804504059490554431319073868536585392761"
             & "499461027105354717875439522251555395721962205454182"
             & "320526933331265975559627069012782771253984502147933"
             & "385329539098391968409079762502579963280862247258278951";
    sint060 : constant string
            := "0.1860551516634466481054383041691617772795158180940"
             & "969962866574421467013568062246697956552897317269569"
             & "331482367632973662863946295951616450809587911975473"
             & "674657386882583528431376745425647102973257265863376"
             & "330875769027412122955079306998487569336277751813046089";
    sint061 : constant string
            := "0.1890686641498062127549978370874602018495823335434"
             & "643580564464916344443836388470574133991029539869581"
             & "228698172865196875361556939155160452506329351211238"
             & "237667909623178989091546264952589338465747745307156"
             & "930684933188324884226393202234579140131982977856517397";
    sint062 : constant string
            := "0.1920803970498924416792882046279485355947208168729"
             & "749228413443541811512053491245791265876135572673897"
             & "645255843973394487309229842331353089630746875740030"
             & "204879999682962117645631317421701876431668531012999"
             & "287670893395925210154842096705118028090189272834116191";
    sint063 : constant string
            := "0.1950903220161282678482848684770222409276916177519"
             & "548077545020894947633187859245802253253092340903817"
             & "309920701055366117589657998348476073888265914344275"
             & "520716094619688003685414063409907786724828762947800"
             & "990803467764437301830279319175379898612584566670257313";
    sint064 : constant string
            := "0.1980984107179535861793249181510733849221204375709"
             & "235023812777319344216088053886274663426277947312680"
             & "136000982196113444302643852171352778723529102710797"
             & "534972444622697128542874157879250586473757295790338"
             & "702862145333883414777964538306558656786340865119197501";
    sint065 : constant string
            := "0.2011046348420919115584435458820671774655273168648"
             & "817352323864815781581188420293117431478902484563306"
             & "608924162163486463145323060741320960191768733194022"
             & "693343927465257251082827506047225322277038145131235"
             & "986162972064445897187782503421976880121387535673838033";
    sint066 : constant string
            := "0.2041089660928168741816969499474189576619453422964"
             & "934729950150565701709863573527480651505149501379071"
             & "287229939319257473383169911550172746049339150215578"
             & "256240664481758564098862074461541420574502274335013"
             & "866844534312833470450091125278671283893203902091348965";
    sint067 : constant string
            := "0.2071113761922185497081160197866674167674139345372"
             & "761775960795601783456612437330322511583414223902359"
             & "857791028401115584741293058922090368444614494912337"
             & "377916361252568235462107901388344718125137502062862"
             & "658535070916193699916785817308362225896522474309158940";
    sint068 : constant string
            := "0.2101118368804696217174899720901250451970695512599"
             & "436652108556451898259504376085813325649770923123670"
             & "942402296479455089230639808555125946624689424916878"
             & "332924930259792731972942926953666559164809329129240"
             & "089170766068038677686990953013294389294030394375801824";
    sint069 : constant string
            := "0.2131103199160913739677575178515008423979589821936"
             & "905244000263807712401308857601312582098440534995357"
             & "671946476457314255418587673490619590764316750261943"
             & "951436586938322126611625397956059690910488794416960"
             & "393861570308218152600602425278478827439736029369896971";
    sint070 : constant string
            := "0.2161067970762195099483851312908295838456946552854"
             & "169793529831931866315990199882245978569335662345437"
             & "291519313656282395853273777053172277024153243943778"
             & "586547984477380395889635753793125098768337450564675"
             & "317863608052684616568245355940651804101894976634690648";
    sint071 : constant string
            := "0.2191012401568697972277375474973577988483607967055"
             & "921085026246329098112435823550704068816745768548328"
             & "961855395039949531274554839545365041418372412606732"
             & "403429514280481257456960184892633828220364417932957"
             & "553370833985986845998388623530175318676566083248427709";
    sint072 : constant string
            := "0.2220936209732035340940947213139774485664983603248"
             & "351303370813517104400706786900738912195230628208055"
             & "526427635552612383947562023820420103411198034583547"
             & "578192138906088112012630650079341764917059914029507"
             & "336875486998363423604546349704585505980769383558228885";
    sint073 : constant string
            := "0.2250839113597928359916421198633534633451336324953"
             & "449929123582830508862734414983307301963032022948069"
             & "632054243872926583990779258320011255789689040176801"
             & "332574006207010964627488323613373781685986245155358"
             & "512601227593755375499905866535272062549470596492496264";
    sint074 : constant string
            := "0.2280720831708857392544573794575372441626994680061"
             & "817717327915962597737135789155123583246830012698831"
             & "926187690978109363291607343971403156220065609660520"
             & "016326734287098594850199715179208562163710799726765"
             & "949507270744941237491153331942184051960213117438883070";
    sint075 : constant string
            := "0.2310581082806711196432360184727066652304106304332"
             & "499034754787627016066610915480229880729291823228214"
             & "413711947833422784072545470803196652838243969480501"
             & "494933870353900912438874340657304242177922777068675"
             & "855022129615151471680667179360611661816518704575613161";
    sint076 : constant string
            := "0.2340419585835434231912420449226212800691004582593"
             & "853317524045774735367698614878370152011238056976716"
             & "759899296341588527702813485902374565522109962745009"
             & "689941964751224090390391236632180229725398455909107"
             & "147772653181924660841205445510482723457141552553801101";
    sint077 : constant string
            := "0.2370236059943672068677359145212644641039880051746"
             & "337859577920504773580535296815838446941170393114248"
             & "363487287172272781497999602026265786719359455676804"
             & "526398599863471191632310017510931722123284771699256"
             & "452562683480845118108207097674548998454002513507941968";
    sint078 : constant string
            := "0.2400030224487414865689223653588895716013486111028"
             & "527542367150976672180677630304581666346936817822285"
             & "389694055052717095204075958405714163449471019194094"
             & "142221257438744223047821530193119615937206840586967"
             & "018664019863403509754195413215624455832494706432754615";
    sint079 : constant string
            := "0.2429801799032638899482741620774711183209907832838"
             & "321260423208417946981663958870487838239414776082706"
             & "567113809212419797439857436275167386038607246447923"
             & "207968504520696899998694668136932500968344527727215"
             & "735712377584170894975493995644857872242814967880556183";
    sint080 : constant string
            := "0.2459550503357946115999247085519682118632786652695"
             & "591815243556320971753901406174334260709052943087037"
             & "003758868440900660697435006974522990692189208064907"
             & "565233738543296985792582449596829373497936309974105"
             & "600364872649854004941068921097489507274661396175553510";
    sint081 : constant string
            := "0.2489276057457201681106828162729877050769317717225"
             & "214148620243492995773755201416183215765857822061422"
             & "872445175780301496531831751370206372035276482653783"
             & "339658756767741896830541474591048434206828052853759"
             & "443838656950864028143655624475924135613748671343975498";
    sint082 : constant string
            := "0.2518978181542169504981066283742714258673799187698"
             & "320950278639112246802345119190755183708115090393417"
             & "316191419001853581146727945127780510161810642567920"
             & "964262383782922375755557454039051787172461172765720"
             & "711435310202169735067846471753860419070753681825078661";
    sint083 : constant string
            := "0.2548656596045145715539807788247035490779248168185"
             & "416691847029015242716092767782953110555934130876439"
             & "041182937402705182787951330575453135979999689649679"
             & "507543581665630041971860901186446298517294697041028"
             & "184852519063728519692632283153661968987935826909267714";
    sint084 : constant string
            := "0.2578311021621590056144712947590814441497126759390"
             & "575943305003334694421679329753588356942119043161244"
             & "808915630639807111173031622253543230624493013092610"
             & "465256488176358703579760725379913467669274913318737"
             & "053915154318330549182704833280228141637368195938334375";
    sint085 : constant string
            := "0.2607941179152755182801865090847883426487919206749"
             & "109574696319510928797532549056279570119960704052461"
             & "438751257556024085926528646163364782215669166121996"
             & "931393470099571189836565746930885190353498615891651"
             & "184707303022415160989947872487287734429277458847807983";
    sint086 : constant string
            := "0.2637546789748313836113493219832235967119934070276"
             & "339532296885820048967988787995909934505185063017054"
             & "696760207946376998182446890249667969547375557195878"
             & "896455639214616776649133409326039693311030804530140"
             & "378028219335307659494522560305535959669516122608738018";
    sint087 : constant string
            := "0.2667127574748983863252865151164363940421169883561"
             & "562081989024345080212831231418753245630552591635424"
             & "722296970341971436103120233338649903679175939660494"
             & "374539713647943961115171968474533409359471590690859"
             & "149803630348613531968833549981618476493350480812814829";
    sint088 : constant string
            := "0.2696683255729151065254644624226725396346525991952"
             & "706398329285603879918784591839717586333477069003868"
             & "284519288108507036379790919840414837440409844927657"
             & "957219584307834148360188416006215608544420001442929"
             & "284583668866578956422292308982732978935504919640144798";
    sint089 : constant string
            := "0.2726213554499489844933474772922102408041755580770"
             & "619026594726427921093737282625702489119160101349605"
             & "901676341992112685083878024714705085895399065584807"
             & "044705396486685788119867419351878225069807273465631"
             & "364612066498592266516260616866110321662058335654719589";
    sint090 : constant string
            := "0.2755718193109581630764251683907659018326603860153"
             & "325642733575371349810524607630141097230611496034073"
             & "543127405812274472377540226685218718471742735953664"
             & "373445332378096592067127043360854033726789992874433"
             & "674946128995773379762237995684735696524028294465247101";
    sint091 : constant string
            := "0.2785196893850531052078485259570194793626584607718"
             & "518665295610233074975007968133888519732729596325275"
             & "491261150862709086548130437666236208319032222344214"
             & "285437035335637193023570176250572621387935350180220"
             & "359045577586420376801747161873969183512023044579356066";
    sint092 : constant string
            := "0.2814649379257579840952310073400376703951576887475"
             & "711932180251708887761310156844315512315350430052172"
             & "224214066170092027372439983515516284687156891265323"
             & "447200593985812641055045575365516557101060318132136"
             & "955164628816373998519684545704382819111056056949385981";
    sint093 : constant string
            := "0.2844075372112718436183106149398162274615807409638"
             & "922373945413262728690190259537764158871113645460130"
             & "107533083835695899356704063070873893120753615774534"
             & "285314141171708031876018012341145181566857779317459"
             & "277173228102898148492944436877432773097970526768057259";
    sint094 : constant string
            := "0.2873474595447295264773318414299190824647101685673"
             & "589622734061895346010495925554047417553950305864075"
             & "616078182247317656111540716627721413082305462594448"
             & "566429522292817727932160330628470925707013973039940"
             & "615683906240542627418515764545777076210348725211918436";
    sint095 : constant string
            := "0.2902846772544623676361923758173952746914762783241"
             & "511114206671131253928982994574468532436606893438468"
             & "421981569266985974726632161936158765975844925431147"
             & "293592719478290267478100097575498632359877236161034"
             & "118169049803847190528699697756632759644662618586927835";
    sint096 : constant string
            := "0.2932191626942586506066085989561645990749759199300"
             & "713653651177739816498514982610098314366431463802979"
             & "988291065963619930004728395109048801084193351040707"
             & "420102400944956011921763431745377249165261013162619"
             & "810023979033442825792044586394435435536991046088316883";
    sint097 : constant string
            := "0.2961508882436238241217861277826588797865038519096"
             & "636779958151129536903429743522448388946435785091030"
             & "929190235355439125946217441183235181498849881735706"
             & "117656099272275882516704680278726845297444047288944"
             & "592892473229338962165210681039673876460102080129589179";
    sint098 : constant string
            := "0.2990798263080404767503369727760827896234515469302"
             & "238094621249196809632193386632939258427223608021178"
             & "241407958947660643493862579352952416641628024837104"
             & "604457714314759067617926524247402688563665800193988"
             & "851176556741485814609570646722019517498660181498694739";
    sint099 : constant string
            := "0.3020059493192280670034632317324243912848169556015"
             & "980256384425145434035973129048750553850354021934377"
             & "521859926303677858550070242817916888960925671089608"
             & "320816500718354751726122314936879149467778599159797"
             & "122113497351556706707621663739744864705676369331343254";
    sint100 : constant string
            := "0.3049292297354024064907286334365223463192519271026"
             & "275918567468129123724603911315960759021201024090157"
             & "123366254791230512527199479340994613222001782536372"
             & "830403246026787867833262254298464656153613126832040"
             & "940690269683651617614468954251785617197263908993002898";
    sint101 : constant string
            := "0.3078496400415348936820636455585202015094184912737"
             & "627553150503322770208324630417182556808257602291804"
             & "930959198088831064734309553913282524033525352458482"
             & "240890319823009148257384171893174835391414044578233"
             & "383924569981911867620822653540022730299159497580666334";
    sint102 : constant string
            := "0.3107671527496114958359972502117632942196132072684"
             & "204428620846281659509548355038269353084251540200658"
             & "468351705514520539825219448442043732174700328541197"
             & "150642183903611387172495182283946097971498710964936"
             & "436121544746259075382083432725211165493532263672499337";
    sint103 : constant string
            := "0.3136817403988914766564788459941003099933775094565"
             & "467851932847388033398137344679026583246149717961566"
             & "235954198732386085078444263815068892657515312999186"
             & "282843991243383176028844313277967411434326092259686"
             & "100712455065313793685805060345543304379345034579842623";
    sint104 : constant string
            := "0.3165933755561658672430470346829517504096969166492"
             & "718883162986204146192975248430663290601076797398273"
             & "531595785673673232792825909766963497280251525934332"
             & "039184840814342829698594436654614045980099458053513"
             & "985773124669232315735065355859851305676091941660960909";
    sint105 : constant string
            := "0.3195020308160156779015182715397565772982056106134"
             & "798799855121117012355370919647825868253578023825912"
             & "881893441156290508573647875448029678151471909056012"
             & "464306249756156734528431239162934545465381087436306"
             & "051310382015314505454745158546160780344503063276045219";
    sint106 : constant string
            := "0.3224076788010698483848074776591158480870902495070"
             & "670106805740766909232272165044162817919786326175238"
             & "291365805299618719579486595649462130206202821498970"
             & "180897827517997424607168554911914575027436970347189"
             & "640897637553680890522296222724722209135292890882744018";
    sint107 : constant string
            := "0.3253102921622629341359547080141967703657466600660"
             & "294553285869535761121159776102343671524645071413736"
             & "559779131711234838468962730339977932466888263798560"
             & "197803752738420522252007697816560858880475209193804"
             & "173441222051189656450613443424349954891397829582958760";
    sint108 : constant string
            := "0.3282098435790925261079168166295056612837898419418"
             & "725944234413214229647379607651698191151528544065564"
             & "816165380066829323599834582312058457523599322425628"
             & "129700121876389729069286264665978486319037177652226"
             & "352149237239808191054062171413651065731683781090562387";
    sint109 : constant string
            := "0.3311063057598764017371907372666501496458861183171"
             & "805551051946784530741512178795861156507461644458143"
             & "553377268515519956600112161840882832345317335804895"
             & "979773646819421527476696825013788373072052050908208"
             & "000645007216697353918181979403669340275843200933730586";
    sint110 : constant string
            := "0.3339996514420094046508654805349790667399323599263"
             & "156946221339881600014597038541584093163038424464117"
             & "165422345402942207689664860133289652187762264869057"
             & "936681671403632647034137528929243768356238759791349"
             & "468080113267526841539545310280498466191712764698165506";
    sint111 : constant string
            := "0.3368898533922200506892532126191475704777667796712"
             & "222840515309421523592933598104109594588644731514875"
             & "566646455342412120020059137247381125522626800604920"
             & "655000913884180437174946336159744752157047743479741"
             & "528451424800326129331200057946292135457048672499657130";
    sint112 : constant string
            := "0.3397768844068268578288258028174093267566821180251"
             & "100390776254190116119713334976331110001903719930996"
             & "126744385025482773285054720339441091322021015242110"
             & "437372129433226458125502597578949378939184477227574"
             & "349439362804542082254507410365172079633569151545790316";
    sint113 : constant string
            := "0.3426607173119943975927819825612868796115591641438"
             & "863646620820902565103339096773237999505132004361568"
             & "708638964487908550150105906557533294004777672240633"
             & "927104908376312475232111021742705253350993011915831"
             & "618864536168239092146755694417593849234483935449952416";
    sint114 : constant string
            := "0.3455413249639890655391917230787393996813489848402"
             & "525087688677092509107784472764210007539829478194774"
             & "765170063260130437397825102845735880846274702067792"
             & "661127484935109139640168434253789540148710219349376"
             & "604193750361927511518695902404460112139943431723278411";
    sint115 : constant string
            := "0.3484186802494345684193085876948116364531209085280"
             & "093224679374582799681831732945144230449272152019415"
             & "834710443118508175920770082552196458983517741342187"
             & "953065797985376580915959088545868454325023334863765"
             & "549736523676681316688674450765320715994373851764874579";
    sint116 : constant string
            := "0.3512927560855671256013076230482732902746930770598"
             & "977253119215370091420174209160787857292082854299676"
             & "870584126632072870259970950803752173976736304610116"
             & "919933889844126837419392934111615943480083574557555"
             & "190706326921945234861497411233060479772878256817873563";
    sint117 : constant string
            := "0.3541635254204903823573957961389289436032296670206"
             & "147338379015563202835449557421950661038042795152235"
             & "634498249296458801321798764353513761013463355676224"
             & "755876316357503135695537511046104818126444661656430"
             & "896806258218828937928591451027672272393596543132213206";
    sint118 : constant string
            := "0.3570309612334300326149540357940367878148326224386"
             & "568290621142200414124655942644050783002348336093305"
             & "425407695875424224171316443481546952504823107767903"
             & "467141354927502940202972866089557383182434104621043"
             & "233250806257035351400555706822672347605920943126325063";
    sint119 : constant string
            := "0.3598950365349881487751045723267564202023174211290"
             & "258497630100778756401389522703121629916592490241555"
             & "351151531079351839697550218630851277756334889268487"
             & "566649028723073289664785663939491949787031415817032"
             & "259848972030576512147653109977094545991390267552669793";
    sint120 : constant string
            := "0.3627557243673972162048544621153241297874011048202"
             & "520155940190663014948554028365186983075440044111335"
             & "827799881434392552984491206801634131858371646854795"
             & "524844938365369340568906744615615541927698763026329"
             & "742041449562206580566382819959960888806222382883699635";
    sint121 : constant string
            := "0.3656129978047738700117459086069781700085140664294"
             & "972528903290581509153243479155962255861874122016102"
             & "017092598724863162982430436751006030357307971349089"
             & "084457710880269669128777370856517378855559772580615"
             & "870254354560164988967385440349088101308780038446204904";
    sint122 : constant string
            := "0.3684668299533723317127462216816982941984166163223"
             & "511021535623999869634686700797144574672689157817468"
             & "955998006689115009732886045956264516413146410978908"
             & "226324905969589166393609331994268237677202037659381"
             & "794052051965950353868947019857382593727594058561154773";
    sint123 : constant string
            := "0.3713171939518375434119349670219232661775725301769"
             & "655376587020886674745864278377441012715035639620971"
             & "063021496578278136708216168846279182766718786129485"
             & "213981365778518731828373379243793492190091159073842"
             & "371387061034998585954009438639955464693849424718071008";
    sint124 : constant string
            := "0.3741640629714579971043930195383235683167971603106"
             & "387167716978305976400246386965264533249658849709630"
             & "537819958693450669225015128400126843381152690200790"
             & "484195366895662142664492182476661442727761707355514"
             & "612937633902544217622876476693693795729803815228344611";
    sint125 : constant string
            := "0.3770074102164182567265678231998572323015378836479"
             & "092748594578595133254809758399906775242029870683486"
             & "807553121426412740061805531281941009270258496256196"
             & "834755034543957318098508846337520980048691290228970"
             & "137058937092121303173402429408342243784221170530792620";
    sint126 : constant string
            := "0.3798472089240511705762811467990666753860053150308"
             & "724931364148298657109440662952902056236223596557586"
             & "600602724894311144766027744437971613620457127339644"
             & "146722897004247871366964698857432464011314000376564"
             & "453552373559288499000969530180196827043253923495298914";
    sint127 : constant string
            := "0.3826834323650897717284599840303988667613445624856"
             & "270414338006356275460339600896922370137853422835471"
             & "484242886614935559007560102009675979208442091777288"
             & "702111639612022873816233516975958506707886601926715"
             & "740526121248262215765061901640578812361838848504530727";
    sint128 : constant string
            := "0.3855160538439188640756079493391946817280674679258"
             & "861679626726611129511902898642133121046856319143705"
             & "318248486199997255852469189307950333625124015310223"
             & "484442734344383073857403243200335825547831823954400"
             & "726381362256631306660846101812373463897045713578725851";
    sint129 : constant string
            := "0.3883450466988262916249935406705281014962833494394"
             & "689315956273329869069584940237001847685357330624351"
             & "347671191899549031927966660847610501348041435309082"
             & "406644023680286871125839960706497538555120977406903"
             & "706266649777002141908953340390952301872584782516217742";
    sint130 : constant string
            := "0.3911703843022538886875129486588618994434759784047"
             & "975256572286445776044968494016656469734591967348379"
             & "205194332085953878914898607410940176466723796537697"
             & "174430765755479847475834303809822193843523652669025"
             & "694364346557132488747198697248971570861555595216846903";
    sint131 : constant string
            := "0.3939920400610481085961886608903134244856522968110"
             & "642025862597880650858688180689520062438121892504346"
             & "348068717173583751488773992105397816400287258478083"
             & "786057334438807339940844005218523128004462106217149"
             & "126150060954107484684615713375286626360722417496013739";
    sint132 : constant string
            := "0.3968099874167103285952909113684558637493841960279"
             & "746464374813170959603329888314450244585686125729661"
             & "948025477930883642129567432813096298359270732102739"
             & "940165520134800978396106814774829997179369279844252"
             & "268501252450074422731901068652822653661301700436789538";
    sint133 : constant string
            := "0.3996241998456468285441170307420208607299235371692"
             & "063270556043397293889007199834800706515987380019972"
             & "789747779575879154006348694283918979140116946187140"
             & "502906882116635199821195647331290433501309423671869"
             & "092829953998007474826531467972203436569015691319301347";
    sint134 : constant string
            := "0.4024346508594184410825339335196296727887695860787"
             & "983330348426745285289846508941450120316989952994853"
             & "837895327016564286535127514480518951059237440075902"
             & "245842106030832411282219186977822783900801280308476"
             & "746722796104814638495852666804747142254798912177024462";
    sint135 : constant string
            := "0.4052413140049898709084813055050524665119477541050"
             & "890138058289893394311023351135448467297414812047645"
             & "030267707081423608221032123397377611618879087991399"
             & "983894792344617677629614466376849019774863447186538"
             & "954747370272130212455165122800285490742686257695427696";
    sint136 : constant string
            := "0.4080441628649786808207474989303193628209072906898"
             & "636067407340343582383381235642012728683126785912433"
             & "359263776235704395158977205530209034970854882544978"
             & "857114091466801361755163442019399600528923014703852"
             & "209766771460117404412578575798879569197188967404282019";
    sint137 : constant string
            := "0.4108431710579039421834666749289471041367041399440"
             & "958550780249699787653543310983875337898182315867885"
             & "892238589329360610218067491164625110717907138488936"
             & "579325325740741289205904108261764476510873955237029"
             & "720194129430061447222752680896551545290840786863633758";
    sint138 : constant string
            := "0.4136383122384345474719443235553664460012225907222"
             & "602846707167034147035999858249255474217540530440558"
             & "651808123523654794896620643981070336838076496120544"
             & "035731250277009130513828411321876400298482615966994"
             & "920235019383690227087718028012519668401363033049035833";
    sint139 : constant string
            := "0.4164295600976371825625989107894802468556074172857"
             & "931058987646494032866402575286675034298735441291868"
             & "556500209651269189808501225573440790067038653625822"
             & "746904457443792064646305114930179953400617795721743"
             & "349909719203889379220899201287564275575535053055385939";
    sint140 : constant string
            := "0.4192168883632239564330100199299511912473338408393"
             & "713432982946760067996254990238688798935469245422415"
             & "359713724822873518408583942948711200645465367464904"
             & "064252948009645921717941911700677685807702749667543"
             & "577265356994355026963044069161211641574090235751628356";
    sint141 : constant string
            := "0.4220002707997996859412879413320511405711334230154"
             & "070671722377080329788737291963477463296638801290977"
             & "153227400580753386771554776748238186498548336771731"
             & "859797164776772112908036152596369455518990453259893"
             & "448699896389337126688918550384873044674863975687925469";
    sint142 : constant string
            := "0.4247796812091088333572261892346645556946332007941"
             & "068757557751999991434480866013883200286352102348576"
             & "607470781492876080640501489970740378210211680063767"
             & "872820248799692557165337245611743558513121940167022"
             & "553179314381737924799902307268398703517694819514985240";
    sint143 : constant string
            := "0.4275550934302820943209668568887985343045786293424"
             & "586393648472046518438249098519379985739495457378037"
             & "624943453661015350997872508849588593069594248458435"
             & "132139985259080075688322788389569584037963740422721"
             & "118685781712369015757158328023873384247065371847091562";
    sint144 : constant string
            := "0.4303264813400826339081990305829805911285178100623"
             & "346043952892612779095579648710350617740521672820193"
             & "384191217738666278541222034291852637642052253186963"
             & "310268373635240834279824187493376776794277106834337"
             & "221682305601011238266142614131484582534277787924075459";
    sint145 : constant string
            := "0.4330938188531519684842226384895773616555246753722"
             & "248816024451496419591000693454963887390401225672107"
             & "089892364853991300388950510119087763080461828777372"
             & "742320119051197899763403952820938025075346468862119"
             & "050919506343864112236348952963176724278702185198664459";
    sint146 : constant string
            := "0.4358570799222554910325440803550755497283539552018"
             & "526669107410880240005586214739967150233125821154632"
             & "754489099702500383748546190017179556733962714043986"
             & "734532533830669665338114458734467198475792597777689"
             & "307991321765947068068055603663632149778841093204576347";
    sint147 : constant string
            := "0.4386162385385276376470257375461343598431014545442"
             & "370392629572239749453814826748644613989040384611131"
             & "123410196502588788254702626969625774843348408799712"
             & "335325981745033929053206324961705045728548293711636"
             & "752807456325607505868609791295015395105368686289744797";
    sint148 : constant string
            := "0.4413712687317166928799889682561164890337329042919"
             & "533093469963042270510760356721922852655307727922950"
             & "686614048011744896885023391417213295622280684976455"
             & "446675822392110109055837824034086275263684140934145"
             & "604781565838193293902148754449370023024922728556223521";
    sint149 : constant string
            := "0.4441221445704292316420694179834663691061023974011"
             & "178866841277090821315679521870047100800232390493591"
             & "703751676799014459898569203824455253477294847969891"
             & "989621441752255967377375802876351265330753027083228"
             & "682967725615819578091477558882135889822220332126777931";
    sint150 : constant string
            := "0.4468688401623741953530443887189265744266000023182"
             & "254650989703709025050177480703385323936744124241248"
             & "871376990986306008179682665447256732513069256941034"
             & "845626973889640708986635593768716834655186276754620"
             & "044405362007417377713940213280312646307043942612513313";
    sint151 : constant string
            := "0.4496113296546066000462945794242270758831870483778"
             & "572623655778812241904757362305784403432955948240634"
             & "275609799906262292838825302304420748030347222668022"
             & "335182490954444623502149504758398142374310098046185"
             & "678061400156689089748674578053536739836841105165726576";
    sint152 : constant string
            := "0.4523495872337708741330267029477754348290775509171"
             & "694613782517666006691883227550741759973290738878598"
             & "811441759203948643191538521683515034197150465402751"
             & "083870536842785755401594175076294390419761041930076"
             & "689613666360746633097999486294541163316935992744360067";
    sint153 : constant string
            := "0.4550835871263438235358692678967193602180057693617"
             & "736026727590793653638793805343698251691131176235850"
             & "745160203728291917576025159159642845407264129570349"
             & "101640040588517859354040305323271412921519400795514"
             & "768677362466518172070908305473789714088601776826359189";
    sint154 : constant string
            := "0.4578133035988772219049611553600256307008096136593"
             & "152933162216188529724108517564474185625293410506884"
             & "199943104202337441852908170850654256723535845547910"
             & "965158853452674405251744307690558887981091434341102"
             & "174146132170327011535054559337391420949442483366838564";
    sint155 : constant string
            := "0.4605387109582400236331814867414884243654576306456"
             & "521936829498960857226391413424386167986485136657837"
             & "965191702727966372051804085781579581756604719443496"
             & "506876040675983178255047799564873949817145373974385"
             & "047213386394513582221578037235243567611525069294564735";
    sint156 : constant string
            := "0.4632597835518601973907196370998114210798243014169"
             & "286606058491378684376946398398267872155880210157891"
             & "011550275174764651373480655709321564286006858887245"
             & "408322752907103572598395066054269020799558531511486"
             & "170543228528612683520608302082508806343855656497246816";
    sint157 : constant string
            := "0.4659764957679661779027560648877787034268950367328"
             & "910379702059564331413700435780236467855527735494041"
             & "939504615996446027026045789171115650755344424379050"
             & "581647999773565124523145997874107120995392996728417"
             & "228499705036849229541891784908196644373802876063314151";
    sint158 : constant string
            := "0.4686888220358279336976178702147423143397742323839"
             & "151459274730391300826704635523173813987940312770929"
             & "397456597377974871993462959538553651352664889663714"
             & "755280342843161215193834898431794129523854326879156"
             & "068734889683022890889738398024537792048042431663226350";
    sint159 : constant string
            := "0.4713967368259976485563876259052543776574603189324"
             & "806214016140310088352216651617529075631968879121177"
             & "453582548889551861368392179583452027771212751033145"
             & "034183892698733708244738941173784662611713889370792"
             & "350055203087704025713804866250608040022940812249952217";
    sint160 : constant string
            := "0.4741002146505500143985800146693645451235605768098"
             & "899666996257525821400268772851507418240100459495701"
             & "337335600033881275919380282470901269082362388249957"
             & "030095149476347508354117837534952977047495899187209"
             & "883056609308508953797740619829184994928954713124162798";
    sint161 : constant string
            := "0.4767992300633221333421581174135755721060696321798"
             & "402275257402669866139120279932334723702748725821837"
             & "837707479032569544807786828871882053164199959542031"
             & "437157737361882395574141347426711897728494371843813"
             & "267032971512556354638720063255408408594261631864673661";
    sint162 : constant string
            := "0.4794937576601530266798397976816903187312904532233"
             & "566297927930487177784233100846998650310625757221028"
             & "380438953111194389195607738922082295955649841715559"
             & "325906886094348040679945168942660371150175467554294"
             & "665815089724929786627891063564538692854834974712635817";
    sint163 : constant string
            := "0.4821837720791227485173444807974086112395892088444"
             & "544430747966272186590561952851698670976509806014910"
             & "801589907839273718328208480332619258661654673289417"
             & "017315950644779562469867555547792451308833403067189"
             & "447031989591762021507709013909084605606092778365175695";
    sint164 : constant string
            := "0.4848692480007911018229516986611111799490454472883"
             & "253494358003366183788613769070560940843463262055170"
             & "566439430812178022139872842741919720303066134013664"
             & "158024974392687930344147191375236474292158653570964"
             & "896822093150234002147136214183285519771829241957510069";
    sint165 : constant string
            := "0.4875501601484359546414850273076499075501098248882"
             & "845855719163327501235851994996111443061515185390876"
             & "195029731963854170254297528986071127060846021974998"
             & "782760713376650994668464058824598291199594545793621"
             & "255965739260237152283262228453638243638079901342878385";
    sint166 : constant string
            := "0.4902264832882911542295984490366087200183623387998"
             & "826273274744020450973640760878845770222081060328612"
             & "474265634926207172254336911739581012679931201281449"
             & "304598153513034096479309556729651051816554381649002"
             & "824635895370417715171395422080421002733956189564320207";
    sint167 : constant string
            := "0.4928981922297840368730266887588092682396873065483"
             & "635811114318744863825867427922036946116984077256090"
             & "878411750930004452956458840869839113436170005989757"
             & "227050317787342276241717671454824315694320067579508"
             & "055643049713148550326968122069775891245436282518828478";
    sint168 : constant string
            := "0.4955652618257725311502666695414869190434385500595"
             & "884044377154166428885268609119064802692431388651958"
             & "958193092953383819754425364895687921601249585092512"
             & "194940567665595255501872070549845181493152705222255"
             & "899225907423649065236666964171454777168894050321787867";
    sint169 : constant string
            := "0.4982276669727818524109838693598428202796319260431"
             & "717415486557398915976030519668069428635918823632385"
             & "358216565213621356169709805662784362369345688036144"
             & "825622401344078450445012372856316964607345893987148"
             & "215000465326784882796304941071719359896635652646825524";
    sint170 : constant string
            := "0.5008853826112407862412850037568218781730189009464"
             & "244529972233464239431259981189002467037347644511193"
             & "338840379929673124018057643141833273147882471650258"
             & "478758052979551107695997222465921481503090594690883"
             & "271209664743738958556672417907305545682581193125041633";
    sint171 : constant string
            := "0.5035383837257175586918670712604959834510947782836"
             & "330000052646722559809208656301981901541887020363397"
             & "998796998482700707057613025220886495807215450775645"
             & "703418784776522919812538536675964614069181195745646"
             & "242272667268465719163254698939054697104418832574760773";
    sint172 : constant string
            := "0.5061866453451552910489423435964958612549020365952"
             & "933947357664755626033671560566358284122682689221931"
             & "802143845120211865269589137355867202539744646628245"
             & "946889971836201866282797771554187541100353119711895"
             & "119285064752368376975671706649881097329467655387344956";
    sint173 : constant string
            := "0.5088301425431070369317493243516528643705928238912"
             & "837603929570591425500949957261576747316699101598933"
             & "636570325570106318731814468263146526485898860684194"
             & "788133222579388069534381393521540235531961655648178"
             & "335280019317000973394620473131687225237672585222236913";
    sint174 : constant string
            := "0.5114688504379703995043910009878180819645876215508"
             & "495736727383302681374994910764658317456346090814265"
             & "081485191413425074293606622319930602026446273847947"
             & "500874295666711144082442765196462018684281728217348"
             & "066674483971696482765748975607870307325342347059960800";
    sint175 : constant string
            := "0.5141027441932217265936938389688157726080491204162"
             & "178042879364539860379272742201183345765668289017940"
             & "589653702775218553661014349955692424614950660832991"
             & "151061170993839429150497385853704975242742537276905"
             & "021806796237038549042085670306980508563543761459364539";
    sint176 : constant string
            := "0.5167317990176498815087538760497663859531390303097"
             & "521104882035920141093319497283077194072311739602127"
             & "232125496375998029131489227838634555602176734112170"
             & "566302994375659012268393189303552613470210063908777"
             & "949381875792424764274072422843513394762225498964452979";
    sint177 : constant string
            := "0.5193559901655895873618299320920460363301387807596"
             & "576524517362128678295448363841925459578449457646464"
             & "884168655477615655466299466985684836967422975133929"
             & "250163253638842686449486296801501446562978335082860"
             & "036646851997798121781268353933474027633790963821683345";
    sint178 : constant string
            := "0.5219752929371543426942583175191106190747666027324"
             & "104558844314629004700707423659574541922168676975124"
             & "514908693267691421782355471703421827198712895225292"
             & "940516144557962048602579310769568932204041107166950"
             & "337886222131504397920176370220105940942050175488755483";
    sint179 : constant string
            := "0.5245896826784689062150984639335456725563206448119"
             & "373588726363875959195963160091506546822732900757611"
             & "815021897390506302939880801036215065561381218662489"
             & "214706401797045989943716337314522241144986857617175"
             & "589077182988899032866331941031069343552647212064804915";
    sint180 : constant string
            := "0.5271991347819013484642745754946705568645495625871"
             & "803817371024059366745883698908225064781859809257946"
             & "743177657178578892368195005972026898866614274425069"
             & "706096914051696429692785432799853349438333859700726"
             & "627157291273177839828787559010724482485677170878576250";
    sint181 : constant string
            := "0.5298036246862946682160546712352691063333056111012"
             & "520340517860983261173999865508431625968047726255931"
             & "363212118239425588202029961273572056067525659638490"
             & "642360362853246829596476243482104699358704216935954"
             & "138576680718778583912914165919907681190668787887239685";
    sint182 : constant string
            := "0.5324031278771979714428052182081180752112809608216"
             & "542138573406891345035927216559117170727126154960716"
             & "978085956499154509764892942151171976222857819121186"
             & "196694173939842747810048154654935395261240744377800"
             & "364034257054141835403528484733361824236446135155094740";
    sint183 : constant string
            := "0.5349976198870972106630769046370179155602656921719"
             & "002684556738510206831412598337502398630470672918217"
             & "663346744546653075097480586334339642430449559792901"
             & "419082865210969055593800980332925818751705832832676"
             & "792167476072342135029664150632155200153362700023351120";
    sint184 : constant string
            := "0.5375870762956454825022149323489381597209519052688"
             & "873708295028039411929637002733699816284603309066279"
             & "852551675360607949876364572790988846475571376965425"
             & "482212858793757847687686320081665195605362722955570"
             & "077214139973241725419087110121672118651162156208991175";
    sint185 : constant string
            := "0.5401714727298928812978454797368391243804773726575"
             & "731328559542076742614251127228693834535014490957262"
             & "471044685523369104702249940763652938638851656514379"
             & "030760393480270717346430701295401401910048329659000"
             & "489098632706235613911151179800661559479344886821055008";
    sint186 : constant string
            := "0.5427507848645159065867686612074833318665104452394"
             & "303870138832970441219322200386415236444326872372056"
             & "700616651039546977208778873417161825165819232663294"
             & "367409788581250851127787264554071339885656300629946"
             & "699297599587686619225376963501243266042690003968630964";
    sint187 : constant string
            := "0.5453249884220464223139873471738261971863535360472"
             & "649722480657047027389676435677651676327386251232354"
             & "422477892582883818991882337412619549802139889527022"
             & "938112208240989211709259534278734439861756380913391"
             & "520267740411356969192838288377746676534726228272692677";
    sint188 : constant string
            := "0.5478940591731001656088205706344303158431835263620"
             & "271505432991519430723608866511100780725112451994454"
             & "175155309410869742519474891911459633158184249777712"
             & "159135720207644178564010278866853944064105119598696"
             & "290435900991550409283117149553680236797453871785147159";
    sint189 : constant string
            := "0.5504579729366048029772898925285391908407957993772"
             & "653460290228409481027739716884170783773450786742276"
             & "397881710743384297383834289429705741080609359850486"
             & "732460112923330928092532902775916387040113386734291"
             & "137094241340681688661709921886736102797091142937368680";
    sint190 : constant string
            := "0.5530167055800275317642269884598102702774565804435"
             & "897678269512669499655247677167205902773101266730008"
             & "728096714408388038337237059491298958630518071912351"
             & "585840827591703202909411003938942366228786122187241"
             & "019872620231111177451203056632533013738605226221327105";
    sint191 : constant string
            := "0.5555702330196022247428308139485328743749371907548"
             & "040459241535282029492475774800683831248803512639343"
             & "555379458249282081463006489743886788062649174003545"
             & "193643593342058942115829145001516024727725570039393"
             & "515637931802982156762049658632383361088879709963540299";
    sint192 : constant string
            := "0.5581185312205561156937029638155638457395303538738"
             & "867342363454889817711309235287972510008872208141190"
             & "248685098635056244014106178472174392967235839413989"
             & "520450995151362141976244092439824757299999880752446"
             & "900609127212240842991367753106626745858513944225231966";
    sint193 : constant string
            := "0.5606615761973360238397102231455323464061827104244"
             & "866484907053386673352487699305395072211804866093502"
             & "257934635032183827706624038579211404216282607364765"
             & "265567095910615325477140211652018982184326452927678"
             & "054085029523081115432213873289629725411691429092580280";
    sint194 : constant string
            := "0.5631993440138341150073637718570158237455984152604"
             & "557338053465946957573353776284664373062059806169926"
             & "968300909939060488541433361950079519774928295983456"
             & "150116970187661993521320457627559607800082579384435"
             & "027546774502899476567190118737881159689931197525864617";
    sint195 : constant string
            := "0.5657318107836131973897650113692660904743014432468"
             & "325846677441006973264595406666136283809542103063533"
             & "494605064479803995331102355444455393368970326780370"
             & "817631851369515515133578406137128091077579872417618"
             & "584996424949954210730099689344805564066875463961666075";
    sint196 : constant string
            := "0.5682589526701315497905484891559202039443500356891"
             & "325755674634753283773618214707633616938651378335280"
             & "668918869772038592305594758396485776304332307982795"
             & "874142095666581244031610331666898933614806404932945"
             & "612089306108173869059379963318637993953146159433570074";
    sint197 : constant string
            := "0.5707807458869672802326528638849623650065504461592"
             & "043109527031644842099082840634183914531040168660657"
             & "726878317569441340161492007543507695867138504065662"
             & "072735513110742565513266925714824913128859776895375"
             & "255430919968575654529812822671132087443081236968254896";
    sint198 : constant string
            := "0.5732971666980422128201712389421399336043776624972"
             & "067178580572285003471764163548099061656682637589078"
             & "026644957261023230214086707748360486449313672393775"
             & "610266725413943045432437550252673402876842583878230"
             & "404133415821431300936154551795410736884829998363648926";
    sint199 : constant string
            := "0.5758081914178453007459724538157308417760084553140"
             & "966022029889262442367885227670337478439943268291258"
             & "751515043273020343810366353419781830506440676659671"
             & "518980652118464684098943308575471213958907573348384"
             & "304353269649828304180965699183686129957249869107781214";
    sint200 : constant string
            := "0.5783137964116555633422450192905771729677952159265"
             & "018988362914963433729100905933673967337558978004133"
             & "322020887956465265115479203698979690676709262350418"
             & "450910827692515362800019470317592992555761326611004"
             & "462516831134521457296146822538284295662860431386886326";
    sint201 : constant string
            := "0.5808139580957645450755952716785138763990615817533"
             & "023272667800751035973385574586595142774992220216564"
             & "099392657515904338935539591361775870307427383175907"
             & "605427131505555344294609408209815186491976336224896"
             & "227649924404260924713154127030120785046573304778948483";
    sint202 : constant string
            := "0.5833086529376982943928309612343084453657295517112"
             & "574064951557779390422051524371047513406771389330365"
             & "941422055182955930529572558471403432675094486967404"
             & "858145504975581535868212595952390301118766459505659"
             & "785658688553522950900985226956531143458022926432008657";
    sint203 : constant string
            := "0.5857978574564388603280808381186622857671368129040"
             & "376588390429865154945995197178191115892068783134378"
             & "389897809831485312898895908113774202896765076021187"
             & "448111199234944620361906441658850060243261152439826"
             & "787730011932446472467457326880890213910881341379006782";
    sint204 : constant string
            := "0.5882815482226453047864398132348808771714497340642"
             & "036716783410379977402712991684314747423022365860932"
             & "521524553563599245086591297387111149268132473028145"
             & "909790636342620474829218395011193956380776772131019"
             & "152116780512862976072416762389213054402623961372167567";
    sint205 : constant string
            := "0.5907597018588742284238879082605695714561338505893"
             & "823218551275930807683842486667611678670213622261462"
             & "234644164803768142010447659695958077037241339037850"
             & "653884091718960108379790080239776847859087575063913"
             & "787465856155341366600245763218391084071548349098382588";
    sint206 : constant string
            := "0.5932322950397998080478094263125274426679727983917"
             & "349942720027926289085759972243915559112345549981837"
             & "045587381074719924616012942706336918588360541129214"
             & "599293559016764362392533709451715201679182866398427"
             & "711953289421707834275277484512436193189159648518445930";
    sint207 : constant string
            := "0.5956993044924333434670365288299698895119263384375"
             & "047868977947345129738936287256664315349110825853769"
             & "003100949115109402206769593983646369984745112052270"
             & "470367907348671295556941207881832549187246403724741"
             & "835730495239655415134415424512611437868624111730332172";
    sint208 : constant string
            := "0.5981607069963423117249586521624986237340592960663"
             & "008127624009663732205155437014182821138648550182198"
             & "041386214761290188460137129959174105597437081328772"
             & "865829961995001177833490987094757323884563441498802"
             & "938757436937358224124095287279421091968688540670192431";
    sint209 : constant string
            := "0.6006164793838689266538758955455595386439545513218"
             & "809953821578919554504540907427737751843393314327215"
             & "945265858934856117472240085206193594067951699686182"
             & "455088949556945850637143749003751259185312054228079"
             & "465873763157372328097313597703440037905129006963253280";
    sint210 : constant string
            := "0.6030665985403482016934306169951154427009111934547"
             & "324024945770349033566458640528168599835260756119262"
             & "962517658281010384013079282803380722831364569429293"
             & "098700971870807368635754704572063442076733353032595"
             & "711399835491027611565213522685257944237661004057844782";
    sint211 : constant string
            := "0.6055110414043255139206269413298796581825983175372"
             & "228999106201976322202011281189863287846441187291491"
             & "144768902529327548649328942384523869387978363799601"
             & "326420206321318552992602076963422942744042909766281"
             & "377707533003594793931758044105307323755737802608745676";
    sint212 : constant string
            := "0.6079497849677736672436426710264256112166043249865"
             & "144038810334510795077714632028975622717157298095742"
             & "720056800440952301363401509647756614414322049635281"
             & "090647914214760476283334776740283608775970986673250"
             & "406950308713970963612042656209758248343741368050020396";
    sint213 : constant string
            := "0.6103828062763094527163521517406882689736974200427"
             & "536464451341914606160555843296894896186928517017990"
             & "988353693387566641824121223213100613749723080635504"
             & "048118257497239988307314607756114861336682523249423"
             & "266386044574523065062969322859783709145541298872684745";
    sint214 : constant string
            := "0.6128100824294097039352119357182669348358081983250"
             & "578624528845973159804129812242222699436589255619098"
             & "191725929940859932847990147660344848133892125946602"
             & "452859016790556720787607880361332452926594290798663"
             & "012021582535632245218513587682230709501227974275751451";
    sint215 : constant string
            := "0.6152315905806268454849135634139842776594300077644"
             & "242274544609298006215692903192269110107483029334541"
             & "586100695585596470941202999304458162281938483573217"
             & "527500616577186174695038742434129133051767196053514"
             & "168041248610859927723311688766616308764435158399419648";
    sint216 : constant string
            := "0.6176473079378039324039794017162291302857021824689"
             & "483307698284270294303693774715731033102023654274543"
             & "958319838812736584804598084988454276662950315835092"
             & "285946962374877471153633831381647745284656863252111"
             & "166655078241950561697507952109187747741157666387802576";
    sint217 : constant string
            := "0.6200572117632891786462681913114239630541487502709"
             & "073619040197912932014748016069110565428674319983132"
             & "222097438528919372298218860694966927203210985132345"
             & "273255632695862848091540478405500303929726336609229"
             & "400872534045475327351747691387202379625193104514419706";
    sint218 : constant string
            := "0.6224612793741499725191667208364860040790529154497"
             & "645824751449366331572160239533105771395867911164168"
             & "771250394018974719867934245588820067335320752002583"
             & "734646553708278987140108083542825027959000990472998"
             & "868981674842172836735398150935592194352787414427934457";
    sint219 : constant string
            := "0.6248594881423863770840728162810527125578830486090"
             & "799052802523184906598474771427390155904186015826647"
             & "617026029911673514394811774103683191113849024921186"
             & "484219358340330388365665141806136399786657344943181"
             & "079919272854696923131788424796975487956856513656739446";
    sint220 : constant string
            := "0.6272518154951441135096225651662877879063712406198"
             & "640212778253309966804512228971185442291205552774555"
             & "657881279188778552260090017495531390796834412686163"
             & "033932350385273609160670981196346173074392269374636"
             & "645172163566830729080189018235252229712444700335019541";
    sint221 : constant string
            := "0.6296382389149270253729813407145820848289792322209"
             & "340930839029102992388619721426669229583297251365555"
             & "110948908165368339932360188977238321854255769101683"
             & "375579667115625566831339058780125686856190650010861"
             & "443422927114772251029283226028776462067637535915806317";
    sint222 : constant string
            := "0.6320187359398090219094037057276920411962095885686"
             & "819279506708617217019092134143227160512764879393447"
             & "329236693861329000262629370333807299618185558429748"
             & "937116183620516991046551511151297421319352544607353"
             & "849675633502861077073561821051857851000341209931205531";
    sint223 : constant string
            := "0.6343932841636454982151716132254933706756870948417"
             & "216064338247186966867261235483941717674035365368534"
             & "697200773200211345837390831366418354602759052326453"
             & "092637223414219185238996320159304210503675542262632"
             & "577892213579543460196716484683219890847267310021826479";
    sint224 : constant string
            := "0.6367618612362842304139434349020790264832412203675"
             & "084467069993495150241594271522426303585879571135169"
             & "367191294864066920773244889594834911217035565775619"
             & "559417372046730572873224138984821301204526738304106"
             & "564387416141977149697204421632859576901242243696734130";
    sint225 : constant string
            := "0.6391244448637757438014881927921738391132953281173"
             & "675848149560865914244967504258663418668696043187461"
             & "478702621533231863461467143009762610058800425620851"
             & "866674244663954215279280957686199967508579382104242"
             & "675129169129752610136795047801708551726713158947941608";
    sint226 : constant string
            := "0.6414810128085831519887398976942527869152617376142"
             & "429124131192750126848703077404111091027113814630904"
             & "504147014763687655094282760651716129085470366369040"
             & "415366464588149723493505349908557131627347668502136"
             & "950880577239571822238926617365301307457790254559293323";
    sint227 : constant string
            := "0.6438315428897914650680860631769553884419570576009"
             & "011868107220178840894113968623573184042765517183682"
             & "450448352411396657497562775897490129876494644964000"
             & "816913186465307245522761472846640093925969921073409"
             & "334841622732148712970263614558542403811869283729383512";
    sint228 : constant string
            := "0.6461760129833163648328022195365852888821737660677"
             & "396609185129137581387590635623783601916288482155253"
             & "193726184521766962974259176327270335371302172596500"
             & "582465130724844714261632131154164953046443971028046"
             & "198947396254655357078476056016656276496908966410524969";
    sint229 : constant string
            := "0.6485144010221124450845605508348932123722212627783"
             & "756061610926447170148487295631545935086705193037931"
             & "445000007908909205541130526176559198014173621508160"
             & "763533858131852650696150971392066442203635401213257"
             & "150891889131225512747679024598541361987495760843166401";
    sint230 : constant string
            := "0.6508466849963809150689755729126486795480597434835"
             & "441866516678657682140513317415725453050144101835682"
             & "623652040969979332109209195256277297234732547788926"
             & "936227541678402409592435335487432111517261634135605"
             & "258772376569517402685339061038266691174023298453219584";
    sint231 : constant string
            := "0.6531728429537767640842030136563054150768600237141"
             & "605475133036582798174067452799172122116779246648060"
             & "015260864579856859481800978000664178278943480117592"
             & "094318091244506544871566473699776328446126145486567"
             & "082375734509356173714763937306477404600040560974473385";
    sint232 : constant string
            := "0.6554928529996153853126797012293059676952069517977"
             & "622662051420064331176687319945239297910626840576399"
             & "522864201321381812401887460807318626824450371896659"
             & "145720645914916792228787665626071250554374520629996"
             & "109794470096138696492690625254558759383245997618072531";
    sint233 : constant string
            := "0.6578066932970786569311822637300037965084167495920"
             & "388765198309430201853225379130424029761764228519113"
             & "297327399596575881098084770108719084505256897724373"
             & "201442274659680856258780304451194133207026716594613"
             & "030023374309323688007018436695040995738297132189493157";
    sint234 : constant string
            := "0.6601143420674204785594907468958086011814235036758"
             & "911913415855637633579784223784688232368289228594303"
             & "332696685255437464966910094400861663964018278153058"
             & "263345643033741900216479725936134665438089942956338"
             & "583132177591265904589403467944643439105101091670060971";
    sint235 : constant string
            := "0.6624157775901717611130698169566881852974425440511"
             & "929825277633930492309617675878204032686382703335051"
             & "241103811789726622294901122683449046583346924057913"
             & "681758342689730173734040003427463561506231192936307"
             & "456296219447900470866974514893738279957738081872698453";
    sint236 : constant string
            := "0.6647109782033448681303249852974479320666974744887"
             & "784859826918528432748115336806696568784741934041869"
             & "177852012879568406387837201809358295485207271013615"
             & "879988540998111761644158035677984258873407390539519"
             & "553865049606529721564629930722444012857646551708442196";
    sint237 : constant string
            := "0.6669999223036375066501542217927266025220265822372"
             & "297847236967136134419351497348575039360969155879917"
             & "812548695728478079612080433103207948244115441275757"
             & "087370137579270255685458854222436746901956207080112"
             & "069124040590963214925059191156365379645373211466982338";
    sint238 : constant string
            := "0.6692825883466360657206963659359292635496952150734"
             & "586104944256345583445023027408017425936581339185043"
             & "079886948720985820335257903354418306733839520813458"
             & "486559794023865522971390482558506574139123857254908"
             & "315763880531298253418309435708939673066420666750498239";
    sint239 : constant string
            := "0.6715589548470184006253768504274218032287506321997"
             & "944988832130731094649118790014151672904931784041606"
             & "634999442618036977416712557804483405115265096832534"
             & "678659042541388268469203423531082217931279143995827"
             & "061540034535933759222130773492608535126292265383202668";
    sint240 : constant string
            := "0.6738290003787560609175683717822752457223626629661"
             & "524063954096971439681915157772764276076471916838191"
             & "306576818726804259491470871459745591261248843121908"
             & "795844543470495001521135960444823170831061510925552"
             & "547707274457457272781930350291184961806226051385736385";
    sint241 : constant string
            := "0.6760927035753159603604192276581525580927731452766"
             & "324608588630133987883403023735467098779767622733026"
             & "503399383041087631843176711454740342176588656697748"
             & "384653027133932164195200434955367179346253803712206"
             & "895830936634935742856018638390500556381228255353936750";
    sint242 : constant string
            := "0.6783500431298614868736550417149607364422033657653"
             & "815716528878847599683946010960859016307172864331593"
             & "347954459439370447413089591940172289229497118424322"
             & "988906249873347043812335270179427238003282568815807"
             & "922145704946523246581753180406960462404916912072473996";
    sint243 : constant string
            := "0.6806009977954530505944304644563996185348546094786"
             & "518119251422058873009963578370745862507508405292100"
             & "544647750739483472645850368808805793462255083137411"
             & "871602170723507378613046224557190211165732744514752"
             & "585877582692122106396056964844805327380883425634413345";
    sint244 : constant string
            := "0.6828455463852480681645961230581124832961157832365"
             & "145886408415108394666365972970657536333661388673231"
             & "852993856399362785345121021406942123697647586624007"
             & "859535276361773049301207600325174918123328996179420"
             & "882845488262639957525881647556621620579026040006627877";
    sint245 : constant string
            := "0.6850836677727003813620525448786688535463171180026"
             & "453736937486985984569486549936815874603508114767262"
             & "435249031549368866246419210743354969037001314237583"
             & "346614301866289686450988210029351890385932382332207"
             & "470173161163489531370965196974484622505339405170996985";
    sint246 : constant string
            := "0.6873153408917591081991869482317427466776116545951"
             & "892451570873902027449124012540074630637094460104829"
             & "544116298456885928283651033453289770738663215570271"
             & "324788603729065524979495516197730986809501172994912"
             & "329232510953329846606484263936064506063901599350195192";
    sint247 : constant string
            := "0.6895405447370669246167306299574847028455368442791"
             & "232258618203924978160457288778645751098290599853791"
             & "044920716402292097131514817222691129798460371270071"
             & "349750121536438396914865545756745628251843769788127"
             & "807239469886936247712890010017031555161552739003932576";
    sint248 : constant string
            := "0.6917592583641577749067341320888287838269067705299"
             & "095911102144465317252155279438740777303584453751094"
             & "186114454525661629646869871671607742040001411735433"
             & "013327733300460427538586219157311474541222161866757"
             & "796223106010909229275848776578160480656388054658992206";
    sint249 : constant string
            := "0.6939714608896540090037343890191158870916147727358"
             & "502188425787927708104317601256810398247547895826050"
             & "658317959047665610448728670515686964997842841618985"
             & "434147995815020117494850961856939556583802252345792"
             & "593578658380341593630576621113695593227872172734145871";
    sint250 : constant string
            := "0.6961771314914629447885825914304004643727746960111"
             & "956668567824901156531946713995050976716898204124082"
             & "610358261379267682119479805636575235936989655810247"
             & "788728923985172219474911976624970265155502571026673"
             & "022614712686166259017377817725161138539928282850986813";
    sint251 : constant string
            := "0.6983762494089728535548135030617225398363276188521"
             & "383524865554105532064124710882403167866541746350997"
             & "289509601103265248176422468895195282044082833735633"
             & "424012052748282162576554149438175344271320497336832"
             & "090033395537495887781753442343432059670124377469781450";
    sint252 : constant string
            := "0.7005687939432483667928663802436673373586637428573"
             & "520206839033491748725997381000686901569091625534354"
             & "810081059523545103522952518474088354656031641208472"
             & "320580851209949653356687626767249607426138827591552"
             & "977454342769824150862106587417938449512481339307064463";
    sint253 : constant string
            := "0.7027547444572253024529144208959899415321728118683"
             & "075951906864810262181316229355224222327794546693269"
             & "444181399226491945322533783088388505865779654358104"
             & "480808121888079402499898415487632077507162772266015"
             & "317228581564582252236911916593183006538963128343486718";
    sint254 : constant string
            := "0.7049340803759049088525237581118148771392938248127"
             & "094351776507713929956213381924396801086866010559808"
             & "095496744752437725362683727308516154922272269289357"
             & "978123894042094293514844424171759663968203183186236"
             & "806275113836634836283206296969661458878881366254579763";
    sint255 : constant string
            := "0.7071067811865475244008443621048490392848359376884"
             & "740365883398689953662392310535194251937671638207863"
             & "675069231154561485124624180279253686063220607485499"
             & "679157066113329637527963778999752505763910302857350"
             & "547799858029851372672984310073642587093204445993047762";

    cost000 : constant string
            := "0.9999952938095761715115801257001198995529876336221"
             & "876541107369867438233011131433523910096238281304075"
             & "117980405044947912524999914655753762623999310202055"
             & "873148096834424091648034152634030721337450075655239"
             & "097995483576270266099534962916496186967123651448640502";
    cost001 : constant string
            := "0.9999811752826011426569904377285677161739172509443"
             & "350919401576950808210694765900655445062260270025721"
             & "684823688320219869651080170042578869600381875222785"
             & "540644987527990380077804258293097521190931630170124"
             & "610298627236283023853431773483415232459094529154950811";
    cost002 : constant string
            := "0.9999576445519638663331209195316303717368608499531"
             & "550935450877445651128811655924872174346546444639073"
             & "767111893453169630997624575282367272537767373291883"
             & "662394375730231532558972573382455406717034670717053"
             & "446087566321159820914512699623749106601381772140292336";
    cost003 : constant string
            := "0.9999247018391445409216464911963832243506064688022"
             & "178302440836730014837213868841334793501741589368376"
             & "047625856712191820356825849276874109509686173466167"
             & "182467053415648433899821190990329806340071793711278"
             & "921743149703858819234363366599750256528858838242978275";
    cost004 : constant string
            := "0.9998823474542125256330496265059156608292573548366"
             & "226134361084596641443980967303008001314853916560931"
             & "653349438935413501465519468406015071132852475665060"
             & "715419680146549357698311269484506073045554247388405"
             & "336495526932037742157105862362282423986724308662684826";
    cost005 : constant string
            := "0.9998305817958234220157222749226655145858957850267"
             & "552231789803669895796702500070420021220880128273906"
             & "036424979212831386200654845653655334494162240918881"
             & "848823504105492343208578125549620468067757749227722"
             & "760786437562316125063336113327466948984714472513810019";
    cost006 : constant string
            := "0.9997694053512153216576170363329930430512377077464"
             & "978564183800381509957301124499121640493646798982950"
             & "476495622091690884331858233452567964065700794539002"
             & "409759436861167339891795612930121810149531567998508"
             & "135987059041410933600340688145626674065845306221800631";
    cost007 : constant string
            := "0.9996988186962042201157656496661721968500610812577"
             & "296246441237879263816687547147750419016055276630732"
             & "636855255741913959553359622019704652270908328616194"
             & "184995415334272627893775125468738404904727378998692"
             & "471103942599802968164076655624156291674831232035803879";
    cost008 : constant string
            := "0.9996188224951785971168306373539819130463131387659"
             & "250077173328177194347814975904941807779361241973360"
             & "657195271241146753220649392753923961106698931614077"
             & "660145573667779461868840006816948195934107266487961"
             & "704192471609038289719261383044167175393693504504955655";
    cost009 : constant string
            := "0.9995294175010931630797033221567409693597944190919"
             & "770211776509087587144040740069696910600031388104767"
             & "199499128371718432551783033492858286378624638537230"
             & "492895163566013474129642671887642097397546208572369"
             & "046717452581083505530858486157153778995041571558056398";
    cost010 : constant string
            := "0.9994306045554617720190083272854685356212209690526"
             & "418374968986454828960440596227698215137450856146726"
             & "631144136593173506955823653212055921772406492571394"
             & "239591050494073275955746642217924625384768522218438"
             & "607331874200800476765273673342919076801934809974878771";
    cost011 : constant string
            := "0.9993223845883495008962210111399103525911368291690"
             & "488414700783462084728541753948606086943047741244546"
             & "247752525150425402250406721323980636886260726142064"
             & "423734164746480005762557448614064486410551870548347"
             & "478871891218134439996444191940357345404537865022998278";
    cost012 : constant string
            := "0.9992047586183638954929500005057034450442347068512"
             & "060424350343548518605611270845924963148470537207840"
             & "039198472117990071370362009854753445537643697419082"
             & "171574745052028733742661565037311881605877507818391"
             & "390488578860442455492733581956340856161661283890311510";
    cost013 : constant string
            := "0.9990777277526453828887819968641261429898648213239"
             & "274790163487757353577643585245814293732138210092255"
             & "607248773365678627962252735812865240556147732684878"
             & "415527759962678155879362023948295858048294089648675"
             & "318658331890643647009958424560797400208276378167032493";
    cost014 : constant string
            := "0.9989412931868568506339302657244654427194216613726"
             & "419750560313946908062727733108435230619910243468558"
             & "963144016761286257850713734910312328796958730775475"
             & "544884985364003473676381459278059348244403419544826"
             & "050562159267428186781689375982410242720065481565546654";
    cost015 : constant string
            := "0.9987954562051723927147716047591006944432036147046"
             & "117943428170736856099059521574038311892442663941710"
             & "411642753012920518949380253048468143861730196515371"
             & "856672868503203328432010716766830495033881033717828"
             & "249563062601422232746759056718718093906839410048508583";
    cost016 : constant string
            := "0.9986402181802652224181990491800163114016797180866"
             & "067141858045441497622818959219527227366797691180266"
             & "339608635264481060322820131985631270884937699434764"
             & "617462383395208255924064741500835698072227919843066"
             & "731206524376774053990045106421677803715717595768269319";
    cost017 : constant string
            := "0.9984755805732947522085590384263713301081889047769"
             & "385719273966753791254822045208650361164603629810427"
             & "519787474875849291463751931448660581845465112182671"
             & "794624607126475797363262090602864229170937363151073"
             & "526366044974240121541572854080998411274394899270792219";
    cost018 : constant string
            := "0.9983015449338928407387821630291092908158931696678"
             & "409510585852725396188756408788652597022456553249756"
             & "536259139950325342203317882486682318273607928146503"
             & "559777876647898383840831972751423987296001079063791"
             & "731660730153521175874736021576357007461460430176397868";
    cost019 : constant string
            := "0.9981181129001492071251558606836062321140632113778"
             & "141992456004556118395173713637433388786880504218344"
             & "134758266407544204991274602293829914522218134807354"
             & "213256770848862236377214171258450254579586390344656"
             & "575613611960009372807594059860454452379706945856233350";
    cost020 : constant string
            := "0.9979252861985960126230254623090194905748672013469"
             & "828081331772000341233114162448742677642034971763837"
             & "781191341880924709057161639297482352084193960696642"
             & "959641129983330214011366861112080784850473350600602"
             & "151416367056510545361626769590979478715888488391869850";
    cost021 : constant string
            := "0.9977230666441916098485467284287755605433413968690"
             & "056141311132658119507200350651581780478796198078540"
             & "327066797917348435717121668413859680351119580935073"
             & "430527355997192778816022107929871559517157323290397"
             & "148328449239019681318152708271062422806024149705125672";
    cost022 : constant string
            := "0.9975114561403034596994483898081441976323805321324"
             & "117935165395169731353813446088763183372462016758985"
             & "638536590586494222438991225955821495722159970377387"
             & "108365192229482466256550835115907020669958592854466"
             & "509366703085357418688916656694225325008539591684354140";
    cost023 : constant string
            := "0.9972904566786902161355971401825678211716886791662"
             & "213732075785240505426808841642150637983837742537777"
             & "780594907648121931286592117752244348771536463150986"
             & "839597913789087926489993079854622590350344157417378"
             & "305485798202634277151208160488580473381227273598949141";
    cost024 : constant string
            := "0.9970600703394829789879899493683807050365624334870"
             & "277165261250437243719240791124805542147053632904292"
             & "330599730055985331381188345307585010791733780719731"
             & "254139982880836470454161425253929894455478597081041"
             & "227319588477862832535065769621199249992364549595593395";
    cost025 : constant string
            := "0.9968202992911657149726293983440371335985977721848"
             & "728384982186920648330428731597999991292167879077734"
             & "762468168290967488196749554676822760362949288937862"
             & "482848974020548400755359313985479067411282542218610"
             & "388965228808206623484909142878102289314521243864972529";
    cost026 : constant string
            := "0.9965711457905548470935669103181862989689003194724"
             & "270469534710706320546528330397481108548428793895736"
             & "022311566723065236243782258042843656218626337656784"
             & "677728795848977366963095610165274861702765387087602"
             & "591296487629135089653446396775912696993714341155449077";
    cost027 : constant string
            := "0.9963126121827780126272261896697316108763234867490"
             & "231742314941322048183052286663034324352828581139530"
             & "266115339499361137423672006257261462429436424245357"
             & "605214765526242324632646393180576402620636976654891"
             & "572764648207738299904366678066248481688962372069837809";
    cost028 : constant string
            := "0.9960447009012519898879448102795669989215824446049"
             & "551206694524513484560542172664028865954502540119396"
             & "265430745718330229660482999919993521378803476796341"
             & "766622444052193559804790467459162643851314480568208"
             & "259396379651638477608993001806997768603783988854982072";
    cost029 : constant string
            := "0.9957674144676597939824956425161862154313075283397"
             & "895108131787511431958335863805779519167633221522049"
             & "342421162073319386192616725265370054096093868393771"
             & "355195867457189410042724626864818781419618015441754"
             & "777315951253576220344172132652227825830685211589232291";
    cost030 : constant string
            := "0.9954807554919269417691716003477196872654293223779"
             & "239383944372761458421190178183284205502344645098952"
             & "402737971598364364472531522370316540935872522019737"
             & "441457535709986641579237030347476737405646790075821"
             & "288722087082764982648052473943087111417461708257599315";
    cost031 : constant string
            := "0.9951847266721968862448369531094799215754748687298"
             & "570618336129657848901668945865379725290842696483902"
             & "877244931182981983461417865099877607189105687761328"
             & "478748617236574343419497441643958171830872294082617"
             & "276700525377249754609992534574515467910511666484225456";
    cost032 : constant string
            := "0.9948793307948056205911661067562028172051981133855"
             & "192006176197000846235633873621879262700855760064836"
             & "581561971022653583756017657107177477441949365146222"
             & "698659706194638029091968687421917383091284445403275"
             & "943867318731205350230046375347716821128864827206736373";
    cost033 : constant string
            := "0.9945645707342554521191062433890623879591963997686"
             & "689929497134733049764975365965519646053125923923732"
             & "021389239686013965974093242507455858116483971953831"
             & "918448215913705899331015679594927428165464877077856"
             & "079615531583643926986453722330291133952759373461538368";
    cost034 : constant string
            := "0.9942404494531879463584134419068988176869989775529"
             & "323156918600807582181213225653677072673736343168314"
             & "268617930921506414506004743990880827454260285572680"
             & "900146734233491110616814739701981924492973841528474"
             & "783674119869794505399399160302440309713984109668790677";
    cost035 : constant string
            := "0.9939069700023560415469228132477998214355953506226"
             & "764763674651868930403405085985075889181078039373099"
             & "271825577816977486659233044089260044331143156447598"
             & "151510259498422159584724731598000671641607656862507"
             & "899261936589386681050206122736891296944969768765531038";
    cost036 : constant string
            := "0.9935641355205953337820216973419474904243639784598"
             & "839099021476841949370539161558293602154930940078513"
             & "184895026119051905690968865168891844911739474676477"
             & "637613462198234566440375509936895997700726869520838"
             & "107987283598348327904086179535767412346016393971895865";
    cost037 : constant string
            := "0.9932119492347945331046010120927827784469943659755"
             & "387491148048822176047274351681217436733955647559801"
             & "150203577147384792449415431093005934236920592093259"
             & "958116388988312667936118291252681596797368647371468"
             & "914192330984992452745966948446718300184854592701480191";
    cost038 : constant string
            := "0.9928504144598650907935633439675960681742676607802"
             & "046356346969945330992396221153599285142673680410430"
             & "885893914298543467143410996471326024756859037727193"
             & "284370105011876364761983977034664283038182743687529"
             & "555147228912495682156221889628653368629094327598042436";
    cost039 : constant string
            := "0.9924795345987099981567672516611178200108206546341"
             & "595464907092780513666698265584048621348108183070855"
             & "256895261610517163417141668727170182189760817570592"
             & "799036145772461882840209092017960980535432014815041"
             & "276253890927639376824527836780717454502602506630611179";
    cost040 : constant string
            := "0.9920993131421917571120854453716869282950854130243"
             & "223527003934927291294853563690001933717867323780960"
             & "189967829459768651360246884263104785409447201970977"
             & "929123197415893603351630293655008832525458900272145"
             & "935754746824485143400908815640345506056228672401606169";
    cost041 : constant string
            := "0.9917097536690995228600499310995717595784264084693"
             & "254023626144153838723120274582293706301438463523444"
             & "900482292754535282253009148337231907881817928311956"
             & "007473742433450199040941204791312641388401707998778"
             & "133031046058857985526959706087265602588600155884442961";
    cost042 : constant string
            := "0.9913108598461154189573497986674833550034585194173"
             & "984391230798458510145478268668717638783054219086403"
             & "442548241938425462697343938183241390617203660937062"
             & "296338470462814423413748520035173718930951026760699"
             & "862384315605132622263898881284034658068604829205103577";
    cost043 : constant string
            := "0.9909026354277800251082370105274337521976237991149"
             & "927256992081100643996553792994128615689906479927039"
             & "920767795591021692784317564445054405244349341666236"
             & "043564363098306066771731501327765160757277874184645"
             & "117836756944083874375890731846253270449254182023331008";
    cost044 : constant string
            := "0.9904850842564570379986822425364353778167451152296"
             & "616438562948984225063638873064009260991606615290055"
             & "476237967462416992014053166584334804817442806539683"
             & "980918594224461313610952392276982657439256949637998"
             & "256509641554444272682284528667622311909845612690198898";
    cost045 : constant string
            := "0.9900582102622971055059064644647793938596186386186"
             & "491597796296552743739230798941892049947506658785750"
             & "631376064947051335215042077090071685325072867424775"
             & "949377154690798529550564273015807806698110938070764"
             & "442324507450652543854595160340391750494939856300840431";
    cost046 : constant string
            := "0.9896220174632008346236944537822198667249656540092"
             & "215756965962402622914860215558207185344718543851355"
             & "413399927592972581569434239086407081160444424448774"
             & "148610748963650561946723857547078076131645540241546"
             & "463589864285185459109780929752622514168509050611444312";
    cost047 : constant string
            := "0.9891765099647809734516737380162430639836895333369"
             & "074010191503644849872398622341144003243385877090478"
             & "889318178176746209443665298299387621976957746370975"
             & "252269809696531388216307934022873523382788238945821"
             & "141332294581887866077168433397519930525255019687180708";
    cost048 : constant string
            := "0.9887216919603237676045164854897910426445959783492"
             & "790177046789209879760118680682887001814430722574416"
             & "927381675121349900082925023049699309086508928992069"
             & "550881301185575782410245546511568736995455512201585"
             & "618034423562581432496255241806898294681614491811868886";
    cost049 : constant string
            := "0.9882575677307494914047925383512565371364538532078"
             & "629899363640808687435572386103632870868405863409365"
             & "007178603879751614794979760897337986868276575329976"
             & "737663714381017952834152073268426850737609058973705"
             & "442492875471590747229238707365730782224682575321893159";
    cost050 : constant string
            := "0.9877841416445721542309690323667278618121977647853"
             & "197604087723990264825241998838387775404939975137203"
             & "179032224245572518557642721076285050165480504853636"
             & "931483704604971496463604456503931967012502755386134"
             & "873614564509402859248119912178433551415108063700205623";
    cost051 : constant string
            := "0.9873014181578583823998158018450177283203725606333"
             & "990723464856282547202161680413296923198214066593567"
             & "839494814302672049182982737111733992398650406808427"
             & "384122445198044599528679247704867483757579042667032"
             & "861899211611777163141205974800816142545534887926466114";
    cost052 : constant string
            := "0.9868094018141854769702359522345500231768165633873"
             & "553340954664046433480249332223007001772523804215731"
             & "333174559781516602153882728947906451694319286080093"
             & "092235961671967485041538832342131653014300291944479"
             & "849071425897887515159837904909999685272245224666686494";
    cost053 : constant string
            := "0.9863080972445986478632975243258948530479327400784"
             & "137071185702432957470187041762660276955538002550988"
             & "534046454887452544769644761909465704276975486168985"
             & "082021562619660758065166564227029113592494615897795"
             & "085862175308547762408225966376728077203449986563106244";
    cost054 : constant string
            := "0.9857975091675674247009949996071128769398239356131"
             & "899318496370830312513571190685468020293868455907788"
             & "213635083950973458550760897373551459099508883020888"
             & "274273454952419606052841346380459631428559527097498"
             & "246746009375148367234436866476456899204705313428947302";
    cost055 : constant string
            := "0.9852776423889412447740184331785477871601291558128"
             & "148744442534259035894601318364005524044913961710172"
             & "645813729490339621374249967889805539357479462643228"
             & "312957306574878983226322587835917287872219356507697"
             & "038403388319405731117602375939758349109890295407604863";
    cost056 : constant string
            := "0.9847485018019042185565531758024046506443781277586"
             & "363605051367748947393685626232759770893447373403575"
             & "237732519398557636531278280699897034853488153506974"
             & "183776302950785823789033918311386199479621261270945"
             & "267626725725103918811513434997070982692081339954121275";
    cost057 : constant string
            := "0.9842100923869290731938743872398332594066332257994"
             & "171627485610797752359154338893194184379871857008731"
             & "683279821884143514920955761479073290483050621126550"
             & "115800996591502479806857033816522174493276572442666"
             & "398212245165655826688672249327197476301184734904066990";
    cost058 : constant string
            := "0.9836624192117302743962377761507951158350349472112"
             & "216953632389302090396314694800686040230697660489981"
             & "458701520749269608038139948505755658680759337483085"
             & "648224595882190349789446356125587572329739028375409"
             & "844770660523773405662017786171998952368436564804813267";
    cost059 : constant string
            := "0.9831054874312163271803011546739453386706019222933"
             & "975526577941595694195922260171330812506503594499562"
             & "047512196040219313048888188959025917036728127809965"
             & "083569161974308555025163703474291359728229501261552"
             & "811484412195292373210866526647683476196839950522359607";
    cost060 : constant string
            := "0.9825393022874412559070403955777269171829204613654"
             & "459166324238584350164968892806025018566800114381842"
             & "599915053934868553996362921116463132372765456762443"
             & "105170212001956364208655876253855483206131498620751"
             & "092087304500499230040836404951972376960231742022548833";
    cost061 : constant string
            := "0.9819638691095552640728481538315649202651768988504"
             & "381096017890537567277169108198098255310590228840645"
             & "490771307225266195515116675838583246726078658549006"
             & "913032376405496924035415232899310622426897276127560"
             & "505748129713137960321634495411026937311299615795024607";
    cost062 : constant string
            := "0.9813791933137545743182241898789480320709097648168"
             & "140978905287205202548896135627687603661635922259410"
             & "372869709676264878779898126423173742012660437686024"
             & "689189920557668203056642538894425203280699473924069"
             & "404517180752654933008039007697762491891103643264703630";
    cost063 : constant string
            := "0.9807852804032304491261822361342390369739337308933"
             & "360950029160885453065135496050639150649858533007632"
             & "598948662798775784681310960848381701091485451909052"
             & "981223580423918286860736338652741318972946739839332"
             & "937486597435047390244869403252433116449612877676624275";
    cost064 : constant string
            := "0.9801821359681173926902100086453528464396877109657"
             & "589618763359728242709151596636564777583396340507656"
             & "320361087734413494625869959685527414592148435028826"
             & "369869078325815188645568369638831148476184420510832"
             & "831409861389690259518042186595104754704800565913931017";
    cost065 : constant string
            := "0.9795697656854405344393261098798955052132344937016"
             & "664832965963308279433648744956009086340504990027179"
             & "548115018951496677702645785702384834566450609502297"
             & "588178843653074065986233446540949745963038081483495"
             & "826543587611063830936072325127925149816563934826945378";
    cost066 : constant string
            := "0.9789481753190621947154801236603956128550203763183"
             & "727240459934244801542948952559810222201621349275438"
             & "624680701677759948708870253317777641747325987002414"
             & "531431297934829583957689250152932673918245874867178"
             & "641310797903310917292143458211305036022189309666964274";
    cost067 : constant string
            := "0.9783173707196276331062400968954989486662500360631"
             & "230816495786015811018698451759744903891812532253407"
             & "320754949995683594864480028801618367471278574333021"
             & "002890242432103674311571288068983190038455922808332"
             & "181403175633644562375654071247405740284565512229427861";
    cost068 : constant string
            := "0.9776773578245099799434047624729313055872572249685"
             & "078999434977059277796082462615007685729231323363141"
             & "915540440286638582988895306281316673512424827272043"
             & "548366774710857734357326551829766724453769358484727"
             & "562162668176235492351067036129956690924500075546505066";
    cost069 : constant string
            := "0.9770281426577543514858662110857144252619956039677"
             & "714046346396481316448798722415793590892244462790233"
             & "410184391327185338324466681844988843168595862489460"
             & "074535355225568175750944958379475374674419506997984"
             & "640521073306083482747739845860335097075500977057255837";
    cost070 : constant string
            := "0.9763697313300211493127321944898351365000802823869"
             & "514605372914800699461585511590539771473127256441116"
             & "325056991560270462048421802762207460402366263261098"
             & "120954053637253702385555796138697575210453281881421"
             & "777242882985991963824749105146234921143151991162805051";
    cost071 : constant string
            := "0.9757021300385285444603957664195279716440122657920"
             & "431654131848601497388340819490924121720223405848031"
             & "097546663046407981482597590735321123648043591111645"
             & "161960194052443483262229571507815075318294476853254"
             & "592103404932590525705973722373493521284304437441109967";
    cost072 : constant string
            := "0.9750253450669941468449134678423469945313388750890"
             & "490217564694191829259000462786791014357315903556766"
             & "970311520154341694557283627586910004799155262365946"
             & "912099856634043206054760249155998366562426621311048"
             & "888711860136518292383720435097888926748293041860603012";
    cost073 : constant string
            := "0.9743393827855758605187216681943645931425572615457"
             & "181085075592055635059417511288773973236613467944281"
             & "869388415754556232372036308846129323885963079968196"
             & "994851206972872770836481957286619140338505299943300"
             & "014336454403986099747183373795066667024767255590432891";
    cost074 : constant string
            := "0.9736442496508119253183839115181956093635632173765"
             & "468164373125622288663310232607292745308186902523029"
             & "637332944977404210854582273598151621885367904820822"
             & "564290259990605567864704084068891324556547448242397"
             & "680839005576593012024654420955191796037855485455973599";
    cost075 : constant string
            := "0.9729399522055601454677201139037965700243897536956"
             & "235170889072514855566027279337974215662072989261852"
             & "619521956867667444228330752269100458433589402847680"
             & "578701128949963089369532954503415808012760960551918"
             & "677781520765892425922601244033455589246176763453690218";
    cost076 : constant string
            := "0.9722264970789363057083211442241431612174246543477"
             & "987177757107756439266440148511065719788243213119213"
             & "086448139107788660366100462035231761667231962259603"
             & "724351800740791986297716938027840693360894608201198"
             & "211077233482781869929446195200894229846973062131890668";
    cost077 : constant string
            := "0.9715038909862517755370996218349531511232991382862"
             & "402580222206632399626060785245148587469449986959431"
             & "076234847968790375495922685552146179777627776215383"
             & "449521847675643355249406639445077430930906896963359"
             & "385604919818369797748706420843260782359897163274125613";
    cost078 : constant string
            := "0.9707721407289503021381696106902808737283328953573"
             & "998065283373776964952818598324852026948909859916561"
             & "538882401569383999909828682004588409285042030238420"
             & "485724950030958500830158589461642181275958069682093"
             & "117683461475146466300951428096220199211204151108273350";
    cost079 : constant string
            := "0.9700312531945439926039842072861002514568659622480"
             & "741009834974506799112392593276800890798848932125239"
             & "392165702750461151864748558447393744389813988141225"
             & "790573092007503538674201800057439002092944570791395"
             & "867972948366990151015060678006900485230960901113956103";
    cost080 : constant string
            := "0.9692812353565484860482907381059832428042708423117"
             & "998036815894241546156347813079270456268886928875822"
             & "272657880043213836783185074221875151196906016147217"
             & "782655102672272472212546353916706837811185778221236"
             & "730593415846334780652978270309233688812179936139120077";
    cost081 : constant string
            := "0.9685220942744173162210883289834431285217750218321"
             & "895003012659367825618090665948528965074262923873720"
             & "341546235384941242411192457498850496205425817698807"
             & "814102010969350374673729599857382743621914910993828"
             & "907311332670134659474282613180889627276534107362771056";
    cost082 : constant string
            := "0.9677538370934754652433919122446032943536820495674"
             & "942035253925875806217946428974688566568351055063505"
             & "094747857299358448159645456991294374610744466910880"
             & "157938678607703178563844898463130834918917754029961"
             & "972238610949120787986398892898157586485446205578788172";
    cost083 : constant string
            := "0.9669764710448521090872202259367862730479839136104"
             & "195989303513041730219611083739525922661999501009065"
             & "522955892289033260477852355688934187463186728482105"
             & "825275670899716069855762198056712812781794267624781"
             & "170712022706415517011303834471675123872847370260037034";
    cost084 : constant string
            := "0.9661900034454125554338329612223419786343840872174"
             & "320922991761673887630240952669599693548990044195387"
             & "761079957248905702662899363490193309108552632882405"
             & "303560185276228722414717945395547013021816990453928"
             & "962390377092688838788303797909945080661213544942095012";
    cost085 : constant string
            := "0.9653944416976893745508438575169030701867730131426"
             & "343205542520119051642586856569821790410914804860161"
             & "691422407432868923619611421339771447521420719177767"
             & "520330600480885803124008837340529427726989062618865"
             & "757794194896676688457820764828726112729880410178431047";
    cost086 : constant string
            := "0.9645897932898127238364321586277055315770498729752"
             & "100975463058254447312938729556496907601426695453131"
             & "598804803610392359981067067043997717501944742253053"
             & "948140043043445656050255072924547876433705411544157"
             & "483536169709179953724524855657173817091820074512536407";
    cost087 : constant string
            := "0.9637760657954398666864643555078351536630838488266"
             & "327043089160414011547947386650792172514930553151595"
             & "866123289104277313575261318315272060684430627529457"
             & "062679572069315752456140983758235747770378107226673"
             & "563977159491716952089339753194460622547439260545987121";
    cost088 : constant string
            := "0.9629532668736838863479214808508748735204565905897"
             & "460682687117760262068253087088095449390223470855219"
             & "494725702861942010088949852782967616529514917874167"
             & "439551236445038730096697298525000836890522957628459"
             & "230511288768290791322806026092946803370144599323915442";
    cost089 : constant string
            := "0.9621214042690415954296043162301533683259328575220"
             & "033671705138807821048031823343814048250043454930134"
             & "788895545768029996974527354458910297296282875199311"
             & "774848672652316948796257213677896561910059685776491"
             & "616424842604141506495129453479641796273795164541385960";
    cost090 : constant string
            := "0.9612804858113206417486596525191235450943555695075"
             & "193414309501472414000188115250119301129932441617631"
             & "827040776129881443195618116971562133632822954049915"
             & "467895033316569890951858729240129179749521031488159"
             & "110387415822115644216040914058665216098321302675378309";
    cost091 : constant string
            := "0.9604305194155658111990351376552656335496268220003"
             & "038114349196996602195860983666547482838877672058391"
             & "302116193311886523342575990980422927759357493161554"
             & "506715385882917009391949902007286328282302944764490"
             & "254525397910873836795280855217925759363905426852519077";
    cost092 : constant string
            := "0.9595715130819845283355281812303626134485965186116"
             & "516972353870707142388981665457541897260610856850318"
             & "891067501494265122392465594008443460070852214300889"
             & "065392844218730597506869709663429923401032058693320"
             & "526862044981366743685496716233037530726083368269893602";
    cost093 : constant string
            := "0.9587034748958715553746457917668690903737111090162"
             & "149144316197753884929659828284853532259509025070689"
             & "421279574207141621812708310802615855474367464078260"
             & "050306846328801791549641415692881578572460975969159"
             & "780371073129516395128634586276904852349287397435226183";
    cost094 : constant string
            := "0.9578264130275328903210370287966757208083863485438"
             & "522062087553172321154698789902443960686854432584119"
             & "800485449246420489158184398767090607894457241012704"
             & "733297445496018610181593674862364535436742577754589"
             & "623166302843922324819682624190247734028333187727646142";
    cost095 : constant string
            := "0.9569403357322088649357978869802699694828492056300"
             & "372613012071998841601453681608249571934246962711819"
             & "922698493296146407927856166538403476343873898787038"
             & "483921127420661559169492775476604492989574586292032"
             & "664974921476368928324844200782923623784640180102102040";
    cost096 : constant string
            := "0.9560452513499964432704798225393117130005746467262"
             & "919480204953574988727611552191831339080422060870827"
             & "081434477470634912163597183392771004520010605784706"
             & "744426813956503270846415023854497197163511408107376"
             & "682454556742282762440560295793537040133734974414858252";
    cost097 : constant string
            := "0.9551411683057707214981577123356394246483804498203"
             & "158529646521887671458409193643991259473977487539849"
             & "435622826524164591301473623416405255290739035812991"
             & "979941060417917073579965311726863340138131803187401"
             & "951535945331539401977390619754475556901648609246182488";
    cost098 : constant string
            := "0.9542280951091056297804307321904861426773098032184"
             & "780401644871465492339712759831325969298147172342623"
             & "947736803642481726033045773183251931480787730797727"
             & "260345750647040372397151473354337647841672983000241"
             & "619444104961010643583051091988319379167320354157023622";
    cost099 : constant string
            := "0.9533060403541938369167403827397938349040969010171"
             & "632867949529339188348896363636206990880373294669276"
             & "509252572913893836171225334338271505790460258061567"
             & "090335585000394782763794730994705128821699930740520"
             & "427772586709627612067802958957482534517874798310139805";
    cost100 : constant string
            := "0.9523750127197658585298936075710087775910962246224"
             & "436027367604950942596168858603174563383596614003172"
             & "550013590984907991364131254913524556327349749967706"
             & "363191238812829653964586099342548666020801181210350"
             & "381706424874581326208812197411296061434881250168058332";
    cost101 : constant string
            := "0.9514350209690083695491755689557821707500772607476"
             & "043079392068176881272993806457738474895188569845765"
             & "008128961530878853247035986429508819934679121008756"
             & "186987956019194001399132880388365990960666265523846"
             & "062677033208585853530945494557324393273965600739207592";
    cost102 : constant string
            := "0.9504860739494817217599261006205415490837005595036"
             & "579598658877148852072280006495223904425141991699141"
             & "856344011234229157694408761497404703246004046906927"
             & "891700598880789661769531929919112645701920030376073"
             & "806881672972858869275726087578389485661033603170216038";
    cost103 : constant string
            := "0.9495281805930366671959360741893450282522241538324"
             & "108524439709653973932083707234524215787997564023190"
             & "217594060381506442301498403297053621440694324110218"
             & "388017345870917938390607630563268737321669903055026"
             & "778902716120434211718259769306127428260089481462897710";
    cost104 : constant string
            := "0.9485613499157302881584948257653042498390321546836"
             & "932216198241182734181066441995003027368498367635506"
             & "693068170640346930117303752068893748922119846083671"
             & "779405011815142166704103844931710442657211880767518"
             & "177717669286944649847069067066909062974824705861762952";
    cost105 : constant string
            := "0.9475855910177411346533873212314649157949189830936"
             & "963643808307761880154723549967661047076255098680343"
             & "782408329870219563520003003446625318713874104284489"
             & "098441143815207664127778082963477197777417498583685"
             & "316833852270463112091946585677069800078383946319935748";
    cost106 : constant string
            := "0.9466009130832835700445998229621097795147633836406"
             & "999140748968915553641624483601833887564938312934209"
             & "814167196816630770799048435887144599691115596487590"
             & "470997106045800588120976104488211300854220203381025"
             & "729758774990820268562686802904041624880088551637881489";
    cost107 : constant string
            := "0.9456073253805213257309453865238450864520238660965"
             & "341455629718085724350291286140933010689868381632885"
             & "393497071967001543043999345728765389608974707499477"
             & "461498313884329941188590708501070240178894866070822"
             & "719577405860104360561333791906580784009385177778944631";
    cost108 : constant string
            := "0.9446048372614802656592654934661592984649847828579"
             & "160291202279954325714120987907086831346001983686749"
             & "156146103124427974532948537295428686786222221556389"
             & "598157774739943662787016838832630007234814837986961"
             & "798716098097205108446061441867879279519411334681274116";
    cost109 : constant string
            := "0.9435934581619603614953014453784386693433755609167"
             & "296540524324871497281597136463378068067271732514632"
             & "983988873848100704043524874787401245106704428865783"
             & "954892413806907311677307994925651029033092059846343"
             & "219823786659142791575791988532067831472835606406377912";
    cost110 : constant string
            := "0.9425731976014468792807587350218082231326560438123"
             & "738601186867913336940993537302886286541050626022328"
             & "707271042244841108527402025371096107592677010749286"
             & "745416226420075338619951474194979929303448080486922"
             & "424300645222993068896756244325572578050686401326313554";
    cost111 : constant string
            := "0.9415440651830207784125094025995023571855897958251"
             & "828675468258789699971203630496978846001890548036847"
             & "188371045666972199780856472480360553974495519123755"
             & "747961700877283280920100428817404371695602453944896"
             & "079388913707311364601113430502903798385826190254253772";
    cost112 : constant string
            := "0.9405060705932683237872913092520213689163721527655"
             & "615618922449578789587957300079664333610328080693847"
             & "265879187455280033723136566760870170901780571025833"
             & "550123692892897136718145130746176429693427071235689"
             & "991361666793584583072499141369409751180059111235039793";
    cost113 : constant string
            := "0.9394592236021899119626692458704222131583983521196"
             & "722128837739195491790912870070251174422800036033831"
             & "138594089180324772370327292565646115271794692284301"
             & "075110568644622677753809130119486526670902322800695"
             & "450657245620894190732879033285745562106203183366206607";
    cost114 : constant string
            := "0.9384035340631081121924207736047628846643663901209"
             & "036726806226122421818740213601685058494060320800211"
             & "099946492216710648648637025282833745645162809297389"
             & "101614149699435181715054971808244448969501095328927"
             & "966241800897373098289619993665028024568526657724201025";
    cost115 : constant string
            := "0.9373390119125749232018995933723808790288308132995"
             & "546963325657458411778002447696447032083161920987636"
             & "106024800686367006129621211582352120543261802601092"
             & "901372092745634897546839576611695940032250600888789"
             & "783110991942164306489276810593863228584842798653185572";
    cost116 : constant string
            := "0.9362656671702782465763109956857752575590818942129"
             & "451192645921308548966316671964632560149229495942432"
             & "466206336423295286475685661589549239608080638624193"
             & "149618097938943783865084294855414243223009253171525"
             & "816405658909646741353437513551259191525123112588369968";
    cost117 : constant string
            := "0.9351835099389475776422074797424988136309455357656"
             & "205372815936361740075933392808785627155405552851652"
             & "442630934265715323015815117931469537084349952485935"
             & "705006433134452639237927792252085160766490280266377"
             & "866607673784331287915278141262381175138148314127001707";
    cost118 : constant string
            := "0.9340925504042589147298778825471074110553903629225"
             & "686445110875327636676881405172620186886075039348755"
             & "511766167491278739005701322617707214672493910674629"
             & "504458570132765845422967347272995998514611973517759"
             & "209385318929108611781687093702898533139576120502220282";
    cost119 : constant string
            := "0.9329927988347388877116602555433024982950155205122"
             & "950488923147712766684426703203108827359068660545814"
             & "659569432560336407190523005268397275774895773985253"
             & "415942843047990721825008468053964009484568309459954"
             & "794751625122685316533052423462853068269878213276902803";
    cost120 : constant string
            := "0.9318842655816681067185571985770264225996390618719"
             & "338789353596792493529343904461903246532570809745433"
             & "789965979393590490867259922579310006773901396822629"
             & "959595753244827475940608248000898115347272134555196"
             & "154884685099871380657752657037559144800643580737013582";
    cost121 : constant string
            := "0.9307669610789837319448723398218081380297515646350"
             & "503346297894971826961242543993839102265110455496911"
             & "467003639029353153841378388064160058413847744218621"
             & "672080783879084233382453092262239037811057000870612"
             & "503807674288574070950435414617672441765002799052798179";
    cost122 : constant string
            := "0.9296408958431812654579180664894332423899170154979"
             & "732937462342707783857564402827043132235152837168917"
             & "667843601635326851408785943104927022830404754987214"
             & "434651216631020467937798291847936592333526301351497"
             & "075714356558362544226284668062289257765736562769965704";
    cost123 : constant string
            := "0.9285060804732155659371673957159456406234109867828"
             & "877479861615510769100582911511828005233458856097702"
             & "717418630343658668283731906080996935255721239666335"
             & "111882702338276043556052893996553257621852117992316"
             & "638989366787426824189487665192446106881512403541381831";
    cost124 : constant string
            := "0.9273625256504010872745369590302413069028138258771"
             & "417244020837490472943330669356551880528131282370273"
             & "941772838448474171009009393309223139221586885216848"
             & "574681891175559939321036281068422201341752246684072"
             & "495272639336320475615242409897177667943014593939454447";
    cost125 : constant string
            := "0.9262102421383113419747933884371432823411412506173"
             & "913153462233814000053246707029808628138042140535213"
             & "919362929011240216327826317891709572055238713450528"
             & "008607985807928876047913711653418083372430667111815"
             & "929297939235964376715323176725778315211134348076945411";
    cost126 : constant string
            := "0.9250492407826775903023718686184477409811358598791"
             & "183891915584516614361127689994280999165891352126804"
             & "051199332561263768040491751050419262693000841448494"
             & "785540228155293522938171762469579699833419955853320"
             & "695312833457818513558776270773403275792801250503663705";
    cost127 : constant string
            := "0.9238795325112867561281831893967882868224166258636"
             & "424861150977312805350075011023587148399348503445960"
             & "979630257822478830308691775799042014275332219995578"
             & "278983938373732927138059433771800144734486056051930"
             & "637452667632039014232534332759013298398643482170954809";
    cost128 : constant string
            := "0.9227011283338785704372642268248986908391320246055"
             & "309606847655752376964394813901218803654534356020026"
             & "270710758250956768690826558758442194204204385066230"
             & "942094894589218939288263671954368935438117947153978"
             & "293587678805261011087603883059515849968497283351116254";
    cost129 : constant string
            := "0.9215140393420419434653963315480622263672465115990"
             & "260347766029679984936339535886574875348797982385112"
             & "818102009493514376737352339371821916229709747935375"
             & "923293331366148832845744118424650900808200980344703"
             & "097793963268105678464045400741192771315883376316334981";
    cost130 : constant string
            := "0.9203182767091105664400765410164427383299927934101"
             & "930909504779609513481251059397152247817000497521287"
             & "604495910517751570283375685078534248854817400476223"
             & "919263791038827409932917237196828498169201722560235"
             & "098640640422925579692974568585855078991509661954387193";
    cost131 : constant string
            := "0.9191138516900577439084777893585916788760254170155"
             & "826476013631182379061388510154227887534912418462755"
             & "963565031760114255395288822124178644437523686130697"
             & "308191511515527590990012311570537057903726373960996"
             & "469006480890319482093257424458564522334234063085697246";
    cost132 : constant string
            := "0.9179007756213904576422762970161218427527792101607"
             & "703185971662225004315031684208421819211416714707249"
             & "659731243077828120808897256767198429723187607886354"
             & "060559829652184408235161627261671874933703365905008"
             & "964879643716179008932099787364983497459078084670023881";
    cost133 : constant string
            := "0.9166790599210426631164570134177923250274471032269"
             & "164055587242622222926817895136578874369584703071366"
             & "161133635270531032986293594408835338321442614595982"
             & "894837224537287817160935498061367171788321202656797"
             & "959905176762126035195389640810930621654318898375230429";
    cost134 : constant string
            := "0.9154487160882678195664312919622163783863060352317"
             & "183605341850852170136655263558567909094958677925637"
             & "054898529874722935142519247384985389511124572545787"
             & "395947212258970075247608081518567965381793095904093"
             & "290825970052209925949711790308157695814500739897082843";
    cost135 : constant string
            := "0.9142097557035306546350148293935774010446911156821"
             & "770013565662410591304975518946726379074291674786932"
             & "682938611970114452390696598450028374233207993449586"
             & "046161934283448895027521883815056943107507623179301"
             & "424908033191600140168195756841642186247079739874094499";
    cost136 : constant string
            := "0.9129621904283981646280182333945881639288239722808"
             & "657867041241923202400392828472730809165261553569508"
             & "411384188554634255113090726752087766639352935909359"
             & "587023054054376398946514893520151699449342249801276"
             & "918904709658040435734208896208845781831390925324925535";
    cost137 : constant string
            := "0.9117060320054298514043973250755404989166873134114"
             & "448679139587726688467450341916215372487880953275350"
             & "734449029108069701741902307232499220158583679788263"
             & "185527926118149266897270931533828742508043376197612"
             & "314979376013500742636324029594594276700980120996801570";
    cost138 : constant string
            := "0.9104412922580671969340953692880071699802655492372"
             & "571006413260539468868768253516712935180550071088981"
             & "642490061003236003411218499297743386388560440866363"
             & "379201043646429385397463366429136248471037291583876"
             & "305213432526166400369381399348403059085735876719180154";
    cost139 : constant string
            := "0.9091679830905223765638847877078063304794860514149"
             & "418013493945486915413433040565043848415036564172092"
             & "068602274258642336248010015988971781204759789307377"
             & "243686337655668966012921420282712896129604165820357"
             & "025700876657113054685872335985537991582820307319377288";
    cost140 : constant string
            := "0.9078861164876662120386814798769818177571351230944"
             & "045028629650232143051085886448449027401107299512114"
             & "804836802691744754718201181112208358112779749186864"
             & "098740127040502166468944521071420550361345585808008"
             & "453739858922197065243813866928792068637735023562566722";
    cost141 : constant string
            := "0.9065957045149153653329605884237134126506369603843"
             & "768155054695993610692483759009919065695462661868140"
             & "788511698218036608866630200951702020083868998448130"
             & "230247058373456042800534737568426330363131849570870"
             & "453981246775358738055039125986277934627474252691905286";
    cost142 : constant string
            := "0.9052967593181187743540483291399726543292442596242"
             & "562214947388471444646149536217812978623617410026166"
             & "949741810274832805908461140320228406725219037503536"
             & "989274690160055849699699266086874950448608875113063"
             & "739565559114787274263748971478381174801933022064352676";
    cost143 : constant string
            := "0.9039892931234433315862002972305370487101320250506"
             & "080496646735759958654405663192393572459493188167709"
             & "686373412858922707859253343703234618218445090546064"
             & "837391160856885997460886509554957672985632690634879"
             & "433558008388897981412849906534900252608730908954296993";
    cost144 : constant string
            := "0.9026733182372588067515023906888894471608558820282"
             & "872818208427247102924789114202515037631176859530826"
             & "863072234237108576692885412297515227043581242041317"
             & "696618685062407714878248429819753500104910357981523"
             & "385424855780630306736662947747683671717529790677891489";
    cost145 : constant string
            := "0.9013488470460220145707460933354787941262707344023"
             & "840430179336531474003903532558241555053161255234923"
             & "638620946759720839918452658395888714256390824053396"
             & "741922549098025931546566306442521299611512591012166"
             & "214482753531296742041390848636310687269043588118102078";
    cost146 : constant string
            := "0.9000158920161602287145352665970513170516412731865"
             & "543879793180050977045286955066038426810044169324302"
             & "716509044047758313084815333473895846061407344841519"
             & "427170780888719375716501263867789680998299554542280"
             & "331935010225902254754719531367305662523132775076590087";
    cost147 : constant string
            := "0.8986744656939538430419767437334850971686375082110"
             & "164316056962924551981783325052487723219296180128524"
             & "264508501372651635826358481686922415569591622209021"
             & "735482262844438686028882632703153078048349292530625"
             & "105259806416656769317969390954537685931544576581531798";
    cost148 : constant string
            := "0.8973245807054182812313918361485724100227964935611"
             & "997250374640428573507988124307865275701413632768739"
             & "634443684821500488585011281460344587293677179773958"
             & "640241516137892248054083373607002240950192253454780"
             & "006386209939066166312996636008653809692361758738579884";
    cost149 : constant string
            := "0.8959662497561851559145602819685074912350183654394"
             & "935322797390060548191757512708320552930334719315116"
             & "383282954855960719239546316518094141468289325540929"
             & "722612003690451502423885179377131350614867746001557"
             & "397723105058682893324807214314373636156198583140185495";
    cost150 : constant string
            := "0.8945994856313826784330721256493119981033722667560"
             & "281534398173861636405663426352853061555401913127539"
             & "125011229424759933894216042650654370393721585120774"
             & "479602027103289447482194192645836643187226799556635"
             & "829563356188326728948760355090187060701457371208354489";
    cost151 : constant string
            := "0.8932243011955153203424164474933979780006255889988"
             & "727896079334615180005880405975123811164719357263991"
             & "195187766929392076457777412648507331833725128916872"
             & "764364733434048927468582172374131820076847095169319"
             & "671562492432511795219289344407335365081675125833365813";
    cost152 : constant string
            := "0.8918407093923427277964786972263580580505297928777"
             & "078756047695710567971745844367413729639262820804910"
             & "986579633662184349547246249321321275924712903683690"
             & "004237603423019183018189510300137716317335761526976"
             & "978223012276489534181511049894006417067315734732451737";
    cost153 : constant string
            := "0.8904487232447578899521505599180370203448336066920"
             & "747170908443080644308771945984098359575822232874916"
             & "336758053355711883600666491396597653030615421281990"
             & "077705411918871497671857100377022438514768004467079"
             & "998272777253460711427952293931155661024810831464436068";
    cost154 : constant string
            := "0.8890483558546645625407777293374767964898556757438"
             & "253030803481144320686554397899980082427075664283992"
             & "961044178603066188631928990097906713501803184482077"
             & "373306870392070004296902979440621556659794194664340"
             & "041912540561578072348818570639858078497405243863409118";
    cost155 : constant string
            := "0.8876396204028539477601816172209067692593291468142"
             & "088200144762682216257218234955710341178697470831918"
             & "157328295087262605310090189786060451396876998788811"
             & "085926118804588487472394117300846026143943984866552"
             & "085978088766141801108965215320608290996568696993990748";
    cost156 : constant string
            := "0.8862225301488806316479908209186334139362462577180"
             & "094050642547035792985516040336626859369658422637591"
             & "303295486232801297491900231883467693215529293305297"
             & "239569905288274424720675357720784633527881160694720"
             & "349052478087307068314423710599563871514739592951611209";
    cost157 : constant string
            := "0.8847970984309377801040070405857408897590863819469"
             & "723775926315124450646871598995194801771557952983844"
             & "030665668740156902220224770101743965130196734257678"
             & "842975234249694755751339313019892802919553765028169"
             & "135966302977097898770442260830743547781315465688728725";
    cost158 : constant string
            := "0.8833633386657315947363080147110533129482287388003"
             & "985211877201370905411934243708024170510672647602797"
             & "983291728415926548255064364204197109382693187094829"
             & "444270035393032979541008255260124021865550743264206"
             & "626827952578673735216144549364379448785324299072160954";
    cost159 : constant string
            := "0.8819212643483550297127568636603883495084426206747"
             & "279806325386167120666470450034970776205819055333733"
             & "858182337103195704203212736540957959522145572052095"
             & "447098687236291489400218618820580374931724583901277"
             & "873316092149643859046996094150297367252323167682812669";
    cost160 : constant string
            := "0.8804708890521607708065429294693790235531168140817"
             & "586696008925956573889244321730609730506078150494142"
             & "579117597149350001483360848768538676238457331918091"
             & "610516743376922818381712065870533423383377502636398"
             & "825869364206028280100909119777315820971320252830546676";
    cost161 : constant string
            := "0.8790122264286334778313237108883723441412213430684"
             & "602518469326246002622438112612767452207407805984759"
             & "724359400005161495820232508892177440625370557720355"
             & "402292896260394818518484682504665071451348149777388"
             & "454634074646125855872219793016816328998460715419780753";
    cost162 : constant string
            := "0.8775452902072612916684707502927493687195095184952"
             & "529155629858048283418734779969913947608768631341774"
             & "315837977144313256104496284145404622478506504866151"
             & "533630085002814715763724473666606074920939868039834"
             & "401133930733256591619157013614246526317900975396364034";
    cost163 : constant string
            := "0.8760700941954066070958442682679904961152886604629"
             & "681197407004228278540027820264873629675303398073130"
             & "325121817520452017768736484533940854332242980139346"
             & "330648741537764962567850726235172334667359026176764"
             & "083252058535604676444310264480787022865366795433663528";
    cost164 : constant string
            := "0.8745866522781761126344318973080041997486952406045"
             & "097017221737198284701982278791467158548106604733500"
             & "243706860395435591659391466244834954015991879172346"
             & "620612129183152466130590176082453774297221520876537"
             & "055833827326266010465294567813857106100634327748400993";
    cost165 : constant string
            := "0.8730949784182900986360859730828656214301581468010"
             & "053820263561488419767988141208542808846027449901986"
             & "621568637702861575167992616473675334434038370663675"
             & "271947536967240691489053917891875821902911153600375"
             & "059769326312512827793734548530001288263026011688040879";
    cost166 : constant string
            := "0.8715950866559510348424814352011506891190164630450"
             & "746728682633457394707186954452361118727204777353454"
             & "175325625503650235158682768758781400007892028624623"
             & "180447275079775496181894225097436126942054038029929"
             & "389459457139532561213668851296095194581721295767739160";
    cost167 : constant string
            := "0.8700869911087114186522924044838488439108277895298"
             & "254871093938228539361812818416287793461042575249436"
             & "405558234224224766132274228313542272339368380284157"
             & "170666653832534860272744930731357304470099851226832"
             & "195536880510906796716558006994974871709641107193534089";
    cost168 : constant string
            := "0.8685707059713408953404498757617203040593764213300"
             & "616888403051165700583584971260315630190900402787686"
             & "606141138178610181801520822263318866364454320785757"
             & "298665627535736235780047108293726603836976250899643"
             & "056997801568618308692013070100588911592992674131202045";
    cost169 : constant string
            := "0.8670462455156926514801956294958485611965851201621"
             & "297660798725181487863996483140811564586056150781468"
             & "485088891407574136942268928652325669606270709248056"
             & "079205755878493072795145456117326239497692572512845"
             & "468865424749118252989569331489948603729479768774770524";
    cost170 : constant string
            := "0.8655136240905690828254883576021393973799296090407"
             & "015002253905303982477455750015399792114601557558816"
             & "311728976300416185433815982839302068394072690010084"
             & "502369510165872832439117703573843618125161962271824"
             & "478620653941588505334762574928364782337630523767131513";
    cost171 : constant string
            := "0.8639728561215867379181470543455525270742850950242"
             & "317811555335678079167517202899177077972415950589790"
             & "376901322859688742463731012624713627807876094611736"
             & "940382184839195327457883207555211599740948091038688"
             & "592978722042498606888085154080046015553234615023831628";
    cost172 : constant string
            := "0.8624239561110405386909338777812499507348754626281"
             & "270153948166998718123396320322869767530182078194690"
             & "535814356515798814791806043933767343433983971726827"
             & "247913537427167701937914563397114492619319307706604"
             & "478392921873143649717990880297048194413537943501821062";
    cost173 : constant string
            := "0.8608669386377672793445838767919507564613366474625"
             & "723387199471186396489662014555571438079539062719344"
             & "711524565345691568293408078626674836309963052773137"
             & "117358184152980202661994382162239945731833900804250"
             & "326737971145035847333178560464957705718105687141018923";
    cost174 : constant string
            := "0.8593018183570084047835821392505594300776504877935"
             & "064551046923612187469093103698470644702895615236871"
             & "627258739117232923418895798382057773699086770035917"
             & "475591656115382175602651793172572324260049316164348"
             & "928428695403833467141830444779365241048700655330036780";
    cost175 : constant string
            := "0.8577286100002720699022699842847701370424907994337"
             & "340186047185423495190873294644473503467443827372978"
             & "075241933570661095405205680771980831058388087056657"
             & "416312480446175849256354721324434157530376262117206"
             & "579783448508375330774264558685873586110592419030421394";
    cost176 : constant string
            := "0.8561473283751944810196307322098535510157021598510"
             & "809975165737872013594207560187965013852508157240295"
             & "969495589085045055356485352241037344562891506790312"
             & "370270095007636429793610626652670314598835123239215"
             & "865820392491926771033599106390060472112090349529599367";
    cost177 : constant string
            := "0.8545579883654005207678622757156412315482088645418"
             & "588110512777579331631426671302978508560211359238233"
             & "785652148547697285290049303292081949746679019154961"
             & "550923500623864924682870214705223278961773562847887"
             & "249825273969097563752539875304531192190378916975094989";
    cost178 : constant string
            := "0.8529606049303636577465880817472955612247429334215"
             & "427919944020427523011457280123176945055587586258150"
             & "967346884391913724395546699436321107760990816403452"
             & "140254557156327076566974473919413269237333137912814"
             & "470258307724446989283162022074884034152031364817333724";
    cost179 : constant string
            := "0.8513551931052651422612903117258743556752017362384"
             & "884190667743706881798048993733629191335831105877273"
             & "749353298792744066798179376960853386224848393073801"
             & "759240725199219549976536232122353001679529257569133"
             & "848218278973143720168303712088212679835881197452940258";
    cost180 : constant string
            := "0.8497417680008524894712683949492727307512283411144"
             & "292886330832639998779886779551721857643521885904581"
             & "879050984940178687228307032853146962157971311277727"
             & "302227335455545023342133018653504686562429460297910"
             & "763948352503991250239727344878028763407970886562974903";
    cost181 : constant string
            := "0.8481203448032972512791335629489720926244788385512"
             & "470982499046933124708185509683309949120128627571915"
             & "130218316031682816771366873017819651088760664528186"
             & "826719757976732231645482199331979699272787293189905"
             & "458460896817189429637127642905165975252974293409056391";
    cost182 : constant string
            := "0.8464909387740520783005444881226811347222124933243"
             & "732803616230977976007739805772335103119046083987670"
             & "200933555266817880758564577811191835561922760351072"
             & "653603752584739296848838850992100678746918635582907"
             & "227799414016068593933737111267109654109603309449727349";
    cost183 : constant string
            := "0.8448535652497070732595712051049570977197859813891"
             & "086266626210143385371786920747277456118093123217801"
             & "672284127815307980250456681241186968318046987479170"
             & "149687666054359960194575439985487245816687792143540"
             & "993304109938224612552989635796214065640098233761574864";
    cost184 : constant string
            := "0.8432082396418454371617438651833876812949431645086"
             & "799259749346880129364719994550891835869905317673165"
             & "238820679271058907453505840236527843443102797466841"
             & "136755188697792716495751133191965836429805505631172"
             & "194093807501829269453354203007408210774523467418099377";
    cost185 : constant string
            := "0.8415549774368984096034995198422270646026529096920"
             & "897822933798480571048569463580896837136732489074750"
             & "928916766346523473312179072962381478955311593518841"
             & "372103856673026586597119022038097749819329155877524"
             & "607286343395862461090787467783310282270846413112137474";
    cost186 : constant string
            := "0.8398937941959995045833839865636283006705890551261"
             & "649383618330564158183111843403396582187935747266575"
             & "296388100943757636812646643192896759091248597200636"
             & "005404555161539609885684638511748885641349212603425"
             & "880128295864851571948394642529787763560691725829150891";
    cost187 : constant string
            := "0.8382247055548380431869968558042860626362709704203"
             & "287236075646985850546882622337805329101311213571178"
             & "587164549514256673277018441313948628286024228673699"
             & "612373820982371128173201484290425327544841332421739"
             & "169381498409444484878665596177169428875540110925182345";
    cost188 : constant string
            := "0.8365477272235119845242857901788464603635954918007"
             & "066490321657778572992257517778545426331240553942111"
             & "453337458115429376977130624425004009409809485689035"
             & "921739374717002881963629639267405913014965888613395"
             & "936407935176279931462901079854378227301558656546568069";
    cost189 : constant string
            := "0.8348628749863800563044013830288509569496169938396"
             & "501691698247621759545705834412655429270889159855589"
             & "619476910317250613553453935567219791610215498667567"
             & "275981229108835244276033092363666761352909720793560"
             & "567112474124966645009255512494352483000864042141580156";
    cost190 : constant string
            := "0.8331701647019131864399159215673905735273240019246"
             & "449123600064517642821648358427136153090160455324920"
             & "141559552077003278017431320396872868721318555891637"
             & "550434894491766596291822787716228639563704717490681"
             & "812535934405033677968712514866677929741763928420364683";
    cost191 : constant string
            := "0.8314696123025452370787883776179057567385608119872"
             & "499634461245902276379201446423439817749190079806503"
             & "360029402384434497288691970476296799430585600237731"
             & "382115302434462470313448879685115835814704970145904"
             & "782763774608865508293008705040105999560660210384482696";
    cost192 : constant string
            := "0.8297612337945230424690237646868315884264976428624"
             & "741896421886786974252694594756932508255874420987909"
             & "567013516541602032822213007055640646341451236946411"
             & "785804663301534385766364308031396109917263068534632"
             & "190912334783384713849398132011969330543714855124662018";
    cost193 : constant string
            := "0.8280450452577557520675275919240992974059373827702"
             & "435013402251673004249376370670847831166593190374881"
             & "776133698714113092717035508578000835834315107376592"
             & "811548869031715516079052019284224792048193392232616"
             & "129472464577175060808771137454201966991824807933343369";
    cost194 : constant string
            := "0.8263210628456634803111954517573735930748846670815"
             & "290526874387859085748608413637484304305595770407856"
             & "863984003157901665804313354600797205983924603387421"
             & "868486276539837205858402066574013346630277107712403"
             & "926767396330154813149686244257470189847480022867739865";
    cost195 : constant string
            := "0.8245893027850252644748037370848904463286468248009"
             & "647017983878967148285245862882494503860782365080786"
             & "744153591441870277785237871272095706568443595470149"
             & "557667723094592542153074672903654750054007184501719"
             & "167903399802858794589839747270709892262668807006190865";
    cost196 : constant string
            := "0.8228497813758263320467800344530491588480961622299"
             & "734493920894871901082223292255327703007795221417197"
             & "598951416093810991185342867924671869863499926815324"
             & "217659677951037804452909651646293832168294980542192"
             & "895916715905261480593044181690915898605806819186042359";
    cost197 : constant string
            := "0.8211025149911046790604308203298999031279085421450"
             & "958038286485691644319044054916388437027171445749242"
             & "223983561181943393282397082563412587060825660418888"
             & "681602099805054872010908744735797397397646708307788"
             & "520199664388495482221235472648030948456853388642631914";
    cost198 : constant string
            := "0.8193475200767969608246896372425308158363671409975"
             & "222751515759479898122243335921307009387027201996769"
             & "810715805870996139514132541894343134371980348305916"
             & "089332125797638654510534011556087403774459513677966"
             & "199616755821689995554297032397148116548307833859069095";
    cost199 : constant string
            := "0.8175848131515836965049208841306338094710425175669"
             & "140941589457011734569854789100989156894233606203180"
             & "184873664026944389114208199220536843189171199820539"
             & "095831801958003370976523598592521049421488957613132"
             & "785742220768964878326283624172444663214162270189987538";
    cost200 : constant string
            := "0.8158144108067337890107726598636859187152307591231"
             & "952886209737757258913335176731035680630728185876238"
             & "524870970146697720946384311332405963558723141532360"
             & "692998610438845589050929565017190873198214272072061"
             & "158495264649060955714005993315963521823856366893415265";
    cost201 : constant string
            := "0.8140363297059483616545166896872007386115265776746"
             & "003774044464669170594417949778338506169691073909138"
             & "906816176203157764999551897059034008881773435794550"
             & "300479938387353008502984923283700952777535650415289"
             & "903802449462323820381860007186266621833220048129889109";
    cost202 : constant string
            := "0.8122505865852039130497441807454391551551357211260"
             & "330985546644175991769320472814227448099419441185352"
             & "713506799518619654331777976069908300935610688145295"
             & "127187009635445067203494385377244887317930795958850"
             & "125460184047955213312665139122257953988021333637804209";
    cost203 : constant string
            := "0.8104571982525947917267034342445271670501183629908"
             & "300867998224525663956923539724574114552269252066104"
             & "086877539263619116755142330548458119971504775022574"
             & "541741728249727638665957277730328640095219420474183"
             & "403969191519009587125898493448259117165639126217073409";
    cost204 : constant string
            := "0.8086561815881749919469681278720381651798162762814"
             & "637093736461951585846591523845537267628945522920492"
             & "806004259693272600790253879857786242279611896846717"
             & "024632206496979428986092861520881846313410331613091"
             & "188534227679637929700045444770086645758283570901543822";
    cost205 : constant string
            := "0.8068475535437992722065143125084730510426323611548"
             & "135326799453325201134607217077195700453616968496817"
             & "170194097655426935719583202649925697054916621689643"
             & "417554499063314902167436485006652852898401149980336"
             & "177482426343885968988869862153713747201801587403547294";
    cost206 : constant string
            := "0.8050313311429635979226592819157852444787123174644"
             & "970190924831606995920885418544645815127031817692053"
             & "367970235907538532236848955656190584681607442292605"
             & "516876043250024420981531163046296576457283043792195"
             & "706257965205197512646947347864144602295261511223660048";
    cost207 : constant string
            := "0.8032075314806449098066765129631419238795694271704"
             & "608349765046543457740533480153272939564002841188636"
             & "971250153578526552530010885661343062055398088116519"
             & "105861847403701040332813708172494993795305252386827"
             & "406924818642259692200865764436703557141101335170138774";
    cost208 : constant string
            := "0.8013761717231402194302477771779672800629108667097"
             & "545508088530614950692727725548326704591481024484530"
             & "711160868288867797317350303472555809945172044613475"
             & "560915258482087179452970979123354924766258770230609"
             & "181511497602290210869727674143625302118757463821508947";
    cost209 : constant string
            := "0.7995372691079050335002462322515077643195102380441"
             & "946207778170017861806312346974341810108662847197345"
             & "623768915238715133112538694287904190134057758494284"
             & "153573434846878978197127125972519771545272992499325"
             & "005678253872204451670996202489653654262371460729798177";
    cost210 : constant string
            := "0.7976908409433911083626627549768352735344843743754"
             & "048329459382966025495650914084700581212228151074172"
             & "464570943941616037990599135563457366014255013219791"
             & "442050942134913044900360283996048223983910612025217"
             & "000730415063680726360755485233491564935173620238685938";
    cost211 : constant string
            := "0.7958369046088835362627919154816773610504779970095"
             & "616551123623411140461660007086782718225784605044266"
             & "041172990870198131525633514907456420893344881277651"
             & "020550192747465540410324933218908567558403959787051"
             & "696096770218060256744382631434851592118865367773240306";
    cost212 : constant string
            := "0.7939754775543371648950837572017825000808051746992"
             & "113525926753850058151174287248220758884883389652453"
             & "035952544896981820160492333033273433833218163143909"
             & "037936159970011636857461290005683027509181676681113"
             & "998045036654100100034898008881537305343014861521000981";
    cost213 : constant string
            := "0.7921065773002123517823428786210334566942030608291"
             & "357291633653650859464030368960855346219609943892685"
             & "144273844792870410830375068049655087346761471843072"
             & "158088416593364384303232708757055802347416690408877"
             & "663091904118582406803283442390574139965373058755186874";
    cost214 : constant string
            := "0.7902302214373100550302171523164021706369885772062"
             & "786870431947462777460459496932499438892559012534068"
             & "801756679515681621336886574275409276943229739710840"
             & "863706545251434242850069512280427532113189922240132"
             & "141292055370783744175381432592117436014697531064309109";
    cost215 : constant string
            := "0.7883464276266062620091647053596892826564931371496"
             & "486506948917380682596142739840278122197078706868444"
             & "118194652198692285559603098321321285698681291877082"
             & "575229653272134898822628299625540053055729020997350"
             & "954469486910307181508923146137551152927336879325011981";
    cost216 : constant string
            := "0.7864552135990857575223194638513642098426265082159"
             & "146577896986477583694630211488069540998437191961291"
             & "223603950270116029744001396408473888470479973323149"
             & "981649742110766605192636941930538209912269492364672"
             & "662099589369967117646606971945905169716133038988593836";
    cost217 : constant string
            := "0.7845565971555752330238925746397839911870155934056"
             & "389736686405662061245531157961897734539388161559437"
             & "799075059399312373791264536599711222251327760093927"
             & "155326267085902114213491498018354470891377225219842"
             & "446865056656227718386802734473782669469171336963571902";
    cost218 : constant string
            := "0.7826505961665757384589493005947525340922521134473"
             & "689623010314942817276545016964798441487213708274377"
             & "967957082220100633201521973210983610987892650128690"
             & "633330498627348559699054140868763539332832312757203"
             & "548526746787948387136400718549609431798920889394888153";
    cost219 : constant string
            := "0.7807372285720944783015884837795336966664392172732"
             & "182218575541571912949186280102646373014695612349156"
             & "112113482611362136388953620037599936843188763244034"
             & "234030134691119565994511941712010300386877671315247"
             & "208124197106369823743591462606368320147741373380347028";
    cost220 : constant string
            := "0.7788165123814759533747243252609644040752852140782"
             & "341266186630414670788747709445718121059506246577092"
             & "029584415733180816147319944829799435477522539652947"
             & "557146458797810873402645116775883944345575025996172"
             & "389606101033174511745119307305561332071337205516555277";
    cost221 : constant string
            := "0.7768884656732324500408279830138537082772169948111"
             & "833586492039017904288027304081652093260699741500773"
             & "302338513387247644368291355576552858279287047437866"
             & "127695876564964853792292522475872202623286537869092"
             & "972947813118387859695346757875430671441813829297436156";
    cost222 : constant string
            := "0.7749531065948738783591292824542866253702789416274"
             & "744528423192510499442521528468115690457394934866699"
             & "693413736913174207723046734513046686499067651910390"
             & "339631629411747035504536509877490331481323114710494"
             & "980791011809611610247774876241578549725031535075584645";
    cost223 : constant string
            := "0.7730104533627369608109066097584698009710412929008"
             & "096093564028966879506053059873020379548140068698393"
             & "742614841345656141239599015933191604613844325890155"
             & "820694712052901845914576238028944169888998549480787"
             & "285186485527254685250048432935860983329847720229695680";
    cost224 : constant string
            := "0.7710605242618137732006057586124016581017888104106"
             & "236957344745001546892119844197616344638790252564130"
             & "343893163048287054135342585468538861934390286780304"
             & "564310775469300896514022975241299896058651434840075"
             & "299615362012000289738298302548310513536923162352092214";
    cost225 : constant string
            := "0.7691033376455796393466260688578576671915689470561"
             & "573767814460706652638146767956577733764554233516231"
             & "180327946495696356196159565122469911571754039372493"
             & "912027922328695647788019704536927590884696674188497"
             & "946834796991254984951163714964423014879696150716732955";
    cost226 : constant string
            := "0.7671389119358203811816945732593211757841089136593"
             & "752564979783986311402782719432800690089805888915770"
             & "964238084871983571454629579077555039203396915710348"
             & "578217942252236805413358012852529009485520248384562"
             & "450416961951271919943102670402087123823553656731935203";
    cost227 : constant string
            := "0.7651672656224589258888159990649059180489680014138"
             & "812013888589144682411874807370646026874414570169458"
             & "384590161001362772241935555694149962299107816061549"
             & "568773434399740026262670070568036577688930918594192"
             & "788926097262645155608151281512449561663417775148539005";
    cost228 : constant string
            := "0.7631884172633812717048382970658654583073195549428"
             & "169546356156115837775787597651957490307989885858463"
             & "266061093738270172833041986570791411085455240872898"
             & "105979922346773437570092113205659244999414585415019"
             & "552618976214967933640684715860294373172157807219142696";
    cost229 : constant string
            := "0.7612023854842618140297098355118759758054588851057"
             & "932668222299793647276679764876301976803922542436070"
             & "157902264777075658697396242589931328292097923392271"
             & "575120217717976197007688744750753766857163590883374"
             & "413538868482537540783435084940451088297526574512487842";
    cost230 : constant string
            := "0.7592091889783880334855254426954076862902140668166"
             & "461023657046837885468643354082282363798126974557708"
             & "643235337239088695687528975677039675099168546663350"
             & "791445101787344298382759416722007868716829443070389"
             & "836988746049288458686201361547968707656409917677795265";
    cost231 : constant string
            := "0.7572088465064845475754640536057844730404337157316"
             & "168500555117689736318985200425229373510021761200010"
             & "792622571275917205513778900219117993025457748532455"
             & "374704743438078624962912617592357978884302442648145"
             & "385417705703397598497186436146516032313931213674900435";
    cost232 : constant string
            := "0.7552013768965365275987107562452439785001995161895"
             & "532408754923372853750076790503787509113637768265447"
             & "350714715893963324234637813109439851455629690502119"
             & "616586675064439897625856098551429523916123620616452"
             & "819732952788976246415324869673886958081723108131644856";
    cost233 : constant string
            := "0.7531867990436124824834304856149119103704206458977"
             & "201688199814238436253236772063638125057986905394760"
             & "682303663068979816634917333308120069765761277335181"
             & "710956995507961208831954744264827125678205098271981"
             & "489731369557249447020870833296745122501701977025954433";
    cost234 : constant string
            := "0.7511651319096864112058194217842731511004673418053"
             & "341235480639161350598565816620550059271121534413610"
             & "077120229376113427063999134073311940308893910124263"
             & "054723679618331000171804793091517519404617633738688"
             & "256600084379819869742039264581985831476483593275662962";
    cost235 : constant string
            := "0.7491363945234593254692032567668843148057429825117"
             & "396427265396568131702532271355663097737255290303954"
             & "302284845481425704745725791008995315191669647628894"
             & "885130591508366035961474706383891342562896416057656"
             & "716362271672798628876836739403921813677691726853154277";
    cost236 : constant string
            := "0.7471006059801801443230788471989288460126118161006"
             & "954948669300443836435085533034929445145991521141590"
             & "539775544396063187312614990823528550362087255118856"
             & "478383848036981390368529693863416641786501677131388"
             & "195709294043255989601534864011666049585006964147779025";
    cost237 : constant string
            := "0.7450577854414659624079073102652695805364073315923"
             & "834408792947418477155037262988420534798000506329271"
             & "671063501378323616602418026806745887622256925022462"
             & "522881450785073645918622703569508257856797097388055"
             & "073417850912713415851443635653055843934355106839301072";
    cost238 : constant string
            := "0.7430079521351216935173622932982347968199697288722"
             & "255960808796217890160970908919259333209551781889205"
             & "644885132582932926955210449387882453646229767194942"
             & "704097783665335818518079438877498764286800368348734"
             & "083987539881716160779801742202359217551235626071514729";
    cost239 : constant string
            := "0.7409511253549590911756168974951627297289553093090"
             & "900457364120438466825355639597576967103438140682232"
             & "065883642126995332937542226601359824256320509819474"
             & "040054556478842439515003835463523859782054944504666"
             & "584627935193780290257696088936571224286543685698587933";
    cost240 : constant string
            := "0.7388873244606151479331165079192798134314610826232"
             & "796435862715351515555688712044074906819612060735618"
             & "750814918894902791475724825261302909442893734288696"
             & "137938294112205945620664964325705913001814404233068"
             & "729637446248937000421855449461358635143706855407562735";
    cost241 : constant string
            := "0.7368165688773698750901325201727469468678844583869"
             & "454714583259896983881579200727962109789520953452696"
             & "376811865263815458028203372031667581338514096314513"
             & "984534943254683775256493610773591653969739832055207"
             & "854187080199309238869850447318267715691104802098199446";
    cost242 : constant string
            := "0.7347388780959634645632236038195365703151066276956"
             & "077451749980610634161053426603278410790659650783109"
             & "111981021484960037020256529977576329814136485864733"
             & "516315789217837438613612726867554556908168481476827"
             & "083420860777681021187827677904134689449367330129218498";
    cost243 : constant string
            := "0.7326542716724128346155466488994934632929541703531"
             & "315202018267179836348103639635123735795806023642575"
             & "919801044736355985436647018667509777873647669258392"
             & "420223802758308616999163773569551734507144448046067"
             & "494555495474227758586600097155009014794862856271981144";
    cost244 : constant string
            := "0.7305627692278275611777588499757241467643701736445"
             & "517067231339860665244616691304998033233471852716385"
             & "818866236733863542060798962517124566870653937734388"
             & "212348170213484202519949454926391367072179550112955"
             & "630319119453441764821640555407437628122399822382582607";
    cost245 : constant string
            := "0.7284643904482251964920354375100570986874932060934"
             & "284710721012965549394854709948367875391766589461991"
             & "270418515452596884368954019335324379030382596624790"
             & "224992085133983588595519627996745822420690832593802"
             & "979608845675485144658646107709003233723873637378709166";
    cost246 : constant string
            := "0.7263591550843459768174943145333992297606624731249"
             & "970521865590593986899205018907678919699763277902271"
             & "624006232812226985254659142758536174503036125514162"
             & "673896428713470182569187757474266116935286490960463"
             & "965421966143565561670875168194501373948066808481936404";
    cost247 : constant string
            := "0.7242470829514669209410692432905531674830930048004"
             & "368801650713787740884271114355639629308132763842398"
             & "050531783730636670194118500994201534533574744390214"
             & "268391804388173744031342520014196237512581120320664"
             & "340285363651637084726321183941634116289375362933925491";
    cost248 : constant string
            := "0.7221281939292153212436071976676250998279170867769"
             & "836496209438781382698583633964719256070148798860019"
             & "224375686390535053150143538824036674811823412636035"
             & "894879789147384680754494633199026081039656540004790"
             & "538395591841852303359362194568455669492192027707158668";
    cost249 : constant string
            := "0.7200025079613816290766829987842242033652881679560"
             & "536078689347656127839871312872838425296149569909680"
             & "334951399313675775948716762484993235014432810392588"
             & "381783942201339572289786022174273610274627170964008"
             & "880358995865732999706149210844915362075812561731354621";
    cost250 : constant string
            := "0.7178700450557317362113253293369297123442181563095"
             & "726148884065016740071136026826328360396420438539133"
             & "177318139657416671738505589733005852735839329063103"
             & "906027998165744554184163108549934941995987308272549"
             & "881246734685373011597046542804747119137365434810557027";
    cost251 : constant string
            := "0.7157308252838186541255326234552024122145960838025"
             & "807150765644749748984228732427893691391169957748675"
             & "480690950464559258239573343440357467530335143757047"
             & "514471952638184238403754950951354728811923599494913"
             & "431166840321368869149650209668159943571020968175771918";
    cost252 : constant string
            := "0.7135848687807935929031250994722950161934721469966"
             & "256124426691340233749164458502113745317734878817331"
             & "891444960841829336511702840732662368863020819875540"
             & "290923208676361541685241670885568900665597601033915"
             & "189848019325376802410251155532916796131301069887591147";
    cost253 : constant string
            := "0.7114321957452164415221302897745546729364570497435"
             & "620758261057567115582716221411872946223730266309968"
             & "851886195702205050546319843337468277673990219056645"
             & "176479905115618243554707136003960504716404623535179"
             & "432824614262977554170802300915323766414275122347334917";
    cost254 : constant string
            := "0.7092728264388656513165337715826299609990577295530"
             & "415148739423231633240870745750155623682388672793876"
             & "320376776854071233970627638237615713485708291508646"
             & "592014686392967610571045744155645841456388360890848"
             & "811333833588022711952337628997328251716422723502646611";
    cost255 : constant string
            := "0.7071067811865475244008443621048490392848359376884"
             & "740365883398689953662392310535194251937671638207863"
             & "675069231154561485124624180279253686063220607485499"
             & "679157066113329637527963778999752505763910302857350"
             & "547799858029851372672984310073642587093204445993047762";

    x : hexa_double;
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

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing hexa double arithmetic ...");
    put_line("  1. test addition and subtraction");
    put_line("  2. test multiplication and division");
    put_line("  3. test reading from a string");
    put_line("  4. write 16 leading doubles of log(10), log(2), exp(1)");
    put_line("  5. write 16 leading double for inverse factorials");
    put_line("  6. input and output");
    put_line("  7. Newton's method for sqrt(2)");
    put_line("  8. test the value of hd_eps");
    put_line("  9. test log(exp(pi)) = pi = exp(log(pi))");
    put_line("  A. write the 16 leading doubles for pi and multiples");
    put_line("  B. write 16 leading doubles for sine and cosine table");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, 9, A, or B to select a test : ");
    Ask_Alternative(ans,"123456789AB");
    case ans is
      when '1' => Test_Add_and_Subtract;
      when '2' => Test_Multiplication_and_Division;
      when '3' => Test_Read;
      when '4' => Log10log2exp1_doubles;
      when '5' => inverse_factorials;
      when '6' => Test_io;
      when '7' => Test_sqrt2;
      when '8' => Test_hd_eps;
      when '9' => Log_exp_of_Pi;
      when 'A' => Write_Pi;
      when 'B' => Sine_Cosine_Table;
      when others => null;
    end case;
  end Main;

end Test_Hexa_Doubles;
