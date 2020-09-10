with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;

package Triple_Double_Constants is

-- DESCRIPTION :
--   This package collects common definitions of triple double constants
--   for use in triple double mathematical functions.

  td_eps : constant double_float := 5.473822126268817e-48; -- 2^(-52*3-1)

-- PI and multiples and factors :

  pi : constant triple_double
     := create( 3.141592653589793116e+00, 1.224646799147353207e-16,
               -2.994769809718339666e-33);
  twopi : constant triple_double -- 2*pi
        := create( 6.283185307179586232e+00, 2.449293598294706414e-16,
                  -5.989539619436679332e-33);
  pi2 : constant triple_double -- pi/2
      := create( 1.570796326794896558e+00, 6.123233995736766036e-17,
                -1.497384904859169833e-33);
  pi4 : constant triple_double -- pi/4
      := create( 7.853981633974482790e-01, 3.061616997868383018e-17,
                -7.486924524295849165e-34);
  threepi4 : constant triple_double -- 3*pi/4
           := create( 2.356194490192344837e+00, 9.1848509936051484375e-17,
                      3.9168984647504003225e-33);
  pi1024 : constant triple_double -- pi/1204
         := create( 3.067961575771282340e-03, 1.195944139792337116e-19,
                   -2.924579892303066080e-36);

-- INVERSE FACTORIALS FOR TAYLOR EXPANSION :

  i_fac0  : constant triple_double
          := Create( 1.66666666666666657e-01, 9.25185853854297066e-18,
                     5.13581318503262866e-34);
  i_fac1  : constant triple_double
          := Create( 4.16666666666666644e-02, 2.31296463463574266e-18,
                     1.28395329625815716e-34);
  i_fac2  : constant triple_double
          := Create( 8.33333333333333322e-03, 1.15648231731787138e-19,
                     1.60494162032269652e-36);
  i_fac3  : constant triple_double
          := Create( 1.38888888888888894e-03,-5.30054395437357706e-20,
                    -1.73868675534958776e-36);
  i_fac4  : constant triple_double
          := Create( 1.98412698412698413e-04, 1.72095582934207053e-22,
                     1.49269123913941271e-40);
  i_fac5  : constant triple_double
          := Create( 2.48015873015873016e-05, 2.15119478667758816e-23,
                     1.86586404892426588e-41);
  i_fac6  : constant triple_double
          := Create( 2.75573192239858925e-06,-1.85839327404647208e-22,
                     8.49175460488199287e-39);
  i_fac7  : constant triple_double
          := Create( 2.75573192239858883e-07, 2.37677146222502973e-23,
                    -3.26318890334088294e-40);
  i_fac8  : constant triple_double
          := Create( 2.50521083854417202e-08,-1.44881407093591197e-24,
                     2.04267351467144546e-41);
  i_fac9  : constant triple_double
          := Create( 2.08767569878681002e-09,-1.20734505911325997e-25,
                     1.70222792889287100e-42);
  i_fac10 : constant triple_double
          := Create( 1.60590438368216133e-10, 1.25852945887520981e-26,
                    -5.31334602762985031e-43);
  i_fac11 : constant triple_double
          := Create( 1.14707455977297245e-11, 2.06555127528307454e-28,
                     6.88907923246664603e-45);
  i_fac12 : constant triple_double
          := Create( 7.64716373181981641e-13, 7.03872877733453001e-30,
                    -7.82753927716258345e-48);
  i_fac13 : constant triple_double
          := Create( 4.77947733238738525e-14, 4.39920548583408126e-31,
                    -4.89221204822661465e-49);
  i_fac14 : constant triple_double
          := Create( 2.81145725434552060e-15, 1.65088427308614326e-31,
                    -2.87777179307447918e-50);

  n_inv_fact : constant natural := 15;
  i_fac : array(0..n_inv_fact-1) of triple_double
        := (i_fac0,i_fac1,i_fac2,i_fac3,i_fac4,i_fac5,i_fac6,i_fac7,
            i_fac8,i_fac9,i_fac10,i_fac11,i_fac12,i_fac13,i_fac14);

-- TABLES of sin(k * pi/1024) and cos(k * pi/1024).

  sin_t0   : constant triple_double
           := create( 3.0679567629659761e-03, 1.2690279085455925e-19,
                      5.2879464245328389e-36);
  sin_t1   : constant triple_double
           := create( 6.1358846491544753e-03, 9.0545257482474933e-20,
                      1.6260113133745320e-37);
  sin_t2   : constant triple_double
           := create( 9.2037547820598194e-03,-1.2136591693535934e-19,
                      5.5696903949425567e-36);
  sin_t3   : constant triple_double
           := create( 1.2271538285719925e-02, 6.9197907640283170e-19,
                     -4.0203726713435555e-36);
  sin_t4   : constant triple_double
           := create( 1.5339206284988102e-02,-8.4462578865401696e-19,
                      4.6535897505058629e-35);
  sin_t5   : constant triple_double
           := create( 1.8406729905804820e-02, 7.4195533812833160e-19,
                      3.9068476486787607e-35);
  sin_t6   : constant triple_double
           := create( 2.1474080275469508e-02,-4.5407960207688566e-19,
                     -2.2031770119723005e-35);
  sin_t7   : constant triple_double
           := create(2.4541228522912288e-02,-9.1868490125778782e-20,
                     4.8706148704467061e-36);
  sin_t8   : constant triple_double
           := create( 2.7608145778965743e-02,-1.5932358831389269e-18,
                     -7.0475416242776030e-35);
  sin_t9   : constant triple_double
           := create( 3.0674803176636626e-02,-1.6936054844107918e-20,
                     -2.0039543064442544e-36);
  sin_t10  : constant triple_double
           := create( 3.3741171851377587e-02,-2.0096074292368340e-18,
                     -1.3548237016537134e-34);
  sin_t11  : constant triple_double
           := create( 3.6807222941358832e-02, 6.1060088803529842e-19,
                     -4.0448721259852727e-35);
  sin_t12  : constant triple_double
           := create( 3.9872927587739811e-02, 4.6657453481183289e-19,
                      3.4119333562288684e-35);
  sin_t13  : constant triple_double
           := create( 4.2938256934940820e-02, 2.8351940588660907e-18,
                      1.6991309601186475e-34);
  sin_t14  : constant triple_double
           := create( 4.6003182130914630e-02,-1.1182813940157788e-18,
                      7.5235020270378946e-35);
  sin_t15  : constant triple_double
           := create( 4.9067674327418015e-02,-6.7961037205182801e-19,
                     -4.4318868124718325e-35);
  sin_t16  : constant triple_double
           := create( 5.2131704680283324e-02,-2.4243695291953779e-18,
                     -1.3675405320092298e-34);
  sin_t17  : constant triple_double
           := create( 5.5195244349689941e-02,-1.3340299860891103e-18,
                     -3.4359574125665608e-35);
  sin_t18  : constant triple_double
           := create( 5.8258264500435759e-02, 2.3299905496077492e-19,
                      1.9376108990628660e-36);
  sin_t19  : constant triple_double
           := create( 6.1320736302208578e-02,-5.1181134064638108e-19,
                     -4.2726335866706313e-35);
  sin_t20  : constant triple_double
           := create( 6.4382630929857465e-02,-4.2325997000052705e-18,
                      3.3260117711855937e-35);
  sin_t21  : constant triple_double
           := create( 6.7443919563664065e-02,-6.9221796556983636e-18,
                      1.5909286358911040e-34);
  sin_t22  : constant triple_double
           := create( 7.0504573389613870e-02,-6.8552791107342883e-18,
                     -1.9961177630841580e-34);
  sin_t23  : constant triple_double
           := create( 7.3564563599667426e-02,-2.7784941506273593e-18,
                     -9.1240375489852821e-35);
  sin_t24  : constant triple_double
           := create( 7.6623861392031492e-02, 2.3253700287958801e-19,
                     -1.3186083921213440e-36);
  sin_t25  : constant triple_double
           := create( 7.9682437971430126e-02,-4.4867664311373041e-18,
                      2.8540789143650264e-34);
  sin_t26  : constant triple_double
           := create( 8.2740264549375692e-02, 1.4735983530877760e-18,
                      3.7284093452233713e-35);
  sin_t27  : constant triple_double
           := create( 8.5797312344439894e-02,-3.3881893830684029e-18,
                     -1.6135529531508258e-34);
  sin_t28  : constant triple_double
           := create( 8.8853552582524600e-02,-3.7501775830290691e-18,
                      3.7543606373911573e-34);
  sin_t29  : constant triple_double
           := create( 9.1908956497132724e-02, 4.7631594854274564e-18,
                      1.5722874642939344e-34);
  sin_t30  : constant triple_double
           := create( 9.4963495329639006e-02,-6.5885886400417564e-18,
                     -2.1371116991641965e-34);
  sin_t31  : constant triple_double
           := create( 9.8017140329560604e-02,-1.6345823622442560e-18,
                     -1.3209238810006454e-35);
  sin_t32  : constant triple_double
           := create( 1.0106986275482782e-01, 3.3164325719308656e-18,
                     -1.2004224885132282e-34);
  sin_t33  : constant triple_double
           := create( 1.0412163387205457e-01, 6.5760254085385100e-18,
                      1.7066246171219214e-34);
  sin_t34  : constant triple_double
           := create( 1.0717242495680884e-01, 6.4424044279026198e-18,
                     -8.3956976499698139e-35);
  sin_t35  : constant triple_double
           := create( 1.1022220729388306e-01,-5.6789503537823233e-19,
                      1.0380274792383233e-35);
  sin_t36  : constant triple_double
           := create( 1.1327095217756435e-01, 2.7100481012132900e-18,
                      1.5323292999491619e-35);
  sin_t37  : constant triple_double
           := create( 1.1631863091190477e-01, 1.0294914877509705e-18,
                     -9.3975734948993038e-35);
  sin_t38  : constant triple_double
           := create( 1.1936521481099137e-01,-3.9500089391898506e-18,
                      3.5317349978227311e-34);
  sin_t39  : constant triple_double
           := create( 1.2241067519921620e-01, 2.8354501489965335e-18,
                      1.8151655751493305e-34);
  sin_t40  : constant triple_double
           := create( 1.2545498341154623e-01, 4.8686751763148235e-18,
                      5.9878105258097936e-35);
  sin_t41  : constant triple_double
           := create( 1.2849811079379317e-01, 3.8198603954988802e-18,
                     -1.8627501455947798e-34);
  sin_t42  : constant triple_double
           := create( 1.3154002870288312e-01,-5.0039708262213813e-18,
                     -1.2983004159245552e-34);
  sin_t43  : constant triple_double
           := create( 1.3458070850712620e-01,-9.1670359171480699e-18,
                      1.5916493007073973e-34);
  sin_t44  : constant triple_double
           := create( 1.3762012158648604e-01, 6.6253255866774482e-18,
                     -2.3746583031401459e-34);
  sin_t45  : constant triple_double
           := create( 1.4065823933284924e-01,-7.9193932965524741e-18,
                      6.0972464202108397e-34);
  sin_t46  : constant triple_double
           := create( 1.4369503315029444e-01, 1.1472723016618666e-17,
                     -5.1884954557576435e-35);
  sin_t47  : constant triple_double
           := create( 1.4673047445536175e-01, 3.7269471470465677e-18,
                      3.7352398151250827e-34);
  sin_t48  : constant triple_double
           := create( 1.4976453467732151e-01, 8.0812114131285151e-18,
                      1.2979142554917325e-34);
  sin_t49  : constant triple_double
           := create( 1.5279718525844344e-01,-7.6313573938416838e-18,
                      5.7714690450284125e-34);
  sin_t50  : constant triple_double
           := create( 1.5582839765426523e-01, 3.0351307187678221e-18,
                     -1.0976942315176184e-34);
  sin_t51  : constant triple_double
           := create( 1.5885814333386145e-01,-4.0163200573859079e-18,
                     -9.2840580257628812e-35);
  sin_t52  : constant triple_double
           := create( 1.6188639378011183e-01, 1.1850519643573528e-17,
                     -5.0440990519162957e-34);
  sin_t53  : constant triple_double
           := create( 1.6491312048996992e-01,-7.0405288319166738e-19,
                      3.3211107491245527e-35);
  sin_t54  : constant triple_double
           := create( 1.6793829497473117e-01, 5.4284533721558139e-18,
                     -3.3263339336181369e-34);
  sin_t55  : constant triple_double
           := create( 1.7096188876030122e-01, 9.1919980181759094e-18,
                     -6.7688743940982606e-34);
  sin_t56  : constant triple_double
           := create( 1.7398387338746382e-01, 5.8151994618107928e-18,
                     -1.6751014298301606e-34);
  sin_t57  : constant triple_double
           := create( 1.7700422041214875e-01, 6.7329300635408167e-18,
                      2.8042736644246623e-34);
  sin_t58  : constant triple_double
           := create( 1.8002290140569951e-01, 7.9701826047392143e-18,
                     -7.0765920110524977e-34);
  sin_t59  : constant triple_double
           := create( 1.8303988795514095e-01, 7.7349918688637383e-18,
                     -4.4803769968145083e-34);
  sin_t60  : constant triple_double
           := create( 1.8605515166344666e-01,-1.2564893007679552e-17,
                      7.5953844248530810e-34);
  sin_t61  : constant triple_double
           := create( 1.8906866414980622e-01,-7.6208955803527778e-18,
                     -4.4792298656662981e-34);
  sin_t62  : constant triple_double
           := create( 1.9208039704989244e-01, 4.3348343941174903e-18,
                     -2.3404121848139937e-34);
  sin_t63  : constant triple_double
           := create( 1.9509032201612828e-01,-7.9910790684617313e-18,
                      6.1846270024220713e-34);
  sin_t64  : constant triple_double
           := create( 1.9809841071795359e-01,-1.8434411800689445e-18,
                      1.4139031318237285e-34);
  sin_t65  : constant triple_double
           := create( 2.0110463484209190e-01, 1.1010032669300739e-17,
                     -3.9123576757413791e-34);
  sin_t66  : constant triple_double
           := create( 2.0410896609281687e-01, 6.0941297773957752e-18,
                     -2.8275409970449641e-34);
  sin_t67  : constant triple_double
           := create( 2.0711137619221856e-01,-1.0613362528971356e-17,
                      2.2456805112690884e-34);
  sin_t68  : constant triple_double
           := create( 2.1011183688046961e-01, 1.1561548476512844e-17,
                      6.0355905610401254e-34);
  sin_t69  : constant triple_double
           := create( 2.1311031991609136e-01, 1.2031873821063860e-17,
                     -3.4142699719695635e-34);
  sin_t70  : constant triple_double
           := create( 2.1610679707621952e-01,-1.0111196082609117e-17,
                      7.2789545335189643e-34);
  sin_t71  : constant triple_double
           := create( 2.1910124015686980e-01,-3.6513812299150776e-19,
                     -2.3359499418606442e-35);
  sin_t72  : constant triple_double
           := create( 2.2209362097320354e-01,-3.0337210995812162e-18,
                      6.6654668033632998e-35);
  sin_t73  : constant triple_double
           := create( 2.2508391135979283e-01, 3.9507040822556510e-18,
                      2.4287993958305375e-35);
  sin_t74  : constant triple_double
           := create( 2.2807208317088573e-01, 8.2361837339258012e-18,
                      6.9786781316397937e-34);
  sin_t75  : constant triple_double
           := create( 2.3105810828067111e-01, 1.0129787149761869e-17,
                     -6.9359234615816044e-34);
  sin_t76  : constant triple_double
           := create( 2.3404195858354343e-01,-6.9922402696101173e-18,
                     -5.7323031922750280e-34);
  sin_t77  : constant triple_double
           := create( 2.3702360599436720e-01, 8.8544852285039918e-18,
                      1.3588480826354134e-34);
  sin_t78  : constant triple_double
           := create( 2.4000302244874150e-01,-1.2137758975632164e-17,
                     -2.6448807731703891e-34);
  sin_t79  : constant triple_double
           := create( 2.4298017990326390e-01,-8.7514315297196632e-18,
                     -6.5723260373079431e-34);
  sin_t80  : constant triple_double
           := create( 2.4595505033579462e-01,-1.1129044052741832e-17,
                      4.3805998202883397e-34);
  sin_t81  : constant triple_double
           := create( 2.4892760574572018e-01,-8.1783436100020990e-18,
                      5.5666875261111840e-34);
  sin_t82  : constant triple_double
           := create( 2.5189781815421697e-01,-1.7591436032517039e-17,
                     -1.0959681232525285e-33);
  sin_t83  : constant triple_double
           := create( 2.5486565960451457e-01,-1.3602299806901461e-19,
                     -6.0073844642762535e-36);
  sin_t84  : constant triple_double
           := create( 2.5783110216215899e-01, 1.8480038630879957e-17,
                      3.3201664714047599e-34);
  sin_t85  : constant triple_double
           := create( 2.6079411791527551e-01, 4.2721420983550075e-18,
                      5.6782126934777920e-35);
  sin_t86  : constant triple_double
           := create( 2.6375467897483140e-01,-1.8837947680038700e-17,
                      1.3720129045754794e-33);
  sin_t87  : constant triple_double
           := create( 2.6671275747489837e-01, 2.0941222578826688e-17,
                     -1.1303466524727989e-33);
  sin_t88  : constant triple_double
           := create( 2.6966832557291509e-01, 1.5765657618133259e-17,
                     -6.9696142173370086e-34);
  sin_t89  : constant triple_double
           := create( 2.7262135544994898e-01, 7.8697166076387850e-18,
                      6.6179388602933372e-35);
  sin_t90  : constant triple_double
           := create( 2.7557181931095814e-01, 1.9320328962556582e-17,
                      1.3932094180100280e-33);
  sin_t91  : constant triple_double
           := create( 2.7851968938505312e-01,-1.0030273719543544e-17,
                      7.2592115325689254e-34);
  sin_t92  : constant triple_double
           := create( 2.8146493792575800e-01,-1.2322299641274009e-17,
                     -1.0564788706386435e-34);
  sin_t93  : constant triple_double
           := create( 2.8440753721127182e-01, 2.2209268510661475e-17,
                     -9.1823095629523708e-34);
  sin_t94  : constant triple_double
           := create( 2.8734745954472951e-01, 1.5461117367645717e-17,
                     -6.3263973663444076e-34);
  sin_t95  : constant triple_double
           := create( 2.9028467725446239e-01,-1.8927978707774251e-17,
                      1.1522953157142315e-33);
  sin_t96  : constant triple_double
           := create( 2.9321916269425863e-01, 2.2385430811901833e-17,
                      1.3662484646539680e-33);
  sin_t97  : constant triple_double
           := create( 2.9615088824362384e-01,-2.0220736360876938e-17,
                     -7.9252212533920413e-35);
  sin_t98  : constant triple_double
           := create( 2.9907982630804048e-01, 1.6701181609219447e-18,
                      8.6091151117316292e-35);
  sin_t99  : constant triple_double
           := create( 3.0200594931922808e-01,-1.7167666235262474e-17,
                      2.3336182149008069e-34);
  sin_t100 : constant triple_double
           := create( 3.0492922973540243e-01,-2.2989033898191262e-17,
                     -1.4598901099661133e-34);
  sin_t101 : constant triple_double
           := create( 3.0784964004153487e-01, 2.7074088527245185e-17,
                      1.2568858206899284e-33);
  sin_t102 : constant triple_double
           := create( 3.1076715274961147e-01, 2.0887076364048513e-17,
                     -3.0130590791065942e-34);
  sin_t103 : constant triple_double
           := create( 3.1368174039889146e-01, 1.4560447299968912e-17,
                      3.6564186898011595e-34);
  sin_t104 : constant triple_double
           := create( 3.1659337555616585e-01, 2.1435292512726283e-17,
                      1.2338169231377316e-33);
  sin_t105 : constant triple_double
           := create( 3.1950203081601569e-01,-1.3981562491096626e-17,
                      8.1730000697411350e-34);
  sin_t106 : constant triple_double
           := create( 3.2240767880106985e-01,-4.0519039937959398e-18,
                      3.7438302780296796e-34);
  sin_t107 : constant triple_double
           := create( 3.2531029216226293e-01, 7.9171249463765892e-18,
                     -6.7576622068146391e-35);
  sin_t108 : constant triple_double
           := create( 3.2820984357909255e-01,-2.6693140719641896e-17,
                      7.8928851447534788e-34);
  sin_t109 : constant triple_double
           := create( 3.3110630575987643e-01,-2.7469465474778694e-17,
                     -1.3401245916610206e-33);
  sin_t110 : constant triple_double
           := create( 3.3399965144200938e-01, 2.2598986806288142e-17,
                      7.8063057192586115e-34);
  sin_t111 : constant triple_double
           := create( 3.3688985339222005e-01,-4.2000940033475092e-19,
                     -2.9178652969985438e-36);
  sin_t112 : constant triple_double
           := create( 3.3977688440682685e-01, 6.6028679499418282e-18,
                      1.2575009988669683e-34);
  sin_t113 : constant triple_double
           := create( 3.4266071731199438e-01, 1.9261518449306319e-17,
                     -9.2754189135990867e-34);
  sin_t114 : constant triple_double
           := create( 3.4554132496398904e-01, 2.7251143672916123e-17,
                      7.0138163601941737e-34);
  sin_t115 : constant triple_double
           := create( 3.4841868024943456e-01, 3.6974420514204918e-18,
                      3.5532146878499996e-34);
  sin_t116 : constant triple_double
           := create( 3.5129275608556715e-01,-2.2670712098795844e-17,
                     -1.6994216673139631e-34);
  sin_t117 : constant triple_double
           := create( 3.5416352542049040e-01,-1.6951763305764860e-17,
                      1.2772331777814617e-33);
  sin_t118 : constant triple_double
           := create( 3.5703096123343003e-01,-4.8218191137919166e-19,
                     -4.1672436994492361e-35);
  sin_t119 : constant triple_double
           := create( 3.5989503653498817e-01,-1.7601687123839282e-17,
                      1.3375125473046791e-33);
  sin_t120 : constant triple_double
           := create( 3.6275572436739723e-01,-9.1668352663749849e-18,
                     -7.4317843956936735e-34);
  sin_t121 : constant triple_double
           := create( 3.6561299780477385e-01, 1.6217898770457546e-17,
                      1.1286970151961055e-33);
  sin_t122 : constant triple_double
           := create( 3.6846682995337232e-01, 1.0463640796159268e-17,
                      2.0554984738517304e-35);
  sin_t123 : constant triple_double
           := create( 3.7131719395183754e-01, 3.4749239648238266e-19,
                     -7.5151053042866671e-37);
  sin_t124 : constant triple_double
           := create( 3.7416406297145799e-01, 8.0114103761962118e-18,
                      5.3429599813406052e-34);
  sin_t125 : constant triple_double
           := create( 3.7700741021641826e-01,-2.7255302041956930e-18,
                      6.3646586445018137e-35);
  sin_t126 : constant triple_double
           := create( 3.7984720892405116e-01, 9.9151305855172370e-18,
                      4.8761409697224886e-34);
  sin_t127 : constant triple_double
           := create( 3.8268343236508978e-01,-1.0050772696461588e-17,
                     -2.0605316302806695e-34);
  sin_t128 : constant triple_double
           := create( 3.8551605384391885e-01, 1.5177665396472313e-17,
                      1.4198230518016535e-33);
  sin_t129 : constant triple_double
           := create( 3.8834504669882630e-01,-1.0053770598398717e-17,
                      7.5942999255057131e-34);
  sin_t130 : constant triple_double
           := create( 3.9117038430225387e-01, 1.7997787858243995e-17,
                     -1.0613482402609856e-33);
  sin_t131 : constant triple_double
           := create( 3.9399204006104810e-01, 9.7649241641239336e-18,
                     -2.1233599441284617e-34);
  sin_t132 : constant triple_double
           := create( 3.9680998741671031e-01, 2.0545063670840126e-17,
                      6.1347058801922842e-34);
  sin_t133 : constant triple_double
           := create( 3.9962419984564684e-01,-1.5065497476189372e-17,
                     -9.9653258881867298e-34);
  sin_t134 : constant triple_double
           := create( 4.0243465085941843e-01, 1.0902619339328270e-17,
                      7.3998528125989765e-34);
  sin_t135 : constant triple_double
           := create( 4.0524131400498986e-01, 9.9111401942899884e-18,
                     -2.5169070895434648e-34);
  sin_t136 : constant triple_double
           := create( 4.0804416286497869e-01,-7.0006015137351311e-18,
                     -1.4108207334268228e-34);
  sin_t137 : constant triple_double
           := create( 4.1084317105790397e-01,-2.4219835190355499e-17,
                     -1.1418902925313314e-33);
  sin_t138 : constant triple_double
           := create( 4.1363831223843456e-01,-1.0393984940597871e-17,
                     -1.1481681174503880e-34);
  sin_t139 : constant triple_double
           := create( 4.1642956009763721e-01,-2.5475580413131732e-17,
                     -3.4482678506112824e-34);
  sin_t140 : constant triple_double
           := create( 4.1921688836322396e-01,-4.2232463750110590e-18,
                     -3.6053023045255790e-34);
  sin_t141 : constant triple_double
           := create( 4.2200027079979968e-01, 4.3543266994128527e-18,
                      3.1734310272251190e-34);
  sin_t142 : constant triple_double
           := create( 4.2477968120910881e-01, 2.7462312204277281e-17,
                     -4.6552847802111948e-34);
  sin_t143 : constant triple_double
           := create( 4.2755509343028208e-01, 9.4111898162954726e-18,
                     -1.7446682426598801e-34);
  sin_t144 : constant triple_double
           := create( 4.3032648134008261e-01, 2.2259686974092690e-17,
                      8.5972591314085075e-34);
  sin_t145 : constant triple_double
           := create( 4.3309381885315196e-01, 1.1224283329847517e-17,
                      5.3223748041075651e-35);
  sin_t146 : constant triple_double
           := create( 4.3585707992225547e-01, 1.6230515450644527e-17,
                     -6.4371449063579431e-35);
  sin_t147 : constant triple_double
           := create( 4.3861623853852766e-01,-2.0883315831075090e-17,
                     -1.4259583540891877e-34);
  sin_t148 : constant triple_double
           := create( 4.4137126873171667e-01, 2.2360783886964969e-17,
                      1.1864769603515770e-34);
  sin_t149 : constant triple_double
           := create( 4.4412214457042926e-01,-2.4218874422178315e-17,
                      2.2205230838703907e-34);
  sin_t150 : constant triple_double
           := create( 4.4686884016237421e-01,-1.9222136142309382e-17,
                     -4.4425678589732049e-35);
  sin_t151 : constant triple_double
           := create( 4.4961132965460660e-01, 4.8831924232035243e-18,
                      2.7151084498191381e-34);
  sin_t152 : constant triple_double
           := create( 4.5234958723377089e-01,-1.4827977472196122e-17,
                     -7.6947501088972324e-34);
  sin_t153 : constant triple_double
           := create( 4.5508358712634384e-01,-1.2379906758116472e-17,
                      5.5289688955542643e-34);
  sin_t154 : constant triple_double
           := create( 4.5781330359887723e-01,-8.4554254922295949e-18,
                     -6.3770394246764263e-34);
  sin_t155 : constant triple_double
           := create( 4.6053871095824001e-01, 1.8488777492177872e-17,
                     -1.0527732154209725e-33);
  sin_t156 : constant triple_double
           := create( 4.6325978355186020e-01,-7.3514924533231707e-18,
                      6.7175396881707035e-34);
  sin_t157 : constant triple_double
           := create( 4.6597649576796618e-01,-3.3023547778235135e-18,
                      3.4904677050476886e-35);
  sin_t158 : constant triple_double
           := create( 4.6868882203582796e-01,-2.2949251681845054e-17,
                     -1.1364757641823658e-33);
  sin_t159 : constant triple_double
           := create( 4.7139673682599764e-01, 6.5166781360690130e-18,
                      2.9457546966235984e-34);
  sin_t160 : constant triple_double
           := create( 4.7410021465055002e-01,-8.1451601548978075e-18,
                     -3.4789448555614422e-34);
  sin_t161 : constant triple_double
           := create( 4.7679923006332214e-01,-1.0293515338305794e-17,
                     -3.6582045008369952e-34);
  sin_t162 : constant triple_double
           := create( 4.7949375766015301e-01, 1.8419999662684771e-17,
                     -1.3040838621273312e-33);
  sin_t163 : constant triple_double
           := create( 4.8218377207912277e-01,-2.5861500925520442e-17,
                     -6.2913197606500007e-36);
  sin_t164 : constant triple_double
           := create( 4.8486924800079112e-01,-1.8034004203262245e-17,
                     -3.5244276906958044e-34);
  sin_t165 : constant triple_double
           := create( 4.8755016014843594e-01, 1.4231090931273653e-17,
                     -1.8277733073262697e-34);
  sin_t166 : constant triple_double
           := create( 4.9022648328829116e-01,-5.1496145643440404e-18,
                     -3.6903027405284104e-34);
  sin_t167 : constant triple_double
           := create( 4.9289819222978404e-01,-1.0257831676562186e-18,
                      6.9520817760885069e-35);
  sin_t168 : constant triple_double
           := create( 4.9556526182577254e-01,-9.4323241942365362e-18,
                      3.1212918657699143e-35);
  sin_t169 : constant triple_double
           := create( 4.9822766697278187e-01,-1.6126383830540798e-17,
                     -1.5092897319298871e-33);
  sin_t170 : constant triple_double
           := create( 5.0088538261124083e-01,-3.9604015147074639e-17,
                     -2.2208395201898007e-33);
  sin_t171 : constant triple_double
           := create( 5.0353838372571758e-01,-1.6731308204967497e-17,
                     -1.0140233644074786e-33);
  sin_t172 : constant triple_double
           := create( 5.0618664534515534e-01,-4.8321592986493711e-17,
                      9.2858107226642252e-34);
  sin_t173 : constant triple_double
           := create( 5.0883014254310699e-01, 4.7836968268014130e-17,
                     -1.0727022928806035e-33);
  sin_t174 : constant triple_double
           := create( 5.1146885043797041e-01,-1.3088001221007579e-17,
                      4.0929033363366899e-34);
  sin_t175 : constant triple_double
           := create( 5.1410274419322177e-01,-4.5712707523615624e-17,
                      1.5488279442238283e-33);
  sin_t176 : constant triple_double
           := create( 5.1673179901764987e-01, 8.3018617233836515e-18,
                      5.8251027467695202e-34);
  sin_t177 : constant triple_double
           := create( 5.1935599016558964e-01,-5.5331248144171145e-17,
                     -3.1628375609769026e-35);
  sin_t178 : constant triple_double
           := create( 5.2197529293715439e-01,-4.6555795692088883e-17,
                      4.6378980936850430e-34);
  sin_t179 : constant triple_double
           := create( 5.2458968267846895e-01,-4.3068869040082345e-17,
                     -4.2013155291932055e-34);
  sin_t180 : constant triple_double
           := create( 5.2719913478190139e-01,-4.2202983480560619e-17,
                      8.5585916184867295e-34);
  sin_t181 : constant triple_double
           := create( 5.2980362468629472e-01,-4.8067841706482342e-17,
                      5.8309721046630296e-34);
  sin_t182 : constant triple_double
           := create( 5.3240312787719801e-01,-4.1020306135800895e-17,
                     -1.9239996374230821e-33);
  sin_t183 : constant triple_double
           := create( 5.3499761988709726e-01,-5.3683132708358134e-17,
                     -1.3900569918838112e-33);
  sin_t184 : constant triple_double
           := create( 5.3758707629564551e-01,-2.2617365388403054e-17,
                     -5.9787279033447075e-34);
  sin_t185 : constant triple_double
           := create( 5.4017147272989285e-01, 2.7072447965935839e-17,
                      1.1698799709213829e-33);
  sin_t186 : constant triple_double
           := create( 5.4275078486451589e-01, 1.7148261004757101e-17,
                     -1.3525905925200870e-33);
  sin_t187 : constant triple_double
           := create( 5.4532498842204646e-01,-4.1517817538384258e-17,
                     -1.5318930219385941e-33);
  sin_t188 : constant triple_double
           := create( 5.4789405917310019e-01,-2.4065878297113363e-17,
                     -3.5639213669362606e-36);
  sin_t189 : constant triple_double
           := create( 5.5045797293660481e-01,-8.3319903015807663e-18,
                     -2.3058454035767633e-34);
  sin_t190 : constant triple_double
           := create( 5.5301670558002758e-01,-4.7061536623798204e-17,
                     -1.0617111545918056e-33);
  sin_t191 : constant triple_double
           := create( 5.5557023301960218e-01, 4.7094109405616768e-17,
                     -2.0640520383682921e-33);
  sin_t192 : constant triple_double
           := create( 5.5811853122055610e-01, 1.3481176324765226e-17,
                     -5.5016743873011438e-34);
  sin_t193 : constant triple_double
           := create( 5.6066157619733603e-01,-7.3956418153476152e-18,
                      3.9680620611731193e-34);
  sin_t194 : constant triple_double
           := create( 5.6319934401383409e-01, 2.3835775146854829e-17,
                      1.3511793173769814e-34);
  sin_t195 : constant triple_double
           := create( 5.6573181078361323e-01,-3.4096079596590466e-17,
                     -1.7073289744303546e-33);
  sin_t196 : constant triple_double
           := create( 5.6825895267013160e-01,-5.0935673642769248e-17,
                     -1.6274356351028249e-33);
  sin_t197 : constant triple_double
           := create( 5.7078074588696726e-01, 2.4568151455566208e-17,
                     -1.2844481247560350e-33);
  sin_t198 : constant triple_double
           := create( 5.7329716669804220e-01, 8.5176611669306400e-18,
                     -6.4443208788026766e-34);
  sin_t199 : constant triple_double
           := create( 5.7580819141784534e-01,-3.7909495458942734e-17,
                     -2.7433738046854309e-33);
  sin_t200 : constant triple_double
           := create( 5.7831379641165559e-01,-2.6237691512372831e-17,
                      1.3679051680738167e-33);
  sin_t201 : constant triple_double
           := create( 5.8081395809576453e-01, 1.8585338586613408e-17,
                      2.7673843114549181e-34);
  sin_t202 : constant triple_double
           := create( 5.8330865293769829e-01, 3.4516601079044858e-18,
                      1.8065977478946306e-34);
  sin_t203 : constant triple_double
           := create( 5.8579785745643886e-01,-3.7485501964311294e-18,
                      2.7965403775536614e-34);
  sin_t204 : constant triple_double
           := create( 5.8828154822264533e-01,-2.9292166725006846e-17,
                     -2.3744954603693934e-33);
  sin_t205 : constant triple_double
           := create( 5.9075970185887428e-01,-4.7013584170659542e-17,
                      2.4808417611768356e-33);
  sin_t206 : constant triple_double
           := create( 5.9323229503979980e-01, 1.2892320944189053e-17,
                      5.3058364776359583e-34);
  sin_t207 : constant triple_double
           := create( 5.9569930449243336e-01,-1.3438641936579467e-17,
                     -6.7877687907721049e-35);
  sin_t208 : constant triple_double
           := create( 5.9816070699634227e-01, 3.8801885783000657e-17,
                     -1.2084165858094663e-33);
  sin_t209 : constant triple_double
           := create( 6.0061647938386897e-01,-4.6398198229461932e-17,
                     -1.6673493003710801e-33);
  sin_t210 : constant triple_double
           := create( 6.0306659854034816e-01, 3.7323357680559650e-17,
                      2.7771920866974305e-33);
  sin_t211 : constant triple_double
           := create( 6.0551104140432555e-01,-3.1202672493305677e-17,
                      1.2761267338680916e-33);
  sin_t212 : constant triple_double
           := create( 6.0794978496777363e-01, 3.5160832362096660e-17,
                     -2.5546242776778394e-34);
  sin_t213 : constant triple_double
           := create( 6.1038280627630948e-01,-2.2563265648229169e-17,
                      1.3185575011226730e-33);
  sin_t214 : constant triple_double
           := create( 6.1281008242940971e-01,-4.2693476568409685e-18,
                      2.5839965886650320e-34);
  sin_t215 : constant triple_double
           := create( 6.1523159058062682e-01, 2.6231417767266950e-17,
                     -1.4095366621106716e-33);
  sin_t216 : constant triple_double
           := create( 6.1764730793780398e-01,-4.7478594510902452e-17,
                     -7.2986558263123996e-34);
  sin_t217 : constant triple_double
           := create( 6.2005721176328921e-01,-2.7983410837681118e-17,
                      1.1649951056138923e-33);
  sin_t218 : constant triple_double
           := create( 6.2246127937414997e-01, 5.2940728606573002e-18,
                     -4.8486411215945827e-35);
  sin_t219 : constant triple_double
           := create( 6.2485948814238634e-01, 3.3671846037243900e-17,
                     -2.7846053391012096e-33);
  sin_t220 : constant triple_double
           := create( 6.2725181549514408e-01, 3.0763585181253225e-17,
                      2.7068930273498138e-34);
  sin_t221 : constant triple_double
           := create( 6.2963823891492698e-01, 4.1115334049626806e-17,
                     -1.9167473580230747e-33);
  sin_t222 : constant triple_double
           := create( 6.3201873593980906e-01,-4.0164942296463612e-17,
                     -7.2208643641736723e-34);
  sin_t223 : constant triple_double
           := create( 6.3439328416364549e-01, 1.0420901929280035e-17,
                      4.1174558929280492e-34);
  sin_t224 : constant triple_double
           := create( 6.3676186123628420e-01, 3.1419048711901611e-17,
                     -2.2693738415126449e-33);
  sin_t225 : constant triple_double
           := create( 6.3912444486377573e-01, 1.2416796312271043e-17,
                     -6.2095419626356605e-34);
  sin_t226 : constant triple_double
           := create( 6.4148101280858316e-01,-9.9883430115943310e-18,
                      4.1969230376730128e-34);
  sin_t227 : constant triple_double
           := create( 6.4383154288979150e-01,-3.2084798795046886e-17,
                     -1.2595311907053305e-33);
  sin_t228 : constant triple_double
           := create( 6.4617601298331639e-01,-2.9756137382280815e-17,
                     -1.0275370077518259e-33);
  sin_t229 : constant triple_double
           := create( 6.4851440102211244e-01, 3.9870270313386831e-18,
                      1.9408388509540788e-34);
  sin_t230 : constant triple_double
           := create( 6.5084668499638088e-01, 3.9714670710500257e-17,
                      2.9178546787002963e-34);
  sin_t231 : constant triple_double
           := create( 6.5317284295377676e-01, 8.5695642060026238e-18,
                     -6.9165322305070633e-34);
  sin_t232 : constant triple_double
           := create( 6.5549285299961535e-01, 3.5638734426385005e-17,
                      1.2695365790889811e-33);
  sin_t233 : constant triple_double
           := create( 6.5780669329707864e-01, 1.9580943058468545e-17,
                     -1.1944272256627192e-33);
  sin_t234 : constant triple_double
           := create( 6.6011434206742048e-01,-1.3960054386823638e-19,
                      6.1515777931494047e-36);
  sin_t235 : constant triple_double
           := create( 6.6241577759017178e-01,-2.2615508885764591e-17,
                      5.0177050318126862e-34);
  sin_t236 : constant triple_double
           := create( 6.6471097820334490e-01,-3.6227793598034367e-17,
                     -9.0607934765540427e-34);
  sin_t237 : constant triple_double
           := create( 6.6699992230363747e-01, 3.5284364997428166e-17,
                     -1.0382057232458238e-33);
  sin_t238 : constant triple_double
           := create( 6.6928258834663612e-01,-5.4592652417447913e-17,
                     -2.5181014709695152e-33);
  sin_t239 : constant triple_double
           := create( 6.7155895484701844e-01,-4.0489037749296692e-17,
                      3.1995835625355681e-34);
  sin_t240 : constant triple_double
           := create( 6.7382900037875604e-01, 2.3091901236161086e-17,
                      5.7428037192881319e-34);
  sin_t241 : constant triple_double
           := create( 6.7609270357531592e-01, 3.7256902248049466e-17,
                      1.7059417895764375e-33);
  sin_t242 : constant triple_double
           := create( 6.7835004312986147e-01, 1.8302093041863122e-17,
                      9.5241675746813072e-34);
  sin_t243 : constant triple_double
           := create( 6.8060099779545302e-01, 2.8473293354522047e-17,
                      4.1331805977270903e-34);
  sin_t244 : constant triple_double
           := create( 6.8284554638524808e-01,-1.2958058061524531e-17,
                      1.8292386959330698e-34);
  sin_t245 : constant triple_double
           := create( 6.8508366777270036e-01, 2.5948135194645137e-17,
                     -8.5030743129500702e-34);
  sin_t246 : constant triple_double
           := create( 6.8731534089175916e-01,-5.5156158714917168e-17,
                      1.1896489854266829e-33);
  sin_t247 : constant triple_double
           := create( 6.8954054473706694e-01,-1.5889323294806790e-17,
                      9.1242356240205712e-34);
  sin_t248 : constant triple_double
           := create( 6.9175925836415775e-01, 2.7406078472410668e-17,
                      1.3286508943202092e-33);
  sin_t249 : constant triple_double
           := create( 6.9397146088965400e-01, 7.4345076956280137e-18,
                      7.5061528388197460e-34);
  sin_t250 : constant triple_double
           := create( 6.9617713149146299e-01,-4.1224081213582889e-17,
                     -3.1838716762083291e-35);
  sin_t251 : constant triple_double
           := create( 6.9837624940897280e-01, 4.8988282435667768e-17,
                      1.9134010413244152e-33);
  sin_t252 : constant triple_double
           := create( 7.0056879394324834e-01, 3.1027960192992922e-17,
                      9.5638250509179997e-34);
  sin_t253 : constant triple_double
           := create( 7.0275474445722530e-01, 2.5278294383629822e-18,
                     -8.6985561210674942e-35);
  sin_t254 : constant triple_double
           := create( 7.0493408037590488e-01, 2.7608725585748502e-17,
                      2.9816599471629137e-34);
  sin_t255 : constant triple_double
           := create( 7.0710678118654757e-01,-4.8336466567264567e-17,
                      2.0693376543497068e-33);

  cos_t0   : constant triple_double
           := create( 9.9999529380957619e-01,-1.9668064285322189e-17,
                     -6.3053955095883481e-34);
  cos_t1   : constant triple_double
           := create( 9.9998117528260111e-01, 3.3568103522895585e-17,
                     -1.4740132559368063e-35);
  cos_t2   : constant triple_double
           := create( 9.9995764455196390e-01,-3.1527836866647287e-17,
                      2.6363251186638437e-33);
  cos_t3   : constant triple_double
           := create( 9.9992470183914450e-01, 3.7931082512668012e-17,
                     -8.5099918660501484e-35);
  cos_t4   : constant triple_double
           := create( 9.9988234745421256e-01,-3.5477814872408538e-17,
                      1.7102001035303974e-33);
  cos_t5   : constant triple_double
           := create( 9.9983058179582340e-01, 1.8825140517551119e-17,
                     -5.1383513457616937e-34);
  cos_t6   : constant triple_double
           := create( 9.9976940535121528e-01, 4.2681177032289012e-17,
                      1.9062302359737099e-33);
  cos_t7   : constant triple_double
           := create( 9.9969881869620425e-01,-2.9851486403799753e-17,
                     -1.9084787370733737e-33);
  cos_t8   : constant triple_double
           := create( 9.9961882249517864e-01,-4.1181965521424734e-17,
                      2.0915365593699916e-33);
  cos_t9   : constant triple_double
           := create( 9.9952941750109314e-01, 2.0517917823755591e-17,
                     -4.7673802585706520e-34);
  cos_t10  : constant triple_double
           := create( 9.9943060455546173e-01, 3.9644497752257798e-17,
                     -2.3757223716722428e-34);
  cos_t11  : constant triple_double
           := create( 9.9932238458834954e-01,-4.2858538440845682e-17,
                      3.3235101605146565e-34);
  cos_t12  : constant triple_double
           := create( 9.9920475861836389e-01, 9.1796317110385693e-18,
                      5.5416208185868570e-34);
  cos_t13  : constant triple_double
           := create( 9.9907772775264536e-01, 2.1419007653587032e-17,
                     -7.9048203318529618e-34);
  cos_t14  : constant triple_double
           := create( 9.9894129318685687e-01,-2.0610641910058638e-17,
                     -1.2546525485913485e-33);
  cos_t15  : constant triple_double
           := create( 9.9879545620517241e-01,-1.2291693337075465e-17,
                      2.4468446786491271e-34);
  cos_t16  : constant triple_double
           := create( 9.9864021818026527e-01,-4.8690254312923302e-17,
                     -2.9470881967909147e-34);
  cos_t17  : constant triple_double
           := create( 9.9847558057329477e-01,-2.2002931182778795e-17,
                     -1.2371509454944992e-33);
  cos_t18  : constant triple_double
           := create( 9.9830154493389289e-01,-5.1869402702792278e-17,
                      1.0480195493633452e-33);
  cos_t19  : constant triple_double
           := create( 9.9811811290014918e-01, 2.7935487558113833e-17,
                      2.4423341255830345e-33);
  cos_t20  : constant triple_double
           := create( 9.9792528619859600e-01, 1.7143659778886362e-17,
                      5.7885840902887460e-34);
  cos_t21  : constant triple_double
           := create( 9.9772306664419164e-01,-2.6394475274898721e-17,
                     -1.6176223087661783e-34);
  cos_t22  : constant triple_double
           := create( 9.9751145614030345e-01, 5.6007205919806937e-18,
                     -5.9477673514685690e-35);
  cos_t23  : constant triple_double
           := create( 9.9729045667869021e-01, 9.1647695371101735e-18,
                      6.7824134309739296e-34);
  cos_t24  : constant triple_double
           := create( 9.9706007033948296e-01, 1.6734093546241963e-17,
                     -1.3169951440780028e-33);
  cos_t25  : constant triple_double
           := create( 9.9682029929116567e-01, 4.7062820708615655e-17,
                      2.8412041076474937e-33);
  cos_t26  : constant triple_double
           := create( 9.9657114579055484e-01, 1.1707179088390986e-17,
                     -7.5934413263024663e-34);
  cos_t27  : constant triple_double
           := create( 9.9631261218277800e-01, 1.1336497891624735e-17,
                      3.4002458674414360e-34);
  cos_t28  : constant triple_double
           := create( 9.9604470090125197e-01, 2.2870031707670695e-17,
                     -3.9184839405013148e-34);
  cos_t29  : constant triple_double
           := create( 9.9576741446765982e-01,-2.3151908323094359e-17,
                     -1.6306512931944591e-34);
  cos_t30  : constant triple_double
           := create( 9.9548075549192694e-01, 3.2084621412226554e-18,
                     -4.9501292146013023e-36);
  cos_t31  : constant triple_double
           := create( 9.9518472667219693e-01,-4.2486913678304410e-17,
                      1.3315510772504614e-33);
  cos_t32  : constant triple_double
           := create( 9.9487933079480562e-01, 4.2130813284943662e-18,
                     -4.2062597488288452e-35);
  cos_t33  : constant triple_double
           := create( 9.9456457073425542e-01, 3.6745069641528058e-17,
                     -3.0603304105471010e-33);
  cos_t34  : constant triple_double
           := create( 9.9424044945318790e-01, 4.4129423472462673e-17,
                     -3.0107231708238066e-33);
  cos_t35  : constant triple_double
           := create( 9.9390697000235606e-01,-1.8964849471123746e-17,
                     -1.5980853777937752e-35);
  cos_t36  : constant triple_double
           := create( 9.9356413552059530e-01, 2.9752309927797428e-17,
                     -4.5066707331134233e-34);
  cos_t37  : constant triple_double
           := create( 9.9321194923479450e-01, 3.3096906261272262e-17,
                      1.5592811973249567e-33);
  cos_t38  : constant triple_double
           := create( 9.9285041445986510e-01,-1.4094517733693302e-17,
                     -1.1954558131616916e-33);
  cos_t39  : constant triple_double
           := create( 9.9247953459870997e-01, 3.1093055095428906e-17,
                     -1.8379594757818019e-33);
  cos_t40  : constant triple_double
           := create( 9.9209931314219180e-01,-3.9431926149588778e-17,
                     -6.2758062911047230e-34);
  cos_t41  : constant triple_double
           := create( 9.9170975366909953e-01,-2.3372891311883661e-18,
                      2.7073298824968591e-35);
  cos_t42  : constant triple_double
           := create( 9.9131085984611544e-01,-2.5192111583372105e-17,
                     -1.2852471567380887e-33);
  cos_t43  : constant triple_double
           := create( 9.9090263542778001e-01, 1.5394565094566704e-17,
                     -1.0799984133184567e-33);
  cos_t44  : constant triple_double
           := create( 9.9048508425645709e-01,-5.5411437553780867e-17,
                     -1.4614017210753585e-33);
  cos_t45  : constant triple_double
           := create( 9.9005821026229712e-01,-1.7055485906233963e-17,
                      1.3454939685758777e-33);
  cos_t46  : constant triple_double
           := create( 9.8962201746320089e-01,-5.2398217968132530e-17,
                      1.3463189211456219e-33);
  cos_t47  : constant triple_double
           := create( 9.8917650996478101e-01,-4.0987309937047111e-17,
                     -4.4857560552048437e-34);
  cos_t48  : constant triple_double
           := create( 9.8872169196032378e-01,-1.0976227206656125e-17,
                      3.2311342577653764e-34);
  cos_t49  : constant triple_double
           := create( 9.8825756773074946e-01, 2.7030607784372632e-17,
                      7.7514866488601377e-35);
  cos_t50  : constant triple_double
           := create( 9.8778414164457218e-01,-2.3600693397159021e-17,
                     -1.2323283769707861e-33);
  cos_t51  : constant triple_double
           := create( 9.8730141815785843e-01,-5.2332261255715652e-17,
                     -2.7937644333152473e-33);
  cos_t52  : constant triple_double
           := create( 9.8680940181418553e-01,-5.0287214351061075e-17,
                     -2.2681526238144461e-33);
  cos_t53  : constant triple_double
           := create( 9.8630809724459867e-01,-2.1520877103013341e-17,
                      1.1866528054187716e-33);
  cos_t54  : constant triple_double
           := create( 9.8579750916756748e-01,-5.1439452979953012e-17,
                      2.6276439309996725e-33);
  cos_t55  : constant triple_double
           := create( 9.8527764238894122e-01, 2.3155637027900207e-17,
                     -7.5275971545764833e-34);
  cos_t56  : constant triple_double
           := create( 9.8474850180190421e-01, 1.0548144061829957e-17,
                      2.8786145266267306e-34);
  cos_t57  : constant triple_double
           := create( 9.8421009238692903e-01, 4.7983922627050691e-17,
                      2.2597419645070588e-34);
  cos_t58  : constant triple_double
           := create( 9.8366241921173025e-01, 1.9864948201635255e-17,
                     -1.0743046281211033e-35);
  cos_t59  : constant triple_double
           := create( 9.8310548743121629e-01, 4.2170007522888628e-17,
                      8.2396265656440904e-34);
  cos_t60  : constant triple_double
           := create( 9.8253930228744124e-01, 1.5149580813777224e-17,
                     -4.1802771422186237e-34);
  cos_t61  : constant triple_double
           := create( 9.8196386910955524e-01, 2.1108443711513084e-17,
                     -1.5253013442896054e-33);
  cos_t62  : constant triple_double
           := create( 9.8137919331375456e-01, 1.3428163260355633e-17,
                     -6.5294290469962986e-34);
  cos_t63  : constant triple_double
           := create( 9.8078528040323043e-01, 1.8546939997825006e-17,
                     -1.0696564445530757e-33);
  cos_t64  : constant triple_double
           := create( 9.8018213596811743e-01,-3.6801786963856159e-17,
                      6.3245171387992842e-34);
  cos_t65  : constant triple_double
           := create( 9.7956976568544052e-01, 1.5573991584990420e-17,
                     -1.3401066029782990e-33);
  cos_t66  : constant triple_double
           := create( 9.7894817531906220e-01,-2.3817727961148053e-18,
                     -1.0694750370381661e-34);
  cos_t67  : constant triple_double
           := create( 9.7831737071962765e-01,-2.1623082233344895e-17,
                      1.0970403012028032e-33);
  cos_t68  : constant triple_double
           := create( 9.7767735782450993e-01, 5.0514136167059628e-17,
                     -1.3254751701428788e-33);
  cos_t69  : constant triple_double
           := create( 9.7702814265775439e-01,-4.3353875751555997e-17,
                      5.4948839831535478e-34);
  cos_t70  : constant triple_double
           := create( 9.7636973133002114e-01, 9.3093931526213780e-18,
                     -4.1184949155685665e-34);
  cos_t71  : constant triple_double
           := create( 9.7570213003852857e-01,-2.5572556081259686e-17,
                     -9.3174244508942223e-34);
  cos_t72  : constant triple_double
           := create( 9.7502534506699412e-01, 2.6642660651899135e-17,
                      1.7819392739353853e-34);
  cos_t73  : constant triple_double
           := create( 9.7433938278557586e-01, 2.3041221476151512e-18,
                      1.0758686005031430e-34);
  cos_t74  : constant triple_double
           := create( 9.7364424965081198e-01,-5.1729808691005871e-17,
                     -1.5508473005989887e-33);
  cos_t75  : constant triple_double
           := create( 9.7293995220556018e-01,-3.1311211122281800e-17,
                     -2.6874087789006141e-33);
  cos_t76  : constant triple_double
           := create( 9.7222649707893627e-01, 3.6461169785938221e-17,
                      3.0309636883883133e-33);
  cos_t77  : constant triple_double
           := create( 9.7150389098625178e-01,-7.9865421122289046e-18,
                     -4.3628417211263380e-34);
  cos_t78  : constant triple_double
           := create( 9.7077214072895035e-01,-4.7992163325114922e-17,
                      3.0347528910975783e-33);
  cos_t79  : constant triple_double
           := create( 9.7003125319454397e-01, 1.8365300348428844e-17,
                     -1.4311097571944918e-33);
  cos_t80  : constant triple_double
           := create( 9.6928123535654853e-01,-4.5663660261927896e-17,
                      9.6147526917239387e-34);
  cos_t81  : constant triple_double
           := create( 9.6852209427441727e-01, 4.9475074918244771e-17,
                      2.8558738351911241e-33);
  cos_t82  : constant triple_double
           := create( 9.6775383709347551e-01,-4.5512132825515820e-17,
                     -1.4127617988719093e-33);
  cos_t83  : constant triple_double
           := create( 9.6697647104485207e-01, 3.8496228837337864e-17,
                     -5.3881631542745647e-34);
  cos_t84  : constant triple_double
           := create( 9.6619000344541250e-01, 5.1298840401665493e-17,
                      1.4564075904769808e-34);
  cos_t85  : constant triple_double
           := create( 9.6539444169768940e-01,-2.3745389918392156e-17,
                      5.9221515590053862e-34);
  cos_t86  : constant triple_double
           := create( 9.6458979328981276e-01,-3.4189470735959786e-17,
                      2.2982074155463522e-33);
  cos_t87  : constant triple_double
           := create( 9.6377606579543984e-01, 2.6463950561220029e-17,
                     -2.9073234590199323e-36);
  cos_t88  : constant triple_double
           := create( 9.6295326687368388e-01, 8.9341960404313634e-18,
                     -3.9071244661020126e-34);
  cos_t89  : constant triple_double
           := create( 9.6212140426904158e-01, 1.5236770453846305e-17,
                     -1.3050173525597142e-33);
  cos_t90  : constant triple_double
           := create( 9.6128048581132064e-01, 2.0933955216674039e-18,
                      1.0768607469015692e-34);
  cos_t91  : constant triple_double
           := create( 9.6043051941556579e-01, 2.4653904815317185e-17,
                     -1.3792169410906322e-33);
  cos_t92  : constant triple_double
           := create( 9.5957151308198452e-01, 1.1000640085000957e-17,
                     -4.2036030828223975e-34);
  cos_t93  : constant triple_double
           := create( 9.5870347489587160e-01,-4.3685014392372053e-17,
                      2.2001800662729131e-33);
  cos_t94  : constant triple_double
           := create( 9.5782641302753291e-01,-1.7696710075371263e-17,
                      1.9164034110382190e-34);
  cos_t95  : constant triple_double
           := create( 9.5694033573220882e-01, 4.0553869861875701e-17,
                     -1.7147013364302149e-33);
  cos_t96  : constant triple_double
           := create( 9.5604525134999641e-01, 3.7705045279589067e-17,
                      1.9678699997347571e-33);
  cos_t97  : constant triple_double
           := create( 9.5514116830577067e-01, 5.0088652955014668e-17,
                     -2.6983181838059211e-33);
  cos_t98  : constant triple_double
           := create( 9.5422809510910567e-01,-3.7545901690626874e-17,
                      1.4951619241257764e-33);
  cos_t99  : constant triple_double
           := create( 9.5330604035419386e-01,-2.5190738779919934e-17,
                     -1.4272239821134379e-33);
  cos_t100 : constant triple_double
           := create( 9.5237501271976588e-01,-2.0269300462299272e-17,
                     -1.0635956887246246e-33);
  cos_t101 : constant triple_double
           := create( 9.5143502096900834e-01, 3.1350584123266695e-17,
                     -2.4824833452737813e-33);
  cos_t102 : constant triple_double
           := create( 9.5048607394948170e-01, 1.9410097562630436e-17,
                     -8.1559393949816789e-34);
  cos_t103 : constant triple_double
           := create( 9.4952818059303667e-01,-7.5544151928043298e-18,
                     -5.1260245024046686e-34);
  cos_t104 : constant triple_double
           := create( 9.4856134991573027e-01, 2.0668262262333232e-17,
                     -5.9440730243667306e-34);
  cos_t105 : constant triple_double
           := create( 9.4758559101774109e-01, 4.3417993852125991e-17,
                     -2.7728667889840373e-34);
  cos_t106 : constant triple_double
           := create( 9.4660091308328353e-01, 3.5056800210680730e-17,
                      9.8578536940318117e-34);
  cos_t107 : constant triple_double
           := create( 9.4560732538052128e-01, 4.6019102478523738e-17,
                     -6.2534384769452059e-34);
  cos_t108 : constant triple_double
           := create( 9.4460483726148026e-01, 8.8100545476641165e-18,
                      5.2291695602757842e-34);
  cos_t109 : constant triple_double
           := create( 9.4359345816196039e-01,-2.4093127844404214e-17,
                      1.0283279856803939e-34);
  cos_t110 : constant triple_double
           := create( 9.4257319760144687e-01, 1.3235564806436886e-17,
                     -5.7048262885386911e-35);
  cos_t111 : constant triple_double
           := create( 9.4154406518302081e-01,-2.7896379547698341e-17,
                      1.6273236356733898e-33);
  cos_t112 : constant triple_double
           := create( 9.4050607059326830e-01, 2.8610421567116268e-17,
                      2.9261501147538827e-33);
  cos_t113 : constant triple_double
           := create( 9.3945922360218992e-01,-7.0152867943098655e-18,
                     -5.6395693818011210e-34);
  cos_t114 : constant triple_double
           := create( 9.3840353406310806e-01, 5.4242545044795490e-17,
                     -1.9039966607859759e-33);
  cos_t115 : constant triple_double
           := create( 9.3733901191257496e-01,-3.6570926284362776e-17,
                     -1.1902940071273247e-33);
  cos_t116 : constant triple_double
           := create( 9.3626566717027826e-01,-1.3013766145497654e-17,
                      5.2229870061990595e-34);
  cos_t117 : constant triple_double
           := create( 9.3518350993894761e-01,-3.2609395302485065e-17,
                     -8.1813015218875245e-34);
  cos_t118 : constant triple_double
           := create( 9.3409255040425887e-01, 4.4662824360767511e-17,
                     -2.5903243047396916e-33);
  cos_t119 : constant triple_double
           := create( 9.3299279883473885e-01, 4.2041415555384355e-17,
                      9.0285896495521276e-34);
  cos_t120 : constant triple_double
           := create( 9.3188426558166815e-01,-4.0785944377318095e-17,
                      1.7631450298754169e-33);
  cos_t121 : constant triple_double
           := create( 9.3076696107898371e-01, 1.9703775102838329e-17,
                      6.5657908718278205e-34);
  cos_t122 : constant triple_double
           := create( 9.2964089584318121e-01, 5.1282530016864107e-17,
                      2.3719739891916261e-34);
  cos_t123 : constant triple_double
           := create( 9.2850608047321559e-01,-2.3306639848485943e-17,
                     -7.7799084333208503e-34);
  cos_t124 : constant triple_double
           := create( 9.2736252565040111e-01,-2.7677111692155437e-17,
                      2.2110293450199576e-34);
  cos_t125 : constant triple_double
           := create( 9.2621024213831138e-01,-3.7303754586099054e-17,
                      2.0464457809993405e-33);
  cos_t126 : constant triple_double
           := create( 9.2504924078267758e-01, 6.0529447412576159e-18,
                     -8.8256517760278541e-35);
  cos_t127 : constant triple_double
           := create( 9.2387953251128674e-01, 1.7645047084336677e-17,
                     -5.0442537321586818e-34);
  cos_t128 : constant triple_double
           := create( 9.2270112833387852e-01, 5.2963798918539814e-17,
                     -5.7135699628876685e-34);
  cos_t129 : constant triple_double
           := create( 9.2151403934204190e-01, 4.1639843390684644e-17,
                      1.1891485604702356e-33);
  cos_t130 : constant triple_double
           := create( 9.2031827670911059e-01,-2.7806888779036837e-17,
                      2.7011013677071274e-33);
  cos_t131 : constant triple_double
           := create( 9.1911385169005777e-01,-2.6496484622344718e-17,
                      6.5403604763461920e-34);
  cos_t132 : constant triple_double
           := create( 9.1790077562139050e-01,-3.9074579680849515e-17,
                      2.3004636541490264e-33);
  cos_t133 : constant triple_double
           := create( 9.1667905992104270e-01,-4.1733978698287568e-17,
                      1.2094444804381172e-33);
  cos_t134 : constant triple_double
           := create( 9.1544871608826783e-01,-1.3591056692900894e-17,
                      5.9923027475594735e-34);
  cos_t135 : constant triple_double
           := create( 9.1420975570353069e-01,-3.6316182527814423e-17,
                     -1.9438819777122554e-33);
  cos_t136 : constant triple_double
           := create( 9.1296219042839821e-01,-4.7932505228039469e-17,
                     -1.7753551889428638e-33);
  cos_t137 : constant triple_double
           := create( 9.1170603200542988e-01,-2.6913273175034130e-17,
                     -5.1928101916162528e-35);
  cos_t138 : constant triple_double
           := create( 9.1044129225806725e-01,-5.0433041673313820e-17,
                      1.0938746257404305e-33);
  cos_t139 : constant triple_double
           := create( 9.0916798309052238e-01,-3.6878564091359894e-18,
                      2.9951330310507693e-34);
  cos_t140 : constant triple_double
           := create( 9.0788611648766626e-01,-4.9459964301225840e-17,
                     -1.6599682707075313e-33);
  cos_t141 : constant triple_double
           := create( 9.0659570451491533e-01, 3.0506718955442023e-17,
                     -1.4478836557141204e-33);
  cos_t142 : constant triple_double
           := create( 9.0529675931811882e-01,-4.1153099826889901e-17,
                      2.9859368705184223e-33);
  cos_t143 : constant triple_double
           := create( 9.0398929312344334e-01,-6.6097544687484308e-18,
                      1.2728013034680357e-34);
  cos_t144 : constant triple_double
           := create( 9.0267331823725883e-01,-1.9250787033961483e-17,
                      1.3242128993244527e-33);
  cos_t145 : constant triple_double
           := create( 9.0134884704602203e-01,-1.3524789367698682e-17,
                      6.3605353115880091e-34);
  cos_t146 : constant triple_double
           := create( 9.0001589201616028e-01,-5.0639618050802273e-17,
                      1.0783525384031576e-33);
  cos_t147 : constant triple_double
           := create( 8.9867446569395382e-01, 2.6316906461033013e-17,
                      3.7003137047796840e-35);
  cos_t148 : constant triple_double
           := create( 8.9732458070541832e-01,-3.6396283314867290e-17,
                     -2.3611649895474815e-33);
  cos_t149 : constant triple_double
           := create( 8.9596624975618511e-01, 4.9025099114811813e-17,
                     -1.9440489814795326e-33);
  cos_t150 : constant triple_double
           := create( 8.9459948563138270e-01,-1.7516226396814919e-17,
                     -1.3200670047246923e-33);
  cos_t151 : constant triple_double
           := create( 8.9322430119551532e-01,-4.1161239151908913e-18,
                      2.5380253805715999e-34);
  cos_t152 : constant triple_double
           := create( 8.9184070939234272e-01, 4.6690228137124547e-18,
                      1.6150254286841982e-34);
  cos_t153 : constant triple_double
           := create( 8.9044872324475788e-01, 1.1781931459051803e-17,
                     -1.3346142209571930e-34);
  cos_t154 : constant triple_double
           := create( 8.8904835585466457e-01,-1.1164514966766675e-17,
                     -3.4797636107798736e-34);
  cos_t155 : constant triple_double
           := create( 8.8763962040285393e-01, 1.2805091918587960e-17,
                      3.9948742059584459e-35);
  cos_t156 : constant triple_double
           := create( 8.8622253014888064e-01,-6.7307369600274315e-18,
                      1.2385593432917413e-34);
  cos_t157 : constant triple_double
           := create( 8.8479709843093779e-01,-9.4331469628972690e-18,
                     -5.7106541478701439e-34);
  cos_t158 : constant triple_double
           := create( 8.8336333866573158e-01, 1.5822643380255127e-17,
                     -7.8921320007588250e-34);
  cos_t159 : constant triple_double
           := create( 8.8192126434835505e-01,-1.9843248405890562e-17,
                     -7.0412114007673834e-34);
  cos_t160 : constant triple_double
           := create( 8.8047088905216075e-01, 1.6311096602996350e-17,
                     -5.7541360594724172e-34);
  cos_t161 : constant triple_double
           := create( 8.7901222642863353e-01,-4.7356837291118011e-17,
                      1.4388771297975192e-33);
  cos_t162 : constant triple_double
           := create( 8.7754529020726124e-01, 5.0113311846499550e-17,
                      2.8382769008739543e-34);
  cos_t163 : constant triple_double
           := create( 8.7607009419540660e-01, 5.8729024235147677e-18,
                      2.7941144391738458e-34);
  cos_t164 : constant triple_double
           := create( 8.7458665227817611e-01,-5.7216617730397065e-19,
                     -2.9705811503689596e-35);
  cos_t165 : constant triple_double
           := create( 8.7309497841829009e-01, 7.8424672990129903e-18,
                     -4.8685015839797165e-34);
  cos_t166 : constant triple_double
           := create( 8.7159508665595109e-01,-5.5272998038551050e-17,
                     -2.2104090204984907e-33);
  cos_t167 : constant triple_double
           := create( 8.7008699110871146e-01,-4.1888510868549968e-17,
                      7.0900185861878415e-34);
  cos_t168 : constant triple_double
           := create( 8.6857070597134090e-01, 2.7192781689782903e-19,
                     -1.6710140396932428e-35);
  cos_t169 : constant triple_double
           := create( 8.6704624551569265e-01, 3.0267859550930567e-18,
                     -1.1559438782171572e-34);
  cos_t170 : constant triple_double
           := create( 8.6551362409056909e-01,-6.3723113549628899e-18,
                      2.3725520321746832e-34);
  cos_t171 : constant triple_double
           := create( 8.6397285612158670e-01, 4.1486355957361607e-17,
                      2.2709976932210266e-33);
  cos_t172 : constant triple_double
           := create( 8.6242395611104050e-01, 3.7008992527383130e-17,
                      5.2128411542701573e-34);
  cos_t173 : constant triple_double
           := create( 8.6086693863776731e-01,-3.0050048898573656e-17,
                     -8.8706183090892111e-34);
  cos_t174 : constant triple_double
           := create( 8.5930181835700836e-01, 4.2435655816850687e-17,
                      7.6181814059912025e-34);
  cos_t175 : constant triple_double
           := create( 8.5772861000027212e-01,-4.8183447936336620e-17,
                     -1.1044130517687532e-33);
  cos_t176 : constant triple_double
           := create( 8.5614732837519447e-01, 9.1806925616606261e-18,
                      5.6328649785951470e-34);
  cos_t177 : constant triple_double
           := create( 8.5455798836540053e-01,-1.2991124236396092e-17,
                      1.2893407722948080e-34);
  cos_t178 : constant triple_double
           := create( 8.5296060493036363e-01, 2.7152984251981370e-17,
                      7.4336483283120719e-34);
  cos_t179 : constant triple_double
           := create( 8.5135519310526520e-01,-5.3279874446016209e-17,
                      2.2281156380919942e-33);
  cos_t180 : constant triple_double
           := create( 8.4974176800085244e-01, 5.1812347659974015e-17,
                      3.0810626087331275e-33);
  cos_t181 : constant triple_double
           := create( 8.4812034480329723e-01, 1.8762563415239981e-17,
                      1.4048773307919617e-33);
  cos_t182 : constant triple_double
           := create( 8.4649093877405213e-01,-4.7969419958569345e-17,
                     -2.7518267097886703e-33);
  cos_t183 : constant triple_double
           := create( 8.4485356524970712e-01,-4.3631360296879637e-17,
                     -2.0307726853367547e-33);
  cos_t184 : constant triple_double
           := create( 8.4320823964184544e-01, 9.6536707005959077e-19,
                      2.8995142431556364e-36);
  cos_t185 : constant triple_double
           := create( 8.4155497743689844e-01,-3.4095465391321557e-17,
                     -8.4130208607579595e-34);
  cos_t186 : constant triple_double
           := create( 8.3989379419599952e-01,-1.6673694881511411e-17,
                     -1.4759184141750289e-33);
  cos_t187 : constant triple_double
           := create( 8.3822470555483808e-01,-3.5560085052855026e-17,
                      1.1689791577022643e-33);
  cos_t188 : constant triple_double
           := create( 8.3654772722351201e-01,-2.0899059027066533e-17,
                     -9.8104097821002585e-35);
  cos_t189 : constant triple_double
           := create( 8.3486287498638001e-01, 4.6048430609159657e-17,
                     -5.1827423265239912e-34);
  cos_t190 : constant triple_double
           := create( 8.3317016470191319e-01, 1.3275129507229764e-18,
                      4.8589164115370863e-35);
  cos_t191 : constant triple_double
           := create( 8.3146961230254524e-01, 1.4073856984728024e-18,
                      4.6951315383980830e-35);
  cos_t192 : constant triple_double
           := create( 8.2976123379452305e-01,-2.9349109376485597e-18,
                      1.1496917934149818e-34);
  cos_t193 : constant triple_double
           := create( 8.2804504525775580e-01,-4.4196593225871532e-17,
                      2.7967864855211251e-33);
  cos_t194 : constant triple_double
           := create( 8.2632106284566353e-01,-5.3957485453612902e-17,
                      6.8976896130138550e-34);
  cos_t195 : constant triple_double
           := create( 8.2458930278502529e-01,-2.6512360488868275e-17,
                      1.6916964350914386e-34);
  cos_t196 : constant triple_double
           := create( 8.2284978137582632e-01, 1.5193019034505495e-17,
                      9.6890547246521685e-34);
  cos_t197 : constant triple_double
           := create( 8.2110251499110465e-01, 3.0715131609697682e-17,
                     -1.7037168325855879e-33);
  cos_t198 : constant triple_double
           := create( 8.1934752007679701e-01,-4.8200736995191133e-17,
                     -1.5574489646672781e-35);
  cos_t199 : constant triple_double
           := create( 8.1758481315158371e-01,-1.4883149812426772e-17,
                     -7.8273262771298917e-34);
  cos_t200 : constant triple_double
           := create( 8.1581441080673378e-01, 8.2652693782130871e-18,
                     -2.3028778135179471e-34);
  cos_t201 : constant triple_double
           := create( 8.1403632970594841e-01,-5.2127351877042624e-17,
                     -1.9047670611316360e-33);
  cos_t202 : constant triple_double
           := create( 8.1225058658520388e-01, 3.1054545609214803e-17,
                      2.2649541922707251e-34);
  cos_t203 : constant triple_double
           := create( 8.1045719825259477e-01, 2.3520367349840499e-17,
                     -7.7530070904846341e-34);
  cos_t204 : constant triple_double
           := create( 8.0865618158817498e-01, 9.3251597879721674e-18,
                     -7.1823301933068394e-34);
  cos_t205 : constant triple_double
           := create( 8.0684755354379922e-01, 4.9220603766095546e-17,
                      2.9796016899903487e-33);
  cos_t206 : constant triple_double
           := create( 8.0503133114296355e-01, 5.1368289568212149e-17,
                      6.3082807402256524e-34);
  cos_t207 : constant triple_double
           := create( 8.0320753148064494e-01,-3.3060609804814910e-17,
                     -1.2242726252420433e-33);
  cos_t208 : constant triple_double
           := create( 8.0137617172314024e-01,-2.0958013413495834e-17,
                     -4.3798162198006931e-34);
  cos_t209 : constant triple_double
           := create( 7.9953726910790501e-01, 2.0356723822005431e-17,
                     -9.7448513696896360e-34);
  cos_t210 : constant triple_double
           := create( 7.9769084094339116e-01,-4.6730759884788944e-17,
                      2.3075897077191757e-33);
  cos_t211 : constant triple_double
           := create( 7.9583690460888357e-01,-3.0062724851910721e-17,
                     -2.2496210832042235e-33);
  cos_t212 : constant triple_double
           := create( 7.9397547755433717e-01,-7.4194631759921416e-18,
                      2.4124341304631069e-34);
  cos_t213 : constant triple_double
           := create( 7.9210657730021239e-01,-3.7087850202326467e-17,
                     -1.4874457267228264e-33);
  cos_t214 : constant triple_double
           := create( 7.9023022143731003e-01, 2.3056905954954492e-17,
                      1.4481080533260193e-33);
  cos_t215 : constant triple_double
           := create( 7.8834642762660623e-01, 3.4396993154059708e-17,
                      1.7710623746737170e-33);
  cos_t216 : constant triple_double
           := create( 7.8645521359908577e-01,-9.7841429939305265e-18,
                      3.3906063272445472e-34);
  cos_t217 : constant triple_double
           := create( 7.8455659715557524e-01,-8.5627965423173476e-18,
                     -2.1106834459001849e-34);
  cos_t218 : constant triple_double
           := create( 7.8265059616657573e-01, 9.0745866975808825e-18,
                      6.7623847404278666e-34);
  cos_t219 : constant triple_double
           := create( 7.8073722857209449e-01,-9.9198782066678806e-18,
                     -2.1265794012162715e-36);
  cos_t220 : constant triple_double
           := create( 7.7881651238147598e-01,-2.4891385579973807e-17,
                      6.7665497024807980e-35);
  cos_t221 : constant triple_double
           := create( 7.7688846567323244e-01, 7.7418602570672864e-18,
                     -5.9986517872157897e-34);
  cos_t222 : constant triple_double
           := create( 7.7495310659487393e-01,-5.2209083189826433e-17,
                     -9.6653593393686612e-34);
  cos_t223 : constant triple_double
           := create( 7.7301045336273699e-01,-3.2565907033649772e-17,
                      1.3860807251523929e-33);
  cos_t224 : constant triple_double
           := create( 7.7106052426181382e-01,-4.4558442347769265e-17,
                     -2.9863565614083783e-33);
  cos_t225 : constant triple_double
           := create( 7.6910333764557959e-01, 5.1546455184564817e-17,
                      2.6142829553524292e-33);
  cos_t226 : constant triple_double
           := create( 7.6713891193582040e-01,-1.8885903683750782e-17,
                     -1.3659359331495433e-33);
  cos_t227 : constant triple_double
           := create( 7.6516726562245896e-01,-3.2707225612534598e-17,
                      1.1177117747079528e-33);
  cos_t228 : constant triple_double
           := create( 7.6318841726338127e-01, 2.6314748416750748e-18,
                      1.4048039063095910e-34);
  cos_t229 : constant triple_double
           := create( 7.6120238548426178e-01, 3.5315510881690551e-17,
                      1.2833566381864357e-33);
  cos_t230 : constant triple_double
           := create( 7.5920918897838807e-01,-3.8558842175523123e-17,
                      2.9720241208332759e-34);
  cos_t231 : constant triple_double
           := create( 7.5720884650648457e-01,-1.9909098777335502e-17,
                      3.9409283266158482e-34);
  cos_t232 : constant triple_double
           := create( 7.5520137689653655e-01,-1.9402238001823017e-17,
                     -3.7756206444727573e-34);
  cos_t233 : constant triple_double
           := create( 7.5318679904361252e-01,-3.7937789838736540e-17,
                     -6.7009539920231559e-34);
  cos_t234 : constant triple_double
           := create( 7.5116513190968637e-01, 4.3499761158645868e-17,
                      2.5227718971102212e-33);
  cos_t235 : constant triple_double
           := create( 7.4913639452345937e-01,-4.4729078447011889e-17,
                     -2.4206025249983768e-33);
  cos_t236 : constant triple_double
           := create( 7.4710060598018013e-01, 1.1874824875965430e-17,
                      2.1992523849833518e-34);
  cos_t237 : constant triple_double
           := create( 7.4505778544146595e-01, 1.5078686911877863e-17,
                      8.0898987212942471e-34);
  cos_t238 : constant triple_double
           := create( 7.4300795213512172e-01,-2.5144629669719265e-17,
                      7.1128989512526157e-34);
  cos_t239 : constant triple_double
           := create( 7.4095112535495911e-01,-1.4708616952297345e-17,
                     -4.9550433827142032e-34);
  cos_t240 : constant triple_double
           := create( 7.3888732446061511e-01, 3.4324874808225091e-17,
                     -1.3706639444717610e-33);
  cos_t241 : constant triple_double
           := create( 7.3681656887736990e-01,-2.8932468101656295e-17,
                     -3.4649887126202378e-34);
  cos_t242 : constant triple_double
           := create( 7.3473887809596350e-01,-3.4507595976263941e-17,
                     -2.3718000676666409e-33);
  cos_t243 : constant triple_double
           := create( 7.3265427167241282e-01, 1.8918673481573520e-17,
                     -1.5123719544119886e-33);
  cos_t244 : constant triple_double
           := create( 7.3056276922782759e-01,-2.9689959904476928e-17,
                     -1.1276871244239744e-33);
  cos_t245 : constant triple_double
           := create( 7.2846439044822520e-01, 1.1924642323370718e-19,
                      5.9001892316611011e-36);
  cos_t246 : constant triple_double
           := create( 7.2635915508434601e-01,-3.1917502443460542e-17,
                      7.7047912412039396e-34);
  cos_t247 : constant triple_double
           := create( 7.2424708295146689e-01, 2.9198471334403004e-17,
                      2.3027324968739464e-33);
  cos_t248 : constant triple_double
           := create( 7.2212819392921535e-01,-2.3871262053452047e-17,
                      1.0636125432862273e-33);
  cos_t249 : constant triple_double
           := create( 7.2000250796138165e-01,-2.5689658854462333e-17,
                     -9.1492566948567925e-34);
  cos_t250 : constant triple_double
           := create( 7.1787004505573171e-01, 2.7006476062511453e-17,
                     -2.2854956580215348e-34);
  cos_t251 : constant triple_double
           := create( 7.1573082528381871e-01,-5.1581018476410262e-17,
                     -1.3736271349300259e-34);
  cos_t252 : constant triple_double
           := create( 7.1358486878079364e-01,-4.2342504403133584e-17,
                     -4.2690366101617268e-34);
  cos_t253 : constant triple_double
           := create( 7.1143219574521643e-01, 7.9643298613856813e-18,
                      2.9488239510721469e-34);
  cos_t254 : constant triple_double
           := create( 7.0927282643886569e-01,-3.7597359110245730e-17,
                      1.0613125954645119e-34);
  cos_t255 : constant triple_double
           := create( 7.0710678118654757e-01,-4.8336466567264567e-17,
                      2.0693376543497068e-33);

  sin_table : array(0..255) of triple_double 
            := (sin_t0  ,sin_t1  ,sin_t2  ,sin_t3  ,sin_t4  ,sin_t5  ,
                sin_t6  ,sin_t7  ,sin_t8  ,sin_t9  ,sin_t10 ,sin_t11 ,
                sin_t12 ,sin_t13 ,sin_t14 ,sin_t15 ,sin_t16 ,sin_t17 ,
                sin_t18 ,sin_t19 ,sin_t20 ,sin_t21 ,sin_t22 ,sin_t23 ,
                sin_t24 ,sin_t25 ,sin_t26 ,sin_t27 ,sin_t28 ,sin_t29 ,
                sin_t30 ,sin_t31 ,sin_t32 ,sin_t33 ,sin_t34 ,sin_t35 ,
                sin_t36 ,sin_t37 ,sin_t38 ,sin_t39 ,sin_t40 ,sin_t41 ,
                sin_t42 ,sin_t43 ,sin_t44 ,sin_t45 ,sin_t46 ,sin_t47 ,
                sin_t48 ,sin_t49 ,sin_t50 ,sin_t51 ,sin_t52 ,sin_t53 ,
                sin_t54 ,sin_t55 ,sin_t56 ,sin_t57 ,sin_t58 ,sin_t59 ,
                sin_t60 ,sin_t61 ,sin_t62 ,sin_t63 ,sin_t64 ,sin_t65 ,
                sin_t66 ,sin_t67 ,sin_t68 ,sin_t69 ,sin_t70 ,sin_t71 ,
                sin_t72 ,sin_t73 ,sin_t74 ,sin_t75 ,sin_t76 ,sin_t77 ,
                sin_t78 ,sin_t79 ,sin_t80 ,sin_t81 ,sin_t82 ,sin_t83 ,
                sin_t84 ,sin_t85 ,sin_t86 ,sin_t87 ,sin_t88 ,sin_t89 ,
                sin_t90 ,sin_t91 ,sin_t92 ,sin_t93 ,sin_t94 ,sin_t95 ,
                sin_t96 ,sin_t97 ,sin_t98 ,sin_t99 ,sin_t100,sin_t101,
                sin_t102,sin_t103,sin_t104,sin_t105,sin_t106,sin_t107,
                sin_t108,sin_t109,sin_t110,sin_t111,sin_t112,sin_t113,
                sin_t114,sin_t115,sin_t116,sin_t117,sin_t118,sin_t119,
                sin_t120,sin_t121,sin_t122,sin_t123,sin_t124,sin_t125,
                sin_t126,sin_t127,sin_t128,sin_t129,sin_t130,sin_t131,
                sin_t132,sin_t133,sin_t134,sin_t135,sin_t136,sin_t137,
                sin_t138,sin_t139,sin_t140,sin_t141,sin_t142,sin_t143,
                sin_t144,sin_t145,sin_t146,sin_t147,sin_t148,sin_t149,
                sin_t150,sin_t151,sin_t152,sin_t153,sin_t154,sin_t155,
                sin_t156,sin_t157,sin_t158,sin_t159,sin_t160,sin_t161,
                sin_t162,sin_t163,sin_t164,sin_t165,sin_t166,sin_t167,
                sin_t168,sin_t169,sin_t170,sin_t171,sin_t172,sin_t173,
                sin_t174,sin_t175,sin_t176,sin_t177,sin_t178,sin_t179,
                sin_t180,sin_t181,sin_t182,sin_t183,sin_t184,sin_t185,
                sin_t186,sin_t187,sin_t188,sin_t189,sin_t190,sin_t191,
                sin_t192,sin_t193,sin_t194,sin_t195,sin_t196,sin_t197,
                sin_t198,sin_t199,sin_t200,sin_t201,sin_t202,sin_t203,
                sin_t204,sin_t205,sin_t206,sin_t207,sin_t208,sin_t209,
                sin_t210,sin_t211,sin_t212,sin_t213,sin_t214,sin_t215,
                sin_t216,sin_t217,sin_t218,sin_t219,sin_t220,sin_t221,
                sin_t222,sin_t223,sin_t224,sin_t225,sin_t226,sin_t227,
                sin_t228,sin_t229,sin_t230,sin_t231,sin_t232,sin_t233,
                sin_t234,sin_t235,sin_t236,sin_t237,sin_t238,sin_t239,
                sin_t240,sin_t241,sin_t242,sin_t243,sin_t244,sin_t245,
                sin_t246,sin_t247,sin_t248,sin_t249,sin_t250,sin_t251,
                sin_t252,sin_t253,sin_t254,sin_t255);

  cos_table : array(0..255) of triple_double
            := (cos_t0  ,cos_t1  ,cos_t2  ,cos_t3  ,cos_t4  ,cos_t5  ,
                cos_t6  ,cos_t7  ,cos_t8  ,cos_t9  ,cos_t10 ,cos_t11 ,
                cos_t12 ,cos_t13 ,cos_t14 ,cos_t15 ,cos_t16 ,cos_t17 ,
                cos_t18 ,cos_t19 ,cos_t20 ,cos_t21 ,cos_t22 ,cos_t23 ,
                cos_t24 ,cos_t25 ,cos_t26 ,cos_t27 ,cos_t28 ,cos_t29 ,
                cos_t30 ,cos_t31 ,cos_t32 ,cos_t33 ,cos_t34 ,cos_t35 ,
                cos_t36 ,cos_t37 ,cos_t38 ,cos_t39 ,cos_t40 ,cos_t41 ,
                cos_t42 ,cos_t43 ,cos_t44 ,cos_t45 ,cos_t46 ,cos_t47 ,
                cos_t48 ,cos_t49 ,cos_t50 ,cos_t51 ,cos_t52 ,cos_t53 ,
                cos_t54 ,cos_t55 ,cos_t56 ,cos_t57 ,cos_t58 ,cos_t59 ,
                cos_t60 ,cos_t61 ,cos_t62 ,cos_t63 ,cos_t64 ,cos_t65 ,
                cos_t66 ,cos_t67 ,cos_t68 ,cos_t69 ,cos_t70 ,cos_t71 ,
                cos_t72 ,cos_t73 ,cos_t74 ,cos_t75 ,cos_t76 ,cos_t77 ,
                cos_t78 ,cos_t79 ,cos_t80 ,cos_t81 ,cos_t82 ,cos_t83 ,
                cos_t84 ,cos_t85 ,cos_t86 ,cos_t87 ,cos_t88 ,cos_t89 ,
                cos_t90 ,cos_t91 ,cos_t92 ,cos_t93 ,cos_t94 ,cos_t95 ,
                cos_t96 ,cos_t97 ,cos_t98 ,cos_t99 ,cos_t100,cos_t101,
                cos_t102,cos_t103,cos_t104,cos_t105,cos_t106,cos_t107,
                cos_t108,cos_t109,cos_t110,cos_t111,cos_t112,cos_t113,
                cos_t114,cos_t115,cos_t116,cos_t117,cos_t118,cos_t119,
                cos_t120,cos_t121,cos_t122,cos_t123,cos_t124,cos_t125,
                cos_t126,cos_t127,cos_t128,cos_t129,cos_t130,cos_t131,
                cos_t132,cos_t133,cos_t134,cos_t135,cos_t136,cos_t137,
                cos_t138,cos_t139,cos_t140,cos_t141,cos_t142,cos_t143,
                cos_t144,cos_t145,cos_t146,cos_t147,cos_t148,cos_t149,
                cos_t150,cos_t151,cos_t152,cos_t153,cos_t154,cos_t155,
                cos_t156,cos_t157,cos_t158,cos_t159,cos_t160,cos_t161,
                cos_t162,cos_t163,cos_t164,cos_t165,cos_t166,cos_t167,
                cos_t168,cos_t169,cos_t170,cos_t171,cos_t172,cos_t173,
                cos_t174,cos_t175,cos_t176,cos_t177,cos_t178,cos_t179,
                cos_t180,cos_t181,cos_t182,cos_t183,cos_t184,cos_t185,
                cos_t186,cos_t187,cos_t188,cos_t189,cos_t190,cos_t191,
                cos_t192,cos_t193,cos_t194,cos_t195,cos_t196,cos_t197,
                cos_t198,cos_t199,cos_t200,cos_t201,cos_t202,cos_t203,
                cos_t204,cos_t205,cos_t206,cos_t207,cos_t208,cos_t209,
                cos_t210,cos_t211,cos_t212,cos_t213,cos_t214,cos_t215,
                cos_t216,cos_t217,cos_t218,cos_t219,cos_t220,cos_t221,
                cos_t222,cos_t223,cos_t224,cos_t225,cos_t226,cos_t227,
                cos_t228,cos_t229,cos_t230,cos_t231,cos_t232,cos_t233,
                cos_t234,cos_t235,cos_t236,cos_t237,cos_t238,cos_t239,
                cos_t240,cos_t241,cos_t242,cos_t243,cos_t244,cos_t245,
                cos_t246,cos_t247,cos_t248,cos_t249,cos_t250,cos_t251,
                cos_t252,cos_t253,cos_t254,cos_t255);

end Triple_Double_Constants;
