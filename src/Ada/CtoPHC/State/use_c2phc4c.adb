with text_io;                           use text_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Greeting_Banners;
with PHCpack_Operations;
with C_to_PHCpack;
with use_syscon,use_syspool;
with use_tabform;
with use_solcon,use_solpool;
with use_scaling;
with use_reduction;
with use_c2pieri,use_c2lrhom;
with use_c2fac,use_c2mbt;
with use_roco;
with use_celcon;
with use_track,use_sweep;
with use_mapcon;
with use_nxtsol;
with unisolve,use_giftwrap;
with use_numbtrop;
with use_series;
with use_padcon;
with use_multip;  -- multiplicity structure
with use_witsols;
with Job_Containers;
with Job_Handlers;
with Multprec_PolySys_Interface;
with Multprec_Solutions_Interface;
with File_Management_Interface;
with Symbol_Table_Interface;
with Continuation_Parameters_Interface;
with Newton_Interface;
with Deflation_Interface;
with Path_Trackers_Interface;
with Witness_Interface;

function use_c2phc4c ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  procedure Write_Welcome is
  begin
    put_line(Greeting_Banners.welcome);
  end Write_Welcome;

  function Handle_Jobs return integer32 is

    use Job_Containers;
    use Job_Handlers;
    use Symbol_Table_Interface;
    use Multprec_PolySys_Interface;
    use Multprec_Solutions_Interface;
    use File_Management_Interface;
    use Newton_Interface;
    use Deflation_Interface;
    use Continuation_Parameters_Interface;
    use Path_Trackers_Interface;
    use Witness_Interface;

  begin
    if vrblvl > 0
     then put_line("-> in use_c2phc4c.Handle_Jobs ...");
    end if;
    case job is
      when 0 => Write_Welcome; return 0;
      when 1 => return Standard_Target_Poly_System_to_Container(vrblvl-1);
      when 2 => return Standard_Container_Poly_System_to_Target(vrblvl-1);
      when 3 => return Standard_Start_Poly_System_to_Container(vrblvl-1);
      when 4 => return Standard_Container_Poly_System_to_Start(vrblvl-1);
      when 5 => return Standard_Target_Solutions_to_Container(vrblvl-1);
      when 6 => return Standard_Container_Solutions_to_Target(vrblvl-1);
      when 7 => return Standard_Start_Solutions_to_Container(vrblvl-1);
      when 8 => return Standard_Container_Solutions_to_Start(vrblvl-1);
      when 9 => return Newton_Standard_Polynomial_Verify(vrblvl-1);
      when 10..15 => return C_to_PHCpack(job-10,0,vrblvl-1);
      when 16 => return Path_Trackers_Standard_Polynomial_Solve(a,vrblvl-1);
      when 17..19 => return C_to_PHCpack(job-10,0,vrblvl-1);
      when 20..29 => return use_syscon(job-20,a,b,c,vrblvl-1);
      when 30..38 => return use_solcon(job-30,a,b,c,vrblvl-1);
      when 39 => return use_c2fac(28,a,b,c,vrblvl-1); -- set state to silent
      when 40..65 => return use_c2fac(job-40,a,b,c,vrblvl-1);
      when 66 => return Witness_Standard_Polynomial_Embed(a,vrblvl-1);
      when 67 => return use_syscon(67,a,b,c,vrblvl-1); -- load poly as string
      when 68 => return use_c2fac(job-42,a,b,c,vrblvl-1); -- return #factors
      when 69 => return use_c2fac(job-42,a,b,c,vrblvl-1); -- return factor
      when 70 => return Continuation_Parameters_Ask_Values(vrblvl-1);
      when 71 => return Continuation_Parameters_Ask_Output_Level(vrblvl-1);
      when 72 => return Continuation_Parameters_Get_All(c,vrblvl-1);
      when 73 => return Continuation_Parameters_Set_All(c,vrblvl-1);
     -- store Laurential as string
      when 74 => return use_syscon(74,a,b,c,vrblvl-1);
      when 75 => return Standard_Laurent_Solver(a,b,vrblvl-1);
     -- store polynomial as string
      when 76 => return use_syscon(76,a,b,c,vrblvl-1);
      when 77 => return Standard_Polynomial_Solver(a,b,vrblvl-1);
      when 78 => return Mixed_Volume(a,vrblvl-1);
      when 79 => return Stable_Mixed_Volume(a,b,vrblvl-1);
      when 80..105 => return use_celcon(job-80,a,b,c,vrblvl-1);
     -- load dobldobl poly as string
      when 106 => return use_syscon(68,a,b,c,vrblvl-1);
     -- load quaddobl poly as string
      when 107 => return use_syscon(69,a,b,c,vrblvl-1);
     -- load multprec poly as string
      when 108 => return use_syscon(70,a,b,c,vrblvl-1);
     -- random system in container
      when 109 => return use_syscon(71,a,b,c,vrblvl-1);
      when 110..118 => return use_roco(job-110,a,b,c);
     -- degree of standard polynomial
      when 119 => return use_syscon(20,a,b,c,vrblvl-1);
     -- operations on Laurent systems :
      when 120..127 => return use_syscon(job-20,a,b,c,vrblvl-1);
      when 128 => return use_syscon(77,a,b,c,vrblvl-1); -- load Laur as string
      when 129 => return Witness_DoblDobl_Polynomial_Embed(a,vrblvl-1);
      when 130..145 => return use_solcon(job-120,a,b,c,vrblvl-1);
     -- drop coordinate by name
      when 146 => return use_solcon(9,a,b,c,vrblvl-1);
      when 147 => return use_syscon(10,a,b,c,vrblvl-1);
      when 148 => return use_syscon(11,a,b,c,vrblvl-1);
      when 149..171 => return use_track(job-150,a,b,c,vrblvl-1);
     -- track operations for double double precision :
      when 172..178 => return use_track(job-150,a,b,c,vrblvl-1);
     -- variable precision Newton step :
      when 179 => return Newton_Varbprec_Step(a,b,vrblvl-1);
     -- double double and quad double L-R homotopies :
      when 180..181 => return use_c2lrhom(job-178,a,b,c,vrblvl-1);
     -- track operations for quad double precision :
      when 182..188 => return use_track(job-150,a,b,c,vrblvl-1);
     -- tuning continuation parameters, deflation, and Newton step
      when 189 => return Continuation_Parameters_Get_Value(a,c,vrblvl-1);
      when 190 => return Continuation_Parameters_Set_Value(a,c,vrblvl-1);
      when 191 => return File_Management_Set_Output(a,b,vrblvl-1);
      when 192 => return File_Management_Close_Output(vrblvl-1);
      when 193 => return Continuation_Parameters_Autotune(a,b,vrblvl-1);
      when 194 => return Continuation_Parameters_Show(vrblvl-1);
      when 195 => return Newton_Multprec_Polynomial_Step(a,vrblvl-1);
      when 196 => return Deflation_Standard_Run(a,b,c,vrblvl-1);
      when 197 => return Newton_QuadDobl_Polynomial_Step(vrblvl-1);
      when 198 => return Newton_DoblDobl_Polynomial_Step(vrblvl-1);
      when 199 => return Newton_Standard_Polynomial_Step(vrblvl-1);
      when 200..209 => return use_solcon(job-170,a,b,c,vrblvl-1);
      when 210..227 => return use_c2pieri(job-210,a,b,c,vrblvl-1);
      when 228..229 => return use_c2lrhom(job-228,a,b,c,vrblvl-1);
      when 230 => return use_track(42,a,b,c,vrblvl-1);
      when 231..235 => return C_to_PHCpack(job-220,0,vrblvl-1);
      when 236 => return Path_Trackers_DoblDobl_Polynomial_Solve(a,vrblvl-1);
      when 237..238 => return C_to_PHCpack(job-220,0,vrblvl-1);
      when 239 => return use_celcon(46,a,b,c,vrblvl-1);
      when 240 => return use_celcon(47,a,b,c,vrblvl-1);
      when 241..245 => return C_to_PHCpack(job-220,0,vrblvl-1);
      when 246 => return Path_Trackers_QuadDobl_Polynomial_Solve(a,vrblvl-1);
      when 247..248 => return C_to_PHCpack(job-220,0,vrblvl-1);
     -- deflation in double double and quad double precision
      when 249 => return Deflation_DoblDobl_Run(a,b,c,vrblvl-1);
      when 250 => return Deflation_QuadDobl_Run(a,b,c,vrblvl-1);
     -- double double versions for jobs 1 to 8
      when 251 => return DoblDobl_Target_Poly_System_to_Container(vrblvl-1);
      when 252 => return DoblDobl_Container_Poly_System_to_Target(vrblvl-1);
      when 253 => return DoblDobl_Start_Poly_System_to_Container(vrblvl-1);
      when 254 => return DoblDobl_Container_Poly_System_to_Start(vrblvl-1);
      when 255 => return DoblDobl_Target_Solutions_to_Container(vrblvl-1);
      when 256 => return DoblDobl_Container_Solutions_to_Target(vrblvl-1);
      when 257 => return DoblDobl_Start_Solutions_to_Container(vrblvl-1);
      when 258 => return DoblDobl_Container_Solutions_to_Start(vrblvl-1);
     -- double double witness set for a hypersurface
      when 259 => return use_track(49,a,b,c,vrblvl-1);
      when 260 => return Witness_QuadDobl_Polynomial_Embed(a,vrblvl-1);
     -- quad double versions for jobs 1 to 8
      when 261 => return QuadDobl_Target_Poly_System_to_Container(vrblvl-1);
      when 262 => return QuadDobl_Container_Poly_System_to_Target(vrblvl-1);
      when 263 => return QuadDobl_Start_Poly_System_to_Container(vrblvl-1);
      when 264 => return QuadDobl_Container_Poly_System_to_Start(vrblvl-1);
      when 265 => return QuadDobl_Target_Solutions_to_Container(vrblvl-1);
      when 266 => return QuadDobl_Container_Solutions_to_Target(vrblvl-1);
      when 267 => return QuadDobl_Start_Solutions_to_Container(vrblvl-1);
      when 268 => return QuadDobl_Container_Solutions_to_Start(vrblvl-1);
     -- quad double witness set for a hypersurface
      when 269 => return use_track(50,a,b,c,vrblvl-1);
     -- interface to diagonal homotopies ...
      when 270 => return use_track(40,a,b,c,vrblvl-1); -- st witset of poly
      when 271 => return use_track(41,a,b,c,vrblvl-1); -- start diagl csc sols
     -- univariate polynomial solvers
      when 272 => return unisolve(1,a,b,c,vrblvl-1); -- double precision
      when 273 => return unisolve(2,a,b,c,vrblvl-1); -- double double precision
      when 274 => return unisolve(3,a,b,c,vrblvl-1); -- quad double precision
      when 275 => return unisolve(4,a,b,c,vrblvl-1); -- multiprecision
     -- read next of solutions
      when 276 => return use_solcon(276,a,b,c,vrblvl-1); -- next standard
      when 277 => return use_solcon(277,a,b,c,vrblvl-1); -- next double double
      when 278 => return use_solcon(278,a,b,c,vrblvl-1); -- next quad double
      when 279 => return use_solcon(279,a,b,c,vrblvl-1); -- next multprec
      when 280 => return use_c2fac(29,a,b,c,vrblvl-1); -- st rnd cmplx nbr
     -- multiprecision versions for jobs 1 to 8
      when 281 => return Multprec_Target_Poly_System_to_Container(vrblvl-1);
      when 282 => return Multprec_Container_Poly_System_to_Target(vrblvl-1);
      when 283 => return Multprec_Start_Poly_System_to_Container(vrblvl-1);
      when 284 => return Multprec_Container_Poly_System_to_Start(vrblvl-1);
      when 285 => return Multprec_Target_Solutions_to_Container(vrblvl-1);
      when 286 => return Multprec_Container_Solutions_to_Target(vrblvl-1);
      when 287 => return Multprec_Start_Solutions_to_Container(vrblvl-1);
      when 288 => return Multprec_Container_Solutions_to_Start(vrblvl-1);
     -- diagonal homotopy in double double and quad double precision
      when 289 => return use_track(43,a,b,c,vrblvl-1); -- dd diagonal homotopy
      when 290 => return use_track(44,a,b,c,vrblvl-1); -- qd diagonal homotopy
     -- manipulation of symbols
      when 291 => return Symbol_Table_Remove_by_Index(a,vrblvl-1);
      when 292 => return Symbol_Table_Sort_Embedded(a,vrblvl-1);
      when 293 => return Symbol_Table_Size(a,vrblvl-1);
      when 294 => return Symbol_Table_Write(vrblvl-1);
      when 295 => return Symbol_Table_String(a,b,vrblvl-1);
      when 296 => return Symbol_Table_Remove_by_Name(a,b,vrblvl-1);
     -- interface to diagonal homotopies continued
      when 297 => return use_track(45,a,b,c,vrblvl-1); -- dd diag startsols
      when 298 => return use_track(46,a,b,c,vrblvl-1); -- qd diag startsols
      when 299 => return use_track(47,a,b,c,vrblvl-1); -- dd collapse diagonal
      when 300..305 => return use_syspool(job-300,a,b,c,vrblvl-1);
      when 306..311 => return use_syscon(job-294,a,b,c,vrblvl-1);
      when 312 => return use_track(48,a,b,c,vrblvl-1); -- qd collapse diagonal
      when 313..317 => return use_syspool(job-307,a,b,c,vrblvl-1);
      when 318 => return use_syspool(11,a,b,c,vrblvl-1); -- init dd sys pool
      when 319 => return use_syspool(12,a,b,c,vrblvl-1); -- init qd sys pool
      when 320..325 => return use_solpool(job-320,a,b,c,vrblvl-1);
     -- one Newton step on Laurent system :
      when 326 => return Newton_Standard_Laurent_Step(vrblvl-1);
      when 327 => return Newton_DoblDobl_Laurent_Step(vrblvl-1);
      when 328 => return Newton_QuadDobl_Laurent_Step(vrblvl-1);
      when 329 => return Newton_Multprec_Laurent_Step(a,vrblvl-1);
     -- operations on double double system container
      when 330..339 => return use_syscon(job-130,a,b,c,vrblvl-1);
     -- operations on double double solution container
      when 340..349 => return use_solcon(job-300,a,b,c,vrblvl-1);
      when 370..371 => return use_solcon(job-300,a,b,c,vrblvl-1);
      when 378 => return use_solcon(job-300,a,b,c,vrblvl-1);
     -- operations on quad double system container
      when 380..389 => return use_syscon(job-170,a,b,c,vrblvl-1);
     -- operations on quad double solution container
      when 390..399 => return use_solcon(job-310,a,b,c,vrblvl-1);
      when 420..421 => return use_solcon(job-310,a,b,c,vrblvl-1);
      when 428 => return use_solcon(job-310,a,b,c,vrblvl-1);
     -- operations on monomial maps as solutions to binomial systems
      when 430..438 => return use_mapcon(job-430,a,b,c,vrblvl-1);
     -- scan for the number of variables
      when 439 => return Symbol_Table_Scan(a,b,vrblvl-1);
     -- operations on multiprecision system container
      when 440..444 => return use_syscon(job-220,a,b,c,vrblvl-1);
      when 447..449 => return use_syscon(job-220,a,b,c,vrblvl-1);
     -- operations on multiprecision solutions :
      when 450..453 => return use_solcon(job-330,a,b,c,vrblvl-1);
      when 457 => return use_solcon(job-330,a,b,c,vrblvl-1);
     -- moving pointer to the current solution
      when 454 => return use_solcon(300,a,b,c,vrblvl-1);
      when 455 => return use_solcon(301,a,b,c,vrblvl-1);
      when 456 => return use_solcon(302,a,b,c,vrblvl-1);
      when 458 => return use_solcon(303,a,b,c,vrblvl-1);
     -- polyhedral homotopies in double double precision :
      when 460..469 => return use_celcon(job-434,a,b,c,vrblvl-1);
     -- polyhedral homotopies in quad double precision :
      when 470..479 => return use_celcon(job-434,a,b,c,vrblvl-1);
     -- string representations of multiprecision solutions :
      when 480..481 => return use_solcon(job-330,a,b,c,vrblvl-1);
      when 488 => return use_solcon(job-330,a,b,c,vrblvl-1);
     -- PHCpack operations for multiprecision arithmetic
      when 491 => return Multprec_PolySys_Prompt_for_Target(a,vrblvl-1);
      when 492 => return Multprec_PolySys_Write_Target(vrblvl-1);
      when 493 => return Multprec_PolySys_Prompt_for_Start(a,vrblvl-1);
      when 494 => return Multprec_PolySys_Write_Start(vrblvl-1);
      when 495 => return Multprec_Solutions_Write_Start(vrblvl-1);
      when 496 => return Path_Trackers_Multprec_Polynomial_Solve(a,vrblvl-1);
      when 497 => return Multprec_Solutions_Write_Target(vrblvl-1);
      when 498 => PHCpack_Operations.Multprec_Clear; return 0;
     -- path trackers with generators :
      when 500..520 => return use_nxtsol(job-500,a,b,c,vrblvl-1);
     -- multiprecision homotopies :
      when 522..524 => return use_track(job-470,a,b,c,vrblvl-1);
     -- get length of current solution string :
      when 525 => return use_solcon(304,a,b,c,vrblvl-1);
      when 526 => return use_solcon(305,a,b,c,vrblvl-1);
      when 527 => return use_solcon(306,a,b,c,vrblvl-1);
      when 528 => return use_solcon(307,a,b,c,vrblvl-1);
     -- multihomogeneous Bezout numbers and start systems
      when 530..532 => return use_roco(job-520,a,b,c);
     -- returns current solution string :
      when 533 => return use_solcon(308,a,b,c,vrblvl-1);
      when 534 => return use_solcon(309,a,b,c,vrblvl-1);
      when 535 => return use_solcon(310,a,b,c,vrblvl-1);
      when 536 => return use_solcon(311,a,b,c,vrblvl-1);
     -- homotopy membership tests
      when 537 => return use_c2mbt(0,a,b,c,vrblvl-1); -- standard membertest
      when 538 => return use_c2mbt(1,a,b,c,vrblvl-1); -- dobldobl membertest
      when 539 => return use_c2mbt(2,a,b,c,vrblvl-1); -- quaddobl membertest
     -- operations to read systems into the containers
      when 540..543 => return use_syscon(job,a,b,c,vrblvl-1);
     -- operations to read systems and solutions into the containers
      when 544..547 => return use_solcon(job,a,b,c);
     -- random dobldobl and quaddobl systems
      when 548 => return use_syscon(78,a,b,c,vrblvl-1);
      when 549 => return use_syscon(79,a,b,c,vrblvl-1);
     -- operations on Laurent container for double doubles :
      when 550..558 => return use_syscon(job-440,a,b,c,vrblvl-1);
      when 559 => return use_syscon(72,a,b,c,vrblvl-1);
     -- operations on Laurent container for quad doubles :
      when 560..568 => return use_syscon(job-440,a,b,c,vrblvl-1);
      when 569 => return use_syscon(73,a,b,c,vrblvl-1);
     -- operations on Laurent container for multiprecision :
      when 570..574 => return use_syscon(job-440,a,b,c,vrblvl-1);
      when 577..579 => return use_syscon(job-440,a,b,c,vrblvl-1);
     -- convex hull via giftwrapping :
      when 580..589 => return use_giftwrap(job-579,a,b,c,vrblvl-1);
     -- scaling systems and solutions :
      when 590..596 => return use_scaling(job-589,a,b,c,vrblvl-1);
     -- copy start solutions from cell container to solutions container
      when 597 => return use_celcon(48,a,b,c,vrblvl-1); -- st start solution
      when 598 => return use_celcon(49,a,b,c,vrblvl-1); -- dd start solution
      when 599 => return use_celcon(50,a,b,c,vrblvl-1); -- qd start solution
     -- size limits of string representations of polynomials
      when 600..607 => return use_syscon(job-520,a,b,c);
     -- make system in the dobldobl and quaddobl systems pool
      when 608 => return use_syspool(16,a,b,c,vrblvl-1); -- k-th dobldobl sys
      when 609 => return use_syspool(17,a,b,c,vrblvl-1); -- k-th quaddobl sys
     -- run the sweep homotopy :
      when 610..621 => return use_sweep(job-610,a,b,c,vrblvl-1);
     -- crude path trackers :
      when 622 => return use_track(55,a,b,c,vrblvl-1); -- sd crude tracker
      when 623 => return use_track(56,a,b,c,vrblvl-1); -- dd crude tracker
      when 624 => return use_track(57,a,b,c,vrblvl-1); -- qd crude tracker
     -- embedding of Laurent systems :
      when 625 => return Witness_Standard_Laurent_Embed(a,vrblvl-1);
      when 626 => return Witness_DoblDobl_Laurent_Embed(a,vrblvl-1);
      when 627 => return Witness_QuadDobl_Laurent_Embed(a,vrblvl-1);
     -- make standard monodromy breakup verbose
      when 630 => return use_c2fac(30,a,b,c,vrblvl-1);
     -- monodromy breakup in double double precision :
      when 631..649 => return use_c2fac(job-600,a,b,c,vrblvl-1);
      when 652..660 => return use_c2fac(job-600,a,b,c,vrblvl-1);
     -- monodromy breakup in quad double precision :
      when 661..679 => return use_c2fac(job-600,a,b,c,vrblvl-1);
      when 682..690 => return use_c2fac(job-600,a,b,c,vrblvl-1);
     -- power series Newton method
      when 691..696 => return use_series(job-690,a,b,c,vrblvl-1);
     -- clear systems pool
      when 697 => return use_syspool(13,a,b,c,vrblvl-1); -- clear st syspool
      when 698 => return use_syspool(14,a,b,c,vrblvl-1); -- clear dd syspool
      when 699 => return use_syspool(15,a,b,c,vrblvl-1); -- clear qd syspool
     -- blackbox solvers in double double and quad double precision
      when 700 => return DoblDobl_Polynomial_Solver(a,b,vrblvl-1);
      when 701 => return DoblDobl_Laurent_Solver(a,b,vrblvl-1);
      when 702 => return QuadDobl_Polynomial_Solver(a,b,vrblvl-1);
      when 703 => return QuadDobl_Laurent_Solver(a,b,vrblvl-1);
     -- Pade approximants
      when 704 => return use_series(7,a,b,c,vrblvl-1); -- in double precision
      when 705 => return use_series(8,a,b,c,vrblvl-1); -- double doubles
      when 706 => return use_series(9,a,b,c,vrblvl-1); -- quad doubles
     -- reduction of polynomial systems
      when 707 => return use_reduction(1,a,b,c,vrblvl-1); -- standard linear
      when 708 => return use_reduction(2,a,b,c,vrblvl-1); -- dobldobl linear
      when 709 => return use_reduction(3,a,b,c,vrblvl-1); -- quaddobl linear
      when 710 => return use_reduction(4,a,b,c,vrblvl-1); -- standard nonlinear
     -- container for numerically computed tropisms
      when 711..731 => return use_numbtrop(job-710,a,b,c,vrblvl-1);
     -- computation of multiplicity structure
      when 732 => return use_multip(0,a,b,c,vrblvl-1); -- double precision
      when 733 => return use_multip(0,a,b,c,vrblvl-1); -- with double doubles
      when 734 => return use_multip(0,a,b,c,vrblvl-1); -- with quad doubles
     -- pade continuation
      when 735 => return use_padcon(0,a,b,c,vrblvl-1); -- set default values
      when 736 => return use_padcon(1,a,b,c,vrblvl-1); -- clear parameter vals
      when 737 => return use_padcon(2,a,b,c,vrblvl-1); -- get a parameter value
      when 738 => return use_padcon(3,a,b,c,vrblvl-1); -- set a parameter value
      when 739 => return use_padcon(4,a,b,c,vrblvl-1); -- track paths
      when 740 => return use_padcon(23,a,b,c,vrblvl-1); -- reset parameters
     -- integer mixed cell configurations
      when 741..758 => return use_celcon(job-690,a,b,c,vrblvl-1);
     -- reading, writing Laurent start and target systems
      when 759..773 => return c_to_phcpack(job-730,0,vrblvl-1);
     -- solve by Laurent homotopy continuation
      when 774 => return Path_Trackers_Standard_Laurent_Solve(a,vrblvl-1);
      when 775 => return Path_Trackers_DoblDobl_Laurent_Solve(a,vrblvl-1);
      when 776 => return Path_Trackers_QuadDobl_Laurent_Solve(a,vrblvl-1);
     -- copying Laurent systems from and into the containers
      when 777 => return Standard_Container_Laur_System_to_Start(vrblvl-1);
      when 778 => return DoblDobl_Container_Laur_System_to_Start(vrblvl-1);
      when 779 => return QuadDobl_Container_Laur_System_to_Start(vrblvl-1);
      when 780 => return Standard_Container_Laur_System_to_Target(vrblvl-1);
      when 781 => return DoblDobl_Container_Laur_System_to_Target(vrblvl-1);
      when 782 => return QuadDobl_Container_Laur_System_to_Target(vrblvl-1);
      when 783 => return Standard_Start_Laur_System_to_Container(vrblvl-1);
      when 784 => return DoblDobl_Start_Laur_System_to_Container(vrblvl-1);
      when 785 => return QuadDobl_Start_Laur_System_to_Container(vrblvl-1);
      when 786 => return Standard_Target_Laur_System_to_Container(vrblvl-1);
      when 787 => return DoblDobl_Target_Laur_System_to_Container(vrblvl-1);
      when 788 => return QuadDobl_Target_Laur_System_to_Container(vrblvl-1);
     -- cascades for Laurent homotopies
      when 789 => return use_track(58,a,b,c,vrblvl-1); -- st csc Laur htpy
      when 790 => return use_track(59,a,b,c,vrblvl-1); -- dd csc Laur htpy
      when 791 => return use_track(60,a,b,c,vrblvl-1); -- qd csc Laur htpy
      when 792 => return Path_Trackers_Standard_Laurent_Homotopy(vrblvl-1);
      when 793 => return Path_Trackers_DoblDobl_Laurent_Homotopy(vrblvl-1);
      when 794 => return Path_Trackers_QuadDobl_Laurent_Homotopy(vrblvl-1);
     -- homotopy membership tests on Laurent systems
      when 795 => return use_c2mbt(3,a,b,c,vrblvl-1); -- Laurent st membertest
      when 796 => return use_c2mbt(4,a,b,c,vrblvl-1); -- Laurent dd membertest
      when 797 => return use_c2mbt(5,a,b,c,vrblvl-1); -- Laurent qd membertest
     -- read a witness set defined by a Laurent polynomial system
      when 798 => return use_c2fac(91,a,b,c,vrblvl-1); -- prompt st Laur witset
      when 799 => return use_c2fac(92,a,b,c,vrblvl-1); -- prompt dd Laur witset
      when 800 => return use_c2fac(93,a,b,c,vrblvl-1); -- prompt qd Laur witset
      when 801 => return use_c2fac(94,a,b,c,vrblvl-1); -- read st Laur witset
      when 802 => return use_c2fac(95,a,b,c,vrblvl-1); -- read dd Laur witset
      when 803 => return use_c2fac(96,a,b,c,vrblvl-1); -- read qd Laur witset
     -- monodromy factorization on witness sets defined by Laurent systems
      when 804 => return use_c2fac(97,a,b,c,vrblvl-1); -- init st Laur sampler
      when 805 => return use_c2fac(98,a,b,c,vrblvl-1); -- init dd Laur sampler
      when 806 => return use_c2fac(99,a,b,c,vrblvl-1); -- init qd Laur sampler
     -- copy embedded system from sampler to container
      when 807 => return use_c2fac(100,a,b,c,vrblvl-1); -- sd Laurent copy
      when 808 => return use_c2fac(101,a,b,c,vrblvl-1); -- dd Laurent copy
      when 809 => return use_c2fac(102,a,b,c,vrblvl-1); -- qd Laurent copy
     -- construct a diagonal Laurent homotopy
      when 810 => return use_track(61,a,b,c,vrblvl-1); -- st diag Laur htp
      when 811 => return use_track(62,a,b,c,vrblvl-1); -- dd diag Laur htp
      when 812 => return use_track(63,a,b,c,vrblvl-1); -- qd diag Laur htp
     -- witness sets for Laurent polynomials
      when 813 => return use_track(64,a,b,c,vrblvl-1); -- st witset Laurent
      when 814 => return use_track(65,a,b,c,vrblvl-1); -- dd witset Laurent
      when 815 => return use_track(66,a,b,c,vrblvl-1); -- qd witset Laurent
     -- swap slack variables to the end
      when 816 => return Witness_Standard_Polynomial_Swap(a,b,vrblvl-1);
      when 817 => return Witness_DoblDobl_Polynomial_Swap(a,b,vrblvl-1);
      when 818 => return Witness_QuadDobl_Polynomial_Swap(a,b,vrblvl-1);
      when 819 => return Witness_Standard_Laurent_Swap(a,b,vrblvl-1);
      when 820 => return Witness_DoblDobl_Laurent_Swap(a,b,vrblvl-1);
      when 821 => return Witness_QuadDobl_Laurent_Swap(a,b,vrblvl-1);
     -- homotopy membership tests with symbolic test points
      when 822 => return use_c2mbt(6,a,b,c,vrblvl-1); -- standard membertest
      when 823 => return use_c2mbt(7,a,b,c,vrblvl-1); -- dobldobl membertest
      when 824 => return use_c2mbt(8,a,b,c,vrblvl-1); -- quaddobl membertest
      when 825 => return use_c2mbt(9,a,b,c,vrblvl-1);  -- st Laur membertest
      when 826 => return use_c2mbt(10,a,b,c,vrblvl-1); -- dd Laur membertest
      when 827 => return use_c2mbt(11,a,b,c,vrblvl-1); -- qd Laur membertest
     -- dropping a variable from Laurent polynomial systems
      when 828 => return use_syscon(22,a,b,c,vrblvl-1); -- st Laurent by idx
      when 829 => return use_syscon(23,a,b,c,vrblvl-1); -- dd Laurent by idx
      when 830 => return use_syscon(24,a,b,c,vrblvl-1); -- qd Laurent by idx
      when 831 => return use_syscon(25,a,b,c,vrblvl-1); -- st Laurent by name
      when 832 => return use_syscon(26,a,b,c,vrblvl-1); -- dd Laurent by name
      when 833 => return use_syscon(27,a,b,c,vrblvl-1); -- qd Laurent by name
     -- extract DEMiCs output data
     -- when 834 => return use_outdata(0,a,b,c); -- allocate memory for lifting
     -- when 835 => return use_outdata(1,a,b,c); -- assign a lifting value
     -- when 836 => return use_outdata(2,a,b,c); -- retrieve a lifting value
     -- when 837 => return use_outdata(3,a,b,c); -- clear lifting values
     -- when 838 => return use_outdata(4,a,b,c); -- append cell indices
     -- when 839 => return use_outdata(5,a,b,c); -- retrieve cell indices
     -- when 840 => return use_outdata(6,a,b,c); -- clear cell indices
     -- when 841 => return use_outdata(7,a,b,c); -- store mixed volume
     -- when 842 => return use_outdata(8,a,b,c); -- retrieve mixed volume
     -- when 843 => return use_outdata(9,a,b,c); -- call DEMiCs for mixed volume
     -- when 844 => return use_outdata(10,a,b,c); -- stable mv by DEMiCs
     -- numerical irreducible decomposition
      when 845..859 => return use_witsols(job-845,a,b,c,vrblvl-1); -- solvers
     -- for the irreducible factor string, use job 993
     -- Pade continuation in a step wise fashion
      when 860 => return use_padcon(5,a,b,c,vrblvl-1); -- init homotopy
      when 861 => return use_padcon(6,a,b,c,vrblvl-1); -- set start solution
      when 862 => return use_padcon(7,a,b,c,vrblvl-1); -- do next step
      when 863 => return use_padcon(8,a,b,c,vrblvl-1); -- get current solution
      when 864 => return use_padcon(9,a,b,c,vrblvl-1); -- clear data
      when 865 => return use_padcon(10,a,b,c,vrblvl-1); -- get pole radius
      when 866 => return use_padcon(11,a,b,c,vrblvl-1); -- get closest radius
      when 867 => return use_padcon(12,a,b,c,vrblvl-1); -- get t value
      when 868 => return use_padcon(13,a,b,c,vrblvl-1); -- get step size
      when 869 => return use_padcon(14,a,b,c,vrblvl-1); -- series coefficient
      when 870 => return use_padcon(15,a,b,c,vrblvl-1); -- get Pade coefficient
      when 871 => return use_padcon(16,a,b,c,vrblvl-1); -- get pole
     -- reading dobldobl and quaddobl target systems without solutions
      when 872 => return use_track(67,a,b,c,vrblvl-1); -- read dobldobl target
      when 873 => return use_track(68,a,b,c,vrblvl-1); -- read quaddobl target
     -- write homotopy continuation paramaters to defined output file
      when 874 => return use_padcon(17,a,b,c,vrblvl-1); 
     -- set value of the continuation parameter to zero
      when 875..877 => return use_solcon(job,a,b,c,vrblvl-1);
     -- initializes natural parameter homotopy in series-Pade tracker
      when 878 => return use_padcon(18,a,b,c,vrblvl-1);
     -- functions for stable mixed cells
      when 879 => return use_celcon(69,a,b,c,vrblvl-1);
      when 880 => return use_celcon(70,a,b,c,vrblvl-1);
      when 881 => return use_celcon(71,a,b,c,vrblvl-1);
      when 882 => return use_celcon(72,a,b,c,vrblvl-1); -- solve st stable
      when 883 => return use_celcon(73,a,b,c,vrblvl-1); -- solve dd stable
      when 884 => return use_celcon(74,a,b,c,vrblvl-1); -- solve qd stable
     -- retrieving step sizes for the step-by-step series-Pade tracker
      when 885 => return use_padcon(19,a,b,c,vrblvl-1); -- get series step
      when 886 => return use_padcon(20,a,b,c,vrblvl-1); -- get pole step
      when 887 => return use_padcon(21,a,b,c,vrblvl-1); -- get distance eta
      when 888 => return use_padcon(22,a,b,c,vrblvl-1); -- get Hessian step
     -- passing the tableau forms to the systems containers
      when 889 => return use_tabform(0,a,b,c); -- store standard tableau form
      when 890 => return use_tabform(1,a,b,c); -- get std tableau dimensions
     -- projective transformations of systems
      when 891..893 => return use_syscon(job,a,b,c,vrblvl-1); -- 1-homogeneous
      when 901..903 => return use_syscon(job,a,b,c,vrblvl-1); -- to affine
     -- multi-homogeneous
      when 904..906 => return use_syscon(job,a,b,c,vrblvl-1);
     -- from m-hom to affine
      when 907..909 => return use_syscon(job,a,b,c,vrblvl-1);
     -- projective transformations of solutions
      when 894..896 => return use_solcon(job,a,b,c,vrblvl-1); -- 1-homogeneous
     -- adding a symbol passed as string
      when 897 => return use_syscon(job,a,b,c,vrblvl-1);
     -- projective transformations of solutions continued ...
      when 898..900 => return use_solcon(job,a,b,c,vrblvl-1); -- to affine
     -- multi-homogeneous
      when 910..912 => return use_solcon(job,a,b,c,vrblvl-1);
     -- from m-hom to affine
      when 913..915 => return use_solcon(job,a,b,c,vrblvl-1);
     -- reading solutions from file with given name
      when 916..918 => return use_solcon(job,a,b,c,vrblvl-1);
      when 919 => return use_padcon(24,a,b,c,vrblvl-1);
      when 920 => return Job_Handlers.Standard_Condition_Report(a,b,c,vrblvl-1);
     -- making and clearing Laurent homotopies
      when 921 => return use_track(69,a,b,c,vrblvl-1);
      when 922 => return use_track(70,a,b,c,vrblvl-1);
      when 923 => return use_track(71,a,b,c,vrblvl-1);
      when 924 => return use_track(72,a,b,c,vrblvl-1);
      when 925 => return use_track(73,a,b,c,vrblvl-1);
      when 926 => return use_track(74,a,b,c,vrblvl-1);
     -- get string representation of irreducible factors
      when 993 => return use_witsols(15,a,b,c,vrblvl-1);
     -- get number of available CPUs
      when 994 => return Job_Handlers.Get_Core_Count(a,vrblvl-1);
     -- get and set the gamma constant
      when 995 => return Job_Handlers.Get_Gamma_Constant(a,c,vrblvl-1);
      when 996 => return Job_Handlers.Set_Gamma_Constant(a,c,vrblvl-1);
     -- get and set the seed and the version string
      when 997 => return Job_Handlers.Get_Seed(a,vrblvl-1);
      when 998 => return Job_Handlers.Set_Seed(a,vrblvl-1);
      when 999 => return Job_Handlers.Version_String(a,b,vrblvl-1);
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_c2phc handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 raise; -- return job;
end use_c2phc4c;
