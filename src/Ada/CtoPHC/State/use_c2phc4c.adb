with Interfaces.C;
with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Parse_Dimensions;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;
with Standard_Root_Refiners;
-- with Standard_Continuation_Data_io;     use Standard_Continuation_Data_io;
with Standard_Deflation_Methods;
with DoblDobl_Deflation_Methods;
with QuadDobl_Deflation_Methods;
with Verification_of_Solutions;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
-- with Assignments_of_Solutions;          use Assignments_of_Solutions;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Multprec_Solutions_Container;
with PHCpack_Operations;
with PHCpack_Operations_io;
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
with Symbol_Table_Interface;
with Newton_Interface;
with Continuation_Parameters_Interface;
with Witness_Interface;

function use_c2phc4c ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer;
                       vrblvl : integer32 := 0 ) return integer32 is

  procedure Write_Menu is
  begin
    new_line;
    put_line("General MENU to use PHCpack from C :");
    put_line("999. returns in the version string of PHCpack;");
    put_line("998. fixes the seed of the random number generator;");
    put_line("  0. display this menu;");
    put_line("MENU to move data from PHCpack from and to containers :");
    put_line("  1. copy target system to systems container;");
    put_line("  2. copy system in container to target system;");
    put_line("  3. copy start system to systems container;");
    put_line("  4. copy system in container to start system;");
    put_line("  5. copy target solutions to solutions container;");
    put_line("  6. copy solutions in container to target solutions;");
    put_line("  7. copy start solutions to solutions container;");
    put_line("  8. copy solutions in container to start solutions;");
    put_line("  9. verify solutions container using systems container.");
    put_line("MENU to the plain C to PHCpack interface :");
    put_line(" 10. displays the menu of available operations;");
    put_line(" 11. read target polynomial system;");
    put_line(" 12. write target polynomial system;");
    put_line(" 13. read start polynomial system;");
    put_line(" 14. write start polynomial system;");
    put_line(" 15. write start solutions;");
    put_line(" 16. solve by homotopy continuation;");
    put_line(" 17. write the target solutions;");
    put_line(" 18. clear the data in PHCpack_Operations;");
    put_line(" 19. define the output file.");
    put_line("191. define the output file from a string.");
    put_line("192. close the defined output file.");
    put_line("199. do one Newton step on container input data.");
    put_line("MENU to the operations in the systems container :");
    put_line(" 20. read polynomial system and put in container;");
    put_line(" 21. write the polynomial system in the container;");
    put_line(" 22. return the dimension of the polynomial system;");
    put_line(" 23. initializes the container with a dimension;"); 
    put_line(" 24. return the number of terms in a polynomial;");
    put_line(" 25. return a term in a polynomial;");
    put_line(" 26. add a term to a polynomial;");
    put_line(" 27. the systems container is cleared;");
    put_line(" 28. returns the total degree of the system;");
    put_line(" 29. clears the symbol table.");
    put_line("120. read Laurent polynomial system and put in container;");
    put_line("121. write the Laurent polynomial system in the container;");
    put_line("MENU to the operations in the solutions container :");
    put_line(" 30. read solutions from file and put in container;");
    put_line(" 31. write solutions in the container;");
    put_line(" 32. return the length of the container;");
    put_line(" 33. return the dimension of the solution vectors;");
    put_line(" 34. return a solution from the container;");
    put_line(" 35. changes a solution in the container;");
    put_line(" 36. append a solution to the container;");
    put_line(" 37. clears all solutions from the container;");
    put_line("MENU to the monodromy factorization :");
    put_line(" 40. displays menu of the factorization operations;");
    put_line(" 41. reads a witness set respresented as embedded system;");
    put_line(" 42. initialize the sampling machine;");
    put_line(" 43. assign coefficient of a slice;");
    put_line(" 44. store a random gamma constant in PHCpack;");
    put_line(" 45. compute a new witness set on the new slices;");
    put_line(" 46. swaps slices and solution sets to turn back;");
    put_line(" 47. copy embedded system from sampler to systems container;");
    put_line(" 48. copy original solutions from sampler to container;");
    put_line(" 49. put solutions from monodromy grid in container;");
    put_line(" 50. initialize maximum number of monodromy loops;");
    put_line(" 51. copy solutions in container to Monodromy_Permutations;");
    put_line(" 52. compute permutation by last stored solution list;");
    put_line(" 53. update the breakup with a new permutation;");
    put_line(" 54. write the current monodromy breakup;");
    put_line(" 55. applies the linear trace to certify the decomposition.");
    put_line(" 56. returns the diagnostics of the trace grid;");
    put_line(" 57. compute the difference in the trace sum for a factor;");
    put_line(" 58. find the index of a solution label in a slice;");
    put_line(" 59. initialize number of slices in Sampling_Operations;");
    put_line(" 60. adds a new slice to Sampling_Operations;");
    put_line(" 61. retrieves a slice from Sampling_Operations;");
    put_line(" 62. sets target slices to a previously stored set of slices;");
    put_line(" 63. completes one loop starting at one solution;");
    put_line(" 64. reads one witness set from file;");
    put_line(" 65. reads one witness set to file.");
    put_line(" 66. replaces system container with embedding of system;");
    put_line(" 67. load a polynomial from the container into a string.");
    put_line("MENU for tuning continuation parameters and output settings :");
    put_line(" 70. interactive tuning of the continuation parameters;");
    put_line(" 71. interactive setting of output level during continuation;");
    put_line(" 72. retrieve the values of the continuation parameters;");
    put_line(" 73. set the continuation parameters with given values.");
    put_line("ACCESS to the blackbox solver of PHCpack :");
    put_line(" 74. pass a Laurential as string to the system container;");
    put_line(" 75. solves the system in the Laurent system container;");
    put_line(" 76. pass a polynomial as string to the system container;");
    put_line(" 77. solves the system in the system container;");
    put_line(" 78. computes the mixed volume of the system in the container;");
    put_line(" 79. computes the mixed volume and the stable mixed volume;");
    put_line("MENU for the operations in the cells container :");
    put_line(" 80. read a mixed-cell configuration, store in container;");
    put_line(" 81. write the mixed-cell configuration in the container;");
    put_line(" 82. return the number of cells in the container;");
    put_line(" 83. return the dimension of the lifted points;");
    put_line(" 84. return number of different supports and occurrences;");
    put_line(" 85. return number of different supports and cardinalities;");
    put_line(" 86. return point of i-th support at position j;");
    put_line(" 87. return inner normal for the i-th cell;");
    put_line(" 88. return the number of points in each support of cell i;");
    put_line(" 89. return k-th point from j-th list of cell i;");
    put_line(" 90. return the mixed volume of the i-th cell;");
    put_line(" 91. sets the number of different supports and occurrences;");
    put_line(" 92. appends a point to the i-th support;");
    put_line(" 93. appends a mixed cell to the cells container;");
    put_line(" 94. clears the cells container;");
    put_line(" 95. retrieves a mixed cell from the container;");
    put_line(" 96. creation of a random coefficient system;");
    put_line(" 97. read a random coefficient system, store in container;");
    put_line(" 98. writes the random coefficient system;");
    put_line(" 99. copy random coefficient system to systems container;");
    put_line("100. copy system in systems container to cells container;");
    put_line("101. create a polyhedral homotopy to solve a random system;");
    put_line("102. solve start system corresponding to a mixed cell;");
    put_line("103. track a path starting at a solution for a mixed cell;");
    put_line("104. copy a target solution to the solutions container;");
    put_line("105. permute a given target system.");
    put_line("109. random polynomial system in the systems container.");
    put_line("MENU for linear-product root counts and start systems :");
    put_line("110. create a supporting set structure;");
    put_line("111. write the created set structure;");
    put_line("112. compute the Bezout bound;");
    put_line("MENU for operations on Laurent poly systems :");
    put_line("120. read system into the container;");
    put_line("121. write system into the container;");
    put_line("122. return dimension of the system in the container;");
    put_line("123. initialize the container with the dimension;");
    put_line("124. return number of terms in a Laurential;");
    put_line("125. return a term of a Laurential;");
    put_line("126. add a term to a Laurential;");
    put_line("127. clear the Laurent system container.");
    put_line("MENU for incremental read/write of solutions :");
    put_line("130. prompt user for input file for solution and open it;");
    put_line("131. prompt user for output file for solution and create it;");
    put_line("132. scan solution input file for banner SOLUTIONS;");
    put_line("133. read from solution input file dimensions of solutions;");
    put_line("134. write to solution output file dimensions of solutions;");
    put_line("135. read the next solution from solution input file;");
    put_line("136. write the next solution to solution output file;");
    put_line("137. close solution input file;");
    put_line("138. close solution output file;");
    put_line("139. write solution banner to the defined output file;");
    put_line("140. write solution dimensions to the defined output file;");
    put_line("141. write next solution to the defined output file;");
    put_line("142. compute next solution of total degree start system;");
    put_line("143. compute next solution of linear product start system;");
    put_line("144. retrieve a solution of linear product start system.");
    put_line("MENU for operations for track with incremental read/write :");
    put_line("147. create an evaluator for the container system;");
    put_line("148. create an evaluator for the Jacobian matrix;");
    put_line("149. refine a solution with Newton on the container system;");
    put_line("150. read target system without its solutions;");
    put_line("151. read start system without start solutions;");
    put_line("152. create homotopy with random gamma constant;");
    put_line("153. create homotopy with gamma given in c;");
    put_line("154. clear the homotopy;");
    put_line("155. track one solution path with a silent homotopy;");
    put_line("156. track one solution path with a reporting homotopy;");
    put_line("157. writes next solution with diagnostics to defined output;");
    put_line("158. writes a string of characters to the defined output file;");
    put_line("159. writes an integer sequence to the defined output file;");
    put_line("160. writes a sequence of doubles to the defined output file;");
    put_line("161. file name given to read target system without solutions;");
    put_line("162. file name given to read start system without solutions;");
    put_line("163. file name given to read linear-product start system;");
    put_line("MENU for diagonal homotopy implementation :");
    put_line("164. create a cascade homotopy from the stored systems;");
    put_line("165. create a diagonal homotopy from the stored systems;");
    put_line("166. reads first or second witness set from file;");
    put_line("167. resets the input file to read in a witness set;");
    put_line("168. returns dimension of the diagonal homotopy at start;");
    put_line("169. create witness set for polynomial k in container;");
    put_line("170. collapse the diagonal after intersecting witness sets;");
    put_line("171. remove the last slack variable when going down cascade.");
    put_line("MENU for double double and quad double tracking :");
    put_line("172. make double double homotopy with random gamma constant;");
    put_line("173. make double double homotopy with given gamma constant;");
    put_line("174. clear double double homotopy;");
    put_line("182. make quad double homotopy with random gamma constant;");
    put_line("183. make quad double homotopy with given gamma constant;");
    put_line("184. clear quad double homotopy;");
    put_line("MENU for writing solutions to strings :");
    put_line("200. returns size of solution string;");
    put_line("201. returns string representation of a solution;");
    put_line("202. returns size of solution introduction string;");
    put_line("203. returns solution introduction string;");
    put_line("204. returns size of solution vector string;");
    put_line("205. returns solution vector string;");
    put_line("206. returns size of solution diagnostics string;");
    put_line("207. returns solution diagnostics string;");
    put_line("208. append a solution string to the solution container;");
    put_line("209. replace a solution string in the solution container.");
    put_line("General MENU to numerical Schubert calculus :");
    put_line("210. display this menu;");
    put_line("211. initialize dimensions (m,p,q);");
    put_line("212. initialize m*p + q*(m+p) input m-planes;");
    put_line("213. initialize m*p + q*(m+p) interpolation points;");
    put_line("214. store pivots of pattern at start solution curve;");
    put_line("215. store pivots of pattern at target solution curve;");
    put_line("216. store coefficients of start solution curve;");
    put_line("217. retrieve coefficients of target solution curve;");
    put_line("218. track solution path without intermediate output;");
    put_line("219. track solution path with output diagnostics;");
    put_line("220. verify intersection conditions without output;");
    put_line("221. verify intersection conditions with extra output;");
    put_line("222. destroy the state machine;");
    put_line("223. compute the combinatorial Pieri root count;");
    put_line("224. return localization poset as string in b;");
    put_line("225. run the Pieri homotopies on random input data;");
    put_line("226. generate real planes osculating a rational curve;");
    put_line("227. put Schubert polynomial system in container.");
    put_line("228. resolve general Schubert intersection condition;");
    put_line("229. Littlewood-Richardson homotopies for Schubert problems.");
    put_line("MENU for double double & quad double plain C to PHCpack :");
    put_line("231. read double double target system;");
    put_line("232. write double double target system;");
    put_line("233. read double double start system;");
    put_line("234. write double double start system;");
    put_line("235. write double double start solutions;");
    put_line("236. solve by double double homotopy continuation;");
    put_line("237. write the double double target solutions;");
    put_line("238. clear the double double data;");
    put_line("241. read quad double target system;");
    put_line("242. write quad double target system;");
    put_line("243. read quad double start system;");
    put_line("244. write quad double start system;");
    put_line("245. write quad double start solutions;");
    put_line("246. solve by quad double homotopy continuation;");
    put_line("247. write the quad double target solutions;");
    put_line("248. clear the quad double data.");
    put_line("MENU to move double double data :");
    put_line("251. copy target system to systems container;");
    put_line("252. copy system in container to target system;");
    put_line("253. copy start system to systems container;");
    put_line("254. copy system in container to start system;");
    put_line("255. copy target solutions to solutions container;");
    put_line("256. copy solutions in container to target solutions;");
    put_line("257. copy start solutions to solutions container;");
    put_line("258. copy solutions in container to start solutions;");
    put_line("259. verify solutions container using systems container;");
    put_line("MENU to move quad double data :");
    put_line("261. copy target system to systems container;");
    put_line("262. copy system in container to target system;");
    put_line("263. copy start system to systems container;");
    put_line("264. copy system in container to start system;");
    put_line("265. copy target solutions to solutions container;");
    put_line("266. copy solutions in container to target solutions;");
    put_line("267. copy start solutions to solutions container;");
    put_line("268. copy solutions in container to start solutions;");
    put_line("269. verify solutions container using systems container;");
    put_line("MENU to move multiprecision data :");
    put_line("281. copy target system to systems container;");
    put_line("282. copy system in container to target system;");
    put_line("283. copy start system to systems container;");
    put_line("284. copy system in container to start system;");
    put_line("285. copy target solutions to solutions container;");
    put_line("286. copy solutions in container to target solutions;");
    put_line("287. copy start solutions to solutions container;");
    put_line("288. copy solutions in container to start solutions;");
    put_line("289. verify solutions container using systems container;");
    put_line("MENU for operations on pool of systems :");
    put_line("300. initializes the pool of systems;");
    put_line("301. returns size of the systems pool;");
    put_line("302. read and create k-th system in the pool;");
    put_line("303. write the k-th system in the pool;");
    put_line("304. create k-th system from the systems container;");
    put_line("305. use k-th system to refine a root.");
    put_line("MENU for operations on pool of solution lists :");
    put_line("320. initializes the pool of solution lists;");
    put_line("321. returns the size of the solution pool;");
    put_line("322. returns the length of a solution list;");
    put_line("323. returns the dimension of a solution list;");
    put_line("324. appends a solution to a list in the pool;");
    put_line("325. retrieves a solution from a list in the pool.");
    put_line("MENU for operations on double double complex poly systems :");
    put_line("330. read system into the container;");
    put_line("331. write system into the container;");
    put_line("332. return dimension of the system in the container;");
    put_line("333. initialize the container with the dimension;");
    put_line("334. return number of terms of a polynomial;");
    put_line("335. return a term in a polynomial;");
    put_line("336. add a term to a polynomial;");
    put_line("337. clear the double double poly system container;");
    put_line("338. store double double polynomial given as a string;");
    put_line("339. returns the degree of a double double polynomial.");
    put_line("MENU for operations on quad double complex poly systems :");
    put_line("380. read system into the container;");
    put_line("381. write system into the container;");
    put_line("382. return dimension of the system in the container;");
    put_line("383. initialize the container with the dimension;");
    put_line("384. return number of terms of a polynomial;");
    put_line("385. return a term in a polynomial;");
    put_line("386. add a term to a polynomial;");
    put_line("387. clear the quad double poly system container;");
    put_line("388. store quad double polynomial given as a string;");
    put_line("389. returns the degree of a quad double polynomial.");
    put_line("MENU for operations on monomial maps :");
    put_line("430. solve binomial system stored in Laurent system container;");
    put_line("431. write monomial maps to screen;");
    put_line("432. clear all monomial maps stored in container;");
    put_line("433. returns top dimension of the monomial maps;");
    put_line("434. returns number of maps of given dimension;");
    put_line("435. returns degree of map, given dimension and index;");
    put_line("436. returns coefficients of map, given dimension and index;");
    put_line("437. returns exponents of map, given dimension and index;");
    put_line("438. returns coefficients and exponents of map.");
    put_line("MENU for operations on multiprecision complex poly systems :");
    put_line("440. read system into the container;");
    put_line("441. write system into the container;");
    put_line("442. return dimension of the system in the container;");
    put_line("443. initialize the container with the dimension;");
    put_line("444. return number of terms of a polynomial;");
    put_line("447. clear the multiprecision systems container;");
    put_line("448. store multiprecision polynomial given as a string;");
    put_line("449. returns the degree of a multiprecision polynomial.");
    put_line("MENU for path tracking with generators :");
    put_line("500. initialize standard homotopy with target and start system;");
    put_line("501. initialize double double homotopy with target and start;");
    put_line("502. initialize quad double homotopy with target and start;");
    put_line("503. initialize start solution in standard double precision;");
    put_line("504. initialize start solution in double double precision;");
    put_line("505. initialize start solution in quad double precision;");
    put_line("506. compute next solution along path with standard doubles;");
    put_line("507. compute next solution along path with double doubles;");
    put_line("508. compute next solution along path with quad doubles;");
    put_line("509. deallocate and reset standard path tracker;");
    put_line("510. deallocate and reset double double path tracker;");
    put_line("511. deallocate and reset quad double path tracker.");
  end Write_Menu;

  function Job9 return integer32 is -- verify the solutions in the container

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Standard_Root_Refiners;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work : Solution_List;
   -- tmp : Solution_List;

  begin
   -- put_line("refining the roots ...");
   -- put(lp.all);
    if lp = null or Is_Null(sols) then
      return 9;             
    else
    -- put_line("Solutions before the refine :");
    -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
      Copy(sols,work);
      declare
        epsxa : constant double_float := 1.0E-12;
        epsfa : constant double_float := 1.0E-12;
        tolsi : constant double_float := 1.0E-8;
        deflate : boolean := false;
        nit : natural32 := 0;
      begin
        if PHCpack_Operations.Is_File_Defined then
          put(PHCpack_Operations.output_file,natural32(lp'last),lp.all);
          Reporting_Root_Refiner
            (PHCpack_Operations.output_file,
             lp.all,work,epsxa,epsfa,tolsi,nit,5,deflate,false);
        else
          Silent_Root_Refiner
            (lp.all,work,epsxa,epsfa,tolsi,nit,5,deflate);
        end if;
        Standard_Solutions_Container.Clear;
        Standard_Solutions_Container.Initialize(work);
       -- put("Refiner Results :");
       -- tmp := work;
       -- while not Is_Null(tmp) loop
       --   put(" ");
       --   put(Head_Of(tmp).m,1);
       --   tmp := Tail_Of(tmp);
       -- end loop;
       -- new_line;
      end;
      return 0;
     -- sols := Standard_Solutions_Container.Retrieve;
     -- put_line("Solutions after the refine :");
     -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
    end if;
  exception
    when others =>
      put_line("Exception 9 at verification of solutions.");
      return 9;
  end Job9;

  function Job179 return integer32 is -- variable precision Newton step

    use Interfaces.C;
    use Verification_of_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(5));
    dim : constant natural32 := natural32(v_a(v_a'first));
    ns : constant natural32 := natural32(v_a(v_a'first+1));
    wanted : constant natural32 := natural32(v_a(v_a'first+2));
    maxitr : constant natural32 := natural32(v_a(v_a'first+3));
    maxprc : constant natural32 := natural32(v_a(v_a'first+4));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(ns-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(ns));
    s : constant String(1..integer(ns)) := C_Integer_Array_to_String(ns,v_b);
    ls : Array_of_Strings(1..integer(dim)) := Split(natural(dim),s,';');
    sols : constant Multprec_Complex_Solutions.Solution_List
         := Multprec_Solutions_Container.Retrieve;
    work : Multprec_Complex_Solutions.Solution_List;

  begin
   -- put("The dimension : "); put(integer32(dim),1); new_line;
   -- put("The number of characters in the system : ");
   -- put(integer32(ns),1); new_line;
   -- put("The number of wanted number of decimal places : ");
   -- put(integer32(wanted),1); new_line;
   -- put("The maximum number of Newton steps : ");
   -- put(integer32(maxitr),1); new_line;
   -- put("The maximum number of decimal places to estimate the loss : ");
   -- put(integer32(maxprc),1); new_line;
   -- put_line("The polynomials : "); put_line(s);
   -- for i in ls'range loop
   --   put("Polynomial "); put(integer32(i),1); put_line(" :");
   --   put_line(ls(i).all);
   -- end loop;
    Multprec_Complex_Solutions.Copy(sols,work);
    Verify_Solutions_of_Laurent_Polynomials
      (ls,work,wanted,maxitr,maxprc);
    Multprec_Solutions_Container.Clear;
    Multprec_Solutions_Container.Initialize(work);
    Clear(ls);
    return 0;
  exception
    when others =>
      put_line("Exception occurred in variable precision Newton step.");
      return 179;
  end Job179;

  function Job439 return integer32 is -- scan for the number of symbols

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbc : constant natural32 := natural32(v_a(v_a'first));
    nbc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nbc-1);
    v_b : constant C_Integer_Array(0..nbc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
    s : constant String(1..integer(nbc)) := C_Integer_Array_to_String(nbc,v_b);
    cnt : constant natural := String_Splitters.Count_Delimiters(s,';');
    ls : Array_of_Strings(1..integer(cnt))
       := String_Splitters.Split(cnt,s,';');
    dim : natural32;
    maxvar : constant natural32 := 2*natural32(cnt); -- better than 1024

  begin
    dim := Parse_Dimensions.Dim(maxvar,ls);
    Assign(integer32(dim),a);
    String_Splitters.Clear(ls);
    Parse_Dimensions.Clear;
    return 0;
  exception
    when others =>
      put_line("Exception occurred when scanning for number of symbols.");
      return 439;
  end Job439;

  function Job196 return integer32 is -- apply deflation as in main driver

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Interfaces.C;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work : Solution_List;
   -- symbolic,output : boolean;
   -- nbdgts : natural32;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    maxitr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    maxdef : constant natural32 := natural32(v_b(v_b'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(3));
    tolerr : constant double_float := double_float(v_c(v_c'first));
    tolres : constant double_float := double_float(v_c(v_c'first+1));
    tolrnk : constant double_float := double_float(v_c(v_c'first+2));

  begin
    Copy(sols,work);
   -- Deflate_Singularities(lp.all,work);
   -- Set_Default_Parameters
   --   (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    Standard_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (lp.all,work,maxitr,maxdef,tolerr,tolres,tolrnk);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(work);
    return 0;
  end Job196;

  function Job249 return integer32 is -- deflation in double double precision

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;
    use Interfaces.C;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    work : Solution_List;
   -- symbolic,output : boolean;
   -- nbdgts : natural32;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    maxitr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    maxdef : constant natural32 := natural32(v_b(v_b'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(3));
    tolerr : constant double_float := double_float(v_c(v_c'first));
    tolres : constant double_float := double_float(v_c(v_c'first+1));
    tolrnk : constant double_float := double_float(v_c(v_c'first+2));

  begin
    Copy(sols,work);
   -- Deflate_Singularities(lp.all,work);
   -- Set_Default_Parameters
   --   (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    DoblDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (lp.all,work,maxitr,maxdef,tolerr,tolres,tolrnk);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(work);
    return 0;
  end Job249;

  function Job250 return integer32 is -- deflation in quad double precision

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;
    use Interfaces.C;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    work : Solution_List;
   -- symbolic,output : boolean;
   -- nbdgts : natural32;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    maxitr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    maxdef : constant natural32 := natural32(v_b(v_b'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(3));
    tolerr : constant double_float := double_float(v_c(v_c'first));
    tolres : constant double_float := double_float(v_c(v_c'first+1));
    tolrnk : constant double_float := double_float(v_c(v_c'first+2));

  begin
    Copy(sols,work);
   -- Deflate_Singularities(lp.all,work);
   -- Set_Default_Parameters
   --   (symbolic,output,maxitr,maxdef,nbdgts,tolerr,tolres,tolrnk);
    QuadDobl_Deflation_Methods.Algorithmic_Deflation_and_Clustering
      (lp.all,work,maxitr,maxdef,tolerr,tolres,tolrnk);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(work);
    return 0;
  end Job250;

  function Job628 return integer32 is -- create standard Laurent cascade
  begin
    return 0;
  exception
    when others => return 628;
  end Job628;

  function Job629 return integer32 is -- create dobldobl Laurent cascade
  begin
    return 0;
  exception
    when others => return 629;
  end Job629;

  function Job191 return integer32 is -- define output file from string

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    s : constant String(1..integer(n)) := C_Integer_Array_to_String(n,v_b);

  begin
    PHCpack_Operations.Define_Output_File(s);
    return 0;
  exception
    when others =>
      put_line("Exception raised when defining output file from string.");
      return 191;
  end Job191;

  function Job192 return integer32 is -- close the defined output file
  begin
    PHCpack_Operations.Close_Output_File;
    return 0;
  exception
    when others =>
      put_line("Exception raised when closing defined output file.");
      return 192;
  end Job192;

  function Job491 return integer32 is -- read multiprecision target system

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    deci : constant natural32 := natural32(v_a(v_a'first));

  begin
    PHCpack_Operations_io.Read_Multprec_Target_System(deci);
    return 0;
  end Job491;

  function Job493 return integer32 is -- read multiprecision start system

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    deci : constant natural32 := natural32(v_a(v_a'first));

  begin
    PHCpack_Operations_io.Read_Multprec_Start_System(deci);
    return 0;
  end Job493;

  function Job16 return integer32 is -- call standard path trackers

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks

  begin
    return PHCpack_Operations.Solve_by_Standard_Homotopy_Continuation(nbt);
  end Job16;

  function Job236 return integer32 is -- call double double path trackers

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks

  begin
    return PHCpack_Operations.Solve_by_DoblDobl_Homotopy_Continuation(nbt);
  end Job236;

  function Job246 return integer32 is -- call quad double path trackers

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks

  begin
    return PHCpack_Operations.Solve_by_QuadDobl_Homotopy_Continuation(nbt);
  end Job246;

  function Job496 return integer32 is -- call multiprecision trackers

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    deci : constant natural32 := natural32(v_a(v_a'first));

  begin
    return PHCpack_Operations.Solve_by_Multprec_Homotopy_Continuation(deci);
  end Job496;

  function Job774 return integer32 is -- solve by standard Laurent continuation

    use PHCpack_Operations;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant natural32 := natural32(v_a(v_a'first));

  begin
    return Solve_by_Standard_Laurent_Homotopy_Continuation(nbr);
  end Job774;

  function Job775 return integer32 is -- solve by dobldobl Laurent continuation

    use PHCpack_Operations;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant natural32 := natural32(v_a(v_a'first));

  begin
    return Solve_by_DoblDobl_Laurent_Homotopy_Continuation(nbr);
  end Job775;

  function Job776 return integer32 is -- solve by quaddobl Laurent continuation

    use PHCpack_Operations;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant natural32 := natural32(v_a(v_a'first));

  begin
    return Solve_by_QuadDobl_Laurent_Homotopy_Continuation(nbr);
  end Job776;

  function Handle_Jobs return integer32 is

    use Job_Containers;
    use Job_Handlers;
    use Symbol_Table_Interface;
    use Newton_Interface;
    use Continuation_Parameters_Interface;
    use Witness_Interface;

  begin
    if vrblvl > 0
     then put_line("-> in use_c2phc4c.Handle_Jobs ...");
    end if;
    case job is
      when 0 => Write_Menu; return 0;
      when 1 => return Standard_Target_Poly_System_to_Container(vrblvl-1);
      when 2 => return Standard_Container_Poly_System_to_Target(vrblvl-1);
      when 3 => return Standard_Start_Poly_System_to_Container(vrblvl-1);
      when 4 => return Standard_Container_Poly_System_to_Start(vrblvl-1);
      when 5 => return Standard_Target_Solutions_to_Container(vrblvl-1);
      when 6 => return Standard_Container_Solutions_to_Target(vrblvl-1);
      when 7 => return Standard_Start_Solutions_to_Container(vrblvl-1);
      when 8 => return Standard_Container_Solutions_to_Start(vrblvl-1);
      when 9 => return Job9; -- verify the solutions in the container
      when 10..15 => return C_to_PHCpack(job-10,0);
      when 16 => return Job16; -- call standard path trackers
      when 17..19 => return C_to_PHCpack(job-10,0);
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
      when 179 => return Job179;
     -- double double and quad double L-R homotopies :
      when 180..181 => return use_c2lrhom(job-178,a,b,c,vrblvl-1);
     -- track operations for quad double precision :
      when 182..188 => return use_track(job-150,a,b,c,vrblvl-1);
     -- tuning continuation parameters, deflation, and Newton step
      when 189 => return Continuation_Parameters_Get_Value(a,c,vrblvl-1);
      when 190 => return Continuation_Parameters_Set_Value(a,c,vrblvl-1);
      when 191 => return Job191; -- define output file from string
      when 192 => return Job192; -- close the defined output file
      when 193 => return Continuation_Parameters_Autotune(a,b,vrblvl-1);
      when 194 => return Continuation_Parameters_Show(vrblvl-1);
      when 195 => return Newton_Multprec_Polynomial_Step(a,vrblvl-1);
      when 196 => return Job196; -- apply deflation
      when 197 => return Newton_QuadDobl_Polynomial_Step(vrblvl-1);
      when 198 => return Newton_DoblDobl_Polynomial_Step(vrblvl-1);
      when 199 => return Newton_Standard_Polynomial_Step(vrblvl-1);
      when 200..209 => return use_solcon(job-170,a,b,c,vrblvl-1);
      when 210..227 => return use_c2pieri(job-210,a,b,c,vrblvl-1);
      when 228..229 => return use_c2lrhom(job-228,a,b,c,vrblvl-1);
      when 230 => return use_track(42,a,b,c,vrblvl-1);
      when 231..235 => return C_to_PHCpack(job-220,0);
      when 236 => return Job236; -- solve by double double path tracking
      when 237..238 => return C_to_PHCpack(job-220,0);
      when 239 => return use_celcon(46,a,b,c,vrblvl-1);
      when 240 => return use_celcon(47,a,b,c,vrblvl-1);
      when 241..245 => return C_to_PHCpack(job-220,0);
      when 246 => return Job246; -- solve by quad double path tracking
      when 247..248 => return C_to_PHCpack(job-220,0);
     -- deflation in double double and quad double precision
      when 249 => return Job249; -- double double deflate
      when 250 => return Job250; -- quad double deflate
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
      when 439 => return Job439;
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
      when 491 => return Job491; -- read multiprecision target system   
      when 492 => PHCpack_Operations_io.Write_Multprec_Target_System; return 0;
      when 493 => return Job493; -- read multiprecision start system
      when 494 => PHCpack_Operations_io.Write_Multprec_Start_System; return 0;
      when 495 =>
        PHCpack_Operations_io.Write_Multprec_Start_Solutions; return 0;
      when 496 => return Job496;
      when 497 =>
        PHCpack_Operations_io.Write_Multprec_Target_Solutions; return 0;
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
     -- cascade of Laurent systems :
      when 628 => return Job628; -- create double Laurent cascade
      when 629 => return Job629; -- create double double Laurent cascade
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
      when 732 => return use_multip(0,a,b,c); -- standard double precision
      when 733 => return use_multip(0,a,b,c); -- double double precision
      when 734 => return use_multip(0,a,b,c); -- quad double precision
     -- pade continuation
      when 735 => return use_padcon(0,a,b,c,vrblvl-1); -- set default values
      when 736 => return use_padcon(1,a,b,c,vrblvl-1); -- clear parameter vals
      when 737 => return use_padcon(2,a,b,c,vrblvl-1); -- get a parameter value
      when 738 => return use_padcon(3,a,b,c,vrblvl-1); -- set a parameter value
      when 739 => return use_padcon(4,a,b,c,vrblvl-1); -- track paths
     -- integer mixed cell configurations
      when 741..758 => return use_celcon(job-690,a,b,c,vrblvl-1);
     -- reading, writing Laurent start and target systems
      when 759..773 => return c_to_phcpack(job-730,0);
     -- solve by Laurent homotopy continuation
      when 774 => return Job774; -- in standard double precision
      when 775 => return Job775; -- in double double precision
      when 776 => return Job776; -- in quad double precision
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
      when 792 =>
        PHCpack_Operations.Create_Standard_Laurent_Homotopy; return 0;
      when 793 =>
        PHCpack_Operations.Create_DoblDobl_Laurent_Homotopy; return 0;
      when 794 =>
        PHCpack_Operations.Create_QuadDobl_Laurent_Homotopy; return 0;
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
      when 898..900 => return use_solcon(job,a,b,c,vrblvl-1); -- to affine
     -- multi-homogeneous
      when 910..912 => return use_solcon(job,a,b,c,vrblvl-1);
     -- from m-hom to affine
      when 913..915 => return use_solcon(job,a,b,c,vrblvl-1);
     -- adding a symbol passed as string
      when 897 => return use_syscon(job,a,b,c,vrblvl-1);
     -- reading solutions from file with given name
      when 916..918 => return use_solcon(job,a,b,c,vrblvl-1);
      when 920 => return Job_Handlers.Standard_Condition_Report(a,b,c,vrblvl-1);
     -- getting, setting the seed and the version string
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
