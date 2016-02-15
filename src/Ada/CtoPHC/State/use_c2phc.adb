with Interfaces.C;
with text_io;                           use text_io;
with String_Splitters;                  use String_Splitters;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
-- with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Random_Numbers;
with Multprec_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_SysFun;
with Multprec_Complex_Laur_JacoMats;
with Standard_Complex_Laur_Systems;
with Standard_Laur_Poly_Convertors;
with DoblDobl_Laur_Poly_Convertors;
with QuadDobl_Laur_Poly_Convertors;
with Parse_Dimensions;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Multprec_Complex_Solutions;
with Continuation_Parameters; 
with Continuation_Parameters_io; 
with Pack_Continuation_Parameters;
-- with Standard_Continuation_Data_io;     use Standard_Continuation_Data_io;
with Drivers_for_Poly_Continuation;     use Drivers_for_Poly_Continuation;
with Standard_Root_Refiners;
with DoblDobl_Root_Refiners;
with QuadDobl_Root_Refiners;
with Multprec_Root_Refiners;
with Drivers_to_Deflate_Singularities;
with Verification_of_Solutions;
with Floating_Mixed_Subdivisions;
with Black_Mixed_Volume_Computations;
with Black_Box_Solvers;
with Witness_Sets_io;
with Square_and_Embed_Systems;
with Greeting_Banners;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
-- with Assignments_of_Solutions;          use Assignments_of_Solutions;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with Multprec_PolySys_Container;
with Multprec_LaurSys_Container;
with Laurent_Systems_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;
with Multprec_Solutions_Container;
with Cells_Container;
with PHCpack_Operations;
with PHCpack_Operations_io;
with C_to_PHCpack;
with use_syscon,use_syspool;
with use_solcon,use_solpool;
with use_scaling;
with use_c2pieri,use_c2lrhom;
with use_c2fac,use_c2mbt;
with use_roco;
with use_celcon;
with use_track,use_sweep;
with use_mapcon;
with use_nxtsol;
with unisolve,use_giftwrap;
with use_numbtrop;

function use_c2phc ( job : integer32;
                     a : C_intarrs.Pointer;
		     b : C_intarrs.Pointer;
                     c : C_dblarrs.Pointer ) return integer32 is

  function Version_String return integer32 is

  -- DESCRIPTION :
  --   Returns the PHCpack version string in b,
  --   with in a the number of characters in the string.

    s : constant string := Greeting_Banners.Version;
    n : constant integer32 := integer32(s'last);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);

  begin
    Assign(n,a);
    Assign(sv,b);
    return 0;
  exception
    when others => return 999;
  end Version_String;

  function Get_Seed return integer32 is -- job 997 to get the seed

  -- DESCRIPTION :
  --   Gets the seed used in the random number generators
  --   and returns its value in the parameter a on return.
  --   Having the seed could be useful for reproducible runs.

    seed : constant integer32 := Standard_Random_Numbers.Get_Seed;

  begin
    Assign(seed,a);
    return 0;
  exception
    when others => return 997;
  end Get_Seed;

  function Set_Seed return integer32 is

  -- DESCRIPTION :
  --   Takes the number stored in a
  --   and uses its value to initialize the seed
  --   for the random number generator.
 
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    fixed_seed : constant natural32 := natural32(v_a(v_a'first));

  begin
    Standard_Random_Numbers.Set_Seed(fixed_seed);
    return 0;
  exception
    when others => return 998;
  end Set_Seed;

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

  function Job1 return integer32 is -- copy target system to container

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 1;
    else
      Standard_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 1 at copy of target system to container.");
      return 1;
  end Job1;

  function Job251 return integer32 is -- copy target system to container

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 251;
    else
      DoblDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 251 at copy of target system to container.");
      return 251;
  end Job251;

  function Job261 return integer32 is -- copy target system to container

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 261;
    else
      QuadDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 261 at copy of target system to container.");
      return 261;
  end Job261;

  function Job281 return integer32 is -- copy target system to container

    use Multprec_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Target_System(lp);
    if lp = null then
      return 281;
    else
      Multprec_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 281 at copy of target system to container.");
      return 281;
  end Job281;

  function Job2 return integer32 is -- copy target system from container

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 2;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 2 at copy of target system from container.");
      return 2;
  end Job2;

  function Job252 return integer32 is -- copy target system from container

    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 252;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 252 at copy of target system from container.");
      return 252;
  end Job252;

  function Job262 return integer32 is -- copy target system from container

    use QuadDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 262;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 262 at copy of target system from container.");
      return 262;
  end Job262;

  function Job282 return integer32 is -- copy target system from container

    use Multprec_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 282;
    else
      PHCpack_Operations.Store_Target_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 282 at copy of target system from container.");
      return 282;
  end Job282;

  function Job3 return integer32 is -- copy start system to container

    use Standard_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 3;
    else
      Standard_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 3 at copy of start system to container.");
      return 3;
  end Job3;

  function Job253 return integer32 is -- copy start system to container

    use DoblDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 253;
    else
      DoblDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 253 at copy of start system to container.");
      return 253;
  end Job253;

  function Job263 return integer32 is -- copy start system to container

    use QuadDobl_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 263;
    else
      QuadDobl_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 263 at copy of start system to container.");
      return 263;
  end Job263;

  function Job283 return integer32 is -- copy start system to container

    use Multprec_Complex_Poly_Systems;
    lp : Link_to_Poly_Sys;

  begin
    PHCpack_Operations.Retrieve_Start_System(lp);
    if lp = null then
      return 283;
    else
      Multprec_PolySys_Container.Initialize(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 283 at copy of start system to container.");
      return 283;
  end Job283;

  function Job4 return integer32 is -- copy start system from container

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 4;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 4 at copy of start system from container.");
      return 4;
  end Job4;

  function Job254 return integer32 is -- copy start system from container

    use DoblDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 254;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 254 at copy of start system from container.");
      return 254;
  end Job254;

  function Job264 return integer32 is -- copy start system from container

    use QuadDobl_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 264;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 264 at copy of start system from container.");
      return 264;
  end Job264;

  function Job284 return integer32 is -- copy start system from container

    use Multprec_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;

  begin
    if lp = null then
      return 284;
    else
      PHCpack_Operations.Store_Start_System(lp.all);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 284 at copy of start system from container.");
      return 284;
  end Job284;

  function Job5 return integer32 is -- copy target solutions to container

    use Standard_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 5;
    else
      Standard_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line("Exception 5 at copy of target solutions to container.");
      return 5;
  end Job5;

  function Job255 return integer32 is -- copy target solutions to container

    use DoblDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 255;
    else
      DoblDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 255 at copy of target solutions system to container.");
      return 255;
  end Job255;

  function Job265 return integer32 is -- copy target solutions to container

    use QuadDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 265;
    else
      QuadDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 265 at copy of target solutions system to container.");
      return 265;
  end Job265;

  function Job285 return integer32 is -- copy target solutions to container

    use Multprec_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Target_Solutions(sols);
    if Is_Null(sols) then 
      return 285;
    else
      Multprec_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 285 at copy of target solutions system to container.");
      return 285;
  end Job285;

  function Job6 return integer32 is -- copy container to target solutions

    use Standard_Complex_Solutions;  
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 6;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 6 at copy of target solutions system from container.");
      return 6;
  end Job6;

  function Job256 return integer32 is -- copy container to target solutions

    use DoblDobl_Complex_Solutions;  
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 256;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 256 at copy of target solutions system from container.");
      return 256;
  end Job256;

  function Job266 return integer32 is -- copy container to target solutions

    use QuadDobl_Complex_Solutions;  
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 266;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 266 at copy of target solutions system from container.");
      return 266;
  end Job266;

  function Job286 return integer32 is -- copy container to target solutions

    use Multprec_Complex_Solutions;  
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 286;
    else
      PHCpack_Operations.Store_Target_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 286 at copy of target solutions system from container.");
      return 286;
  end Job286;

  function Job7 return integer32 is -- copy start solutions to container

    use Standard_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 7;
    else
      Standard_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 7 at copy of start solutions system to container.");
      return 7;
  end Job7;

  function Job257 return integer32 is -- copy start solutions to container

    use DoblDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 257;
    else
      DoblDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 257 at copy of start solutions system to container.");
      return 257;
  end Job257;

  function Job267 return integer32 is -- copy start solutions to container

    use QuadDobl_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 267;
    else
      QuadDobl_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 267 at copy of start solutions system to container.");
      return 267;
  end Job267;

  function Job287 return integer32 is -- copy start solutions to container

    use Multprec_Complex_Solutions;
    sols : Solution_List;

  begin
    PHCpack_Operations.Retrieve_Start_Solutions(sols);
    if Is_Null(sols) then
      return 287;
    else
      Multprec_Solutions_Container.Initialize(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 287 at copy of start solutions system to container.");
      return 287;
  end Job287;

  function Job8 return integer32 is -- copy container to start solutions

    use Standard_Complex_Solutions;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 8;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 8 at copy of start solutions system from container.");
      return 8;
  end Job8;

  function Job258 return integer32 is -- copy container to start solutions

    use DoblDobl_Complex_Solutions;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 258;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 258 at copy of start solutions system from container.");
      return 258;
  end Job258;

  function Job268 return integer32 is -- copy container to start solutions

    use QuadDobl_Complex_Solutions;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 268;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 268 at copy of start solutions system from container.");
      return 268;
  end Job268;

  function Job288 return integer32 is -- copy container to start solutions

    use Multprec_Complex_Solutions;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;

  begin
    if Is_Null(sols) then
      return 288;
    else
      PHCpack_Operations.Store_Start_Solutions(sols);
      return 0;
    end if;
  exception
    when others =>
      put_line
        ("Exception 288 at copy of start solutions system from container.");
      return 288;
  end Job288;

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

  function Job439 return integer32 is -- scan for the nubmer of symbols

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbc : constant natural32 := natural32(v_a(v_a'first));
    nbc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nbc-1);
    v_b : constant C_Integer_Array(0..nbc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc));
    s : constant String(1..integer(nbc)) := C_Integer_Array_to_String(nbc,v_b);
    cnt : constant natural := String_Splitters.Count_Delimiters(s,';');
    ls : Array_of_Strings(1..integer(cnt)) := String_Splitters.Split(cnt,s,';');
    dim : natural32;

  begin
    dim := Parse_Dimensions.Dim(1024,ls);
    Assign(integer32(dim),a);
    return 0;
  exception
    when others =>
      put_line("Exception occurred when scanning for number of symbols.");
      return 439;
  end Job439;

  function Job195 return integer32 is -- 1 Newton step on multprec containers

    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Poly_SysFun;
    use Multprec_Complex_Jaco_Matrices;
    use Multprec_Complex_Solutions;
    use Multprec_Root_Refiners;
    lp : constant Link_to_Poly_Sys := Multprec_PolySys_Container.Retrieve;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;
    work : Solution_List;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    deci : constant natural32 := natural32(v_a(v_a'first));
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(deci);

  begin
    if lp = null or Is_Null(sols) then
      return 195;             
    else
      Copy(sols,work);
      Set_Size(work,size);
      declare
        f : Eval_Poly_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          Multprec_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        Multprec_Solutions_Container.Clear;
        Multprec_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  end Job195;

  function Job329 return integer32 is -- one multprec Newton step on Laurent 

    use Multprec_Complex_Laur_Systems;
    use Multprec_Complex_Laur_SysFun;
    use Multprec_Complex_Laur_JacoMats;
    use Multprec_Complex_Solutions;
    use Multprec_Root_Refiners;
    lp : constant Link_to_Laur_Sys := Multprec_LaurSys_Container.Retrieve;
    sols : constant Solution_List := Multprec_Solutions_Container.Retrieve;
    work : Solution_List;
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    deci : constant natural32 := natural32(v_a(v_a'first));
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(deci);

  begin
    if lp = null or Is_Null(sols) then
      return 329;             
    else
      Copy(sols,work);
      Set_Size(work,size);
      declare
        f : Eval_Laur_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          Multprec_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        Multprec_Solutions_Container.Clear;
        Multprec_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  end Job329;

  function Job197 return integer32 is -- 1 Newton step on quaddobl containers

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Root_Refiners;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if lp = null or Is_Null(sols) then
      return 197;             
    else
      Copy(sols,work);
      declare
        f : Eval_Poly_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          QuadDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        QuadDobl_Solutions_Container.Clear;
        QuadDobl_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  end Job197;

  function Job328 return integer32 is -- one quaddobl Newton step on Laurent

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Root_Refiners;
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if lp = null or Is_Null(sols) then
      return 328;             
    else
      Copy(sols,work);
      declare
        f : Eval_Laur_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          QuadDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        QuadDobl_Solutions_Container.Clear;
        QuadDobl_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  end Job328;

  function Job198 return integer32 is -- 1 Newton step on dobldobl containers

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Root_Refiners;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if lp = null or Is_Null(sols) then
      return 198;             
    else
      Copy(sols,work);
      declare
        f : Eval_Poly_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          DoblDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        DoblDobl_Solutions_Container.Clear;
        DoblDobl_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  end Job198;

  function Job327 return integer32 is -- one dobldobl Newton step on Laurent

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Root_Refiners;
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if lp = null or Is_Null(sols) then
      return 327;             
    else
      Copy(sols,work);
      declare
        f : Eval_Laur_Sys(lp'range) := Create(lp.all);
        jm : Jaco_Mat(lp'range,lp'range) := Create(lp.all);
        jf : Eval_Jaco_Mat(lp'range,lp'range) := Create(jm);
        tmp : Solution_List := work;
        ls : Link_to_Solution;
      begin
        while not Is_Null(tmp) loop
          ls := Head_Of(tmp);
          DoblDobl_Newton_Step(f,jf,ls.v,ls.err,ls.rco,ls.res);
          tmp := Tail_Of(tmp);
        end loop;
        Clear(f); Clear(jm); Clear(jf);
        DoblDobl_Solutions_Container.Clear;
        DoblDobl_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  end Job327;

  function Job199 return integer32 is -- 1 Newton step on standard containers

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Standard_Root_Refiners;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    if lp = null or Is_Null(sols) then
      return 199;             
    else
      Copy(sols,work);
      declare
        epsxa : constant double_float := 1.0E-12;
        epsfa : constant double_float := 1.0E-12;
        tolsi : constant double_float := 1.0E-8;
        deflate : boolean := false;
        nit : natural32 := 0;
      begin
        Silent_Root_Refiner(lp.all,work,epsxa,epsfa,tolsi,nit,1,deflate);
        Standard_Solutions_Container.Clear;
        Standard_Solutions_Container.Initialize(work);
      end;
      return 0;
    end if;
  end Job199;

  function Job326 return integer32 is -- one standard Newton step on Laurent

    use Standard_Complex_Laur_Systems,Standard_Complex_Solutions;
    use Standard_Root_Refiners;
    lp : constant Link_to_Laur_Sys := Laurent_Systems_Container.Retrieve;
    sols : Solution_List := Standard_Solutions_Container.Retrieve;
    refsols : Solution_List;

  begin
    if lp = null or Is_Null(sols) then
      return 326;             
    else
      declare
        epsxa : constant double_float := 1.0E-12;
        epsfa : constant double_float := 1.0E-12;
        tolsi : constant double_float := 1.0E-8;
        nit : natural32 := 0;
      begin
        Silent_Root_Refiner(lp.all,sols,refsols,epsxa,epsfa,tolsi,nit,1);
        Standard_Solutions_Container.Clear;
        Standard_Solutions_Container.Initialize(sols);
      end;
      return 0;
    end if;
  end Job326;

  function Job196 return integer32 is -- apply deflation as in main driver

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Drivers_to_Deflate_Singularities;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    Copy(sols,work);
    Deflate_Singularities(lp.all,work);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(work);
    return 0;
  end Job196;

  function Job249 return integer32 is -- deflation in double double precision

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;
    use Drivers_to_Deflate_Singularities;

    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    Copy(sols,work);
    Deflate_Singularities(lp.all,work);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(work);
    return 0;
  end Job249;

  function Job250 return integer32 is -- deflation in quad double precision

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;
    use Drivers_to_Deflate_Singularities;

    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;
    work : Solution_List;

  begin
    Copy(sols,work);
    Deflate_Singularities(lp.all,work);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(work);
    return 0;
  end Job250;

  function Job66 return integer32 is -- embed standard double system

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    ep : Link_to_Poly_Sys;

  begin
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(ep.all);
    Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => return 66;
  end Job66;

  function Job129 return integer32 is -- embed double double system

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    ep : Link_to_Poly_Sys;

  begin
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(ep.all);
    Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => return 129;
  end Job129;

  function Job260 return integer32 is -- embed quad double system

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    ep : Link_to_Poly_Sys;

  begin
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(ep.all);
    Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => return 260;
  end Job260;

  function Job70 return integer32 is -- interactive tuning of parameters
  begin
    if PHCpack_Operations.Is_File_Defined then
      Driver_for_Continuation_Parameters
        (PHCpack_Operations.output_file);
    else
      Driver_for_Continuation_Parameters;
    end if;
    return 0;
  end Job70;

  function Job71 return integer32 is -- interactive setting of output level

    oc : natural32;

  begin
    if PHCpack_Operations.Is_File_Defined then
      Driver_for_Process_io
        (PHCpack_Operations.output_file,oc);
    else
      Driver_for_Process_io(standard_output,oc);
    end if;
    return 0;
  end Job71;

  function Job72 return integer32 is -- retrieve continuation parameters

    v : constant Standard_Floating_Vectors.Vector(1..34)
      := Pack_Continuation_Parameters.Get;

  begin
    Assign(v,c);
    return 0;
  end Job72;

  function Job73 return integer32 is -- set continuation parameters

    v : Standard_Floating_Vectors.Vector(1..34);

  begin
    Assign(34,c,v);
    Pack_Continuation_Parameters.Set(v);
    return 0;
  end Job73;

  function Job193 return integer32 is -- autotune continuation parameters

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    level : constant natural32  := natural32(v_a(v_a'first));
    nbdgt : constant natural32  := natural32(v_b(v_b'first));

  begin
    Continuation_Parameters.Tune(level,nbdgt);
    return 0;
  end Job193;

  function Job194 return integer32 is -- show continuation parameters
  begin
    Continuation_Parameters_io.put;
    return 0;
  end Job194;

  function Job75 return integer32 is -- solve Laurent system in container

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
    use Standard_Complex_Laur_Systems;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Laur_Sys := Laurent_Systems_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc : natural32;
    sols : Solution_List;

  begin
   -- put("nv = "); put(integer32(nv),1);
   -- put("  nq = "); put(integer32(nq),1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 75;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 75;
    end if;
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
      Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        p : constant Poly_Sys := Positive_Laurent_Polynomial_System(lp.all);
      begin
        Black_Box_Solvers.Solve(ntasks,p,silent,rc,sols);
      end;
    end if;
    Assign(integer32(rc),a);
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  end Job75;

  function Job77 return integer32 is -- calls the blackbox solver

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;
   -- n : constant natural := Standard_PolySys_Container.Dimension;
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    nt : constant natural32  := natural32(v_b(v_b'first));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc : natural32;
    sols : Solution_List;

  begin
   -- put("Dimension of the system in the container : "); put(n,1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 77;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 77;
    end if;
    Black_Box_Solvers.Solve(nt,lp.all,false,rc,sols); -- not silent by default
    Assign(integer32(rc),a);
    Standard_Solutions_Container.Initialize(sols);
    return 0;
  end Job77;

  function Job700 return integer32 is -- dobldobl polynomial blackbox solver

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    nt : constant natural32  := natural32(v_b(v_b'first));
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc : natural32;
    sols : Solution_List;

  begin
   -- put("Dimension of the system in the container : "); put(n,1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 700;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 700;
    end if;
    Black_Box_Solvers.Solve(nt,lp.all,false,rc,sols); -- not silent by default
    Assign(integer32(rc),a);
    DoblDobl_Solutions_Container.Initialize(sols);
    return 0;
  end Job700;

  function Job701 return integer32 is -- dobldobl Laurent blackbox solver

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Laur_Systems;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc : natural32;
    sols : Solution_List;

  begin
   -- put("nv = "); put(integer32(nv),1);
   -- put("  nq = "); put(integer32(nq),1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 701;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 701;
    end if;
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
      Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
    else
      declare
        use DoblDobl_Laur_Poly_Convertors;
        p : constant Poly_Sys := Positive_Laurent_Polynomial_System(lp.all);
      begin
        Black_Box_Solvers.Solve(ntasks,p,silent,rc,sols);
      end;
    end if;
    Assign(integer32(rc),a);
    DoblDobl_Solutions_Container.Initialize(sols);
    return 0;
  end Job701;

  function Job702 return integer32 is -- quaddobl polynomial blackbox solver

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    nt : constant natural32  := natural32(v_b(v_b'first));
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc : natural32;
    sols : Solution_List;

  begin
   -- put("Dimension of the system in the container : "); put(n,1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 702;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 702;
    end if;
    Black_Box_Solvers.Solve(nt,lp.all,false,rc,sols); -- not silent by default
    Assign(integer32(rc),a);
    QuadDobl_Solutions_Container.Initialize(sols);
    return 0;
  end Job702;

  function Job703 return integer32 is -- quaddobl Laurent blackbox solver

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Laur_Systems;
    use Interfaces.C;

    v_b : constant C_Integer_Array(0..1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    silval : constant natural32 := natural32(v_b(v_b'first));
    silent : constant boolean := (silval = 1);
    ntasks : constant natural32 := natural32(v_b(v_b'first+1));
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    nv : constant natural32 := Size_of_Support(lp.all);
    nq : constant natural32 := natural32(lp'last);
    rc : natural32;
    sols : Solution_List;

  begin
   -- put("nv = "); put(integer32(nv),1);
   -- put("  nq = "); put(integer32(nq),1); new_line;
    if nv < nq then
      put_line("The system is overdetermined, add slack variables.");
      return 703;
    elsif nv > nq then
      put_line("The system is underdetermined, add linear equations.");
      return 703;
    end if;
    if QuadDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
      Black_Box_Solvers.Solve(ntasks,lp.all,silent,rc,sols);
    else
      declare
        use QuadDobl_Laur_Poly_Convertors;
        p : constant Poly_Sys := Positive_Laurent_Polynomial_System(lp.all);
      begin
        Black_Box_Solvers.Solve(ntasks,p,silent,rc,sols);
      end;
    end if;
    Assign(integer32(rc),a);
    QuadDobl_Solutions_Container.Initialize(sols);
    return 0;
  end Job703;

  function Job78 return integer32 is -- computes mixed volume

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Black_Mixed_Volume_Computations;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv : natural32;

  begin
    if lp /= null then
      Black_Box_Mixed_Volume_Computation
        (lp.all,mix,perm,iprm,lifsup,mixsub,mv);
    else
      declare
        lq : constant Link_to_Laur_Sys := Laurent_Systems_Container.Retrieve;
      begin
        Black_Box_Mixed_Volume_Computation
          (lq.all,mix,perm,iprm,lifsup,mixsub,mv);
      end;
    end if;
    Assign(integer32(mv),a);
    Cells_Container.Initialize(mix,lifsup,mixsub);
    return 0;
  end Job78;

  function Job79 return integer32 is -- computes mixed and stable volume

    use Standard_Complex_Poly_Systems;
    use Black_Mixed_Volume_Computations;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    mix,perm,iprm : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    mixsub,mcc1,mcc2 : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv,smv,tot,c1,c2 : natural32;
    stlb : double_float;

  begin
    Black_Box_Mixed_Volume_Computation
      (lp.all,mix,perm,iprm,stlb,lifsup,mixsub,mcc1,mcc2,mv,smv,tot,c1,c2);
    Assign(integer32(mv),a);
    Assign(integer32(smv),b);
    Cells_Container.Initialize(mix,lifsup,mixsub);
    return 0;
  end Job79;

  function Job189 return integer32 is -- get value of continuation parameter

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    res : Standard_Floating_Vectors.Vector(1..1);

  begin
    if k = 0 or k > 34 then
      return 189;
    else
      res(1) := Pack_Continuation_Parameters.Get_Value(k);
      Assign(res,c);
    end if;
    return 0;
  exception
    when others => return 189;
  end Job189;

  function Job190 return integer32 is -- set value of continuation parameter

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    k : constant natural32 := natural32(v_a(v_a'first));
    v_c : constant C_Double_Array := C_dblarrs.Value(c);
    v : constant double_float := double_float(v_c(v_c'first));

  begin
    if k = 0 or k > 34
     then return 190;
     else Pack_Continuation_Parameters.Set_Value(k,v);
    end if;
    return 0;
  end Job190;

  function Job191 return integer32 is -- define output file from string

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
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

  function Job291 return integer32 is -- remove a symbol from the table

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    n : constant natural32 := natural32(v_a(v_a'first));

  begin
   -- put("removing symbol "); put(n,1); put(" from the table... ");
    Symbol_Table.Remove(n);
    Symbol_Table.Downsize(1);
   -- put_line(" ... done");
    return 0;
  exception
    when others =>
      put_line("Exception raised when removing a symbol from the table");
      return 291;
  end Job291;

  function Job296 return integer32 is -- remove a symbol by name

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nc : constant natural32 := natural32(v_a(v_a'first));
    vb : constant C_Integer_Array(0..Interfaces.C.size_t(nc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc+1));
    sv : constant String(1..integer(nc)) := C_Integer_Array_to_String(nc,vb);   
    sb : Symbol_Table.Symbol;
    ind : natural32;

  begin
    for i in 1..integer(nc) loop
      sb(i) := sv(i);
    end loop;
    for i in integer(nc)+1..sb'last loop
      sb(i) := ' ';
    end loop;
    ind := Symbol_Table.Get(sb);
    Symbol_Table.Remove(ind);
    Symbol_Table.Downsize(1);
    return 0;
  exception
    when others =>
      put_line("Exception raised when removing a symbol name");
      return 296;
  end Job296;

  function Job292 return integer32 is -- sort embed symbols

    use Standard_Complex_Poly_Systems;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    dim : natural32;

  begin
    if lp /= null then
      dim := Witness_Sets_io.Count_Embed_Symbols(natural32(lp'last),"zz");
      if dim > 0 then
        Witness_Sets_io.Swap_Symbols_to_End(natural32(lp'last),dim,"zz",lp.all);
        if dim > 1 then
          Witness_Sets_io.Sort_Embed_Symbols
            (natural32(lp'last),natural32(lp'last)-dim,dim,lp.all);
        end if;
      end if;
      Assign(integer32(dim),a);
    else
      Assign(-1,a);
    end if;
    return 0;
  exception
    when others =>
      put_line("Exception raised when sorting the embed symbols");
      return 292;
  end Job292;

  function Job293 return integer32 is -- number of symbols
  begin
    Assign(integer32(Symbol_Table.Number),a); 
    return 0;
  exception
    when others => return 293;
  end Job293;

  function Job294 return integer32 is -- write symbols to screen
  begin
    for i in 1..Symbol_Table.Number loop
      put(" "); Symbol_Table_io.put(Symbol_Table.Get(i));
    end loop;
    return 0;
  exception
    when others => return 294;
  end Job294;

  function Characters_of_Symbol ( i : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns the non-blank characters of the i-th symbol.

    sb : constant Symbol_Table.Symbol := Symbol_Table.Get(i);
    res : string(sb'range);
    cnt : natural32 := 0;

  begin
    for i in sb'range loop
      exit when (sb(i) = ' ');
      cnt := cnt + 1;
      res(i) := sb(i);
    end loop;
    return res(1..integer(cnt));
  end Characters_of_Symbol;

  function String_of_Symbols ( i : natural32 ) return string is

  -- DESCRIPTION :
  --   Returns a string of all symbols in the symbol table,
  --   separated by one space.

  begin
    if i > Symbol_Table.Number then
      return "";
    elsif i = Symbol_Table.Number then
      return Characters_of_Symbol(i);
    else
      return Characters_of_Symbol(i) & " " & String_of_Symbols(i+1);
    end if;  
  end String_of_Symbols;

  function Job295 return integer32 is -- return string of symbols

    s : constant string := String_of_Symbols(1);
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(s);

  begin
    Assign(integer32(s'last),a);
    Assign(sv,b);
    return 0;
  exception
    when others => return 295;
  end Job295;

  function Job491 return integer32 is -- read multiprecision target system

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    deci : constant natural32 := natural32(v_a(v_a'first));

  begin
    PHCpack_Operations_io.Read_Multprec_Target_System(deci);
    return 0;
  end Job491;

  function Job493 return integer32 is -- read multiprecision start system

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    deci : constant natural32 := natural32(v_a(v_a'first));

  begin
    PHCpack_Operations_io.Read_Multprec_Start_System(deci);
    return 0;
  end Job493;

  function Job16 return integer32 is -- call standard path trackers

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks

  begin
    return PHCpack_Operations.Solve_by_Standard_Homotopy_Continuation(nbt);
  end Job16;

  function Job236 return integer32 is -- call double double path trackers

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks

  begin
    return PHCpack_Operations.Solve_by_DoblDobl_Homotopy_Continuation(nbt);
  end Job236;

  function Job246 return integer32 is -- call quad double path trackers

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks

  begin
    return PHCpack_Operations.Solve_by_QuadDobl_Homotopy_Continuation(nbt);
  end Job246;

  function Job496 return integer32 is -- call multiprecision trackers

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    deci : constant natural32 := natural32(v_a(v_a'first));

  begin
    return PHCpack_Operations.Solve_by_Multprec_Homotopy_Continuation(deci);
  end Job496;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => Write_Menu; return 0;
      when 1 => return Job1; -- copy target system to container
      when 2 => return Job2; -- copy target system from container
      when 3 => return Job3; -- copy start system to container
      when 4 => return Job4; -- copy start system from container
      when 5 => return Job5; -- copy target solutions to container
      when 6 => return Job6; -- copy container to target solutions
      when 7 => return Job7; -- copy start solutions to container
      when 8 => return Job8; -- copy container to start solutions
      when 9 => return Job9; -- verify the solutions in the container
      when 10..15 => return C_to_PHCpack(job-10,0);
      when 16 => return Job16; -- call standard path trackers
      when 17..19 => return C_to_PHCpack(job-10,0);
      when 20..29 => return use_syscon(job-20,a,b,c);
      when 30..38 => return use_solcon(job-30,a,b,c);
      when 39 => return use_c2fac(28,a,b,c); -- set state to silent
      when 40..65 => return use_c2fac(job-40,a,b,c);
      when 66 => return Job66; -- return embedded standard double system 
      when 67 => return use_syscon(67,a,b,c); -- load polynomial as string
      when 68 => return use_c2fac(job-42,a,b,c); -- return #factors
      when 69 => return use_c2fac(job-42,a,b,c); -- return irreducible factor
      when 70 => return Job70; -- interactive tuning of parameters
      when 71 => return Job71; -- interactive setting of output level
      when 72 => return Job72; -- retrieve values of continuation parameters
      when 73 => return Job73; -- set values of continuation parameters
      when 74 => return use_syscon(74,a,b,c); -- store Laurential as string
      when 75 => return Job75; -- calls blackbox solver on Laurent system
      when 76 => return use_syscon(76,a,b,c); -- store polynomial as string
      when 77 => return Job77; -- calls blackbox solver on polynomial system
      when 78 => return Job78; -- computes mixed volume
      when 79 => return Job79; -- for mixed volume and stable mixed volume
      when 80..105 => return use_celcon(job-80,a,b,c);
      when 106 => return use_syscon(68,a,b,c); -- load dobldobl poly as string
      when 107 => return use_syscon(69,a,b,c); -- load quaddobl poly as string
      when 108 => return use_syscon(70,a,b,c); -- load multprec poly as string
      when 109 => return use_syscon(71,a,b,c); -- random system in container
      when 110..118 => return use_roco(job-110,a,b,c);
      when 119 => return use_syscon(20,a,b,c); -- degree of standard polynomial
     -- operations on Laurent systems :
      when 120..127 => return use_syscon(job-20,a,b,c);
      when 128 => return use_syscon(77,a,b,c); -- load standard Laur as string
      when 129 => return Job129; -- embed double double system
      when 130..145 => return use_solcon(job-120,a,b,c);
      when 146 => return use_solcon(9,a,b,c); -- drop coordinate by name
      when 147 => return use_syscon(10,a,b,c);
      when 148 => return use_syscon(11,a,b,c);
      when 149..171 => return use_track(job-150,a,b,c);
     -- track operations for double double precision :
      when 172..178 => return use_track(job-150,a,b,c);
     -- variable precision Newton step :
      when 179 => return Job179;
     -- double double and quad double L-R homotopies :
      when 180..181 => return use_c2lrhom(job-178,a,b,c);
     -- track operations for quad double precision :
      when 182..188 => return use_track(job-150,a,b,c);
     -- tuning continuation parameters, deflation, and Newton step
      when 189 => return Job189; -- get value of a continuation parameter
      when 190 => return Job190; -- set value of a continuation parameter
      when 191 => return Job191; -- define output file from string
      when 192 => return Job192; -- close the defined output file
      when 193 => return Job193; -- autotune continuation parameters
      when 194 => return Job194; -- print continuation parameters to screen
      when 195 => return Job195; -- 1 Newton step on multiprecision containers
      when 196 => return Job196; -- apply deflation
      when 197 => return Job197; -- 1 Newton step on quaddobl containers
      when 198 => return Job198; -- 1 Newton step on dobldobl containers
      when 199 => return Job199; -- 1 Newton step on standard containers
      when 200..209 => return use_solcon(job-170,a,b,c);
      when 210..227 => return use_c2pieri(job-210,a,b,c);
      when 228..229 => return use_c2lrhom(job-228,a,b,c);
      when 230 => return use_track(42,a,b,c);
      when 231..235 => return C_to_PHCpack(job-220,0);
      when 236 => return Job236; -- solve by double double path tracking
      when 237..238 => return C_to_PHCpack(job-220,0);
      when 239 => return use_celcon(46,a,b,c);
      when 240 => return use_celcon(47,a,b,c);
      when 241..245 => return C_to_PHCpack(job-220,0);
      when 246 => return Job246; -- solve by quad double path tracking
      when 247..248 => return C_to_PHCpack(job-220,0);
     -- deflation in double double and quad double precision
      when 249 => return Job249; -- double double deflate
      when 250 => return Job250; -- quad double deflate
     -- double double versions for jobs 1 to 8
      when 251 => return Job251; -- copy target system to container
      when 252 => return Job252; -- copy target system from container
      when 253 => return Job253; -- copy start system to container
      when 254 => return Job254; -- copy start system from container
      when 255 => return Job255; -- copy target solutions to container
      when 256 => return Job256; -- copy container to target solutions
      when 257 => return Job257; -- copy start solutions to container
      when 258 => return Job258; -- copy container to start solutions
     -- double double witness set for a hypersurface
      when 259 => return use_track(49,a,b,c);
      when 260 => return Job260; -- embed quad double system
     -- quad double versions for jobs 1 to 8
      when 261 => return Job261; -- copy target system to container
      when 262 => return Job262; -- copy target system from container
      when 263 => return Job263; -- copy start system to container
      when 264 => return Job264; -- copy start system from container
      when 265 => return Job265; -- copy target solutions to container
      when 266 => return Job266; -- copy container to target solutions
      when 267 => return Job267; -- copy start solutions to container
      when 268 => return Job268; -- copy container to start solutions
     -- quad double witness set for a hypersurface
      when 269 => return use_track(50,a,b,c);
     -- interface to diagonal homotopies ...
      when 270 => return use_track(40,a,b,c); -- standard witset of hypersurface
      when 271 => return use_track(41,a,b,c); -- start diagonal cascade sols
     -- univariate polynomial solvers
      when 272 => return unisolve(1,a,b,c); -- standard double precision
      when 273 => return unisolve(2,a,b,c); -- double double precision
      when 274 => return unisolve(3,a,b,c); -- quad double precision
      when 275 => return unisolve(4,a,b,c); -- arbitrary multiprecision
     -- read next of solutions
      when 276 => return use_solcon(276,a,b,c); -- next standard solution
      when 277 => return use_solcon(277,a,b,c); -- next double double solution
      when 278 => return use_solcon(278,a,b,c); -- next quad double solution
      when 279 => return use_solcon(279,a,b,c); -- next multprec initialize
      when 280 => return use_c2fac(29,a,b,c); -- standard random complex number
     -- multiprecision versions for jobs 1 to 8
      when 281 => return Job281; -- copy target system to container
      when 282 => return Job282; -- copy target system from container
      when 283 => return Job283; -- copy start system to container
      when 284 => return Job284; -- copy start system from container
      when 285 => return Job285; -- copy target solutions to container
      when 286 => return Job286; -- copy container to target solutions
      when 287 => return Job287; -- copy start solutions to container
      when 288 => return Job288; -- copy container to start solutions
     -- diagonal homotopy in double double and quad double precision
      when 289 => return use_track(43,a,b,c); -- dobldobl diagonal homotopy
      when 290 => return use_track(44,a,b,c); -- quaddobl diagonal homotopy
     -- manipulation of symbols
      when 291 => return Job291; -- remove a symbol from the table
      when 292 => return Job292; -- sort the embed symbols and permute
      when 293 => return Job293; -- number of symbols in the table
      when 294 => return Job294; -- writes the symbols to screen
      when 295 => return Job295; -- returns string of symbols
      when 296 => return Job296; -- removes symbol by name
     -- interface to diagonal homotopies continued
      when 297 => return use_track(45,a,b,c); -- dobldobl diagonal startsols
      when 298 => return use_track(46,a,b,c); -- quaddobl diagonal startsols
      when 299 => return use_track(47,a,b,c); -- dobldobl collapse diagonal
      when 300..305 => return use_syspool(job-300,a,b,c);
      when 306..311 => return use_syscon(job-294,a,b,c);
      when 312 => return use_track(48,a,b,c); -- quaddobl collapse diagonal
      when 320..325 => return use_solpool(job-320,a,b,c);
     -- one Newton step on Laurent system :
      when 326 => return Job326; -- standard Newton step on Laurent
      when 327 => return Job327; -- dobldobl Newton step on Laurent
      when 328 => return Job328; -- quaddobl Newton step on Laurent
      when 329 => return Job329; -- multprec Newton step on Laurent
     -- operations on double double system container
      when 330..339 => return use_syscon(job-130,a,b,c);
     -- operations on double double solution container
      when 340..349 => return use_solcon(job-300,a,b,c);
      when 370..371 => return use_solcon(job-300,a,b,c);
      when 378 => return use_solcon(job-300,a,b,c);
     -- operations on quad double system container
      when 380..389 => return use_syscon(job-170,a,b,c);
     -- operations on quad double solution container
      when 390..399 => return use_solcon(job-310,a,b,c);
      when 420..421 => return use_solcon(job-310,a,b,c);
      when 428 => return use_solcon(job-310,a,b,c);
     -- operations on monomial maps as solutions to binomial systems
      when 430..438 => return use_mapcon(job-430,a,b,c);
     -- scan for the number of variables
      when 439 => return Job439;
     -- operations on multiprecision system container
      when 440..444 => return use_syscon(job-220,a,b,c);
      when 447..449 => return use_syscon(job-220,a,b,c);
     -- operations on multiprecision solutions :
      when 450..453 => return use_solcon(job-330,a,b,c);
      when 457 => return use_solcon(job-330,a,b,c);
     -- moving pointer to the current solution
      when 454 => return use_solcon(300,a,b,c);
      when 455 => return use_solcon(301,a,b,c);
      when 456 => return use_solcon(302,a,b,c);
      when 458 => return use_solcon(303,a,b,c);
     -- polyhedral homotopies in double double precision :
      when 460..469 => return use_celcon(job-435,a,b,c);
     -- polyhedral homotopies in quad double precision :
      when 470..479 => return use_celcon(job-435,a,b,c);
     -- string representations of multiprecision solutions :
      when 480..481 => return use_solcon(job-330,a,b,c);
      when 488 => return use_solcon(job-330,a,b,c);
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
      when 500..520 => return use_nxtsol(job-500,a,b,c);
     -- multiprecision homotopies :
      when 522..524 => return use_track(job-470,a,b,c);
     -- get length of current solution string :
      when 525 => return use_solcon(304,a,b,c);
      when 526 => return use_solcon(305,a,b,c);
      when 527 => return use_solcon(306,a,b,c);
      when 528 => return use_solcon(307,a,b,c);
     -- multihomogeneous Bezout numbers and start systems
      when 530..532 => return use_roco(job-520,a,b,c);
     -- returns current solution string :
      when 533 => return use_solcon(308,a,b,c);
      when 534 => return use_solcon(309,a,b,c);
      when 535 => return use_solcon(310,a,b,c);
      when 536 => return use_solcon(311,a,b,c);
     -- homotopy membership tests
      when 537 => return use_c2mbt(0,a,b,c); -- standard membership test
      when 538 => return use_c2mbt(1,a,b,c); -- dobldobl membership test
      when 539 => return use_c2mbt(2,a,b,c); -- quaddobl membership test
     -- operations to read systems into the containers
      when 540..543 => return use_syscon(job,a,b,c);
     -- operations to read systems and solutions into the containers
      when 544..547 => return use_solcon(job,a,b,c);
     -- operations on Laurent container for double doubles :
      when 550..558 => return use_syscon(job-440,a,b,c);
      when 559 => return use_syscon(72,a,b,c);
     -- operations on Laurent container for quad doubles :
      when 560..568 => return use_syscon(job-440,a,b,c);
      when 569 => return use_syscon(73,a,b,c);
     -- operations on Laurent container for multiprecision :
      when 570..574 => return use_syscon(job-440,a,b,c);
      when 577..579 => return use_syscon(job-440,a,b,c);
     -- convex hull via giftwrapping :
      when 580..589 => return use_giftwrap(job-579,a,b,c);
     -- scaling systems and solutions :
      when 590..596 => return use_scaling(job-589,a,b,c);
     -- size limits of string representations of polynomials
      when 600..607 => return use_syscon(job-520,a,b,c);
     -- run the sweep homotopy :
      when 610..621 => return use_sweep(job-610,a,b,c);
     -- make standard monodromy breakup verbose
      when 630 => return use_c2fac(30,a,b,c);
     -- monodromy breakup in double double precision :
      when 631..649 => return use_c2fac(job-600,a,b,c);
      when 652..660 => return use_c2fac(job-600,a,b,c);
     -- monodromy breakup in quad double precision :
      when 661..679 => return use_c2fac(job-600,a,b,c);
      when 682..690 => return use_c2fac(job-600,a,b,c);
     -- blackbox solvers in double double and quad double precision
      when 700 => return Job700; -- dobldobl poly system blackbox solver
      when 701 => return Job701; -- dobldobl Laurent poly blackbox solver
      when 702 => return Job702; -- dobldobl poly system blackbox solver
      when 703 => return Job703; -- dobldobl Laurent poly blackbox solver
     -- container for numerically computed tropisms
      when 711..731 => return use_numbtrop(job-710,a,b,c);
     -- getting, setting the seed and the version string
      when 997 => return Get_Seed;
      when 998 => return Set_Seed;
      when 999 => return Version_String;
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
                 return job;
end use_c2phc;
