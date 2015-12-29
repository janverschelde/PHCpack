About PHCpack
=============

This section points to phcpy, a scripting interface to PHCpack.
The development has a long history and owes a lot to many people.

Project History
---------------

The software originated in the development of new homotopy algorithms
to solve polynomial systems.  The main novelty of the first release
of the sources was the application of polyhedral homotopies in the
blackbox solver.  Polyhedral homotopies are generically optimal for
sparse polynomial systems.  Although the number of solutions may grow
exponentially in the number of equations, variables, and degrees,
for systems where the coefficients are sufficiently generic,
every solution path defined by a polyhedral homotopy will lead
to one isolated solution.

Version 2.0 of the code implemented SAGBI and Pieri homotopies
to solve problem in enumerative geometry.  A classical problem
in Schubert calculus is the problem of the two lines that meet
four general lines in 3-space.  Pieri homotopies are generically
optimal to compute all solutions to such geometric problems.
They solve the output pole placement problem in linear systems control.
With message passing, parallel versions of the Pieri homotopies
lead to good speedups on parallel distributed memory computers.

Starting with version 2.0 was the gradual introduction of new
homotopies to deal with positive dimensional solution sets.
Cascades of homotopies provide generic points on every solution set,
at every dimension.  After the application of cascade homotopies
to compute generic points on all equidimensional components,
the application of monodromy loops with the linear trace stop test
classifies the generic points on the equidimensional component
into irreducible components.  This leads to a numerical irreducible
decomposition of the solution set of a polynomial system.
Cascade of homotopies are the top down method.
A bottom up method applies diagonal homotopies to intersect
positive dimensional solution sets in an equation-by-equation solver.

To deal with singular solutions of polynomial systems,
the deflation method was added in version 2.3.
Version 2.3 was quickly followed by a bug release 2.3.01
and subsequently by many more quick releases.
The introduction of the fast mixed volume calculator MixedVol in 2.3.13
was followed by capabilities to compute stable mixed volumes in 2.3.31,
and an upgrade of the blackbox solver in version 2.3.34.

Shared memory multitasking provided the option -t,
followed by the number of tasks, to speedup the path tracking.
Our main motivation of parallelism is to offset the extra cost
of multiprecision arithmetic, in particular double double and quad
double arithmetic.
Marking a milestone after one hundred quick releases,
version 2.4 provided path tracking methods on graphics processing units.
A collection of Python scripts defines a simple web interface to the
blackbox solver and the path trackers,
enabling the solution of polynomial systems in the cloud.

phcpy: An Application Programming Interface to PHCpack
------------------------------------------------------

Because code development on PHCpack has taken a very long time,
looking at the code may be a bit too overwhelming at first.
A good starting point could be the Python interface
and in particular phcpy, with documentation at
<http://www.math.uic.edu/~jan/phcpy_doc_html/index.html>.

The main executable ``phc`` built by the code in PHCpack 
is called at the command line with options to invoke specific tools
and with file names as arguments in which the input and output data goes.
In contrast, the scripting interface replaces the files with persistent
objects and instead of selecting options from menus, the user runs scripts.

References
----------

PHCpack relies for its fast mixed volume computation on MixedVol
and on QDlib for its double double and quad double arithmetic.
Pointers to the literature are mentioned below.


1. N. Bliss, J. Sommars, J. Verschelde and X. Yu:
   **Solving polynomial systems in the cloud with polynomial
   homotopy continuation.**
   In *Computer Algebra in Scientific Computing, 17th International 
   Workshop, CASC 2015, Aachen, Germany*,
   edited by V.P. Gerdt, W. Koepf, E.W. Mayr, and E.V. Vorozhtsov.
   Volume 9301 of *Lecture Notes in Computer Science*, pages 87-100,
   Springer-Verlag, 2015.

2. T. Gao, T. Y. Li, M. Wu:
   **Algorithm 846: MixedVol: a software package for mixed-volume 
   computation.**
   *ACM Transactions on Mathematical Software*, 31(4):555-560, 2005.

3. E. Gross, S. Petrovic, and J. Verschelde: **PHCpack in Macaulay2.**
   *The Journal of Software for Algebra and Geometry: Macaulay2*,
   5:20-25, 2013.

4. Y. Guan and J. Verschelde: 
   **PHClab: A MATLAB/Octave interface to PHCpack.**
   In *IMA Volume 148: Software for Algebraic Geometry*,
   edited by M. E. Stillman, N. Takayama, and J. Verschelde,
   pages 15-32, Springer-Verlag, 2008. 

5. Y. Hida, X.S. Li, and D.H. Bailey:
   **Algorithms for quad-double precision floating point arithmetic.**
   In *15th IEEE Symposium on Computer Arithmetic (Arith-15 2001)*,
   11-17 June 2001, Vail, CO, USA, pages 155-162.
   IEEE Computer Society, 2001.
   Shortened version of Technical Report LBNL-46996.

6. A. Leykin and J. Verschelde: 
   **PHCmaple: A Maple Interface to the Numerical Homotopy Algorithms
   in PHCpack.**
   In the *Proceedings of the Tenth International Conference 
   on Applications of Computer Algebra (ACA'2004)*,
   edited by Q. N. Tran, pages 139-147, 2004.

7. A. Leykin and J. Verschelde: 
   **Interfacing with the Numerical Homotopy Algorithms in PHCpack.**
   In *proceedings of ICMS 2006, LNCS 4151*,
   edited by A. Iglesias and N. Takayama,
   pages 354-360, Springer-Verlag, 2006. 

8. M. Lu. and B. He and Q. Luo
   **Supporting extended precision on graphics processors.**
   In *Proceedings of the Sixth International Workshop on Data 
   Management on New Hardware (DaMoN 2010), 
   June 7, 2010, Indianapolis, Indiana*, edited by
   A. Ailamaki and P.A. Boncz, pages 19-26, 2010.

9. K. Piret and J. Verschelde:
   **Sweeping Algebraic Curves for Singular Solutions.**
   *Journal of Computational and Applied Mathematics*,
   234(4): 1228-1237, 2010. 

10. A. J. Sommese, J. Verschelde, and C. W. Wampler.
    **Numerical irreducible decomposition using PHCpack.**
    In *Algebra, Geometry, and Software Systems*, 
    edited by M. Joswig and N. Takayama,
    pages 109-130. Springer-Verlag, 2003.

11. J. Verschelde:
    **Algorithm 795: PHCpack: A general-purpose solver for polynomial
    systems by homotopy continuation.**
    *ACM Transactions on Mathematical Software*, 25(2):251--276, 1999.

12. J. Verschelde:
    **Polynomial homotopy continuation with PHCpack.**
    *ACM Communications in Computer Algebra*, 44(4):217-220, 2010.

13. J. Verschelde:
    **Modernizing PHCpack through phcpy.**
    In Proceedings of the 6th European Conference on Python in Science
    (EuroSciPy 2013), edited by Pierre de Buyl and Nelle Varoquaux,
    pages 71-76, 2014, available at
    <http://arxiv.org/abs/1310.0056>.

14. J. Verschelde and G. Yoffe.
    **Polynomial homotopies on multicore workstations.**
    In M.M. Maza and J.-L. Roch, editors, *Proceedings of the 4th
    International Workshop on Parallel Symbolic Computation (PASCO 2010),
    July 21-23 2010, Grenoble, France*, pages 131--140. ACM, 2010.

15. J. Verschelde and X. Yu:
    **Polynomial Homotopy Continuation on GPUs.**
    *ACM Communications in Computer Algebra*, to appear.

Users
-----

To demonstrate the relevance of the software, the first version
of the software was released with a collection of about eighty 
different polynomial systems, collected from the literature. 
This section points to a different collection of problems,
problems that have been solved by users of the software,
without intervention of its developers.

The papers listed below report the use of PHCpack in the fields of
algebraic statistics, communication networks,
geometric constraint solving, real algebraic geometry,
computation of Nash equilibria, signal processing, magnetism,
mechanical design, computational geometry, computer vision,
optimal control, image processing, pattern recognition,
global optimization, and computational physics:

1. M. Abdullahi, B.I. Mshelia, and S. Hamma:
   **Solution of polynomial system using PHCpack**.
   *Journal of Physical Sciences and Innovation*, 4:44-53, 2012.

2. Min-Ho Ahn, Dong-Oh Nam and Chung-Nim Lee:
   **Self-Calibration with Varying Focal Lengths Using 
   the Infinity Homography**. In *Proceedings of the 
   4th Asian Conference on Computer Vision* (ACCV2000),
   pages 140-145, 2000.

3. Gianni Amisano and Oreste Tristani:
   **Exact likelihood computation for nonlinear DSGE models with
   heteroskedastic innovations**.
   *Journal of Economic Dynamics and Control* 35:2167-2185, 2011.

4. D. Arzelier, C. Louembet, A. Rondepierre, and M. Kara-Zaitri:
   **A New Mixed Iterative Algorithm to Solve the Fuel-Optimal Linear 
   Impulsive Rendezvous Problem.**
   *Journal of Optimization Theory and Applications*, 2013.

5. Bassi, I.G., Abdullahi Mohammed, and Okechukwu C.E.:
   **Analysis Of Solving Polynomial Equations Using Homotopy Continuation
   Method**. *International Journal of Engineering Research &
   Technology (IJERT)* 2(8):1401-1411, 2013.

6. Daniel J. Bates and Frank Sottile:
   **Khovanskii-Rolle Continuation for Real Solutions**.
   *Foundations of Computational Mathematics* 11:563-587, 2011.

7. Jahan Bayat and Carl D. Crane III:
   **Closed-Form Equilibrium Analysis of Planar Tensegrity Mechanisms**.
   In *2006 Florida Conference on Recent Advances in Robotics*, FCRAR 2006.

8. Genevieve Belanger, Kristjan Kannike, Alexander Pukhov, and Martti Raidal:
   **Minimal semi-annihilating Z_n scalar dark matter**.
   *Journal of Cosmology and Astroparticle Physics*, June 2014 (Open Access).

9. Ivo W.M. Bleylevens, Michiel E. Hostenbach, and Ralf L.M. Peeters:
   **Polynomial Optimization and a Jacobi-Davidson type method for
   commuting matrices**,
   *Applied Mathematics and Computation* 224(1): 564-580, 2013.

10. Guy Bresler, Dustin Cartwright, David Tse:
    **Feasibility of Interference Alignment for the MIMO interference
    channel**.
    *IEEE Transactions on Information Theory* 60(9):5573-5586, 2014.

11. M.-L. G. Buot and D. St. P. Richards:
    **Counting and Locating the Solutions of Polynomial Systems of
    Maximum Likelihood Equations I**.
    *Journal of Symbolic Computation* 41(2): 234-244, 2005.

12. Max-Louis G. Buot, Serkan Hosten and Donald St. P. Richards:   
    **Counting and locating the solutions of polynomial systems of maximum 
    likelihood equations, II: The Behrens-Fisher problem**.
    *Statistica Sinica* 17(4):1343-1354, 2007.

13. Enric Celaya, Tom Creemers, Lluis Ros:
    **Exact interval propagation for the efficient solution of position
    analysis problems on planar linkages**.
    *Mechanism and Machine Theory* 54: 116-131, 2012.

14. Tom Creemers, Josep M. Porta, Lluis Ros, and Federico Thomas:
    **Fast Multiresolutive Approximations of Planar Linkage Configuration
    Spaces**. *IEEE 2006 International Conference on Robotics and Automation.*

15. R.S. Datta:
    **Using Computer Algebra To Compute Nash Equilibria**.
    In *Proceedings of the 2003 International Symposium on Symbolic and
    Algebraic Computation (ISSAC 2003)*, pages 74-79, ACM 2003.

16. R.S. Datta:
    **Finding all Nash equilibria of a finite game using
    polynomial algebra**.  *Economic Theory* 42(1):55-96, 2009.

17. B.H. Dayton:
    **Numerical Local Rings and Local Solution of Nonlinear
    Systems**.  In *Proceedings of the 2007 International Workshop on
    Symbolic-Numeric Computation (SNC'07)*, pages 79-86, ACM 2007.

18. Max Demenkov:
    **Estimating region of attraction for polynomial vector fields
    by homotopy methods**.
    *ACM Communications in Computer Algebra* 46(3):84-85, 2012.

19. Max Demenkov:
    A Matlab Tool for Regions of Attraction Estimation
    via Numerical Algebraic Geometry.</B>
    In the *2015 International Conference on Mechanics - Seventh
    Polyakhov's Reading*, February 2-6, 2015, Russia,
    Saint Petersburg State University,
    Proceedings Edited by A.A. Tikhonov. IEEE 2015.

20. Ian H. Dinwoodie, Emily Gamundi, and Ed Mosteig:
    **Multiple Solutions for Blocking Probabilities in Asymmetric Networks**.
    *Open Systems and Information Dynamics* 12(3):273-288, 2005.

21. Csaba Domokos and Zoltan Kato: 
    **Parametric Estimation of Affine Deformations of Planar Shapes**.
    *Pattern Recognition*, 2009. In press.

22. C. Durand and C.M. Hoffmann:
    **Variational Constraints in 3D**.
    In *Proceedings of the International Conference on Shape Modeling 
    and Applications*, Aizu-Wakamatsu, Japan, pages 90-98, IEEE Computer
    Society, 1999.

23. C. Durand and C.M. Hoffmann:
    **A systematic framework for solving
    geometric constraints analytically**.
    *Journal of Symbolic Computation* 30(5):493-520, 2000.

24. I.Z. Emiris, E. Tsigaridas, G. Tzoumas:
    **The predicates for the Voronoi diagram of ellipses**. 
    In *Proc. ACM Symp. Comput. Geom.* 2006. 

25. Jonathan P. Epperlein and Bassam Bamieh:
    **A Frequency Domain Method for Optimal Periodic Control**.
    2012 American Control Conference (ACC), pages 5501-5506, IEEE 2012.

26. F. Ferrari:
    **On the geometry of super Yang-Mills theories: phases and 
    irreducible polynomials**.
    *Journal of High Energy Physics* 1, paper 26, 2009.

27. Jaime Gallardo-Alvarado:
    **A simple method to solve the forward displacement analysis of
    the general six-legged parallel manipulator**.
    *Robotics and Computer-Integrated Manufacturing* 30:55-61, 2014.

28. Jaime Gallardo-Alvarado and Juan-de-Dios Posadas-Garcia:
    **Mobility analysis and kinematics of the semi-general 2(3-RPS)
    series-parallel manipulator**.
    *Robotics and Computer-Integrated Manufactoring* 29(6): 463-472, 2013.

29. Bertrand Haas:
    **A Simple Counterexample to Kouchnirenko's Conjecture**.
    *Beitraege zur Algebra und Geometrie/Contributions to Algebra
    and Geometry* 43(1):1-8, 2002.

30. Adlane Habed and Boubakeur Boufama:
    **Camera self-calibration from bivariate polynomial equations and
    the coplanarity constraint**.
    *Image and Vision Computing* 24(5):498-514, 2006.

31. Marshall Hampton and Richard Moeckel:
    **Finiteness of stationary configurations of the four-vortex problem**.
    *Transactions of the American Mathematical Society* 361(3): 1317-1332,
    2009.

32. Jonathan Hauenstein, Jose Israel Rodriguez, and Bernd Sturmfels:
    **Maximum Likelihood for Matrices with Rank Constraints**.
    *Journal of Algebraic Statistics* 5(1): 18-38, 2014.

33. Christoph Hellings, David A. Schmidt, and Wolfgang Utschick:
    **Optimized beamforming for the two stream MIMO interference channel
    at high SNR**. In 2009 Internatial ITG Workshop on Smart Antennas
    (WSA 2009), February 16-19, Berlin, Germany, pages 88-95.

34. Gabor Horvath:
    **Moment Matching-Based Distribution Fitting with Generalized
    Hyper-Erlang Distributions**.
    In *Analytical and Stochastic Modeling Techniques and Applications*,
    Lecture Notes in Computer Science, Volume 7984, pages 232-246, 2013.

35. X.G. Huang:
    **Forward Kinematics for a Parallel Platform Robot**.
    *Communications in Computer and Information Sciences* 86:529-532, 2011.

36. Xiguang Huang, Qizheng Liao, Shimin Wei, and Qiang Xu:
    **Five precision point-path synthesis of planar four-bar linkage
    using algebraic method**.
    *Frontiers of Electrical and Electronic Engineering in China*
    3(4):470-474, 2008.

37. Xiguang Huang, Qizheng Liao, Shimin Wei, Qiang Xu, and Shuguang Huang:
    **The 4SPS-2CCS generalized Stewart-Gough Platform mechanisms and its
    direct kinematics**.
    In *Proceedings of the 2007 IEEE International Conference on
    Mechatronics and Automation*, August 5-8, 2007, Harbin, China.
    Pages 2472-2477, 2007.

38. Hamadi Jamali, Tokunbo Ogunfunmi:
    **Stationary points of the finite length constant modulus optimization**.
    *Signal Processing* 82(4): 625-641, 2002.

39. Bjorn Johansson, Magnus Oskarsson, and Kalle Astrom:
    **Structure and motion estimation from complex features
    in three views**.
    In the Online ICVGIP-2002 Proceedings
    (Indian Conference on Computer Vision, Graphics and Image Processing).

40. M. Kara-Zaitri, D. Arzelier, and C. Louembet:
    **Mixed iterative algorithm for solving optimal implusive time-fixed
    rendezvous problem**.
    *American Institute of Aeronautics and Astronautics Guidance, Navigation,
    and Control Conference*, Toronto, Canada, 02-05 August 2010.

41. P.U. Lamalle, A. Messiaen, P. Dumortier, F. Durodie, M. Evrard, F. Louche:
    **Study of mutual coupling effects in the antenna array of the ICRH
    plug-in for ITER**. 
    *Fusion Engineering and Design* 74:359-365, 2005.

42. E. Lee and C. Mavroidis:
    **Solving the Geometric Design Problem of Spatial
    3R Robot Manipulators Using Polynomial Continuation**.
    *Journal of Mechanical Design, Transactions of the ASME* 124(4):652-661,
    2002.

43. E. Lee and C. Mavroidis:
    **Four Precision Points Geometric Design of Spatial 3R Manipulators**.
    In the *Proceedings of the 11th World Congress in Mechanism and Machine 
    Sciences*, August 18-21, 2003, Tianjin, China.
    China Machinery Press, edited by Tian Huang.

44. E. Lee and C. Mavroidis:
    **Geometric Design of 3R Manipulators for
    Reaching Four End-Effector Spatial Poses**.
    *International Journal for Robotics Research*, 23(3):247-254, 2004.

45. E. Lee, C. Mavroidis, and J. Morman:
    **Geometric Design of Spatial 3R Manipulators**.
    In *Proceedings of the 2002 NSF Design, Service, and
    Manufacturing Grantees and Research Conference*, San Juan, Puerto Rico,
    January 7-10, 2002.

46. Dimitri Leggas and Oleg V. Tsodikov:
    **Determination of small crystal structures from a minimum set of
    diffraction intensities by homotopy continuation**.
    *Acta Crystallographica Section A* 71(3): 319-324, 2015.

47. Dawei Leng and Weidong Sun:
    **Finding All the Solutions of PnP Problem**.
    In *IST 2009 - International Workshop on Imaging Systems and Techniques*,
    Shenzhen, China, May 11-12, 2009.  Pages 348-352, IEEE, 2009.

48. Anton Leykin:
    **Numerical Primary Decomposition**.
    In *Proceedings of ISSAC 2008*,
    edited by David Jeffrey, pages 165-164, ACM 2008.

49. Anton Leykin and Frank Sottile:
    **Computing Monodromy via Parallel Homotopy Continuation**.
    In *Proceedings of the 2007 International
    Workshop on Parallel Symbolic Computation (PASCO'07)*, 
    pages 97-98, ACM 2007. (on CDROM)

50. Anton Leykin and Frank Sottile:
    **Galois groups of Schubert problems via homotopy computation**.
    *Mathematics of Computation* 78: 1749-1765, 2009.

51. Shaobai Li, Srinandan Dasmahapatra, and Koushik Maharatna:
    **Dynamical System Approach for Edge Detection Using Coupled
    FitzHugh-Naguma Neurons**.
    *IEEE Transactions on Image Processing* 24(12), 5206-5219, 2015.

52. Ross A. Lippert:
    **Fixing multiple eigenvalues by a minimal perturbation**.
    *Linear Algebra Appl.* 432(7): 1785-1817, 2010.

53. M. Maniatis and O. Nachtmann:
    **Stability and symmetry breaking in the general three-Higgs-double
    model**.
    *Journal of High Energy Physics* 2015:58, February 2015.

54. Hyosang Moon and Nina P. Robson:
    **Design of spatial non-anthropomorphic articulated systems based on
    arm joint constraint kinematic data for human interactive robotics
    applications**. DETC2015-46530.  In the *Proceedings of the ASME 2015
    International Design Engineering Technical Conferences & Computers
    and Information in Engineering Conference*. IDETC/CIE 2015.
    August 2-5, 2015, Boston Massachusetts.

55. Marc Moreno Maza, Greg Reid, Robin Scott, and Wenyuan Wu:
    **On Approximate Triangular Decompositions I. Dimension Zero**.
    In the *SNC 2005 Proceedings*.
    International Workshop on Symbolic-Numeric Computation.
    Xi'an, China, July 19-21, 2005.
    Edited by Dongming Wang and Lihong Zhi.
    Pages 250-275, 2005.

56. Andrew J. Newell:
    **Transition to supermagnetism in chains of magnetosome crystals**.
    *Geochemistry Geophysics Geosystems* 10(11):1-19, 2009.

57. M. Oskarsson, A. Zisserman and K. Astrom:
    **Minimal Projective Reconstruction for combinations of Points
    and Lines in Three Views**.
    In the *Electronic Proceedings of BMVC2002 - The 13th British Machine
    Vision Conference 2002*, pages 63 - 72.

58. P.A. Parrilo and B. Sturmfels.
    **Minimizing polynomial functions**.
    In S. Basu and L. Gonzalez-Vega, editors,
    *Algorithmic and quantitative real algebraic geometry*,
    volume 60 of *DIMACS Series in Discrete Mathematics and 
    Theoretical Computer Science*, pages 83-99. AMS, 2003.

59. Alba Perez and J.M. McCarthy:
    **Dual Quaternion Synthesis of Constrained Robotic Systems**.
    *Journal of Mechanical Design* 126(3): 425-435, 2004.

60. Nina Patarinsky-Robson, J. Michael McCarthy, and Irem Y. Tumer:
    **The algebraic synthesis of a spatial TS chain for a prescribed
    acceleration task**.
    *Mechanism and Machine Theory* 43(10): 1268-1280, 2008.

61. Nina Patarinsky-Robson, J. Michael McCarthy, and Irem Y. Tumer:
    **Failure Recovery Planning for an Arm Mounted on an
    Exploratory Rover**.
    *IEEE Transactions on Robotics* 25(6):1448-1453, 2009.

62. Jose Israel Rodriguez:
    **Combinatorial excess intersection**.
    *Journal of Symbolic Computation* 68(2): 297-307, 2015.

63. Roger E. Sanchez-Alonso, Jose-Joel Gonzalez-Barbosa, Eduardo
    Castilo-Castaneda, and Jaime Gallardo-Alvarado:
    **Kinematic analysis of a novel 2(3-RUS) parallel manipulator**.
    *Robotica*, available on CJO2015.

64. H. Schreiber, K. Meer, and B.J. Schmitt:
    **Dimensional synthesis of planar Stephenson mechanisms for motion
    generation using circlepoint search and homotopy methods**.
    *Mechanism and Machine Theory* 37(7):717-737, 2002.

65. Ben Shirt-Ediss, Ricard V. Sole, and Kepa Ruiz-Mirazo:
    **Emergent Chemical Behavior in Variable-Volume Protocells**.
    *Life* 5: 181-121, 2015.

66. Frank Sottile:
    **Real Schubert Calculus: Polynomial systems and a conjecture
    of Shapiro and Shapiro**.
    *Experimental Mathematics* 9(2): 161-182, 2000.

67. H. Stewenius and K. Astrom:
    **Structure and Motion Problems for Multiple Rigidly Moving Cameras**.
    In *Computer Vision - ECCV 2004: 8th European Conference on
    Computer Vision, Prague, Czech Republic, May 11-14, 2004. 
    Proceedings, Part III*.  Edited by T. Pajdla and J. Matas.
    Lecture Notes in Computer Science 3023, pages 252-263, Springer, 2004.

68. H.-J. Su and J.M. McCarthy:
    **Kinematic Synthesis of RPS Serial Chains**.
    In the *Proceedings of the ASME Design Engineering Technical
    Conferences* (CDROM).
    Paper DETC03/DAC-48813.  Chicago, IL, Sept. 02-06, 2003.

69. H.-J. Su and J.M. McCarthy:
    **Synthesis of Compliant Mechanisms with Specified Equilibrium 
    Positions**. In the *Proceedings of the ASME International
    Design Engineering Technical Conferences*.
    Paper DETC 2005-85085.  Long Beach, CA, Sept. 24-28 2005.

70. H.-J. Su and J.M. McCarthy:
    **Kinematic Synthesis of RPS Serial Chains for a Given Set of 
    Task Positions**.
    *Mechanism and Machine Theory* 40(7):757-775, 2005

71. H.-J. Su and J.M. McCarthy:
    **A Polynomial Homotopy Formulation of the Inverse Static Analysis of
    Planar Compliant Mechanisms**.
    *ASME Journal of Mechanical Design* 128(4): 776-786, 2006.

72. H.-J. Su, C.W. Wampler, and J.M. McCarthy:
    **Geometric Design of Cylindric PRS Serial Chains**.
    *ASME Design Engineering Technical Conferences*,
    Chicago, IL, Sep 2-6, 2003.

73. Attila Tanács and Joakim Lindblad and Nataša Sladoje and Zoltan Ka:
    **Estimation of linear deformations of 2D and 3D fuzzy objects**.
    *Pattern Recognition* 48(4):1391-1403, 2015.

74. N. Trawny, X.S. Zhou, K.X. Zhou, S.I. Roumeliotis:
    **3D Relative Pose Estimation from Distance-Only Measurements**.
    In the *Proceedings of the 2007/IEEE/RSJ International Conference
    on intelligent Robots and Systems*. San Diego, CA, Oct 29-Nov 2, 2007,
    pages 1071-1078, IEEE, 2007.

75. T. Turocy:
    **Towards a black-box solver for finite games: Computing all equilibria
    with Gambit and PHCpack**.
    In *Software for Algebraic Geometry*, volume 148 of the IMA
    volumes in Mathematics and its Applications, edited by M.E. Stillman,
    N. Takayama, and J. Verschelde, pages 133-148, Springer-Verlag, 2008.

76. Konstantin Usevich and Ivan Markovsky:
    **Structured low-rank approximation as a rational function
    minimization**.
    In 16th IFAC Symposium on System Identification Brussels, 
    11-13 Jul 2012, pages 722-727.

77. C.W. Wampler:
    **Isotropic coordinates, circularity and Bezout numbers:
    planar kinematics from a new perspective**.
    In the *Proceedings of the 1996 ASME Design Engineering Technical
    Conference*. Irvine, CA, Aug 18-22, 1996. Available on CD-ROM.

78. Wenyuan Wu and Greg Reid:
    **Symbolic-numeric computation of implicit Riquier bases for PDE**.
    In the *Proceedings of the 2007 International Symposium on Symbolic and
    Algebraic Computation*, edited by C.W. Brown, pages 377-385, ACM 2007.

79. Jonathan Widger and Daniel Grosu:
    **Parallel Computation of Nash Equilibria in N-Player Games**.
    In the *Proceedings of the 12th IEEE International Conference
    on Computational Science and Engineering (CSE 2009)*,
    August 29-31, 2009, Vancouver, Canada, pages 209-215.

80. F. Xie, G. Reid, and S. Valluri:
    **A numerical method for the
    one dimensional action functional for FBG structures**.
    *Can J. Phys.* 76: 1-21, 2002.

81. Hong Bing Xin, Qiang Huang, and Yueqing Yu:
    **Position and Orientation Analyses of Mechanism by PHCpack Solver
    of Homotopy Continuation**.
    *Applied Mechanics and Materials* 152-254: 1779-1784, 2012.

82. K. Yang and R. Orsi:
    **Static output feedback pole placement via a trust region approach**.
    *IEEE Transactions on Automatic Control* 52(11): 2146-2150, 2007.

83. Yan Yang, Yao Zhang, Fangxing Li, and Haoyong Chen:
    **Computing All Nash Equilibria of Multiplayer Games in Electricity
    Markets by Solving Polynomial Equations**.
    *IEEE Transactions on Power Systems* 27(1): 81-91, 2012.

84. Jun Zhang and Mohan Sarovar:
    **Identification of open quantum systems from observable time traces**.
    *Physical Review A* 91, 052121, 2015.

85. Xun S. Zhou and Stergios I. Roumeliotis:
    **Determining 3-D Relative Transformations for Any Combination of
    Range and Bearing Measurements.**
    *IEEE Transactions on Robotics* 29(2):458-474, 2013.

In addition to the publications listed above, PHCpack was used as a
benchmark to measure the progress of new algorithms in the following papers:

86. T. Gao and T.Y. Li:
    **Mixed volume computation via linear programming**.
    *Taiwanese Journal of Mathematics* 4(4): 599-619, 2000.

87. T. Gao and T.Y. Li:
    **Mixed volume computation for semi-mixed systems**.
    *Discrete Comput. Geom.* 29(2):257-277, 2003.

88. L. Granvilliers:
    **On the Combination of Interval Constraint Solvers**.
    *Reliable Computing* 7(6): 467-483, 2001.

89. S. Kim and M. Kojima:
    **Numerical Stability of Path Tracing in Polyhedral Homotopy 
    Continuation Methods**.
    *Computing* 73(4): 329-348, 2004.

90. Y. Lebbah, C. Michel, M. Rueher, D. Daney, and J.P. Merlet:
    **Efficient and safe global constraints for handling numerical
    constraint systems**.
    *SIAM J. Numer. Anal.* 42(5):2076-2097, 2005.

91. T.L. Lee, T.Y. Li, and C.H. Tsai:
    **HOM4PS-2.0: a software package for solving polynomial systems
    by the polyhedral homotopy continuation method**.
    *Computing* 83(2-3): 109-133, 2008.

92. Anton Leykin:
    **Numerical Algebraic Geometry**.
    *The Journal of Software for Algebra and Geometry*
    volume 3, pages 5-10, 2011. 

93. T.Y. Li and X. Li:
    **Finding Mixed Cells in the Mixed Volume Computation**.
    *Foundations of Computational Mathematics* 1(2): 161-181, 2001.

94. T.Y. Li, X. Wang, and M. Wu:
    **Numerical Schubert Calculus by the Pieri Homotopy Algorithm**.
    *SIAM J. Numer Anal.* 40(2): 578-600, 2002.

95. J.M. Porta, L. Ros, T. Creemers, and F. Thomas:
    **Box approximations of planar linkage configuration spaces**.
    *Journal of Mechanical Design* 129(4):397-405, 2007.

96. Laurent Sorber, Marc Van Barel, and Lieven De Lathauwer:
    **Numerical solution of bivariate and polyanalytic polynomial systems**.
    *SIAM J. Numer. Anal.* 52(4):1551-1572, 2014.

97. Yang Sun, Yu-Hui Tao, Feng-Shan Bai:
    **Incomplete Groebner basis as a preconditioner for polynomial systems**.
    *Journal of Computational and Applied Mathematics* 226(1):2-9, 2009.

PHCpack was used to develop new homotopy algorithms:

98. Bo Dong, Bo Yu, and Yan Yu:
    **A symmetric and hybrid polynomial system solving method for mixed
    trigonometric polynomial systems**.
    *Mathematics of Computation* 83(288): 1847-1868, 2014.

99. Bo Yu and Bo Dong:
    **A hybrid polynomial system solving method for mixed
    trigonometric polynomial systems**.
    *SIAM J. Numer. Anal.* 46(3): 1503-1518, 2008.

100. Xuping Zhang, Jintao Zhang, and Bo Yu:
     **Eigenfunction expansion method for multiple solutions
     of semilinear elliptic equations with polynomial nonlinearity**>
     *SIAM J. Numer. Anal.* 51(5): 2680-2699, 2013.

Last, but certainly not least, there is the wonderful book of
Bernd Sturmfels which contains a section on computing Nash
equilibria with PHCpack.

101. B. Sturmfels:
     **Solving Systems of Polynomial Equations**.
     CBMS Regional Conference Series of the AMS, Number 97, 2002.

So we have to end quoting Bernd Sturmfels:
*polynomial systems are for everyone.*

Acknowledgments
---------------

This material is based upon work supported by the 
National Science Foundation under Grants No. 9804846, 0105739, 0134611,
0410036, 0713018, 1115777, and 1440534.
Any opinions, findings, and conclusions or recommendations expressed 
in this material are those of the author(s) and do not necessarily 
reflect the views of the National Science Foundation. 

Since 2001, the code in PHCpack improved thanks to the contributions
of many PhD students at the University of Illinois at Chicago.
Their names, titles of PhD dissertation, and year of PhD are listed below:

1. Yusong Wang: 
   *Computing Dynamic Output Feedback Laws with Pieri Homotopies on a 
   Parallel Computer*, 2005.

2. Ailing Zhao:
   *Newton's Method with Deflation for Isolated Singularities
   of Polynomial Systems*, 2007.

3. Yan Zhuang:
   *Parallel Implementation of Polyhedral Homotopy Methods*, 2007.

4. Kathy Piret:
   *Computing Critical Points of Polynomial Systems
   using PHCpack and Python*, 2008.

5. Yun Guan:
   *Numerical Homotopies for Algebraic Sets on a Parallel Computer*, 2010.

6. Genady Yoffe:
   *Using Parallelism to compensate for Extended Precision in Path 
   Tracking for Polynomial System Solving*, 2012.

7. Danko Adrovic:
   *Solving Polynomial Systems with Tropical Methods*, 2012.

8. Xiangcheng Yu:
   *Accelerating Polynomial Homotopy Continuation
   on Graphics Processing Units*, 2015.

Anton Leykin contributed to the application of message passing 
in a parallel implementation of monodromy to decompose an equidimensional
solution set into irreducible components.  
The Maple interface ``PHCmaple`` was written jointly with Anton Leykin.
The work of Anton Leykin also paved the way for the Macaulay2 interface,
which was further developed into ``PHCpack.m2`` in joint work with
Elizabeth Gross and Sonja Petrovic.

The software has been developed with GNAT GPL, the gnu-ada compiler.
