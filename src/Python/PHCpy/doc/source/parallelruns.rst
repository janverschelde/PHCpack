Parallel Runs
=============

Almost all computers are parallel and have multiple cores available.
In a polynomial homotopy, all solution paths can be computed
independently from each other.

::

    from phcpy.dimension import get_core_count
    nbcores = get_core_count()
    nbcores

In the experiment we use the cyclic 7-roots problem.

::

    from phcpy.families import cyclic
    c7 = cyclic(7)
    for pol in c7:
        print(pol)

which shows the polynomials

::

    x0 + x1 + x2 + x3 + x4 + x5 + x6;
    x0*x1 + x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6 + x6*x0;
    x0*x1*x2 + x1*x2*x3 + x2*x3*x4 + x3*x4*x5 + x4*x5*x6 + x5*x6*x0 + x6*x0*x1;
    x0*x1*x2*x3 + x1*x2*x3*x4 + x2*x3*x4*x5 + x3*x4*x5*x6 + x4*x5*x6*x0 + x5*x6*x0*x1 + x6*x0*x1*x2;
    x0*x1*x2*x3*x4 + x1*x2*x3*x4*x5 + x2*x3*x4*x5*x6 + x3*x4*x5*x6*x0 + x4*x5*x6*x0*x1 + x5*x6*x0*x1*x2 + x6*x0*x1*x2*x3;
    x0*x1*x2*x3*x4*x5 + x1*x2*x3*x4*x5*x6 + x2*x3*x4*x5*x6*x0 + x3*x4*x5*x6*x0*x1 + x4*x5*x6*x0*x1*x2 + x5*x6*x0*x1*x2*x3 + x6*x0*x1*x2*x3*x4;
    x0*x1*x2*x3*x4*x5*x6 - 1;

The mixed volume provides a generically sharp upper bound
on the number of isolated solutions.

::

    from phcpy.volumes import mixed_volume
    mixed_volume(c7)

For this problem, there are as many solutions as the mixed volume.
So, there are ``924`` paths to track.

::

    from phcpy.solver import solve

speeding up the solver on multiple cores
----------------------------------------

To measure the speedup, the elapsed time between the start 
and the end of the run has to be computed.  
The most honest time measurement is the *wall clock time* 
which as suggested uses the time on the wall clock.  
The timers provided by Python do not measure the CPU time 
of compiled code that is executed by the solver.

::

    from datetime import datetime

Then the code cell below does a run on one core
and prints the elapsed wall clock time.

::

    timestart = datetime.now()
    s = solve(c7)
    timestop = datetime.now()
    elapsed_onecore = timestop - timestart
    print('elapsed wall clock time on 1 core :', elapsed_onecore)

We check on ``len(s)`` 
whether we have as many solutions as the mixed volume.

Now we solve again, using all available cores,
computed earlier in ``nbcores``.

::

    timestart = datetime.now()
    s = solve(c7, tasks=nbcores)
    timestop = datetime.now()
    elapsed_manycores = timestop - timestart
    print('elapsed wall clock time on', nbcores, 'cores:', elapsed_manycores)

We observe the reduction in the elapsed wall clock time
and compute the speedup as follows.

::

    speedup = elapsed_onecore/elapsed_manycores
    speedup

quality up
----------

Can multithreading compensate for the overhead of double double arithmetic?
If we can afford the time for a sequential run, by how much can we increase
the precision in a multithreaded run in the same time or less?

The code cell below computes all solutions in double double precision,
on all avaible cores.

::

    timestart = datetime.now()
    s = solve(c7, tasks=nbcores, precision='dd')
    timestop = datetime.now()
    elapsed = timestop - timestart
    print('elasped wall clock time on', nbcores, 'cores :', elapsed)

Again, we check whether we have as many solutions as the mixed volume.

If ``elapsed < elapsed_onecore`` evaluates to ``True``,
then we have achieved quality up in the multicore run.
With the multicore run, we compensate for the cost overhead 
of double double arithmetic, if the elapsed wall clock time 
on many cores in double double precision is less than the run 
on one core in double precision."
