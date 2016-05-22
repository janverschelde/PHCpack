"""
A Schubert problem on k-planes in n-space is defined by a sequence of m
brackets that specify at what dimensions the k-planes intersect m flags.
This script prompts the user for n and k and then generates 3 random
intersection conditions on k-planes in n-space.
Another possibility is to enumerate all possible triple conditions.
"""

from random import randint
from phcpy.schubert import resolve_schubert_conditions as resolve

def partition2bracket(n, k, p):
    """
    Given in p a list of k numbers representing a partition
    for conditions of a k-plane in n-space,
    returns the corresponding bracket notation for the conditions.
    """
    result = []
    for i in range(k):
        result.append(n - k + i + 1 - p[i])
    return result

def random_partition(n, k):
    """
    Returns a decreasing list of k numbers, entries of at most n-k,
    generated uniformly with the first number at least 1.
    """
    result = [randint(1, n-k-1)]
    for i in range(k-1):
        result.append(randint(0, result[i]))
    return result

def generate_all_partitions2(n, k, L, fun=None):
    """
    Prints the partition if len(L) == k (or calls fun if not None)
    and then continues the generation.
    The argument fun is a callback function, called with n, k, and L.
    """
    if(len(L) == k):
        if(fun == None):
            print L
        else:
            fun(n, k, L)
    else:
        lastval = L[len(L)-1]+1
        if(len(L) == k-1):
            if(lastval == n-k):
                lastval = lastval-1
        for j in range(0, lastval):
            generate_all_partitions2(n, k, L+[j], fun)

def generate_all_partitions(n, k, fun=None):
    """
    Prints all partitions for a k-plane in n-space,
    following the specifications of the above generate_partition(n, k).
    The last argument is a callback function, called with n, k,
    and the list that represents the generated partition.
    """
    for i in range(1, n-k):
        # print 'i =', i
        generate_all_partitions2(n, k, [i], fun)

def random_complement(n, k, p):
    """
    Returns an increasing list q of k numbers, with the restriction
    that p[i] + q[i] < n - k.  The last element is at least 1.
    """
    result = [randint(0, n-k-p[0]-1)]
    for i in range(k-1):
        result.append(randint(result[i],n-k-p[i+1]-1))
    if(result[k-1] == 0):
        result[k-1] = 1
    return result

def generate_all_complements2(n, k, p, L, fun=None):
    """
    Prints the complement if len(L) == k (or calls fun in not None)
    and then continues the generation.
    The argument fun is a callback function, called with n, k, p, and L.
    """
    if(len(L) == k):
        if(fun == None):
            print 'pair :', p, L
        else:
            fun(n, k, p, L)
    else:
        idx = len(L)-1
        firstval = L[idx]
        if(len(L) == k-1):
            if(firstval == 0):
                firstval = 1
        for j in range(firstval, n-k-p[idx+1]):
            generate_all_complements2(n, k, p, L+[j], fun)

def generate_all_complements(n, k, p, fun=None):
    """
    Prints all complements of the partition p for a k-plane in n-space,
    following the specifications of the above generate_complement(n, k, p).
    The callback function takes as arguments n, k, p, and the complement L.
    """
    for i in range(n-k-p[0]):
        # print 'i =', i
        generate_all_complements2(n, k, p, [i], fun)

def difference(n, k, p, q):
    """
    Returns the list of elements n - k - p[i] - q[i].
    """
    result = []
    for i in range(len(p)):
        result.append(n - k - p[i] - q[i])
    return result

def random_difference(n, k, p, q):
    """
    Returns a difference that has the same sum as
    the difference between p and q, but with a redistribution.
    """
    diff = difference(n, k, p, q)
    print 'the original difference :', diff
    rootcount = 0
    while(rootcount == 0):
        S = sum(diff)
        result = []
        for i in range(k-1):
            if(S > n-k):
                r = randint(0, n-k)
            else:
                r = randint(0, S)
            result.append(r)
            S = S - r
        result.append(S)
        result.sort(reverse=True)
        print 'original:', diff, 'redistributed:', result
        ind = 0
        while((ind < k-1) and (result[ind] > n-k)):
            d = result[ind] - n + k
            result[ind] = n-k
            result[ind+1] = result[ind+1] + d
            ind = ind + 1
        result.sort(reverse=True)
        q.sort(reverse=True)
        bp = partition2bracket(n, k, p)
        bq = partition2bracket(n, k, q)
        br = partition2bracket(n, k, result)
        rootcount = resolve(n, k, [bp, bq, br], False)
    return result

def generate_all_differences2(n, k, p, q, val, r, fun=None):
    """
    Recursive enumerator for all differences,
    as called by generate_all_differences.
    """
    if(len(r) == k):
        br = partition2bracket(n, k, r)
        rc = resolve(n, k, [p, q, br], False)
        if(rc > 2):
            if(fun != None):
                fun(p,q,br,rc)
            else:
                print 'brackets :', [p, q, br], 'root count :', rc
    elif(len(r) == k - 1):
        if(val <= r[k-2]):
            generate_all_differences2(n, k, p, q, 0, r+[val], fun)
    else:
        upper = min([val, n-k, r[len(r)-1]])
        lower = val/(k - len(r))
        for i in range(lower, upper+1):
            generate_all_differences2(n, k, p, q, val-i, r+[i], fun)

def generate_all_differences(n, k, p, q, fun=None):
    """
    Computes the difference between p and q, given as partitions,
    generates all possible redistributions of the conditions
    imposed by that difference.  Each time a new triplet with
    nonzero root count is generated, the brackets are printed.
    """
    diff = difference(n, k, p, q)
    S = sum(diff)
    lower = S/k # need to redistribute S into a partition
    upper = min([S, n-k])
    q.sort(reverse=True)
    bp = partition2bracket(n, k, p)
    bq = partition2bracket(n, k, q)
    for i in range(lower,upper+1):
        generate_all_differences2(n, k, bp, bq, S-i, [i], fun)

def check_sum(n, k, p, q, r):
    """
    Given three partitions in p, q, and r,
    checks whether the sum of their entries equals k*(n-k).
    """
    S = sum(p) + sum(q) + sum(r)
    c = k*(n-k)
    ps = ', %d*(%d-%d) = ' % (n, n, k)
    result = (S == c)
    print 'sum : ' + str(S) + ps + str(c) + ' ' + str(result)

def random_triplet(n, k):
    """
    Generates a triple intersection condition on k-planes in n-space,
    returned as three partitions.
    """
    p = random_partition(n, k)
    print 'the first random partition :', p
    q = random_complement(n, k, p)
    print 'the second random partition :', q
    print 'the second random partition :', q
    r = random_difference(n, k, p, q)
    print 'the difference :', r
    q.sort(reverse=True)
    r.sort(reverse=True)
    print 'the partitions, sorted :', p, q, r
    print 'checking the sum ...'
    check_sum(n, k, p, q, r)
    return p, q, r

def write_triplet(n, k, p, q):
    """
    Given a partition p and its complement q, computes the difference,
    and checks the sum.
    """
    r = difference(n, k, p, q)
    q.sort(reverse=True)
    r.sort(reverse=True)
    # print 'the triplet :', p, q, r
    bp = partition2bracket(n, k, p)
    bq = partition2bracket(n, k, q)
    br = partition2bracket(n, k, r)
    rc = resolve(n, k, [bp, bq, br], verbose=False)
    print 'the brackets :', bp, bq, br, 'root count :', rc

def write(n, k, L):
    """
    Simple callback function for generate_all_partitions.
    """
    # print 'n =', n, 'k =', k, 'L =', L
    # generate_all_complements(n, k, L, write_triplet)
    generate_all_complements(n, k, L, generate_all_differences)

def main():
    """
    Prompts the user for n, k, and m.
    """
    n = input('Give the ambient dimension : ')
    k = input('Give the dimension of the planes : ')
    m = 3 # input('Give the number of conditions : ')
    print '-> the parameters :', n, k, m
    ans = raw_input('generate all partitions ? (y/n) ')
    if(ans == 'y'):
        generate_all_partitions(n, k, write)
    else:
        p, q, r = random_triplet(n, k)
        bp = partition2bracket(n, k, p)
        bq = partition2bracket(n, k, q)
        br = partition2bracket(n, k, r)
        print 'the bracket representation :', bp, bq, br
        rc = resolve(n,k, [bp, bq, br], verbose=False)
        print 'the root count :', rc

if __name__=="__main__": main()
