

# w in an inverse Lyndon word if s <' w for each proper suffix s of w
# x <' y if x is a proper prefix of y or x = ras and y = rbt

# m1, ..., mk is an inverse Lyndon factorization if:
# w = m1...mk; mi is an inverse Lyndon word; m1 << m2 << ... << mk
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[i] for i in reversed(seq))

def duval(data):
    if type(data) == str:
        result = duval_(data)
    else:
        result = []
        for word in data:
            result.append(duval_(word))
    return result


def duval_(s):
    """
    input: a string s
    output: the Lyndon factorization (also CFL) of the string s
    """
    n = len(s);
    i = 0
    res = []

    while i < n:
        j = i + 1;
        k = i
        while j < n and s[k] <= s[j]:
            if s[k] < s[j]:
                k = i
            else:
                k += 1
            j += 1

        res.append(s[i: i + j - k])
        i += j - k

    return res


# Pref_bre(w) = {(p, p') | p in an inverse Lyndon word and a proper prefix of w}
def find_prefix(w):
    """
    input: a string w
    output: (x, y) where x = w0, y = '' if w in an inverse Lyndon word
        w = xy, x = pp' where (p, p') ∈ Pref_bre(w), otherwise.
        p is an inverse Lyndon word which is a proper prefix of w = pv;
        p' is the bounded right extension of p in w.
        A bounder right extension is a proper prefix of v such that:
            - p' is an inverse Lyndon word
            - pz' is an inverse Lyndon word for each proper prefix z' of p'
            - pp' is not an inverse Lyndon word
            - p << p' (p < p' and p is not a proper prefix of p')
        Pref_bre(w) = {(p, p') | p is an inverse Lyndon word which is a non
            empty proper prefix of w }
    """
    n = len(w)
    if n == 1:
        return w + '0', ''

    i = 0;
    j = 1
    while j < n - 1 and w[j] <= w[i]:
        if w[j] < w[i]:
            i = 0
        else:
            i += 1
        j += 1

    if j == n - 1:
        if w[j] <= w[i]:
            return w + '0', ''
    return w[:j + 1], w[j + 1:]


def find_bre(x, y):
    """
    input: (x, y) where w = xy is not an inverse Lyndon word;
        x = pp' = raurb, (p, p') ∈ Pref_bre(w)
    output: (p, p', y, last) = (rau, rb, y, |r|)
    """
    w = x + y
    n = len(x) - 1
    f = get_failure_function(x[:-1])  # Border(raur)
    i = n - 1
    last = n

    while i >= 0:
        if w[f[i]] < x[-1]:
            last = f[i] - 1
        i = f[i] - 1

    return (w[:n - last - 1],
            w[n - last - 1: n + 1],
            y,
            last + 1)


def get_failure_function(s):
    f = [0 for _ in s]
    i = 1;
    j = 0;
    m = len(s)
    while i < m:
        if s[j] == s[i]:
            f[i] = j + 1
            i += 1
            j += 1
        elif j > 0:
            j = f[j - 1]
        else:
            f[i] = 0
            i += 1
    return f


def icfl(data):
    if type(data) == str:
        result = icfl_(data)
    else:
        result = []
        for word in data:
            result.append(icfl_(word))
    return result


def icfl_(w):
    """
    input: a string w
    output: the inverse factorization of w obtained with the algorithm ICFL
        If w is an inverse lyndon word, ICFL(w) = w otherwise we have w=pv
        and ICFL(v) = (m1', ..., mk') and p' bounded right extension of p in w.
        ICFL(w) = (p) + ICFL(v)         if p' = rb <= m1'
                  (pm1', m2', ..., mk') if m1' <= r
    """
    x, y = find_prefix(w)
    if x == w + '0':
        return [w]
    p, bre, y, last = find_bre(x, y)
    l = icfl(bre + y)
    if len(l[0]) > last:  # |m1'| > |r|
        l.insert(0, p)
    else:
        l[0] = p + l[0]
    return l


def cfl_icfl(data, cfl_max=30, sep=False):
    if type(data) == str:
        result = cfl_icfl_(data, cfl_max, sep)
    else:
        result = []
        for word in data:
            result.append(cfl_icfl_(word, cfl_max, sep))
    return result


def cfl_icfl_(w, cfl_max=30, sep=False):
    result = []
    cfl_fact = cfl(w)
    for factor in cfl_fact:
        if len(factor) > cfl_max:
            icfl_fact = icfl(factor)
            if sep:
                icfl_fact = ['<<'] + icfl_fact + ['>>']
            result += icfl_fact
        else:
            result.append(factor)
    return result


def d_duval(data):
    if type(data) == str:
        result = d_duval_(data, duval, k=None)
    else:
        result = []
        for word in data:
            result.append(d_duval_(word, duval))

    return result


def d_cfl(data,T):
    return d_duval(data)


def d_icfl(data,T):
    if type(data) == str:
        result = d_duval_(data, icfl, k = None)
    else:
        result = []
        for word in data:
            result.append(d_duval_(word, icfl))
    return result


def d_cfl_icfl(data,k):
    if type(data) == str:
        result = d_duval_(data, cfl_icfl, k=k)
    else:
        result = []
        for word in data:
            result.append(d_duval_(word, cfl_icfl))
    return result


def d_duval_(seq, alg, k):
    factors1 = None
    if k == None:
        factors1 = [len(i) for i in alg(seq)]
    else:
        factors1 = [len(i) for i in alg(seq,k)]

    complement = reverse_complement(seq)
    factors2 = [len(i) for i in reversed(alg(complement))]

    rest = seq;
    result = []
    while factors1 and factors2:
        if factors1[0] < factors2[0]:
            n = factors1.pop(0)
            factors2[0] = factors2[0] - n
            if factors2[0] == 0: factors2.pop(0)
        else:
            n = factors2.pop(0)
            factors1[0] = factors1[0] - n
            if factors1[0] == 0: factors1.pop(0)
        f, rest = rest[:n], rest[n:]
        result.append(f)

    while factors1:
        n = factors1.pop(0)
        f, rest = rest[:n], rest[n:]
        result.append(f)
    while factors2:
        n = factors2.pop(0)
        f, rest = rest[:n], rest[n:]
        result.append(f)

    return result


cfl = duval
cfl_comb = d_cfl
icfl_comb = d_icfl
cfl_icfl_comb = d_cfl_icfl
algs_dict = {
    'cfl': cfl,
    'icfl': icfl,
    'cfl_icfl': cfl_icfl,
    'd_cfl': d_cfl,
    'd_icfl': d_icfl,
    'd_cfl_icfl': d_cfl_icfl,
}