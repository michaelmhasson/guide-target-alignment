start = -5
extend = -1


# determine if there is a match or mismatch at given position
def sim(s1, s2, i, j):
    if s1[i-1] == s2[j-1]:
        return 1
    else:
        return -1


# initialization fxns for the three matrices
def inx(i, j):
    if j > 0 and i == 0:
        return [-float('inf'), 0]
    elif j == 0:
        return [start + extend * i, 2]
    else:
        return [0, 0]


def iny(i, j):
    if i > 0 and j == 0:
        return [-float('inf'), 0]
    elif i == 0:
        return [start + extend * j, 3]
    else:
        return [0, 0]


def inm(i, j):
    if i == 0 and j == 0:
        return [0, 0]
    elif i == 0 or j == 0:
        return [-float('inf'), 0]
    else:
        return [0, 0]


# using the fxns to make alignment matrix
def alignmat(s1, s2):
    l1 = len(s1) + 1
    l2 = len(s2) + 1
    X = [[inx(i, j) for j in range(0, l2)] for i in range(0, l1)]
    Y = [[iny(i, j) for j in range(0, l2)] for i in range(0, l1)]
    M = [[inm(i, j) for j in range(0, l2)] for i in range(0, l1)]
    for j in range(1, l2):
        for i in range(1, l1):
            if (start + M[i-1][j][0]) > X[i-1][j][0]:
                X[i][j] = [extend + start + M[i-1][j][0],1]
            else:
                X[i][j] = [extend + X[i-1][j][0],2]
            if (start + M[i][j-1][0]) > Y[i][j-1][0]:
                Y[i][j] = [extend + start + M[i][j-1][0],1]
            else:
                Y[i][j] = [extend + Y[i][j-1][0],3]
            m = max(M[i-1][j-1][0], X[i-1][j-1][0], Y[i-1][j-1][0])
            if m == M[i-1][j-1][0]:
                M[i][j] = [sim(s1, s2, i, j) + M[i-1][j-1][0],1]
            elif m == X[i-1][j-1][0]:
                M[i][j] = [sim(s1, s2, i, j) + X[i-1][j-1][0],2]
            else:
                M[i][j] = [sim(s1, s2, i, j) + Y[i-1][j-1][0],3]
    s1a = ''
    s2a = ''
    adata = ''
    insstr = ''
    ins = 0
    istart = -10
    delete = 0
    dstart = -10
    mismatch = 0
    i = len(s1)
    j = len(s2)
    inv = [M[i][j][0], X[i][j][0], Y[i][j][0]]
    inp = [M[i][j][1], X[i][j][1], Y[i][j][1]]
    p = inp[inv.index(max(inv))]
    score = max(inv)

    def newptr(p, i, j):
        inp = [M[i][j][1], X[i][j][1], Y[i][j][1]]
        return inp[p-1]
    while i > 0 or j > 0:
        if p == 1:
            s1a = s1[i-1] + s1a
            s2a = s2[j-1] + s2a
            if s1[i-1] == s2[j-1]:
                adata = 'm' + adata
            else:
                adata = 'x' + adata
                mismatch += 1
            p = newptr(p, i, j)
            i -= 1
            j -= 1
        elif p == 2:
            s1a = s1[i-1] + s1a
            s2a = '-' + s2a
            adata = 'd' + adata
            delete += 1
            if dstart == -10:
                dstart = i-1
            p = newptr(p, i, j)
            i -= 1
        else:
            s1a = '-' + s1a
            s2a = s2[j-1] + s2a
            insstr = s2[j-1] + insstr
            adata = 'i' + adata
            ins += 1
            if istart == -10:
                istart = j-1
            p = newptr(p, i, j)
            j -= 1
    return [s1a, s2a, adata, str(ins), str(delete), str(mismatch), str(score), str(dstart), str(istart), insstr]
