def global_alignment(query, target, match, mismatch, gap):
    n, m = len(query), len(target)
    dp = [[0]*(m+1) for _ in range(n+1)]
    pointer = [[(0, 0)]*(m+1) for _ in range(n+1)]

    # initialization
    for i in range(n+1):
        dp[i][0] = gap * i
        pointer[i][0] = (i-1, 0)
    for j in range(m+1):
        dp[0][j] = gap * j
        pointer[0][j] = (0, j-1)

    # main loop
    for i in range(1, n+1):
        for j in range(1, m+1):
            if query[i-1] == target[j-1]:
                score = match
            else:
                score = mismatch

            dp[i][j] = max(
                dp[i-1][j-1] + score,  # diagonal (match or mismatch)
                dp[i-1][j] + gap,  # vertical (insertion in target)
                dp[i][j-1] + gap  # horizontal (deletion from target)
            )
            if dp[i][j] == dp[i-1][j-1] + score:
                pointer[i][j] = (i-1, j-1)
            elif dp[i][j] == dp[i-1][j] + gap:
                pointer[i][j] = (i-1, j)
            else:
                pointer[i][j] = (i, j-1)

    # backtracking
    cigar = []
    i, j = n, m
    while (i, j) != (0, 0):
        if pointer[i][j] == (i-1, j-1):
            if query[i-1] == target[j-1]:
                cigar.append('M')  # match
            else:
                cigar.append('X')  # mismatch
        elif pointer[i][j] == (i-1, j):
            cigar.append('D')  # deletion
        else:
            cigar.append('I')  # insertion
        i, j = pointer[i][j]
    cigar = ''.join(cigar[::-1])  # reverse and join

    return dp[n][m], cigar
