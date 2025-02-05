import time

# Recursive solution without memoization
def f_rec(n):
    if n == 0:
        return 1
    if n < 0:
        return 0
    return f_rec(n-1) + f_rec(n-2) + f_rec(n-5)

# Recursive solution with memoization
def f_memo(n, memo=None):
    if memo is None:
        memo = {}
    if n in memo:
        return memo[n]
    if n == 0:
        return 1
    if n < 0:
        return 0
    memo[n] = f_memo(n-1, memo) + f_memo(n-2, memo) + f_memo(n-5, memo)
    return memo[n]

# Iterative solution
def f_it(n):
    dp = [0] * (n + 1)
    dp[0] = 1
    for i in range(1, n + 1):
        dp[i] += dp[i - 1] if i - 1 >= 0 else 0
        dp[i] += dp[i - 2] if i - 2 >= 0 else 0
        dp[i] += dp[i - 5] if i - 5 >= 0 else 0
    return dp[n]

# Timing for f_rec
n_values = [10, 25]
rec_times = {}
for n in n_values:
    start_time = time.time()
    result = f_rec(n)
    end_time = time.time()
    rec_times[n] = (result, end_time - start_time)

# Timing for f_memo
n_values = [10, 25, 50, 100]
memo_times = {}
for n in n_values:
    start_time = time.time()
    result = f_memo(n)
    end_time = time.time()
    memo_times[n] = (result, end_time - start_time)

# Timing for f_it
it_times = {}
for n in n_values:
    start_time = time.time()
    result = f_it(n)
    end_time = time.time()
    it_times[n] = (result, end_time - start_time)

rec_times, memo_times, it_times

