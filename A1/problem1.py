import time

# Recursive solution without memoization
def f_rec(n):
    if n < 0:
        return 0
    if n == 0:
        return 1
    return f_rec(n - 1) + f_rec(n - 2) + f_rec(n - 5)

# Recursive solution with memoization
def f_memo(n, memo=None):
    if memo is None:
        memo = {}
    if n < 0:
        return 0
    if n == 0:
        return 1
    if n in memo:
        return memo[n]
    memo[n] = f_memo(n - 1, memo) + f_memo(n - 2, memo) + f_memo(n - 5, memo)
    return memo[n]

# Iterative solution
def f_it(n):
    if n < 0:
        return 0
    dp = [0] * (n + 1)
    dp[0] = 1
    for i in range(1, n + 1):
        dp[i] = (dp[i - 1] if i >= 1 else 0) + \
                (dp[i - 2] if i >= 2 else 0) + \
                (dp[i - 5] if i >= 5 else 0)
    return dp[n]

# Timing the functions
test_values = [10, 25, 50, 100]

# Timing recursive without memoization (only for small values)
for n in [10, 25]:
    start_time = time.time()
    result = f_rec(n)
    elapsed_time = time.time() - start_time
    print(f"f_rec({n}) = {result}, time: {elapsed_time:.6f} sec")

# Timing recursive with memoization
for n in test_values:
    start_time = time.time()
    result = f_memo(n)
    elapsed_time = time.time() - start_time
    print(f"f_memo({n}) = {result}, time: {elapsed_time:.6f} sec")

# Timing iterative solution
for n in test_values:
    start_time = time.time()
    result = f_it(n)
    elapsed_time = time.time() - start_time
    print(f"f_it({n}) = {result}, time: {elapsed_time:.6f} sec")
