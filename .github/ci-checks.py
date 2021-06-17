

# can vary by up to 0.01 it seems
planck_like_expected = -6288.687206797994

# demo 2 check
planck_like = float(open("output/demo2/likelihoods/values.txt").read().split()[-1])

print("Planck like expected: ", planck_like_expected)
print("Planck like found: ", planck_like)
assert abs(planck_like - planck_like_expected) < 0.02