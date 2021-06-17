

# can vary by up to 1e-3 I believe
planck_like_expected = -6288.687206797994

# sets
planck_like = float(open("output/demo2/likelihoods/values.txt").read().split()[-1])


assert abs(planck_like - planck_like_expected) < 1e-3