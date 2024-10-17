from scipy.stats import binom
import matplotlib.pyplot as plt

'''
Binomial is used for yes, no data 

k = desired value (y in FCBS)
n = sample size ()
p = parameters (θ in FCBS)
'''

y = 1
num_samples = 1000
θ = p = 0.4 #% true


# Probability Mass Function = sum all the way up until the value you are concerned about
pmf = binom.ppf(0.01, num_samples, p)
print(pmf)

# Sample random values from the binomial distribution
# Syntax: binom.rvs(max range, percentile to sample from, num samples)

for i in range(1,10):
    print(i)
    samples = binom.rvs(1, i, size=1000)
    #print(samples)
    plt.hist(samples)
    plt.show()
    plt.clf()