import pandas
import numpy as np
import matplotlib.pyplot as plt

sa=1
sb=1
r=0.7
nd=sa*sb*r
cov = [[sa**2, nd], [nd, sb**2]]


measurements=np.random.multivariate_normal(mean, cov, 10)
for measurement in measurements:
    mean = measurement
    x, y = np.random.multivariate_normal(mean, cov, 5000).T
    plt.plot(x, y, 'x')
    plt.axis('equal')
    plt.show()

 # diagonal covariance
