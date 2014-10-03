import random
#import matplotlib.pyplot as plt
#http://matplotlib.org/api/pyplot_api.html


import numpy as np

#
# this modified version of the scipy.stats.kde.gaussian_kde() relies
#   only on numpy. simply using "import numpy as np" will enable any script
#   to use this class.
#
class gaussian_kde(object):
    """
    Representation of a kernel-density estimate using Gaussian kernels.



    Attributes
    ----------
    d : int
        number of dimensions
    n : int
        number of datapoints

    Methods
    -------
    kde.evaluate(points) : array
        evaluate the estimated pdf on a provided set of points
    kde(points) : array
        same as kde.evaluate(points)
    kde.resample(size=None) : array
        randomly sample a dataset from the estimated pdf.
    kde.covariance_factor() : float
        computes the coefficient that multiplies the data covariance matrix to
        obtain the kernel covariance matrix. Set this method to
        kde.scotts_factor or kde.silverman_factor (or subclass to provide your
        own). The default is scotts_factor.
    *For integration methods, see scipy distribution: scipy.stats.kde.gaussian_kde()

    Parameters
    ----------
    dataset : (# of dims, # of data)-array
        datapoints to estimate from

    """

    def __init__(self, dataset, bw_method="scott"):
        self.dataset = np.atleast_2d(dataset)
        
        if bw_method == "scott":
            self.covariance_factor = self.scotts_factor
        elif bw_method == "silverman":
            self.covariance_factor = self.silverman_factor
        else:
            self.covariance_factor = self.arbitrary_factor
            self.bw_method = float(bw_method)
            
        self.d, self.n = self.dataset.shape

        self._compute_covariance()
        
        
            
    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError if the dimensionality of the input points is different than
        the dimensionality of the KDE.
        """

        points = np.atleast_2d(points).astype(self.dataset.dtype)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d,
                    self.d)
                raise ValueError(msg)

        result = np.zeros((m,), points.dtype)

        if m >= self.n:
            # there are more points than data, so loop over data
            for i in range(self.n):
                diff = self.dataset[:,i,np.newaxis] - points
                tdiff = np.dot(self.inv_cov, diff)
                energy = np.sum(diff*tdiff,axis=0)/2.0
                result += np.exp(-energy)
        else:
            # loop over points
            for i in range(m):
                diff = self.dataset - points[:,i,np.newaxis]
                tdiff = np.dot(self.inv_cov, diff)
                energy = np.sum(diff*tdiff,axis=0)/2.0
                result[i] = np.sum(np.exp(-energy),axis=0)

        result /= self._norm_factor

        return result

    __call__ = evaluate

    def resample(self, size=None):
        """Randomly sample a dataset from the estimated pdf.

        Parameters
        ----------
        size : int, optional
            The number of samples to draw.
            If not provided, then the size is the same as the underlying
            dataset.

        Returns
        -------
        dataset : (self.d, size)-array
            sampled dataset
        """

        if size is None:
            size = self.n

        norm = transpose(multivariate_normal(zeros((self.d,), float),
            self.covariance, size=size))
        indices = randint(0, self.n, size=size)
        means = self.dataset[:,indices]

        return means + norm


    def scotts_factor(self):
        return np.power(self.n, -1./(self.d+4))

    def silverman_factor(self):
        return power(self.n*(self.d+2.0)/4.0, -1./(self.d+4))
    
    def arbitrary_factor(self):
        return self.bw_method

    # This can be replaced with silverman_factor if one wants to use Silverman's
    # rule for choosing the bandwidth of the kernels.
    

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor
        """
        self.factor = self.covariance_factor()
        self.covariance = np.atleast_2d(np.cov(self.dataset, rowvar=1, bias=False) *
            self.factor * self.factor)
        self.inv_cov = np.linalg.inv(self.covariance)
        self._norm_factor = np.sqrt(np.linalg.det(2*np.pi*self.covariance)) * self.n
