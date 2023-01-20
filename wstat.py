from __future__ import print_function
import sys

import numpy as np
# import doctest; doctest.run_docstring_examples(wstat.wsem, globals())

einsum_bug = tuple(map(int, np.__version__.split(".")[:2])) < (1, 7)
# Bugs (due to np.einsum) for moments>1 and dim not () and np._version <= 1.6.1
"""
      >>> np.__version__
      '1.6.1'
      >>> x = np.arange(10).reshape(1,10)
      >>> print(np.einsum('ij,ij,ij->j',x, x,x))
      [  0   1   8  27  64 125 216 343   0  81]
      >>> print((x*x*x).sum(axis=0))
      [  0   1   8  27  64 125 216 343 512 729]
      #https://github.com/numpy/numpy/issues/3142
      #http://stackoverflow.com/questions/18365073/why-is-numpys-einsum-faster-than-numpys-built-in-functions/19612080
"""

def wmom(y, w=None, moment=1, axis=None, e=None, dim=(), keepdims=False):
   """
   Weighted moments.

   Returns sum(w*y^moment). np.einsum is used for efficient computation.
   Einsum can nicely handle dimensions.

   Parameters
   ----------
   y : array_like
      Sample values.
   w : array_like, optional
      An array of weights associated with the values in y.

   Examples
   --------
   >>> a = np.arange(12.).reshape(3,4)
   >>> wstat.wmom(a, moment=0)
   12.0
   >>> wstat.wmom(a, w=a*1., moment=(0,1,2))
   [66.0, 506.0, 4356.0]
   >>> wstat.wmom(a, w=a*1., moment=(0,2), dim=0)
   [array([  6.,  22.,  38.]), array([   36.,   748.,  3572.])]

   >>> a = np.arange(10.).reshape(10,1)
   >>> wstat.wmom(a, w=a*1., moment=2, axis=1)
   array([   0.,    1.,    8.,   27.,   64.,  125.,  216.,  343.,  512.,  729.])
   >>> (a*a*a).sum(axis=1)
   array([   0.,    1.,    8.,   27.,   64.,  125.,  216.,  343.,  512.,  729.])

   """
   y = np.array(y, dtype=np.float)

   d = range(y.ndim)
   if axis is not None:
      # convert axis to dim; create complement
      if isinstance(axis, int): axis = (axis,)
      axis = [d[a] for a in axis]   # ensure positive index
      dim = set(d) - set(axis)
   else:
      # projection onto axis
      if isinstance(dim, int): dim = (dim,)
      dim = [d[a] for a in dim]   # positive index

   if w is None:
      w = np.ones_like(y) if e is None else 1./e**2

   scalar = isinstance(moment, int)
   if scalar: moment = [moment]
   # compute (w*y^m).sum(axis)
   if einsum_bug and dim:
      m = [np.einsum(w, d, y**i, d, dim) for i in moment]
   else:
      m = [np.einsum(w, d, *(y,d)*i+(dim,)) for i in moment]

   if keepdims:
      # insert one at those position
      kdim = [(ni if i in dim else 1) for i,ni in enumerate(y.shape)]
      m = [x.reshape(kdim) for x in m]
      ## wstat.wmom(range(10), keepdims=True)  does not work for dim=(), i.e. scalar moments
      #s = [(slice(None) if i in dim else None) for i in d]
      #m = [x[s] for x in m]

   return m[0] if scalar else m

def wmean(y, w=None, axis=None, dim=None):
   """
   Weighted mean.

   Parameters
   ----------
   y : array_like
      Sample values.
   w : array_like, optional
      An array of weights associated with the values in y.
   axis : None or int or tuple of ints, optional
      Axes over which to average y. These axis are collapsed.
   dim : None or int or tuple of ints, optional
      Axis along which to average y. These axis will remain. The dim keyword is complementary to the axis keyword.
      It is ignored if axis is given.

   Notes
   -----
   Due to use of np.dot and/or np.einsum should be faster than np.average.
   Also it support the dim keyword.

   Examples
   --------
   >>> a = np.arange(30).reshape(2,5,3)
   >>> wstat.wmean(a); np.mean(a)
   14.5
   14.5
   >>> wstat.wmean(a,axis=1); np.mean(a,axis=1)
   array([[  6.,   7.,   8.],
          [ 21.,  22.,  23.]])
   array([[  6.,   7.,   8.],
          [ 21.,  22.,  23.]])
   >>> wstat.wmean(a,w=a); np.average(a, weights=a)
   19.666666666666668
   19.666666666666668
   >>> wstat.wmean(a, w=a, axis=1); np.average(a, weights=a, axis=1)
   array([[  9.        ,   9.57142857,  10.25      ],
          [ 21.85714286,  22.81818182,  23.7826087 ]])
   array([[  9.        ,   9.57142857,  10.25      ],
          [ 21.85714286,  22.81818182,  23.7826087 ]])

   Using the dim keyword:

   >>> wstat.wmean(a, dim=0); wstat.wmean(a, w=a*0+1, dim=0); np.mean(a, axis=1).mean(axis=1)
   array([  7.,  22.])
   array([  7.,  22.])
   array([  7.,  22.])
   >>> wstat.wmean(a, dim=(0,2)); wstat.wmean(a, w=a*0+1, dim=(0,2)); np.mean(a, axis=1)
   array([[  6.,   7.,   8.],
          [ 21.,  22.,  23.]])
   array([[  6.,   7.,   8.],
          [ 21.,  22.,  23.]])
   array([[  6.,   7.,   8.],
          [ 21.,  22.,  23.]])

   """
   if w is None and dim is None:
      return np.mean(y, axis=axis)

   if dim is None and axis is None:
      # use the dot product over all axis
      wysum = np.dot(w.ravel(), y.ravel())
      wsum = float(w.sum())
   else:
      d = range(y.ndim)
      if axis is not None:
         # convert axis to dim; create complement
         if isinstance(axis, int): axis = (axis,)
         axis = [d[a] for a in axis]    # ensure positive index
         dim = set(d) - set(axis)   # [i for i in d if i not in axis]   # create complement; convert axis to dim
         #wsum = np.sum(w, axis=axis).astype(float) # only np.vxx supports tuple for axis
         #return np.average(y, weights=w, axis=axis)
      else:
         # projection onto axis
         if isinstance(dim, int): dim = (dim,)
         dim = [d[a] for a in dim]   # positive index

      if w is None:
         wysum = np.einsum(y, d, dim)
         wsum = float(y.size / wysum.size)
      else:
         wysum = np.einsum(w, d, y, d, dim)
         wsum = np.einsum(w, d, dim).astype(float)

   return wysum / wsum

def wsem(y, mean=None, rescale=True, ddof=1, keepdims=False, **kwargs):
   """
   Standard error of weighted mean.

   Returns the the weighted mean and its uncertainty as a tuple

   Parameters
   ----------
   w : Weights
   e : One sigma error estimates.

   Examples
   --------
   Compare with results from scipy.stats.sem

   >>> a = np.arange(20).reshape(5,4)
   >>> wstat.wsem(a, ddof=0)
   (9.5, 1.2893796958227628)
   >>> wstat.wsem(a, axis=0)
   (array([  8.,   9.,  10.,  11.]), array([ 2.82842712,  2.82842712,  2.82842712,  2.82842712]))
   >>> a = np.array([55])
   >>> wstat.wsem(a, e=0.2*a)
   (55.0, 11.0)

   """
   kwargs['keepdims'] = keepdims or kwargs.get('dim') or kwargs.get('axis')  # for broadcasting of mean

   wsum, wy = wmom(y, moment=(0,1), **kwargs)
   mean = wy / wsum
   var_mean = 1. / wsum

   if rescale:
      # Rescale for under or overdispersion
      dof = float(y.size / mean.size)
      if dof > 1:
         if ddof: dof -= ddof
         # to broadcast the mean we need to keep the dimensions
         var_mean = var_mean * wmom(y-mean, moment=2, **kwargs) / dof

   if kwargs['keepdims'] and not keepdims:
   #if not keepdims:
      mean = mean.squeeze()   # np.float64(190.).reshape((1,1,)).squeeze()
                              # returns an array of zero dimension: array(190.0)
      var_mean = var_mean.squeeze()

   return mean, np.sqrt(var_mean)


def wrms(y, w=None):
   """
   root mean square (= quadratic mean)

   rms = sqrt(mean(x**2)) (not centered)

   Parameters
   ----------
   w is the weight = 1./yerr**2

   Examples
   --------
   >>> x = np.random.normal(size=1000)
   >>> wstat.wrms(x), wstat.wrms(x, w=x)

   """
   W, quadsum = (len(y), np.dot(y,y)) if w is None else (
                 np.sum(w), np.einsum('i,i,i', w,y,y))
   return np.sqrt(quadsum/W)

quadmean = rms = wrms   # quadratic mean

def wstd(y, e, axis=None, dim=(), ret_err=False):
   """
   Compute the standard deviation along the specified axis.

   Returns the weighted standard deviation (centered, biased) and weighted mean
   wstd = sqrt(wmean((x - x.wmean())**2))

   Returns
   -------
   weighted_standard_deviation: tuple(wstd, wmean [, ret_err])

   Examples
   --------
   >>> import numpy as np; import wstat
   >>> x = np.random.normal(size=(100,50))
   >>> wstat.wstd(x, x*0+1, axis=0)
   >>> wstat.wstd(x, x*0+1,dim=1, ret_err=True)

   """
   w = np.zeros_like(e, dtype=np.float)
   with np.errstate(invalid='ignore'):
       ind = e > 0
   w[ind] = 1. / e[ind]**2

   d = range(y.ndim)
   if axis is not None:
      if isinstance(axis, int): axis = (axis,)
      axis = [d[a] for a in axis]    # ensure positive index
      dim = [i for i in d if i not in axis]   # create complement

   if isinstance(dim, int): dim = (dim,)
   s = None
   if dim:
      dim = [d[a] for a in dim]   # positive index
      s = tuple(slice(None) if (a in dim) else None for a in d)   # new shape (to broadcast mean)

   with np.errstate(divide='ignore'):
      nsum = np.einsum(ind.astype(float), d, dim)
      wsum = np.einsum(w, d, dim).astype(float)
      wmean =  np.einsum(w, d, y, d, dim) / wsum
      res = y - (wmean[s] if s else wmean)       # centering data
      #wstd1 = (np.einsum(w, d, res, d, res, d, dim)/wsum)**.5   # does not work in 1.6.1 due to einsum bug
      wstd1 = (np.einsum(w, d, res*res, d, dim) / wsum)**.5
      out = (wstd1, wmean)

   if ret_err:  # append mean error
      out += ((nsum/wsum)**.5,)

   return out

def wstd_new(y, mean=None, ddof=1, keepdims=False, **kwargs):
   """
   Weighted standard deviation.

   Returns the the weighted mean and its uncertainty as a tuple

   Parameters
   ----------
   w : Weights
   e : One sigma error estimates

   Examples
   --------
   Compare with results from scipy.stats.sem

   >>> a = np.arange(20).reshape(5,4)
   >>> wstat.nanwstd(a);  np.std(a, ddof=1)
   5.9160797830996161
   5.9160797830996161
   >>> wstat.wsem(a, axis=0)
   (array([  8.,   9.,  10.,  11.]), array([ 2.82842712,  2.82842712,  2.82842712,  2.82842712]))
   >>> a = np.array([55])
   >>> wstat.wsem(a, e=0.2*a)
   (55.0, 11.0)

   """
   kwargs['keepdims'] = keepdims or kwargs.get('dim') or kwargs.get('axis')  # for broadcasting of mean

   wsum, wy = wmom(y, moment=(0,1), **kwargs)
   mean = wy / wsum
   var_mean = 1. / wsum

   dof = wmom(y, moment=0, **kwargs)
   if dof > 1:
      # to broadcast the mean we need to keep the dimensions
      var_mean = var_mean * wmom(y-mean, moment=2, **kwargs) * dof / (dof-ddof)

   if kwargs['keepdims'] and not keepdims:
      mean = mean.squeeze()   # np.float64(190.).reshape((1,1,)).squeeze()
                              # returns an array of zero dimension: array(190.0)
      var_mean = var_mean.squeeze()

   return np.sqrt(var_mean)


def wstd_v00(y, e, axis=None, ret_err=False):
   """ Compute the standard deviation along the specified axis.
   Returns the weighted standard deviation (centered, biased) and weighted mean
   wstd = sqrt(wmean((x - x.wmean())**2))

   Returns
   -------
   weighted_standard_deviation: tuple(wstd, wmean [, ret_err])
   """
   if y.ndim==1:
      e = np.array(e)
      w = np.zeros_like(e, dtype=np.float)
      ind = e > 0
      w[ind] = 1. / e[ind]**2
      #if len(y)==1: return y,e
      W = np.sum(w)
      if W>0:
         wmean = np.dot(w,y) / W
         wstd1 = (np.dot(w,(y-wmean)**2)/W)**.5, wmean
         if ret_err:  wstd1 += ((np.sum(ind)/W)**.5,)  # append mean error
      else:
         wstd1 = (0, 0, 0) if ret_err else (0, 0)
   elif axis==None:
      wstd1 = wstd_v00(y.reshape(-1),e.reshape(-1),ret_err=ret_err)
   elif axis==0:
      #wstd1 = map(np.array, zip(*(wstd_v00(yi,ei,ret_err=ret_err) for yi,ei in zip(y,e))))
      wstd1 = map(np.array, zip(*map(wstd_v00, y,e)))
      #pause()
   elif axis==1:
      wstd1 = map(np.array, zip(*(wstd_v00(yi,ei,ret_err=ret_err) for yi,ei in zip(y.T,e.T))))

   return wstd1

def quantile(x, p, w=None, middle=False):
   """
   Quantile.

   Weights are possible. Returned values are discrete, meaning a member of the sample.

   Parameters
   ----------
   x : array_like
      The sample values.
   p : array_like
      The fractions in the CDF.
   w : array_like, optional
      Weights for sample values.

   Examples
   --------
   >>> x = np.arange(10)
   >>> p = np.arange(1000)/10.
   >>> print np.median(x), quantile(x, p=[0, 0.5, 1.0])

   Visualise the CDF and the difference between quantile (returning discrete values)
   and np.percentile (returning )

   >>> gplot(quantile(x, p=p/100),p, 'w lp,', np.percentile(x, list(p)),p)

   """
   ii = np.argsort(x)

   if w is None:
      i = np.multiply(p, len(x)).astype(int)
   else:
      cdf = np.cumsum(w[ii])         # CDF of the re-arranged weights
      cdf = cdf / float(cdf[-1])       # normalise, alternatively p could be scaled
      i = np.searchsorted(cdf, p, side='right') # if p==cdf_k the lower value is returned
                               #; so in even case the upper value is returned  by default
   scalar = np.isscalar(i)
   if scalar: i = [i] # workaround for np.take in numpy<1.8.0

   quantile = x[ii.take(i, mode='clip')]

   if middle:
      edge = i%1.0 == 0. if w is None else cdf[ix-1] == p

   if scalar:   # workaround for numpy<1.8.0
      quantile = quantile[0]
   return quantile

def iqr(x, w=None, sigma=False):
   """
   Interquartile range.

   Parameters
   ----------
   x : array_like
      The sample values.
   w : array_like, optional
      Weights for sample values.
   sigma : boolean
      Convert to 1-sigma estimate.

   """
   q = quantile(x, [0.25,0.75], w=w)
   iqr = q[1] - q[0]
   if sigma:
      iqr /= 1.349
   return iqr

def mad(data, axis=None, sigma=False):
   """
   Median absolute deviation.

   Parameters
   ----------
   x : array_like
      The sample values.
   w : array_like, optional
      Weights for sample values.
   sigma : boolean
      Convert to 1-sigma estimate.

   """
   mad = np.median(np.absolute(data - np.median(data, axis)), axis)
   if sigma:
      mad *= 1.4826
   return mad


def wnan_to_num(y, w=None, e=None):
   """
   Preparation for nan functions.

   Detects nans in y and w (or e). Returns nan replaced by zero and weights set to zero.
   Negative weights are not tested. Only errors greater than zero will have non-zero weights.

   Examples
   --------
   >>> y = np.array([np.nan, 0, 1, 2, 3])
   >>> e = np.array([np.nan, np.nan, -1, 0, 2])
   >>> y, e
   (array([ nan,   0.,   1.,   2.,   3.]), array([ nan,  nan,  -1.,   0.,   2.]))
   >>> wstat.wnan_to_num(y, e=e)
   (array([ 0.,  0.,  1.,  2.,  3.]), array([ 0.  ,  0.  ,  0.  ,  0.  ,  0.25]))

   """
   if w is None:
      if e is not None:
         w = np.zeros(y.shape)
         ind = e > 0
         w[ind] = 1. / e[ind]**2
   else:
      w = np.nan_to_num(w)

   w = np.isfinite(y) * (1. if w is None else w)
   y = np.nan_to_num(y)
   return y, w

def nanwsem(y, w=None, e=None, **kwargs):
   """
   Standard error of weighted mean.

   Returns the the weighted mean and its uncertainty as a tuple

   Parameters
   ----------
   w : Weights
   e : One sigma error estimates. Only positive errors are used.

   Examples
   --------
   Compare with results from scipy.stats.sem

   >>> a = np.arange(20).reshape(5,4)
   >>> wstat.wsem(a, ddof=0)
   (array(9.5), 1.2893796958227628)
   >>> wstat.wsem(a, axis=0)
   (array([  8.,   9.,  10.,  11.]), array([ 2.82842712,  2.82842712,  2.82842712,  2.82842712]))
   >>> a = np.array([55])
   >>> wstat.wsem(a, e=0.2*a)
   (55.0, 11.0)

   """
   y, w = wnan_to_num(y, w=w, e=e)
   return wsem(y, w=w, **kwargs)

def nanwstd(y, w=None, e=None, **kwargs):
   """
   Compute the standard deviation along the specified axis.
   Same as wstd but with wmom

   Returns the weighted standard deviation (centered, biased) and weighted mean
   wstd = sqrt(wmean((x - x.wmean())**2))

   Returns
   -------
   weighted_standard_deviation: tuple(wstd, wmean [, ret_err])

   Examples
   --------
   >>> y = np.array([np.nan, 0, 1, 2, 3, 4, 5])
   >>> e = np.array([np.nan, np.nan, -1, 0, 2, 2, 2])
   >>> wstat.nanwstd(y, e=e)
   1.1547005383792515

   """
   y, w = wnan_to_num(y, w=w, e=e)
   return wstd_new(y, w=w, **kwargs)

def naniqr(x, w=None, e=None, **kwargs):
   """
   Interquartile range.

   Parameters
   ----------
   x : array_like
      The sample values.
   w : array_like, optional
      Weights for sample values.
   sigma : boolean
      Convert to 1-sigma estimate.

   Examples
   --------
   >>> y = np.array([np.nan, 0, 1, 2, 3, 4, 5])
   >>> e = np.array([np.nan, np.nan, -1, 0, 2, 2, 2])
   >>> wstat.naniqr(y, e=e, sigma=True)
   1.4825796886582654

   """
   x, w = wnan_to_num(x, w=w, e=e)
   return iqr(x, w=w, **kwargs)

def mlrms(y, e, s=0., verbose=False, ml=True, ret_mean=False):
   '''
   Maximum likelihood rms.

   Parameters
   ----------
   y :  array_like
      The sample values.
   e :  array_like
      One sigma error estimates.
   s : float, optional
      Start guess for jitter.
   ml : boolean, optional
      If true, iterate towards maximum likelihood. If false, iterate for chi_red=1.
   ret_mean : boolean, optional
      If true, the reweighted mean is returned.

   Returns
   -------
   The reweighted rms and the unknown variance.

   '''
   n = y.size
   i = 0
   eps = .000001
   while True:
      i += 1
      w = 1 / (e**2+s**2)          # updated weights
      W = w.sum()
      q = 1 / np.sqrt(w.mean())    # "mean total jitter"
      Y = np.dot(w, y) / W         # the weighted mean of y
      r = y - Y

      chi2 = np.sum(w*r**2)
      wrms = np.sqrt(chi2 / W)
      if ml:
         # iterate for maximum lnL
         wwrr = np.sum(w*w*r**2)
         #print(np.sum(w),np.sum(r),np.sum(y))
         #s = np.sqrt(np.clip(np.sum(w*w*(r**2-e**2)) / np.sum(w*w),0,None))
         s = np.sqrt((np.sum(w*w*(r**2-e**2)) / np.sum(w*w)).clip(min=0))
         rr = wwrr / W
      else:
         # iterate for chi_red=1
         s = np.sqrt((s**2+wrms**2-q**2).clip(min=0)) # add to the jitter the remaining difference between external and internal scatter
         rr = wrms / q
      #print(-eps<rr-1<eps , s==0, wwrr, np.sum(w*w*(r**2-e**2)) , np.sum(w*w))
      lnL = -0.5 * np.sum(np.log(2*np.pi/w)) - 0.5 * chi2

      if verbose: print('mean %.5g' %Y,' err', q, ' mlrms', wrms, lnL, ' rchi', chi2/n, rr, ' jit', s, s/wrms, wwrr,W)
      #if verbose: print('s',s,'wwrr', wwrr, 'r',rr, 'W', W)
      if -eps<rr-1<eps or s==0 or i>20:
         if ret_mean:
            return wrms, s, Y
         return wrms, s


def speed_comparison():
   import timeit
   setup = '''
import numpy as np
import wstat
a = np.arange(3000).reshape(20,5,-1)
'''
   print("wstat.wmean",min(timeit.Timer('wstat.wmean(a, w=a, axis=1)', setup=setup).repeat(7, 1000)))
   print("np.average ",min(timeit.Timer('np.average(a, weights=a, axis=1)', setup=setup).repeat(7, 1000)))

