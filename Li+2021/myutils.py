""" These utility functions may or may not be related to the oPDF method. Some of them are just general-purpose plotting or monitoring functions.
"""
import sys
import numpy as np
from scipy.stats import gaussian_kde,norm,chi2
import matplotlib
#matplotlib.user('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import MaxNLocator # added 

def create31fig(sharex=True, sharey=False, figsize=(8,8)):
    '''create a figure with 3 tightly packed subplots'''
    f,ax = plt.subplots(3, sharex=sharex, sharey=sharey, figsize=figsize)
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    nbins = len(ax[0].get_yticklabels())
    plt.setp([a.yaxis for a in ax[:-1]], major_locator=MaxNLocator(nbins=nbins, prune='lower',symmetric=True))
    return ax

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    
      cmap : The matplotlib colormap to be altered
      
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
      
      Credit: http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap
  
def plot_circle(cen=[0,0], r=1, **kwargs):
    '''plot a circle'''
    phi=np.arange(0,2*np.pi+0.11,0.1)
    x=cen[0]+r*np.cos(phi)
    y=cen[1]+r*np.sin(phi)
    h=plt.plot(x,y,**kwargs)
    return h

def lnADPDF(lnAD):
  gpar=[[0.569,   -0.570,    0.511], [ 0.431,    0.227,    0.569]] #w, mu, sigma, bi-normal fit
  p=gpar[0][0]*norm.pdf(lnAD,gpar[0][1],gpar[0][2])+gpar[1][0]*norm.pdf(lnAD, gpar[1][1], gpar[1][2])
  return p
  
def ADSurvFunc(AD):
  '''survival function for AD test on uniform distributions'''
  gpar=[[0.569,   -0.570,    0.511], [ 0.431,    0.227,    0.569]] #w, mu, sigma, bi-normal fit
  lnAD=np.log(AD)
  p=gpar[0][0]*norm.sf(lnAD,gpar[0][1],gpar[0][2])+gpar[1][0]*norm.sf(lnAD, gpar[1][1], gpar[1][2])
  return p
  #return norm.sf(np.log(AD),loc=-0.22,scale=0.66);

def P2Sig(pval):
  """convert pval to sigma"""
  #return norm.ppf(1.0-pval/2)
  return -norm.ppf(pval/2) #equivalent to the above, but more accurate

def AD2Sig(AD):
  """convert AndersonDarling TS to sigma"""
  AD=np.array(AD)
  sig=P2Sig(ADSurvFunc(AD))
  #sig[AD>5]=(np.log(AD[AD>5])+0.22)/0.66
  return sig

def Sig2TS(sig,dof=1):
  '''convert a sigma value to the 2*likelihood ratio'''
  return chi2.ppf(chi2.cdf(np.array(sig)**2,1),dof)

def Chi2Sig(x, dof):
  '''convert chi-square value to significance level, for dof degrees of freedom'''
  return P2Sig(chi2.sf(x,dof))
  
class ProgressMonitor:
    """monitor progress of your loops"""
    def __init__(self,total_steps, total_show=100, init_show=0):
        """init_show: initial value progress percentage, set to 0 if no reason
        total_steps: maximum iteration steps to monitor
        total_show: number of revealing times
        """
        self.current_show=int(init_show)
        self.total_steps=total_steps
        self.total_show=total_show
        print(" %02d%%"%(self.current_show*100/self.total_show), end=' ')
        sys.stdout.flush()
        
    def monitor_progress(self,current_step):
        """put this inside loop to monitor progress, to print the percent of
        job finished."""
        #print when current_step first exceeds current show percent
        if current_step>=self.total_steps*self.current_show/self.total_show:
            print("\b\b\b\b\b %02d%%"%(self.current_show*100/self.total_show), end=' ')
            sys.stdout.flush()
            self.current_show+=1
            
def percent2level(p,z):
    ''' convert percentiles to levels '''
    try:
      p=list(p)
    except: #single number
      p=[p]
    x=np.sort(np.array(z.ravel()))
    x=x[::-1]
    frac=x.cumsum()/x.sum()
    l=[x[abs(frac-pi).argmin()] for pi in p]
    return l

def contour_handle(color, linestyle='solid'):
  '''return a patch object to be used for labelling patch objects in legends'''
  return Ellipse((0,0),0,0,fill=False, color=color, linestyle=linestyle)  

def density_of_points(data, bins=100, method='kde', weights=None):
  ''' estimate density of points with kde or histogram2d 
  
  data: should be shape [2,n] array
  
  bins: can be an integer or [nx,ny] for number of bins, an ndarray or a list of two arrays  for bin edges
  
  method: 'kde' or 'hist', kernel-density-estimate or 2d-histogram estimate
  
  weights: whether to use weights or not. currently only supports hist method.
  
  return: (X,Y,Z)
      ready to be used for contour plots as contour(X, Y, Z). 
      X and Y are mid points of the bins on which Z is calculated.
  '''
  if data.shape[0]!=2 and data.shape[1]==2:
    data=data.T
  l=data.min(axis=1)
  r=data.max(axis=1)
  if isinstance(bins, int):  
      x=np.linspace(l[0], r[0], bins+1)
      y=np.linspace(l[1], r[1], bins+1)
  elif isinstance(bins, np.ndarray):
      x=bins
      y=bins
  elif isinstance(bins, list):
    if isinstance(bins[0], int):
      x=np.linspace(l[0], r[0], bins[0]+1)
      y=np.linspace(l[1], r[1], bins[1]+1)
    else:  
      x=bins[0]
      y=bins[1]
  
  X,Y=np.meshgrid((x[:-1]+x[1:])/2, (y[:-1]+y[1:])/2) #mid points
  if method=='kde':
      positions =np.vstack([X.ravel(), Y.ravel()])
      kernel = gaussian_kde(data)
      Z = np.reshape(kernel(positions).T, X.shape)
  else:
      Z=np.histogram2d(data[0], data[1], bins=[x, y], weights=weights)
      Z=Z[0].T
  return X,Y,Z
      
def percentile_contour(X,Y,Z, percents=0.683, colors=None, fill=False, linestyles='solid', **kwargs):
  """
  plot contour at specific percentile levels

      X,Y can be both 2-d arrays as Z, or 1-d array specifying the column(horizontally varying) and row coordinates for Z.
  
      percents can be a list, specify the contour percentile levels
  
      colors should be a tuple, e.g, (r,) 
  
      fill: bool, whether to plot filled contours
  
      kwargs specify linestyles
  
  return:
      a handle artist of the same linestyle (but not the contour object) to be used in legends
  """    
  #if type is 'image':
    #extent=(x.min()-(x[1]-x[0])/2, x.max()+(x[1]-x[0])/2, y.min()-(y[1]-y[0])/2, y.max()+(y[1]-y[0])/2)
    #if colors is None:
      #colors=plt.cm.summer
    #h0=plt.imshow(Z, extent=extent, cmap=colors)
  #else:
  lvls=percent2level(percents,Z)
  if fill:
    h0=plt.contourf(X,Y,Z,lvls, colors=colors, linestyles=linestyles, **kwargs)
  else:
    h0=plt.contour(X,Y,Z,lvls, colors=colors, linestyles=linestyles, **kwargs)

  try:
    color=list(colors)[0]
  except:
    color=colors
  h=Ellipse((0,0),0,0,fill=fill, color=color, linestyle=linestyles, **kwargs)
  
  return h,h0,lvls

def get_extent(X,Y):
  ''' get extent for X,Y vectors or meshgrids. 
  
  the output is (xmin,xmax,ymin,ymax), the edge-padded boudaries,
      assuming X,Y specifies the mid points of bins and uniformly spaced.
      
  can be used to specify extent for imshow()'''
  if X.squeeze().ndim==2:
    dx=X[0,1]-X[0,0]
    dy=Y[1,0]-Y[0,0]
  else:
    dx=X[1]-X[0]
    dy=Y[1]-Y[0]
  extent=(X.ravel().min()-dx/2, X.ravel().max()+dx/2, Y.ravel().min()-dy/2, Y.ravel().max()+dy/2)
  return extent

def plot_cov_ellipse(cov, pos, nstd=1, fill=False, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    
        cov : The 2x2 covariance matrix to base the ellipse on
        
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 1 standard deviations.
        
        ax : The axis that the ellipse will be plotted on. Defaults to the current axis.
        
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    TS={1:2.3, 2:6.18, 3:11.8} #the TS=x'* at the specific nstd, for a 2-d
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * np.sqrt(TS[nstd]*vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, fill=fill, **kwargs)
    #ellip.set_facecolor('none')
    ax.add_artist(ellip)
    return ellip
  
def skeleton(x,y,nbin=10,alpha=0.683,weights=None):
    """
    to divide x into bins and give estimation of center and variance of y inside each bin
    
    input:
        x,y: column vectors to extract skeleton from
        
        nbin: number of bins or bin edges for x
        
        alpha: confidence level for boundary estimation
    """

    x=np.array(x)
    y=np.array(y)
    
    count,xbin=np.histogram(x,nbin,weights=weights)
    nbin=len(xbin)-1
    bin=np.digitize(x,xbin)-1
    
    xm=np.empty(nbin)
    #ym=xm[:]     #this is wrong! even though id(ym)!=id(xm), and id(ym[0])!=id(xm[0])
    ym=np.empty_like(xm)
    ysig=np.empty_like(xm)
    xmed=np.empty_like(xm)
    ymed=np.empty_like(xm)
    ylim=np.empty([2,nbin]);
    alpha=(1-alpha)/2;
    
    for i in range(nbin):
        if weights is not None:
          xm[i]=np.sum(x[bin==i]*weights[bin==i])/np.sum(weights[bin==i])
          ym[i]=np.sum(y[bin==i]*weights[bin==i])/np.sum(weights[bin==i])
        xm[i]=np.mean(x[bin==i])
        xmed[i]=np.median(x[bin==i])
        ym[i]=np.mean(y[bin==i])
        ymed[i]=np.median(y[bin==i])
        ysig[i]=np.std(y[bin==i])
        if np.sum(bin==i):
            ylim[:,i]=np.percentile(y[bin==i], [alpha*100, (1-alpha)*100])
        else:
            ylim[:,i]=[np.NaN,np.NaN]
            
    return {'x':{'median':xmed,'mean':xm,'bin':xbin,'hist':count},'y':{'median':ymed,'mean':ym,'std':ysig,'CI':ylim}}

class Enum(object):
  def __init__(self, names):
    for number, name in enumerate(names.split()):
      setattr(self, name, number)
      
class NamedValues(object):
  def __init__(self, value, name):
    self.value=value
    self.name=name
  def __repr__(self):
    return self.name
  def __str__(self):
    return self.name  

class NamedEnum(object):
  def __init__(self, names):
    '''create a named enum object from a list of empty-space separated names'''
    for number, name in enumerate(names.split()):
      setattr(self, name, NamedValues(number, name))

#import ctypesGsl as cgsl
#def fmin_gsl(func, x0, args=[], xtol=1e-3, ftolabs=0.01, xstep=1.0, maxiter=1000, full_output=False):
    #'''
    #minimize function with gsl_simplex method
    #func(x [,args]): function to be minimized
    #x0: [a,b,c...], initial parameter
    #args: list of additional parameter if any
    #xstep: initial simplex size
    #mimics scipy.optimize.fmin() interface
    #this fmin() is faster and more accurate than the scipy.optimize.fmin(), also better than fmin_powell() in scipy.
    #'''
    #args=list(args)
    #if args is []:
      #myfunc=lambda x,arg: func(x)
    #else:
      #myfunc=lambda x,arg: func(x, *arg)
    #F = cgsl.gsl_multimin_function(myfunc, len(x0), args)
    #x = cgsl.vector(x0)
    #T = cgsl.multimin_fminimizer_nmsimplex
    #s = cgsl.multimin_fminimizer(T, F)
    #s.init(x, cgsl.vector([xstep] * F.n))
    #it = 0
    #f1=F(x)
    #while True:
        #it += 1
        #s.iterate()
        #f0=f1
        #f1=s.minimum()
        #status = s.test_size(xtol)
        #xx = s.x()
        #if status and abs(f1-f0)<ftolabs:
            #print "Optimization terminated successfully."
            #print "\t Current function value: ", f1
            #print "\t Iterations: ", it
            #print "\t x abs err: ", s.size()
            #print "\t", xx
            #status=0
            #break
        #if it >= maxiter:
            #print "Maximum number of %d iterations reached"%maxiter
            #print "Failed to converge"
            #status=1
            #break
    #x=np.array([xx[i] for i in xrange(F.n)])
    #if full_output:
      #return x,f1,it,status
    #else:
      #return x      
