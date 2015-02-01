import scipy.stats
import math

def cdf_multi(l1, l2, u1, u2, mu, cov):
	stdev = math.sqrt(cov)
	l1 = (l1-mu)/stdev
	l2 = (l2-mu)/stdev
	u1 = (u1-mu)/stdev
	u2 = (u2-mu)/stdev
	correl = [0.0]
	infin = [2,2]
	error, cdfvalue, inform = scipy.stats.kde.mvn.mvndst([l1,l2], [u1,u2], infin, correl)
	return cdfvalue