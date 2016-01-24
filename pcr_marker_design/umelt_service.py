##interrogate umelt service at Univ of Utah
import sys
import urllib.request, urllib.parse, urllib.error
import urllib.request, urllib.error, urllib.parse
import xml.etree.ElementTree as ET
import numpy as np
import scipy as sp
from scipy import interpolate
from scipy.interpolate import interp1d


url='https://www.dna.utah.edu/db/services/cgi-bin/udesign.cgi'
timeout_sec=500 ## default for timeout

# Query UW melt prediction service for a single sequence, returning default array for helicity, assuming within temperature range of 65-95
def getmelt(input_seq):
    values = {'seq' : input_seq, 'rs':0, 'dmso':0,'cation': 20 ,'mg': 2} # Note the buffer conditions
    data = urllib.parse.urlencode(values)
    req = urllib.request.Request(url, data.encode('utf-8'))
    try:
        response = urllib.request.urlopen(req,timeout=timeout_sec)
    except urllib.error.HTTPError as e:
        print('The server couldn\'t fulfill the request.', file=sys.stderr)
        print('Error code: ', e.code, file=sys.stderr)
    except urllib.error.URLError as e:
        print('We failed to reach a server.', file=sys.stderr)
        print('Reason: ', e.reason, file=sys.stderr)
    else:
        melt_data = response.read()
        tree = ET.fromstring(melt_data)
        helicity = [amp.find('helicity').text.split() for amp in tree.findall('amplicon')]
        hels = np.array(helicity[0], dtype=np.float32).transpose()
        return hels

# helicity[0] used because the default retreived data is a list of 3 lists (for WT, mu and hets). The 3 lists are identical for our data (no IUPAC) so only [0] is used.

def getTm(hel_array):
	temps=np.arange(65,100.5,0.5) # Temperature range of 65-100.5. Step of 0.5 NEEDS TO BE CHECKED..REcent change?
	tck = interpolate.splrep(temps,hel_array,s=0)
	xnew =np.arange(65,95.5,0.05)
	yder = interpolate.splev(xnew,tck,der=1) # der=1, first derivative
	return xnew[yder.argmin()] # Returns the x value corresponding to the minimum (peak) y value -> Tm
