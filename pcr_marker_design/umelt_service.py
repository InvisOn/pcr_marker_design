# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:29:48 2016

@author: cflbxd

This is a rewrite of the umelt_service module
for calculating melting temperatures of DNA
sequences using the online uMelt service at
the University of Utah. This rewrite is meant
to be easier to test and more object-oriented.
"""

from __future__ import print_function
import requests
import xml.etree.ElementTree as ET
import numpy as np
from scipy import interpolate

# Silence InsecureRequestWarning
requests.packages.urllib3.disable_warnings()


class MeltSeq:
    """A DNA melting experiment

    This class knows all about melting DNA
    sequences, including all
    attributes they might have and
    how to convert between different
    formats for attributes that
    other classes might want to use.
    """

    def __init__(self, sequence, resolution=0, dmso_percent=0, cations=20, free_mg=2):
        self.sequence = sequence
        self.resolution = resolution
        self.dmso_percent = dmso_percent
        self.cations = cations
        self.free_mg = free_mg


class HelicityInfo:
    """Information about a sample's helicity.

    This class knows about helicity
    charts and the data in them. If
    you give it x and y data, it
    will help you to interpret the
    data.
    """

    def __init__(self, helicity_array, temperature_range=None, min_temp=65, max_temp=100.5):
        """Create a HelicityInfo object.

        helicity_array is a numpy array of
        y data in a temperature vs. helicity
        % chart.
        temperature_range is an optional
        numpy array of x data in a
        temperature vs. helicity % chart.
        If temperature_range is not given,
        one will be
        """

        self.helicity_data = helicity_array
        self.min_temp = min_temp
        self.max_temp = max_temp

        if temperature_range is None:
            temperature_range = np.linspace(self.min_temp, self.max_temp,
                                            self.helicity_data.size)

        self.temperature_range = temperature_range

    def get_melting_temp(self):

        # interpolate as a spline
        # to increase resolution
        # s = 0 means no smoothing
        # just straight interpolation
        tck = interpolate.splrep(self.temperature_range, self.helicity_data, s=0)

        # make some new x points with 10 times
        # the resolution of our original points
        xnew = np.linspace(self.min_temp, self.max_temp, self.helicity_data.size * 10)

        # first derivative of the line
        ynew_derivative = interpolate.splev(xnew, tck, der=1)

        # return the x value corresponding to the
        # point with the steepest downward slope
        return xnew[ynew_derivative.argmin()]


class UmeltService:
    """An API for the uMelt service.

    This class knows everything about how to
    use the uMelt web service, so no other
    classes have to.
    """

    def __init__(self):
        self.url = 'https://www.dna.utah.edu/db/services/cgi-bin/udesign.cgi'
        self.timeout = 500

    def get_response(self, sequence):
        """Send a sequence to uMelt and return the response.

        sequence is a MeltSeq object.

        Returns a requests response object.
        """

        values = {'seq': sequence.sequence,
                  'rs': sequence.resolution,
                  'dmso': sequence.dmso_percent,
                  'cation': sequence.cations,
                  'mg': sequence.free_mg}

        # TODO: add in error handling
        response = requests.get(self.url, params=values, timeout=self.timeout, verify=False)
        # TODO: replace 'verify = False' with something
        # that's not a massive security hole

        return response

    def get_helicity_info(self, response):

        melt_data = response.text
        tree = ET.fromstring(melt_data)
        helicity = [amp.find('helicity').text.split() for amp in tree.findall('amplicon')]

        # helicity is a list of 3 lists, which for our
        # purposes are identical
        helicity_array = np.array(helicity[0], dtype=np.float32).transpose()
        temperature_range = np.arange(65, 100.5, 0.5)

        helicity_info = HelicityInfo(helicity_array, temperature_range)

        return helicity_info


def getmelt(input_seq):
    """A copy of the original getmelt function.

    This function takes an input sequence
    string, queries the online umelt service
    at UoU and returns the helicity array
    that results.
    """

    sequence = MeltSeq(input_seq)
    umelt = UmeltService()
    response = umelt.get_response(sequence)
    helicity = umelt.get_helicity_info(response)

    return helicity.helicity_data
