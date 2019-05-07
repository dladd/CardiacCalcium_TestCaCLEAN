#!/usr/bin/env python
'''
Author: David Ladd

Brief: Convert finite element FCa field results to simulated confocal
       microscopy data to be processed by CaCLEAN.

Copyright 2019 David Ladd, University of Melbourne

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import re
import shutil
import requests

#=================================================================
# F u n c t i o n s
#=================================================================

def download_file(local_file, targetDir):
    """Download the local_file from the URL"""
    resourceLinks = {
        'psf_50x50x140_53nmPx_fwhm_410xy_1800z.tif' : "https://melbourne.figshare.com/ndownloader/files/15081458",
        'psf_50x50x140_53nmPx_fwhm_205xy_900z.tif' : "https://melbourne.figshare.com/ndownloader/files/15081464",
        'psf_50x50x140_53nmPx_fwhm_820xy_3600z.tif' : "https://melbourne.figshare.com/ndownloader/files/15081461",
        'Combined_8Sarc_503kNodes_surfaceFixed.stl' : "https://melbourne.figshare.com/ndownloader/files/15045725",
        'Combined_8Sarc_1319kNodes_node.h5' : "https://melbourne.figshare.com/ndownloader/files/15045695"
    }
    url = resourceLinks[local_file]
    r = requests.get(url, stream=True)
    with open(targetDir + local_file, 'wb') as f:
        shutil.copyfileobj(r.raw, f)


def findBetween( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


def readExnode(filename,info,nodeData,totalNumberOfNodes,dependentFieldNumber):
    """
    Reads fields from data produced by OpenCMISS in the .exnode format
    """
    try:
        with open(filename):
            f = open(filename,"r")

            #Read header
            line=f.readline()
            if '/' in line:
                line=f.readline()
                info.region=findBetween(line, ' Region: ', '\n')
            info.group=findBetween(line, ' Group name: ', '\n')
            numberOfFieldComponents = []
            numberOfFields = 0
            readExnodeHeader(f,numberOfFieldComponents)
            #print(numberOfFieldComponents)
            numberOfFields = len(numberOfFieldComponents)

            #Read node data
            endOfFile = False
            while endOfFile == False:
                previousPosition = f.tell()
                line=f.readline()
                #print(line)
                line = line.strip()
                if line:
                    if 'Node:' in line:
                        s = re.findall(r'\d+', line)
                        node = int(s[0])
                        if node > totalNumberOfNodes:
                            endOfFile = True
                            break
                        c = -1
                        for field in range(numberOfFields):
                            for component in range(numberOfFieldComponents[field]):                                
                                line=f.readline()
                                line = line.strip()
                                value = float(line)
                                c+=1
                                nodeData[node-1,c] = value

                    elif 'Fields' in line:
                        f.seek(previousPosition)
                        numberOfFieldComponents = []
                        numberOfFields = 0
                        readExnodeHeader(f,numberOfFieldComponents)
                        numberOfFields = len(numberOfFieldComponents)

                else:
                    endOfFile = True
                    f.close()
    except IOError:
       print ('Could not open file: ' + filename)



def readExnodeHeader(f,numberOfFieldComponents):
    #Read header info
    line=f.readline()
    if '/' in line:
        line=f.readline()
    s = re.findall(r'\d+', line)
    numberOfFields = int(s[0])
    for field in range(numberOfFields):
        line=f.readline()
        fieldName = findBetween(line, str(field + 1) + ') ', ',')
        numberOfComponents = int(findBetween(line, '#Components=', '\n'))
        numberOfFieldComponents.append(numberOfComponents)
        for skip in range(numberOfComponents):
            line=f.readline()
       
def inHull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0


def isNumber(s):
    # Check if is a real number
    try:
        int(s)
    except ValueError:
        return False
    return True

