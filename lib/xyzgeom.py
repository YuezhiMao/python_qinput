import os, glob, re, sys
import subprocess as sp
import numpy as np

def parse_xyz_file(filename):
   AtomList = []
   CoordList = []
   fr = open(filename, 'r')
   for line in fr.readlines():
      l = re.search('(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line)
      if l!=None:
         AtomList.append(l.group(1))
         coord = float(l.group(2)),float(l.group(3)),float(l.group(4))
         CoordList.append(coord)
   fr.close()
   Coords = np.array(CoordList)
   return AtomList, Coords

def compute_diheral(Coords, idx1, idx2, idx3, idx4):
   #the input index starts from 1
   b1 = Coords[idx2-1] - Coords[idx1-1]
   b2 = Coords[idx3-1] - Coords[idx2-1]
   b3 = Coords[idx4-1] - Coords[idx3-1]
   n1 = np.cross(b1, b2)
   n1 /= np.linalg.norm(n1)
   n2 = np.cross(b2, b3)
   n2 /= np.linalg.norm(n2)
   b2 /= np.linalg.norm(b2)
   m1 = np.cross(n1, b2)
   x = np.dot(n1, n2)
   y = np.dot(m1, n2)
   theta = np.arctan2(y,x) * 180.0/np.pi
   return theta
