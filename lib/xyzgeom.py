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

def compute_distance(Coords, idx1, idx2):
   vec_12 = Coords[idx2-1] - Coords[idx1-1]
   return np.linalg.norm(vec_12)

def compute_angle(Coords, idx1, idx2, idx3):
   vec_21 = Coords[idx1-1] - Coords[idx2-1]
   vec_21 /= np.linalg.norm(vec_21)
   vec_23 = Coords[idx3-1] - Coords[idx2-1]
   vec_23 /= np.linalg.norm(vec_23)
   theta = np.arccos(np.dot(vec_21, vec_23))
   return theta

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

def bohr_to_angs(Coords):
   Coords *= 0.529177
   return Coords

def write_xyz_file(outfile, AtomList, Coords):
   fw = open(outfile, 'w')
   fw.write("%d\n" %len(AtomList))
   fw.write("\n")
   for iAtom in range(0, len(AtomList)):
      print Coords[iAtom]
      x, y, z = Coords[iAtom][0], Coords[iAtom][1], Coords[iAtom][2]
      fw.write("%-3s %15.10f %15.10f %15.10f\n" %(AtomList[iAtom], x, y, z))
   fw.close()

def set_origin(Coords, idx_orig):
   print "Set the position of Atom %d as the origin" %idx_orig
   origin = np.copy(Coords[idx_orig-1])
   for i in range(0, Coords.shape[0]):
      Coords[i] -= origin

def set_xy_plane(Coords, idx1, idx2, idx3):
   print "Set the plane determined by atoms %d, %d, %d as the xy-plane" %(idx1, idx2, idx3) 
   A0, A1, A2 = Coords[idx1-1], Coords[idx2-1], Coords[idx3-1]
   dx = A1 - A0
   dy = A2 - A0
   xprime = dx/np.linalg.norm(dx)
   zprime = np.cross(dx, dy)
   zprime /= np.linalg.norm(zprime)
   yprime = np.cross(zprime, xprime)
   rotmat = np.zeros((3,3))
   rotmat[0] = xprime
   rotmat[1] = yprime
   rotmat[2] = zprime
   rotmat = rotmat.transpose()
   Coords = np.dot(Coords, rotmat)
   return Coords

def axis_translate(Coords, distance, axis='z'): 
   print "Translate atoms along the %s axis by %.4f Angstrom" %(axis.lower(), distance)
   if axis.lower() == 'x':
      Coords[:, 0] += distance
   elif axis.lower() == 'y':
      Coords[:, 1] += distance
   elif axis.lower() == 'z':
      Coords[:, 2] += distance
   return Coords

def xy_rotate(Coords, angle):
   print "rotate in the xy-plane by %f degrees" %angle
   rotmat = np.zeros((3,3))
   rotmat[0] = np.array([np.cos(angle/180.0*np.pi), -np.sin(angle/180.0*np.pi), 0.0])
   rotmat[1] = np.array([np.sin(angle/180.0*np.pi), np.cos(angle/180.0*np.pi), 0.0])
   rotmat[2] = np.array([0.0, 0.0, 1.0])
   Coords = np.dot(Coords, rotmat)
   return Coords

def xy_rotate_nfold(Coords, fold):
   print "rotate in the xy-plane by 2pi/%d" %fold
   n = float(fold) #n-fold rotational axis
   rotmat = np.zeros((3,3))
   rotmat[0] = np.array([np.cos(2*np.pi/n), -np.sin(2*np.pi/n), 0.0])
   rotmat[1] = np.array([np.sin(2*np.pi/n), np.cos(2*np.pi/n), 0.0])
   rotmat[2] = np.array([0.0, 0.0, 1.0])
   Coords = np.dot(Coords, rotmat)
   return Coords

def reflection(Coords, mirror_pos, normal_axis='z'):
   if normal_axis.lower() == 'z':
      print "reflection with respect to the mirror at z = %.2f" %mirror_pos
      Coords[:,2] = 2*mirror_pos - Coords[:,2]
   elif normal_axis.lower() == 'x':
      print "reflection with respect to the mirror at x = %.2f" %mirror_pos
      Coords[:,0] = 2*mirror_pos - Coords[:,0]
   elif normal_axis.lower() == 'y':
      print "reflection with respect to the mirror at y = %.2f" %mirror_pos
      Coords[:,0] = 2*mirror_pos - Coords[:,1]
   return Coords 

def pivot(Coords_orig, reverse=False):
   Coords = np.zeros(Coords_orig.shape)
   if reverse:
      print "pivot the axes: x->z, y->x, z->y"
      Coords[:,0] = np.copy(Coords_orig[:,1])
      Coords[:,1] = np.copy(Coords_orig[:,2])
      Coords[:,2] = np.copy(Coords_orig[:,0])
   else:
      print "pivot the axes: x->y, y->z, z->x"
      Coords[:,0] = np.copy(Coords_orig[:,2])
      Coords[:,1] = np.copy(Coords_orig[:,0])
      Coords[:,2] = np.copy(Coords_orig[:,1])
   return Coords
