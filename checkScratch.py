#!/user/bin/python
import os
max_host_number = 102
host_list = ["n0{:03d}.mhg".format(x) for x in range(30, max_host_number)]

print "there are {} MHG nodes".format(max_host_number)
for host in host_list:
   os.system("echo {}".format(host))
   #os.system("ssh {} \"df -h \" ".format(host))
   os.system("ssh {} \"du -sh /scratch/* \" ".format(host))
