from subprocess import Popen, PIPE, check_call, call
import shlex
import sys

filename = "out_2n_1_1_10"
command = "ffmpeg -i imgs/"+str(filename)+"_%d.png "+str(filename)+".mp4"
args = shlex.split(command)
print args
Popen(args)
print "Fin"