from subprocess import Popen, PIPE, check_call, call
import shlex
#import pretend as p
#import pexpect
#child = pexpect.spawn('ssh csp256@193.157.210.92')
#child.expect('Password:')
#child.sendline('')
#proc.stdin.flush()  
#stdout,stderr = proc.communicate('\n')  

#subprocess.call()
#subprocess.call(['cd', 'Documents/meanders/src'])
#subprocess.call(["make", "run"])
#subprocess.call(["ffmpeg", '-y', '-r', '20', '-start_number', '0', '-filter_complex', 'scale=360:-1,tile=2x2:margin=10:padding=4', '-b:v', '320k', '-i', '../imgs/img%10d.bmp', '-vcodec', 'libx264', './../imgs/output.avi'])
#subprocess.call(["mplayer", "../imgs/output.avi"])


import sys

#call("./monte_carlo","out","<","std")

#aa = ['a','b','c'].join('\n')
#print aa
import csv
with open('std') as file:
	reader = csv.reader(file,delimiter=' ')
	aa = ''
	for row in reader:
		aa = aa + str(row[0]) + '\n'
		#print row
		#print aa

command = "./monte_carlo out 1"
args = shlex.split(command)
print args
p = Popen(args, stdout=PIPE, stdin=PIPE)
q = Popen(args, stdout=PIPE, stdin=PIPE)
std_out1 = p.communicate(input=aa)[0]
print std_out1
std_out2 = q.communicate(input=aa)[0]
print std_out2