#/bin/python
import sys
#import time
import os.path
import subprocess

#script to disect the contents of a HMM search table and disssect on the first match to a marker
class Hit(object):
	acc = ""
	e_val = 0.0
	leng = 0
	def __init__(self, a, e, l):
		self.acc = a
		self.e_val = float(e)
		self.leng = l

file_name = sys.argv[-1]

fhandle = open(file_name, "r")

table = {}

for line in fhandle:
	if line[0] != "#":
		entry = line[0:-1].split(" ")
		entry = filter(None,entry)
		if (entry[3] in table.keys()) == False:
			table[entry[3]] = Hit(entry[0],entry[6], (int(entry[2])/int(entry[5])))
		
		else:
			if table[entry[3]].e_val >  float(entry[6]):
				table[entry[3]] = Hit(entry[0],entry[6], (int(entry[2])/int(entry[5])))
fhandle.close()

#now you have the best match for each marker in the 
db_name = file_name.split(".")
db_name = db_name[0:-2]
db_name.append('fa')
db_name.append('cidx')
db_name = '.'.join(db_name)

#for each key in hash ie maker
for i in table.keys():
#get the seq ID table[i].acc
	seq_id = table[i].acc
#i will be the aligned_marker_numbers
	marker = i.split("_")[1]
	
#syscall to do cdbyank of the peptide from and append to marker file
#	print "cdbyank "+seq_id+" >>"+marker+".fasta"
	#p.subProcess
	p = subprocess.Popen(["echo '\n' >>../marker_DB/marker_files/"+marker+".fasta"], shell=True)
	p = subprocess.Popen(["cdbyank -a "+seq_id+" /tmp/"+db_name+ " >>../marker_DB/marker_files/"+marker+".fasta"], shell=True)
