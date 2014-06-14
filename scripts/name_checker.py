#! /usr/bin/python
import subprocess
import os.path
import sys

file_handle = sys.argv[1]

if os.path.exists(file_handle):
	p_grep = subprocess.Popen(["grep -m 1 'No name' "+file_handle+" | wc -l"], shell=True, stdout=subprocess.PIPE)
	no_name = p_grep.stdout.read()

	fasta = open(file_handle,"r")
	fix = open(file_handle[0:-6]+".fx"+file_handle[-6:len(file_handle)], "w")
	i = 1
	for line in fasta:
		new_line = ""
		if line[0] == '>':
			if int(no_name) == 1:
				new_line = ">Unknown_"+`i`+"\n"
				i+=1
			else:
				new_line = line
			new_line.replace(":","_")
			fix.write(new_line)
		else:
			fix.write(line)
	fasta.close()
	fix.close()
else:
	sys.exit("ERROR: Invalid file "+ file_handle)
