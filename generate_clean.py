f_1 = "assembly_data/SRR065390.dna"
f_2 = "assembly_data/SRR065390_align.txt"
f_3 = "tempfile"
f_out = "assembly_data/SRR065390.clean"
numreads = 67617092
readlen = 100

file_1 = open(f_1,'r')
file_2 = open(f_2,'r')
file_3 = open(f_3,'r')
file_out = open(f_out,'w')

i = 0
read_align = file_2.readline()
read_align = read_align[:readlen]
for i in range(numreads):
	read_original = file_1.readline().rstrip('\n')
	if(read_original == read_align):
		read_clean = file_3.readline().rstrip('\n')
		if(read_clean == '-'*readlen):
			file_out.write(read_align+'\n')
		else:
			file_out.write(read_clean+'\n')	
		read_align = file_2.readline()
		read_align = read_align[:readlen]		
	else:
		file_out.write('#'*readlen+'\n')
