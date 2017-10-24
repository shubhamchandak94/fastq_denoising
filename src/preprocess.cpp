#include <iostream>
#include <fstream>
#include <string>

std::string infile;
std::string outfile;
std::string outfilequality;
std::string outfileid;
std::string outfilenumreads;

int readlen;


int preprocess();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[2]);
	infile = std::string(argv[1]);
	readlen = atoi(argv[3]);
	outfile = basedir + "/output/input.dna";
	outfilequality = basedir + "/output/input.quality";
	outfileid = basedir + "/output/output.id";
	outfilenumreads = basedir + "/output/numreads.bin";
	int status = preprocess();
	if(status != 0)
		return -1;
	std::cout << "Preprocessing Done!\n";
	return 0;
}

int preprocess()
{
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	std::ofstream f(outfile);
	std::ofstream f_quality(outfilequality);
	std::ofstream f_id(outfileid);
	
	int i = 0;
	uint64_t readnum = 0;
	while(std::getline(myfile, line))
	{
		switch(i)
		{
			case 0:	f_id << line << "\n";
				break;
			case 1: 
				if(line.length() != readlen)
				{	
					std::cout << "Read length not fixed. Found two different read lengths: "<< 
						  readlen << " and " << line.length() << "\n";
					return -1;
				}
				f << line << "\n";
				break;
			case 2: break;
			case 3: f_quality << line << "\n";
					readnum++;
					break;
		}
		i = (i+1)%4;
	}
	if(readnum > 4294967290)
	{
		std::cout << "Too many reads. HARC supports at most 4294967290 reads\n";
		return -1;
	}
	else
	{
		std::ofstream f_numreads(outfilenumreads,std::ios::binary);
		uint32_t readnum_32 = readnum;
		f_numreads.write((char*)&readnum_32, sizeof(uint32_t));
		std::cout << "Read length: " << readlen << "\n";
		std::cout << "Total number of reads: " << readnum <<"\n";
		f_numreads.close();
	}	
	return 0;	
}
