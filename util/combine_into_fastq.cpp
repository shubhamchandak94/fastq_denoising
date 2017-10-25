#include <fstream>
#include <iostream>
#include <string>

std::string basename;
std::string infile_reads,infile_quality,infile_id,outfile_fastq;

int main(int argc, char** argv)
{
	basename = std::string(argv[1]);
	infile_reads = basename + ".dna.d";
	infile_quality = basename + ".quality";
	infile_id = basename + ".id";
	outfile_fastq = basename + ".output.fastq";

	std::ifstream fin_reads(infile_reads);
	std::ifstream fin_id(infile_id);
	std::ifstream fin_quality(infile_quality);
	std::ofstream fout_fastq(outfile_fastq);


	std::string read,id,quality;
	while(std::getline(fin_reads,read))
	{
		std::getline(fin_id,id);
		std::getline(fin_quality,quality);
		fout_fastq << id << "\n";
		fout_fastq << read << "\n";	
		fout_fastq << "+\n";
		fout_fastq << quality << "\n";
	}
	fin_reads.close();
	fin_quality.close();
	fin_id.close();
	fout_fastq.close();
	return 0;
}
