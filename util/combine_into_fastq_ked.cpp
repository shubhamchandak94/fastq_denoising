#include <fstream>
#include <iostream>
#include <string>

std::string basename;
std::string suffix;
std::string infile_reads,infile_quality,infile_id,outfile_fastq_1,outfile_fastq_2;

int main(int argc, char** argv)
{
	basename = std::string(argv[1]);
	infile_reads = basename + ".dna.d";
	infile_quality = basename + ".quality";
	infile_id = basename + ".id";
	outfile_fastq_1 = basename + "_1.harc.fastq";
	outfile_fastq_2 = basename + "_2.harc.fastq";

	std::ifstream fin_reads(infile_reads);
	std::ifstream fin_id(infile_id);
	std::ifstream fin_quality(infile_quality);
	std::ofstream fout_fastq_1(outfile_fastq_1);
	std::ofstream fout_fastq_2(outfile_fastq_2);

	std::string read,id,quality;
	int number_of_reads = 0;
	while(std::getline(fin_reads, read))
        	++number_of_reads;
	int num_reads_per_file = number_of_reads/2;
	fin_reads.clear();
	fin_reads.seekg(0, std::ios::beg);

	for(int i=0; i < num_reads_per_file; i++) 
	{
		std::getline(fin_reads, read);	
		std::getline(fin_id,id);
		std::getline(fin_quality,quality);
		fout_fastq_1 << id << "\n";
		fout_fastq_1 << read << "\n";	
		fout_fastq_1 << "+\n";
		fout_fastq_1 << quality << "\n";
	}
	for(int i=0; i < num_reads_per_file; i++) 
	{
		std::getline(fin_reads, read);	
		std::getline(fin_id,id);
		std::getline(fin_quality,quality);
		fout_fastq_2 << id << "\n";
		fout_fastq_2 << read << "\n";	
		fout_fastq_2 << "+\n";
		fout_fastq_2 << quality << "\n";
	}

	fin_reads.close();
	fin_quality.close();
	fin_id.close();
	fout_fastq_1.close();
	fout_fastq_2.close();
	return 0;
}
