#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdio>
#include <omp.h>
#include "config.h"


typedef std::bitset<3*readlen> bitset;

uint32_t numreads = 0;
std::string mode;
double param1, param2, param3;

std::string infilenumreads;

std::string infile;
std::string infile_quality;
std::string outfile;
std::string outfile_quality;
std::string outdir;


//Some global arrays (some initialized in setglobalarrays())
char revinttochar[8] = {'A','N','G','#','C','#','T','#'};//used in bitsettostring
char inttochar[] = {'A','C','G','T','N'};
char chartorevchar[128];//A-T etc for reverse complement
bitset basemask[readlen][128];//bitset for A,G,C,T,N at each position 
//used in stringtobitset, chartobitset and bitsettostring
bitset positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
bitset mask63;//bitset with 63 bits set to 1 (used in bitsettostring for conversion to ullong)
char longtochar[] = {'A','C','G','T','N'};
long chartolong[128];


int *dict_start;
int *dict_end; 
bitset mask[maxshift];
bitset revmask[maxshift];
bitset indexmask[numdict];


double quality_to_p_pbar[42];
double log_1_p[42];
double log_p_3[42];
double one_minus_p[42];

char(*read_char)[readlen+1];
char(*read_char_RC)[readlen+1];


void denoise(bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict);

void denoise_read(char *current_read, char *current_quality, char *denoised_read, char *denoised_quality ,std::vector<char*> &overlap_reads, std::vector<int> &overlap_shift, std::vector<char*> &overlap_quality, std::vector<bool> &overlap_RC);

void find_overlapping_reads(uint32_t &current_rid, bitset &current_bitset, char *current_quality, std::vector<char*> &overlap_reads, std::vector<int> &overlap_shift, std::vector<char*> &overlap_quality, std::vector<bool> &overlap_RC, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict);

void find_overlapping_reads_at_shift(int shift, bool rev, bitset &current_bitset, char *current_quality, std::vector<char*> &overlap_reads, std::vector<int> &overlap_shift, std::vector<char*> &overlap_quality, std::vector<bool> &overlap_RC, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict, uint32_t &current_rid);

void readDnaFile(bitset *read);

void readQualityFile(char(*quality)[readlen+1]);

void constructdictionary(bitset *read, std::unordered_map<uint64_t,uint32_t*> *dict);

void findMajority(std::vector<std::array<double,5>> &count, char *denoised_read);

int string_hamming(char *a, char *b, int &startposa, int &startposb, int &matchlen);
//compute hamming distance between a[startposa:startposa+matchlen] and b[startposb:startposb+matchlen]

bitset chartobitset(char *s);

std::string findMajority(std::vector<std::array<double,5>> &count);



bitset stringtobitset(std::string s);

void generateindexmasks(bitset *indexmask);

void generatemasks(bitset *mask,bitset *revmask);

void bitsettostring(bitset b,char *s);

void reverse_complement(char* s, char* s1);

std::string bitsettostring(bitset b);

std::string reverse_complement(std::string s);

void setglobalarrays();


int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outdir = basedir + "/output/";
	infile = basedir + "/output/input.dna";
	infile_quality = basedir + "/output/input.quality";
	outfile = basedir + "/output/output.dna";
	outfile_quality = basedir + "/output/output.quality";
	infilenumreads = basedir + "/output/numreads.bin";
	
	mode = std::string(argv[2]);
	param1 = atof(argv[3]);
	param2 = atof(argv[4]);
	param3 = atof(argv[5]);
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));	
	omp_set_num_threads(num_thr);	
	setglobalarrays();

	std::cout << "Reading file: " << infile << std::endl;

	read_char = new char [numreads][readlen+1];
	read_char_RC = new char [numreads][readlen+1];
	bitset *read = new bitset [numreads];
	readDnaFile(read);

	char(*quality)[readlen+1] = new char [numreads][readlen+1];
	readQualityFile(quality);

	std::unordered_map<uint64_t,uint32_t*> dict[numdict];
	std::cout << "Constructing dictionaries\n";
	constructdictionary(read,dict);

	std::cout << "Denoising reads\n";
	denoise(read,quality,dict);

	delete[] read;
	delete[] quality;
	delete[] read_char;
 	delete[] read_char_RC;
	std::cout << "Done!\n";
	return 0;
}



void denoise(bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	std::ofstream fout(outfile+'.'+std::to_string(tid));
	std::ofstream fout_quality(outfile_quality+'.'+std::to_string(tid));
	std::vector<char*> overlap_reads;//reads overlapping with current read
	std::vector<int> overlap_shift;//shift of overlapping reads wrt current read - 0 for current read, -ve for reads on left
	std::vector<char*> overlap_quality;//qualities of overlapping reads (reversed if read was reversed)
	std::vector<bool> overlap_RC;//reverse_complemented wrt current read (true if so), used for handling quality

	bitset current_bitset;
	char current_quality[readlen+1], current_read[readlen+1], denoised_read[readlen+1], denoised_quality[readlen+1];
	while(i < stop)
	{
		current_bitset = read[i];

		strcpy(current_quality,quality[i]);

		find_overlapping_reads(i, current_bitset,current_quality,overlap_reads,overlap_shift,overlap_quality,overlap_RC,read,quality,dict);
		strcpy(current_read, read_char[i]);

		strcpy(denoised_quality,current_quality);
		strcpy(denoised_read,current_read);
		
		denoise_read(current_read, current_quality,denoised_read,denoised_quality,overlap_reads,overlap_shift,overlap_quality,overlap_RC);
		fout << denoised_read << "\n";
		fout_quality << denoised_quality << "\n";
		i++;
		overlap_reads.clear();
		overlap_shift.clear();	
		overlap_quality.clear();
		overlap_RC.clear();	
	}
	fout.close();
	fout_quality.close();
	}//parallel end

	//combine files produced by the threads
	std::ofstream fout(outfile);
	std::ofstream fout_quality(outfile_quality);

	for(int tid = 0; tid < num_thr; tid++)
	{
		std::ifstream fin(outfile+'.'+std::to_string(tid));
		std::ifstream fin_quality(outfile_quality+'.'+std::to_string(tid));

		fout << fin.rdbuf();
		fout_quality << fin_quality.rdbuf();
		
		fout.clear();
		fout_quality.clear();

		remove((outfile+'.'+std::to_string(tid)).c_str());
		remove((outfile_quality+'.'+std::to_string(tid)).c_str());
	}
	return;
}

void denoise_read(char *current_read, char *current_quality, char *denoised_read, char *denoised_quality ,std::vector<char*> &overlap_reads, std::vector<int> &overlap_shift, std::vector<char*> &overlap_quality, std::vector<bool> &overlap_RC)
{        
	std::vector<std::array<double,5>> read_count(readlen,{0,0,0,0,0});
	auto reads_it = overlap_reads.begin();
	auto shift_it = overlap_shift.begin();
	auto quality_it = overlap_quality.begin();
	auto RC_it = overlap_RC.begin();
	double read_weight;
	int startposa, startposb, matchlen;//variables for string_hamming
	int quality_pos; // used for handling reverse complementation
	// accumulate counts
	while(reads_it != overlap_reads.end())
	{
		if(*shift_it <= 0)
		{
		    startposa = 0;
		    matchlen = readlen + *shift_it;
		    startposb = readlen -  matchlen;
		}
		else
		{
		    startposb = 0;
		    matchlen = readlen - *shift_it;
		    startposa = readlen -  matchlen;
		}
		read_weight = 1 - sqrt((double)string_hamming(current_read,*reads_it,startposa,startposb,matchlen)/(matchlen*0.25));
		//	read_weight = 1;
		if(read_weight < 0)
			read_weight = 0;
		for(int i = 0; i < matchlen; i++)
		{
			if(*RC_it == true)
				quality_pos = readlen - (startposb+i) - 1;
			else
				quality_pos = (startposb+i);
			read_count[startposa+i][chartolong[(*reads_it)[startposb+i]]] += read_weight*one_minus_p[(uint32_t)((*quality_it)[quality_pos])-33];
		}
		++shift_it;
		++reads_it;
		++quality_it;
		++RC_it;
	}
	findMajority(read_count, denoised_read);
	for(int j = 0; j < readlen; j++)
	{
		if(read_count[j][chartolong[current_read[j]]] > param2)
			denoised_read[j] = current_read[j];
	}
	return;
}

void find_overlapping_reads(uint32_t &current_rid, bitset &current_bitset, char *current_quality, std::vector<char*> &overlap_reads, std::vector<int> &overlap_shift, std::vector<char*> &overlap_quality, std::vector<bool> &overlap_RC, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict)
{
	bitset forward_bitset = current_bitset;
	bitset reverse_bitset = chartobitset(read_char_RC[current_rid]);

	//find matches at shift 0
	find_overlapping_reads_at_shift(0, false, forward_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, overlap_RC, read, quality, dict, current_rid);
	find_overlapping_reads_at_shift(0, true, reverse_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, overlap_RC, read, quality, dict, current_rid);

	//find matches to right
	for(int shift = 1; shift < maxshift; shift++)
	{	
		forward_bitset >>= 3;
		reverse_bitset <<= 3;
		find_overlapping_reads_at_shift(shift, false, forward_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, overlap_RC, read, quality, dict, current_rid);
		find_overlapping_reads_at_shift(shift, true, reverse_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, overlap_RC, read, quality, dict, current_rid);
	}

	//find matches to left
	forward_bitset = current_bitset;
	reverse_bitset = chartobitset(read_char_RC[current_rid]);
	for(int shift = -1; shift > -maxshift; shift--)
	{	
		forward_bitset <<= 3;
		reverse_bitset >>= 3;
		find_overlapping_reads_at_shift(shift, false, forward_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, overlap_RC, read, quality, dict, current_rid);
		find_overlapping_reads_at_shift(shift, true, reverse_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, overlap_RC, read, quality, dict, current_rid);
	}
	return;
}

void find_overlapping_reads_at_shift(int shift, bool rev, bitset &current_bitset, char *current_quality, std::vector<char*> &overlap_reads, std::vector<int> &overlap_shift, std::vector<char*> &overlap_quality, std::vector<bool> &overlap_RC, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict, uint32_t &current_rid)
{
	std::unordered_set<uint32_t> added_rids;//to make sure rids are not inserted multiple times into the vectors
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	uint64_t ull;//for indexing dict

	//if shift is 0 and rev is false, add the current read
	if(shift == 0 && rev == false)
	{
		added_rids.insert(current_rid);
		overlap_reads.push_back(read_char[current_rid]);
		overlap_quality.push_back(quality[current_rid]);
		overlap_shift.push_back(0);	
		overlap_RC.push_back(false);
	}

	bitset readmask;
	if(rev == false && shift >= 0)
		readmask = mask[shift];
	else if(rev == true && shift >= 0)
		readmask = revmask[shift];
	else if(rev == false && shift < 0)
		readmask = revmask[-shift];
	else if(rev == true && shift < 0)
		readmask = mask[-shift];

	for(int l = 0; l < numdict; l++)
	{
		//check if dictionary endpoints fit in the shifted read
		if(shift >= 0)
		{	
			if(!rev)
			{
				if(dict_end[l]+shift >= readlen)
					continue;
			}
			else
			{
				if(dict_start[l] < shift)
					continue;
			}
		}
		else
		{	
			if(!rev)
			{
				if(dict_start[l] - (-shift) < 0)
					continue;
			}
			else
			{
				if(dict_end[l] + (-shift) >= readlen)
					continue;
			}
		}
		bitset b = current_bitset&indexmask[l];
		ull = (b>>3*dict_start[l]).to_ullong();
		if(dict[l].count(ull) == 1)
		{
		auto s = dict[l][ull];
			for (int64_t i = s[0] - 1; i >=1 && i >= long(s[0]) - maxsearch ; i--)
			{
				auto rid = s[i];
				if((current_bitset^(read[rid]&readmask)).count()<=param1)
				{	
					if(added_rids.find(rid) == added_rids.end())//not already present
					{
						added_rids.insert(rid);
						overlap_shift.push_back(shift);	
						if(!rev)	
						{
							overlap_reads.push_back(read_char[rid]);
							overlap_quality.push_back(quality[rid]);
							overlap_RC.push_back(false);
						}
						else
						{
							overlap_reads.push_back(read_char_RC[rid]);
							overlap_quality.push_back(quality[rid]);
							overlap_RC.push_back(true);
						}	
					}
				}
			}
		}
	}
	return;
}

void readDnaFile(bitset *read)
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	std::ifstream f(infile, std::ifstream::in);
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	std::string s;
	while(i < stop)
	{
		f.getline(read_char[i],readlen+1);
		read[i] = chartobitset(read_char[i]);
		reverse_complement(read_char[i],read_char_RC[i]);
		read_char_RC[i][readlen] = '\0';		
		i++;
	}
	f.close();
	}
	return;
}


void readQualityFile(char(*quality)[readlen+1])
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;
	std::ifstream f(infile_quality, std::ifstream::in);
	f.seekg(uint64_t(i)*(readlen+1), f.beg);
	std::string s;
	while(i < stop)
	{
		f.getline(quality[i],readlen+1);
		i++;
	}
	f.close();
	}
	return;
}

void constructdictionary(bitset *read, std::unordered_map<uint64_t,uint32_t*> *dict)
{

	//Parallelizing construction of the multiple dictionaries
	#pragma omp parallel num_threads(numdict)
	{
	#pragma omp for
	for(int j = 0; j < numdict; j++)
	{	
		bitset b;
		uint64_t ull;
		//find number of times each key occurs
		for(uint32_t i = 0; i < numreads; i++)
		{
			b = read[i]&indexmask[j];
			ull = (b>>3*dict_start[j]).to_ullong();
			if(dict[j].count(ull) == 1)
				(*dict[j][ull])++;
			else
			{
				dict[j][ull]=new uint32_t;
				(*dict[j][ull]) = 1;
			}
		}
		//allocate memory for each bin (number of reads with given key + 1) 1 for storing the length
		for(auto it = dict[j].begin(); it !=  dict[j].end(); ++it)
		{
			uint32_t binsize = *(it->second);
			delete it->second;
			dict[j][it->first] =  new uint32_t[binsize+1];
			dict[j][it->first][0] = 1;
		}
		//fill in the read ids in each bin, dict[j][ull][0] stores the position where next id is put - at the
		//end it stores the size of the array
		for(uint32_t i = 0; i < numreads; i++)
		{
			b = read[i]&indexmask[j];
			ull = (b>>3*dict_start[j]).to_ullong();
			dict[j][ull][dict[j][ull][0]++] = i;
		}

	}
	}
	return;
}


std::string findMajority(std::vector<std::array<double,5>> &count)
{
    std::string ref(count.size(),'A');
    for(long i = 0; i < count.size(); i++)
    {
        double max = 0;
        long indmax = 0;
        for(long j = 0; j < 5; j++)
            if(count[i][j]>max)
            {
                max = count[i][j];
                indmax = j;
            }
        ref[i] = longtochar[indmax];
    }
    return ref;
}

void findMajority(std::vector<std::array<double,5>> &count, char *denoised_read)
{
    for(long i = 0; i < count.size(); i++)
    {
        double max = 0;
        long indmax = 0;
        for(long j = 0; j < 5; j++)
            if(count[i][j]>max)
            {
                max = count[i][j];
                indmax = j;
            }
        denoised_read[i] = longtochar[indmax];
    }
    return;
}
/*
int string_hamming(std::string &a, std::string &b, int &startposa, int &startposb, int &matchlen)
{
    int dist = 0;
    for(int k = 0; k < matchlen; k++)
        if(a[startposa+k] != b[startposb+k])
            dist++;
    return dist;
}*/
int string_hamming(char *a, char *b, int &startposa, int &startposb, int &matchlen)
{
    int dist = 0;
    for(int k = 0; k < matchlen; k++)
        if(a[startposa+k] != b[startposb+k])
            dist++;
    return dist;
}

bitset chartobitset(char *s)
{
	bitset b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

bitset stringtobitset(std::string s)
{
	bitset b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void generateindexmasks(bitset *indexmask)
//masks for dictionary positions
{
	for(int j = 0; j < numdict; j++)
		indexmask[j].reset();
	for(int j = 0; j < numdict; j++)
		for(int i = 3*dict_start[j]; i < 3*(dict_end[j]+1); i++)
			indexmask[j][i] = 1;
	return;
}

void generatemasks(bitset *mask,bitset *revmask)
{
	for(int i = 0; i < maxshift; i++)
	{	
		mask[i].reset();
		revmask[i].reset();
		for(int j = 0; j < 3*readlen - 3*i; j++)
			mask[i][j] = 1;
		for(int j = 3*i; j < 3*readlen; j++)
			revmask[i][j] = 1; 	
	}
	return;
}

void bitsettostring(bitset b,char *s)
{
	unsigned long long ull,rem;
	for(int i = 0; i < 3*readlen/63+1; i++)
	{	
		ull = (b&mask63).to_ullong();
		b>>=63;
		for(int j = 21*i  ; j < 21*i+21 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%8];	
			ull/=8;
		}
	}
	s[readlen] = '\0';
	return;
}


std::string bitsettostring(bitset b)
{
	char s[readlen+1];
	s[readlen] = '\0';
	unsigned long long ull,rem;
	for(int i = 0; i < 3*readlen/63+1; i++)
	{	
		ull = (b&mask63).to_ullong();
		b>>=63;
		for(int j = 21*i  ; j < 21*i+21 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%8];	
			ull/=8;
		}
	}
	std::string s1 = s;
	return s1;
}


std::string reverse_complement(std::string s)
{
	std::string s1(s);
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return s1;
}

void reverse_complement(char* s, char* s1)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return;
}

void setglobalarrays()
{
	chartorevchar['A'] = 'T';
	chartorevchar['C'] = 'G';
	chartorevchar['G'] = 'C';
	chartorevchar['T'] = 'A';
	chartorevchar['N'] = 'N';

	chartolong['A'] = 0;
	chartolong['C'] = 1;
	chartolong['G'] = 2;
	chartolong['T'] = 3;
	chartolong['N'] = 4;
	#if numdict == 1
	{
		dict_start = new int[1];
		dict_end = new int[1];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
	}
	#elif numdict == 2
	{
		dict_start = new int[2];
		dict_end = new int[2];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
	}
	#elif numdict == 3
	{
		dict_start = new int[3];
		dict_end = new int[3];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
		dict_start[2] = dict3_start;
		dict_end[2] = dict3_end;
	}
	#elif numdict == 4
	{
		dict_start = new int[4];
		dict_end = new int[4];
		dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;
		dict_start[2] = dict3_start;
		dict_end[2] = dict3_end;
		dict_start[3] = dict4_start;
		dict_end[3] = dict4_end;
	}
	#endif
	for(int i = 0; i < 63; i++)
		mask63[i] = 1;
	for(int i = 0; i < readlen; i++)
	{
		basemask[i]['A'][3*i] = 0;
		basemask[i]['A'][3*i+1] = 0;
		basemask[i]['A'][3*i+2] = 0;
		basemask[i]['C'][3*i] = 0;
		basemask[i]['C'][3*i+1] = 0;
		basemask[i]['C'][3*i+2] = 1;
		basemask[i]['G'][3*i] = 0;
		basemask[i]['G'][3*i+1] = 1;
		basemask[i]['G'][3*i+2] = 0;
		basemask[i]['T'][3*i] = 0;
		basemask[i]['T'][3*i+1] = 1;
		basemask[i]['T'][3*i+2] = 1;
		basemask[i]['N'][3*i] = 1;
		basemask[i]['N'][3*i+1] = 0;
		basemask[i]['N'][3*i+2] = 0;
		positionmask[i][3*i] = 1;
		positionmask[i][3*i+1] = 1;
		positionmask[i][3*i+2] = 1;
	}		
	generatemasks(mask,revmask);
	generateindexmasks(indexmask);

	for(int i=0; i<42;i++) 
	{	
		double _prob = pow(10.0,-((double)i/10.0));
		quality_to_p_pbar[i] = _prob/(1.0-_prob);
		log_1_p[i] = log10(1-_prob);
		log_p_3[i] = log10(_prob/3);
		one_minus_p[i] = 1-_prob;
	}
	return;
}

