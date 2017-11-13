#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdio>
#include <omp.h>
#include "config.h"


typedef std::bitset<3*readlen> bitset;

uint32_t numreads = 0;
std::string mode;
double param1, param2, param3, param4, param5;

std::string infilenumreads;

std::string infile;
std::string infile_quality;
std::string outfile;
std::string outfile_quality;
std::string outdir;

struct overlap 
{
	char *read;
	int shift;
	uint32_t hamming;
	double weight;
	char *quality;
	uint32_t rid;
	bool RC;
};	 


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

void denoise_read(char *current_read, char *current_quality, char *denoised_read, char *denoised_quality ,std::vector<overlap> &overlap_vec);

void find_overlapping_reads(uint32_t &current_rid, bitset &current_bitset, char *current_quality, std::vector<overlap> &overlap_vec, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict, bool *added_rids);

void find_overlapping_reads_at_shift(int shift, bool rev, bitset &current_bitset, char *current_quality, std::vector<overlap> &overlap_vec, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict, uint32_t &current_rid, bool *added_rids);

void readDnaFile(bitset *read);

void readQualityFile(char(*quality)[readlen+1]);

void constructdictionary(bitset *read, std::unordered_map<uint64_t,uint32_t*> *dict);

void findMajority(std::vector<std::array<double,5>> &count, char *denoised_read);

int string_hamming(const char *a, const char *b, int startposa, int startposb, int matchlen);
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
	param1 = atof(argv[3]); //max hamming thresh
	param2 = atof(argv[4]); //denoising param
	param3 = atof(argv[5]); //alpha
	param4 = atof(argv[6]); //max number thresh
	param5 = atof(argv[7]); //max number thresh
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
	std::vector<overlap> overlap_vec;
	std::ofstream fout_overlap(outfile+".overlap."+std::to_string(tid));
	std::ifstream fin_clean("/srv/data/shubham/assembly/fastq_denoising/assembly_data/c_elegans_c25_100_180.upper.clean");
	fin_clean.seekg(uint64_t(i)*(readlen+1), fin_clean.beg);
	char clean_read[readlen+1];

	bitset current_bitset;
	char current_quality[readlen+1], current_read[readlen+1], denoised_read[readlen+1], denoised_quality[readlen+1];
	bool *added_rids = new bool [numreads];
	std::fill(added_rids, added_rids+numreads,false);
	bool exact_match_found;
	while(i < stop)
	{
		fin_clean.getline(clean_read,readlen+1);
		current_bitset = read[i];

		strcpy(current_quality,quality[i]);

		find_overlapping_reads(i, current_bitset,current_quality,overlap_vec,read,quality,dict,added_rids);
		strcpy(current_read, read_char[i]);

		strcpy(denoised_quality,current_quality);
		strcpy(denoised_read,current_read);
		denoise_read(current_read, current_quality,denoised_read,denoised_quality,overlap_vec);
		if(strcmp(clean_read,denoised_read) != 0 && overlap_vec.size() < 100)
		{
			////////// DEBUG mode
			fout_overlap << std::string(100, '#') << "\n";
			std::sort(overlap_vec.begin(), overlap_vec.end(), [](overlap read1, overlap read2) {
				return read1.shift < read2.shift;
			});
			fout_overlap << clean_read << "\n";
			fout_overlap << current_read << "\n";
			fout_overlap << denoised_read << "\n";
			fout_overlap << std::string(100, '#') << "\n";
			fout_overlap << clean_read << "\n";
			for( int j = 0; j < readlen; j++)
			{
				if(current_read[j] != clean_read[j])
					fout_overlap << current_read[j];
				else
					fout_overlap << "-";
			}
			fout_overlap << "\n";
			for( int j = 0; j < readlen; j++)
			{
				if(denoised_read[j] != clean_read[j])
					fout_overlap << denoised_read[j];
				else
					fout_overlap << "-";
			}
			fout_overlap << "\n";
				
			for(std::vector<overlap>::iterator it = overlap_vec.begin() ; it != overlap_vec.end(); ++it)
			{
				fout_overlap << std::string(maxshift+(*it).shift, '*') << (*it).read << "\n";
			}
		}
		fout << denoised_read << "\n";
		fout_quality << denoised_quality << "\n";
		i++;
		overlap_vec.clear();
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


void denoise_read(char *current_read, char *current_quality, char *denoised_read, char *denoised_quality ,std::vector<overlap> &overlap_vec)
{       
	int startposa, startposb, matchlen;//variables for string_hamming
	//find hamming distances of all reads to current read
	
	for(auto it = overlap_vec.begin(); it != overlap_vec.end(); ++it)
	{
		if((*it).shift <= 0)
		{
			startposa = 0;
			matchlen = readlen + (*it).shift;
			startposb = readlen -  matchlen;
		}
		else
		{
			startposb = 0;
			matchlen = readlen - (*it).shift;
			startposa = readlen -  matchlen;
		}
		(*it).hamming = string_hamming(current_read,(*it).read,startposa,startposb,matchlen);
	}
	//sort overlap_vec in increasing hamming distance. For same hamming distance, higher priority to less shift
	std::sort(overlap_vec.begin(), overlap_vec.end(), [](overlap read1, overlap read2) {
					return (read1.hamming==read2.hamming)?(std::abs(read1.shift)<std::abs(read2.shift)):(read1.hamming<read2.hamming);
	});	
	std::vector<std::array<double,5>> read_count(readlen,{0,0,0,0,0});
	std::vector<double> total_count(readlen,0);
	int num_bases_denoised = 0;
	std::vector<int> denoised_bases(readlen,0);
	for(auto it = overlap_vec.begin(); it != overlap_vec.end(); ++it)
	{
		if((*it).shift <= 0)
		{
			startposa = 0;
			matchlen = readlen + (*it).shift;
			startposb = readlen -  matchlen;
		}
		else
		{
			startposb = 0;
			matchlen = readlen - (*it).shift;
			startposa = readlen -  matchlen;
		}
		for(int i = 0; i < matchlen; i++)
		{
			if(denoised_bases[startposa+i])
				continue;
			read_count[startposa+i][chartolong[((*it).read)[startposb+i]]] += 1;
			total_count[startposa+i]++;
			if((read_count[startposa+i][chartolong[((*it).read)[startposb+i]]]+0.25)/(total_count[startposa+i]+1) >= param4 && total_count[startposa+i]>param5)
			{
				denoised_read[startposa+i] = ((*it).read)[startposb+i];	
				denoised_bases[startposa+i] = 1;
				num_bases_denoised++;	
			}
		}
		if(num_bases_denoised == readlen)
			break;
	}
	return;
}

void find_overlapping_reads(uint32_t &current_rid, bitset &current_bitset, char *current_quality, std::vector<overlap> &overlap_vec, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict, bool *added_rids)
{
	bitset forward_bitset = current_bitset;
	bitset reverse_bitset = chartobitset(read_char_RC[current_rid]);

	//find matches at shift 0
	find_overlapping_reads_at_shift(0, false, forward_bitset, current_quality, overlap_vec, read, quality, dict, current_rid, added_rids);
	find_overlapping_reads_at_shift(0, true, reverse_bitset, current_quality, overlap_vec, read, quality, dict, current_rid, added_rids);

	for(auto it = overlap_vec.begin(); it != overlap_vec.end(); ++it)
		added_rids[(*it).rid] = false;
        bitset forward_bitset_left = current_bitset;
        bitset reverse_bitset_left = chartobitset(read_char_RC[current_rid]);
        //find matches to right/left
        for(int shift = 1; shift < maxshift; shift++)
        {
                forward_bitset >>= 3;
                reverse_bitset <<= 3;
                forward_bitset_left <<= 3;
                reverse_bitset_left >>= 3;
                find_overlapping_reads_at_shift(shift, false, forward_bitset, current_quality, overlap_vec, read, quality, dict, current_rid, added_rids);
                find_overlapping_reads_at_shift(shift, true, reverse_bitset, current_quality, overlap_vec, read, quality, dict, current_rid, added_rids);
                find_overlapping_reads_at_shift(-shift, false, forward_bitset_left, current_quality, overlap_vec, read, quality, dict, current_rid, added_rids);
                find_overlapping_reads_at_shift(-shift, true, reverse_bitset_left, current_quality, overlap_vec, read, quality, dict, current_rid, added_rids);
	for(auto it = overlap_vec.begin(); it != overlap_vec.end(); ++it)
		added_rids[(*it).rid] = false;
        }

	//clear added_rids array now
	for(auto it = overlap_vec.begin(); it != overlap_vec.end(); ++it)
		added_rids[(*it).rid] = false;

	return;
}

void find_overlapping_reads_at_shift(int shift, bool rev, bitset &current_bitset, char *current_quality, std::vector<overlap> &overlap_vec, bitset *read, char(*quality)[readlen+1], std::unordered_map<uint64_t,uint32_t*> *dict, uint32_t &current_rid, bool *added_rids)
{
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	uint64_t ull;//for indexing dict
	//if shift is 0 and rev is false, add the current read
	if(shift == 0 && rev == false)
	{
		added_rids[current_rid] = true;
		overlap_vec.push_back({read_char[current_rid],0,0,0,quality[current_rid],current_rid,false});
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
				auto read_hamming_distance = (current_bitset^(read[rid]&readmask)).count();
				double read_weight = param3*read_hamming_distance + (1.0-param3)*std::abs(shift);
				if(read_hamming_distance <= param1)
				{	
					//std::cout << read_weight << std::endl;
					if(added_rids[rid] == false)//not already present
					{
						added_rids[rid] = true;
						if(!rev)	
						{
							overlap_vec.push_back({read_char[rid],shift,0,read_weight,quality[rid],rid,false});
						}
						else
						{
							overlap_vec.push_back({read_char_RC[rid],shift,0,read_weight,quality[rid],rid,true});
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
int string_hamming(const char *a, const char *b, int startposa, int startposb, int matchlen)
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

