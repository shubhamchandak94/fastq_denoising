#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <array>
#include <unordered_set>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdio>
#include <omp.h>
#include "BooPHF.h"
#include "config.h"


typedef std::bitset<3*readlen> bitset;
typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t> boophf_t;

class bbhashdict
{
	public:
	boophf_t * bphf;
	uint32_t numkeys;
	uint32_t *startpos;
	uint32_t *read_id;
	bool *empty_bin;
	void findpos(int64_t *dictidx, uint32_t &startposidx);
	void remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current);
	bbhashdict()
	{
		bphf = NULL;
		startpos = NULL;
		read_id = NULL;
	}
	~bbhashdict()
	{
		delete[] startpos;
		delete[] read_id;
		delete bphf;
	}	
};

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


void denoise(bitset *read, char(*quality)[readlen+1], bbhashdict *dict);

void find_overlapping_reads(uint32_t &current_rid, bitset &current_bitset, std::string &current_quality, std::vector<std::string> &overlap_reads, std::vector<int> &overlap_shift, 
							std::vector<std::string> &overlap_quality, bitset *read, char(*quality)[readlen+1], bbhashdict *dict);

void find_overlapping_reads_at_shift(int shift, bool rev, bitset &current_bitset, std::string &current_quality, std::vector<std::string> &overlap_reads, std::vector<int> &overlap_shift, 
							std::vector<std::string> &overlap_quality, bitset *read, char(*quality)[readlen+1], bbhashdict *dict, uint32_t &current_rid);

void find_overlapping_reads(uint32_t &rid, bitset &current_bitset, std::string &current_quality, std::vector<std::string> &overlap_reads, std::vector<int> &overlap_shift, 
							std::vector<std::string> &overlap_quality, bitset *read, char(*quality)[readlen+1], bbhashdict *dict);


void denoise_read(std::string &current_read, std::string &current_quality, std::string &denoised_read, std::string &denoised_quality ,std::vector<std::string> &overlap_reads, 
		std::vector<int> &overlap_shift, std::vector<std::string> &overlap_quality);

void readDnaFile(bitset *read);

void readQualityFile(char(*quality)[readlen+1]);

void constructdictionary(bitset *read, bbhashdict *dict);

std::string findMajority(std::vector<std::array<double,5>> &count);

int string_hamming(std::string &a, std::string &b, int &startposa, int &startposb, int &matchlen);
//compute hamming distance between a[startposa:startposa+matchlen] and b[startposb:startposb+matchlen]

bitset stringtobitset(std::string s);

void generateindexmasks(bitset *indexmask);

void generatemasks(bitset *mask,bitset *revmask);

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

	bitset *read = new bitset [numreads];
	readDnaFile(read);

	char(*quality)[readlen+1] = new char [numreads][readlen+1];
	readQualityFile(quality);

	bbhashdict dict[numdict];
	std::cout << "Constructing dictionaries\n";
	constructdictionary(read,dict);

	std::cout << "Denoising reads\n";
	denoise(read,quality,dict);

	delete[] read;
	delete[] quality;
	std::cout << "Done!\n";
	return 0;
}



void denoise(bitset *read, char(*quality)[readlen+1], bbhashdict *dict)
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
		std::vector<std::string> overlap_reads;//reads overlapping with current read
		std::vector<int> overlap_shift;//shift of overlapping reads wrt current read - 0 for current read, -ve for reads on left
		std::vector<std::string> overlap_quality;//qualities of overlapping reads (reversed if read was reversed)
	while(i < stop)
	{
		bitset current_bitset = read[i];
		std::string current_quality = quality[i];

		find_overlapping_reads(i, current_bitset,current_quality,overlap_reads,overlap_shift,overlap_quality,read,quality,dict);
		std::string current_read = bitsettostring(current_bitset);
		std::string denoised_read = current_read, denoised_quality = current_quality;
		denoise_read(current_read, current_quality,denoised_read,denoised_quality,overlap_reads,overlap_shift,overlap_quality);
		fout << denoised_read << "\n";
		fout_quality << denoised_quality << "\n";
		i++;
		overlap_reads.clear();
		overlap_shift.clear();	
		overlap_quality.clear();	
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

void denoise_read(std::string &current_read, std::string &current_quality, std::string &denoised_read, std::string &denoised_quality ,std::vector<std::string> &overlap_reads, 
		std::vector<int> &overlap_shift, std::vector<std::string> &overlap_quality)
{        
	std::vector<std::array<double,5>> read_count(readlen,{0,0,0,0,0});
	auto reads_it = overlap_reads.begin();
	auto shift_it = overlap_shift.begin();
	auto quality_it = overlap_quality.begin();
    double read_weight;
    int startposa, startposb, matchlen;//variables for string_hamming

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
                read_count[startposa+i][chartolong[(*reads_it)[startposb+i]]] += read_weight*one_minus_p[(uint32_t)((*quality_it)[startposb+i])-33];
        }
        ++shift_it;
        ++reads_it;
        ++quality_it;
    }
    denoised_read = findMajority(read_count);
    for(int j = 0; j < readlen; j++)
        if(read_count[j][chartolong[current_read[j]]] > param2)
                denoised_read[j] = current_read[j];
	return;
}

void find_overlapping_reads(uint32_t &current_rid, bitset &current_bitset, std::string &current_quality, std::vector<std::string> &overlap_reads, std::vector<int> &overlap_shift, 
							std::vector<std::string> &overlap_quality, bitset *read, char(*quality)[readlen+1], bbhashdict *dict)
{
	bitset forward_bitset = current_bitset;
	bitset reverse_bitset = stringtobitset(reverse_complement(bitsettostring(forward_bitset)));

	//find matches at shift 0
	find_overlapping_reads_at_shift(0, false, forward_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, read, quality, dict, current_rid);
	find_overlapping_reads_at_shift(0, true, reverse_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, read, quality, dict, current_rid);

	//find matches to right
	for(int shift = 1; shift < maxshift; shift++)
	{	
		forward_bitset >>= 3;
		reverse_bitset <<= 3;
		find_overlapping_reads_at_shift(shift, false, forward_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, read, quality, dict, current_rid);
		find_overlapping_reads_at_shift(shift, true, reverse_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, read, quality, dict, current_rid);
	}

	//find matches to left
	forward_bitset = current_bitset;
	reverse_bitset = stringtobitset(reverse_complement(bitsettostring(forward_bitset)));
	for(int shift = -1; shift > -maxshift; shift--)
	{	
		forward_bitset <<= 3;
		reverse_bitset >>= 3;
		find_overlapping_reads_at_shift(shift, false, forward_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, read, quality, dict, current_rid);
		find_overlapping_reads_at_shift(shift, true, reverse_bitset, current_quality, overlap_reads, overlap_shift, overlap_quality, read, quality, dict, current_rid);
	}
	return;
}

void find_overlapping_reads_at_shift(int shift, bool rev, bitset &current_bitset, std::string &current_quality, std::vector<std::string> &overlap_reads, std::vector<int> &overlap_shift, std::vector<std::string> &overlap_quality, bitset *read, char(*quality)[readlen+1], bbhashdict *dict, uint32_t &current_rid)
{
	std::unordered_set<uint32_t> added_rids;//to make sure rids are not inserted multiple times into the vectors
	int64_t dictidx[2];//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	uint64_t ull;//for indexing dict

	//if shift is 0 and rev is false, add the current read
	if(shift == 0 && rev == false)
	{
		added_rids.insert(current_rid);
		overlap_reads.push_back(bitsettostring(read[current_rid]));
		overlap_quality.push_back(std::string(quality[current_rid]));
		overlap_shift.push_back(0);	
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
		startposidx = dict[l].bphf->lookup(ull);
		if(startposidx >= dict[l].numkeys)//not found
			continue;
		dict[l].findpos(dictidx,startposidx);
		if(dict[l].empty_bin[startposidx])//bin is empty
			continue;
		uint64_t ull1 = ((read[dict[l].read_id[dictidx[0]]]&indexmask[l])>>3*dict_start[l]).to_ullong();
		if(ull == ull1)//checking if ull is actually the key for this bin
		{	
			for (int64_t i = dictidx[0]; i < dictidx[1] && i < dictidx[0] + maxsearch; i++)
			{
				auto rid = dict[l].read_id[i];
				if((current_bitset^(read[rid]&readmask)).count()<=param1)
				{	
					if(added_rids.find(rid) == added_rids.end())//not already present
					{
						added_rids.insert(rid);
						overlap_shift.push_back(shift);	
						if(!rev)	
						{
							overlap_reads.push_back(bitsettostring(read[rid]));
							overlap_quality.push_back(std::string(quality[rid]));
						}
						else
						{
							overlap_reads.push_back(reverse_complement(bitsettostring(read[rid])));
							std::string reverse_quality = quality[rid];
							std::reverse(reverse_quality.begin(),reverse_quality.end());
							overlap_quality.push_back(reverse_quality);
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
		std::getline(f,s);
		read[i] = stringtobitset(s);
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

void constructdictionary(bitset *read, bbhashdict *dict)
{
	for(int j = 0; j < numdict; j++)
	{
		uint64_t *ull = new uint64_t[numreads];
		#pragma omp parallel
		{
		bitset b;
		int tid = omp_get_thread_num();
		std::ofstream foutkey(outdir+std::string("keys.bin.")+std::to_string(tid),std::ios::binary);
		uint32_t i, stop;
		i = uint64_t(tid)*numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads;
		//compute keys and write to file and store in ull
		for(; i < stop; i++)
		{
			b = read[i]&indexmask[j];
			ull[i] = (b>>3*dict_start[j]).to_ullong();
			foutkey.write((char*)&ull[i], sizeof(uint64_t));
		}
		foutkey.close();
		}//parallel end
		
		//deduplicating ull
		std::sort(ull,ull+numreads);
		uint32_t k = 0;
		for (uint32_t i = 1; i < numreads; i++) 
		    if (ull[i] != ull[k])         
				ull[++k] = ull[i];
		dict[j].numkeys = k+1;
		//construct mphf
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(ull), static_cast<const u_int64_t*>(ull+dict[j].numkeys));
		double gammaFactor = 5.0;//balance between speed and memory
		dict[j].bphf = new boomphf::mphf<u_int64_t,hasher_t>(dict[j].numkeys,data_iterator,num_thr,gammaFactor,true,false);
	
		delete[] ull;

		//compute hashes for all reads
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();	
		std::ifstream finkey(outdir+std::string("keys.bin.")+std::to_string(tid),std::ios::binary);
		std::ofstream fouthash(outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
		uint64_t currentkey,currenthash;
		uint32_t i, stop;
		i = uint64_t(tid)*numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads;
		for(; i < stop; i++)
		{
			finkey.read((char*)&currentkey, sizeof(uint64_t));
			currenthash = (dict[j].bphf)->lookup(currentkey);
			fouthash.write((char*)&currenthash, sizeof(uint64_t));
		}
		finkey.close();
		remove((outdir+std::string("keys.bin.")+std::to_string(tid)).c_str());
		fouthash.close();
		}//parallel end
		
	}
		
	//for rest of the function, use numdict threads to parallelize
	omp_set_num_threads(std::min(numdict,num_thr));
	#pragma omp parallel
	{
	#pragma omp for	
	for(int j = 0; j < numdict; j++)
	{
		//fill startpos by first storing numbers and then doing cumulative sum
		dict[j].startpos = new uint32_t[dict[j].numkeys+1]();//1 extra to store end pos of last key
		uint64_t currenthash;
		for(int tid = 0; tid < num_thr; tid++)
		{
			std::ifstream finhash(outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				dict[j].startpos[currenthash+1]++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
		}
	
		dict[j].empty_bin = new bool[dict[j].numkeys]();
		for(uint32_t i = 1; i < dict[j].numkeys; i++)
			dict[j].startpos[i] =  dict[j].startpos[i] +  dict[j].startpos[i-1];

		//insert elements in the dict array
		dict[j].read_id = new uint32_t[numreads];
		uint32_t i = 0;
		for(int tid = 0; tid < num_thr; tid++)
		{
			std::ifstream finhash(outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				dict[j].read_id[dict[j].startpos[currenthash]++] = i;
				i++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
			remove((outdir+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j)).c_str());
		}
		
		//correcting startpos array modified during insertion
		for(int64_t i = dict[j].numkeys; i >= 1 ; i--)	
			dict[j].startpos[i] = dict[j].startpos[i-1];
		dict[j].startpos[0] = 0;
	}//for end
	}//parallel end
	omp_set_num_threads(num_thr);
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

int string_hamming(std::string &a, std::string &b, int &startposa, int &startposb, int &matchlen)
{
    int dist = 0;
    for(int k = 0; k < matchlen; k++)
        if(a[startposa+k] != b[startposb+k])
            dist++;
    return dist;
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

void bbhashdict::findpos(int64_t *dictidx, uint32_t &startposidx)
{
	dictidx[0] = startpos[startposidx];
	auto endidx = startpos[startposidx+1];
	if(read_id[endidx-1] == numreads)//means exactly one read has been removed
		dictidx[1] = endidx-1;
	else if(read_id[endidx-1] == numreads+1)//means two or more reads have been removed (in this case second last entry stores the number of reads left)
		dictidx[1] = dictidx[0] + read_id[endidx-2];
	else
		dictidx[1] = endidx;//no read deleted
	return;
}

void bbhashdict::remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current)
{
	auto size = dictidx[1] - dictidx[0];
	if(size == 1)//just one read left in bin
	{
		empty_bin[startposidx] = 1;
		return; //need to keep one read to check during matching
	}
	uint32_t pos = std::lower_bound(read_id+dictidx[0],read_id+dictidx[1],current)-(read_id+dictidx[0]);
	
	std::move(read_id+dictidx[0]+pos+1,read_id+dictidx[1],read_id+dictidx[0]+pos);
	auto endidx = startpos[startposidx+1];
	if(dictidx[1] == endidx)//this is first read to be deleted
		read_id[endidx-1] = numreads;
	else if(read_id[endidx-1] == numreads)//exactly one read has been deleted till now
	{
		read_id[endidx-1] = numreads + 1;
		read_id[endidx-2] = size - 1;//number of reads left in bin
	}
	else//more than two reads have been deleted
		read_id[endidx-2]--;

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

