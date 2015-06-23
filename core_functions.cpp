#include <ctime>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <memory>
#include <queue>
#include <array>
#include <ctime>
#include <omp.h>

#include "core_functions.h"

using namespace std;

map<unsigned int, shared_ptr<Multi_Seq_Obj>> pos_samples_map;
Parameters params;

int printhelp(){
    // Width 12345678911234567892123456789312345678941234567895123456789612345678971234567898
    cout << "MultiGeMS 1.0\n";
    cout << "\n";
    cout << "Overview: Multi-sample Genotype Model Selection (MultiGeMS) is a multiple sample\n";
    cout << "          single nucleotide polymorphism (SNP) caller that works with alignment\n";
    cout << "          files of high-throughput sequencing (HTS) data. MultiGeMS calls SNPs\n";
    cout << "          based on a statistical model selection procedure and accounts for\n";
    cout << "          enzymatic substitution sequencing errors.\n";
    cout << "\n";
    cout << "Input:    MultiGeMS accepts a text file, listing on seperate lines, paths of\n";
    cout << "          SAMtools pileup format files. To convert a SAM/BAM alignment file into\n";
    cout << "          the pileup format, users can use the SAMtools mpileup procedure with\n";
    cout << "          option -s.\n";
    cout << "\n";
    cout << "Filter:   Alignment file reads with undesirable characteristics can be filtered\n";
    cout << "          before running MultiGeMS. For an explanation why filtering may be\n";
    cout << "          desirable and for a brief tutorial on how to filter SAM/BAM alignment\n";
    cout << "          files using SAMtools view, please see the PDF document entitled\n";
    cout << "          \"Pre-Filtering Alignment Files\" available at\n";
    cout << "          https://github.com/cui-lab/multigems.\n";
    cout << "\n";
    cout << "Usage:    multigems -i pileuplist.txt -o multigems.out [OPTIONS]\n";
    cout << "\n";
    cout << "Options:  -b INT   minimum base-calling quality score considered, default is 17\n";
    cout << "          -m INT   minimum mapping quality score considered, default is 20\n";
    cout << "          -s FLOAT maximum likelihood computing steps float value, smaller is\n";
    cout << "                   slower yet more precise, default is 0.01\n";
    cout << "          -e FLOAT EM algorithm convergence threshold, smaller is slower yet\n";
    cout << "                   more precise, default is 0.001\n";
    cout << "          -M INT   maximum number of bases to be considered from each sample of\n";
    cout << "                   each site, 0 indicates unbounded, default is 255\n";
    cout << "          -f FLOAT lFDR SNP threshold, between 0 and 1, default is 0.1\n";
    cout << "          -C INT   number of sites to be analyzed per analysis cycle, smaller\n";
    cout << "                   is slower in general and uses less RAM, default is 200\n";
    cout << "          -t INT   number of threads, can be used in conjunction with -C option\n";
    cout << "                   for multi-thread efficiency, default is 1\n";
    cout << "          -n FLOAT site non-reference allele proportion filter (sites where all\n";
    cout << "                   samples have a non-reference proportion less than this filter\n";
    cout << "                   are not analyzed), smaller is slower and more susceptible to\n";
    cout << "                   false positive SNP calls, higher is faster and more\n";
    cout << "                   susceptible to false negative SNP calls, between 0 and 1,\n";
    cout << "                   default is 0.2\n";
    cout << "          -l FLOAT sample deletion placeholder proportion filter (samples with\n";
    cout << "                   pileup deletion placeholder proportions greater than filter\n";
    cout << "                   are not analyzed), between 0 and 1, default is 0.05\n";
    cout << "\n";
    cout << "Output:   The MultiGeMS output is similar to that of the Variant Call Format\n";
    cout << "          (VCF) file format. Meta-information lines are provided at the\n";
    cout << "          beginning of each output. Only sites less than the user-selected lFDR\n";
    cout << "          SNP threshold are output and are considered SNP calls.\n";
    cout << "\n";
    cout << "Contact:  Xinping Cui <xinping.cui@ucr.edu>\n";
    cout << "\n";
    cout << "          https://sites.google.com/a/bioinformatics.ucr.edu/xinping-cui/\n";
   // Width 12345678911234567892123456789312345678941234567895123456789612345678971234567898
    exit(0);
}

int Get_Name_List(const string &listname, vector<string> &infilename)
{
    ifstream infile(listname, ifstream::in);
    if (!infile)
    {
    	cerr <<"Open list error : " << listname << endl;
    	exit(0);
    }

    string buffer;
    int count = 0;

    while(getline(infile, buffer))
    {
    	if (!buffer.empty())
    	{
    		infilename.push_back(buffer);
    		count++;
    	}
    }

    infile.close();
    return count;
}

int calculate_values(double end)
{
	int count = 0;
	for (map<unsigned int, shared_ptr<Multi_Seq_Obj>>::iterator it = pos_samples_map.begin(); it != pos_samples_map.end(); it++)
	{
		count++;
		it->second.get()->Calc_EM(end, params.step, params.eps);
		it->second.get()->Calc_W(2, 200);
	}
	return count;
}

int calculate_values_omp(double end, int thread)
{
	omp_set_num_threads(thread);
	unsigned int pos_list[pos_samples_map.size()];
	int count = 0;
	for (map<unsigned int, shared_ptr<Multi_Seq_Obj>>::iterator it = pos_samples_map.begin(); it != pos_samples_map.end(); it++)
	{
		pos_list[count] = it->first;
		count++;
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < count; i++)
	{
		if (params.debug) {
			cout << pos_list[i] << endl << endl;
		}
		pos_samples_map[pos_list[i]].get()->Calc_EM(end, params.step, params.eps);
		pos_samples_map[pos_list[i]].get()->Calc_W(2, 200);
	}

	return count;
}

int String_Split(const string &buffer, array<string, 7> &obj, int n)
{
	stringstream strin(buffer);
	for (int i = 0; i < n; i++)
		strin >> obj[i];
	return 0;
}

void calculate_preprocess(const vector<string> &infilename, string &outfilename)
{
	ifstream* ifstream_array = new ifstream[params.sample_count];
	vector<queue<string>> buffer_queue(params.sample_count);

	//open files
	for (int i = 0; i < params.sample_count; i++)
	{
		ifstream_array[i].open(infilename[i], ios::in);
		if (!ifstream_array[i])
		{
			cerr <<"Open infile error : " << infilename[i] << endl;
			exit(0);
		}
	}

	ofstream output_file(outfilename, ios::out);
	if (!output_file)
	{
		cerr << "Open outfile error : " << outfilename << endl;
		exit(0);
	}

	//Let us circle
	core_calculate(ifstream_array, buffer_queue, output_file);

	//close files
	for (int i = 0; i < params.sample_count; i++)
		ifstream_array[i].close();
	output_file.close();

	cout << "GeMS finished, please check the results at " << outfilename << endl;
}

void output_header(ofstream &out)
{
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [9];
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer, 9, "%Y%m%d", timeinfo);
	out << "##fileDate=" << buffer << endl;
	out << "##source=multiGeMSV2.0" << endl;
	out << "##reference=UNKNOWN" << endl;
	out << "##contig=<ID=NA,length=unknown,assembly=NA,md5=NA,species=NA,taxonomy=unknown>" << endl;
	out << "##phasing=partial" << endl;
	out << "##INFO=<ID=NA,Number=NA,Type=NA,Description=\"INFO is not applicable in this version.\"" << endl;
	out << "##FILTER=<ID=q10,Description=\"Quality below 10\"" << endl;
	out << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t"
	    << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << endl;
}

void core_calculate(ifstream* ifstream_array, vector<queue<string>> &buffer_queue, ofstream &output_file)
{
	srand((int)time(NULL));
	vector<unsigned int> count_vector(params.sample_count, 0);
	cout << "Sample number : " << params.sample_count << endl;
	vector<bool> queue_flag(params.sample_count, true);
	int true_flag_count = params.sample_count;

	output_header(output_file);

	int circle_count = 0;
	while (true_flag_count > 0)
	{
		for (int i = 0; i < params.sample_count; i++)
		{
			if (queue_flag[i])
			{
				new_read(ifstream_array[i], buffer_queue[i], buffer_queue[i].size());
				if (buffer_queue[i].size() < params.one_circle_limit)
				{
					queue_flag[i] = false;
					true_flag_count--;
				}
			}
		}

		int checkin_limit = MAX_NUM;

		if (true_flag_count != 0)
			checkin_limit = min_last_element(buffer_queue);
		cout << "checkin_limit = " << checkin_limit << endl;
		clock_t start = clock();
		clock_t end = start;
		for (int i = 0; i < params.sample_count; i++)
		{
			if (buffer_queue[i].size() > 0)
			{
				//cout << i << endl;
				data_checkin(buffer_queue[i], count_vector, checkin_limit, i);
			}
			end = clock();
			cout << end - start << endl;
			start = end;
		}

		cout << endl << ++circle_count << " Loading finished" << endl;
		cout << pos_samples_map.size() << " Pos loaded" << endl;
		cout << position_reduce() << " Pos Removed" << endl;
		cout << pos_samples_map.size() << " Pos remained"<< endl;
		end = clock();
		cout << end - start << endl;
		start = end;
		//calculate_values(params.end_condition);
		calculate_values_omp(params.end_condition, params.thread);
		end = clock();
		cout << end - start << endl;
		start = end;
		cout << "Writing results.." << endl;

		output_values(output_file);
		cout << clock() - start << endl;
	}
}

int new_read(ifstream &in, queue<string> &buffer, int len)
{
	string line;
	while((len < params.one_circle_limit) && (getline(in, line)))
	{
		buffer.push(line);
		len++;
	}
	return len;
}

int min_last_element(vector<queue<string>> &buffer_queue)
{
	int min = MAX_NUM;
	vector<queue<string>>::iterator it = buffer_queue.begin();
	while(it != buffer_queue.end())
	{
		if (!it->empty())
		{
			string line = it->back();
			array<string, 7> obj;
			String_Split(line, obj, 7);
			int pos = stoi(obj[1]);
			if (pos < min)
				min = pos;
		}
		it++;
	}
	return min;
}

void position_add(vector<unsigned int> &count_vector, shared_ptr<Seq_Obj> &seq_obj, int sample)
{
	unsigned int pos = seq_obj.get()->Get_Pos();

	if (pos_samples_map.count(pos) == 0)
	{
		pos_samples_map[pos].reset(new Multi_Seq_Obj(params.sample_count, params.type));
	}
	pos_samples_map[pos].get()->Insert(seq_obj, sample);
	count_vector[sample]++;
}

void data_checkin(queue<string> &buffer, vector<unsigned int> &count_vector, int checkin_limit, int sample)
{
	while (!buffer.empty())
	{
		array<string, 7> obj;
		String_Split(buffer.front(), obj, 7);

		if (stoi(obj[1]) > checkin_limit)
			break;
		//Check if qual_bq_length = qual_mq_length
		if (obj[5].size() != obj[6].size()) cout <<"Qual Seq error : " << buffer.front()<< endl;
		shared_ptr<Seq_Obj> seq_obj(new Seq_Obj(obj, params.type, params.step, params.end_condition));

		if (seq_obj.get()->Get_Ref() == "N") {
			buffer.pop();
			continue;
		}
		if (seq_obj.get()->Seq_Init_Filter() == 1)
			cout << "this is the one with plus error : " << buffer.front() << endl;
		if (seq_obj.get()->Get_Ratio_del() < params.ratio_del)
		{
			if (seq_obj.get()->Get_Ratio_nchar() >= params.ratio_nchar)
			{
				unsigned int temp_pos = seq_obj.get()->Get_Pos();
				seq_obj.get()->Seq_Qual_Filter(params.bp, params.mp);
				seq_obj.get()->Seq_Max_Filter(params.max_count);
				position_add(count_vector, seq_obj, sample);
				pos_samples_map[temp_pos].get()->Enable();
			}
			else
			{
				seq_obj.get()->Seq_Qual_Filter(params.bp, params.mp);
				seq_obj.get()->Seq_Max_Filter(params.max_count);
				position_add(count_vector, seq_obj, sample);
			}
		}
		buffer.pop();
	}
}

unsigned int position_reduce()
{
	unsigned int reduce_count = 0;

	//Remove disabled
	map<unsigned int, shared_ptr<Multi_Seq_Obj>>::iterator it = pos_samples_map.begin();
	while (it != pos_samples_map.end())
	{
		if (!it->second.get()->Get_Is_Qual())
		{
			map<unsigned int, shared_ptr<Multi_Seq_Obj>>::iterator toErase = it;
			it++;
			reduce_count++;
			pos_samples_map.erase(toErase);
		}
		else
		{
			it++;
		}
	}

	return reduce_count;
}

void output_values(ofstream &out)
{
	//char Consensus_letter[11] = "ACGTMRWSYK";

	map<unsigned int, shared_ptr<Multi_Seq_Obj>>::iterator it = pos_samples_map.begin();

	while (it != pos_samples_map.end())
	{
		unsigned int pos = it->first;

		//Let the filter effect
		if (pos_samples_map.at(pos).get()->Get_Is_Qual()
			&& (pos_samples_map.at(pos).get()->Get_W() < params.result_filter))
		{
			out << pos_samples_map.at(pos).get()->Get_Chrom() << "\t";
			out << pos << "\tNA\t" << pos_samples_map.at(pos).get()->Get_Ref() << "\t";
			out << pos_samples_map.at(pos).get()->Get_Max_Allele() << "\t";

			if ((pos_samples_map.at(pos).get()->Get_Sample_Count() >= 1)
				&& (pos_samples_map.at(pos).get()->Get_W() >= 0))
			{
				if (pos_samples_map.at(pos).get()->Get_W() < pow(10, -100))
				{
					out << 999.999 << "\tPASS\t";
				}
				else
				{
					out <<  -10 * log(pos_samples_map.at(pos).get()->Get_W()) << "\tPASS\t";
				}
			}
			else
				out << "NA\tNA\t";
			for (unsigned int i = 0; i < params.sample_count; i++)
			{
				if (pos_samples_map.at(pos).get()->Get_Is_Sample(i)) {
					out << i << ":" << pos_samples_map.at(pos).get()->Get_E_Value_Max(i) + 1 << ",";
				} else
					out << i << ":NA,";
			}

			out << "\t" << "Sample Count: " << pos_samples_map.at(pos).get()->Get_Sample_Count();
			out << "\t" << pos_samples_map.at(pos).get()->Get_P(0);
			out << "\t" << pos_samples_map.at(pos).get()->Get_P(1);
			out << "\t" << pos_samples_map.at(pos).get()->Get_Value(0);
			out << "\t" << pos_samples_map.at(pos).get()->Get_Value(1);
			out << "\t" << pos_samples_map.at(pos).get()->Get_Value(2);
			out << "\t" << pos_samples_map.at(pos).get()->Get_W();
			out << endl;


		}
		it++;
	}
	pos_samples_map.clear();
}


void test()
{
}
