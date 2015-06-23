#include <iostream>
#include <string>
#include <vector>

#include "core_functions.h"

using namespace std;

int main(int argc, char * argv[])
{
    string listname = "list.txt";
    string outfilename = "test.vcf";
    
    vector<string> infilename;

    //Parameters params;
    params.sample_count = 0;
    params.type = 3; //polidy
    params.max_count = 255;//Max allele count
    params.thread = 1;
    params.ratio_nchar = 0.2;
    params.ratio_del = 0.05;
    
    params.step = 0.01;//Steps
    params.eps = 0.001;
    params.p_snp = 0.1;
    
    //lFDR output filter with a default of 0.1.  In this way, only sites with an estimated lFDR less than 0.1 should be provided in the final output.  If a user decides to output all analyzed sites, he can change this parameter setting to 1
    params.result_filter = 0.1;
    params.one_circle_limit = 200;
    params.end_condition = 0.5;

    params.bp = 17;
    params.mp = 20;
    
    params.min = 1;
    params.max = 200;
    
    params.start_position = 16691745;

    int arg_pos = 1;

    if (argc == 1)
	    {
		    params.eps = 0.001;
		    //test();
		    printhelp();
	    }

    params.debug = false;
      
    while(arg_pos < argc)
    {
    	if (argv[arg_pos][0] != '-')
    	{
    		cerr << "Argument " << arg_pos << " error : Arguments must start with -" << endl;
    		exit(0);
    	}
    	if (argv[arg_pos][1] == '\0')
    	{
    		cerr << "Argument " << arg_pos << " error : No option found" << endl;
    		exit(0);
    	}

    	int option_pos = arg_pos + 1;
    	switch(argv[arg_pos][1])
    	{
                            case 'i':
                            	listname = argv[option_pos];
                            	break;
                            case 'o':
                            	outfilename = argv[option_pos];
                            	break;
                            case 'n' :
                            	params.ratio_nchar = stof(argv[option_pos]);
                            	break;
                            case 'l' :
                            	params.ratio_del = stof(argv[option_pos]);
                            	break;
                            case 'b' :
                            	params.bp = stoi(argv[option_pos]);
                            	break;
                            case 'm' :
                            	params.mp = stoi(argv[option_pos]);
                            	break;
                            case 's':
                            	params.step = stof(argv[option_pos]);
                            	break;
                            case 'M':
                            	params.max_count = stoi(argv[option_pos]);
                            	break;
                            case 'e':
                            	params.eps = stof(argv[option_pos]);
                            	break;
                            case 'h':
                            	printhelp();
                            	break;
                            case 'f':
                            	params.result_filter = stof(argv[option_pos]);
                            	break;
                            case 't':
                            	params.thread = stoi(argv[option_pos]);
                            	break;
                            case 'C':
                                params.one_circle_limit = stoi(argv[option_pos]);
                                break;
                            default :
                            	cerr<<"Unrec argument: " << argv[arg_pos] << endl;
                            	printhelp();
                            	break;
    	}
    	arg_pos += 2;
    }
    
    params.result_filter = (params.result_filter < 0.0) ? 0.0 : ((params.result_filter > 1.0) ? 1.0 : params.result_filter);
    
    params.sample_count = Get_Name_List(listname.c_str(), infilename);

    if (params.sample_count == 0)
    {
    	cerr << "Sample number must be at least 1" << endl;
    	exit(0);
    }
    
    if (params.max_count < 0)
    {
    	cerr << "Error : Max allele count must larger than 0!" << endl;
    	exit(0);
    }

    calculate_preprocess(infilename, outfilename);

    return 0;
}

