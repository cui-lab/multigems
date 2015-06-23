#include <iostream>
#include <cmath>
#include "seq_obj.h"

float Seq_Obj::Get_Value_Result_Max() {
	if (this->Value.empty()) {
		cerr << "No value\t" << ID << endl;
		return 0.0;
	}
	if (this->Value.size() < Type)
	{
		cerr << "out of range\t" << ID << endl;
		return 0.0;
	}
	float max = Value[0].result;

	for (unsigned int i = 1; i < Type; i++)
		if (Value[i].result > max)
			max = Value[i].result;

	return max;
}

int Seq_Obj::Get_Value_Result_Max_Index()
{
	if (this->Value.empty())
	{
		cerr << "No value\t" << ID << endl;
		return 0;
	}
	if (this->Value.size() < Type)
	{
		cerr << "out of range\t" << ID << endl;
		return 0;
	}

	float max = Value[0].result;
	int max_index = 0;

	for (unsigned int i = 1; i < Type; i++)
		if (Value[i].result > max)
		{
			max = Value[i].result;
			max_index = i;
		}

	return max_index;
}

int Seq_Obj::Seq_Init_Filter()
{
	int iter = 0;
	int valid = 0;
	int length = this->Ref_Info.size();
	string temp_ref_info = "";
	while (iter < length) {
		char test = this->Ref_Info[iter];
		string temp = "";
		switch (test)
		{
			case '+':
				iter++;
				temp = "";
				while ((this->Ref_Info[iter] >= '0') && (this->Ref_Info[iter] <= '9')) {
					temp += Ref_Info[iter];
					iter++;
				}
				iter += stoi(temp);
				break;
			case '-':
				iter++;
				temp = "";
				while ((this->Ref_Info[iter] >= '0') && (this->Ref_Info[iter] <= '9')) {
					temp += Ref_Info[iter];
					iter++;
				}
				iter += atoi(temp.c_str());
				break;
			case '$':
				iter++;
				break;
			case '^':
				iter += 2;
				break;
			default:
				temp_ref_info += test;
				iter++;
				break;
		}
	}
	this->Ref_Info = temp_ref_info;
	if (this->Ref_Info.size() != this->Seq_Qual_1.size())
	{
		valid = 1;
		cout << "diff: " << this->Ref_Info.size() << " " << this->Seq_Qual_1.size() << endl;
	}
	return valid;
}

int Seq_Obj::Seq_Qual_Filter(int bq, int mq)
{
	//Quality
	int iter = 0;
	int length = this->Seq_Qual_1.size();
	string temp_qual_1 = "";
	string temp_qual_2 = "";
	string temp_ref_info = "";

	int qual_bq = bq + 33;
	int qual_mq = mq + 33;

	array<int, 4> allele_count = {{0, 0, 0, 0}};
	array<char, 4> allele_array = {{'A', 'C', 'G', 'T'}};
	while (iter < length)
	{
		if ((this->Ref_Info[iter] != 'N') && (this->Ref_Info[iter] != 'n') && (this->Ref_Info[iter] != '*'))
		{
			if ((this->Seq_Qual_1[iter] >= qual_bq) && (this->Seq_Qual_2[iter] >= qual_mq))
			{
				temp_qual_1 += this->Seq_Qual_1[iter];
				temp_qual_2 += this->Seq_Qual_2[iter];

				if ((this->Ref_Info[iter] != '.') && (this->Ref_Info[iter] != ','))
				{
					char this_allele = toupper(this->Ref_Info[iter]);
					temp_ref_info += this_allele;
					for (int i = 0; i < 4; i++)
						if (this_allele == allele_array[i])
						{
							allele_count[i]++;
							break;
						}
				}
				else
					temp_ref_info += this->Ref;
			}
		}
		iter++;
	}

	this->Seq_Qual_1 = temp_qual_1;
	this->Seq_Qual_2 = temp_qual_2;
	this->Ref_Info = temp_ref_info;

	for (int i = 0; i < 4; i++)
	{
		this->classCounter.at(i) = allele_count[i];
	}

	//Get MAX N
	int max_allele_count = 0;
	char max_allele = 'N';
	for (int i = 0; i < 4; i++)
		if (allele_count[i] > max_allele_count)
		{
			max_allele_count = allele_count[i];
			max_allele = allele_array[i];
		}
	this->Max_allele = max_allele;
	this->Max_allele_count = max_allele_count;

	//RN
	iter = 0;
	length = this->Seq_Qual_1.size();
	temp_qual_1 = "";
	temp_qual_2 = "";
	temp_ref_info = "";

	while (iter < length)
	{
		if (this->Ref_Info[iter] == this->Ref[0])
		{
			temp_ref_info += 'R';
			temp_qual_1 += this->Seq_Qual_1[iter];
			temp_qual_2 += this->Seq_Qual_2[iter];
		}
		else if (this->Ref_Info[iter] == max_allele)
		{
			temp_ref_info += 'N';
			temp_qual_1 += this->Seq_Qual_1[iter];
			temp_qual_2 += this->Seq_Qual_2[iter];
		}
		iter++;
	}

	this->Seq_Qual_1 = temp_qual_1;
	this->Seq_Qual_2 = temp_qual_2;
	this->Ref_Info = temp_ref_info;
	return this->Seq_Qual_1.size();
}

int Seq_Obj::Seq_Max_Filter(const unsigned int max_count)
{
	if ((max_count == 0) || (max_count >= this->Ref_Info.size()))
		return this->Ref_Info.size();
	string temp_qual_1 = "";
	string temp_qual_2 = "";
	string temp_ref_info = "";

	int order[max_count];
	int check = Get_Random(this->Ref_Info.size(), max_count, order);
	for (unsigned int i = 0; i < max_count; i++)
	{
		temp_qual_1 += this->Seq_Qual_1[order[i]];
		temp_qual_2 += this->Seq_Qual_2[order[i]];
		temp_ref_info += this->Ref_Info[order[i]];
	}

	this->Seq_Qual_1 = temp_qual_1;
	this->Seq_Qual_2 = temp_qual_2;
	this->Ref_Info = temp_ref_info;

	return check;

}

float Seq_Obj::Get_Ratio_nchar()
{
	int iter = 0;
	int length = this->Ref_Info.size();
	int valid = 0;
	while (iter < length)
	{
		char test = this->Ref_Info[iter];
		if ((test != '.') && (test != ','))
			valid++;
		iter++;
	}
	return (float) valid / (float) length;
}

float Seq_Obj::Get_Ratio_del()
{
	int iter = 0;
	int length = this->Ref_Info.size();
	int star = 0;
	while (iter < length) {
		char test = this->Ref_Info[iter];
		if (test == '*')
			star++;
		iter++;
	}
	return (float) star / (float) length;
}

void Seq_Obj::Calc_W()
{
	int length = this->Ref_Info.size();

	for (int i = 0; i < length; i++)
		W[i] = 1 - pow(10.0, (((Seq_Qual_1[i] < Seq_Qual_2[i]) ? Seq_Qual_1[i] : Seq_Qual_2[i]) - 33) * (-0.1));
}

int Seq_Obj::Calc_Value(float end, float step)
{
	end = end - step / 10.0;
	Calc_W();

	//Function table
	float (*func_ptr[3][2])(float, float);
	Func_Init(func_ptr);

	// For each genotype of RR and NN
	for (unsigned int i = 0; i < Type - 1; i++)
	{
		this->Value[i].result = MIN;
		this->Value[i].p = 0.0;
		this->Value[i].p2 = 0.0;
		float test_p = step;

		while (test_p < end)
		{
			float test_result = this->Calc_Value(test_p, 0, i);
			if (test_result > this->Value[i].result)
			{
				this->Value[i].result = test_result;
				this->Value[i].p = test_p;
				this->Value[i].p2 = 0;
			}
			test_p += step;
		}
	}

	//Genotype NR and RN

	this->Value[2].result = this->Calc_Value(this->Value[0].p, this->Value[1].p, 2);
	this->Value[2].p = this->Value[0].p;
	this->Value[2].p2 = this->Value[1].p;

	for (unsigned int i = 0; i < this->Type; i++)
	{
		this->valuesVector.at(i * 3 + 0) = this->Value.at(i).result;
		this->valuesVector.at(i * 3 + 1) = this->Value.at(i).p;
		this->valuesVector.at(i * 3 + 2) = this->Value.at(i).p2;
	}

	return 0;
}

float Seq_Obj::Calc_Value(float test_p, float test_p_2, int type)
{
	int length = this->Ref_Info.size();
	//Function table
	float (*func_ptr[3][2])(float, float);
	Func_Init(func_ptr);



	float test_result = 0.0;
	//For each point in the Ref_Info to get the sum of Ref_Info
	for (int j = 0; j < length; j++)
	{
		int num = 0;
		switch (this->Ref_Info[j])
		{
			case 'R':
				num = 0;
				break;
			case 'N':
				num = 1;
				break;
			default:
				break;
		}
		//Get the sum of each row of the table
		float row_sum = 0.0;
		for (int k = 0; k < 2; k++)
			if (k == num)
				row_sum += (*func_ptr[type][k])(test_p, test_p_2) * W[j];
			else
				row_sum += (*func_ptr[type][k])(test_p, test_p_2) * (1.0 - W[j]);                                  ///3.0;

		test_result = test_result + log(row_sum);
	}

	//cout << test_p << "\t" << test_p_2 << "\t" << type << "\t" << test_result << endl;
	return test_result;
}

void Seq_Obj::pre_Calc_Value()
{
	for (int i = 0; i < typeoneVec.size(); i++)
	{
		typeoneVec[i] = this->Calc_Value((i + 1) * this->step, 0, 0);
		typetwoVec[i] = this->Calc_Value((i + 1) * this->step, 0, 1);
		for (int j = 0; j < typeoneVec.size(); j++)
		{
			typethreeVec[i * typeoneVec.size() + j] = this->Calc_Value((i + 1) * this->step, (j + 1) * this->step, 2);
		}
	}
}

float Seq_Obj::get_Calc_Value(float step0, float step1, int type, float step_length)
{
	if (type == 0)
	{
		return this->typeoneVec[floor((step0 - step_length / 10) / step_length)];
	}
	if (type == 1)
	{
		return this->typetwoVec[floor((step0 - step_length / 10) / step_length)];
	}
	if (type == 2)
	{
		//cout << typeoneVec.size() * floor((step0 - step_length / 10) / step_length) + floor((step1 - step_length / 10) / step_length) << "\t";
		return this->typethreeVec[typeoneVec.size() * floor((step0 - step_length / 10) / step_length) + floor((step1 - step_length / 10) / step_length)];
	}
	cerr << "Wrong Type" << endl;
	return  this->typeoneVec[floor((step0 - step_length / 10) / step_length)];
}

int Get_Random(const unsigned int total, const unsigned int n, int * order)
{
	int table[total];
	int order_table[total];

	for (unsigned int i = 0; i < total; i++)
	{
		table[i] = i;
		order_table[i] = 0;
	}
	for (unsigned int i = 0; i < n; i++)
	{
		int j = (int) ((float) (total - 1 - i) * rand() / (RAND_MAX + 1.0));
		order_table[table[j]] = 1;
		table[j] = table[total - 1 - i];
		order[i] = 0;
	}

	int j = 0;

	for (unsigned int i = 0; i < total; i++)
		if (order_table[i] == 1)
		{
			order[j] = i;
			j++;
		}

	return j;
}
