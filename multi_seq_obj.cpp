#include "multi_seq_obj.h"

#include "core_functions.h"

using namespace std;

char Multi_Seq_Obj::Get_Max_Allele()
{
	char max_allele = 'N';
	int max_count = 0;
	array<int, 4> allele_counter = {{0, 0, 0, 0}};
	for (unsigned int i = 0; i < Sample; i++)
	{
		if (Seq_obj_s[i])
		{
			switch (Seq_obj_s[i].get()->Get_Max_Allele())
			{
				case 'A':
					allele_counter[0]++;
					break;
				case 'C':
					allele_counter[1]++;
					break;
				case 'G':
					allele_counter[2]++;
					break;
				case 'T':
					allele_counter[3]++;
					break;
				default:
					break;
			}
		}
	}
	if (allele_counter[0] > max_count)
	{
		max_allele = 'A';
		max_count = allele_counter[0];
	}
	if (allele_counter[1] > max_count)
	{
		max_allele = 'C';
		max_count = allele_counter[1];
	}
	if (allele_counter[2] > max_count)
	{
		max_allele = 'G';
		max_count = allele_counter[2];
	}
	if (allele_counter[3] > max_count)
	{
		max_allele = 'T';
	}

	return max_allele;
}

int Multi_Seq_Obj::Insert(shared_ptr<Seq_Obj> &seq_obj, int n)
{
	if (n > Sample)
		return -1;
	if (Sample_Count == 0)
	{
		this->Ref = seq_obj.get()->Get_Ref();
		this->Chrom = seq_obj.get()->Get_ID();
	}
	if (!Seq_obj_s[n])
		Sample_Count++;
	Seq_obj_s[n].swap(seq_obj);
	return 0;
}

int Multi_Seq_Obj::Get_Value_Max()
{
	float max = Value[0];
	int max_index = 0;
	for (int i = 1; i < Type; i++)
		if (Value[i] > max) {
			max = Value[i];
			max_index = i;
		}
	return max_index;

}

int Multi_Seq_Obj::Get_E_Value_Max(int sample)
{
	if ((sample >= Sample) || (!Seq_obj_s[sample]))
		return -1;
	float max = E_Value[sample * Type + 0];
	int max_index = 0;
	for (int i = 1; i < Type; i++)
		if (E_Value[sample * Type + i] > max)
		{
			max = E_Value[sample * Type + i];
			max_index = i;
		}
	return max_index;
}

int Multi_Seq_Obj::Calc_EM(float end, float step, float eps)
{
	vector<float> E_value(Sample * Type, 0.0);
	vector<float> FS_value(Sample * Type, 0.0);

	vector<float> Init_value(Type, 0.0); //p0
	vector<float> Calc_value(Type, 0.0); //p1

	float Init_p = 0.0;   //q0
	float Init_p_2 = 0.0; //q0
	float Calc_p = 0.0; //q1
	float Calc_p_2 = 0.0; //q1

	//For if only 1 sample
	int Single_sample_index = 0;
	unsigned int Max_value_index = 0;

	for (unsigned int i = 0; i < Sample; i++)
	{
		if (Seq_obj_s[i])
		{
			if (Seq_obj_s[i].get()->Get_Ref_Length() == 0)
			{
				Seq_obj_s[i].reset();
				Sample_Count--;
			}
			else
			{
				Seq_obj_s[i].get()->Calc_Value(end, step);
				for (unsigned int j = 0; j < Type; j++)
					E_value[i * Type + j] = Seq_obj_s[i].get()->Get_Value_Result(j);
				Seq_obj_s[i].get()->pre_Calc_Value();
			}
		}
	}

	if (Sample_Count == 0) {
		Is_Qual = false;
		return 0;
	}

	if (params.debug)
	{
		cout << "Init E_value" << endl;
		for (int i = 0; i < Sample * Type; i++)
			cout << E_value[i] << "\t";
		cout << endl;
	}

	for (int i = 0; i < Sample; i++)
	{
		if (Seq_obj_s[i])
		{
			Single_sample_index = i;
			Max_value_index = Seq_obj_s[i].get()->Get_Value_Result_Max_Index();
			float max = Seq_obj_s[i].get()->Get_Value_Result(Max_value_index);
			for (int j = 0; j < Type; j++)
				E_value[i * Type + j] = exp(Seq_obj_s[i].get()->Get_Value_Result(j) - max);

		}
	}

	Matrix_Norm(E_value, Type, Sample);
	Matrix_Ave(Init_value, E_value, Type, Sample, Sample_Count);

	if (params.debug)
	{
		cout << endl;
		cout << "Before first step of EM" << endl;
		cout << "E_value" << endl;
		for (int i = 0; i < Sample * Type; i++)
			cout << E_value[i] << "\t";
		cout << endl;
		cout << "FS_value" <<endl;
		for (int i = 0; i < Sample * Type; i++)
			cout << FS_value[i] << "\t";
		cout << endl;
		cout << "Init_value" << endl;
		for (int i = 0; i < Type; i++)
			cout << Init_value[i] << endl;
		cout << endl;
		cout << "Calc_value" << endl;
		for (int i = 0; i < Type; i++)
			cout << Calc_value[i] << endl;
		cout << endl;
	}

	//Single Sample
	if (Sample_Count == 1)
	{
		P = Seq_obj_s[Single_sample_index].get()->Get_Value_P(Max_value_index, 0);
		P_2 = Seq_obj_s[Single_sample_index].get()->Get_Value_P(Max_value_index, 1);
		Value = Calc_value;
		E_Value = E_value;
		return 0; //Single GeMS
	}

	Basic_EM(FS_value, E_value, end, step, Init_p, Init_p_2);

	if (params.debug)
	{
		cout << endl;
		cout << "After first step of EM" << endl;
		cout << "E_value" << endl;
		for (int i = 0; i < Sample * Type; i++)
			cout << E_value[i] << "\t";
		cout << endl;
		cout << "FS_value" <<endl;
		for (int i = 0; i < Sample * Type; i++)
			cout << FS_value[i] << "\t";
		cout << endl;
		cout << "Init_value" << endl;
		for (int i = 0; i < Type; i++)
			cout << Init_value[i] << endl;
		cout << endl;
		cout << "Calc_value" << endl;
		for (int i = 0; i < Type; i++)
			cout << Calc_value[i] << endl;
		cout << endl;
	}

	//Loop

	int Loop = 0;
	float diff = MAX;

	while ((diff > eps) && (Loop < MAX_LOOP))
	{
		for (int i = 0; i < Sample; i++)
			if (Seq_obj_s[i])
				for (int j = 0; j < Type; j++)
					E_value[i * Type + j] = exp(FS_value[i * Type + j]) * Init_value[j];


		Matrix_Norm(E_value, Type, Sample);
		Matrix_Ave(Calc_value, E_value, Type, Sample, Sample_Count);

		Basic_EM(FS_value, E_value, end, step, Calc_p, Calc_p_2);

		if (params.debug)
		{
			cout << endl;
			cout << "After Loop " << Loop + 1 << " of EM" << endl;
			cout << "E_value" << endl;
			for (int i = 0; i < Sample * Type; i++)
				cout << E_value[i] << "\t";
			cout << endl;
			cout << "FS_value" <<endl;
			for (int i = 0; i < Sample * Type; i++)
				cout << FS_value[i] << "\t";
			cout << endl;
			cout << "Init_value" << endl;
			for (int i = 0; i < Type; i++)
				cout << Init_value[i] << endl;
			cout << endl;
			cout << "Calc_value" << endl;
			for (int i = 0; i < Type; i++)
				cout << Calc_value[i] << endl;
			cout << endl;
			cout << "Calc_p" << endl;
			cout << Calc_p << endl;
			cout << "Calc_p_2" << endl;
			cout << Calc_p_2 << endl;
		}

		//Check Diff

		diff = 0.0;
		float temp = 0.0;
		for (int i = 0; i < Type; i++)
		{
			temp = fabs(Init_value[i] - Calc_value[i]);
			diff = (diff > temp) ? diff : temp;
		}

		temp = fabs(Init_p - Calc_p);
		diff = (diff > temp) ? diff : temp;
		temp = fabs(Init_p_2 - Calc_p_2);
		diff = (diff > temp) ? diff : temp;

		//Iteratation
		Init_value = Calc_value;
		Init_p = Calc_p;
		Init_p_2 = Calc_p_2;
		Loop++;
	}
	Value = Calc_value;
	E_Value = E_value;

	//min of E_Values
	float min_E_RR = MAX;
	float min_E_NN = MAX;
	for (int i = 0; i < Sample; i++)
	{
		if (E_Value[i * Type] < min_E_RR)
			min_E_RR = E_Value[i * Type];
		if (E_Value[i * Type + 1] < min_E_NN)
			min_E_NN = E_Value[i * Type + 1];
	}

	P = Calc_p; //RN condition 1
	P_2 = Calc_p_2; //NR contiditon 2
	return Loop;
}

void Multi_Seq_Obj::Basic_EM(vector<float> &FS_value, vector<float> &E_value, float end,
		float step, float &p, float &p_2)
{
	end = end - step / 10.0;
	float test_p = step;
	float Sum_max = MIN;

	while (test_p < end)
	{
		float test_p_2 = step;
		while (test_p_2 < end)
		{
			float Sum_temp = 0;
			vector<float> FS_temp(Sample * Type, 0.0);

			for (int i = 0; i < Sample; i++)
				if (Seq_obj_s[i])
				{
					FS_temp[i * Type] = Seq_obj_s[i].get()->get_Calc_Value(test_p, 0, 0, step);
					//cout << test_p << "\t" << test_p_2 << "\t" << "should equal " << FS_temp[i * Type] << "\t" << Seq_obj_s[i].get()->Calc_Value(test_p, 0, 0) << endl;
					FS_temp[i * Type + 1] = Seq_obj_s[i].get()->get_Calc_Value(test_p_2, 0, 1, step);
					//cout << test_p << "\t" << test_p_2 << "\t" << "should equal " << FS_temp[i * Type + 1] << "\t" << Seq_obj_s[i].get()->Calc_Value(test_p_2, 0, 1) << endl;
					FS_temp[i * Type + 2] = Seq_obj_s[i].get()->get_Calc_Value(test_p, test_p_2, 2, step);
					//cout << test_p << "\t" << test_p_2 << "\t" << 	"should equal " << FS_temp[i * Type + 2] << "\t" << Seq_obj_s[i].get()->Calc_Value(test_p, test_p_2, 2) << endl;
					Sum_temp += (FS_temp[i * Type] * E_value[i * Type]);
					Sum_temp += (FS_temp[i * Type + 1] * E_value[i * Type + 1]);
					Sum_temp += (FS_temp[i * Type + 2] * E_value[i * Type + 2]);
				}

			if (Sum_temp * 100.0 >= Sum_max * 100.0)
			{
				p = test_p;
				p_2 = test_p_2;
				FS_value = FS_temp;
				Sum_max = Sum_temp;
			}
			test_p_2 += step;
		}
		test_p += step;
	}
}

float Multi_Seq_Obj::Calc_W(int min, int max) {

	int w = 0;

	float sum = 0.0;

	if (Sample_Count == 1) { //When only 1 sample
		for (int i = 0; i < Sample; i++)
			if (Seq_obj_s[i]) {
				W = E_Value[i * Type + 0];
				return W;
			}
	}

	for (int i = 0; i < Sample; i++)
		if (Seq_obj_s[i])
		{
			int temp_w;
			float temp_e;

			temp_w = (Seq_obj_s[i].get()->Get_Ref_Length() >= min) ? Seq_obj_s[i].get()->Get_Ref_Length() : 0;
			temp_w = (temp_w <= max) ? temp_w : max;
			temp_e = (E_Value[i * Type + 0] * 1000 > 0) ? E_Value[i * Type + 0] : 1e-323;
			w += temp_w;
			sum += (float) temp_w * log(temp_e);
		}

	if (w != 0)
		W = exp(sum / (float) w);

	return W;
}

int Multi_Seq_Obj::Get_Load()
{
	int load = 0;

	for (int i = 0; i < Sample; i++)
		if (Seq_obj_s[i])
			load += Seq_obj_s[i].get()->Get_Ref_Length();

	return load;
}

int Multi_Seq_Obj::Matrix_Norm(vector<float> &m, int w, int h)
{
	for (int i = 0; i < h; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < w; j++)
			sum += m[i * w + j];
		if (sum == 0.0)
			continue;
		for (int j = 0; j < w; j++)
			m[i * w + j] /= sum;
	}

	return 0;
}

int Multi_Seq_Obj::Matrix_Ave(vector<float> &result, vector<float> &m, int w, int h, int count)
{
	for (int i = 0; i < w; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < h; j++)
			sum += m[j * w + i];
		result[i] = sum / count;
	}

	return 0;
}
