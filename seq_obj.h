#include <string>
#include <array>
#include <vector>
#ifndef SEQ_OBJ_H
#define SEQ_OBJ_H

#define MIN -1000000000.0
#define MAX 1000000000.0

using namespace std;

typedef struct _value {
	float result;
	float p;
	float p2;
} internal_value;

class Seq_Obj {
public:

	friend class Multi_Seq_Obj;

	Seq_Obj(string &_ID, unsigned int _Pos, string &_Ref, int _Num_Two,
			string &_Ref_Info, string &_Seq, string &_Seq_Two, unsigned int type, float step, float end)
	{
		this->ID = _ID;
		this->Pos = _Pos;
		this->Ref = _Ref;
		this->Num_Two = _Num_Two;
		this->Ref_Info = _Ref_Info;
		this->Seq_Qual_1 = _Seq;
		this->Seq_Qual_2 = _Seq_Two;
		this->Type = type;
		initVectors(type, step, end);

		this->Max_allele = 'N';
		this->Max_allele_count = 0;
	}

	Seq_Obj()
	{
		this->ID = "";
		this->Pos = 0;
		this->Ref = "";
		this->Num_Two = -1;
		this->Ref_Info = "";
		this->Seq_Qual_1 = "";
		this->Seq_Qual_2 = "";
		this->Type = 0;

		this->Max_allele = 'N';
		this->Max_allele_count = 0;
	}

	Seq_Obj(vector<string> &obj, unsigned int type, float step, float end)
	{
		this->ID = obj[0];
		this->Pos = stoi(obj[1]);
		this->Ref = obj[2];
		this->Num_Two = stoi(obj[3]);
		this->Ref_Info = obj[4];
		this->Seq_Qual_1 = obj[5];
		this->Seq_Qual_2 = obj[6];
		this->Type = type;
		initVectors(type, step, end);

		this->Max_allele = 'N';
		this->Max_allele_count = 0;
	}

	Seq_Obj(array<string, 7> &obj, unsigned int type, float step, float end)
	{
		this->ID = obj[0];
		this->Pos = stoi(obj[1]);
		this->Ref = obj[2];
		this->Num_Two = stoi(obj[3]);
		this->Ref_Info = obj[4];
		this->Seq_Qual_1 = obj[5];
		this->Seq_Qual_2 = obj[6];
		this->Type = type;
		initVectors(type, step, end);

		this->Max_allele = 'N';
		this->Max_allele_count = 0;
	}

	inline int Get_Max_Allele_Count()
	{
		return this->Max_allele_count;
	}

	inline char Get_Max_Allele()
	{
		return this->Max_allele;
	}

	inline string Get_ID()
	{
		return this->ID;
	}

	inline unsigned int Get_Pos()
	{
		return this->Pos;
	}

	inline int Get_Ref_Length()
	{
		return Ref_Info.size();
	}

	inline string Get_Ref()
	{
		return this->Ref;
	}

	inline string Get_Ref_Info()
	{
		return this->Ref_Info;
	}

	inline string Get_Seq_Qual(int n)
	{
		return (n == 0) ? this->Seq_Qual_1 : this->Seq_Qual_2;
	}

	inline float Get_Value_Result(const unsigned int n)
	{
		if (this->Value.empty())
		{
			cerr << "No value\t" << ID << endl;
			return 0.0;
		}
		if (this->Value.size() <= n)
		{
			cerr << "out of range\t" << ID << endl;
			return 0.0;
		}
		return this->Value[n].result;
	}

	inline float Get_Value_P(const unsigned int n, const int i)
	{
		if (this->Value.empty())
		{
			cerr << "No value\t" << ID << endl;
			return 0.0;
		}
		if (this->Value.size() <= n)
		{
			cerr << "out of range\t" << ID << endl;
			return 0.0;
		}
		if (i == 0)
			return this->Value[n].p;
		else
			return this->Value[n].p2;
	}

	inline vector<int> getClassCounter()
	{
		return this->classCounter;
	}

	inline vector<float> getValuesVector()
	{
		return this->valuesVector;
	}


	int Seq_Init_Filter();
	int Seq_Qual_Filter(int bq, int mq);
	int Seq_Max_Filter(const unsigned int max_count);
	float Get_Ratio_nchar();
	float Get_Ratio_del();
	void Calc_W();
	int Calc_Value(float end, float step);float Calc_Value(float test_p, float test_p_2, int type);
	float Get_Value_Result_Max();
	int Get_Value_Result_Max_Index();

	void pre_Calc_Value();
	float get_Calc_Value(float step0, float step1, int type, float step_length);

private:
	int Num_Two;
	unsigned int Type;
	unsigned int Pos;
	string ID;
	string Ref;
	string Ref_Info;
	string Seq_Qual_1;
	string Seq_Qual_2;
	vector<float> W;
	vector<internal_value> Value;

	char Max_allele;
	int Max_allele_count;

	vector<int> classCounter;
	vector<float> valuesVector;

	//To save the Calc_Values
	vector<float> typeoneVec;
	vector<float> typetwoVec;
	vector<float> typethreeVec;
	float step;
	float end;

	inline void initVectors(int type, float step, float end)
	{
		this->step = step;
		this->end = end;
		this->classCounter.resize(4, 0);
		this->valuesVector.resize(3 * type, 0.0);
		this->W.resize(this->Ref_Info.size(), 0.0);
		this->Value.resize(type);
		this->typeoneVec.resize(floor((end - step / 10) / step));
		this->typetwoVec.resize(floor((end - step / 10) / step));
		this->typethreeVec.resize(floor((end - step / 10) / step) * floor((end - step / 10) / step));
	}
};



inline float func_1(float p1, float p2)
{
	return 1.0 - p1;
}

inline float func_2(float p1, float p2)
{
	return p1;
}

inline float func_3(float p1, float p2)
{
	return 0.5 * (1.0 - p1) + 0.5 * p2;
}

inline float func_4(float p1, float p2)
{
	return 0.5 * (1.0 - p2) + 0.5 * p1;
}

inline int Func_Init(float (* (&func_ptr)[3][2])(float, float))
{
	func_ptr[0][0] = func_1;
    func_ptr[0][1] = func_2;

    func_ptr[1][0] = func_2;
    func_ptr[1][1] = func_1;

    func_ptr[2][0] = func_3;
    func_ptr[2][1] = func_4;

    return 0;
}

int Get_Random(const unsigned int total, const unsigned int n, int * order);

#endif
