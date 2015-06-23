#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include "seq_obj.h"

using namespace std;

#ifndef MULTI_SEQ_OBJ_H
#define MULTI_SEQ_OBJ_H

#define MAX_LOOP 300

class Multi_Seq_Obj {

public:
	Multi_Seq_Obj() {
		Sample = 0;
		Sample_Count = 0;
		Type = 0;
		P = 0;
		P_2 = 0;
		W = -1;
		Is_Qual = false;
		Chrom = "NA";
		Ref = "NA";
	}

	Multi_Seq_Obj(unsigned int n, unsigned int type) {
		Sample = n;
		Sample_Count = 0;
		Type = type;
		Seq_obj_s.resize(Sample);
		Value.resize(Type, 0.0);
		P = 0;
		P_2 = 0;
		E_Value.resize(Sample * Type, 0.0);
		W = -1;
		Is_Qual = false;
		Chrom = "NA";
		Ref = "NA";
	}

	inline int Get_Sample_Count()
	{
		return Sample_Count;
	}

	inline int Get_Type()
	{
		return Type;
	}

	inline float Get_Value(int n)
	{
		return (n < Type) ? Value[n] : 0;
	}

	inline float Get_P(int n)
	{
		return (n == 0) ? P : P_2;
	}

	inline float Get_E_Value(int sample, int type)
	{
		return ((sample < Sample) && (type < Type)) ? E_Value[sample * Type + type] : 0;
	}

	inline string Get_Ref()
	{
		return this->Ref;
	}

	inline bool Get_Is_Sample(int sample)
	{
		return ((sample >= Sample) || (!Seq_obj_s[sample])) ? false : true;
	}

	inline int Get_Sample_Ref_Length(int sample)
	{
		return ((sample >= Sample) || (!Seq_obj_s[sample])) ? 0 : Seq_obj_s[sample].get()->Get_Ref_Length();
	}

	inline float Get_W()
	{
		return this->W;
	}

	inline int Enable()
	{
		if (Is_Qual)
			return 0;
		else {
			Is_Qual = true;
			return 1;
		}
	}

	inline bool Get_Is_Qual()
	{
		return Is_Qual;
	}

	inline void Display(int sample)
	{
		if (!Seq_obj_s[sample])
			cout << "NULL" << endl;
		else {
			cout << Seq_obj_s[sample].get()->Get_Ref_Info() << endl;
			cout << Seq_obj_s[sample].get()->Get_Seq_Qual(0) << endl;
			cout << Seq_obj_s[sample].get()->Get_Seq_Qual(1) << endl;
		}
	}

	inline string Get_Sample(int sample)
	{
		if (sample < Seq_obj_s.size() && Seq_obj_s[sample])
			return Seq_obj_s[sample].get()->Get_Ref_Info();
		else
			return "";
	}

	inline string Get_Chrom()
	{
		return Chrom;
	}

	char Get_Max_Allele();
	int Get_Load(); //Sample * Coverage
	int Insert(shared_ptr<Seq_Obj> &seq_obj, int n);
	int Calc_EM(float end, float step, float eps);
	float Calc_W(int min, int max);
	int Get_Value_Max();
	int Get_E_Value_Max(int sample);

private:
	unsigned int Sample;
	unsigned int Sample_Count;
	unsigned int Type;
	float P;
	float P_2;
	float W;
	string Chrom;
	string Ref;
	vector<shared_ptr<Seq_Obj>> Seq_obj_s;
	vector<float> Value;
	vector<float> E_Value;
	bool Is_Qual;

	void Basic_EM(vector<float> &FS_value, vector<float> &E_value, float end, float step, float &p, float &p_2);
	int Matrix_Norm(vector<float> &m, int w, int h);
	int Matrix_Ave(vector<float> &result, vector<float> &m, int w, int h, int count);
};



#endif
