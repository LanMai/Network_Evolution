/*
 * Network.h
 *
 *  Created on: 2010/05/27
 *      Author: Amit Lal & Yoshitaka Oku
 */

#ifndef NETWORK_H_
#define NETWORK_H_

# define NoN 82 // total no. of neurons within NG
# define NoNe 80 // total excitatory neurons
# define NoNi 2 // total inhibitory neurons
# define NoPMe 15 // no. of excitatory PM neurons
# define NoPMi 1 // no. of inhibitory PM neurons
# define NoP 80 // no. of neurons to be used for performance analysis of NG
# define synGrade 4 // no. of synaptic grade levels
# define synGradeEE 4 // no. of synaptic grade levels
# define synGradeEI 4 // no. of synaptic grade levels
# define synGradeIE 4 // no. of synaptic grade levels
# define synGradeII 4 // no. of synaptic grade levels


# define binLength 1500 // no of time bins of size equal to interSpikeThreshold
# define NoEB 100 // no. of expected population bursts
# define PopBurstCountMax NoEB
# define PI 3.141592
# define CCxbins 50 // no. of bins into which x-axis (-1 to +1) should be divided
# define Error_Max 5000 // 5000000 (for EvalPopNeuronBurstCC02)
# define SetPopBurstFreq 7 // was initially set at 15

#include "OneNeuron.h"

typedef struct
{
	double synS; // synaptic strength
	double synSinh; // inhibitory synaptic strength
	double SynFrac; // fraction of all-to-all connectivity
	double SynFracEI;
	double SynFracIE;
	double SynFracII;
	double mutFrac; // fraction of synaptic connections to be mutated
	int synDelete; // no. of synaptic connections to be deleted
	int synIncrease;  // no. of synaptic connections to be added
	double tonicConductance; // constant value for gtonic_e
	double tonicConductanceI; // for inhibitory neurons
} netconfig_param;

typedef struct
{
	int NetTag;
	double GNaP[NoN];
	double GLeak[NoN];
	int conMat[NoN][NoN];
} netconfig_cond;

typedef struct
{
	double eval;
	double eval2;
  	double CCx[CCxbins+1];
  	double CCy[CCxbins+1];
} netconfig_result;

namespace pBC {

class Network {
	
	double synS; // synaptic strength
	double synSinh; // inhibitory synaptic strength
	double SynFrac; // fraction of all-to-all connectivity
	double SynFracEI;
	double SynFracIE;
	double SynFracII;
	double mutFrac; // fraction of synaptic connections to be mutated
	int synDelete; // no. of synaptic connections to be deleted
	int synIncrease;  // no. of synaptic connections to be added
	double tonicConductance; // constant value for gtonic_e
	double tonicConductanceI; // constant value for gtonic_e
	double popThresholdFrac; // fraction of NoN that must burst simultaneously to constitute a population burst

	pBC::OneNeuron cp[NoN];

	double V[NoN][2];
	double n[NoN][2];
	double h[NoN][2];
	double s[NoN][2];
	double gtonic_e[NoN];
	double gsyn_e[NoN][NoN];
	int RN1_copy;
	int RN2_copy;
	double dt;

	double t[binLength];
	int BurstData[binLength][NoN+2];
	int PopBurstData[NoEB][NoP];
	int PopBurstDataA[NoEB][NoN]; //only for intentional mutation
	int EctopicBurstData[NoN];
	int BLcount;
	int NoPoPburst;
	double dummyC, Dummy3;
	double CCx[CCxbins+1];
	double CCy[CCxbins+1];
	double CCy2[CCxbins+1];
	double Ey[CCxbins+1];
	int SynAdd, SynDel;
	int EEsum, IIsum, EIsum, IEsum;

	double PopBurstInterval;
	int PopActBeginEnd[PopBurstCountMax][2]; //[PopActBegin End]
	int NeuronActBegin[NoP][PopBurstCountMax]; // actBegin timings of neurons at each popBurst
	int NeuronActEnd[NoP][PopBurstCountMax]; // actEnd timings of neurons at each popBurst
	int NeuronPeriActBegin[NoP][PopBurstCountMax]; //actBegin relative to PopActBegin
	float PeriActBeginVariability[NoP][2]; // [mean & abs(mean_Deviation)/std]
	int NeuronPopBurstData[NoP][PopBurstCountMax]; // matrix representing participation of neurons in PopAct
	float NeuronPopBurstDataVariability[NoP][2]; // [mean & std] of no. of neurons bursting with each popAct

	void AllocatePMprop(double gNaP[], double gLeak[], int N);
	void AllocateNPMprop (double gNaP[], double gLeak[], int N);
	void InitializeNeuronProp();
	void ReadNeuronProp();
	void ReadConMat();
	void CopyEE();
	void NetworkTopology(int MatTag);
	void InitializeNetworkTopology();
	void DescendSort(double RanMat[], int RanIndex[], int NoSyn);
	void DescendSortInteger(int RanMat[], int RanIndex[], int NoSyn);
	void ReadExptCharacteristic();
	double EvalErrorSum();
	double EvalErrorSum2();	
	void BurstTimeOnsetSporadic();
	void BurstTimeOnset();
	void ComputePopBurstData();
	void ComputePopBurstDataA();
	void ComputePopBurstData00();
	void ComputePopBurstData01();

	double SumExcitatoryPop();
	double SumExcitatorySporadic();
	double EvalPopNeuronCC0();
	double EvalPopNeuronCC01();
	double EvalPopNeuronCC02();
	double EvalPopNeuronCC02A();
	double EvalPopNeuronCC02B();
	double EvalPopNeuronCC02C();
	double EvalPopNeuronCC02D();

	void InitializeNetworkTopologyG();
	void NetworkTopologyHybrid(int MatTag);
	void NetworkTopologyG(int MatTag);
	void ConnectMutationSegG(int iter);	

	void MutationEERowSeg(int neuron);
	void MutationEIRowSeg(int neuron);
	void MutationIERowSeg(int neuron);
	void MutationIIRowSeg(int neuron);
	void MutationRowSeg();
	void MutationEEColSeg(int neuron);
	void MutationEIColSeg(int neuron);
	void MutationIEColSeg(int neuron);
	void MutationIIColSeg(int neuron);
	void MutationColSeg();	
	void MutationRowColSeg(int iter);
	void ConnectMutationSeg(int iter);

	void MutationRowSwitch(int iter);
	void MutationColumnSwitch(int iter);
	void MutationNeuronSwitch(int iter);
	void MutationRowTranslocate(int iter);
	void MutationColumnTranslocate(int iter);
	void MutationPartialRowTranslocate(int iter);
	void MutationPartialColumnTranslocate(int iter);

	void MutationColumnSwitchEI(int iter);
	void MutationRowSwitchEI(int iter);
	void MutationColumnSwitchIE(int iter);
	void MutationRowSwitchIE(int iter);		
	void MutationNeuronSwitchING(int iter);
	void MutationColumnSwitchING(int iter);
	void MutationRowSwitchING(int iter);
	void MutationNeuronSwitchENG(int iter);
	void MutationColumnSwitchENG(int iter);
	void MutationPartialColTranslocateENG(int iter);
	void MutationPartialColTranslocateENG1(int iter);
	void MutationPartialColTranslocateEI1(int iter);
	void MutationPartialColTranslocateIE1(int iter);
	void MutationPartialColTranslocateING1(int iter);
	void MutationRowSwitchENG(int iter);
	void MutationPartialRowTranslocateENG(int iter);
	void MutationPartialRowTranslocateENG1(int iter);
	void MutationPartialRowTranslocateEI1(int iter);
	void MutationPartialRowTranslocateIE1(int iter);
	void MutationPartialRowTranslocateING1(int iter);
	void MutationSynSwap(int iter, int EE);
	void MutationSynDeleteEE(int iter);
	void MutationSynDeleteEI(int iter);
	void MutationSynDeleteIE(int iter);
	void MutationSynDeleteII(int iter);
	void MutationSynAddEE(int iter);
	void MutationSynAddEI(int iter);
	void MutationSynAddIE(int iter);
	void MutationSynAddII(int iter);

	void MutationSynColAddENG1(int iter);
	void MutationSynRowAddENG1(int iter);
	void MutationSynColDeleteENG1(int iter);
	void MutationSynRowDeleteENG1(int iter);
	void ColAddII(int iter);
	void RowAddII(int iter);
	void ColDeleteII(int iter);
	void RowDeleteII(int iter);
	void ColAddEI(int iter);
	void RowAddEI(int iter);
	void ColDeleteEI(int iter);
	void RowDeleteEI(int iter);
	void ColAddIE(int iter);
	void RowAddIE(int iter);
	void ColDeleteIE(int iter);
	void RowDeleteIE(int iter);

	double V_now( int icell );
	void Integrate();

	void SmallWorldEE(int n0, double p);
	void SmallWorldEE2(int n0, double p);
	

public:
	Network(netconfig_param netparam, int id, int f_id, int J_id);
	virtual ~Network();

	int process_id;
	int file_id;
	int job_id;
	//double synSinh; // inhibitory synaptic strength
	int NetTag;
	double gNaPe[NoNe];
	double gLeake[NoNe];
	double gNaPi[NoNi];
	double gLeaki[NoNi];
	int conMatEE[NoNe][NoNe];
	int conMatEI[NoNe][NoNi];
	int conMatIE[NoNi][NoNe];
	int conMatII[NoNi][NoNi];
	int conMatEEcopy[NoNe][NoNe];
	int inSynV[NoN];
	int outSynV[NoN];
	
	netconfig_cond Initialize(int z);	
	void SumConMat(int MatTag, bool Pre_Mutation);
	void Execute(double T, double t_dispStep, double t_saveStep, double Tws, const char *filename);
	void ExecuteNoWrite(double T, double t_dispStep, double t_saveStep, double Tws);
	void CopyNetworkCondition(netconfig_cond nc, int z);
	netconfig_cond MutatateNetwork(int iter, int PsuedoMutate);
	netconfig_result EvalPopNeuronBurstCC();
	netconfig_result EvalPopNeuronBurstCC01();
	netconfig_result EvalPopNeuronBurstCC02();
	netconfig_result EvalCC1();
	void writeNetworkCondition(netconfig_cond nc, int a);
	void DoubleSort(double E1double[], double E2[], int Index[], int IndexSorted[], int NoG, int shift);
	netconfig_cond AddInhibitoryConnection();
	netconfig_result EvalVariability();
	
	netconfig_cond InitializeG(int z);
	netconfig_cond MutatateNetworkG(int iter, int PsuedoMutate);
	netconfig_cond MutatateNetworkHybrid(int iter, int PsuedoMutate);
	netconfig_cond MutatateNetworkHybrid0(int iter, int PsuedoMutate);
	netconfig_cond MutatateNetworkHybrid1(int iter, int PsuedoMutate);

	netconfig_cond InitializeSmallWorldEE(int z, int n0, double q);
};

}

#endif /* NETWORK_H_ */
