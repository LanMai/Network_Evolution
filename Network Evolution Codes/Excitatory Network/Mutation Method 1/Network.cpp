/*
 * Network.cpp
 *
 *  Created on: 2010/05/27
 *      Author: Amit Lal & Yoshitaka Oku
 */

#include "Network.h"
#include "OneNeuron.h"
#include "math.h"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;
using std::cerr;
using std::cout;
using std::endl;

namespace pBC {

Network::Network(netconfig_param netparam, int id, int f_id, int J_id) {
	// TODO Auto-generated constructor stub
	process_id = id;
	file_id = f_id;
	job_id = J_id;
	synS = netparam.synS; // synaptic strength
	synSinh = netparam.synSinh; // inhibitory synaptic strength
	SynFrac = netparam.SynFrac; // fraction of all-to-all connectivity
	SynFracEI = netparam.SynFracEI;
	SynFracIE = netparam.SynFracIE;
	SynFracII = netparam.SynFracII;
	mutFrac = netparam.mutFrac; // fraction of synaptic connections to be mutated
	synDelete = netparam.synDelete; // no. of synaptic connections to be deleted
	synIncrease = netparam.synIncrease;  // no. of synaptic connections to be added	
	tonicConductance = netparam.tonicConductance;
	tonicConductanceI = netparam.tonicConductanceI;
	ReadExptCharacteristic();
	popThresholdFrac = 0.3; // fraction of NoN that must burst simultaneously to constitute a population burst
}

Network::~Network() {
	// TODO Auto-generated destructor stub
}


void Network::ReadExptCharacteristic()
{
	// reading experimental data from file
	//=============================================
	ostringstream stream;
	stream << "ExptSyn" << job_id << ".txt";
	ifstream inFile1;
	inFile1.open(stream.str().c_str());
	if (!inFile1)
	{
		cout << "Unable to open file Expt.txt" << endl;
		//exit(1); // terminate with error
	}
	double dummy;
	for (int i=0; i<=CCxbins; i++)
	{
		inFile1 >> dummy;
		Ey[i] = dummy;
	}
	inFile1.close();
	//=================================================
}


double Network::EvalErrorSum()
{

	double ErrorSum = 0;
	/*for (int i=0; i<=CCxbins; i++)
	{
		if ( (CCy[i]-Ey[i])<0 ) ErrorSum = ErrorSum - exp((50-i)/10)*(CCy[i]-Ey[i]);
		else ErrorSum = ErrorSum + exp((50-i)/10)*(CCy[i]-Ey[i]);
		// cout << CCy[i] << "\t" << Ey[i] << endl;
	}*/
	for (int i=0; i<=CCxbins; i++) ErrorSum = ErrorSum + (CCy[i]-Ey[i])*(CCy[i]-Ey[i]);	
	return ErrorSum;
}

double Network::EvalErrorSum2()
{

	double ErrorSum = 0;
	for (int i=0; i<=CCxbins; i++)
	{
		if ( (CCy2[i]-Ey[i])<0 ) ErrorSum = ErrorSum - exp((50-i)/10)*(CCy2[i]-Ey[i]);
		else ErrorSum = ErrorSum + exp((50-i)/10)*(CCy2[i]-Ey[i]);
	}
	return ErrorSum;
}


void Network::BurstTimeOnsetSporadic()
{
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCountMax; j++) 
		{
			NeuronPopBurstData[i][j] = 0;
			NeuronActBegin[i][j] = 0;
			NeuronActEnd[i][j] = 0;
		}
	}

	for (int i=0; i<PopBurstCountMax; i++)
	{
		for (int j=0; j<2; j++) PopActBeginEnd[i][j] = 0; 
	}
	
	int* TonicNeuron;
	TonicNeuron = new int[NoP];
	int ActSum;
	for (int i=0; i<NoP; i++)
	{
		ActSum = 0;
		for (int j=10; j<BLcount; j++)
		{
			ActSum = ActSum + BurstData[j][i];
		}
		if ( ActSum > BLcount-20 ) TonicNeuron[i] = 1;
		else TonicNeuron[i] = 0;
	}

	cout << "process_id: " << process_id  << " BLcount = " << BLcount << endl;
		
	//float popThresholdFrac = 0.3;
	int PopBurstState = 0;
	int PopBurstCount = 0;
	for (int i=10; i<BLcount; i++)
	{
		if ( PopBurstState==0 )
		{
			if ( BurstData[i][NoN+1]>=floor(popThresholdFrac*NoNe) && BurstData[i-1][NoN+1]<floor(popThresholdFrac*NoNe) )
			{
				PopBurstCount = PopBurstCount + 1;
				PopActBeginEnd[PopBurstCount-1][0] = i;
				PopBurstState = 1;
			}
		}

		if ( PopBurstState==1 )
		{
			if ( BurstData[i-1][NoN+1]>=floor(popThresholdFrac*NoNe) && BurstData[i][NoN+1]<floor(popThresholdFrac*NoNe) )
			{
				PopActBeginEnd[PopBurstCount-1][1] = i-1;
				PopBurstState = 0;
			}

			for (int j=0; j<NoP; j++)
			{
				if ( PopBurstState==1 && NeuronPopBurstData[j][PopBurstCount-1]==0 && BurstData[i][j]==1 ) // tonic neurons are not excluded as in BurstTimeOnset()
					NeuronPopBurstData[j][PopBurstCount-1] = 1;
			}
		}		
	}

	NoPoPburst = PopBurstCount;
	if (NoPoPburst>3)
	{
		PopBurstInterval = double( PopActBeginEnd[NoPoPburst-1][0]-PopActBeginEnd[0][0] )/double(NoPoPburst-1);
	}

	//cout << "popThresholdFrac*NoNe = " << floor(popThresholdFrac*NoNe) << endl;
	cout << "process_id: " << process_id  << " PopBurstCount = " << PopBurstCount 
		<< "; PopBurstInterval = " << PopBurstInterval << "; PopBurstInterval(s) = " << PopBurstInterval*160/1000 << endl;
	//cout << "process_id: " << process_id  << " LastPopBurstBegin = " << PopActBeginEnd[NoPoPburst-1][0] << "; FirstPopBurstBegin = " << PopActBeginEnd[0][0] << endl;

	/*ofstream myfile0p;
	ostringstream stream0p;
	stream0p << "NeuronPopBurstData" << file_id << ".txt"; // time & Vsum/NoNe
	myfile0p.open (stream0p.str().c_str());
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCount; j++)
		{
			if (j==PopBurstCount-1) myfile0p << NeuronPopBurstData[i][j] << endl;
			else myfile0p <<  NeuronPopBurstData[i][j] << "\t";
		}
	}
	myfile0p.close();*/

	/*ofstream myfileAp;
	ostringstream streamAp;
	streamAp << "PopActBeginEnd" << file_id << ".txt"; // time & Vsum/NoNe
	myfileAp.open (streamAp.str().c_str());
	for (int i=0; i<PopBurstCount; i++)
	{
		for (int j=0; j<2; j++)
		{
			if (j==1) myfileAp << PopActBeginEnd[i][j] << endl;
			else myfileAp <<  PopActBeginEnd[i][j] << "\t";
		}
	}
	myfileAp.close();*/

	// identifying NeuronAct"Begin"
	int initialT, k;
	for (int i=0; i<PopBurstCount; i++)
	{
		for (int j=0; j<NoP; j++)
		{
			if ( NeuronPopBurstData[j][i]>0 )
			{
				initialT = PopActBeginEnd[i][0];
				k = 0;
				if ( BurstData[initialT][j]>0 )
				{
					if ( TonicNeuron[j]==1 ) 
					{
						if (initialT-10>0) NeuronActBegin[j][i] = initialT - 10;
						else NeuronActBegin[j][i] = 1;
					}
					else
					{
						while ( initialT-k>0 && BurstData[initialT-k][j]>0 ) k = k + 1;
						if ( initialT-k>0 ) NeuronActBegin[j][i] = initialT - (k-1);
						else NeuronActBegin[j][i] = 1;
					}
				}
				else // ( BurstData[initialT][j]==0 )
				{
					while ( initialT+k<BLcount-1 && BurstData[initialT+k][j]<1 ) k = k + 1;
					if ( initialT+k<BLcount-1 ) NeuronActBegin[j][i] = initialT + k;
					else NeuronActBegin[j][i] = BLcount-1;				
				}
			}
		}
	}

	// identifying NeuronAct"End"
	for (int i=0; i<PopBurstCount; i++)
	{
		for (int j=0; j<NoP; j++)
		{
			if ( NeuronPopBurstData[j][i]>0 )
			{
				/*initialT = PopActBeginEnd[i][1];
				k = 0;
				if ( BurstData[initialT][j]>0 )
				{
					if ( TonicNeuron[j]==1 ) 
					{
						if (initialT+10<BLcount-1) NeuronActEnd[j][i] = initialT + 10;
						else NeuronActEnd[j][i] = BLcount;
					}
					else
					{
						while ( initialT+k<BLcount-1 && BurstData[initialT+k][j]>0 ) k = k + 1;
						if ( initialT+k<BLcount-1 ) NeuronActEnd[j][i] = initialT + (k-1);
						else NeuronActEnd[j][i] = BLcount;
					}
				}
				else // ( BurstData[initialT][j]==0 )
				{
					while ( initialT-k>0 && BurstData[initialT+k][j]<1 ) k = k + 1;
					if ( initialT-k>0 ) NeuronActEnd[j][i] = initialT - k; // this condition is always satisfied
				}*/
				initialT = NeuronActBegin[j][i];
				k = 0;
				if ( TonicNeuron[j]==1 ) 
				{
					if (PopActBeginEnd[i][1]+10<BLcount-1) NeuronActEnd[j][i] = PopActBeginEnd[i][1] + 10;
					else NeuronActEnd[j][i] = BLcount;
				}
				else
				{
					while ( initialT+k<BLcount-1 && BurstData[initialT+k][j]>0 ) k = k + 1;
					if ( initialT+k<BLcount-1 ) NeuronActEnd[j][i] = initialT + (k-1);
					else NeuronActEnd[j][i] = BLcount;
				}
			}
		}
	}

	int SporadicBurstPossible = 0;
	int SporadicBurstCount = 0;
	for (int i=0; i<PopBurstCount-1; i++)
	{
		for (int j=0; j<NoP; j++)
		{
			if ( NeuronPopBurstData[j][i]>0 )
			{
				for (int k=NeuronActEnd[j][i]+1; k<NeuronActBegin[j][i+1]; k++ )
				{
					SporadicBurstCount = SporadicBurstCount + BurstData[k][j];
					SporadicBurstPossible = SporadicBurstPossible + 1;
				}
			}
			else
			{
				for (int k=PopActBeginEnd[i][1]; k<PopActBeginEnd[i+1][0]; k++ )
				{	
					SporadicBurstCount = SporadicBurstCount + BurstData[k][j];
					SporadicBurstPossible = SporadicBurstPossible + 1;
				}
			}
		}
	}

	cout << "process_id: " << process_id  << "\t SporadicBurstPossible = " << SporadicBurstPossible 
		<< "\t SporadicBurstCount = " << SporadicBurstCount 
		<< "\t SporadicFrac = " << 100*double(SporadicBurstCount)/double(SporadicBurstPossible) << endl;



	delete [] TonicNeuron;

	/*ofstream myfileBp;
	ostringstream streamBp;
	streamBp << "NeuronActBegin" << file_id << ".txt";
	myfileBp.open (streamBp.str().c_str());
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCount; j++)
		{
			if (j==PopBurstCount-1) myfileBp << NeuronActBegin[i][j] << endl;
			else myfileBp <<  NeuronActBegin[i][j] << "\t";
		}
	}
	myfileBp.close();*/

	/*ofstream myfileCp;
	ostringstream streamCp;
	streamCp << "NeuronActEnd" << file_id << ".txt";
	myfileCp.open (streamCp.str().c_str());
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCount; j++)
		{
			if (j==PopBurstCount-1) myfileCp << NeuronActEnd[i][j] << endl;
			else myfileCp <<  NeuronActEnd[i][j] << "\t";
		}
	}
	myfileCp.close();*/

}


void Network::BurstTimeOnset()
{
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCountMax; j++) NeuronPopBurstData[i][j] = 0;
	}

	for (int i=0; i<PopBurstCountMax; i++)
	{
		for (int j=0; j<2; j++) PopActBeginEnd[i][j] = 0; 
	}

	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCountMax; j++) NeuronActBegin[i][j] = 0; 
	}

	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<2; j++) PeriActBeginVariability[i][j] = 0; 
	}
	
	int* TonicNeuron;
	TonicNeuron = new int[NoP];
	int ActSum;
	for (int i=0; i<NoP; i++)
	{
		ActSum = 0;
		for (int j=10; j<BLcount; j++)
		{
			ActSum = ActSum + BurstData[j][i];
		}
		if ( ActSum > BLcount-20 ) TonicNeuron[i] = 1;
		else TonicNeuron[i] = 0;
	}
	
	//float popThresholdFrac = 0.3;
	int PopBurstState = 0;
	int PopBurstCount = 0;
	for (int i=10; i<BLcount; i++)
	{
		if ( PopBurstState==0 )
		{
			if ( BurstData[i][NoN+1]>=floor(popThresholdFrac*NoNe) && BurstData[i-1][NoN+1]<floor(popThresholdFrac*NoNe) )
			{
				PopBurstCount = PopBurstCount + 1;
				PopActBeginEnd[PopBurstCount-1][0] = i;
				PopBurstState = 1;
			}
		}

		if ( PopBurstState==1 )
		{
			if ( BurstData[i-1][NoN+1]>=floor(popThresholdFrac*NoNe) && BurstData[i][NoN+1]<floor(popThresholdFrac*NoNe) )
			{
				PopActBeginEnd[PopBurstCount-1][1] = i-1;
				PopBurstState = 0;
			}

			for (int j=0; j<NoP; j++)
			{
				if ( PopBurstState==1 && NeuronPopBurstData[j][PopBurstCount-1]==0 && BurstData[i][j]==1 && TonicNeuron[j]==0) 
					NeuronPopBurstData[j][PopBurstCount-1] = 1;
			}
		}		
	}

	NoPoPburst = PopBurstCount;

	//cout << "popThresholdFrac*NoNe = " << floor(popThresholdFrac*NoNe) << endl;
	//cout << "PopBurstCount = " << PopBurstCount << endl;

	/*ofstream myfile0p;
	ostringstream stream0p;
	stream0p << "NeuronPopBurstData" << file_id << ".txt"; // time & Vsum/NoNe
	myfile0p.open (stream0p.str().c_str());
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCount; j++)
		{
			if (j==PopBurstCount-1) myfile0p << NeuronPopBurstData[i][j] << endl;
			else myfile0p <<  NeuronPopBurstData[i][j] << "\t";
		}
	}
	myfile0p.close();*/

	/*ofstream myfileAp;
	myfileAp.open ("PopActBeginEnd.txt");
	for (int i=0; i<PopBurstCount; i++)
	{
		for (int j=0; j<2; j++)
		{
			if (j==1) myfileAp << PopActBeginEnd[i][j] << endl;
			else myfileAp <<  PopActBeginEnd[i][j] << "\t";
		}
	}
	myfileAp.close();*/

	int dummyA;
	for (int i=0; i<NoP; i++)
	{
		dummyA = 0;
		for (int j=0; j<PopBurstCount; j++) 
		{
			dummyA = dummyA + NeuronPopBurstData[i][j];
		}
		NeuronPopBurstDataVariability[i][0] = float(dummyA)/float(PopBurstCount); 
		// mean no. of neurons participating in each popAct
	}

	float dummyB;
	dummyC = 0;
	for (int i=0; i<NoP; i++)
	{
		dummyB = 0;
		for (int j=0; j<PopBurstCount; j++) 
		{
			dummyB = dummyB + pow(float(NeuronPopBurstData[i][j])-NeuronPopBurstDataVariability[i][0],2);
		}
		NeuronPopBurstDataVariability[i][1] = sqrt(dummyB/float(PopBurstCount));
		dummyC = dummyC + double(pow(NeuronPopBurstDataVariability[i][1],2));
		// total variation of no. of neurons participating in each popAct
	}
	//dummyC = sqrt(dummyC/double(NoP));
	cout << "process_id: " << process_id << " dummyC = " << dummyC << endl;


	int initialT, k;
	for (int i=0; i<PopBurstCount; i++)
	{
		for (int j=0; j<NoP; j++)
		{
			if ( NeuronPopBurstData[j][i]>0 )
			{
				initialT = PopActBeginEnd[i][0];
				k = 0;
				if ( BurstData[initialT][j]>0 )
				{
					while ( initialT-k>0 && BurstData[initialT-k][j]>0 ) k = k + 1;
					if ( initialT-k>0 ) NeuronActBegin[j][i] = initialT - (k-1);
				}
				else // ( BurstData[initialT][j]==0 )
				{
					while ( initialT+k<BLcount-1 && BurstData[initialT+k][j]<1 ) k = k + 1;
					if ( initialT+k<BLcount-1 ) NeuronActBegin[j][i] = initialT + k;
				
				}
			}
		}
	}

	ofstream myfileBp;
	ostringstream streamBp;
	streamBp << "NeuronActBegin" << file_id << ".txt";
	myfileBp.open (streamBp.str().c_str());
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCount; j++)
		{
			if (j==PopBurstCount-1) myfileBp << NeuronActBegin[i][j] << endl;
			else myfileBp <<  NeuronActBegin[i][j] << "\t";
		}
	}
	myfileBp.close();

	
	int dummy1, dummy2;
	for (int i=0; i<NoP; i++)
	{
		dummy1 = 0;
		dummy2 = 0;
		for (int j=0; j<PopBurstCount; j++) 
		{
			if (NeuronPopBurstData[i][j]>0) 
			{
				NeuronPeriActBegin[i][j] = NeuronActBegin[i][j] - PopActBeginEnd[j][0];
				dummy1 = dummy1 + NeuronPeriActBegin[i][j];
				dummy2 = dummy2 + NeuronPopBurstData[i][j];
			}
		}
		if (dummy2>0) PeriActBeginVariability[i][0] = float(dummy1)/float(dummy2);
	}

	float Dummy1, Dummy2; 
	Dummy3 = 0;
	// cout << "Dummy1 = " << Dummy1 << " Dummy3 = " << Dummy3 << endl;
	for (int i=0; i<NoP; i++)
	{
		Dummy1 = 0;
		Dummy2 = 0;
		for (int j=0; j<PopBurstCount; j++) 
		{
			if (NeuronPopBurstData[i][j]>0)
			{
				/*if (float(NeuronPeriActBegin[i][j])>PeriActBeginVariability[i][0]) Dummy1 = Dummy1 + float(NeuronPeriActBegin[i][j]) - PeriActBeginVariability[i][0];
				else Dummy1 = Dummy1 - float(NeuronPeriActBegin[i][j]) + PeriActBeginVariability[i][0];*/
				Dummy1 = Dummy1 + pow(float(NeuronPeriActBegin[i][j]) - PeriActBeginVariability[i][0],2);
				Dummy2 = Dummy2 + float(NeuronPopBurstData[i][j]);
			}
		}
		if (Dummy2>0)
		{
			PeriActBeginVariability[i][1] = sqrt(Dummy1/Dummy2);
			Dummy3 = Dummy3 + double(Dummy1/Dummy2);
		}
		else PeriActBeginVariability[i][1] = 0; // redundant!
	}

	cout << "process_id: " << process_id << " Dummy3 = " << Dummy3 << endl;

	delete [] TonicNeuron;

	ofstream myfileCp;
	ostringstream streamCp;
	streamCp << "NeuronPeriActBegin" << file_id << ".txt";
	myfileCp.open (streamCp.str().c_str());
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<PopBurstCount; j++)
		{
			if (j==PopBurstCount-1) myfileCp << NeuronPeriActBegin[i][j] << endl;
			else myfileCp <<  NeuronPeriActBegin[i][j] << "\t";
		}
	}
	myfileCp.close();

	/*ofstream myfileDp;
	myfileDp.open ("PeriActBeginVariability.txt");
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<2; j++)
		{
			if (j==1) myfileDp << PeriActBeginVariability[i][j] << endl;
			else myfileDp <<  PeriActBeginVariability[i][j] << "\t";
		}
	}
	myfileDp.close();*/

	/*ofstream myfileEp;
	myfileEp.open ("NeuronPopBurstDataVariability.txt");
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<2; j++)
		{
			if (j==1) myfileEp << NeuronPopBurstDataVariability[i][j] << endl;
			else myfileEp <<  NeuronPopBurstDataVariability[i][j] << "\t";
		}
	}
	myfileEp.close();*/

}


netconfig_result Network::EvalVariability()
{
	netconfig_result nr;

	BurstTimeOnset();

	if (NoPoPburst <= 7)
	{
		nr.eval = Error_Max;
		nr.eval2 = dummyC;
		//nr.eval2 = Dummy3;
		//nr.eval2 = double(NoNe*NoNe*synGrade-EEsum);
		//nr.eval2 = 0;
		for (int i=0; i<=CCxbins; i++)
		{
			nr.CCx[i] = CCx[i];
			nr.CCy[i] = CCy[i];
		}
	}
	else // (NoPoPburst > 7)
	{
		nr.eval = Error_Max - Dummy3;
		nr.eval2 = dummyC;
		// dummyC = total variation of no. of neurons participating in each popAct
		//nr.eval2 = Dummy3;
		for (int i=0; i<=CCxbins; i++)
		{
			nr.CCx[i] = 0;
			if ( i==0 ) nr.CCy[i] = double(NoPoPburst);
			else nr.CCy[i] = 0;
		}
	}
	
	return nr;
}


void Network::ComputePopBurstData()
{
	cout << "process_id: " << process_id
			<< "  BLcount = " << BLcount << endl;
	NoPoPburst = 0;
	int actBegin = 0;
	int actEnd = 0;
	int j;
	int popThreshold = int(round(popThresholdFrac*NoP));
	int Bstart[NoEB];
	int Bend[NoEB];

	//ofstream file1a;
	//file1a.open ("actTiming.txt");
	for (int i=0; i<BLcount-1; i++)
	{
		if (BurstData[i][NoN]>popThreshold) //(BurstData[i][NoN] contains RasterSum of all neurons
		{
			actBegin = i;
		}
		if (BurstData[i][NoN]>=popThreshold)
		{
			NoPoPburst = NoPoPburst + 1;
			j=0;
			while ( BurstData[i+j][NoN]>=popThreshold && (i+j)<BLcount-1 )
			{
				j = j + 1;
				actEnd = i+j;
			}

			int NNpopLength = actEnd-actBegin+1;
			int* NNpop;
			NNpop = new int[NNpopLength];
			int* NNpopIndex;
			NNpopIndex = new int[NNpopLength];
			for (int k=actBegin; k<=actEnd; k++)
			{
				NNpop[k-actBegin] = BurstData[k][NoN];
				NNpopIndex[k-actBegin] = k;
			}
			DescendSortInteger(NNpop, NNpopIndex, NNpopLength); // after sorting t[NNpopIndex[0]] is the timing when population burst is maximum
			//file1a << NoPoPburst << "\t" << t[NNpopIndex[0]] << "\t" << actBegin+1 << "\t" << actEnd-1 << "\t" << t[actBegin+1] << "\t" << t[actEnd-1] << endl;
			Bstart[NoPoPburst-1] = actBegin+1;
			Bend[NoPoPburst-1] = actEnd-1;
			delete NNpop;
			delete NNpopIndex;

			i = actEnd;
		}
	}
	//file1a.close();

	//int PopBurstData[NoEB][NoP];
	for (int i=0; i<NoEB; i++)
	{
		for (int j=0; j<NoP; j++)
		{
			PopBurstData[i][j] = 0;
		}
	}

	// the following is to discard incomplete population bursts
	// FirstBurstStart is index of first complete population burst
	// LastBurstEnd is index of last complete population burst
	int FirstBurstStart = 0;
	int LastBurstEnd = NoPoPburst;
	if (Bstart[0]==0) FirstBurstStart = 1;
	if (Bstart[NoPoPburst-1]==BLcount-1) LastBurstEnd = NoPoPburst-1;
	NoPoPburst = LastBurstEnd - FirstBurstStart;
	cout << "process_id : " << process_id << "\t NoPoPburst = " << NoPoPburst << endl;

	int Tcount;
	for (int j=0; j<NoP; j++)
	{
		for (int i=0; i<NoPoPburst; i++)
		{
			Tcount = Bstart[i+FirstBurstStart];
			while ( Tcount<=Bend[i+FirstBurstStart] && PopBurstData[i][j]==0)
			{
				if (BurstData[Tcount][j] == 1) PopBurstData[i][j] = 1;
				Tcount = Tcount + 1;
			}
		}
	}

	/*ofstream myfile0q;
	ostringstream stream2;
	stream2 << "PopBurstData" << process_id << ".txt";
	myfile0q.open (stream2.str().c_str());
	for (int i=0; i<NoEB; i++)
	{
		for (int j=0; j<NoP; j++)
		{
			if (j==NoP-1) myfile0q << PopBurstData[i][j] << endl;
			else myfile0q <<  PopBurstData[i][j] << "\t";
		}
	}
	myfile0q.close();*/

}


void Network::ComputePopBurstData00()
{
	/* same as ComputePopBurstData except that popThreshold is 
	computed taking only excitatory neurons into account */

	cout << "process_id: " << process_id
			<< "  BLcount = " << BLcount << endl;
	NoPoPburst = 0;
	int actBegin = 0;
	int actEnd = 0;
	int j;
	int popThreshold = int(round(popThresholdFrac*NoNe));
	int Bstart[NoEB];
	int Bend[NoEB];

	//ofstream file1a;
	//file1a.open ("actTiming.txt");
	for (int i=0; i<BLcount-1; i++)
	{
		if (BurstData[i][NoN+1]>popThreshold) //(BurstData[i][NoN] contains RasterSum of all neurons
		{
			actBegin = i;
		}
		if (BurstData[i][NoN+1]>=popThreshold)
		{
			NoPoPburst = NoPoPburst + 1;
			j=0;
			while ( BurstData[i+j][NoN+1]>=popThreshold && (i+j)<BLcount-1 )
			{
				j = j + 1;
				actEnd = i+j;
			}

			int NNpopLength = actEnd-actBegin+1;
			int* NNpop;
			NNpop = new int[NNpopLength];
			int* NNpopIndex;
			NNpopIndex = new int[NNpopLength];
			for (int k=actBegin; k<=actEnd; k++)
			{
				NNpop[k-actBegin] = BurstData[k][NoN+1];
				NNpopIndex[k-actBegin] = k;
			}
			DescendSortInteger(NNpop, NNpopIndex, NNpopLength); // after sorting t[NNpopIndex[0]] is the timing when population burst is maximum
			//file1a << NoPoPburst << "\t" << t[NNpopIndex[0]] << "\t" << actBegin+1 << "\t" << actEnd-1 << "\t" << t[actBegin+1] << "\t" << t[actEnd-1] << endl;
			Bstart[NoPoPburst-1] = actBegin+1;
			Bend[NoPoPburst-1] = actEnd-1;
			delete NNpop;
			delete NNpopIndex;

			i = actEnd;
		}
	}
	//file1a.close();

	//int PopBurstData[NoEB][NoP];
	for (int i=0; i<NoEB; i++)
	{
		for (int j=0; j<NoNe; j++)
		{
			PopBurstData[i][j] = 0;
		}
	}

	// the following is to discard incomplete population bursts
	// FirstBurstStart is index of first complete population burst
	// LastBurstEnd is index of last complete population burst
	int FirstBurstStart = 0;
	int LastBurstEnd = NoPoPburst;
	if (Bstart[0]==0) FirstBurstStart = 1;
	if (Bstart[NoPoPburst-1]==BLcount-1) LastBurstEnd = NoPoPburst-1;
	NoPoPburst = LastBurstEnd - FirstBurstStart;
	cout << "process_id : " << process_id << "\t NoPoPburst = " << NoPoPburst << endl;

	int Tcount;
	for (int j=0; j<NoNe; j++)
	{
		for (int i=0; i<NoPoPburst; i++)
		{
			Tcount = Bstart[i+FirstBurstStart];
			while ( Tcount<=Bend[i+FirstBurstStart] && PopBurstData[i][j]==0)
			{
				if (BurstData[Tcount][j] == 1) PopBurstData[i][j] = 1;
				Tcount = Tcount + 1;
			}
		}
	}

	/*ofstream myfile0q;
	ostringstream stream2;
	stream2 << "PopBurstData" << process_id << ".txt";
	myfile0q.open (stream2.str().c_str());
	for (int i=0; i<NoEB; i++)
	{
		for (int j=0; j<NoNe; j++)
		{
			if (j==NoNe-1) myfile0q << PopBurstData[i][j] << endl;
			else myfile0q <<  PopBurstData[i][j] << "\t";
		}
	}
	myfile0q.close();*/
	

}


void Network::ComputePopBurstData01()
{
	/* same as ComputePopBurstData except that popThreshold is 
	computed taking only excitatory neurons into account */

	cout << "process_id: " << process_id
			<< "  BLcount = " << BLcount << endl;
	NoPoPburst = 0;
	int actEnd = 0;
	int j;
	int popThreshold = int(round(popThresholdFrac*NoNe));
	
	for (int i=0; i<BLcount-1; i++)
	{		
		if (BurstData[i][NoN+1]>=popThreshold)
		{
			NoPoPburst = NoPoPburst + 1;
			j=0;
			while ( BurstData[i+j][NoN+1]>=popThreshold && (i+j)<BLcount-1 )
			{
				j = j + 1;
				actEnd = i+j;
			}

			i = actEnd;
		}
	}	

}


double Network::SumExcitatoryPop()
{
	int popThreshold = int(round(popThresholdFrac*NoNe));
	int Ci;
	double CrossCorrSum = 0;
	/*for (int i=0; i<NoNe; i++)
	{
		Ci = 0;
		for (int k=0; k<BLcount; k++)
		{
			if ( BurstData[k][NoN+1]>=popThreshold )
			{
				Ci = Ci + BurstData[k][i];
			}
		}
		CrossCorrSum = CrossCorrSum + double(Ci);
	}*/

	int BLpopBurst = 0;
	for (int k=0; k<BLcount; k++)
	{
		if ( BurstData[k][NoN+1]>=popThreshold ) 
		{
			CrossCorrSum = CrossCorrSum + (NoNe - BurstData[k][NoN+1]);
			BLpopBurst = BLpopBurst + 1;
		}
	}
	CrossCorrSum = CrossCorrSum/double(BLpopBurst);

	return CrossCorrSum;
}


double Network::SumExcitatorySporadic()
{
	int popThreshold = int(round(popThresholdFrac*NoNe));
	int Ci;
	double CrossCorrSum = 0;
	for (int i=0; i<NoNe; i++)
	{
		Ci = 0;
		for (int k=0; k<BLcount; k++)
		{
			if ( BurstData[k][NoN+1]<popThreshold )
			{
				Ci = Ci + BurstData[k][i];
			}
		}
		CrossCorrSum = CrossCorrSum + double(Ci);
	}

	return CrossCorrSum;
}


double Network::EvalPopNeuronCC0()
{
	int popThreshold = int(round(popThresholdFrac*NoP));
	int Ci;
	double CrossCorrSum = 0;
	for (int i=0; i<NoP; i++)
	{
		Ci = 0;
		for (int k=0; k<BLcount; k++)
		{
			Ci = Ci + BurstData[k][i];
			/*if ( (BurstData[k][NoN+1]>=popThreshold) && (BurstData[k][i]==1) )
			{
				Ci = Ci + 1;
			}
			else if ( (BurstData[k][NoNe]<popThreshold) && (BurstData[k][i]==1) )
			{
				Ci = Ci - 1;
			}
			else if ( (BurstData[k][NoNe]>=popThreshold) && (BurstData[k][i]==0) )
			{
				Ci = Ci - 1;
			}*/
		}
		CrossCorrSum = CrossCorrSum + double(Ci);
	}

	return CrossCorrSum;
}



double Network::EvalPopNeuronCC01()
{
	int popThreshold = int(round(popThresholdFrac*NoP));
	double popActMean = 0;
	int* popAct;
	popAct = new int[BLcount];
	for (int i=0; i<BLcount-1; i++)
	{
		if ( BurstData[i][NoP]>=popThreshold )	popAct[i] = 1; 
		else popAct[i] = 0;
		popActMean = popActMean + double(popAct[i]);
	}
	popActMean = popActMean/double(BLcount);
		
	double* CCvector;
	CCvector = new double[NoP];
	int* CCindex;
	CCindex = new int[NoP];
	double CCmean;
	double CCi, CCj;
	for (int i=0; i<NoP; i++)
	{
		CCmean = 0;
		for (int j=0; j<BLcount; j++)
		{
			CCmean = CCmean + BurstData[j][i];
		}
		CCmean = double(CCmean)/double(BLcount);
		if (CCmean > 0.999999) CCmean = 0.999999;
		if (CCmean < 0.000001) CCmean = 0.000001;

		CCvector[i] = 0;
		CCi = 0;
		CCj = 0;
		for (int j=0; j<BLcount; j++)
		{
			CCvector[i] = CCvector[i] + (BurstData[j][i]-CCmean)*(popAct[j]-popActMean);
			CCi = CCi + pow( (BurstData[j][i]-CCmean),2 );
			CCj = CCj + pow( (popAct[j]-popActMean),2 );
		}
		CCvector[i] = CCvector[i]/sqrt(CCi*CCj);
		CCindex[i] = i;
	}

	DescendSort(CCvector,CCindex,NoP);

	//cout << "process_id: " << process_id << "  NoPoPburst = " << NoPoPburst << endl;
	for (int i=0; i<=CCxbins; i++)
	{
		CCx[i] = -1.0 + double(i)*(2.0/double(CCxbins));
	}

	int count = 0;
	for (int i=0; i<=CCxbins; i++)
	{
		if (NoP-1-count>=0)
		{
			while (CCvector[NoP-1-count]<=CCx[i])
			{
				count = count + 1;
				if (NoP-1-count<0) break;
			}
			if (count==0) CCy2[i] = 0;
			else CCy2[i] = double(count)/double(NoP);			
		}
		else CCy2[i] = CCy2[i-1];
	}

	delete[] CCvector;
	delete[] CCindex;


	double eval2 = EvalErrorSum2();
	
	return eval2;
}

double Network::EvalPopNeuronCC02C()
{
	// to effect ectopic bursts
	/* Sums up total neuronal activity;
	slight variation from EvalPopNeuronCC02B(); no threshold limit
	only excitatory neurons are considered here */

	int popThreshold = int(round(popThresholdFrac*NoNe));
	double NeuronActSum = 0;
	int* popAct;
	popAct = new int[BLcount];
	for (int i=0; i<BLcount-1; i++)
	{
		popAct[i] = BurstData[i][NoN+1];
		NeuronActSum = NeuronActSum + double(popAct[i]);
	}
	
	delete[] popAct;

	/*cout << "process_id: " << process_id << " popThreshold = " << popThreshold 
		<< "  NeuronActSum = " << NeuronActSum << endl;*/
	
	return NeuronActSum;
}

double Network::EvalPopNeuronCC02B()
{
	// to effect ectopic bursts
	/* Sums up total neuronal activity;
	during population burst total active neurons = threshold neurons
	only excitatory neurons are considered here */

	int popThreshold = int(round(popThresholdFrac*NoNe));
	double NeuronActSum = 0;
	//double Weight = 0;
	int* popAct;
	popAct = new int[BLcount];
	for (int i=0; i<BLcount-1; i++)
	{
		if ( BurstData[i][NoN+1]>popThreshold ) popAct[i] = 2*popThreshold + (BurstData[i][NoN+1]-popThreshold);	
		else 
		{
			//popAct[i] = popThreshold;
			popAct[i] = 2*BurstData[i][NoN+1];
		}
		NeuronActSum = NeuronActSum + double(popAct[i]);
	}
	
	delete[] popAct;

	/*cout << "process_id: " << process_id << " popThreshold = " << popThreshold 
		<< "  NeuronActSum = " << NeuronActSum << endl;*/
	
	return NeuronActSum;
}


double Network::EvalPopNeuronCC02D()
{
	// to effect ectopic bursts
	/* Sums up longest stretch of neuronal activity;
	during population burst; only excitatory neurons 
	are considered here */

	int popThreshold = int(round(popThresholdFrac*NoNe));
	double NeuronActSum0 = 0;
	double NeuronActSum1 = 0;
	int* popAct;
	popAct = new int[BLcount];
	for (int i=0; i<BLcount-1; i++)
	{
		if ( BurstData[i][NoN+1]>popThreshold )
		{
			popAct[i] = 0.0*popThreshold + 1.0*(BurstData[i][NoN+1]-popThreshold);
			NeuronActSum0 = NeuronActSum0 + double(popAct[i]);
		}
		else 
		{
			popAct[i] = 0;
			if ( NeuronActSum0>NeuronActSum1 )
			{
				NeuronActSum1 = NeuronActSum0;
			}
			NeuronActSum0 = 0;
		}
	}
	
	delete[] popAct;

	/*cout << "process_id: " << process_id << " popThreshold = " << popThreshold 
		<< "  NeuronActSum = " << NeuronActSum << endl;*/
	
	return NeuronActSum1;
}


double Network::EvalPopNeuronCC02A()
{
	/* popThresholdFrac is assumed to be zero/nonzero
	only excitatory neurons are considered here */

	int popThreshold = int(round(popThresholdFrac*NoNe));
	double NeuronActSum = 0;
	int* popAct;
	popAct = new int[BLcount];
	for (int i=0; i<BLcount-1; i++)
	{
		if ( BurstData[i][NoN+1]>=popThreshold )	popAct[i] = BurstData[i][NoN+1]; //used for Result1
		//if ( BurstData[i][NoN+1]>=1 )	popAct[i] = 1;
		else popAct[i] = 0;
		NeuronActSum = NeuronActSum + double(popAct[i]);
	}
	
	delete[] popAct;	
	
	return NeuronActSum;
}


double Network::EvalPopNeuronCC02()
{
	/* same as EvalPopNeuronCC01 except that cross-corelation among
	only excitatory neurons are considered here */

	int popThreshold = int(round(popThresholdFrac*NoNe));
	double popActMean = 0;
	int* popAct;
	popAct = new int[BLcount];
	for (int i=0; i<BLcount-1; i++)
	{
		//if ( BurstData[i][NoN+1]>=popThreshold )	popAct[i] = 1; 
		if ( BurstData[i][NoN+1]>=1 )	popAct[i] = 1;
		else popAct[i] = 0;
		popActMean = popActMean + double(popAct[i]);
	}
	popActMean = popActMean/double(BLcount);

	//double alpha;
	double* CCvector;
	CCvector = new double[NoNe];
	int* CCindex;
	CCindex = new int[NoNe];
	double CCmean;
	double CCi, CCj;
	for (int i=0; i<NoNe; i++)
	{
		CCmean = 0;
		for (int j=0; j<BLcount; j++)
		{
			CCmean = CCmean + BurstData[j][i];
		}
		CCmean = double(CCmean)/double(BLcount);
		if (CCmean > 0.999999) CCmean = 0.999999;
		if (CCmean < 0.000001) CCmean = 0.000001;

		CCvector[i] = 0;
		CCi = 0;
		CCj = 0;
		for (int j=0; j<BLcount; j++)
		{
			//if ((BurstData[j][i]-CCmean)*(popAct[j]-popActMean)<0) alpha = 10.0;
			//else alpha = 0.0;
			//CCvector[i] = CCvector[i] + alpha*(BurstData[j][i]-CCmean)*(popAct[j]-popActMean);
			//CCi = CCi + pow( (BurstData[j][i]-CCmean),2 );
			//CCj = CCj + pow( (popAct[j]-popActMean),2 );
			CCvector[i] = CCvector[i] + (BurstData[j][i]-CCmean)*(popAct[j]-popActMean);
			CCi = CCi + pow( (BurstData[j][i]-CCmean),2 );
			CCj = CCj + pow( (popAct[j]-popActMean),2 );
		}
		CCvector[i] = CCvector[i]/sqrt(CCi*CCj);
		CCindex[i] = i;
	}

	DescendSort(CCvector,CCindex,NoNe);

	//cout << "process_id: " << process_id << "  NoPoPburst = " << NoPoPburst << endl;
	for (int i=0; i<=CCxbins; i++)
	{
		CCx[i] = -1.0 + double(i)*(2.0/double(CCxbins));
	}

	int count = 0;
	for (int i=0; i<=CCxbins; i++)
	{
		if (NoNe-1-count>=0)
		{
			while (CCvector[NoNe-1-count]<=CCx[i])
			{
				count = count + 1;
				if (NoNe-1-count<0) break;
			}
			if (count==0) CCy2[i] = 0;
			else CCy2[i] = double(count)/double(NoNe);			
		}
		else CCy2[i] = CCy2[i-1];
	}

	delete[] CCvector;
	delete[] CCindex;
	delete[] popAct;

	double eval2 = EvalErrorSum2();
	
	return eval2;
}


netconfig_result Network::EvalPopNeuronBurstCC()
{
	netconfig_result nr;

	//ComputePopBurstData(); //uses first NoP neurons 
	//ComputePopBurstData00(); // uses only NoNe neurons
	BurstTimeOnsetSporadic();

	/*PopBurstData[NoEB][NoP]
	NeuronPopBurstData[NoP][PopBurstCountMax]
	for ( int i=0; i<NoP; i++ )
	{
		for ( int j=0; j<NoEB; j++) PopBurstData[j][i] = 0;
	}*/

	for ( int i=0; i<NoP; i++ )
	{
		for ( int j=0; j<NoEB; j++) PopBurstData[j][i] = NeuronPopBurstData[i][j];
	}

	double Error1, Error2, Error3, Error4;
	int popThreshold = int(round(popThresholdFrac*NoNe));
	if (NoPoPburst <= 3)
	{
		Error1 = Error_Max;
		Error2 = pow(PopBurstInterval*160/1000-SetPopBurstFreq,2);
		Error3 = 0;
		Error4 = Error1 + Error2 + Error3;
		cout << "process_id: " << process_id << "  Error1 = " << Error1 
			<< "  Error2 = " << Error2 << "  Error3 = " << Error3 << endl;
		cout << "process_id: " << process_id << "  Error1F = " << Error1/Error4 
			<< "  Error2F = " << Error2/Error4 << "  Error3F = " << Error3/Error4 << endl;

		nr.eval = Error_Max;
		//nr.eval2 = EvalPopNeuronCC0(); // NeuronActSum (uses NoP)
		//nr.eval2 = double(NoNe*NoNe*synGrade-EEsum);
		//nr.eval2 = double(BLcount*NoNe)-SumExcitatorySporadic();
		//nr.eval2 = pow(PopBurstInterval*160/1000-SetPopBurstFreq,2);
		nr.eval2 = double(NoNe*NoNe-EEsum);		
		for (int i=0; i<=CCxbins; i++)
		{
			nr.CCx[i] = CCx[i];
			nr.CCy[i] = CCy[i];
		}		
	}
	else // (NoPoPburst > 3)
	{
		double* CCvector;
		CCvector = new double[NoP];
		int* CCindex;
		CCindex = new int[NoP];
		int CCsum;
		for (int i=0; i<NoP; i++)
		{
			CCsum = 0;
			for (int j=0; j<NoPoPburst; j++)
			{
				CCsum = CCsum + PopBurstData[j][i];
			}
			CCvector[i] = double(CCsum)/double(NoPoPburst);
			CCindex[i] = i;
		}

		DescendSort(CCvector,CCindex,NoP);

		//cout << "process_id: " << process_id << "  NoPoPburst = " << NoPoPburst << endl;

		for (int i=0; i<=CCxbins; i++)
		{
			CCx[i] = 0.0 + double(i)*(1.0/double(CCxbins));
		}

		int count = 0;
		for (int i=0; i<=CCxbins; i++)
		{
			if (NoP-1-count>=0)
			{
				while (CCvector[NoP-1-count]<=CCx[i])
				{
					count = count + 1;
					if (NoP-1-count<0) break;
				}
				if (count==0) CCy[i] = 0;
				else CCy[i] = double(count)/double(NoP);			
			}
			else CCy[i] = CCy[i-1];
			//cout << "CCy[" << i << "] = " << CCy[i] << endl;
		}

		delete[] CCvector;
		delete[] CCindex;


		Error1 = EvalErrorSum();
		Error2 = pow(PopBurstInterval*160/1000-SetPopBurstFreq,2)/10;
		Error3 = 0;
		Error4 = Error1 + Error2 + Error3;
		cout << "process_id: " << process_id << "  Error1 = " << Error1 
			<< "  Error2 = " << Error2 << "  Error3 = " << Error3 << endl;
		cout << "process_id: " << process_id << "  Error1F = " << Error1/Error4 
			<< "  Error2F = " << Error2/Error4 << "  Error3F = " << Error3/Error4 << endl;


		nr.eval = Error1 + Error2;
		//nr.eval = SumExcitatoryPop();//double(NoPoPburst);
		//nr.eval2 = EvalPopNeuronCC0(); // NeuronActSum (uses NoP)
		//nr.eval2 = double(NoNe*NoNe*synGrade-EEsum);
		//nr.eval = EvalPopNeuronCC01();
		//nr.eval2 = double(BLcount*NoNe)-SumExcitatorySporadic();
		//nr.eval2 = pow(PopBurstInterval*160/1000-SetPopBurstFreq,2);
		nr.eval2 = double(NoNe*NoNe-EEsum);
		for (int i=0; i<=CCxbins; i++)
		{
			nr.CCx[i] = CCx[i];
			nr.CCy[i] = CCy[i];
			//nr.CCy[i] = CCy2[i];
		}
	}
	
	return nr;
}


netconfig_result Network::EvalPopNeuronBurstCC01()
{
	netconfig_result nr;

	//ComputePopBurstData();
	//ComputePopBurstData00();
	ComputePopBurstData01();

	if (NoPoPburst <= 3)
	{
		nr.eval = Error_Max;
		nr.eval2 = EvalPopNeuronCC02A(); //NeuronActSum
		for (int i=0; i<=CCxbins; i++)
		{
			nr.CCx[i] = CCx[i];
			nr.CCy[i] = CCy2[i];
		}
	}
	else // (NoPoPburst > 3)
	{		
		//nr.eval = EvalPopNeuronCC01();
		nr.eval = EvalPopNeuronCC02();
		nr.eval2 = EvalPopNeuronCC02A(); //NeuronActSum
		for (int i=0; i<=CCxbins; i++)
		{
			nr.CCx[i] = CCx[i];
			nr.CCy[i] = CCy2[i];
		}
	}
	
	return nr;
}

netconfig_result Network::EvalPopNeuronBurstCC02()
{
	// same a EvalPopNeuronBurstCC01() except that NoP=NoNe
	netconfig_result nr;
	
	//ComputePopBurstData01();
	
	nr.eval = -EvalPopNeuronCC02D();
	nr.eval2 = EvalPopNeuronCC0(); // NeuronActSum
	//nr.eval2 = double(NoNe*NoNe-EEsum);
	//nr.eval2 = double(NoN*NoN*synGrade-EEsum-EIsum-IEsum-IIsum);
	for (int i=0; i<=CCxbins; i++)
	{
		nr.CCx[i] = 0;
		nr.CCy[i] = 0;
	}
	/*cout << "process_id: " << process_id
		<< "  EvalPopNeuronCC02B() = " << EvalPopNeuronCC02B()
			<< "  BLcount = " << BLcount << endl;*/
	
	return nr;
}


netconfig_result Network::EvalCC1()
{
	netconfig_result nr;

	ComputePopBurstData();

	cout << "process_id: " << process_id
			<< "  BLcount = " << BLcount << endl;
	int Mi, Mj;
	double Cii, Cjj, Cij, iMean, jMean;
	double CrossCorr[NoP][NoP];
	for (int i=0; i<NoP; i++)
	{
		for (int j=i; j<NoP; j++)
		{
			Mi = 0;
			for (int k=0; k<BLcount; k++)
			{
				Mi = Mi + BurstData[k][i];
			}
			iMean = double(Mi)/double(BLcount);
			if (iMean > 0.99999) iMean = 0.99999;
			if (iMean < 0.000001) iMean = 0.000001;

			Cii = 0;
			for (int k=0; k<BLcount; k++)
			{
				Cii = Cii + (double(BurstData[k][i])-iMean)*(double(BurstData[k][i])-iMean);
			}

			if (i==j) CrossCorr[i][j] = 1;
			else
			{
				Mj = 0;
				for (int k=0; k<BLcount; k++)
				{
					Mj = Mj + BurstData[k][j];
				}
				jMean = double(Mj)/double(BLcount);
				if (jMean > 0.99999) jMean = 0.99999;
				if (jMean < 0.000001) jMean = 0.000001;

				Cij = 0;
				Cjj = 0;
				for (int k=0; k<BLcount; k++)
				{
					Cjj = Cjj + (double(BurstData[k][j])-jMean)*(double(BurstData[k][j])-jMean);
					Cij = Cij + (double(BurstData[k][i])-iMean)*(double(BurstData[k][j])-jMean);
				}
				CrossCorr[i][j] = Cij/sqrt(Cii*Cjj);
				CrossCorr[j][i] = CrossCorr[i][j];
			}
			//cout << "CrossCorr[" << i << "][" << j << "] = " << CrossCorr[i][j] << endl;
		}
	}


	double* CCvector;
	CCvector = new double[NoP*NoP];
	int* CCindex;
	CCindex = new int[NoP*NoP];
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<NoP; j++)
		{
			CCvector[i*NoP+j] = CrossCorr[i][j];
			CCindex[i*NoP+j] = i*NoP+j;
			//cout << "process_id: " << process_id << "\t CCvector[" << i*NoP+j << "] = " << CCvector[i*NoP+j] << endl;
		}
	}

	DescendSort(CCvector,CCindex,NoP*NoP);

	for (int i=0; i<=CCxbins; i++)
	{
		CCx[i] = -1.0 + double(i)*(2.0/double(CCxbins));
	}

	int count = 0;
	for (int i=0; i<=CCxbins; i++)
	{
		if (NoP*NoP-1-count>=0)
		{
			while (CCvector[NoP*NoP-1-count]<=CCx[i])
			{
				count = count + 1;
				if (NoP*NoP-1-count<0) break;
			}
			if (count==0) CCy[i] = 0;
			else CCy[i] = double(count)/double(NoP*NoP);
			//cout << "CCy[" << i << "] = " << CCy[i] << endl;
		}
		else CCy[i] = CCy[i-1];
	}

	delete[] CCvector;
	delete[] CCindex;

	if ( NoPoPburst<3 )
	{
		nr.eval = Error_Max;
	}
	else nr.eval = EvalErrorSum();
	nr.eval2 = EvalPopNeuronCC0();
	for (int i=0; i<=CCxbins; i++)
	{
		nr.CCx[i] = CCx[i];
		nr.CCy[i] = CCy[i];
	}

	return nr;
}



void Network::ReadConMat()
{
	ostringstream stream;
	//stream << "FinalConMatA.txt";
	//stream << "FinalConMatB.txt";
	stream << "FinalConMatA" << job_id << ".txt";
	//char filename0[] = "FinalConMat0.txt";
	ifstream inFile1;
	inFile1.open(stream.str().c_str());
	if (!inFile1)
	{
		cout << "Unable to open file FinalConMat.txt" << endl;
		//exit(1); // terminate with error
	}
	int dummy;
	for (int i=0; i<NoN; i++)
	{
		for (int j=0; j<NoN; j++)
		{
			if ( i<NoNe )
			{
				if ( j<NoNe )
				{
					inFile1 >> dummy;
					conMatEE[i][j] = dummy;
				}
				else
				{
					inFile1 >> dummy;
					conMatEI[i][j-NoNe] = dummy;
				}
			}
			else
			{
				if ( j<NoNe )
				{
					inFile1 >> dummy;
					conMatIE[i-NoNe][j] = dummy;
				}
				else
				{
					inFile1 >> dummy;
					conMatII[i-NoNe][j-NoNe] = dummy;
				}
			}
		}
	}
	inFile1.close();
	//StoreInitialConMat();
}

void Network::CopyEE()
{
	int EErow = 42; 
	int EEcol = 42;
	ostringstream stream;
	stream << "FinalConMatA.txt";
	//stream << "FinalConMatB.txt";
	//stream << "FinalConMatA" << process_id << ".txt";
	//char filename0[] = "FinalConMat0.txt";
	ifstream inFile1;
	inFile1.open(stream.str().c_str());
	if (!inFile1)
	{
		cout << "Unable to open file FinalConMat.txt" << endl;
		//exit(1); // terminate with error
	}
	int dummy;
	for (int i=0; i<EErow; i++)
	{
		for (int j=0; j<EEcol; j++)
		{
			if ( i<NoNe )
			{
				if ( j<NoNe )
				{
					inFile1 >> dummy;
					conMatEE[i][j] = dummy;
				}
				else
				{
					inFile1 >> dummy;
					//conMatEI[i][j-NoNe] = dummy;
				}
			}
			else
			{
				if ( j<NoNe )
				{
					inFile1 >> dummy;
					//conMatIE[i-NoNe][j] = dummy;
				}
				else
				{
					inFile1 >> dummy;
					//conMatII[i-NoNe][j-NoNe] = dummy;
				}
			}
		}
	}
	inFile1.close();
	//StoreInitialConMat();
}


void Network::ReadNeuronProp()
{
	ostringstream stream;
	stream << "NeuronProp82A" << job_id << ".txt";
	//stream << "examplePP1.txt";
	//stream << "NeuronProp60e.txt";
	//stream << "NeuronProp80c.txt";
	//stream << "Neuron40Shrink.txt";
	//stream << "NeuronPropUni100.txt";
	//stream << "NeuronProp82.txt";
	//stream << "NeuronProp160.txt";
	ifstream inFile1;
	inFile1.open(stream.str().c_str());
	if (!inFile1)
	{
		cout << "Unable to open file NeuronProp.txt" << endl;
		//exit(1); // terminate with error
	}
	double dummy;
	for (int i=0; i<NoN; i++)
	{	
		if ( i<NoNe )
		{
			inFile1 >> dummy;
			gLeake[i] = dummy;
			inFile1 >> dummy;
			gNaPe[i] = dummy;
		}
		else
		{
			inFile1 >> dummy;
			gLeaki[i-NoNe] = dummy;
			inFile1 >> dummy;
			gNaPi[i-NoNe] = dummy;
		}
	}
	inFile1.close();
}


void Network::DoubleSort(double E1double[], double E2[], int Index[], int IndexSorted[], int NoG, int shift) {
	//int E1[NoG], E1Copy[NoG], Index[NoG], IndexSorted[NoG];	
	//double E2[NoG];	
	
	int* E1;
	E1 = new int[NoG];	
	for (int i=0; i<NoG; i++)	
	{		
		E1[i] = round( 100000*double(E1double[i]) );
	}

			
	/*for (int i=0; i<NoG; i++)	
	{		
		cout << Index[i] << "\t" << E1[i] << "\t" << E2[i] << endl;	
	}*/	
		
	DescendSortInteger(E1, Index, NoG);	
		
	/*for (int i=0; i<NoG; i++)	
	{		
		cout << Index[i] << "\t" << E1[i] << "\t" << E2[Index[i]-shift] << endl;	
	}*/	
		
	int Count = 0;		
	int CountLast = 0;	
	//int Es[NoG], En[NoG]; // s: same E1 value; n: no. of same E1 values	
	int* Es; // s: same EvalArray value	
	Es = new int[NoG];	
	int* En; // n: no. of same EvalArray values	
	En = new int[NoG];	
	int Eg = E1[Count];	
	for (int i=0; i<NoG; i++)	
	{		
		if (Count<NoG)		
		{			
			if ( Eg != E1[i] )			
			{				
				En[Count] = i - CountLast;				
				Es[Count] = Eg;				
				Count = Count + 1;				
				CountLast = i;				
				Eg = E1[i];			
			}
			
			if (i==NoG-1)			
			{				
				En[Count] = i - CountLast + 1;	
				Es[Count] = Eg;			
			}		
		}	
	}	
	
		
	/*for (int i=0; i<=Count; i++)	
	{		
		cout << Es[i] << "\t" << En[i] << endl;	
	}*/
			
	double* E2sub;	
	int* E2subIndex;	
	int counter = 0;	
	for (int i=0; i<=Count; i++)	
	{		
		//cout << "\n En[" << i << "]=" << En[i] << endl;
		E2sub = new double[En[i]];		
		E2subIndex = new int[En[i]];		
		for (int j=0; j<En[i]; j++)		
		{			
			E2sub[j] = E2[Index[counter+j]-shift];
			E2subIndex[j] = Index[counter+j];
			//cout << E2sub[j] << "\t" << E2subIndex[j] << endl;
		}		
			
		DescendSort(E2sub, E2subIndex, En[i]);		
			
		for (int j=0; j<En[i]; j++)		
		{			
			//IndexSorted[counter+j] = E2subIndex[j];
			IndexSorted[counter+j] = E2subIndex[En[i]-1-j];
			//cout << counter+j << "\t" << E2subIndex[j] << "\t" << E2[IndexSorted[counter+j]] << endl;
		}
		
		counter = counter + En[i];	
	}
	
	delete [] E2sub;	
	delete [] E2subIndex;	
	delete [] Es;	
	delete [] En;	
	delete [] E1;	
	
	/*cout << "\n \n " << "Final array" << endl;	
	for (int i=0; i<NoG; i++)	
	{		
		cout << IndexSorted[i] << "\t" << E2[IndexSorted[i]-shift] << endl;
	}*/
}



void Network::DescendSort(double RanMat[], int RanIndex[], int NoSyn)
{
	int flipCount;
	do
	{
		flipCount = 0;
		for (int i=0; i<NoSyn-1; i++)
		{
			double dumN;
			int dumIn;
			if (RanMat[i]<RanMat[i+1])
			{
				dumN = RanMat[i];
				RanMat[i] = RanMat[i+1];
				RanMat[i+1] = dumN;
				dumIn = RanIndex[i];
				RanIndex[i] = RanIndex[i+1];
				RanIndex[i+1] = dumIn;
				flipCount = flipCount + 1;
			}
		}
	} while (flipCount > 0);
}

void Network::DescendSortInteger(int RanMat[], int RanIndex[], int NoSyn)
{
	int flipCount;
	do
	{
		flipCount = 0;
		for (int i=0; i<NoSyn-1; i++)
		{
			int dumN;
			int dumIn;
			if (RanMat[i]<RanMat[i+1])
			{
				dumN = RanMat[i];
				RanMat[i] = RanMat[i+1];
				RanMat[i+1] = dumN;
				dumIn = RanIndex[i];
				RanIndex[i] = RanIndex[i+1];
				RanIndex[i+1] = dumIn;
				flipCount = flipCount + 1;
			}
		}
	} while (flipCount > 0);
}


void Network::AllocatePMprop (double gNaP[], double gLeak[], int N)
{
	double gNaP_mean, gNaP_std, gLeak_mean, gLeak_std;
	gNaP_mean = 0.84;
	gNaP_std = 107.0/100*gNaP_mean;
	gLeak_mean = 6.84;
	gLeak_std = 25.0/100*gLeak_mean;

	//srand ( time(NULL)+process_id );	
	double U1, U2, U3, U4;
	double gNaP_temp, gLeak_temp;
	int n=0;
	while ( n<N )
	{
		U1 = (double) rand()/RAND_MAX;
		U2 = (double) rand()/RAND_MAX;
		U3 = (double) rand()/RAND_MAX;
		U4 = (double) rand()/RAND_MAX;
		gNaP_temp = gNaP_mean + gNaP_std*sqrt(-2*log(U1))*cos(2*PI*U2);
		gLeak_temp = gLeak_mean + gLeak_std*sqrt(-2*log(U3))*cos(2*PI*U4);
		if ( (4.3*gLeak_temp-0.9*gNaP_temp-4.3*0.4+0.7*0.9)>0 && (3.8*gLeak_temp-4.7*gNaP_temp+0.2*0.9)<0 )
		{
			gNaP[n] = gNaP_temp;
			gLeak[n] = gLeak_temp;		
			n = n + 1;
		} 		
	}

}


void Network::AllocateNPMprop (double gNaP[], double gLeak[], int N)
{
	double gNaP_mean, gNaP_std, gLeak_mean, gLeak_std;
	gNaP_mean = 1.13;
	gNaP_std = 26.0/100*gNaP_mean;
	gLeak_mean = 1.27;
	gLeak_std = 107.0/100*gLeak_mean;

	//srand ( time(NULL)+100.0*process_id );
	double U1, U2, U3, U4;
	double gNaP_temp, gLeak_temp;
	int n=0;
	while ( n<N )
	{
		U1 = (double) rand()/RAND_MAX;
		U2 = (double) rand()/RAND_MAX;
		U3 = (double) rand()/RAND_MAX;
		U4 = (double) rand()/RAND_MAX;
		gNaP_temp = gNaP_mean + gNaP_std*sqrt(-2*log(U1))*cos(2*PI*U2);
		gLeak_temp = gLeak_mean + gLeak_std*sqrt(-2*log(U3))*cos(2*PI*U4);
		if ( gNaP_temp>0 && (3.4*gLeak_temp-4.7*gNaP_temp-0.3*3.4)>0 )
		{
			gNaP[n] = gNaP_temp;
			gLeak[n] = gLeak_temp;		
			n = n + 1;
		} 		
	}  	

}


void Network::InitializeNeuronProp()
{
	int NoNPMe = NoNe - NoPMe; // no. of excitatory NPM neurons
	int NoNPMi = NoNi - NoPMi; // no. of inhibitory NPM neurons
	
	double gNaP1[NoPMe], gLeak1[NoPMe];
	//double gNaP2[NoNPMe], gLeak2[NoNPMe];
	double gNaP3[NoPMi], gLeak3[NoPMi];
	//double gNaP4[NoNPMi], gLeak4[NoNPMi];

	double* gNaP2;
	gNaP2 = new double[NoNPMe];
	double* gLeak2;
	gLeak2 = new double[NoNPMe];
	double* gNaP4;
	gNaP4 = new double[NoNPMi];
	double* gLeak4;
	gLeak4 = new double[NoNPMi];

	AllocatePMprop(gNaP1, gLeak1, NoPMe);
	AllocateNPMprop(gNaP2, gLeak2, NoNPMe);
	AllocatePMprop(gNaP3, gLeak3, NoPMi);
	AllocateNPMprop(gNaP4, gLeak4, NoNPMi);
	
	for ( int i=0 ; i<NoNe ; i++ )
	{
		if ( i<NoPMe )	
		{
			gNaPe[i] = gNaP1[i];
			gLeake[i] = gLeak1[i];
		}
		else
		{
			gNaPe[i] = gNaP2[i-NoPMe];
			gLeake[i] = gLeak2[i-NoPMe];
		}
	}
	
	for ( int i=0 ; i<NoNi ; i++ )
	{
		if ( i<NoPMi )	
		{
			gNaPi[i] = gNaP3[i];
			gLeaki[i] = gLeak3[i];
		}
		else
		{
			gNaPi[i] = gNaP4[i-NoPMi];
			gLeaki[i] = gLeak4[i-NoPMi];
		}
	}
	
	delete [] gNaP2;
	delete [] gLeak2;
	delete [] gNaP4;
	delete [] gLeak4;

}


void Network::NetworkTopology(int MatTag)
{
	int M, N;
	if ( MatTag==1 || MatTag==2 ) M = NoNe; 
	else M = NoNi;

	if ( MatTag==1 || MatTag==3 ) N = NoNe; 
	else N = NoNi;

	double SynFrac0;
	if (MatTag==1) SynFrac0 = SynFrac;
	if (MatTag==2) SynFrac0 = SynFracEI;
	if (MatTag==3) SynFrac0 = SynFracIE;
	if (MatTag==4) SynFrac0 = SynFracII;

	int NoSyn = int( round(SynFrac0*double(M*N)) );	
	int RN1_copy = 0;
	int RN2_copy = 0;

	double* RanMat;
	RanMat = new double[M*N];
	double* RanMatCopy;
	RanMatCopy = new double[M*N];
	int* RanIndex;
	RanIndex = new int[M*N];
	double RN1;
	//srand(time(NULL)+200.0*process_id);
	for (int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			do
			{
				RN1 = (double)rand( ) / RAND_MAX;
				RN2_copy = (int) (RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;
			RanMat[i*N+j] = RN1;
			RanMatCopy[i*N+j] = RN1;
			RanIndex[i*N+j] = i*N+j;			
		}
	}

	DescendSort(RanMat,RanIndex,M*N);
	

	double threshold;
	int thresholdIndex = NoSyn-1;
	if ( thresholdIndex < 0 ) threshold = 1.0;
	else threshold = RanMat[thresholdIndex];
	delete [] RanMat;
	delete [] RanIndex;
	
	if ( MatTag == 1 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) conMatEE[i][j] = 1;
				else conMatEE[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 2 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) conMatEI[i][j] = 1;
				else conMatEI[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 3 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) conMatIE[i][j] = 1;
				else conMatIE[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 4 )
	{			
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) conMatII[i][j] = 1;
				else conMatII[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}


	//ConnectIsolatedNeurons(conMat, M);
	//StoreInitialConMat();
}

void Network::InitializeNetworkTopology()
{
	for ( int i=1 ; i<=4; i++ )
	{		
		NetworkTopology(i);
	}
}


netconfig_cond Network::Initialize(int z) {
	
	srand ( time(NULL)+process_id );
	netconfig_cond nc;
	
	if ( z==0 )
	{
		NetTag = process_id;
		InitializeNeuronProp();
		InitializeNetworkTopology();
	}
	
	if ( z==1 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		InitializeNetworkTopology();
	}

	if ( z==2 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		ReadConMat();
	}

	if ( z==3 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		InitializeNetworkTopology();
		CopyEE();
	}

	
	for (int i=0; i<NoN; i++)
	{
		if ( i<NoNe )
		{
			nc.GNaP[i] = gNaPe[i];
			nc.GLeak[i] = gLeake[i];
		}
		else
		{
			nc.GNaP[i] = gNaPi[i-NoNe];
			nc.GLeak[i] = gLeaki[i-NoNe];
		}		
	}

	
	for (int m = 0; m < NoN; m++)
	{
		//gtonic_e[m] = tonicConductance;
		for(int mm = 0; mm < NoN; mm++)
		{
			gsyn_e[m][mm] = synS; // index m is for pre-synaptic neuron
		}
		cp[m].Gtonic_e(gtonic_e[m]);
		if ( m<NoNe )
		{
			cp[m].GNaPo_set(gNaPe[m]);
			cp[m].GLeako_set(gLeake[m]);
			gtonic_e[m] = tonicConductance;
		}
		else
		{
			cp[m].GNaPo_set(gNaPi[m-NoNe]);
			cp[m].GLeako_set(gLeaki[m-NoNe]);
			gtonic_e[m] = tonicConductanceI;
		}	
		
		
	}

	
	for (int i=0; i<NoN; i++)
	{
		for (int j=0; j<NoN; j++)
		{
			if ( i<NoNe )
			{
				if ( j<NoNe )	nc.conMat[i][j] = conMatEE[i][j];
				else nc.conMat[i][j] = conMatEI[i][j-NoNe];
			}
			else
			{
				if ( j<NoNe )	nc.conMat[i][j] = conMatIE[i-NoNe][j];
				else nc.conMat[i][j] = conMatII[i-NoNe][j-NoNe];
			}
		}
	}
	
	nc.NetTag = NetTag;

	return nc;
}



netconfig_cond Network::InitializeSmallWorldEE(int z, int n0, double q) {
	
	srand ( time(NULL)+process_id );
	netconfig_cond nc;
	
	if ( z==0 )
	{
		NetTag = process_id;
		InitializeNeuronProp();
		InitializeNetworkTopology(); // initializes conMatEE, conMatEI...
		SmallWorldEE( n0, q); // initializes only conMatEE
		
	}
	
	if ( z==1 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		InitializeNetworkTopology();// initializes conMatEE, conMatEI...
		SmallWorldEE( n0, q); // initializes only conMatEE
	}

	if ( z==2 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		ReadConMat();
	}

	
	for (int i=0; i<NoN; i++)
	{
		if ( i<NoNe )
		{
			nc.GNaP[i] = gNaPe[i];
			nc.GLeak[i] = gLeake[i];
		}
		else
		{
			nc.GNaP[i] = gNaPi[i-NoNe];
			nc.GLeak[i] = gLeaki[i-NoNe];
		}
	}

	
	for (int m = 0; m < NoN; m++)
	{
		//gtonic_e[m] = tonicConductance;
		for(int mm = 0; mm < NoN; mm++)
		{
			gsyn_e[m][mm] = synS; // index m is for pre-synaptic neuron
		}
		cp[m].Gtonic_e(gtonic_e[m]);
		if ( m<NoNe )
		{
			cp[m].GNaPo_set(gNaPe[m]);
			cp[m].GLeako_set(gLeake[m]);
			gtonic_e[m] = tonicConductance;
		}
		else
		{
			cp[m].GNaPo_set(gNaPi[m-NoNe]);
			cp[m].GLeako_set(gLeaki[m-NoNe]);
			gtonic_e[m] = tonicConductanceI;
		}	
		
		
	}

	
	for (int i=0; i<NoN; i++)
	{
		for (int j=0; j<NoN; j++)
		{
			if ( i<NoNe )
			{
				if ( j<NoNe )	nc.conMat[i][j] = conMatEE[i][j];
				else nc.conMat[i][j] = conMatEI[i][j-NoNe];
			}
			else
			{
				if ( j<NoNe )	nc.conMat[i][j] = conMatIE[i-NoNe][j];
				else nc.conMat[i][j] = conMatII[i-NoNe][j-NoNe];
			}
		}
	}
	
	nc.NetTag = NetTag;

	return nc;
}


void Network::Integrate() {
	dt = 0.05; //default 0.05(ms), for variable step methods, set 0.5(ms).
	double sum_gsyn_eXs; // for excitatory syns
	double sum_gsyn_iXs; // for inhibitory syns
		
	for (int iN = 0; iN < NoN; iN++)
	{
		sum_gsyn_eXs = 0; 
		sum_gsyn_iXs = 0;
		for (int jN = 0; jN < NoN; jN++)
		{
			if ( jN<NoNe )
			{	
				if ( iN<NoNe ) sum_gsyn_eXs = sum_gsyn_eXs + synS*conMatEE[jN][iN]*s[jN][0];
				else sum_gsyn_eXs = sum_gsyn_eXs + synS*conMatEI[jN][iN-NoNe]*s[jN][0];
			}
			else
			{
				if ( iN<NoNe ) sum_gsyn_iXs = sum_gsyn_iXs + synSinh*conMatIE[jN-NoNe][iN]*s[jN][0];
				else sum_gsyn_iXs = sum_gsyn_iXs + synSinh*conMatII[jN-NoNe][iN-NoNe]*s[jN][0];
			}
		}	
		cp[iN].RK4(V[iN][0], n[iN][0], h[iN][0], s[iN][0], sum_gsyn_eXs, sum_gsyn_iXs, dt);
		V[iN][1] = cp[iN].V_Next();
		n[iN][1] = cp[iN].n_Next();
		h[iN][1] = cp[iN].h_Next();
		s[iN][1] = cp[iN].s_Next();
	}

	for (int k = 0; k < NoN; k++)
	{
		V[k][0] = V[k][1];
		n[k][0] = n[k][1];
		h[k][0] = h[k][1];
		s[k][0] = s[k][1];
	}
}


void Network::Execute(double T, double t_dispStep,
		double t_saveStep, double Tws, const char *filename)
{
	//synSinh = 0.006*(process_id-1);
	//cout << "process_id = " << process_id << "\t synSinh = " << synSinh << endl;
	//===============================
	for (int i=0; i<binLength; i++)
	{
		t[i] = 0;
		for (int j=0; j<=NoN+1; j++)
		{
			BurstData[i][j] = 0;
		}
	}
	//===============================

	//===============================
	for (int m = 0; m < NoN; m++)
	{
		V[m][0] = -61.5; //mV
		n[m][0] = 0;
		h[m][0] = 0;
		s[m][0] = 0;
	}
	//===============================

	double time = 0;
	double t_disp = 0.0;
	int length = int((T-Tws)/t_saveStep);
	int interSpikeThreshold = 160; //ms
	double burstThreshold = -20; //mV
	BLcount = 0;
	int RowCount = 0;
	int kcount = 0;

	ofstream myfile2;
	myfile2.open (filename); // time & Bursting_Neuron_#
	ofstream myfile01p;
	ostringstream stream01;
	stream01 << "NeuronBurstData" << file_id << ".txt"; // time & Vsum/NoP
	myfile01p.open (stream01.str().c_str());
	ofstream myfile02p;
	ostringstream stream02;
	stream02 << "ExcitatoryBurstData" << file_id << ".txt"; // time & Vsum/NoNe
	myfile02p.open (stream02.str().c_str());
	double Vsum;
	while ( time <= T )
	{
		// for display on console window
		if ( time>=t_disp )
		{
			cout << "process_id: " << process_id
					<< "  t = " << time << endl;
			t_disp = t_disp + t_dispStep;
		}

		Integrate();
		time = time + dt;


		if ( time>=Tws && int(round(time*100))%int(round(t_saveStep*100))==0 )
		{
			Vsum = 0;
			for (int iN = 0; iN < NoN; iN++)
			{
				
				Vsum = Vsum + V_now(iN);
				/*if (iN == 0) myfile01p << time << "\t" << V_now(iN);
				if (iN>0 && iN<NoN-1) myfile01p << "\t" << V_now(iN);
				if (iN == NoN-1) myfile01p << "\t" << V_now(iN)
							<< "\t" << (double) Vsum/NoP << endl;*/
				if (iN == NoN-1) myfile01p << time << "\t" << (double) Vsum/NoN << endl;
				if (iN == NoNe-1) myfile02p << time << "\t" << (double) Vsum/NoNe << endl;
				
				if (V_now(iN)>burstThreshold)
				{
					myfile2 << time << "\t" << iN << endl;
				}
			}



			if ( RowCount<(length-interSpikeThreshold) )
			{
				if ( kcount==round(double(interSpikeThreshold)/2) ) t[BLcount] = time;
				for (int i=0; i<NoN; i++)
				{
					if (BurstData[BLcount][i]==0)
					{
						if (V_now(i)>burstThreshold) BurstData[BLcount][i] = 1;
					}
				}

				kcount = kcount + 1;
				if (kcount==interSpikeThreshold)
				{
					kcount = 0;
					BLcount = BLcount + 1;
				}
			}
			RowCount = RowCount + 1;
		}
	}
	myfile01p.close();
	myfile02p.close();
	myfile2.close();


	int Bsum;
	int BsumE; // accounting only excitatory neurons
	for (int i=0; i<BLcount-1; i++)
	{
		Bsum = 0;
		BsumE = 0;
		for (int j=0; j<NoN; j++)
		{ 
			Bsum = Bsum + BurstData[i][j];
			if ( j<NoNe ) BsumE = BsumE + BurstData[i][j];
		}
		BurstData[i][NoN] = Bsum;
		BurstData[i][NoN+1] = BsumE;
	}

	ofstream myfile0p;
	//myfile0p.open ("BurstData.txt");
	ostringstream stream1;
	stream1 << "BurstData" << file_id << ".txt"; /* time & raster_data_of_(NoP)_neurons & SUM_of_raster_data_of_(NoP)_neurons & SUM_of_raster_data_of_(NoNe)_neurons */
	myfile0p.open (stream1.str().c_str());
	for (int i=0; i<binLength; i++)
	{
		myfile0p <<  t[i]  << "\t" ;
		for (int j=0; j<NoN+2; j++)
		{
			if (j==NoN+1) myfile0p << BurstData[i][j] << endl;
			else myfile0p <<  BurstData[i][j] << "\t";
		}
	}
	myfile0p.close();
}



void Network::ExecuteNoWrite(double T, double t_dispStep,
		double t_saveStep, double Tws)
{
	//===============================
	for (int i=0; i<binLength; i++)
	{
		t[i] = 0;
		for (int j=0; j<=NoN+1; j++)
		{
			BurstData[i][j] = 0;
		}
	}
	//===============================

	//===============================
	for (int m = 0; m < NoN; m++)
	{
		V[m][0] = -61.5; //mV
		n[m][0] = 0;
		h[m][0] = 0;
		s[m][0] = 0;
	}
	//===============================

	double time = 0;
	double t_disp = 0.0;
	int length = int((T-Tws)/t_saveStep);
	int interSpikeThreshold = 160; //ms
	double burstThreshold = -20; //mV
	BLcount = 0;
	int RowCount = 0;
	int kcount = 0;
	
	while ( time <= T )
	{
		// for display on console window
		if ( time>=t_disp )
		{
			cout << "process_id: " << process_id
					<< "  t = " << time << endl;
			t_disp = t_disp + t_dispStep;
		}

		Integrate();
		time = time + dt;

		if ( time>=Tws && int(round(time*100))%int(round(t_saveStep*100))==0 )
		{
			if ( RowCount<(length-interSpikeThreshold) )
			{
				if ( kcount==round(double(interSpikeThreshold)/2) ) t[BLcount] = time;
				for (int i=0; i<NoN; i++)
				{
					if (BurstData[BLcount][i]==0)
					{
						if (V_now(i)>burstThreshold) BurstData[BLcount][i] = 1;
					}
				}

				kcount = kcount + 1;
				if (kcount==interSpikeThreshold)
				{
					kcount = 0;
					BLcount = BLcount + 1;
				}
			}
			RowCount = RowCount + 1;
		}
	}
	

	int Bsum;
	int BsumE; // accounting only excitatory neurons
	for (int i=0; i<BLcount-1; i++)
	{
		Bsum = 0;
		BsumE = 0;
		for (int j=0; j<NoN; j++)
		{ 
			Bsum = Bsum + BurstData[i][j];
			if ( j<NoNe ) BsumE = BsumE + BurstData[i][j];
		}
		BurstData[i][NoN] = Bsum;
		BurstData[i][NoN+1] = BsumE;
	}

	
}


double Network::V_now( int icell ) {
	return V[icell][0];
}


void Network::CopyNetworkCondition(netconfig_cond nc, int z) 
{
	for (int i=0; i<NoN; i++)
	{
		if ( i<NoNe )
		{
			gNaPe[i] = nc.GNaP[i];
			gLeake[i] = nc.GLeak[i];
		}
		else
		{
			gNaPi[i-NoNe] = nc.GNaP[i];
			gLeaki[i-NoNe] = nc.GLeak[i];
		}
	}

	cout << "process_id : " << process_id << " z = " << z << endl;
	if (z!=0)
	{
		NetTag = nc.NetTag;

		for (int i=0; i<NoN; i++)
		{
			for (int j=0; j<NoN; j++)
			{
				if ( i<NoNe )
				{
					if ( j<NoNe )	conMatEE[i][j] = nc.conMat[i][j];
					else conMatEI[i][j-NoNe] = nc.conMat[i][j];
				}
				else
				{
					if ( j<NoNe )	conMatIE[i-NoNe][j] = nc.conMat[i][j];
					else conMatII[i-NoNe][j-NoNe] = nc.conMat[i][j];
				}
			}
		}
	}
	

}



void Network::ConnectMutationSeg(int iter)
{
	/* pick 4 neurons as centers and define their neighbors (N)
	 * delete connections from neurons with center N[0] to neurons with center N[1]
	 * Add connections from neurons with center N[2] to neurons with center N[3] */	
	
	srand ( time(NULL)+process_id+iter*1000);

	//int NoN = NoNe + NoNi;
	int* neuron;
	neuron = new int[4];
	int* N;
	N = new int[4]; // no. of neurons in neighborhood of neuron[4]	
	int* Ng;
	Ng = new int[4]; // no. of neurons in network containing neuron[4]
	int NeighborMaxE = round(0.2*double(NoNe));
	int NeighborMaxI = round(0.2*double(NoNi));
	int RN1_copy = 0;
	int RN2_copy = 0;
	
	double RN1;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	if (RN1<0.5) // small mutation/neighbor
	{
		NeighborMaxE = round(RN1*0.2*double(NoNe));
		NeighborMaxI = round(RN1*0.2*double(NoNi));
	}
	// else // large mutation/neighbor

	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if (RN1>0.5)
		{
			neuron[2*k] = floor(double(NoNe)*RN1);
			//neuron[2*k] = floor(double(NoN)*RN1);
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			//neuron[2*k+1] = NoNe + floor(double(NoNi)*RN1);
			neuron[2*k+1] = floor(double(NoNe)*RN1);
			/*if ( neuron[2*k]<NoNe ) neuron[2*k+1] = floor(double(NoNe)*RN1);
			else neuron[2*k+1] = NoNe + floor(double(NoNi)*RN1);*/
		}
		else
		{
			//neuron[2*k] = NoNe + floor(double(NoNi)*RN1);
			neuron[2*k] = floor(double(NoNe)*RN1);
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			neuron[2*k+1] = floor(double(NoNe)*RN1);
			/*if ( neuron[2*k]<NoNe ) neuron[2*k+1] = floor(double(NoNe)*RN1);
			else neuron[2*k+1] = NoNe + floor(double(NoNi)*RN1);*/
		}	
		
	}

	/* selecting neurons so that during each mutation
	    addition and deletion of synapses happen in the same 
		sub matrix [ EE, EI, IE or II ]*/
	/*for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoN)*RN1);
	}

	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if ( neuron[k]<NoNe ) neuron[2+k] = floor(double(NoNe)*RN1);
			else neuron[2+k] = NoNe + floor(double(NoNi)*RN1);
	}*/


	for (int k=0; k<4; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if ( neuron[k]<NoNe )	
		{
			N[k] = ceil(double(NeighborMaxE)*RN1);
			Ng[k] = NoNe;
		}
		else 
		{
			N[k] = ceil(double(NeighborMaxI)*RN1);
			Ng[k] = NoNi;
		}
	}


	// choosing nature of mutation by redefining neuron[k] and N[k]
	// *************************************************************
	/*do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	if ( RN1<0.2 ) // flip local inSyn to outSyn 
	{
		neuron[2] = neuron[1];
		neuron[3] = neuron[0];
		N[2] = N[1];
		N[3] = N[0];
		Ng[2] = Ng[1];
		Ng[3] = Ng[0];
	}
	if ( RN1>=0.2 && RN1<0.4 ) // displace outSyn from neighborhood of neuron[0] 
	{
		if (neuron[0]<NoNe)
		{
			if (neuron[2]<NoNe)
			{
				neuron[2] = neuron[0];
				N[2] = N[0];
				Ng[2] = Ng[0];
			}
			else
			{
				neuron[2] = neuron[1];
				N[2] = N[1];
				Ng[2] = Ng[1];
			}
		}
		else // if (neuron[0]>NoNe)
		{
			if (neuron[2]>NoNe)
			{
				neuron[2] = neuron[0];
				N[2] = N[0];
				Ng[2] = Ng[0];
			}
			else
			{
				neuron[2] = neuron[1];
				N[2] = N[1];
				Ng[2] = Ng[1];
			}			
		}
	}
	if ( RN1>=0.4 && RN1<0.6 ) // displace inSyn from neighborhood of neuron[1] 
	{
		if (neuron[1]<NoNe)
		{
			if (neuron[3]<NoNe)
			{
				neuron[3] = neuron[1];
				N[3] = N[1];
				Ng[3] = Ng[1];
			}
			else
			{
				neuron[3] = neuron[0];
				N[3] = N[0];
				Ng[3] = Ng[0];
			}			
		}
		else // if (neuron[1]>NoNe)
		{
			if (neuron[3]>NoNe)
			{
				neuron[3] = neuron[1];
				N[3] = N[1];
				Ng[3] = Ng[1];
			}
			else
			{
				neuron[3] = neuron[0];
				N[3] = N[0];
				Ng[3] = Ng[0];
			}			
		}
	}	*/	
	// if RN1>=0.6, then usual mutation comences
	// **************************************************************** 

	
	for (int k=0; k<4; k++)
	{
		cout << "process_id: " << process_id << "\t neuron[" << k << "] = " << neuron[k] << "\t N[" << k << "] = " <<  N[k] << "\t Ng[" << k << "] = " << Ng[k] << endl;		
	}

	double* dist0;
	int* dist0index;
	dist0 = new double[Ng[0]];	
	dist0index = new int[Ng[0]];

	double* dist1;
	int* dist1index;
	dist1 = new double[Ng[1]];	
	dist1index = new int[Ng[1]];

	double* dist2;
	int* dist2index;
	dist2 = new double[Ng[2]];	
	dist2index = new int[Ng[2]];

	double* dist3;
	int* dist3index;
	dist3 = new double[Ng[3]];	
	dist3index = new int[Ng[3]];

		
	for (int k=0; k<Ng[0]; k++)
	{
		if (neuron[0]<NoNe) dist0[k] = pow( (gLeake[k]-gLeake[neuron[0]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[0]]),2 );
		else dist0[k] = pow( (gLeaki[k]-gLeaki[neuron[0]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[0]-NoNe]),2 );
		dist0index[k] = k;
	}
	

	for (int k=0; k<Ng[1]; k++)
	{
		if (neuron[1]<NoNe) dist1[k] = pow( (gLeake[k]-gLeake[neuron[1]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[1]]),2 );
		else dist1[k] = pow( (gLeaki[k]-gLeaki[neuron[1]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[1]-NoNe]),2 );
		dist1index[k] = k;		
	}
	

	for (int k=0; k<Ng[2]; k++)
	{
		if (neuron[2]<NoNe) dist2[k] = pow( (gLeake[k]-gLeake[neuron[2]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[2]]),2 );
		else dist2[k] = pow( (gLeaki[k]-gLeaki[neuron[2]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[2]-NoNe]),2 );
		dist2index[k] = k;
	}
	

	for (int k=0; k<Ng[3]; k++)
	{
		if (neuron[3]<NoNe) dist3[k] = pow( (gLeake[k]-gLeake[neuron[3]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[3]]),2 );
		else dist3[k] = pow( (gLeaki[k]-gLeaki[neuron[3]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[3]-NoNe]),2 );
		dist3index[k] = k;			
	}
	

	DescendSort(dist0, dist0index, Ng[0]);
	DescendSort(dist1, dist1index, Ng[1]);
	DescendSort(dist2, dist2index, Ng[2]);
	DescendSort(dist3, dist3index, Ng[3]);
	

	int SynDeletePossible = 0;
	double* SynDeleteArray;
	SynDeleteArray = new double[N[0]*N[1]];
	double* SynDeleteArrayCopy;
	SynDeleteArrayCopy = new double[N[0]*N[1]];
	int* SynDeleteArrayIndex;
	SynDeleteArrayIndex = new int[N[0]*N[1]];
	int countD = 0;	
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (conMatII[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]]==1)
					{
						SynDeletePossible = SynDeletePossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynDeleteArray[countD] = RN1;
						SynDeleteArrayCopy[countD] = RN1;
						SynDeleteArrayIndex[countD] = countD;
					}
					else
					{
						SynDeleteArray[countD] = 0;
						SynDeleteArrayCopy[countD] = 0;
						SynDeleteArrayIndex[countD] = countD;
					}					
					countD = countD + 1;
				}
			}			
		}
		else 
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (conMatIE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]]==1)
					{
						SynDeletePossible = SynDeletePossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynDeleteArray[countD] = RN1;
						SynDeleteArrayCopy[countD] = RN1;
						SynDeleteArrayIndex[countD] = countD;
					}
					else
					{
						SynDeleteArray[countD] = 0;
						SynDeleteArrayCopy[countD] = 0;
						SynDeleteArrayIndex[countD] = countD;
					}					
					countD = countD + 1;
				}
			}
		}
	}
	else
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (conMatEI[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]]==1)
					{
						SynDeletePossible = SynDeletePossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynDeleteArray[countD] = RN1;
						SynDeleteArrayCopy[countD] = RN1;
						SynDeleteArrayIndex[countD] = countD;
					}
					else
					{
						SynDeleteArray[countD] = 0;
						SynDeleteArrayCopy[countD] = 0;
						SynDeleteArrayIndex[countD] = countD;
					}					
					countD = countD + 1;
				}
			}
		}
		else
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (conMatEE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]]==1)
					{
						SynDeletePossible = SynDeletePossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynDeleteArray[countD] = RN1;
						SynDeleteArrayCopy[countD] = RN1;
						SynDeleteArrayIndex[countD] = countD;
					}
					else
					{
						SynDeleteArray[countD] = 0;
						SynDeleteArrayCopy[countD] = 0;
						SynDeleteArrayIndex[countD] = countD;
					}					
					countD = countD + 1;
				}
			}
		}
	}
	
	cout << "process_id: " << process_id << "\t SynDeletePossible = " << SynDeletePossible << endl;	

	
	int SynAdditionPossible = 0;
	double* SynAddArray;
	SynAddArray = new double[N[2]*N[3]];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[N[2]*N[3]];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[N[2]*N[3]];
	int countA = 0;
	if ( neuron[2]>=NoNe )
	{
		if ( neuron[3]>=NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (conMatII[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]==0)
					{
						SynAdditionPossible = SynAdditionPossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynAddArray[countA] = RN1;
						SynAddArrayCopy[countA] = RN1;
						SynAddArrayIndex[countA] = countA;
					}
					else
					{
						SynAddArray[countA] = 0;
						SynAddArrayCopy[countA] = 0;
						SynAddArrayIndex[countA] = countA;
					}
					countA = countA + 1;
				}
			}		
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (conMatIE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]==0)
					{
						SynAdditionPossible = SynAdditionPossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynAddArray[countA] = RN1;
						SynAddArrayCopy[countA] = RN1;
						SynAddArrayIndex[countA] = countA;
					}
					else
					{
						SynAddArray[countA] = 0;
						SynAddArrayCopy[countA] = 0;
						SynAddArrayIndex[countA] = countA;
					}
					countA = countA + 1;
				}
			}
		}
	}
	else	
	{
		if ( neuron[3]>=NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (conMatEI[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]==0)
					{
						SynAdditionPossible = SynAdditionPossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynAddArray[countA] = RN1;
						SynAddArrayCopy[countA] = RN1;
						SynAddArrayIndex[countA] = countA;
					}
					else
					{
						SynAddArray[countA] = 0;
						SynAddArrayCopy[countA] = 0;
						SynAddArrayIndex[countA] = countA;
					}
					countA = countA + 1;
				}
			}
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (conMatEE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]==0)
					{
						SynAdditionPossible = SynAdditionPossible + 1;
						do
						{
							RN1 = double(rand()) / double(RAND_MAX);
							RN2_copy = int(RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;

						SynAddArray[countA] = RN1;
						SynAddArrayCopy[countA] = RN1;
						SynAddArrayIndex[countA] = countA;
					}
					else
					{
						SynAddArray[countA] = 0;
						SynAddArrayCopy[countA] = 0;
						SynAddArrayIndex[countA] = countA;
					}
					countA = countA + 1;
				}
			}
		}
	}
	
	cout << "process_id: " << process_id << "\t SynAdditionPossible = " << SynAdditionPossible << endl;

	int SynMutatePossible = SynDeletePossible;
	if (SynAdditionPossible<SynDeletePossible) SynMutatePossible = SynAdditionPossible;
	cout << "process_id: " << process_id << "\t SynMutatePossible = " << SynMutatePossible << endl;

	int MinSynMutate = 5;
	int SynMutate = 0;
	if (SynMutatePossible<=MinSynMutate) SynMutate = SynMutatePossible;
	else // (SynMutatePossible>MinSynMutate)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		SynMutate = round(RN1*double(SynMutatePossible));
		if (SynMutate<=MinSynMutate) SynMutate = MinSynMutate;
	}
	if ( SynMutate==0 ) SynMutate = SynMutatePossible;

	cout << "process_id: " << process_id << "\t SynMutate = " << SynMutate << endl;

	DescendSort(SynDeleteArray, SynDeleteArrayIndex, countD);
	DescendSort(SynAddArray, SynAddArrayIndex, countA);
	
	
	double SynDeleteThreshold = SynDeleteArray[SynMutate-1];
	double SynAddThreshold = SynAddArray[SynMutate-1];

if ( SynMutate > 0 )
{
	int countDD = 0;
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
					{						
						conMatII[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = 0;
					}
					countDD = countDD + 1;
				}
			}
		}
		else
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
					{
						conMatIE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = 0;
					}
					countDD = countDD + 1;
				}
			}
		}
	}
	else
	{
		if ( neuron[1]>NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
					{
						conMatEI[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = 0;
					}
					countDD = countDD + 1;
				}
			}
		}
		else
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
					{
						conMatEE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = 0;
					}
					countDD = countDD + 1;
				}
			}
		}
	}

	
	int countAA = 0;
	if ( neuron[2]>=NoNe )
	{
		if ( neuron[3]>=NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (SynAddArrayCopy[countAA]>=SynAddThreshold)
					{
						conMatII[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = 1;
					}
					countAA = countAA + 1;
				}
			}
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (SynAddArrayCopy[countAA]>=SynAddThreshold)
					{
						conMatIE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = 1;
					}
					countAA = countAA + 1;
				}
			}
		}
	}
	else
	{
		if ( neuron[3]>=NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (SynAddArrayCopy[countAA]>=SynAddThreshold)
					{
						conMatEI[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = 1;
					}
					countAA = countAA + 1;
				}
			}
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					if (SynAddArrayCopy[countAA]>=SynAddThreshold)
					{
						conMatEE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = 1;
					}
					countAA = countAA + 1;
				}
			}
		}
	}
}
	
	
	delete [] SynDeleteArray;
	delete [] SynDeleteArrayCopy;
	delete [] SynDeleteArrayIndex;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

	delete [] dist0index;
	delete [] dist0;
	delete [] dist1index;
	delete [] dist1;
	delete [] dist2index;
	delete [] dist2;
	delete [] dist3index;
	delete [] dist3;

	delete [] neuron;
	delete [] N;
	delete [] Ng;	

}


void Network::MutationEERowSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationEERowSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNe];	
	dist0index = new int[NoNe];
	
	for ( int i=0 ; i<NoNe ; i++ )
	{
		dist0[i] = pow( (gLeake[i]-gLeake[neuron]),2 ) + pow( (gNaPe[i]-gNaPe[neuron]),2 );
		dist0index[i] = i;		
	}	
	
	DescendSort(dist0, dist0index, NoNe);	

	for ( int i=0 ; i<NoNe ; i++ )
	{
		if ( conMatEE[neuron][i] == 0 )
		{
			int j = 0;
			while ( j < NoNe-1 && conMatEE[neuron][i] == 0 )
			{
				if ( conMatEE[dist0index[NoNe-2-j]][i] == 1 )
				{
					conMatEE[neuron][i] = 1;
					conMatEE[dist0index[NoNe-2-j]][i] = 0;
				}
				j = j + 1;
			}		
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}


void Network::MutationEIRowSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationEIRowSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNe];
	dist0index = new int[NoNe];	

	for ( int i=0 ; i<NoNe ; i++ )
	{
		dist0[i] = pow( (gLeake[i]-gLeake[neuron]),2 ) + pow( (gNaPe[i]-gNaPe[neuron]),2 );
		dist0index[i] = i;		
	}
	
	DescendSort(dist0, dist0index, NoNe);

	for ( int i=0 ; i<NoNi ; i++ )
	{
		if ( conMatEI[neuron][i] == 0 )
		{
			int j = 0;
			while ( j < NoNe-1 && conMatEI[neuron][i] == 0 )
			{
				if ( conMatEI[dist0index[NoNe-2-j]][i] == 1 )
				{
					conMatEI[neuron][i] = 1;
					conMatEI[dist0index[NoNe-2-j]][i] = 0;
				}
				j = j + 1;
			}		
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}


void Network::MutationIERowSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationIERowSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNi];	
	dist0index = new int[NoNi];

	for ( int i=0 ; i<NoNi ; i++ )
	{
		dist0[i] = pow( (gLeaki[i]-gLeaki[neuron]),2 ) + pow( (gNaPi[i]-gNaPi[neuron]),2 );
		dist0index[i] = i;		
	}
	
	
	DescendSort(dist0, dist0index, NoNi);

	for ( int i=0 ; i<NoNe ; i++ )
	{
		if ( conMatIE[neuron][i] == 0 )
		{
			int j = 0;
			while ( j < NoNi-1 && conMatIE[neuron][i] == 0 )
			{
				if ( conMatIE[dist0index[NoNi-2-j]][i] == 1 )
				{
					conMatIE[neuron][i] = 1;
					conMatIE[dist0index[NoNi-2-j]][i] = 0;
				}
				j = j + 1;
			}	
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}


void Network::MutationIIRowSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationIIRowSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNi];	
	dist0index = new int[NoNi];
	
	for ( int i=0 ; i<NoNi ; i++ )
	{
		dist0[i] = pow( (gLeaki[i]-gLeaki[neuron]),2 ) + pow( (gNaPi[i]-gNaPi[neuron]),2 );
		dist0index[i] = i;		
	}	
	
	DescendSort(dist0, dist0index, NoNi);	

	for ( int i=0 ; i<NoNi ; i++ )
	{
		if ( conMatII[neuron][i] == 0 )
		{
			int j = 0;
			while ( j < NoNi-1 && conMatII[neuron][i] == 0 )
			{
				if ( conMatII[dist0index[NoNi-2-j]][i] == 1 )
				{
					conMatII[neuron][i] = 1;
					conMatII[dist0index[NoNi-2-j]][i] = 0;
				}				
				j = j + 1;
			}		
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}


void Network::MutationRowSeg()
{
	/* pick one neuron at random and increase its outgoing synapses
	 * by swaping the outgoing synapse of it's neighbors; the number
	 * of outgoing synapse of its neighbor decreases as result */	

	int neuron1, neuron2;
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if ( k==0 ) neuron1 = floor(double(NoP)*RN1);
		else neuron2 = floor(double(NoP)*RN1);
	}
	cout << "precess_id: " << process_id << "\t neuron1 = " << neuron1 << "\t neuron2 = " << neuron2 << endl;

	if ( neuron1<NoNe )
	{
		if ( neuron2<NoNe ) MutationEERowSeg(neuron1);
		else MutationEIRowSeg(neuron1);
	}
	else
	{
		if ( neuron2<NoNe ) MutationIERowSeg(neuron1-NoNe);
		else MutationIIRowSeg(neuron1-NoNe);
	}	

}



void Network::MutationEEColSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationEEColSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNe];	
	dist0index = new int[NoNe];
	
	for ( int i=0 ; i<NoNe ; i++ )
	{
		dist0[i] = pow( (gLeake[i]-gLeake[neuron]),2 ) + pow( (gNaPe[i]-gNaPe[neuron]),2 );
		dist0index[i] = i;		
	}	
	
	DescendSort(dist0, dist0index, NoNe);	

	for ( int i=0 ; i<NoNe ; i++ )
	{
		if ( conMatEE[i][neuron] == 0 )
		{
			int j = 0;
			while ( j < NoNe-1 && conMatEE[i][neuron] == 0 )
			{
				if ( conMatEE[i][dist0index[NoNe-2-j]] == 1 )
				{
					conMatEE[i][neuron] = 1;
					conMatEE[i][dist0index[NoNe-2-j]] = 0;
				}
				j = j + 1;
			}	
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}


void Network::MutationEIColSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationEIColSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNi];
	dist0index = new int[NoNi];
	
	for ( int i=0 ; i<NoNi ; i++ )
	{
		dist0[i] = pow( (gLeaki[i]-gLeaki[neuron]),2 ) + pow( (gNaPi[i]-gNaPi[neuron]),2 );
		dist0index[i] = i;		
	}	
	
	DescendSort(dist0, dist0index, NoNi);	

	for ( int i=0 ; i<NoNe ; i++ )
	{
		if ( conMatEI[i][neuron] == 0 )
		{
			int j = 0;
			while ( j < NoNi-1 && conMatEI[i][neuron] == 0 )
			{
				if ( conMatEI[i][dist0index[NoNi-2-j]] == 1 )
				{
					conMatEI[i][neuron] = 1;
					conMatEI[i][dist0index[NoNi-2-j]] = 0;
				}
				j = j + 1;
			}		
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}


void Network::MutationIEColSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationIEColSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNe];	
	dist0index = new int[NoNe];
	
	for ( int i=0 ; i<NoNe ; i++ )
	{
		dist0[i] = pow( (gLeake[i]-gLeake[neuron]),2 ) + pow( (gNaPe[i]-gNaPe[neuron]),2 );
		dist0index[i] = i;
	}	
	
	DescendSort(dist0, dist0index, NoNe);
	
	for ( int i=0 ; i<NoNi ; i++ )
	{
		if ( conMatIE[i][neuron] == 0 )
		{
			int j = 0;
			while ( j < NoNe-1 && conMatIE[i][neuron] == 0 )
			{
				if ( conMatIE[i][dist0index[NoNe-2-j]] == 1 )
				{
					conMatIE[i][neuron] = 1;
					conMatIE[i][dist0index[NoNe-2-j]] = 0;
				}
				j = j + 1;
			}		
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}


void Network::MutationIIColSeg(int neuron)
{	
	cout << "precess_id: " << process_id << "\t MutationIIColSeg: neuron = " << neuron << endl;	

	double* dist0;
	int* dist0index;
	dist0 = new double[NoNi];	
	dist0index = new int[NoNi];	

	for ( int i=0 ; i<NoNi ; i++ )
	{
		dist0[i] = pow( (gLeaki[i]-gLeaki[neuron]),2 ) + pow( (gNaPi[i]-gNaPi[neuron]),2 );
		dist0index[i] = i;		
	}
	
	
	DescendSort(dist0, dist0index, NoNi);	

	for ( int i=0 ; i<NoNi ; i++ )
	{
		if ( conMatII[i][neuron] == 0 )
		{
			int j = 0;
			while ( j < NoNi-1 && conMatII[i][neuron] == 0 )
			{
				if ( conMatII[i][dist0index[NoNi-2-j]] == 1 )
				{
					conMatII[i][neuron] = 1;
					conMatII[i][dist0index[NoNi-2-j]] = 0;
				}				
				j = j + 1;
			}		
		}
	}
	
	delete [] dist0index;
	delete [] dist0;
}



void Network::MutationColSeg()
{
	/* pick one neuron at random and increase its incoming synapses
	 * by swaping the incoming synapse of it's neighbors; the number
	 * of incoming synapse of its neighbor decreases as result */	

	int neuron1, neuron2;
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if ( k==0 ) neuron1 = floor(double(NoP)*RN1);
		else neuron2 = floor(double(NoP)*RN1);
	}
	cout << "precess_id: " << process_id << "\t neuron1 = " << neuron1 << "\t neuron2 = " << neuron2 << endl;

	if ( neuron1<NoNe )
	{
		if ( neuron2<NoNe ) MutationEEColSeg(neuron1);
		else MutationIEColSeg(neuron1);
	}
	else
	{
		if ( neuron2<NoNe ) MutationEIColSeg(neuron1-NoNe);
		else MutationIIColSeg(neuron1-NoNe);
	}	

}


void Network::MutationRowColSeg(int iter)
{
	srand ( time(NULL)+process_id+iter*1000);
	double RN = double(rand()) / double(RAND_MAX);	
	if ( RN > 0.5 ) MutationColSeg();
	else MutationRowSeg();
}

netconfig_cond Network::MutatateNetwork(int iter, int PsuedoMutate)
{
	srand ( time(NULL)+process_id+iter*1010);
	double RN = double(rand()) / double(RAND_MAX);
	RN = 0.3;
	if (PsuedoMutate==0)
	{
		if ( RN > 0.5 ) MutationRowColSeg(iter);
		else ConnectMutationSeg(iter);
	}

	netconfig_cond nc;
	nc = Initialize(11); // any integer > 2 would suffice here	

	return nc;
}


void Network::writeNetworkCondition(netconfig_cond nc, int a) 
{
	ofstream file1b;
	if ( a==0 ) file1b.open ("NeuronPropInitial.txt");
	else file1b.open ("NeuronPropFinal.txt");
	for (int i=0; i<NoP; i++)
	{	
		file1b << nc.GLeak[i] << "\t" << nc.GNaP[i] << endl;
	}
	file1b.close();


	ostringstream stream;
	if ( a==0 )	stream << "InitialConMat.txt";
	else stream << "FinalConMat.txt";
	ofstream file0b;
	file0b.open (stream.str().c_str());
	for (int i=0; i<NoP; i++)
	{
		for (int j=0; j<NoP; j++)
		{			
			if (j==NoN-1) file0b << nc.conMat[i][j] << endl;
			else file0b << nc.conMat[i][j] << "\t" ;
		}
	}
	file0b.close();	
	
}

netconfig_cond Network::AddInhibitoryConnection()
{
	
	
	ComputePopBurstDataA();

	
	int* Nepop;
	
	Nepop = new int[NoNe];
	
	int* NepopIndex;
	
	NepopIndex = new int[NoNe];
	
	for (int k=0; k<NoNe; k++)
	
	{
		
		Nepop[k] = EctopicBurstData[k];
		
		NepopIndex[k] = k;
	
	}
	
	
	DescendSortInteger(Nepop, NepopIndex, NoNe);
	

	
	int* Nipop;
	
	Nipop = new int[NoNi];
	
	int* NipopIndex;
	
	NipopIndex = new int[NoNi];
	
		
	for (int k=NoNe; k<NoN; k++)
	
	{
		
		Nipop[k-NoNe] = EctopicBurstData[k];
		
		NipopIndex[k-NoNe] = k-NoNe;
	
	}
	
		
	DescendSortInteger(Nipop, NipopIndex, NoNi);
	

	
	int NoIC = SynFrac * NoNi;
	
	for (int i=0; i<NoNi; i++)		
	{
		
		for (int j=0; j<NoNe; j++)
		
		{
			
			conMatIE[i][j] = 0;
		
		}
		
	}

		
	for (int i=0; i<NoIC; i++)
	
	{
		
		int ii = NipopIndex[i];
		
		for (int j=0; j<NoIC; j++)
		
		{
			
			int jj = NepopIndex[j];
			
			conMatIE[ii][jj] = 1;
		
		}
	
	}
	
		
	cout << "process_id : " << process_id << " AddInhibitoryConnection ready " << endl;

	
		
	netconfig_cond nc;
	
	nc = Initialize(3); // any integer > 2 would suffice here

	
	return nc;	
	
}


void Network::ComputePopBurstDataA()
{
	
	
	cout << "process_id: " << process_id
			
		<< "  BLcount = " << BLcount << endl;
	
		
	NoPoPburst = 0;
		
	int actBegin = 0;		
	int actEnd = 0;
	
	int j;
	
	int popThreshold = int(round(popThresholdFrac*NoN));
	
	int Bstart[NoEB];
	
	int Bend[NoEB];

	
		
	//ofstream file1a;
	
	//file1a.open ("actTiming.txt");
	
	for (int i=0; i<BLcount-1; i++)
	
	{
		
		if (BurstData[i][NoN+1]>popThreshold)
		
		{
			
			actBegin = i;
		
		}
		
		if (BurstData[i][NoN+1]>=popThreshold)
		
		{
			
			NoPoPburst = NoPoPburst + 1;
			
			//cout << "NoPoPburst = " << NoPoPburst << endl;
			
			j=0;
			
			while ( BurstData[i+j][NoN+1]>=popThreshold && (i+j)<BLcount-1 )
			
			{
				
				j = j + 1;
				
				actEnd = i+j;
			
			}

			
			int NNpopLength = actEnd-actBegin+1;
			
			int* NNpop;
		 
			NNpop = new int[NNpopLength];
			
			int* NNpopIndex;
			
			NNpopIndex = new int[NNpopLength];
			
			for (int k=actBegin; k<=actEnd; k++)
			
			{
				
				NNpop[k-actBegin] = BurstData[k][NoN+1];
				
				NNpopIndex[k-actBegin] = k;
			
			}
			
				
			DescendSortInteger(NNpop, NNpopIndex, NNpopLength); // after sorting t[NNpopIndex[0]] is the timing when population burst is maximum
			
			//file1a << NoPoPburst << "\t" << t[NNpopIndex[0]] << "\t" << actBegin+1 << "\t" << actEnd-1 << "\t" << t[actBegin+1] << "\t" << t[actEnd-1] << endl;
			
			Bstart[NoPoPburst-1] = actBegin+1;
			
			Bend[NoPoPburst-1] = actEnd-1;
			
			delete NNpop;
			
			delete NNpopIndex;

			
			i = actEnd;
		
				
		}
		
			
	}
	
	
	//cout << "process_id: " << process_id << " NoPoPburst = " << NoPoPburst << endl;
	
	//file1a.close();

	//int PopBurstData[NoEB][NoP];
	
	
	for (int i=0; i<NoEB; i++)
	
	{
		
		for (int j=0; j<NoN; j++)
		
		{
			
			PopBurstDataA[i][j] = 0;
		
		}
	
	}
	
	for (int j=0; j<NoN; j++)
	
	{
		
		EctopicBurstData[j] = 0;
	
	}

	
	int FirstBurstStart = 0;
	
	int LastBurstEnd = NoPoPburst;
	
	if (Bstart[0]==0) FirstBurstStart = 1;
	
	if (Bstart[NoPoPburst-1]==BLcount-1) 
	LastBurstEnd = NoPoPburst-1;
	
	NoPoPburst = LastBurstEnd - FirstBurstStart;
	
	cout << "process_id: " << process_id << " NoPoPburst = " << NoPoPburst << endl;
	

	
	int Tcount;
	
	for (int j=0; j<NoN; j++)
	
	{
		
		for (int i=0; i<NoPoPburst; i++)
		
		{
			
			Tcount = Bstart[i+FirstBurstStart];
			
			while ( Tcount<=Bend[i+FirstBurstStart] && PopBurstDataA[i][j]==0)
			
			{
				
				if (BurstData[Tcount][j] == 1) PopBurstDataA[i][j] = 1;
				
				Tcount = Tcount + 1;
			
			}
		
		}
	
	}
	
		
	cout << "process_id: " << process_id << " PoPBurstData ready " << endl;
	
		
	for (int j=0; j<NoN; j++)
	
	{
		
		for (int i=0; i<NoPoPburst-1; i++)
		
		{
			
			Tcount = Bend[i+FirstBurstStart] + 1;
			
			while (Tcount<Bstart[i+FirstBurstStart+1])
			
			{
				
				if (BurstData[Tcount][j] == 1) EctopicBurstData[j] += 1;
				
				//cout << "EctopicBurstData[" << j << "]= " << EctopicBurstData[j] << endl;
				
				Tcount = Tcount + 1;
			
			}
		
		}
	
	}

	
					
	ostringstream stream2;
	
	ofstream myfile0q;
	
	stream2 << "EctopicBurstData" << process_id << ".txt";
	
	myfile0q.open (stream2.str().c_str());
	
	for (int j=0; j<NoN; j++)
	
	{
		
		myfile0q << EctopicBurstData[j] << endl;
	
	}

	myfile0q.close();
	
		
	cout << "process_id: " << process_id << " EctopicBurstData ready " << endl;


		
}






netconfig_cond Network::InitializeG(int z) {
	
	srand ( time(NULL)+process_id );
	netconfig_cond nc;
	
	if ( z==0 )
	{
		NetTag = process_id;
		InitializeNeuronProp();
		InitializeNetworkTopologyG();
	}
	
	if ( z==1 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		InitializeNetworkTopologyG();
	}

	if ( z==2 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		ReadConMat();
	}

	if ( z==3 )
	{
		NetTag = process_id;
		ReadNeuronProp();
		InitializeNetworkTopologyG();
		CopyEE();
	}

	
	for (int i=0; i<NoN; i++)
	{
		if ( i<NoNe )
		{
			nc.GNaP[i] = gNaPe[i];
			nc.GLeak[i] = gLeake[i];
		}
		else
		{
			nc.GNaP[i] = gNaPi[i-NoNe];
			nc.GLeak[i] = gLeaki[i-NoNe];
		}		
	}

	
	for (int m = 0; m < NoN; m++)
	{
		//gtonic_e[m] = tonicConductance;
		for(int mm = 0; mm < NoN; mm++)
		{
			gsyn_e[m][mm] = synS; // index m is for pre-synaptic neuron
		}
		cp[m].Gtonic_e(gtonic_e[m]);
		if ( m<NoNe )
		{
			cp[m].GNaPo_set(gNaPe[m]);
			cp[m].GLeako_set(gLeake[m]);
			gtonic_e[m] = tonicConductance;
		}
		else
		{
			cp[m].GNaPo_set(gNaPi[m-NoNe]);
			cp[m].GLeako_set(gLeaki[m-NoNe]);
			gtonic_e[m] = tonicConductanceI;
		}	
		
		
	}

	
	for (int i=0; i<NoN; i++)
	{
		for (int j=0; j<NoN; j++)
		{
			if ( i<NoNe )
			{
				if ( j<NoNe )	nc.conMat[i][j] = conMatEE[i][j];
				else nc.conMat[i][j] = conMatEI[i][j-NoNe];
			}
			else
			{
				if ( j<NoNe )	nc.conMat[i][j] = conMatIE[i-NoNe][j];
				else nc.conMat[i][j] = conMatII[i-NoNe][j-NoNe];
			}
		}
	}
	
	nc.NetTag = NetTag;

	return nc;
}

void Network::InitializeNetworkTopologyG()
{
	for ( int i=1 ; i<=4; i++ )
	{		
		// NetworkTopologyG(i);
		NetworkTopologyHybrid(i);
	}
}

void Network::NetworkTopologyHybrid(int MatTag)
{
	srand ( time(NULL)+process_id+101*MatTag );

	int M, N;
	if ( MatTag==1 || MatTag==2 ) M = NoNe; 
	else M = NoNi;

	if ( MatTag==1 || MatTag==3 ) N = NoNe; 
	else N = NoNi;	

	double SynFrac0;
	if (MatTag==1) SynFrac0 = SynFrac;
	if (MatTag==2) SynFrac0 = SynFracEI;
	if (MatTag==3) SynFrac0 = SynFracIE;
	if (MatTag==4) SynFrac0 = SynFracII;

	int NoSyn = int( round(SynFrac0*double(M*N)) );	
	int RN1_copy = 0;
	int RN2_copy = 0;

	double* RanMat;
	RanMat = new double[M*N];
	double* RanMatCopy;
	RanMatCopy = new double[M*N];
	int* RanIndex;
	RanIndex = new int[M*N];
	double RN1;
	srand(time(NULL));
	for (int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			do
			{
				RN1 = (double)rand( ) / RAND_MAX;
				RN2_copy = (int) (RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;
			RanMat[i*N+j] = RN1;
			RanMatCopy[i*N+j] = RN1;
			RanIndex[i*N+j] = i*N+j;			
		}
	}

	DescendSort(RanMat,RanIndex,M*N);
	

	double threshold;
	int thresholdIndex = NoSyn-1;
	if ( thresholdIndex < 0 ) threshold = 1.0;
	else threshold = RanMat[thresholdIndex];
	delete [] RanMat;
	delete [] RanIndex;
	
	if ( MatTag == 1 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) 
				{
					do
					{
						RN1 = (double)rand( ) / RAND_MAX;
						RN2_copy = (int) (RN1*100000);
					} while (RN1_copy == RN2_copy);
					RN1_copy = RN2_copy;
					conMatEE[i][j] = ceil(RN1*synGradeEE);
				}
				else conMatEE[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 2 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) 
				{
					do
					{
						RN1 = (double)rand( ) / RAND_MAX;
						RN2_copy = (int) (RN1*100000);
					} while (RN1_copy == RN2_copy);
					RN1_copy = RN2_copy;
					conMatEI[i][j] = ceil(RN1*synGradeEI);
				}
				else conMatEI[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 3 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) 
				{
					do
					{
						RN1 = (double)rand( ) / RAND_MAX;
						RN2_copy = (int) (RN1*100000);
					} while (RN1_copy == RN2_copy);
					RN1_copy = RN2_copy;
					conMatIE[i][j] = ceil(RN1*synGradeIE);
				}
				else conMatIE[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 4 )
	{			
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				if (RanMatCopy[i*N+j]>=threshold) 
				{
					do
					{
						RN1 = (double)rand( ) / RAND_MAX;
						RN2_copy = (int) (RN1*100000);
					} while (RN1_copy == RN2_copy);
					RN1_copy = RN2_copy;
					conMatII[i][j] = ceil(RN1*synGradeII);
				}
				else conMatII[i][j] = 0;
			}
		}
		delete [] RanMatCopy;
	}




	//ConnectIsolatedNeurons(conMat, M);
	//StoreInitialConMat();
}


void Network::NetworkTopologyG(int MatTag)
{
	int M, N;
	if ( MatTag==1 || MatTag==2 ) M = NoNe; 
	else M = NoNi;

	if ( MatTag==1 || MatTag==3 ) N = NoNe; 
	else N = NoNi;

	double SynFrac0;
	if (MatTag==1) SynFrac0 = SynFrac;
	if (MatTag==2) SynFrac0 = SynFracEI;
	if (MatTag==3) SynFrac0 = SynFracIE;
	if (MatTag==4) SynFrac0 = SynFracII;

	int NoSyn = int( round(SynFrac0*double(M*N*synGrade)) );	
	int RN1_copy = 0;
	int RN2_copy = 0;

	double* RanMat;
	RanMat = new double[M*N*synGrade];
	double* RanMatCopy;
	RanMatCopy = new double[M*N*synGrade];
	int* RanIndex;
	RanIndex = new int[M*N*synGrade];
	double RN1;
	srand(time(NULL));
	for (int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			for (int g=0; g<synGrade; g++)
			{
				do
				{
					RN1 = (double)rand( ) / RAND_MAX;
					RN2_copy = (int) (RN1*100000);
				} while (RN1_copy == RN2_copy);
				RN1_copy = RN2_copy;
				RanMat[(i*N+j)*synGrade+g] = RN1;
				RanMatCopy[(i*N+j)*synGrade+g] = RN1;
				RanIndex[(i*N+j)*synGrade+g] = (i*N+j)*synGrade+g;
			}			
		}
	}

	DescendSort(RanMat,RanIndex,M*N*synGrade);
	

	double threshold;
	int thresholdIndex = NoSyn-1;
	if ( thresholdIndex < 0 ) threshold = 1.0;
	else threshold = RanMat[thresholdIndex];
	delete [] RanMat;
	delete [] RanIndex;
	
	if ( MatTag == 1 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				conMatEE[i][j] = 0;
				for (int g=0; g<synGrade; g++)
				{
					if (RanMatCopy[(i*N+j)*synGrade+g]>=threshold) conMatEE[i][j] = conMatEE[i][j] + 1;
				}
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 2 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				conMatEI[i][j] = 0;
				for (int g=0; g<synGrade; g++)
				{
					if (RanMatCopy[(i*N+j)*synGrade+g]>=threshold) conMatEI[i][j] = conMatEI[i][j] +1;
				}
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 3 )
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				conMatIE[i][j] = 0;
				for (int g=0; g<synGrade; g++)
				{
					if (RanMatCopy[(i*N+j)*synGrade+g]>=threshold) conMatIE[i][j] = conMatIE[i][j] + 1;
				}
			}
		}
		delete [] RanMatCopy;
	}

	if ( MatTag == 4 )
	{			
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				conMatII[i][j] = 0;
				for (int g=0; g<synGrade; g++)
				{
					if (RanMatCopy[(i*N+j)*synGrade+g]>=threshold) conMatII[i][j] = conMatII[i][j] + 1;
				}
			}
		}
		delete [] RanMatCopy;
	}

	//ConnectIsolatedNeurons(conMat, M);
	//StoreInitialConMat();
}

netconfig_cond Network::MutatateNetworkG(int iter, int PsuedoMutate)
{
	srand ( time(NULL)+process_id+iter*1010);
	double RN = double(rand()) / double(RAND_MAX);
	if (PsuedoMutate==0)
	{
		if ( RN < 0.2 ) MutationPartialRowTranslocate(iter);
		if ( RN >= 0.2 && RN < 0.4 ) MutationPartialColumnTranslocate(iter);
		if ( RN >= 0.4 && RN < 0.6 ) MutationRowTranslocate(iter);
		if ( RN >= 0.6 && RN < 0.7 ) MutationColumnTranslocate(iter);
		if ( RN >= 0.7 && RN < 0.8 ) MutationRowSwitch(iter);
		if ( RN >= 0.8 && RN < 0.9 ) MutationColumnSwitch(iter);
		if ( RN > 0.9 ) MutationNeuronSwitch(iter);
	}

	netconfig_cond nc;
	nc = Initialize(11); // any integer > 3 would suffice here	
	//nc = InitializeG(11); // any integer > 3 would suffice here	

	return nc;
}

netconfig_cond Network::MutatateNetworkHybrid(int iter, int PsuedoMutate)
{
	srand ( time(NULL)+process_id+iter*1010);
	double RN = double(rand()) / double(RAND_MAX);
	if (PsuedoMutate==0)
	{
		// if ( RN < 0.2 ) MutationColumnSwitchEI(iter);
		// if ( RN >= 0.2 && RN < 0.4 ) MutationRowSwitchEI(iter);
		// if ( RN >= 0.4 && RN < 0.6 ) MutationColumnSwitchIE(iter);
		// if ( RN >= 0.6 && RN < 0.7 ) MutationRowSwitchIE(iter);		
		// if ( RN >= 0.7 && RN < 0.8 ) MutationNeuronSwitchING(iter);
		// if ( RN >= 0.8 && RN < 0.9 ) MutationColumnSwitchING(iter);
		// if ( RN >= 0.8 && RN < 0.9 ) MutationRowSwitchING(iter);
		//if ( RN < 0.1 ) MutationNeuronSwitchENG(iter);
		//if ( RN >= 0.1 && RN < 0.2 ) MutationColumnSwitchENG(iter);
		//if ( RN >= 0.2 && RN < 0.3 ) MutationRowSwitchENG(iter);
		//if ( RN < 0.34 ) MutationSynAddEI(iter); 
		//if ( RN >= 0.34 && RN < 0.67 ) MutationSynAddIE(iter);
		//if ( RN >= 0.67 ) MutationSynAddII(iter); 
		//if ( RN < 0.34 ) MutationSynDeleteEI(iter); 
		//if ( RN >= 0.34 && RN < 0.67 ) MutationSynDeleteIE(iter);
		//if ( RN >= 0.67 ) MutationSynDeleteII(iter);
		
		//if ( RN >= 0.0 ) MutationSynAddEE(iter);
		//if ( RN >= 0.0 ) MutationSynDeleteEE(iter);

		//if ( RN >= 0.0 ) MutationNeuronSwitchENG(iter);
		
		if ( RN < 0.25 ) MutationPartialColTranslocateENG1(iter);
		if ( RN >= 0.25 && RN < 0.5 ) MutationPartialRowTranslocateENG1(iter);
		if ( RN >= 0.5 && RN < 0.75 ) MutationSynColAddENG1(iter);
		if ( RN >= 0.75 ) MutationSynRowAddENG1(iter);
		//if ( RN >= 0.5 && RN < 0.75 ) MutationSynColDeleteENG1(iter);
		//if ( RN >= 0.75 ) MutationSynRowDeleteENG1(iter);
		
		/*double w = 0.1666;
		if ( RN < w ) MutationPartialColTranslocateENG1(iter);
		if ( RN >= w && RN < 2*w ) MutationPartialRowTranslocateENG1(iter);
		if ( RN >= 2*w && RN < 3*w ) MutationSynColAddENG1(iter);
		if ( RN >= 3*w && RN < 4*w ) MutationSynRowAddENG1(iter);
		if ( RN >= 4*w && RN < 5*w ) MutationSynColDeleteENG1(iter);
		if ( RN >= 5*w ) MutationSynRowDeleteENG1(iter);*/

		//if ( RN < 0.5 ) ColAddIE(iter);
		//if ( RN >= 0.5 ) RowAddIE(iter);
		//if ( RN < 0.5 ) ColDeleteIE(iter);
		//if ( RN >= 0.5 ) RowDeleteIE(iter);

		/*if ( RN < 0.125 ) MutationPartialColTranslocateENG1(iter);
		if ( RN >= 0.125 && RN < 0.25 ) MutationPartialColTranslocateEI1(iter);
		if ( RN >= 0.25 && RN < 0.375 ) MutationPartialColTranslocateIE1(iter);
		if ( RN >= 0.375 && RN < 0.5 ) MutationPartialColTranslocateING1(iter);
		if ( RN >= 0.5 && RN < 0.625 )MutationPartialRowTranslocateENG1(iter);
		if ( RN >= 0.625 && RN < 0.75 )MutationPartialRowTranslocateEI1(iter);
		if ( RN >= 0.75 && RN < 0.875 )MutationPartialRowTranslocateIE1(iter);
		if ( RN >= 0.875 )MutationPartialRowTranslocateING1(iter);*/

		/*if ( RN < 0.25 ) MutationPartialColTranslocateENG1(iter);
		if ( RN >= 0.25 && RN < 0.5 ) MutationPartialColTranslocateENG(iter);
		if ( RN >= 0.5 && RN < 0.75 )MutationPartialRowTranslocateENG1(iter);
		if ( RN >= 0.75 )MutationPartialRowTranslocateENG(iter);*/
	}

	netconfig_cond nc;
	//nc = Initialize(11); // any integer > 3 would suffice here	
	nc = InitializeG(11); // any integer > 3 would suffice here	

	return nc;
}

netconfig_cond Network::MutatateNetworkHybrid0(int iter, int PsuedoMutate)
{
	srand ( time(NULL)+process_id+iter*1010);
	double RN = double(rand()) / double(RAND_MAX);

	double w = 1/6;
	if (PsuedoMutate==0)
	{
		
		if ( RN < w ) MutationPartialColTranslocateENG1(iter);
		if ( RN >= w && RN < 2*w ) MutationPartialRowTranslocateENG1(iter);
		if ( RN >= 2*w && RN < 3*w ) MutationSynColAddENG1(iter);
		if ( RN >= 3*w && RN < 4*w ) MutationSynRowAddENG1(iter);
		if ( RN >= 4*w && RN < 5*w ) MutationSynColDeleteENG1(iter);
		if ( RN >= 5*w ) MutationSynRowDeleteENG1(iter);
				
	}

	netconfig_cond nc;
	//nc = Initialize(11); // any integer > 3 would suffice here	
	nc = InitializeG(11); // any integer > 3 would suffice here	

	return nc;
}

netconfig_cond Network::MutatateNetworkHybrid1(int iter, int PsuedoMutate)
{
	srand ( time(NULL)+process_id+iter*1010);
	double RN = double(rand()) / double(RAND_MAX);
	if (PsuedoMutate==0)
	{		
		if ( RN < 0.25 ) MutationPartialColTranslocateENG1(iter);
		if ( RN >= 0.25 && RN < 0.5 ) MutationPartialRowTranslocateENG1(iter);
		if ( RN >= 0.5 && RN < 0.75 ) MutationSynColAddENG1(iter);
		if ( RN >= 0.75 ) MutationSynRowAddENG1(iter);
		//if ( RN >= 0.5 && RN < 0.75 ) MutationSynColDeleteENG1(iter);
		//if ( RN >= 0.75 ) MutationSynRowDeleteENG1(iter);
				
	}

	netconfig_cond nc;
	//nc = Initialize(11); // any integer > 3 would suffice here	
	nc = InitializeG(11); // any integer > 3 would suffice here	

	return nc;
}

void Network::MutationColumnSwitchEI(int iter)
{
	/* pick 2 neurons at random in Inhibitory NG and sawp 
	 * their incoming excitatory synapses from Excitatory NG;
	 * this is equivalent to swaping their corresponding 
	 * columns in EI connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = NoNe + floor(double(NoNi)*RN1);
	}


	int SynTemp;

	// Column swaping (begin)
	for (int i=0; i<NoNe; i++)
	{
		SynTemp = conMatEI[i][neuron[0]-NoNe];
		conMatEI[i][neuron[0]-NoNe] = conMatEI[i][neuron[1]-NoNe];
		conMatEI[i][neuron[1]-NoNe] = SynTemp;
	}
	

	cout << "process_id: " << process_id << "\t EI-column of neuron " << neuron[0] << " swaped by that neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationRowSwitchEI(int iter)
{
	/* pick 2 neurons at random in Excitatory NG and sawp 
	 * their outgoing excitatory synapses to Inhibitory NG;
	 * this is equivalent to swaping their corresponding 
	 * columns in EI connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}


	int SynTemp;

	// Row swaping (begin)
	for (int i=0; i<NoNi; i++)
	{
		SynTemp = conMatEI[neuron[0]][i];
		conMatEI[neuron[0]][i] = conMatEI[neuron[1]][i];
		conMatEI[neuron[1]][i] = SynTemp;
	}
	

	cout << "process_id: " << process_id << "\t EI-row of neuron " << neuron[0] << " swaped by that neuron " << neuron[1] << endl;	
	
	delete [] neuron;

}


void Network::MutationColumnSwitchIE(int iter)
{
	/* pick 2 neurons at random in Excitatory NG and sawp 
	 * their incoming inhibitory synapses from Inhibitory NG;
	 * this is equivalent to swaping their corresponding 
	 * columns in IE connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}


	int SynTemp;

	// Column swaping (begin)
	for (int i=0; i<NoNi; i++)
	{
		SynTemp = conMatIE[i][neuron[0]];
		conMatIE[i][neuron[0]] = conMatIE[i][neuron[1]];
		conMatIE[i][neuron[1]] = SynTemp;
	}
	

	cout << "process_id: " << process_id << "\t IE-column of neuron " << neuron[0] << " swaped by that neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationRowSwitchIE(int iter)
{
	/* pick 2 neurons at random in Inhibitory NG and sawp 
	 * their outgoing inhibitory synapses to Excitatory NG;
	 * this is equivalent to swaping their corresponding 
	 * rows in IE connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = NoNe + floor(double(NoNi)*RN1);
	}


	int SynTemp;

	// Row swaping (begin)
	for (int i=0; i<NoNe; i++)
	{
		SynTemp = conMatIE[neuron[0]-NoNe][i];
		conMatIE[neuron[0]-NoNe][i] = conMatIE[neuron[1]-NoNe][i];
		conMatIE[neuron[1]-NoNe][i] = SynTemp;
	}
	

	cout << "process_id: " << process_id << "\t IE-row of neuron " << neuron[0] << " swaped by that neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationNeuronSwitchENG(int iter)
{
	/* pick 2 neurons within excitatory NG and swap them; 
	 * equivalently, swap the incoming and outgoing synapses
	 * of first neuron with that of second, and vice versa */

	srand ( time(NULL)+process_id+iter*1000);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}


	int SynTemp;

	// Row swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		SynTemp = conMatEE[neuron[0]][i];
		conMatEE[neuron[0]][i] = conMatEE[neuron[1]][i];
		conMatEE[neuron[1]][i] = SynTemp;
	}				
	

	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		SynTemp = conMatEE[i][neuron[0]];
		conMatEE[i][neuron[0]] = conMatEE[i][neuron[1]];
		conMatEE[i][neuron[1]] = SynTemp;
	}

	cout << "process_id: " << process_id << "\t neuron " << neuron[0] << " in ENG is swaped by neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationColumnSwitchENG(int iter)
{
	/* pick 2 neurons within excitatory NG and swap their incoming synapses; 
	 * this is equivalent to swaping their corresponding columns in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}


	int SynTemp;	

	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		SynTemp = conMatEE[i][neuron[0]];
		conMatEE[i][neuron[0]] = conMatEE[i][neuron[1]];
		conMatEE[i][neuron[1]] = SynTemp;
	}

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " in ENG is swaped by that of neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationSynColDeleteENG1(int iter)
{
	/* picks 1 neuron within excitatory NG and "deletes" existing incoming synapses to it; 
	 * this is done by switching column elements from 1 to 0 in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;
	int synCount = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		while ( synCount==0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;
		
			neuron[k] = floor(double(NoNe)*RN1);

			synCount = 0;
			for (int i=0; i<NoNe; i++) 
			{
				synCount = synCount + conMatEE[i][neuron[0]];
			}
		}
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNe];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNe];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEE[i][neuron[0]]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNe);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) <= countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t Existing syn deleted from Column of neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatEE[i][neuron[0]] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}


void Network::MutationSynRowDeleteENG1(int iter)
{
	/* picks 1 neuron within excitatory NG and "deletes" existing outgoing synapses to it; 
	 * this is done by switching row elements from 1 to 0 in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;
	int synCount = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		while ( synCount==0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;
		
			neuron[k] = floor(double(NoNe)*RN1);

			synCount = 0;
			for (int i=0; i<NoNe; i++) 
			{
				synCount = synCount + conMatEE[neuron[0]][i];
			}
		}
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNe];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNe];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEE[neuron[0]][i]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNe);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) <= countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t Existing syn deleted from Row of neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatEE[neuron[0]][i] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}


void Network::MutationSynColAddENG1(int iter)
{
	/* picks 1 neuron within excitatory NG and "adds" new incoming synapses to it; 
	 * this is done by switching column elements from 0 to 1 in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNe];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNe];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEE[i][neuron[0]]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNe);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t New syn added to Column of neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatEE[i][neuron[0]] = 1;
		}
	}

	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}


void Network::MutationSynRowAddENG1(int iter)
{
	/* picks 1 neuron within excitatory NG and "adds" new outgoing synapses to it; 
	 * this is done by switching Row elements from 0 to 1 in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNe];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNe];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEE[neuron[0]][i]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNe);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t New syn added to Row of neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatEE[neuron[0]][i] = 1;
		}
	}

	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}


void Network::ColDeleteII(int iter)
{
	/* picks 1 neuron within ihibitory NG and "deletes" existing incoming synapses to it; 
	 * this is done by switching column elements from 1 to 0 in ConMatII */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNi];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNi];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatII[i][neuron[0]]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNi);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) <= countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t ColDeleteII, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within conMatII
	for (int i=0; i<NoNi; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatII[i][neuron[0]] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}

void Network::RowDeleteII(int iter)
{
	/* picks 1 neuron within ihibitory NG and "deletes" existing outgoing synapses to it; 
	 * this is done by switching row elements from 1 to 0 in ConMatII */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNi];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNi];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNi];
	int countD = 0;
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatII[neuron[0]][i]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNi);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) <= countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t RowDeleteII, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within conMatII
	for (int i=0; i<NoNi; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatII[neuron[0]][i] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}


void Network::ColAddII(int iter)
{
	/* picks 1 neuron within inhibitory NG and "adds" new incoming synapses to it; 
	 * this is done by switching column elements from 0 to 1 in ConMatII */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNi];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNi];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatII[i][neuron[0]]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNi);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t ColAddII, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within conMatII
	for (int i=0; i<NoNi; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatII[i][neuron[0]] = 1;
		}
	}

	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}

void Network::RowAddII(int iter)
{
	/* picks 1 neuron within inhibitory NG and "adds" new outgoing synapses to it; 
	 * this is done by switching row elements from 0 to 1 in ConMatII */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNi];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNi];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatII[neuron[0]][i]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNi);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t RowAddII, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within conMatII
	for (int i=0; i<NoNi; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatII[neuron[0]][i] = 1;
		}
	}

	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}


void Network::ColDeleteEI(int iter)
{
	/* picks 1 neuron within inhibitory NG and "deletes" existing incoming E-synapses in it; 
	 * this is done by switching column elements from 1 to 0 in ConMatEI */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNe];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNe];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEI[i][neuron[0]]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNe);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) < countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t ColDeleteEI, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within ConMatEI
	for (int i=0; i<NoNe; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatEI[i][neuron[0]] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}


void Network::RowDeleteEI(int iter)
{
	/* picks 1 neuron within E-NG and "deletes" existing outgoing E-synapses in it; 
	 * this is done by switching row elements from 1 to 0 in ConMatEI */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNi];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNi];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatEI[neuron[0]][i]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNi);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) < countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t RowDeleteEI, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within ConMatEI
	for (int i=0; i<NoNi; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatEI[neuron[0]][i] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}

void Network::ColAddEI(int iter)
{
	/* picks 1 neuron within inhibitory NG and "adds" new incoming E-synapses to it; 
	 * this is done by switching column elements from 0 to 1 in ConMatEI */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNe];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNe];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEI[i][neuron[0]]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNe);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t ColAddEI, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within ConMatEI
	for (int i=0; i<NoNe; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatEI[i][neuron[0]] = 1;
		}
	}

	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}

void Network::RowAddEI(int iter)
{
	/* picks 1 neuron within excitatory NG and "adds" new outgoing E-synapses to it; 
	 * this is done by switching row elements from 0 to 1 in ConMatEI */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNi];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNi];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatEI[neuron[0]][i]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNi);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t RowAddEI, E-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within ConMatEI
	for (int i=0; i<NoNi; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatEI[neuron[0]][i] = 1;
		}
	}
	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}

void Network::ColDeleteIE(int iter)
{
	/* picks 1 neuron within excitatory NG and "deletes" existing incoming I-synapses to it; 
	 * this is done by switching column elements from 1 to 0 in ConMatIE */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNi];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNi];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatIE[i][neuron[0]]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNi);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) < countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t ColDeleteIE, E-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within ConMatIE
	for (int i=0; i<NoNi; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatIE[i][neuron[0]] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}

void Network::RowDeleteIE(int iter)
{
	/* picks 1 neuron within inhibitory NG and "deletes" existing outgoing I-synapses to it; 
	 * this is done by switching column elements from 1 to 0 in ConMatIE */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynDelArray;
	SynDelArray = new double[NoNe];
	double* SynDelArrayCopy;
	SynDelArrayCopy = new double[NoNe];
	int* SynDelArrayIndex;
	SynDelArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatIE[neuron[0]][i]>0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDelArray[i] = RN1;
			SynDelArrayCopy[i] = RN1;
			SynDelArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynDelArray[i] = 0;
			SynDelArrayCopy[i] = 0;
			SynDelArrayIndex[i] = i;
		}
	}

	DescendSort(SynDelArray,SynDelArrayIndex,NoNe);

	int SynDelPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynDelPossible)) < countD ) SynDelPossible = ceil(RN1*double(SynDelPossible));
	else SynDelPossible = 0;

	double threshold = SynDelArray[SynDelPossible];

	cout << "process_id: " << process_id << "\t ColDeleteIE, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynDelPossible = " << SynDelPossible << endl;	
		
	// Deleting synpases within ConMatIE
	for (int i=0; i<NoNe; i++)
	{
		if ( SynDelArrayCopy[i]>threshold )
		{
			conMatIE[neuron[0]][i] = 0;
		}
	}

	
	delete [] neuron;
	delete [] SynDelArray;
	delete [] SynDelArrayCopy;
	delete [] SynDelArrayIndex;

}

void Network::ColAddIE(int iter)
{
	/* picks 1 neuron within excitatory NG and "adds" new incoming I-synapses to it; 
	 * this is done by switching column elements from 0 to 1 in ConMatIE */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNi];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNi];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatIE[i][neuron[0]]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNi);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t ColAddIE, E-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within ConMatIE
	for (int i=0; i<NoNi; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatIE[i][neuron[0]] = 1;
		}
	}

	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}

void Network::RowAddIE(int iter)
{
	/* picks 1 neuron within inhibitory NG and "adds" new outgoing I-synapses to it; 
	 * this is done by switching column elements from 0 to 1 in ConMatIE */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[1];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<1; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		neuron[k] = floor(double(NoNi)*RN1);
	}

	
	double* SynAddArray;
	SynAddArray = new double[NoNe];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[NoNe];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatIE[neuron[0]][i]<1 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAddArray[i] = RN1;
			SynAddArrayCopy[i] = RN1;
			SynAddArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynAddArray[i] = 0;
			SynAddArrayCopy[i] = 0;
			SynAddArrayIndex[i] = i;
		}
	}

	DescendSort(SynAddArray,SynAddArrayIndex,NoNe);

	int SynAddPossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	if ( ceil(RN1*double(SynAddPossible)) < countD ) SynAddPossible = ceil(RN1*double(SynAddPossible));
	else SynAddPossible = 0;

	double threshold = SynAddArray[SynAddPossible];

	cout << "process_id: " << process_id << "\t RowAddIE, I-neuron " << neuron[0] << endl;	
	cout << "process_id: " << process_id << "\t SynAddPossible = " << SynAddPossible << endl;	
		
	// Adding synpases within ConMatIE
	for (int i=0; i<NoNe; i++)
	{
		if ( SynAddArrayCopy[i]>threshold )
		{
			conMatIE[neuron[0]][i] = 1;
		}
	}

	
	delete [] neuron;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

}

void Network::MutationPartialColTranslocateENG1(int iter)
{
	/* pick 2 neurons within excitatory NG and partially "translocates" their incoming synapses; 
	 * this is equivalent to partially translocating their corresponding columns in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;		
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNe];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNe];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEE[i][neuron[0]]>conMatEE[i][neuron[1]] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNe);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/	
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " in ENG is swaped by that of neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;	
	
	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatEE[i][neuron[0]];
			conMatEE[i][neuron[0]] = conMatEE[i][neuron[1]];
			conMatEE[i][neuron[1]] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[i][neuron[0]] = " << conMatEE[i][neuron[0]] << "\t conMatEE[i][neuron[1]] = " << conMatEE[i][neuron[1]] << endl;	
		}
	}

	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationPartialColTranslocateEI1(int iter)
{
	/* pick 2 neurons within inhibitory NG and partially "translocates" their EI incoming synapses; 
	 * this is equivalent to partially translocating their corresponding columns in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;		
		
		neuron[k] = NoNe + floor(double(NoNi)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNe];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNe];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEI[i][neuron[0]-NoNe]>conMatEI[i][neuron[1]-NoNe] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNe);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/	
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " in EI is translocated to neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;	
	
	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatEI[i][neuron[0]-NoNe];
			conMatEI[i][neuron[0]-NoNe] = conMatEI[i][neuron[1]-NoNe];
			conMatEI[i][neuron[1]-NoNe] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[i][neuron[0]] = " << conMatEE[i][neuron[0]] << "\t conMatEE[i][neuron[1]] = " << conMatEE[i][neuron[1]] << endl;	
		}
	}

	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationPartialColTranslocateIE1(int iter)
{
	/* pick 2 neurons within excitatory NG and partially "translocates" their IE incoming synapses; 
	 * this is equivalent to partially translocating their corresponding columns in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;		
		
		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNi];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNi];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatIE[i][neuron[0]]>conMatIE[i][neuron[1]] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNi);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/	
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " in IE is translocated to neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;	
	
	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNi; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatIE[i][neuron[0]];
			conMatIE[i][neuron[0]] = conMatIE[i][neuron[1]];
			conMatIE[i][neuron[1]] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[i][neuron[0]] = " << conMatEE[i][neuron[0]] << "\t conMatEE[i][neuron[1]] = " << conMatEE[i][neuron[1]] << endl;	
		}
	}

	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationPartialColTranslocateING1(int iter)
{
	/* pick 2 neurons within excitatory NG and partially "translocates" their II incoming synapses; 
	 * this is equivalent to partially translocating their corresponding columns in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;		
		
		neuron[k] = NoNe + floor(double(NoNi)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNi];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNi];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatII[i][neuron[0]]>conMatII[i][neuron[1]] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNi);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/	
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " in II is translocated to neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;	
	
	int SynTemp;
	// Column swaping (begin) within Inhibitory NG
	for (int i=0; i<NoNi; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatII[i][neuron[0]-NoNe];
			conMatII[i][neuron[0]-NoNe] = conMatII[i][neuron[1]-NoNe];
			conMatII[i][neuron[1]-NoNe] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[i][neuron[0]] = " << conMatEE[i][neuron[0]] << "\t conMatEE[i][neuron[1]] = " << conMatEE[i][neuron[1]] << endl;	
		}
	}

	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}



void Network::MutationPartialColTranslocateENG(int iter)
{
	/* pick 2 neurons within excitatory NG and partially "swaps" their incoming synapses; 
	 * this is equivalent to partially swaping their corresponding columns in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNe];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNe];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( pow(conMatEE[i][neuron[0]]-conMatEE[i][neuron[1]],2) > 0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNe);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/	
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " in ENG is swaped by that of neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;	
	
	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatEE[i][neuron[0]];
			conMatEE[i][neuron[0]] = conMatEE[i][neuron[1]];
			conMatEE[i][neuron[1]] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[i][neuron[0]] = " << conMatEE[i][neuron[0]] << "\t conMatEE[i][neuron[1]] = " << conMatEE[i][neuron[1]] << endl;	
		}
	}

	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationRowSwitchENG(int iter)
{
	/* pick 2 neurons within excitatory NG and swap their outgoing synapses; 
	 * this is equivalent to swaping their corresponding rows in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}


	int SynTemp;	

	// Row swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		SynTemp = conMatEE[neuron[0]][i];
		conMatEE[neuron[0]][i] = conMatEE[neuron[1]][i];
		conMatEE[neuron[1]][i] = SynTemp;
	}

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " in ENG is swaped by that of neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationPartialRowTranslocateENG1(int iter)
{
	/* pick 2 neurons within excitatory NG and partially "translocates" their outgoing synapses; 
	 * this is equivalent to partially translocates their corresponding rows in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNe];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNe];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatEE[neuron[0]][i]>conMatEE[neuron[1]][i] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNe);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " in ENG is swaped by that of neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;

	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatEE[neuron[0]][i];
			conMatEE[neuron[0]][i] = conMatEE[neuron[1]][i];
			conMatEE[neuron[1]][i] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[neuron[0]][i] = " << conMatEE[neuron[0]][i] << "\t conMatEE[neuron[1]][i] = " << conMatEE[neuron[1]][i] << endl;	
		}
	}
	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationPartialRowTranslocateEI1(int iter)
{
	/* pick 2 neurons within excitatory NG and partially "translocates" their EI outgoing synapses; 
	 * this is equivalent to partially translocates their corresponding rows in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNi];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNi];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatEI[neuron[0]][i]>conMatEI[neuron[1]][i] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNi);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " in EI is translocated to neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;

	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNi; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatEI[neuron[0]][i];
			conMatEI[neuron[0]][i] = conMatEI[neuron[1]][i];
			conMatEI[neuron[1]][i] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[neuron[0]][i] = " << conMatEE[neuron[0]][i] << "\t conMatEE[neuron[1]][i] = " << conMatEE[neuron[1]][i] << endl;	
		}
	}
	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationPartialRowTranslocateIE1(int iter)
{
	/* pick 2 neurons within inhibitory NG and partially "translocates" their IE outgoing synapses; 
	 * this is equivalent to partially translocates their corresponding rows in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = NoNe + floor(double(NoNi)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNe];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNe];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( conMatIE[neuron[0]-NoNe][i]>conMatIE[neuron[1]-NoNe][i] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNe);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " in IE is translocated to neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;

	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatIE[neuron[0]-NoNe][i];
			conMatIE[neuron[0]-NoNe][i] = conMatIE[neuron[1]-NoNe][i];
			conMatIE[neuron[1]-NoNe][i] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[neuron[0]][i] = " << conMatEE[neuron[0]][i] << "\t conMatEE[neuron[1]][i] = " << conMatEE[neuron[1]][i] << endl;	
		}
	}
	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationPartialRowTranslocateING1(int iter)
{
	/* pick 2 neurons within inhibitory NG and partially "translocates" their II outgoing synapses; 
	 * this is equivalent to partially translocates their corresponding rows in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = NoNe + floor(double(NoNi)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNi];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNi];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNi];
	int countD = 0;	
	
	for (int i=0; i<NoNi; i++)
	{
		if ( conMatII[neuron[0]-NoNe][i]>conMatII[neuron[1]-NoNe][i] )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNi);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " in II is translocated to neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;

	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNi; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatII[neuron[0]-NoNe][i];
			conMatII[neuron[0]-NoNe][i] = conMatII[neuron[1]-NoNe][i];
			conMatII[neuron[1]-NoNe][i] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[neuron[0]][i] = " << conMatEE[neuron[0]][i] << "\t conMatEE[neuron[1]][i] = " << conMatEE[neuron[1]][i] << endl;	
		}
	}
	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}



void Network::MutationPartialRowTranslocateENG(int iter)
{
	/* pick 2 neurons within excitatory NG and partially "swaps" their outgoing synapses; 
	 * this is equivalent to partially swaping their corresponding rows in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}

	
	double* SynTranslocateArray;
	SynTranslocateArray = new double[NoNe];
	double* SynTranslocateArrayCopy;
	SynTranslocateArrayCopy = new double[NoNe];
	int* SynTranslocateArrayIndex;
	SynTranslocateArrayIndex = new int[NoNe];
	int countD = 0;	
	
	for (int i=0; i<NoNe; i++)
	{
		if ( pow(conMatEE[neuron[0]][i]-conMatEE[neuron[1]][i],2) > 0 )
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynTranslocateArray[i] = RN1;
			SynTranslocateArrayCopy[i] = RN1;
			SynTranslocateArrayIndex[i] = i;
			countD = countD + 1;
		}
		else
		{
			SynTranslocateArray[i] = 0;
			SynTranslocateArrayCopy[i] = 0;
			SynTranslocateArrayIndex[i] = i;
		}
	}

	DescendSort(SynTranslocateArray,SynTranslocateArrayIndex,NoNe);

	int SynTranslocatePossible = 5;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	/*if ( countD > SynTranslocatePossible  )
	{
		if ( round(RN1*double(countD)) > SynTranslocatePossible ) //< was also tried, but didn't work 
			SynTranslocatePossible = round(RN1*double(countD));
	}
	else SynTranslocatePossible = countD;*/
	SynTranslocatePossible = ceil(RN1*double(countD));

	double threshold = SynTranslocateArray[SynTranslocatePossible];

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " in ENG is swaped by that of neuron " << neuron[1] << endl;	
	cout << "process_id: " << process_id << "\t SynTranslocatePossible = " << SynTranslocatePossible << endl;

	int SynTemp;
	// Column swaping (begin) within Excitatory NG
	for (int i=0; i<NoNe; i++)
	{
		if ( SynTranslocateArrayCopy[i]>threshold )
		{
			SynTemp = conMatEE[neuron[0]][i];
			conMatEE[neuron[0]][i] = conMatEE[neuron[1]][i];
			conMatEE[neuron[1]][i] = SynTemp;
			//cout << "process_id: " << process_id << "\t conMatEE[neuron[0]][i] = " << conMatEE[neuron[0]][i] << "\t conMatEE[neuron[1]][i] = " << conMatEE[neuron[1]][i] << endl;	
		}
	}
	
	delete [] neuron;
	delete [] SynTranslocateArray;
	delete [] SynTranslocateArrayCopy;
	delete [] SynTranslocateArrayIndex;

}


void Network::MutationNeuronSwitchING(int iter)
{
	/* pick 2 neurons within inhibitory NG and swap them; 
	 * equivalently, swap the incoming and outgoing synapses
	 * of first neuron with that of second, and vice versa */

	srand ( time(NULL)+process_id+iter*1000);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = NoNe+floor(double(NoNi)*RN1);
	}


	int SynTemp;

	// Row swaping (begin) within Inhibitory NG
	for (int i=0; i<NoNi; i++)
	{
		SynTemp = conMatII[neuron[0]-NoNe][i];
		conMatII[neuron[0]-NoNe][i] = conMatII[neuron[1]-NoNe][i];
		conMatII[neuron[1]-NoNe][i] = SynTemp;
	}				
	

	// Column swaping (begin) within Inhibitory NG
	for (int i=0; i<NoNi; i++)
	{
		SynTemp = conMatII[i][neuron[0]-NoNe];
		conMatII[i][neuron[0]-NoNe] = conMatII[i][neuron[1]-NoNe];
		conMatII[i][neuron[1]-NoNe] = SynTemp;
	}

	cout << "process_id: " << process_id << "\t neuron " << neuron[0] << " in ING is swaped by neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationColumnSwitchING(int iter)
{
	/* pick 2 neurons within inhibitory NG and swap their incoming synapses; 
	 * this is equivalent to swaping their corresponding columns in connectivity matrix */
	
	srand ( time(NULL)+process_id+iter*1000);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = NoNe+floor(double(NoNi)*RN1);
	}


	int SynTemp;

	// Column swaping (begin) within Inhibitory NG
	for (int i=0; i<NoNi; i++)
	{
		SynTemp = conMatII[i][neuron[0]-NoNe];
		conMatII[i][neuron[0]-NoNe] = conMatII[i][neuron[1]-NoNe];
		conMatII[i][neuron[1]-NoNe] = SynTemp;
	}

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " in ING is swaped by neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationRowSwitchING(int iter)
{
	/* pick 2 neurons within inhibitory NG and swap their outgoing synapses; 
	 * this is equivalent to swaping their corresponding rows in connectivity matrix */	
	
	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = NoNe+floor(double(NoNi)*RN1);
	}


	int SynTemp;

	// Row swaping (begin) within Inhibitory NG
	for (int i=0; i<NoNi; i++)
	{
		SynTemp = conMatII[neuron[0]-NoNe][i];
		conMatII[neuron[0]-NoNe][i] = conMatII[neuron[1]-NoNe][i];
		conMatII[neuron[1]-NoNe][i] = SynTemp;
	}
	
	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " in ING is swaped by neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}

void Network::MutationSynSwap(int iter, int EE)
{
	/* pick 2 synapses and swap them; neuron[0] and neuron[1] define the 
	 * first synapse and neuron[2] and neuron[3] define the second */

	srand ( time(NULL)+process_id+iter*1000);
	
	int* neuron;
	neuron = new int[4];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<4; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		if ( EE==1 ) neuron[k] = floor(double(NoNe)*RN1);
		else neuron[k] = floor(double(NoN)*RN1);
	}

	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	// structured mutation
	if ( RN1<0.1 )
	{
		neuron[2] = neuron[1];
		neuron[3] = neuron[0];
	}
	if ( RN1>=0.1 && RN1<0.2 ) neuron[3] = neuron[1];
	if ( RN1>=0.2 && RN1<0.3 ) neuron[2] = neuron[0];


	int SynTemp;

	// Row swaping (begin)
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			SynTemp = conMatII[neuron[0]-NoNe][neuron[1]-NoNe];
			if ( neuron[2]>=NoNe )
			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatII[neuron[0]-NoNe][neuron[1]-NoNe] = conMatII[neuron[2]-NoNe][neuron[3]-NoNe];
					conMatII[neuron[2]-NoNe][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatII[neuron[0]-NoNe][neuron[1]-NoNe] = conMatIE[neuron[2]-NoNe][neuron[3]];
					conMatIE[neuron[2]-NoNe][neuron[3]] = SynTemp;
				}
			}
			else // ( neuron[2]<NoNe )
			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatII[neuron[0]-NoNe][neuron[1]-NoNe] = conMatEI[neuron[2]][neuron[3]-NoNe];
					conMatEI[neuron[2]][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatII[neuron[0]-NoNe][neuron[1]-NoNe] = conMatEE[neuron[2]][neuron[3]];
					conMatEE[neuron[2]][neuron[3]] = SynTemp;
				}
			}
		}
		else // ( neuron[1]<NoNe )
		{
			SynTemp = conMatIE[neuron[0]-NoNe][neuron[1]];
			if ( neuron[2]>=NoNe )
			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatIE[neuron[0]-NoNe][neuron[1]] = conMatII[neuron[2]-NoNe][neuron[3]-NoNe];
					conMatII[neuron[2]-NoNe][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatIE[neuron[0]-NoNe][neuron[1]] = conMatIE[neuron[2]-NoNe][neuron[3]];
					conMatIE[neuron[2]-NoNe][neuron[3]] = SynTemp;
				}
			}
			else // ( neuron[2]<NoNe )
			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatIE[neuron[0]-NoNe][neuron[1]] = conMatEI[neuron[2]][neuron[3]-NoNe];
					conMatEI[neuron[2]][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatIE[neuron[0]-NoNe][neuron[1]] = conMatEE[neuron[2]][neuron[3]];
					conMatEE[neuron[2]][neuron[3]] = SynTemp;
				}
			}
		}		
	}
	else // ( neuron[0]<NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			SynTemp = conMatEI[neuron[0]][neuron[1]-NoNe];
			if ( neuron[2]>=NoNe )
			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatEI[neuron[0]][neuron[1]-NoNe] = conMatII[neuron[2]-NoNe][neuron[3]-NoNe];
					conMatII[neuron[2]-NoNe][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatEI[neuron[0]][neuron[1]-NoNe] = conMatIE[neuron[2]-NoNe][neuron[3]];
					conMatIE[neuron[2]-NoNe][neuron[3]] = SynTemp;
				}
			}
			else // ( neuron[2]<NoNe )

			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatEI[neuron[0]][neuron[1]-NoNe] = conMatEI[neuron[2]][neuron[3]-NoNe];
					conMatEI[neuron[2]][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatEI[neuron[0]][neuron[1]-NoNe] = conMatEE[neuron[2]][neuron[3]];
					conMatEE[neuron[2]][neuron[3]] = SynTemp;
				}
			}
		}
		else // ( neuron[1]<NoNe )
		{
			SynTemp = conMatEE[neuron[0]][neuron[1]];
			if ( neuron[2]>=NoNe )
			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatEE[neuron[0]][neuron[1]] = conMatII[neuron[2]-NoNe][neuron[3]-NoNe];
					conMatII[neuron[2]-NoNe][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatEE[neuron[0]][neuron[1]] = conMatIE[neuron[2]-NoNe][neuron[3]];
					conMatIE[neuron[2]-NoNe][neuron[3]] = SynTemp;
				}
			}
			else // ( neuron[2]<NoNe )
			{
				if ( neuron[3]>=NoNe ) 
				{
					conMatEE[neuron[0]][neuron[1]] = conMatEI[neuron[2]][neuron[3]-NoNe];

					conMatEI[neuron[2]][neuron[3]-NoNe] = SynTemp;
				}
				else
				{
					conMatEE[neuron[0]][neuron[1]] = conMatEE[neuron[2]][neuron[3]];
					conMatEE[neuron[2]][neuron[3]] = SynTemp;
				}
			}
		}		
	}

	

	cout << "process_id: " << process_id << "\t syn " << neuron[0] << "-" << neuron[1] << " swaped by syn " << neuron[2] << "-" << neuron[3] << endl;	
	
	delete [] neuron;	

}


void Network::MutationSynDeleteEE(int iter)
{
	/* pick 1 synapse randomly from conMatEE and delete it; 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1100);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = 0;
	while ( SynTemp == 0 )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			neuron[k] = floor(double(NoNe)*RN1);
		}
		
		SynTemp = conMatEE[neuron[0]][neuron[1]];
		if ( SynTemp>0 ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDel = ceil(RN1*double(SynTemp));
			conMatEE[neuron[0]][neuron[1]] = SynTemp-SynDel;				
		}
	}

	
	cout << "process_id: " << process_id << "\t Syn " << neuron[0] << "-" << neuron[1] << " reduced from " << SynTemp << " to " << SynTemp - SynDel << endl;	
	
	delete [] neuron;
	

}

void Network::MutationSynDeleteEI(int iter)
{
	/* pick 1 synapse randomly from conMatEI and delete it; 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1200);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = 0;
	while ( SynTemp == 0 )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			if ( k==0 ) neuron[k] = floor(double(NoNe)*RN1);
			else neuron[k] = floor(double(NoNi)*RN1);
		}
		
		SynTemp = conMatEI[neuron[0]][neuron[1]];
		if ( SynTemp>0 ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDel = ceil(RN1*double(SynTemp));
			conMatEI[neuron[0]][neuron[1]] = SynTemp-SynDel;				
		}
	}

	
	cout << "process_id: " << process_id << "\t Syn " << neuron[0] << "-" << neuron[1]+NoNe << " reduced from " << SynTemp << " to " << SynTemp - SynDel << endl;	
	
	delete [] neuron;
	

}


void Network::MutationSynDeleteIE(int iter)
{
	/* pick 1 synapse randomly from conMatIE and delete it; 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1300);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = 0;
	while ( SynTemp == 0 )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			if ( k==0 ) neuron[k] = floor(double(NoNi)*RN1);
			else neuron[k] = floor(double(NoNe)*RN1);
		}
		
		SynTemp = conMatIE[neuron[0]][neuron[1]];
		if ( SynTemp>0 ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDel = ceil(RN1*double(SynTemp));
			conMatIE[neuron[0]][neuron[1]] = SynTemp-SynDel;				
		}
	}

	
	cout << "process_id: " << process_id << "\t Syn " << neuron[0]+NoNe << "-" << neuron[1] << " reduced from " << SynTemp << " to " << SynTemp - SynDel << endl;	
	
	delete [] neuron;
	

}

void Network::MutationSynDeleteII(int iter)
{
	/* pick 1 synapse randomly from conMatII and delete it; 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1400);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = 0;
	while ( SynTemp == 0 )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			neuron[k] = floor(double(NoNi)*RN1);
		}
		
		SynTemp = conMatII[neuron[0]][neuron[1]];
		if ( SynTemp>0 ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynDel = ceil(RN1*double(SynTemp));
			conMatII[neuron[0]][neuron[1]] = SynTemp-SynDel;				
		}
	}

	
	cout << "process_id: " << process_id << "\t Syn " << neuron[0]+NoNe << "-" << neuron[1]+NoNe << " reduced from " << SynTemp << " to " << SynTemp - SynDel << endl;	
	
	delete [] neuron;
	

}


void Network::MutationSynAddEE(int iter)
{
	/* pick 1 synapse randomly from conMatEE and increase its synaptic 
	 * strength by random amount (MaxSynStrength is fixed); 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1010);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = synGradeEE;
	while ( SynTemp == synGradeEE )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			neuron[k] = floor(double(NoNe)*RN1);			
		}
		
		SynTemp = conMatEE[neuron[0]][neuron[1]];
		if ( SynTemp<synGradeEE ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAdd = ceil(RN1*double(synGradeEE-SynTemp)); 
			conMatEE[neuron[0]][neuron[1]] = SynTemp + SynAdd;
		}
	}

	
	cout << "process_id: " << process_id << "\t syn " << neuron[0] << "-" << neuron[1] << " increased from " << SynTemp << " to " << SynTemp + SynAdd << endl;	
	
	delete [] neuron;
	
	

}


void Network::MutationSynAddEI(int iter)
{
	/* pick 1 synapse randomly from conMatEI and increase its synaptic 
	 * strength by random amount (MaxSynStrength is fixed); 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1020);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = synGradeEI;
	while ( SynTemp == synGradeEI )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			if ( k==0 ) neuron[k] = floor(double(NoNe)*RN1);
			else neuron[k] = floor(double(NoNi)*RN1);			
		}
		
		SynTemp = conMatEI[neuron[0]][neuron[1]];
		if ( SynTemp<synGradeEI ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAdd = ceil(RN1*double(synGradeEI-SynTemp)); 
			conMatEI[neuron[0]][neuron[1]] = SynTemp + SynAdd;
		}
	}

	
	cout << "process_id: " << process_id << "\t syn " << neuron[0] << "-" << neuron[1]+NoNe << " increased from " << SynTemp << " to " << SynTemp + SynAdd << endl;	
	
	delete [] neuron;
	
	

}


void Network::MutationSynAddIE(int iter)
{
	/* pick 1 synapse randomly from conMatIE and increase its synaptic 
	 * strength by random amount (MaxSynStrength is fixed); 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1030);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = synGradeIE;
	while ( SynTemp == synGradeIE )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			if ( k==0 ) neuron[k] = floor(double(NoNi)*RN1);
			else neuron[k] = floor(double(NoNe)*RN1);			
		}
		
		SynTemp = conMatIE[neuron[0]][neuron[1]];
		if ( SynTemp<synGradeIE ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAdd = ceil(RN1*double(synGradeIE-SynTemp)); 
			conMatIE[neuron[0]][neuron[1]] = SynTemp + SynAdd;
		}
	}

	
	cout << "process_id: " << process_id << "\t syn " << neuron[0]+NoNe << "-" << neuron[1] << " increased from " << SynTemp << " to " << SynTemp + SynAdd << endl;	
	
	delete [] neuron;
	
	

}


void Network::MutationSynAddII(int iter)
{
	/* pick 1 synapse randomly from conMatII and increase its synaptic 
	 * strength by random amount (MaxSynStrength is fixed); 
	 * neuron[0] and neuron[1] define the synapse */

	srand ( time(NULL)+process_id+iter*1040);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	int SynTemp = synGradeII;
	while ( SynTemp == synGradeII )
	{
		for (int k=0; k<2; k++)
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;
			
			neuron[k] = floor(double(NoNi)*RN1);			
		}
		
		SynTemp = conMatII[neuron[0]][neuron[1]];
		if ( SynTemp<synGradeII ) 
		{
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			SynAdd = ceil(RN1*double(synGradeII-SynTemp)); 
			conMatII[neuron[0]][neuron[1]] = SynTemp + SynAdd;
		}
	}

	
	cout << "process_id: " << process_id << "\t syn " << neuron[0]+NoNe << "-" << neuron[1]+NoNe << " increased from " << SynTemp << " to " << SynTemp + SynAdd << endl;	
	
	delete [] neuron;
	
	

}

void Network::SumConMat(int MatTag, bool Pre_Mutation)
{
	int M, N;
	if ( MatTag==1 || MatTag==2 ) M = NoNe; 
	else M = NoNi;

	if ( MatTag==1 || MatTag==3 ) N = NoNe; 
	else N = NoNi;

	
	if ( MatTag == 1 )
	{
		EEsum = 0;
		for ( int i=0; i<M ; i++ ) 
		{
			for (int j=0; j<N ; j++) 
			{			
				EEsum = EEsum + conMatEE[i][j];
			}
		}
		if (Pre_Mutation) cout << "process_id: " << process_id << "\t EEsum_before =  " << EEsum << endl;
		else cout << "process_id: " << process_id << "\t EEsum_after =  " << EEsum << endl;
	}

	if ( MatTag == 2 )
	{
		EIsum = 0;
		for ( int i=0; i<M ; i++ ) 
		{
			for (int j=0; j<N ; j++) 
			{			
				EIsum = EIsum + conMatEI[i][j];
			}
		}
		if (Pre_Mutation) cout << "process_id: " << process_id << "\t EIsum_before =  " << EIsum << endl;
		else cout << "process_id: " << process_id << "\t EIsum_after =  " << EIsum << endl;
	}

	if ( MatTag == 3 )
	{
		IEsum = 0;
		for ( int i=0; i<M ; i++ ) 
		{
			for (int j=0; j<N ; j++) 
			{			
				IEsum = IEsum + conMatIE[i][j];
			}
		}
		if (Pre_Mutation) cout << "process_id: " << process_id << "\t IEsum_before =  " << IEsum << endl;
		else cout << "process_id: " << process_id << "\t IEsum_after =  " << IEsum << endl;
	}

	if ( MatTag == 4 )
	{
		IIsum = 0;
		for ( int i=0; i<M ; i++ ) 
		{
			for (int j=0; j<N ; j++) 
			{			
				IIsum = IIsum + conMatII[i][j];
			}
		}
		if (Pre_Mutation) cout << "process_id: " << process_id << "\t IIsum_before =  " << IIsum << endl;
		else cout << "process_id: " << process_id << "\t IIsum_after =  " << IIsum << endl;
	}
	
}


void Network::MutationPartialRowTranslocate(int iter)
{
	/* pick 2 neurons at random and (partially) translocates outgoing 
	 * synapse of neuron1 to neuron0; this is equivalent to 
	 * adding row elements of neuron1 to that of neuron0 in connectivity  
	 * matrix; row elements of neuron1 reduces to 0 (if P=1) */	
	// synapse from neuron1 to first M of neuron[2..4]
	// to be translocated to neuron0

	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[5];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<5; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);

		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		//neuron[k] = floor(double(NoN)*RN1);
		neuron[k] = floor(double(NoNe)*RN1); // for altering only conMatEE
		//if ( k>1 ) neuron[k] = NoNe + floor(double(NoNi)*RN1); // for preserving conMatEE
	}

	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);

	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	int M; // no. of neighbours of neuron1		
	if ( RN1<0.34 )	M = 1;
	else if ( RN1>0.66 ) M = 3;
	else M = 2;

	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);

	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	double P = 1.0; // level of synapse translocation (0.5 or 1)
	if ( RN1<0.5 ) P = 0.5;
	int SynTemp;

	
	// Row translocation (begin)
	for ( int k=0; k<M; k++ )
	{	
		if ( neuron[0]>=NoNe )
		{
			if ( neuron[1]>=NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatII[neuron[1]-NoNe][neuron[2+k]-NoNe]));
					conMatII[neuron[0]-NoNe][neuron[2+k]-NoNe] = conMatII[neuron[0]-NoNe][neuron[2+k]-NoNe] + SynTemp;
					conMatII[neuron[1]-NoNe][neuron[2+k]-NoNe] = conMatII[neuron[1]-NoNe][neuron[2+k]-NoNe] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatIE[neuron[1]-NoNe][neuron[2+k]]));				
					conMatIE[neuron[0]-NoNe][neuron[2+k]] = conMatIE[neuron[0]-NoNe][neuron[2+k]] + SynTemp;
					conMatIE[neuron[1]-NoNe][neuron[2+k]] = conMatIE[neuron[1]-NoNe][neuron[2+k]] - SynTemp;
				}
			}
			else // ( neuron[1]<NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatEI[neuron[1]][neuron[2+k]-NoNe]));
					conMatII[neuron[0]-NoNe][neuron[2+k]-NoNe] = conMatII[neuron[0]-NoNe][neuron[2+k]-NoNe] + SynTemp;
					conMatEI[neuron[1]][neuron[2+k]-NoNe] = conMatII[neuron[1]][neuron[2+k]-NoNe] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatEE[neuron[1]][neuron[2+k]]));				
					conMatIE[neuron[0]-NoNe][neuron[2+k]] = conMatIE[neuron[0]-NoNe][neuron[2+k]] + SynTemp;
					conMatEE[neuron[1]][neuron[2+k]] = conMatEE[neuron[1]][neuron[2+k]] - SynTemp;
				}
			}		
		}
		else // ( neuron[0]<NoNe )
		{
			if ( neuron[1]>=NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatII[neuron[1]-NoNe][neuron[2+k]-NoNe]));
					conMatEI[neuron[0]][neuron[2+k]-NoNe] = conMatEI[neuron[0]][neuron[2+k]-NoNe] + SynTemp;
					conMatII[neuron[1]-NoNe][neuron[2+k]-NoNe] = conMatII[neuron[1]-NoNe][neuron[2+k]-NoNe] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatIE[neuron[1]-NoNe][neuron[2+k]]));				
					conMatEE[neuron[0]][neuron[2+k]] = conMatEE[neuron[0]][neuron[2+k]] + SynTemp;
					conMatIE[neuron[1]-NoNe][neuron[2+k]] = conMatIE[neuron[1]-NoNe][neuron[2+k]] - SynTemp;
				}
			}
			else // ( neuron[1]<NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatEI[neuron[1]][neuron[2+k]-NoNe]));
					conMatEI[neuron[0]][neuron[2+k]-NoNe] = conMatEI[neuron[0]][neuron[2+k]-NoNe] + SynTemp;
					conMatEI[neuron[1]][neuron[2+k]-NoNe] = conMatEI[neuron[1]][neuron[2+k]-NoNe] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatEE[neuron[1]][neuron[2+k]]));				
					conMatEE[neuron[0]][neuron[2+k]] = conMatEE[neuron[0]][neuron[2+k]] + SynTemp;
					conMatEE[neuron[1]][neuron[2+k]] = conMatEE[neuron[1]][neuron[2+k]] - SynTemp;
				}
			}		
		}
	}	

	cout << "process_id: " << process_id << "\t Synapse from neuron " << neuron[1] << " to neurons " << neuron[2] << ", " << neuron[3] << " and " << neuron[4] << " translocated to neuron " << neuron[0] << "; M = " << M << " and P = " << P << endl;	
	
	delete [] neuron;	

}

void Network::MutationPartialColumnTranslocate(int iter)
{
	/* pick 2 neurons at random and (partially) translocate incoming 
	 * synapses of neuron1 to neuron0; this is equivalent to 
	 * adding column elements of neuron1 to that of neuron0 in connectivity  
	 * matrix; column elements of neuron1 reduces to 0 (if P=1) */	
	// incoming synapse to neuron1 from first M of neuron[2..4]
	// to be translocated to neuron0

	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[5];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<5; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);

		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		//neuron[k] = floor(double(NoN)*RN1);
		neuron[k] = floor(double(NoNe)*RN1); // for altering only conMatEE
		// if ( k<2 ) neuron[k] = NoNe + floor(double(NoNi)*RN1); // for preserving conMatEE
	}

	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);

	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	int M; // no. of neighbours of neuron1		
	if ( RN1<0.34 )	M = 1;
	else if ( RN1>0.66 ) M = 3;
	else M = 2;

	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);

	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	double P = 1.0; // level of synapse translocation (0.5 or 1)
	if ( RN1<0.5 ) P = 0.5;
	int SynTemp;

	
	// Column translocation (begin)
	for ( int k=0; k<M; k++ )
	{	
		if ( neuron[0]>=NoNe )
		{
			if ( neuron[1]>=NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatII[neuron[2+k]-NoNe][neuron[1]-NoNe]));
					conMatII[neuron[2+k]-NoNe][neuron[0]-NoNe] = conMatII[neuron[2+k]-NoNe][neuron[0]-NoNe] + SynTemp;
					conMatII[neuron[2+k]-NoNe][neuron[1]-NoNe] = conMatII[neuron[2+k]-NoNe][neuron[1]-NoNe] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatEI[neuron[2+k]][neuron[1]-NoNe]));	
					conMatEI[neuron[2+k]][neuron[0]-NoNe] = conMatEI[neuron[2+k]][neuron[0]-NoNe] + SynTemp;
					conMatEI[neuron[2+k]][neuron[1]-NoNe] = conMatEI[neuron[2+k]][neuron[1]-NoNe] - SynTemp;
				}		
			}
			else // ( neuron[1]<NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatIE[neuron[2+k]-NoNe][neuron[1]]));
					conMatII[neuron[2+k]-NoNe][neuron[0]-NoNe] = conMatII[neuron[2+k]-NoNe][neuron[0]-NoNe] + SynTemp;
					conMatIE[neuron[2+k]-NoNe][neuron[1]] = conMatIE[neuron[2+k]-NoNe][neuron[1]] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatEE[neuron[2+k]][neuron[1]]));
					conMatEI[neuron[2+k]][neuron[0]-NoNe] = conMatEI[neuron[2+k]][neuron[0]-NoNe] + SynTemp;
					conMatEE[neuron[2+k]][neuron[1]] = conMatEE[neuron[2+k]][neuron[1]] - SynTemp;
				}
			}		
		}
		else // ( neuron[0]<NoNe )
		{
			if ( neuron[1]>=NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatII[neuron[2+k]-NoNe][neuron[1]-NoNe]));
					conMatIE[neuron[2+k]-NoNe][neuron[0]] = conMatIE[neuron[2+k]-NoNe][neuron[0]] + SynTemp;
					conMatII[neuron[2+k]-NoNe][neuron[1]-NoNe] = conMatII[neuron[2+k]-NoNe][neuron[1]-NoNe] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatEI[neuron[2+k]][neuron[1]-NoNe]));
					conMatEE[neuron[2+k]][neuron[0]] = conMatEE[neuron[2+k]][neuron[0]] + SynTemp;
					conMatEI[neuron[2+k]][neuron[1]-NoNe] = conMatEI[neuron[2+k]][neuron[1]-NoNe] - SynTemp;
				}
			}
			else // ( neuron[1]<NoNe )
			{
				if ( neuron[2+k]>=NoNe )
				{	
					SynTemp = ceil(P*double(conMatIE[neuron[2+k]-NoNe][neuron[1]]));
					conMatIE[neuron[2+k]-NoNe][neuron[0]] = conMatIE[neuron[2+k]-NoNe][neuron[0]] + SynTemp;
					conMatIE[neuron[2+k]-NoNe][neuron[1]] = conMatIE[neuron[2+k]-NoNe][neuron[1]] - SynTemp;
				}
				else
				{
					SynTemp = ceil(P*double(conMatEE[neuron[2+k]][neuron[1]]));
					conMatEE[neuron[2+k]][neuron[0]] = conMatEE[neuron[2+k]][neuron[0]] + SynTemp;
					conMatEE[neuron[2+k]][neuron[1]] = conMatEE[neuron[2+k]][neuron[1]] - SynTemp;
				}
			}		
		}
	}	

	cout << "process_id: " << process_id << "\t Incoming synapse to neuron " << neuron[1] << " from neurons " << neuron[2] << ", " << neuron[3] << " and " << neuron[4] << " translocated to neuron " << neuron[0] << "; M = " << M << " and P = " << P << endl;	
	
	delete [] neuron;	

}

void Network::MutationRowTranslocate(int iter)
{
	/* pick 2 neurons at random and translocate outgoing synapses
	 * of neuron1 to neuron0; this is equivalent to 
	 * adding row elements of neuron1 to that of neuron0 in connectivity  
	 * matrix; row elements of neuron1 reduces to 0 */
	 
	 srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;


	double RN1;
	for (int k=0; k<2; k++)
	{

		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);

		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		//neuron[k] = floor(double(NoN)*RN1);
		neuron[k] = floor(double(NoNe)*RN1);// for altering only conMatEE
	}

	
	// Row translocation (begin)
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{				
				conMatIE[neuron[0]-NoNe][i] = conMatIE[neuron[0]-NoNe][i] + conMatIE[neuron[1]-NoNe][i];
				conMatIE[neuron[1]-NoNe][i] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{				
				conMatII[neuron[0]-NoNe][i] = conMatII[neuron[0]-NoNe][i] + conMatII[neuron[1]-NoNe][i];
				conMatII[neuron[1]-NoNe][i] = 0;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				conMatIE[neuron[0]-NoNe][i] = conMatIE[neuron[0]-NoNe][i] + conMatEE[neuron[1]][i];
				conMatEE[neuron[1]][i] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{				
				conMatII[neuron[0]-NoNe][i] = conMatII[neuron[0]-NoNe][i] + conMatEI[neuron[1]][i];
				conMatEI[neuron[1]][i] = 0;
			}
		}		
	}
	else // ( neuron[0]<NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{	
				conMatEE[neuron[0]][i] = conMatEE[neuron[0]][i] + conMatIE[neuron[1]-NoNe][i];
				conMatIE[neuron[1]-NoNe][i] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{
				conMatEI[neuron[0]][i] = conMatEI[neuron[0]][i] + conMatII[neuron[1]-NoNe][i];
				conMatII[neuron[1]-NoNe][i] = 0;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{				
				conMatEE[neuron[0]][i] = conMatEE[neuron[0]][i] + conMatEE[neuron[1]][i];
				conMatEE[neuron[1]][i] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{
				conMatEI[neuron[0]][i] = conMatEI[neuron[0]][i] + conMatEI[neuron[1]][i];
				conMatEI[neuron[1]][i] = 0;
			}
		}		

	}	

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[1] << " translocated to neuron " << neuron[0] << endl;	

	
	delete [] neuron;	

}


void Network::MutationColumnTranslocate(int iter)
{
	/* pick 2 neurons at random and translocate incoming synapses
	 * of neuron1 to neuron0; this is equivalent to 
	 * adding column elements of neuron1 to that of neuron0 in connectivity  
	 * matrix; column elements of neuron1 reduces to 0 */
	 
	 srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		//neuron[k] = floor(double(NoN)*RN1);
		neuron[k] = floor(double(NoNe)*RN1);// for altering only conMatEE

	}

	
	// Column translocation (begin)
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				conMatEI[i][neuron[0]-NoNe] = conMatEI[i][neuron[0]-NoNe] + conMatEI[i][neuron[1]-NoNe];
				conMatEI[i][neuron[1]-NoNe] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{
				conMatII[i][neuron[0]-NoNe] = conMatII[i][neuron[0]-NoNe] + conMatII[i][neuron[1]-NoNe];
				conMatII[i][neuron[1]-NoNe] = 0;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				conMatEI[i][neuron[0]-NoNe] = conMatEI[i][neuron[0]-NoNe] + conMatEE[i][neuron[1]];
				conMatEE[i][neuron[1]] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{
				conMatII[i][neuron[0]-NoNe] = conMatII[i][neuron[0]-NoNe] + conMatIE[i][neuron[1]];
				conMatIE[i][neuron[1]] = 0;
			}
		}		
	}
	else // ( neuron[0]<NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				conMatEE[i][neuron[0]] = conMatEE[i][neuron[0]] + conMatEI[i][neuron[1]-NoNe];
				conMatEI[i][neuron[1]-NoNe] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{
				conMatIE[i][neuron[0]] = conMatIE[i][neuron[0]] + conMatII[i][neuron[1]-NoNe];
				conMatII[i][neuron[1]-NoNe] = 0;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				conMatEE[i][neuron[0]] = conMatEE[i][neuron[0]] + conMatEE[i][neuron[1]];
				conMatEE[i][neuron[1]] = 0;
			}
			for (int i=0; i<NoNi; i++)
			{
				conMatIE[i][neuron[0]] = conMatIE[i][neuron[0]] + conMatIE[i][neuron[1]];
				conMatIE[i][neuron[1]] = 0;
			}
		}		
	}
		

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[1] << " translocated to neuron " << neuron[0] << endl;	
	
	delete [] neuron;	

}


void Network::MutationRowSwitch(int iter)
{
	/* pick 2 neurons at random and sawp their outgoing synapses;
	 * this is equivalent to swaping their corresponding rows in connectivity matrix */	

	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		//neuron[k] = floor(double(NoN)*RN1);
		neuron[k] = floor(double(NoNe)*RN1); // for altering only conMatEE
	}


	int SynTemp;

	// Row swaping (begin)
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatIE[neuron[0]-NoNe][i];
				conMatIE[neuron[0]-NoNe][i] = conMatIE[neuron[1]-NoNe][i];
				conMatIE[neuron[1]-NoNe][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[neuron[0]-NoNe][i];
				conMatII[neuron[0]-NoNe][i] = conMatII[neuron[1]-NoNe][i];
				conMatII[neuron[1]-NoNe][i] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatIE[neuron[0]-NoNe][i];
				conMatIE[neuron[0]-NoNe][i] = conMatEE[neuron[1]][i];
				conMatEE[neuron[1]][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[neuron[0]-NoNe][i];
				conMatII[neuron[0]-NoNe][i] = conMatEI[neuron[1]][i];
				conMatEI[neuron[1]][i] = SynTemp;
			}
		}		
	}
	else // ( neuron[0]<NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[neuron[0]][i];
				conMatEE[neuron[0]][i] = conMatIE[neuron[1]-NoNe][i];
				conMatIE[neuron[1]-NoNe][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatEI[neuron[0]][i];
				conMatEI[neuron[0]][i] = conMatII[neuron[1]-NoNe][i];
				conMatII[neuron[1]-NoNe][i] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[neuron[0]][i];
				conMatEE[neuron[0]][i] = conMatEE[neuron[1]][i];
				conMatEE[neuron[1]][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatEI[neuron[0]][i];
				conMatEI[neuron[0]][i] = conMatEI[neuron[1]][i];
				conMatEI[neuron[1]][i] = SynTemp;
			}
		}		
	}	

	cout << "process_id: " << process_id << "\t Row of neuron " << neuron[0] << " swaped by that of neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationColumnSwitch(int iter)
{
	/* pick 2 neurons at random and sawp their incoming synapses;
	 * this is equivalent to swaping their corresponding columns in connectivity matrix */	

	srand ( time(NULL)+process_id+iter*1000);

	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;


	double RN1;
	for (int k=0; k<2; k++)
	{

		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);

		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		//neuron[k] = floor(double(NoN)*RN1);
		neuron[k] = floor(double(NoNe)*RN1); // for altering only conMatEE

	}


	int SynTemp;

	// Column swaping (begin)
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEI[i][neuron[0]-NoNe];
				conMatEI[i][neuron[0]-NoNe] = conMatEI[i][neuron[1]-NoNe];
				conMatEI[i][neuron[1]-NoNe] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[i][neuron[0]-NoNe];
				conMatII[i][neuron[0]-NoNe] = conMatII[i][neuron[1]-NoNe];
				conMatII[i][neuron[1]-NoNe] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEI[i][neuron[0]-NoNe];
				conMatEI[i][neuron[0]-NoNe] = conMatEE[i][neuron[1]];
				conMatEE[i][neuron[1]] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[i][neuron[0]-NoNe];
				conMatII[i][neuron[0]-NoNe] = conMatIE[i][neuron[1]];
				conMatIE[i][neuron[1]] = SynTemp;
			}
		}		
	}
	else // ( neuron[0]<NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[i][neuron[0]];
				conMatEE[i][neuron[0]] = conMatEI[i][neuron[1]-NoNe];
				conMatEI[i][neuron[1]-NoNe] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatIE[i][neuron[0]];
				conMatIE[i][neuron[0]] = conMatII[i][neuron[1]-NoNe];
				conMatII[i][neuron[1]-NoNe] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[i][neuron[0]];
				conMatEE[i][neuron[0]] = conMatEE[i][neuron[1]];
				conMatEE[i][neuron[1]] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatIE[i][neuron[0]];
				conMatIE[i][neuron[0]] = conMatIE[i][neuron[1]];
				conMatIE[i][neuron[1]] = SynTemp;
			}
		}		
	}		

	cout << "process_id: " << process_id << "\t Column of neuron " << neuron[0] << " swaped by that of neuron " << neuron[1] << endl;	
	
	delete [] neuron;	

}


void Network::MutationNeuronSwitch(int iter)
{
	/* pick 2 neurons and swap them; equivalently, swap the incoming and 
	 * outgoing synapses of first neuron with that of second, and vice versa */

	srand ( time(NULL)+process_id+iter*1000);
	
	int* neuron;
	neuron = new int[2];
	int RN1_copy = 0;
	int RN2_copy = 0;

	double RN1;
	for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		//neuron[k] = floor(double(NoN)*RN1);
		neuron[k] = floor(double(NoNe)*RN1); // for altering only conMatEE
	}


	int SynTemp;

	// Row swaping (begin)
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatIE[neuron[0]-NoNe][i];
				conMatIE[neuron[0]-NoNe][i] = conMatIE[neuron[1]-NoNe][i];
				conMatIE[neuron[1]-NoNe][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[neuron[0]-NoNe][i];
				conMatII[neuron[0]-NoNe][i] = conMatII[neuron[1]-NoNe][i];
				conMatII[neuron[1]-NoNe][i] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatIE[neuron[0]-NoNe][i];
				conMatIE[neuron[0]-NoNe][i] = conMatEE[neuron[1]][i];
				conMatEE[neuron[1]][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[neuron[0]-NoNe][i];
				conMatII[neuron[0]-NoNe][i] = conMatEI[neuron[1]][i];
				conMatEI[neuron[1]][i] = SynTemp;
			}
		}		
	}
	else // ( neuron[0]<NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[neuron[0]][i];
				conMatEE[neuron[0]][i] = conMatIE[neuron[1]-NoNe][i];
				conMatIE[neuron[1]-NoNe][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatEI[neuron[0]][i];
				conMatEI[neuron[0]][i] = conMatII[neuron[1]-NoNe][i];
				conMatII[neuron[1]-NoNe][i] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[neuron[0]][i];
				conMatEE[neuron[0]][i] = conMatEE[neuron[1]][i];
				conMatEE[neuron[1]][i] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatEI[neuron[0]][i];
				conMatEI[neuron[0]][i] = conMatEI[neuron[1]][i];
				conMatEI[neuron[1]][i] = SynTemp;
			}
		}		
	}

	// Column swaping (begin)
	if ( neuron[0]>=NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEI[i][neuron[0]-NoNe];
				conMatEI[i][neuron[0]-NoNe] = conMatEI[i][neuron[1]-NoNe];
				conMatEI[i][neuron[1]-NoNe] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[i][neuron[0]-NoNe];
				conMatII[i][neuron[0]-NoNe] = conMatII[i][neuron[1]-NoNe];
				conMatII[i][neuron[1]-NoNe] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEI[i][neuron[0]-NoNe];
				conMatEI[i][neuron[0]-NoNe] = conMatEE[i][neuron[1]];
				conMatEE[i][neuron[1]] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatII[i][neuron[0]-NoNe];
				conMatII[i][neuron[0]-NoNe] = conMatIE[i][neuron[1]];
				conMatIE[i][neuron[1]] = SynTemp;
			}
		}		
	}
	else // ( neuron[0]<NoNe )
	{
		if ( neuron[1]>=NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[i][neuron[0]];
				conMatEE[i][neuron[0]] = conMatEI[i][neuron[1]-NoNe];
				conMatEI[i][neuron[1]-NoNe] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatIE[i][neuron[0]];
				conMatIE[i][neuron[0]] = conMatII[i][neuron[1]-NoNe];
				conMatII[i][neuron[1]-NoNe] = SynTemp;
			}
		}
		else // ( neuron[1]<NoNe )
		{
			for (int i=0; i<NoNe; i++)
			{
				SynTemp = conMatEE[i][neuron[0]];
				conMatEE[i][neuron[0]] = conMatEE[i][neuron[1]];
				conMatEE[i][neuron[1]] = SynTemp;
			}
			for (int i=0; i<NoNi; i++)
			{
				SynTemp = conMatIE[i][neuron[0]];
				conMatIE[i][neuron[0]] = conMatIE[i][neuron[1]];
				conMatIE[i][neuron[1]] = SynTemp;
			}
		}		
	}

	cout << "process_id: " << process_id << "\t neuron " << neuron[0] << " swaped by neuron " << neuron[1] << endl;	
	
	delete [] neuron;

}


void Network::ConnectMutationSegG(int iter)
{
	/* pick 4 neurons as centers and define their neighbors (N)
	 * delete connections from neurons with center N[0] to neurons with center N[1]
	 * Add connections from neurons with center N[2] to neurons with center N[3] */	
	
	srand ( time(NULL)+process_id+iter*1000);

	//int NoN = NoNe + NoNi;
	int* neuron;
	neuron = new int[4];
	int* N;
	N = new int[4]; // no. of neurons in neighborhood of neuron[4]	
	int* Ng;
	Ng = new int[4]; // no. of neurons in network containing neuron[4]
	int NeighborMaxE = round(0.10*double(NoNe));
	int NeighborMaxI = round(0.10*double(NoNi));
	int RN1_copy = 0;
	int RN2_copy = 0;
	
	double RN1;
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;

	NeighborMaxE = round(RN1*NeighborMaxE);
	NeighborMaxI = round(RN1*NeighborMaxI);
	/*if (RN1<0.5) // small mutation/neighbor
	{
		NeighborMaxE = round(RN1*NeighborMaxE);
		NeighborMaxI = round(RN1*NeighborMaxI);
	}*/
	// else // large mutation/neighbor

	// for effecting mutations in EI and IE synapses
	/*for (int k=0; k<2; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if (RN1>0.5)
		{
			neuron[2*k] = floor(double(NoNe)*RN1);
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			neuron[2*k+1] = NoNe + floor(double(NoNi)*RN1);
		}
		else
		{
			neuron[2*k] = NoNe + floor(double(NoNi)*RN1);
			do
			{
				RN1 = double(rand()) / double(RAND_MAX);
				RN2_copy = int(RN1*100000);
			} while (RN1_copy == RN2_copy);
			RN1_copy = RN2_copy;

			neuron[2*k+1] = floor(double(NoNe)*RN1);
		}	
		
	}


    for (int k=0; k<4; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if ( neuron[k]<NoNe )
		{
			N[k] = ceil(double(NeighborMaxE)*RN1);
			Ng[k] = NoNe;
		}
		else 
		{
			N[k] = ceil(double(NeighborMaxI)*RN1);
			Ng[k] = NoNi;
		}
	}


	// choosing nature of mutation by redefining neuron[k] and N[k]
	// *************************************************************
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	if ( RN1<0.2 ) // flip local inSyn to outSyn 
	{
		neuron[2] = neuron[1];
		neuron[3] = neuron[0];
		N[2] = N[1];
		N[3] = N[0];
		Ng[2] = Ng[1];
		Ng[3] = Ng[0];
	}
	if ( RN1>=0.2 && RN1<0.4 ) // displace outSyn from neighborhood of neuron[0] 
	{
		if (neuron[0]<NoNe)
		{
			if (neuron[2]<NoNe)
			{
				neuron[2] = neuron[0];
				N[2] = N[0];
				Ng[2] = Ng[0];
			}
			else
			{
				neuron[2] = neuron[1];
				N[2] = N[1];
				Ng[2] = Ng[1];
			}			
		}
		else // if (neuron[0]>NoNe)
		{
			if (neuron[2]>NoNe)
			{
				neuron[2] = neuron[0];
				N[2] = N[0];
				Ng[2] = Ng[0];
			}
			else
			{
				neuron[2] = neuron[1];
				N[2] = N[1];
				Ng[2] = Ng[1];
			}			
		}
	}
	if ( RN1>=0.4 && RN1<0.6 ) // displace inSyn from neighborhood of neuron[1] 
	{
		if (neuron[1]<NoNe)
		{
			if (neuron[3]<NoNe)
			{
				neuron[3] = neuron[1];
				N[3] = N[1];
				Ng[3] = Ng[1];
			}
			else
			{
				neuron[3] = neuron[0];
				N[3] = N[0];
				Ng[3] = Ng[0];
			}			
		}
		else // if (neuron[1]>NoNe)
		{
			if (neuron[3]>NoNe)
			{
				neuron[3] = neuron[1];
				N[3] = N[1];
				Ng[3] = Ng[1];
			}
			else
			{
				neuron[3] = neuron[0];
				N[3] = N[0];
				Ng[3] = Ng[0];
			}			
		}
	}		
	// if RN1>=0.6, then usual mutation comences
	// **************************************************************** 
	*/

	//for effecting mutations in EE synapse only
	for (int k=0; k<4; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		neuron[k] = floor(double(NoNe)*RN1);
	}

	for (int k=0; k<4; k++)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;

		if ( neuron[k]<NoNe )	
		{
			N[k] = 1;
			//N[k] = ceil(double(NeighborMaxE)*RN1);
			Ng[k] = NoNe;
		}
		else 
		{
			N[k] = 1;
			//N[k] = ceil(double(NeighborMaxI)*RN1);
			Ng[k] = NoNi;
		}
	}


	// choosing nature of mutation by redefining neuron[k] and N[k]
	// *************************************************************
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	if ( RN1<0.2 ) // flip local inSyn to outSyn 
	{
		neuron[2] = neuron[1];
		neuron[3] = neuron[0];
		N[2] = N[1];
		N[3] = N[0];
		Ng[2] = Ng[1];
		Ng[3] = Ng[0];
	}
	if ( RN1>=0.2 && RN1<0.4 ) // displace outSyn from neighborhood of neuron[0] 
	{
		neuron[2] = neuron[0];
		N[2] = N[0];
		Ng[2] = Ng[0];
	}
	if ( RN1>=0.4 && RN1<0.6 ) // displace inSyn from neighborhood of neuron[1] 
	{
		neuron[3] = neuron[1];
		N[3] = N[1];
		Ng[3] = Ng[1];
	}
	// if RN1>=0.6, then usual mutation comences
	// **************************************************************** 

	
	for (int k=0; k<4; k++)
	{
		cout << "process_id: " << process_id << "\t neuron[" << k << "] = " << neuron[k] << "\t N[" << k << "] = " <<  N[k] << "\t Ng[" << k << "] = " << Ng[k] << endl;		
	}

	double* dist0;
	int* dist0index;
	dist0 = new double[Ng[0]];	
	dist0index = new int[Ng[0]];

	double* dist1;
	int* dist1index;
	dist1 = new double[Ng[1]];	
	dist1index = new int[Ng[1]];

	double* dist2;
	int* dist2index;
	dist2 = new double[Ng[2]];	
	dist2index = new int[Ng[2]];

	double* dist3;
	int* dist3index;
	dist3 = new double[Ng[3]];	
	dist3index = new int[Ng[3]];

		
	for (int k=0; k<Ng[0]; k++)
	{
		if (Ng[0]==NoNe) dist0[k] = pow( (gLeake[k]-gLeake[neuron[0]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[0]]),2 );
		else dist0[k] = pow( (gLeaki[k]-gLeaki[neuron[0]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[0]-NoNe]),2 );
		dist0index[k] = k;
	}
	

	for (int k=0; k<Ng[1]; k++)
	{
		if (Ng[1]==NoNe) dist1[k] = pow( (gLeake[k]-gLeake[neuron[1]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[1]]),2 );
		else dist1[k] = pow( (gLeaki[k]-gLeaki[neuron[1]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[1]-NoNe]),2 );
		dist1index[k] = k;		
	}
	

	for (int k=0; k<Ng[2]; k++)
	{
		if (Ng[2]==NoNe) dist2[k] = pow( (gLeake[k]-gLeake[neuron[2]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[2]]),2 );
		else dist2[k] = pow( (gLeaki[k]-gLeaki[neuron[2]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[2]-NoNe]),2 );
		dist2index[k] = k;
	}
	

	for (int k=0; k<Ng[3]; k++)
	{
		if (Ng[3]==NoNe) dist3[k] = pow( (gLeake[k]-gLeake[neuron[3]]),2 ) + pow( (gNaPe[k]-gNaPe[neuron[3]]),2 );
		else dist3[k] = pow( (gLeaki[k]-gLeaki[neuron[3]-NoNe]),2 ) + pow( (gNaPi[k]-gNaPi[neuron[3]-NoNe]),2 );
		dist3index[k] = k;			
	}
	

	DescendSort(dist0, dist0index, Ng[0]);
	DescendSort(dist1, dist1index, Ng[1]);
	DescendSort(dist2, dist2index, Ng[2]);
	DescendSort(dist3, dist3index, Ng[3]);
	

	int SynDeletePossible = 0;
	double* SynDeleteArray;
	SynDeleteArray = new double[N[0]*N[1]*(synGrade+1)];
	double* SynDeleteArrayCopy;
	SynDeleteArrayCopy = new double[N[0]*N[1]*(synGrade+1)];
	int* SynDeleteArrayIndex;
	SynDeleteArrayIndex = new int[N[0]*N[1]*(synGrade+1)];
	int countD = 0;	
	if ( neuron[0]>NoNe )
	{
		if ( neuron[1]>NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					SynDeletePossible = SynDeletePossible + conMatII[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]];
					for (int g=0; g<=synGrade; g++)
					{
						if (g<conMatII[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynDeleteArray[countD] = RN1;
							SynDeleteArrayCopy[countD] = RN1;
							SynDeleteArrayIndex[countD] = countD;
						}
						else
						{
							SynDeleteArray[countD] = 0;
							SynDeleteArrayCopy[countD] = 0;
							SynDeleteArrayIndex[countD] = countD;
						}					
						countD = countD + 1;
					}
				}
			}			
		}
		else
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					SynDeletePossible = SynDeletePossible + conMatIE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]];
					for (int g=0; g<=synGrade; g++)
					{
						if (g<conMatIE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynDeleteArray[countD] = RN1;
							SynDeleteArrayCopy[countD] = RN1;
							SynDeleteArrayIndex[countD] = countD;
						}
						else
						{
							SynDeleteArray[countD] = 0;
							SynDeleteArrayCopy[countD] = 0;
							SynDeleteArrayIndex[countD] = countD;
						}					
						countD = countD + 1;
					}
				}
			}
		}
	}
	else
	{
		if ( neuron[1]>NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					SynDeletePossible = SynDeletePossible + conMatEI[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]];
					for (int g=0; g<=synGrade; g++)
					{
						if (g<conMatEI[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynDeleteArray[countD] = RN1;
							SynDeleteArrayCopy[countD] = RN1;
							SynDeleteArrayIndex[countD] = countD;
						}
						else
						{
							SynDeleteArray[countD] = 0;
							SynDeleteArrayCopy[countD] = 0;
							SynDeleteArrayIndex[countD] = countD;
						}					
						countD = countD + 1;
					}
				}
			}
		}
		else
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					SynDeletePossible = SynDeletePossible + conMatEE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]];
					for (int g=0; g<=synGrade; g++)
					{
						if (g<conMatEE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynDeleteArray[countD] = RN1;
							SynDeleteArrayCopy[countD] = RN1;
							SynDeleteArrayIndex[countD] = countD;
						}
						else
						{
							SynDeleteArray[countD] = 0;
							SynDeleteArrayCopy[countD] = 0;
							SynDeleteArrayIndex[countD] = countD;
						}					
						countD = countD + 1;
					}
				}
			}
		}
	}
	
	cout << "process_id: " << process_id << "\t SynDeletePossible = " << SynDeletePossible << endl;	

	
	int SynAdditionPossible = 0;
	double* SynAddArray;
	SynAddArray = new double[N[2]*N[3]*(synGrade+1)];
	double* SynAddArrayCopy;
	SynAddArrayCopy = new double[N[2]*N[3]*(synGrade+1)];
	int* SynAddArrayIndex;
	SynAddArrayIndex = new int[N[2]*N[3]*(synGrade+1)];
	int countA = 0;
	if ( neuron[2]>NoNe )
	{
		if ( neuron[3]>NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					SynAdditionPossible = SynAdditionPossible + (synGrade-conMatII[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]);
					for (int g=0; g<=synGrade; g++)
					{
						if (g>conMatII[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynAddArray[countA] = RN1;
							SynAddArrayCopy[countA] = RN1;
							SynAddArrayIndex[countA] = countA;
						}
						else
						{
							SynAddArray[countA] = 0;
							SynAddArrayCopy[countA] = 0;
							SynAddArrayIndex[countA] = countA;
						}
						countA = countA + 1;
					}
				}
			}		
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					SynAdditionPossible = SynAdditionPossible + (synGrade-conMatIE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]);
					for (int g=0; g<=synGrade; g++)
					{
						if (g>conMatIE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynAddArray[countA] = RN1;
							SynAddArrayCopy[countA] = RN1;
							SynAddArrayIndex[countA] = countA;
						}
						else
						{
							SynAddArray[countA] = 0;
							SynAddArrayCopy[countA] = 0;
							SynAddArrayIndex[countA] = countA;
						}
						countA = countA + 1;
					}
				}
			}
		}
	}
	else	
	{
		if ( neuron[3]>NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					SynAdditionPossible = SynAdditionPossible + (synGrade-conMatEI[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]);
					for (int g=0; g<=synGrade; g++)
					{
						if (g>conMatEI[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynAddArray[countA] = RN1;
							SynAddArrayCopy[countA] = RN1;
							SynAddArrayIndex[countA] = countA;
						}
						else
						{
							SynAddArray[countA] = 0;
							SynAddArrayCopy[countA] = 0;
							SynAddArrayIndex[countA] = countA;
						}
						countA = countA + 1;
					}
				}
			}
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					SynAdditionPossible = SynAdditionPossible + (synGrade-conMatEE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]]);
					for (int g=0; g<=synGrade; g++)
					{
						if (g>conMatEE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]])
						{
							do
							{
								RN1 = double(rand()) / double(RAND_MAX);
								RN2_copy = int(RN1*100000);
							} while (RN1_copy == RN2_copy);
							RN1_copy = RN2_copy;

							SynAddArray[countA] = RN1;
							SynAddArrayCopy[countA] = RN1;
							SynAddArrayIndex[countA] = countA;
						}
						else
						{
							SynAddArray[countA] = 0;
							SynAddArrayCopy[countA] = 0;
							SynAddArrayIndex[countA] = countA;
						}
						countA = countA + 1;
					}
				}
			}
		}
	}
	
	cout << "process_id: " << process_id << "\t SynAdditionPossible = " << SynAdditionPossible << endl;

	int SynMutatePossible = SynDeletePossible;
	if (SynAdditionPossible<SynDeletePossible) SynMutatePossible = SynAdditionPossible;
	cout << "process_id: " << process_id << "\t SynMutatePossible = " << SynMutatePossible << endl;

	int MinSynMutate = 5*(synGrade);
	do
	{
		RN1 = double(rand()) / double(RAND_MAX);
		RN2_copy = int(RN1*100000);
	} while (RN1_copy == RN2_copy);
	RN1_copy = RN2_copy;
	if (RN1<0.5) // small mutation/depth
	{
		MinSynMutate = 5*(0.4*synGrade);
	}
	// else // large mutation/depth
	
	int SynMutate;
	if (SynMutatePossible<=MinSynMutate) SynMutate = SynMutatePossible;
	else // (SynMutatePossible>MinSynMutate)
	{
		do
		{
			RN1 = double(rand()) / double(RAND_MAX);
			RN2_copy = int(RN1*100000);
		} while (RN1_copy == RN2_copy);
		RN1_copy = RN2_copy;
		
		SynMutate = round(RN1*double(SynMutatePossible));
		if (SynMutate<=MinSynMutate) SynMutate = MinSynMutate;
	}
	cout << "process_id: " << process_id << "\t SynMutate = " << SynMutate << endl;

	DescendSort(SynDeleteArray, SynDeleteArrayIndex, countD);
	DescendSort(SynAddArray, SynAddArrayIndex, countA);
	
	
	double SynDeleteThreshold = SynDeleteArray[SynMutate-1];
	double SynAddThreshold = SynAddArray[SynMutate-1];

	int countDD = 0;
	if ( neuron[0]>NoNe )
	{
		if ( neuron[1]>NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
						{
							conMatII[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = conMatII[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] - 1;
						}
						countDD = countDD + 1;
					}					
				}
			}
		}
		else
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
						{
							conMatIE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = conMatIE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] - 1;
						}
						countDD = countDD + 1;
					}
				}
			}
		}
	}
	else
	{
		if ( neuron[1]>NoNe )
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
						{
							conMatEI[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = conMatEI[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] - 1;
						}
						countDD = countDD + 1;
					}
				}
			}
		}
		else
		{
			for (int i=0; i<N[0]; i++)
			{
				for (int j=0; j<N[1]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynDeleteArrayCopy[countDD]>=SynDeleteThreshold)
						{
							conMatEE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] = conMatEE[dist0index[Ng[0]-1-i]][dist1index[Ng[1]-1-j]] - 1;
						}
						countDD = countDD + 1;
					}
				}
			}
		}
	}

	
	int countAA = 0;
	if ( neuron[2]>NoNe )
	{
		if ( neuron[3]>NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynAddArrayCopy[countAA]>=SynAddThreshold)
						{
							conMatII[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = conMatII[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] + 1;
						}
						countAA = countAA + 1;
					}
				}
			}
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynAddArrayCopy[countAA]>=SynAddThreshold)
						{
							conMatIE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = conMatIE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] + 1;
						}
						countAA = countAA + 1;
					}
				}
			}
		}
	}
	else
	{
		if ( neuron[3]>NoNe )
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynAddArrayCopy[countAA]>=SynAddThreshold)
						{
							conMatEI[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = conMatEI[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] + 1;
						}
						countAA = countAA + 1;
					}
				}
			}
		}
		else
		{
			for (int i=0; i<N[2]; i++)
			{
				for (int j=0; j<N[3]; j++)
				{
					for (int g=0; g<=synGrade; g++)
					{
						if (SynAddArrayCopy[countAA]>=SynAddThreshold)
						{
							conMatEE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] = conMatEE[dist2index[Ng[2]-1-i]][dist3index[Ng[3]-1-j]] + 1;
						}
						countAA = countAA + 1;
					}
				}
			}
		}
	}	

	
	
	delete [] SynDeleteArray;
	delete [] SynDeleteArrayCopy;
	delete [] SynDeleteArrayIndex;
	delete [] SynAddArray;
	delete [] SynAddArrayCopy;
	delete [] SynAddArrayIndex;

	delete [] dist0index;
	delete [] dist0;
	delete [] dist1index;
	delete [] dist1;
	delete [] dist2index;
	delete [] dist2;
	delete [] dist3index;
	delete [] dist3;

	delete [] neuron;
	delete [] N;
	delete [] Ng;	

}


void Network::SmallWorldEE(int n0, double p)
{
	// n0 = no. of neighbours of each node
	// p = probability with which each edge is randomly altered
	int nLeft = int(floor(double(n0)/2));
	int nRight = n0-nLeft;
	cout << "nLeft = " << nLeft << endl;
	cout << "nRight = " << nRight << endl;
	
	int M = NoNe;
	int N = NoNe;	
		
	for ( int i=0; i<M ; i++ ) 
	{
		if ( i<nLeft ) 
		{
			for (int j=0; j<N ; j++)
			{
				if ( j>i+nRight && j<N+(i-nLeft) ) conMatEE[i][j] = 0;
				else conMatEE[i][j] = 1;

				if ( j==i ) conMatEE[i][j] = 0;
			}				
		}
		else if ( i+nRight>N-1 ) 
		{
			for (int j=0; j<N ; j++)
			{
				if ( j>i+nRight-(N) && j<i-nLeft ) conMatEE[i][j] = 0;
				else conMatEE[i][j] = 1;

				if ( j==i ) conMatEE[i][j] = 0;
			}				
		}
		else
		{
			for (int j=0; j<N ; j++)
			{
				if ( j>=i-nLeft  && j<=i+nRight ) conMatEE[i][j] = 1;
				else conMatEE[i][j] = 0;

				if ( j==i ) conMatEE[i][j] = 0;
			}
		}
	}
	
	SmallWorldEE2(n0,p);

}


void Network::SmallWorldEE2(int n0, double p)
{
	
	int RN1_copy = 0;
	int RN2_copy = 0;
	
	double* RanMatD;
	RanMatD = new double[NoNe*NoNe];
	double* RanMatDCopy;
	RanMatDCopy = new double[NoNe*NoNe];
	int* RanDIndex;
	RanDIndex = new int[NoNe*NoNe];

	double RN1;
	//srand( time(NULL)+10101 );
	for (int i=0; i<NoNe; i++)
	{
		for (int j=0; j<NoNe; j++)
		{
			if (conMatEE[i][j]==1)
			{
				do
				{
					RN1 = (double)rand( ) / RAND_MAX;
					RN2_copy = (int) (RN1*100000);
				} while (RN1_copy == RN2_copy);
				RN1_copy = RN2_copy;
				RanMatD[i*NoNe+j] = RN1;
				RanMatDCopy[i*NoNe+j] = RN1;
			}
			else
			{
				RanMatD[i*NoNe+j] = 0;
				RanMatDCopy[i*NoNe+j] = 0;
			}
			RanDIndex[i*NoNe+j] = i*NoNe+j;
		}
	}
	
	DescendSort(RanMatD,RanDIndex,NoNe*NoNe);
	int SynMutate = int(  round(p*double(NoNe*n0)) );
	cout << "SynMutate(SW) = " << SynMutate << endl;

	/*for ( int i=0; i<NoNe*NoNe; i++)
		cout << i << "\t" << RanMatD[i] << "\t" << RanMatDCopy[i] << "\t" << RanDIndex[i] << endl;*/

	//cout << RanMatD[SynMutate] << endl;

	for (int i=0; i<NoNe; i++)
	{
		for (int j=0; j<NoNe; j++) 
		{	
			conMatEEcopy[i][j] = conMatEE[i][j];
			//cout << "conMatEEcopy[" << i << "][" << j << "] = " << conMatEEcopy[i][j] << endl;
		}
	}
	

	double* RanMatA;
	RanMatA = new double[NoNe];
	double* RanMatACopy;
	RanMatACopy = new double[NoNe];
	int* RanAIndex;
	RanAIndex = new int[NoNe];
	for (int i=0; i<NoNe; i++)
	{
		for (int j=0; j<NoNe; j++)
		{
			//cout << conMatEEcopy[i][j] << "\t" << RanMatDCopy[i*NoNe+j] << endl;
			if (conMatEEcopy[i][j]==1 && RanMatDCopy[i*NoNe+j]>RanMatD[SynMutate])
			{
				//cout << i << "," << j << "\t" << RanMatDCopy[i*NoNe+j] << endl;				
				for (int k=0; k<NoNe; k++)
				{
					if ( conMatEE[i][k]==1 || k==i ) //conMatEE is used instead of conMatEEcopy because conMatEE is the the updated during synaptic reconfiguration 
					{
						RanMatA[k] = 0;
						RanMatACopy[k] = 0;
					}
					else
					{	
						do
						{
							RN1 = (double)rand( ) / RAND_MAX;
							RN2_copy = (int) (RN1*100000);
						} while (RN1_copy == RN2_copy);
						RN1_copy = RN2_copy;
						RanMatA[k] = RN1;
						RanMatACopy[k] = RN1;
					}
					RanAIndex[k] = k;
				}
				
				DescendSort(RanMatA,RanAIndex,NoNe);
				
				for (int k=0; k<NoNe; k++)
				{
					if (RanMatACopy[k]>RanMatA[1]) 
					{
						conMatEE[i][k]=1;
						//cout << "conMatEEcopy[" << i << "][" << j << "] =" <<  conMatEEcopy[i][j] << "\t conMatEE[" << i << "][" << k << "] = " << conMatEE[i][k] << endl;
					}
				}
				conMatEE[i][j]=0;
			}			
		}
	}
	delete [] RanMatA;
	delete [] RanMatACopy;
	delete [] RanAIndex;
	
	delete [] RanMatD;
	delete [] RanMatDCopy;
	delete [] RanDIndex;
	
}






}
