#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iomanip>
using namespace ::std;
//This structure is created to help decode RM codes with List
struct InformationBits {
	vector <int> InfoBits;
	vector < vector <double> > ScratchPad;
	double LogProbability;
};
/////////////////////////////////
//Basic Functions (In Common with previous designs)
void showD(vector <double> v);
void show(vector <int> v);
void show2D(vector <vector<int> > v);
void showD2D(vector <vector<double> > v);
double RandGenAWGN(double sigma);
vector <int> NumtoSign(vector <double> input);
vector < vector <int> > NumtoSign2D(vector< vector <double> >input);
double comparetwo(vector <int> a,vector <int> b,double M);
double comparetwo2D(vector< vector <int> >a,vector <vector <int> >b,double N, double M);//N*M
////////////////////////////////
vector <int> Information_Creator(int k, double condition);

vector <vector <int> > Parity_Check_Creator(int k,int d,vector <vector <int> > information_bits);
vector <double> EPS_Information_Bits(int k,vector <int> information_bits, double Sigma);
vector <vector <double> > EPS_Parity_Check_Bits(int k, int d,vector <vector <int> > parity_check_bits, double Sigma);
vector <vector <double> > SDDecoder2(int k, int d, vector <vector <double> > eps_information_bits,vector <vector <double> > eps_parity_check_bits, int num_iteration);
vector <vector <double> > SDDecoder(int k, int d, vector <vector <double> > eps_information_bits,vector <vector <double> > eps_parity_check_bits, int num_iteration); 
/////////////////////////////
//Polar Precoding 
vector <vector <int> > RMListDecoder(vector <double> EpsRec,int Log2M, int r,int Log2L,vector <int> FrozenBitsPosition);
vector <vector <double> > RightCollapse(vector <vector<double> > input,int Log2M,int row,int StartPoint);
vector <int> PrefixDetector(int NInfoBitDecodedSofar,int Log2M);
vector <int> RMPlotkinEncoder(vector <int> Information,int M,int r);
vector <vector <double> > SortedArray(vector <vector <double> >Input,int log2M );
vector <int> BinaryRep(int N,int Digits);
int BinaryToint(vector <int> v);
int RMInfoBitCounter(int M,int R);
vector <int> RMTableChecker(vector <vector <int> > ListofEstimates,int Log2M, int r,int Log2L,vector <double> EpsRec);
vector <int> FrozenBitsPolar(int k1);
vector <int> Polar_information_bit_creator(int k1,vector <int> Frozen_bits_polar, double condition);
vector <int> RMTableChecker_fullcode(vector <vector <int> > ListofEstimates,int Log2M, int r,int Log2L,vector <double> EpsRec
,vector <double> PCEpsRec,vector < vector <int> > ParityCheckPositions);
int precoding_info_bits(int k1,vector <int> Frozen_bits_polar);

/////////////////////////////
vector <int> Encoder1D (vector <int> InfoBit,vector < vector <int> > ParityCheckPositions);
vector < vector <int> > PC_Position_Creator(int k,int d_2, int d_4);
vector <double> PCEps1D (vector <int> InfoBit,vector <int> paritychecks, vector < vector <int> > ParityCheckPositions, double Sigma);
vector <double> BPSDDecoder(int k,int d_2, int d_4,vector <double> eps_info_bit,vector <double> paritychecks_eps, int num_iteration,
	vector < vector <int> > ParityCheckPositions);

int main()
{
srand (time(NULL));	
freopen ("BPSD plus polar precoding with parity checks of weight 2 and 4 k=64 d=138 n=4096 1-2-1-2022.txt","a",stdout);

int k=64;
int d_2=100;//64
int d_4=32;//96
int n=k*d_4/4+k*d_2/2;
int log2L=5;
int num_rep_info_bits=6;
double condition=0.50;
int num_iterations=6*log2(k);//2*log2(k);
int Trial=100;
vector < vector <int> > pc_pos= PC_Position_Creator(k,d_2,d_4);
//cout<<"\n"<<pc_pos.size()<<"	"<<pc_pos[0].size()<<"\n";
//show2D(pc_pos);
vector <int> weight_PC (k,0);
for(int i=0;i<n;i++)
{
	int size_w=pc_pos[i].size();
	for(int j=0;j<size_w;j++)
	{
		weight_PC[pc_pos[i][j]]++;
	}
}
vector <int> Frozen_Bits_Polar = FrozenBitsPolar(k);
cout<<"\nInformation Bits positions are:\n";
show(Frozen_Bits_Polar);
int precoding_info_bits1=precoding_info_bits(k,Frozen_Bits_Polar);
double rate2=precoding_info_bits1/double(k);
double rate1=k/double(n+k*num_rep_info_bits);
double Initial_SNR=0.0000;
double Final_SNR=3.6000;
double SNR_Step=0.5;//in dB
double SNR=pow(10,Initial_SNR/10);//Initial SNR
vector <vector <int> > Histogram1(100,vector <int> (k+1,0));
vector <vector <int> > Histogram2(100,vector <int> (k+1,0));
vector <vector <int> > Histogram3(100,vector <int> (k+1,0));
cout<<"\nBPSD Error (Gallager) C(n,k)";
cout<<"\nNumber of information bits of precoding (k) : "<<precoding_info_bits1;
cout<<"\nNumber of information bits of gallager design (k2) : "<<k;
cout<<"\nSize of the entire code(n) : "<<n+k*num_rep_info_bits;
cout<<"\nDistance : "<<d_2+d_4+1*num_rep_info_bits<<"	w1="<<num_rep_info_bits<<"	w2="<<d_2<<"	w4="<<d_4;
cout<<"\nprecoding rate : "<<rate2;
cout<<"\nBP rate = "<<rate1;
cout<<"\nrate = "<<rate1*rate2;
cout<<"\nNumberof iterations : "<<num_iterations;
cout<<"\nNumber of times information bits were repeated (r) : "<<num_rep_info_bits;
cout<<"\nNumber of trials : "<<Trial;

int Histindex=0;
while(SNR<pow(10,(Final_SNR/10))+0.0001)
{
	double SNRindB=floor(100*log10(SNR)+0.01)/10;
	double output_BER=0;
	double input_BER=0;
	double polar_BER=0;
	double Sigma=1/sqrt(2*SNR*rate1*rate2);
	for(int t=0;t<Trial;t++)
	{	
		vector <int> polar_info_bit1=Polar_information_bit_creator(k,Frozen_Bits_Polar,condition);
		vector <int> information_bits=RMPlotkinEncoder(polar_info_bit1,log2(k),log2(k));		
		//vector <int> information_bits = Information_Creator (k,condition);
		vector <int> paritychecks = Encoder1D(information_bits,pc_pos);
		vector <double> eps_info_bit = EPS_Information_Bits(k,information_bits,Sigma/sqrt(num_rep_info_bits));
		vector <double> paritychecks_eps = PCEps1D(information_bits,paritychecks,pc_pos,Sigma);
		vector <double> eps_decoded = BPSDDecoder(k,d_2,d_4,eps_info_bit,paritychecks_eps,num_iterations,pc_pos);
		vector <vector <int> > list_of_decoded_polar_bits=RMListDecoder(eps_decoded,log2(k),log2(k),log2L,Frozen_Bits_Polar);
		vector <int> decoded_bits=RMTableChecker_fullcode(list_of_decoded_polar_bits,log2(k),log2(k),log2L,eps_info_bit,
		paritychecks_eps,pc_pos);
		//cout<<"\n";
		//showD(eps_info_bit);
		//showD(eps_decoded);
		double error_0=comparetwo(information_bits,NumtoSign(eps_info_bit),k);
		double error_1=comparetwo(information_bits,NumtoSign(eps_decoded),k);
		double error_2=comparetwo(information_bits,decoded_bits,k);
		
		input_BER+=error_0/double(k);
		output_BER+=error_1/double(k);
		polar_BER+=error_2/double(k);
		Histogram1[Histindex][error_0]++;
		Histogram2[Histindex][error_1]++;
		Histogram3[Histindex][error_2]++;
	}
	cout<<"\nSNR is :	"<<setprecision(3)<<SNRindB;
	cout<<"	Initial BER :	"<<setprecision(3)<<input_BER/Trial;
	cout<<"	WER :	"<<setprecision(3)<<1-Histogram3[Histindex][0]/double(Trial);	
	cout<<"	BP BER :	"<<setprecision(3)<<output_BER/Trial;
	cout<<"	Polar BER :	"<<setprecision(3)<<polar_BER/Trial;
	cout<<"		Sigma = "<<setprecision(3)<<Sigma;
	Histindex++;
	SNR=SNR*pow(10,SNR_Step/10);///SNR Steps
}
cout<<"\n\n";
for(int i=0;i<Histindex;i++)
{
	cout<<"\nSNR is : "<<Initial_SNR+i*SNR_Step;
	cout<<"\nInitial Error :\n";
	show(Histogram1[i]);
	cout<<"\nBPSD Error :\n";
	show(Histogram2[i]);
	cout<<"\nPolar Error :\n";
	show(Histogram3[i]);	
}

//show(weight_PC);

//showD(paritychecks_eps);
//showD(eps_info_bit);
//showD(eps_decoded);
//show(infromation_bits);
//show(paritychecks);

fclose (stdout);
return 0;
}


vector <vector <double> > SDDecoder(int k, int d, vector <vector <double> > eps_information_bits,vector <vector <double> > eps_parity_check_bits, int num_iteration)	 
{
	int k1=k/2;
	vector < vector < vector <double> > > eps_current_stage (d+1,vector<vector<double> > (2,vector <double> (k1,0)));
	for(int i=0;i<d+1;i++)
	{
		eps_current_stage[i]=eps_information_bits;
	}
	for (int iteration=0;iteration<num_iteration;iteration++ )
	{
		vector < vector <vector <double> > > LLR_next_stage (d+1,vector<vector<double> > (2,vector <double> (k1,0)));
		for(int i1=0;i1<d;i1++)
		{
			for(int i2=0;i2<k1;i2++)
			{
				double Temp1=eps_current_stage[i1][0][i2]*eps_parity_check_bits[i1][i2];
				double Temp2=eps_current_stage[i1][1][(i1+i2)%k1]*eps_parity_check_bits[i1][i2];
				LLR_next_stage[d][1][(i1+i2)%k1]+=log((1+Temp1)/(1-Temp1));
				LLR_next_stage[d][0][i2]+=log((1+Temp2)/(1-Temp2));
				
				LLR_next_stage[i1][1][(i1+i2)%k1]=log((1+Temp1)/(1-Temp1));
				LLR_next_stage[i1][0][i2]=log((1+Temp2)/(1-Temp2));
				//output[i1][i2]=information_bits[0][i2]*information_bits[1][(i1+i2)%k1];	
			}
		}
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<k1;i2++)
			{
				LLR_next_stage[d][i1][i2]+=log((1+eps_information_bits[i1][i2])/(1-eps_information_bits[i1][i2]));
			}
		}
		for(int j=0;j<d+1;j++)
		{
			for(int i1=0;i1<2;i1++)
			{
				for(int i2=0;i2<k1;i2++)
				{
					double Temp1=LLR_next_stage[d][i1][i2]-LLR_next_stage[j][i1][i2];
					eps_current_stage[j][i1][i2]=(exp(Temp1)-1)/(exp(Temp1)+1);
					if(j==d)
					{
						double Temp2=LLR_next_stage[d][i1][i2];
						eps_current_stage[d][i1][i2]=(exp(Temp2)-1)/(exp(Temp2)+1);
					}
				}
			}
		}
	}
	//showD2D(eps_current_stage);
	return eps_current_stage[d];
}

double comparetwo2D(vector< vector <int> >a,vector <vector <int> >b,double N, double M)//N*M
{
	double difference=0.0;
	for(int i=0;i<N;i++)
	{
		difference+=comparetwo(a[i],b[i],M);
	}
	//difference=difference;
	return difference;
}

vector < vector <int> > NumtoSign2D(vector< vector <double> >input)
{
	int row=input.size();
	int col=input[0].size();
	vector< vector <int> >output(row,vector <int> (col,0));
	for(int i=0;i<row;i++)
	{
		output[i]=NumtoSign(input[i]);
	}
	return output;
}
double comparetwo(vector <int> a,vector <int> b,double M)
{
	double difference=0.0;
	for(int i=0;i<M;i++)
	{
		if(a[i]!=b[i])
		{
			difference+=1;
		}
	}
	//difference=difference;
	return difference;
}
vector <int> NumtoSign(vector <double> input)
{
	vector <int> output(input.size(),0);
	for(int i=0;i<output.size();i++)
	{
		if(input[i]>0)
		{
			output[i]=1;
		}
		else
		{
			output[i]=-1;
		}
	}
	return output;
}
double RandGenAWGN(double sigma)
{
	double a,b,c;
	a=rand();
	a=(a+0.0001)/(0.0002+RAND_MAX);
	b=rand();
	b=(b/RAND_MAX);
	c=sqrt((-2)*log(a))*cos(2*(3.14159265)*b);
	if(c>10)
	{
		c=10;
		cout<<"\n+infinity";
	}
	if(c<-10)
	{
		c=-10;
		cout<<"\n-infinity";		
	}
	c=c*sigma;
	return c;
}
void showD(vector <double> v)
{
	int size1=v.size();
	for(int i=0;i<size1;i++)
	{
		cout<<v[i]<<" ";
	}
	cout<<endl;
}
void show(vector <int> v)
{
	int size1=v.size();
	for(int i=0;i<size1;i++)
	{
		cout<<v[i]<<" ";
	}
	cout<<endl;
}
void show2D(vector <vector<int> > v)
{
	cout<<"\n";
	for(int i=0;i<v.size();i++)
	{
		show(v[i]);
	}
}
void showD2D(vector <vector<double> > v)
{
	cout<<"\n";
	for(int i=0;i<v.size();i++)
	{
		showD(v[i]);
	}
}
/////////////////////////////////////////////
vector <vector <double> > SDDecoder2(int k, int d, vector <vector <double> > eps_information_bits,
	 vector <vector <double> > eps_parity_check_bits, int num_iteration)	 
{
	int k1=k/2;
	vector <vector <double> > eps_current_stage=eps_information_bits;
	for (int iteration=0;iteration<num_iteration;iteration++ )
	{
		vector <vector <double> > LLR_next_stage(2,vector <double> (k1,0));
		for(int i1=0;i1<d;i1++)
		{
			for(int i2=0;i2<k1;i2++)
			{
				double Temp1=eps_current_stage[0][i2]*eps_parity_check_bits[i1][i2];
				double Temp2=eps_current_stage[1][(i1+i2)%k1]*eps_parity_check_bits[i1][i2];
				LLR_next_stage[1][(i1+i2)%k1]+=log((1+Temp1)/(1-Temp1));
				LLR_next_stage[0][i2]+=log((1+Temp2)/(1-Temp2));
				//output[i1][i2]=information_bits[0][i2]*information_bits[1][(i1+i2)%k1];	
			}
		}
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<k1;i2++)
			{
				LLR_next_stage[i1][i2]+=log((1+eps_information_bits[i1][i2])/(1-eps_information_bits[i1][i2]));
				eps_current_stage[i1][i2]=(exp(LLR_next_stage[i1][i2])-1)/(exp(LLR_next_stage[i1][i2])+1);
			}
		}
		
	}
	//showD2D(eps_current_stage);
	return eps_current_stage;
}
vector <int> Information_Creator(int k, double condition)
{

	vector <int> output (k,1);
	for(int i1=0;i1<k;i1++)
	{
		double a=(rand()+0.001)/(RAND_MAX+0.001);
		if(a<condition)
		{
				output[i1]=-1;
		}
	}
	return output;
}

vector <vector <int> > Parity_Check_Creator(int k,int d,vector <vector <int> > information_bits)
{
	int k1=k/2;
	vector <vector <int> > output (d,vector <int> (k1,0));
	for(int i1=0;i1<d;i1++)
	{
		for(int i2=0;i2<k1;i2++)
		{
			
			output[i1][i2]=information_bits[0][i2]*information_bits[1][(i1+i2)%k1];	
		}
	}
	return output;
		
}

vector <double>  EPS_Information_Bits(int k, vector <int>  information_bits, double Sigma)
{
	vector <double>  output (k,0.0);
	for(int i1=0;i1<k;i1++)
	{
		double a=RandGenAWGN(Sigma) ;
		double Temp1=2*(information_bits[i1]+a)/(Sigma*Sigma);
		double Temp2=(exp(Temp1)-1)/(exp(Temp1)+1);
		output[i1]=Temp2;
	}
	return output;
}

vector <vector <double> > EPS_Parity_Check_Bits(int k, int d,vector <vector <int> > parity_check_bits, double Sigma)
{
	int k1=k/2;
	vector <vector <double> > output (d,vector <double> (k1,0));
	for(int i1=0;i1<d;i1++)
	{
		for(int i2=0;i2<k1;i2++)
		{	
			double a=RandGenAWGN(Sigma);
			double Temp1=2*(parity_check_bits[i1][i2]+a)/(Sigma*Sigma);
			double Temp2=(exp(Temp1)-1)/(exp(Temp1)+1);
			output[i1][i2]=Temp2;		
		}
	}
	return output;
}
vector <int> Encoder1D (vector <int> InfoBit,vector < vector <int> > ParityCheckPositions)
{
	int N = ParityCheckPositions.size();
	vector <int> output (N,1);
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<ParityCheckPositions[i].size();j++)
		{
			output[i]=output[i]*InfoBit[ParityCheckPositions[i][j]];
		}
	}
	return output;
}
vector < vector <int> > PC_Position_Creator(int k,int d_2, int d_4)
{
	int n=k*d_4/4+k*d_2/2;
	vector < vector <int> > pc_pos (n,vector <int> (2,0));
	int index=0;
	for(int i=0;i<k;i++)
	{
		for(int j=0;j<d_2/2;j++)
		{
			pc_pos[index][0]=i;
			pc_pos[index][1]=(i+1+j)%k;
			index++;
		}
	}
	vector <int> Temp_PC_4 (4,0);
	for(int i=0;i<k/4;i++)
	{
		
		for(int j=0;j<d_4;j++)
		{
			int dis_index=(j+1)/(k/4);
			Temp_PC_4[0]=4*i;
			Temp_PC_4[1]=(4*(i+j)+1)%k;
			Temp_PC_4[2]=(4*(i+(3+4*dis_index)*(j+1))+2)%k;
			Temp_PC_4[3]=(4*(i+(5+4*dis_index)*(j+1))+3)%k;
			pc_pos[index]=Temp_PC_4;
			index++;
		}	
	}
	return pc_pos;
}

vector <double> PCEps1D (vector <int> InfoBit,vector <int> paritychecks, vector < vector <int> > ParityCheckPositions, double Sigma)
{
	int n= paritychecks.size();
	vector <double> output (n,1);
	for(int i1=0;i1<n;i1++)
	{
		double a=RandGenAWGN(Sigma) ;
		double Temp1=2*(paritychecks[i1]+a)/(Sigma*Sigma);
		double Temp2=(exp(Temp1)-1)/(exp(Temp1)+1);
		output[i1]=Temp2;
	}
	return output;
}
vector <double> BPSDDecoder(int k,int d_2, int d_4,vector <double> eps_info_bit,vector <double> paritychecks_eps, int num_iteration,
	vector < vector <int> > ParityCheckPositions)
{
	//int n=k*d_4/4+k*d_2/2;
	int n=paritychecks_eps.size();
	vector <double> output (k,0);
	vector < vector <double> > Eps_Next(k,vector <double> (k,0));
	
	for(int i=0;i<k;i++)
	{
		Eps_Next[i]=eps_info_bit;
	}
	for(int iteration=0;iteration<num_iteration;iteration++)
	{
		vector < vector <double> > Partial_llr(k,vector <double> (k,0));
		for(int i=0;i<n;i++)
		{
			if(ParityCheckPositions[i].size()==2)
			{
				double temp1=0;
				double temp2=0;
				temp1=paritychecks_eps[i]*Eps_Next[ParityCheckPositions[i][0]][ParityCheckPositions[i][1]];
				temp2=log((1+temp1)/(1-temp1));
				Partial_llr[ParityCheckPositions[i][1]][ParityCheckPositions[i][0]]+=temp2;
				Partial_llr[ParityCheckPositions[i][1]][ParityCheckPositions[i][1]]+=temp2;

				temp1=paritychecks_eps[i]*Eps_Next[ParityCheckPositions[i][1]][ParityCheckPositions[i][0]];
				temp2=log((1+temp1)/(1-temp1));
				Partial_llr[ParityCheckPositions[i][0]][ParityCheckPositions[i][1]]+=temp2;
				Partial_llr[ParityCheckPositions[i][0]][ParityCheckPositions[i][0]]+=temp2;		
			}
			if(ParityCheckPositions[i].size()==4)
			{
				double temp1=0;
				double temp2=0;
				temp1=paritychecks_eps[i]*Eps_Next[ParityCheckPositions[i][1]][ParityCheckPositions[i][0]]*
				Eps_Next[ParityCheckPositions[i][2]][ParityCheckPositions[i][0]]*Eps_Next[ParityCheckPositions[i][3]][ParityCheckPositions[i][0]];
				temp2=log((1+temp1)/(1-temp1));				
				Partial_llr[ParityCheckPositions[i][0]][ParityCheckPositions[i][1]]+=temp2;
				Partial_llr[ParityCheckPositions[i][0]][ParityCheckPositions[i][2]]+=temp2;
				Partial_llr[ParityCheckPositions[i][0]][ParityCheckPositions[i][3]]+=temp2;				
				Partial_llr[ParityCheckPositions[i][0]][ParityCheckPositions[i][0]]+=temp2;

				temp1=paritychecks_eps[i]*Eps_Next[ParityCheckPositions[i][0]][ParityCheckPositions[i][1]]*
				Eps_Next[ParityCheckPositions[i][2]][ParityCheckPositions[i][1]]*Eps_Next[ParityCheckPositions[i][3]][ParityCheckPositions[i][1]];
				temp2=log((1+temp1)/(1-temp1));				
				Partial_llr[ParityCheckPositions[i][1]][ParityCheckPositions[i][0]]+=temp2;
				Partial_llr[ParityCheckPositions[i][1]][ParityCheckPositions[i][2]]+=temp2;
				Partial_llr[ParityCheckPositions[i][1]][ParityCheckPositions[i][3]]+=temp2;				
				Partial_llr[ParityCheckPositions[i][1]][ParityCheckPositions[i][1]]+=temp2;

				temp1=paritychecks_eps[i]*Eps_Next[ParityCheckPositions[i][0]][ParityCheckPositions[i][2]]*
				Eps_Next[ParityCheckPositions[i][1]][ParityCheckPositions[i][2]]*Eps_Next[ParityCheckPositions[i][3]][ParityCheckPositions[i][2]];
				temp2=log((1+temp1)/(1-temp1));				
				Partial_llr[ParityCheckPositions[i][2]][ParityCheckPositions[i][0]]+=temp2;
				Partial_llr[ParityCheckPositions[i][2]][ParityCheckPositions[i][1]]+=temp2;
				Partial_llr[ParityCheckPositions[i][2]][ParityCheckPositions[i][3]]+=temp2;				
				Partial_llr[ParityCheckPositions[i][2]][ParityCheckPositions[i][2]]+=temp2;

				temp1=paritychecks_eps[i]*Eps_Next[ParityCheckPositions[i][0]][ParityCheckPositions[i][3]]*
				Eps_Next[ParityCheckPositions[i][1]][ParityCheckPositions[i][3]]*Eps_Next[ParityCheckPositions[i][2]][ParityCheckPositions[i][3]];
				temp2=log((1+temp1)/(1-temp1));				
				Partial_llr[ParityCheckPositions[i][3]][ParityCheckPositions[i][0]]+=temp2;
				Partial_llr[ParityCheckPositions[i][3]][ParityCheckPositions[i][1]]+=temp2;
				Partial_llr[ParityCheckPositions[i][3]][ParityCheckPositions[i][2]]+=temp2;				
				Partial_llr[ParityCheckPositions[i][3]][ParityCheckPositions[i][3]]+=temp2;

			}
		}
		for(int i=0;i<k;i++)
		{
			Partial_llr[i][i]+=log((1+eps_info_bit[i])/(1-eps_info_bit[i]));
			Eps_Next[i][i]=(exp(Partial_llr[i][i])-1)/(exp(Partial_llr[i][i])+1);
			for(int j=0;j<k;j++)
			{
				if(j!=i)
				{
					double temp1=Partial_llr[i][i]-Partial_llr[i][j];
					Eps_Next[i][j]=(exp(temp1)-1)/(exp(temp1)+1);
				}
				
			}
		}
	}
	//showD2D(Partial_llr);
	for(int i=0;i<k;i++)
	{
		output[i]=Eps_Next[i][i];
	}
	return output;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Polar Codes
vector <vector <int> > RMListDecoder(vector <double> EpsRec,int Log2M, int r,int Log2L,vector <int> FrozenBitsPosition)
{
	int List=pow(2,Log2L);
	//
	//int List=pow(2,Log2L+1);
	//int List=pow(2,Log2L+2);
	//Log2L++;
	//Log2L++;
	//
	int M=pow(2,Log2M);
	vector <vector <double> > Mat1 (Log2M+1,vector <double> (M,0));
	Mat1[0]=EpsRec;
	Mat1=RightCollapse(Mat1,Log2M,0,0);
	vector <int> InitialInfo (M,0);
	InitialInfo[M-1]=1;
	double InitialProb=-10000;
	InformationBits InitialCase;
	InitialCase.ScratchPad=Mat1;
	InitialCase.InfoBits=InitialInfo;
	InitialCase.LogProbability=InitialProb;
	vector <InformationBits> output (List,InitialCase);
	if(FrozenBitsPosition[M-1]==1)
	{
		output[0].LogProbability=0;
	}
	else
	{
		double EpsilonFirstBit=Mat1[Log2M][M-1];
		output[0].LogProbability=log((1+EpsilonFirstBit)/2);
		output[1].LogProbability=log((1-EpsilonFirstBit)/2);
		output[1].InfoBits[M-1]=-1;
	}	
	for(int i=1;i<M;i++)
	{
		vector <int> Prefix=PrefixDetector(i,Log2M);
		int CommonPrefix=Prefix[0];					//out[0]=CommonPrefix;
		int Prfix=Prefix[1];						//out[1]=Prefix;
		for(int j=0;j<List;j++)
		{
			int SizeofBlock=pow(2,Log2M-CommonPrefix);
			int StartPoint=(pow(2,CommonPrefix)-1-Prfix)*SizeofBlock;
			int EndPoint=StartPoint+SizeofBlock-1;
			vector <int> TempInfoBit(SizeofBlock/2);
			for(int k=0;k<SizeofBlock/2;k++)
			{
				TempInfoBit[k]=output[j].InfoBits[StartPoint+k+SizeofBlock/2];
			}
			vector <int> TempCodeWord=RMPlotkinEncoder(TempInfoBit,Log2M-CommonPrefix-1,Log2M-CommonPrefix-1);
			for(int k=StartPoint;k<StartPoint+SizeofBlock/2;k++)
			{
				output[j].ScratchPad[CommonPrefix+1][k]=(output[j].ScratchPad[CommonPrefix][k]+output[j].ScratchPad[CommonPrefix][k+SizeofBlock/2]*TempCodeWord[k-StartPoint])/(1+(output[j].ScratchPad[CommonPrefix][k]*output[j].ScratchPad[CommonPrefix][k+SizeofBlock/2]*TempCodeWord[k-StartPoint]));	
			}
			int row=CommonPrefix+1;
			int initSP=StartPoint+pow(2,Log2M-row-1);
			if(row!=Log2M)//CommonPrefix+1
			{
				for(int k=row+1;k<Log2M+1;k++)
				{
					for(int q=initSP;q<initSP+pow(2,Log2M-k);q++)
					{
						output[j].ScratchPad[k][q]=output[j].ScratchPad[k-1][q]*output[j].ScratchPad[k-1][q-pow(2,Log2M-k)];
					}
					initSP=initSP+pow(2,Log2M-k-1);
				}
			}
		}
		vector <vector <double> > sortPositions(2,vector <double> (2*List));
		for(int k=0;k<List;k++)
		{
			double Epsilon=output[k].ScratchPad[Log2M][M-1-i];
			double LLRP=0;
			double LLRN=0;
			if(Epsilon==1)
			{
				
				LLRP=log((1+Epsilon)/2);
				LLRN=-20;
			}
			if(Epsilon==-1)
			{
				LLRP=-20;
				LLRN=log((1-Epsilon)/2);
			}	
			if((-1<Epsilon)&&(Epsilon<1))
			{
				LLRP=log((1+Epsilon)/2);
				LLRN=log((1-Epsilon)/2);
			}	
			if(FrozenBitsPosition[M-1-i]==1)
			{
				sortPositions[0][k]=output[k].LogProbability+LLRP;
				sortPositions[1][k]=k;
				sortPositions[0][k+List]=output[k].LogProbability-200000;
				sortPositions[1][k+List]=k+List;
			}
			if(FrozenBitsPosition[M-1-i]!=1)//in this case if FrozenBitsPosition[M-1-i]!=1 it should be decoded
			{
				sortPositions[0][k]=output[k].LogProbability+LLRP;
				sortPositions[1][k]=k;
				sortPositions[0][k+List]=output[k].LogProbability+LLRN;
				sortPositions[1][k+List]=k+List;
			}
			
		}
		sortPositions=SortedArray(sortPositions,Log2L+1);
		//cout<<"\nThis is Sorted Position\n";
		//showD2D(sortPositions);
		vector <int> NumofOcc(List,-1);
		
		for(int k=0;k<List;k++)
		{
			int counter=0;
			for(int h=0;h<List;h++)
			{
				if(sortPositions[1][h]==k)
				{
					counter++;
				}
				if(sortPositions[1][h]==k+List)
				{
					counter=counter+10;
				}
			}
			NumofOcc[k]=counter;
		}
		//cout<<"\nThis is Number of occurences : ";
		//show(NumofOcc);
		vector < vector <int> > HashArray1(2,vector <int> (List/2,-1));
		int index1=0;
		int index2=0;
		for(int k=0;k<List;k++)
		{
			if(NumofOcc[k]==0)
			{
				HashArray1[0][index1]=k;
				index1++;
			}
			if(NumofOcc[k]==11)
			{
				HashArray1[1][index2]=k;
				index2++;
			}
		}
		//cout<<"\nThis is Hash Array :\n";
		//show2D(HashArray1);
		for(int k=0;k<List;k++)
		{
			if(NumofOcc[k]==1)
			{
				output[k].InfoBits[M-1-i]=1;
				double Epsilon=output[k].ScratchPad[Log2M][M-1-i];
				double LLRP=0;
				double LLRN=0;
				if(Epsilon==1)
				{
					
					LLRP=log((1+Epsilon)/2);
					LLRN=-20;
				}
				if(Epsilon==-1)
				{
					LLRP=-20;
					LLRN=log((1-Epsilon)/2);
				}	
				if((-1<Epsilon)&&(Epsilon<1))
				{
					LLRP=log((1+Epsilon)/2);
					LLRN=log((1-Epsilon)/2);
				}	
				output[k].LogProbability=output[k].LogProbability+LLRP;			
			}
			if(NumofOcc[k]==10)
			{
				output[k].InfoBits[M-1-i]=-1;
				double Epsilon=output[k].ScratchPad[Log2M][M-1-i];
				double LLRP=0;
				double LLRN=0;
				if(Epsilon==1)
				{
					
					LLRP=log((1+Epsilon)/2);
					LLRN=-20;
				}
				if(Epsilon==-1)
				{
					LLRP=-20;
					LLRN=log((1-Epsilon)/2);
				}	
				if((-1<Epsilon)&&(Epsilon<1))
				{
					LLRP=log((1+Epsilon)/2);
					LLRN=log((1-Epsilon)/2);
				}	
				output[k].LogProbability=output[k].LogProbability+LLRN;			
			}
			
		}
		//cout<<"\nIndex1 : "<<index1;
		for(int k=0;k<index1;k++)
		{
			output[HashArray1[0][k]]=output[HashArray1[1][k]];
			output[HashArray1[0][k]].InfoBits[M-1-i]=1;
			output[HashArray1[1][k]].InfoBits[M-1-i]=-1;
			double Epsilon=output[HashArray1[1][k]].ScratchPad[Log2M][M-1-i];
			double LLRP=0;
			double LLRN=0;
			if(Epsilon==1)
			{
				
				LLRP=log((1+Epsilon)/2);
				LLRN=-20;
			}
			if(Epsilon==-1)
			{
				LLRP=-20;
				LLRN=log((1-Epsilon)/2);
			}	
			if((-1<Epsilon)&&(Epsilon<1))
			{
				LLRP=log((1+Epsilon)/2);
				LLRN=log((1-Epsilon)/2);
			}
			output[HashArray1[0][k]].LogProbability=output[HashArray1[0][k]].LogProbability+LLRP;
			output[HashArray1[1][k]].LogProbability=output[HashArray1[1][k]].LogProbability+LLRN;
		}	
	}
	vector <vector <int> > out1(List);
	vector <vector <double> > sortPositions(2,vector <double> (List));
	for(int k=0;k<List;k++)
	{
		sortPositions[0][k]=output[k].LogProbability;
		sortPositions[1][k]=k;
	}
	sortPositions=SortedArray(sortPositions,Log2L);
	for(int j=0;j<List;j++)
	{
		out1[j]=output[sortPositions[1][j]].InfoBits;
	}

	return out1;
}

vector <vector <double> > RightCollapse(vector <vector<double> > input,int Log2M,int row,int StartPoint)
{
	//vector <vector <double> > output=input;output==>   input????FAST???
	int initSP=StartPoint+pow(2,Log2M-row-1);
	if(row==Log2M)
	{
		return input;//I don't Know if it works or not?
	}
	for(int i=row+1;i<Log2M+1;i++)
	{
		for(int j=initSP;j<initSP+pow(2,Log2M-i);j++)
		{
			input[i][j]=input[i-1][j]*input[i-1][j-pow(2,Log2M-i)];
		}
		initSP=initSP+pow(2,Log2M-i-1);
	}
	return input;
}
vector <int> PrefixDetector(int NInfoBitDecodedSofar,int Log2M)
{
	vector <int> BRFirst=BinaryRep(NInfoBitDecodedSofar-1,Log2M);
	vector <int> BRSecond=BinaryRep(NInfoBitDecodedSofar,Log2M);
	int CommonPrefix=0;
	for(int i=0;i<Log2M;i++)
	{
		if(BRFirst[i]==BRSecond[i])
		{
			CommonPrefix++;
		}
		else
		{
			break;
		}
	}
	vector <int> PrefixBinary(CommonPrefix,0);
	for(int i=0;i<CommonPrefix;i++)
	{
		PrefixBinary[i]=BRFirst[i];
	}
	int Prefix=BinaryToint(PrefixBinary);
	vector <int> out(2,0);
	out[0]=CommonPrefix;
	out[1]=Prefix;
	return out;
}
vector <int> RMPlotkinEncoder(vector <int> Information,int M,int r)
{
	if(r==0)
	{
		int Length=pow(2,M);
		vector <int> out(Length);
		for(int i=0;i<Length;i++)
		{
			out[i]=Information[0];
		}
		return out;
	}
	if((M==r)&&(r==1))
	{
		vector <int> out(2);
		out[0]=Information[0];
		out[1]=Information[1]*Information[0];
		return out;
	}
	
	if(M==r)//this should be changed if you want to use real plotkin
	{
		int Length=pow(2,M);
		vector <int> out(Length);
		vector <int> Leftin(Length/2);
		vector <int> Rightin(Length/2);
		for(int i=0;i<Length/2;i++)
		{
			Leftin[i]=Information[i];
			Rightin[i]=Information[i+Length/2];
		}
		vector <int> Leftout=RMPlotkinEncoder(Leftin,M-1,r-1);
		vector <int> Rightout=RMPlotkinEncoder(Rightin,M-1,r-1);
		for(int i=0;i<Length/2;i++)
		{
			out[i]=Leftout[i];
			out[i+Length/2]=Leftout[i]*Rightout[i];
		}
		
		return out;
	}

	int k=RMInfoBitCounter(M,r);
	int k1=RMInfoBitCounter(M-1,r);
	int k2=RMInfoBitCounter(M-1,r-1);
	vector <int> InfL(k1);
	vector <int> InfR(k2);
	for(int i=0;i<k1;i++)
	{
		InfL[i]=Information[i];
	}
	for(int i=0;i<k2;i++)
	{
		InfR[i]=Information[i+k1];
	}
	int Length=pow(2,M);
	vector <int> CWL=RMPlotkinEncoder(InfL,M-1,r);
	vector <int> CWR=RMPlotkinEncoder(InfR,M-1,r-1);
	vector <int> out(Length);
	for(int i=0;i<Length/2;i++)
	{
		out[i]=CWL[i];
		out[i+Length/2]=CWL[i]*CWR[i];
	}
	return out;
}
vector <vector <double> > SortedArray(vector <vector <double> >Input,int log2M )
{
	int M=pow(2,log2M);
	if(log2M==0)
	{
		return Input;
	}
	else
	{
		vector <vector <double> >LeftHalf(2);
		LeftHalf[0].assign(Input[0].begin(),Input[0].begin()+M/2);
		LeftHalf[1].assign(Input[1].begin(),Input[1].begin()+M/2);
		vector <vector <double> >RightHalf(2);
		RightHalf[0].assign(Input[0].begin()+M/2,Input[0].end());
		RightHalf[1].assign(Input[1].begin()+M/2,Input[1].end());
		LeftHalf=SortedArray(LeftHalf,log2M-1);
		RightHalf=SortedArray(RightHalf,log2M-1);
		vector <vector <double> > out (2,vector <double> (M,0));
		int index1=0;
		int index2=0;
		for(int i=0;i<M;i++)
		{
			if(index1==M/2)
			{
				for(int j=i;j<M;j++)
				{
					out[0][j]=LeftHalf[0][index2];
					out[1][j]=LeftHalf[1][index2];
					index2++;
				}
				return out;
			}
			if(index2==M/2)
			{
				for(int j=i;j<M;j++)
				{
					out[0][j]=RightHalf[0][index1];
					out[1][j]=RightHalf[1][index1];
					index1++;
				}
				return out;
			}
			if(RightHalf[0][index1]>LeftHalf[0][index2])
			{
				out[0][i]=RightHalf[0][index1];
				out[1][i]=RightHalf[1][index1];
				index1++;
			}
			else
			{
				out[0][i]=LeftHalf[0][index2];
				out[1][i]=LeftHalf[1][index2];
				index2++;
			}
		}
		return out;
	}
}
vector <int> BinaryRep(int N,int Digits)
{
	vector <int> out(Digits);
	int temp=N;
	for(int i=0;i<Digits;i++)
	{
		out[Digits-1-i]=temp%2;
		temp=temp/2;
	}
	return out;
}
int BinaryToint(vector <int> v)
{
	int size1=v.size();
	int out=0;
	for(int i=0;i<size1;i++)
	{
		out=out*2+v[i];
	}
	return out;
}
int RMInfoBitCounter(int M,int R)
{
	if(M==R)
	{
		return pow(2,M);
	}
	if(R==0)
	{
		return 1;
	}
	return RMInfoBitCounter(M-1,R)+RMInfoBitCounter(M-1,R-1);
}
vector <int> RMTableChecker_fullcode(vector <vector <int> > ListofEstimates,int Log2M, int r,int Log2L,vector <double> EpsRec
,vector <double> PCEpsRec,vector < vector <int> > ParityCheckPositions)
{
	int List=pow(2,Log2L);
	int M=EpsRec.size();
	int n=PCEpsRec.size();
	vector <double> ProbMeas (List,0);
	for(int i=0;i<List;i++)
	{
		vector <int> TempCW=RMPlotkinEncoder(ListofEstimates[i],Log2M,r);
		for(int j=0;j<M;j++)
		{
			if((1+TempCW[j]*EpsRec[j])!=0)
			{
				ProbMeas[i]=ProbMeas[i]+log((1+TempCW[j]*EpsRec[j])/2);
			}
			else if((1+TempCW[j]*EpsRec[j])==0)
			{
				ProbMeas[i]=ProbMeas[i]-20;
			}
		}
		vector <int> temp2=Encoder1D(TempCW,ParityCheckPositions);
		for(int j=0;j<n;j++)
		{
			ProbMeas[i]+=log((1+temp2[j]*PCEpsRec[j])/(1-temp2[j]*PCEpsRec[j]));
		}
	}
	int Position=0;
	double MeasureMax=ProbMeas[0];
	for(int i=0;i<List;i++)
	{
		if(MeasureMax<ProbMeas[i])
		{
			MeasureMax=ProbMeas[i];
			Position=i;
		}
	}
	//cout<<"\nPosition of the Maximum Estimate is : "<<Position;
	vector <int> out=RMPlotkinEncoder(ListofEstimates[Position],Log2M,r);
	return out;
}

vector <int> RMTableChecker(vector <vector <int> > ListofEstimates,int Log2M, int r,int Log2L,vector <double> EpsRec)
{
	int List=pow(2,Log2L);
	int M=pow(2,Log2M);
	vector <double> ProbMeas (List,0);
	for(int i=0;i<List;i++)
	{
		vector <int> TempCW=RMPlotkinEncoder(ListofEstimates[i],Log2M,r);
		for(int j=0;j<M;j++)
		{
			if((1+TempCW[j]*EpsRec[j])!=0)
			{
				ProbMeas[i]=ProbMeas[i]+log((1+TempCW[j]*EpsRec[j])/2);
			}
			else if((1+TempCW[j]*EpsRec[j])==0)
			{
				ProbMeas[i]=ProbMeas[i]-20;
			}
		}
	}
	int Position=0;
	double MeasureMax=ProbMeas[0];
	for(int i=0;i<List;i++)
	{
		if(MeasureMax<ProbMeas[i])
		{
			MeasureMax=ProbMeas[i];
			Position=i;
		}
	}
	//cout<<"\nPosition of the Maximum Estimate is : "<<Position;
	vector <int> out=ListofEstimates[Position];
	return out;
}
vector <int> FrozenBitsPolar(int k1)
{
	
	vector <int> FrozenBits
	{
	//2047 , 2046 , 2045 , 2043 , 2039 , 2031 , 2015 , 1983 , 1919 , 1791 , 1535 , 1023 , 2044 , 2042 , 2041 , 2038 , 2037 , 2035 , 2030 , 2029 , 2027 , 2023 , 2014 , 2013 , 2011 , 2007 , 1999 , 1982 , 1981 , 1979 , 1975 , 1967 , 1951 , 1918 , 1917 , 1915 , 1911 , 1903 , 1887 , 1855 , 1790 , 1789 , 1787 , 1783 , 1775 , 1759 , 1727 , 1534 , 1533 , 1531 , 1527 , 1519 , 1503 , 1663 , 1471 , 1022 , 1021 , 1019 , 1015 , 1007 , 1407 , 2040 , 2036 , 2034 , 2033 , 2028 , 2026 , 2025 , 2022 , 2021 , 2019 , 2012 , 2010 , 2009 , 2006 , 2005 , 2003 , 1998 , 1997 , 1995 , 991 , 1991 , 1980 , 1978 , 1977 , 1974 , 1973 , 1971 , 1966 , 1965 , 1963 , 1959 , 1950 , 1949 , 1947 , 1943 , 1935 , 1916 , 1914 , 1913 , 1910 , 1909 , 1907 , 1902 , 1901 , 1899 , 1895 , 1886 , 1885 , 1883 , 1879 , 1871 , 1854 , 1853 , 1851 , 1847 , 1839 , 1788 , 1786 , 1785 , 1782 , 1781 , 1779 , 1774 , 1773 , 1771 , 1767 , 959 , 1758 , 1757 , 1755 , 1751 , 1743 , 1823 , 1726 , 1725 , 1723 , 1719 , 1279 , 1711 , 1532 , 1530 , 1529 , 1526 , 1525 , 1523 , 1518 , 1517 , 1515 , 1511 , 1695 , 1502 , 1501 , 1499 , 1495 , 895 , 1662 , 1661 , 1659 , 1655 , 1487 , 1647 , 1470 , 1469 , 1467 , 1463 , 1455 , 1631 , 1020 , 1018 , 1017 , 1014 , 1013 , 1011 , 1006 , 1005 , 1003 , 999 , 1406 , 1405 , 1403 , 1439 , 767 , 1399 , 2032 , 2024 , 2020 , 2018 , 2017 , 2008 , 2004 , 2002 , 2001 , 1996 , 1994 , 1993 , 990 , 1990 , 989 , 1989 , 987 , 1987 , 1976 , 1972 , 1970 , 1969 , 1964 , 1962 , 1961 , 1958 , 1957 , 1955 , 1948 , 1946 , 1945 , 1942 , 1941 , 1939 , 1934 , 1933 , 983 , 1931 , 1912 , 1908 , 1906 , 1905 , 1900 , 1898 , 1897 , 1894 , 1893 , 1891 , 1884 , 1882 , 1881 , 1878 , 1877 , 1927 , 1875 , 1870 , 1869 , 1867 , 1863 , 1391 , 1852 , 1850 , 1849 , 1846 , 1845 , 1843 , 1599 , 1838 , 1837 , 1835 , 975 , 1831 , 1784 , 1780 , 1778 , 1777 , 1772 , 1770 , 1769 , 1766 , 1765 , 1763 , 958 , 957 , 955 , 1756 , 1754 , 1753 , 1750 , 1749 , 1747 , 1742 , 1741 , 1822 , 1821 , 1739 , 1819 , 951 , 1735 , 1815 , 1724 , 1722 , 1721 , 1718 , 1717 , 1715 , 1278 , 1277 , 1275 , 1710 , 1709 , 1707 , 1375 , 1271 , 1703 , 943 , 511 , 1528 , 1524 , 1522 , 1521 , 1516 , 1514 , 1513 , 1510 , 1807 , 1509 , 1507 , 1694 , 1693 , 1500 , 1498 , 1497 , 1494 , 1691 , 1493 , 1491 , 894 , 893 , 891 , 1660 , 1658 , 1657 , 1654 , 1653 , 1651 , 1486 , 1485 , 1483 , 1687 , 1263 , 1646 , 1645 , 1643 , 887 , 1468 , 1466 , 1465 , 1462 , 1461 , 1479 , 1459 , 1343 , 927 , 1639 , 1454 , 1453 , 1451 , 1679 , 1630 , 1629 , 1627 , 879 , 1447 , 1016 , 1012 , 1010 , 1247 , 1009 , 1004 , 1002 , 1001 , 998 , 997 , 1404 , 1402 , 1438 , 1401 , 1437 , 766 , 995 , 765 , 1398 , 1397 , 2016 , 2000 , 1992 , 988 , 1988 , 986 , 1986 , 1968 , 1960 , 1956 , 1954 , 985 , 1985 , 1953 , 1944 , 1940 , 1938 , 1435 , 1937 , 1623 , 763 , 1932 , 982 , 1930 , 981 , 1929 , 1395 , 1904 , 1896 , 1892 , 1890 , 1889 , 1880 , 1876 , 1926 , 1874 , 1925 , 1873 , 1868 , 1866 , 979 , 1865 , 1862 , 1923 , 1861 , 1390 , 1848 , 1844 , 1389 , 1842 , 1841 , 1598 , 1597 , 1836 , 1834 , 1833 , 974 , 1859 , 973 , 1387 , 1830 , 1431 , 1829 , 759 , 1595 , 1776 , 1768 , 863 , 1764 , 1762 , 956 , 971 , 954 , 1761 , 1752 , 1748 , 953 , 1746 , 1745 , 1827 , 1740 , 1820 , 1738 , 1818 , 950 , 1737 , 1817 , 1615 , 949 , 1734 , 1814 , 1733 , 1813 , 1720 , 1716 , 1215 , 1714 , 1713 , 1276 , 1274 , 947 , 1708 , 1273 , 1706 , 1374 , 1705 , 1383 , 1373 , 1731 , 1811 , 1591 , 1270 , 1269 , 1702 , 967 , 1701 , 942 , 1423 , 941 , 751 , 1371 , 510 , 1520 , 1512 , 509 , 1806 , 1508 , 1506 , 1805 , 1692 , 1267 , 1496 , 1505 , 1690 , 1492 , 1490 , 1699 , 892 , 1689 , 890 , 1489 , 1656 , 939 , 1652 , 1650 , 889 , 1484 , 1482 , 507 , 1649 , 1686 , 1262 , 1481 , 1803 , 1644 , 1685 , 831 , 1642 , 1261 , 886 , 1464 , 1641 , 885 , 1460 , 1583 , 1478 , 1458 , 1367 , 1342 , 1477 , 926 , 1457 , 1638 , 1341 , 1452 , 925 , 1683 , 1450 , 1259 , 1637 , 935 , 883 , 1449 , 1151 , 1678 , 1628 , 503 , 1475 , 1799 , 1626 , 735 , 878 , 1677 , 1446 , 1339 , 1246 , 1008 , 923 , 1000 , 1625 , 877 , 996 , 1635 , 1445 , 1245 , 1400 , 1436 , 994 , 764 , 1396 , 984 , 1984 , 1952 , 1434 , 1936 , 1622 , 762 , 980 , 1928 , 1255 , 1394 , 993 , 1359 , 1888 , 1924 , 1872 , 1675 , 978 , 1864 , 1433 , 1621 , 761 , 1567 , 1922 , 1860 , 1393 , 875 , 1388 , 1840 , 1596 , 977 , 1443 , 1832 , 1858 , 1243 , 972 , 1386 , 1921 , 1430 , 495 , 1828 , 758 , 1335 , 1594 , 919 , 862 , 970 , 1760 , 952 , 1744 , 1857 , 1826 , 1385 , 1736 , 1429 , 1816 , 1614 , 1619 , 948 , 757 , 1593 , 861 , 969 , 1732 , 1812 , 1671 , 1214 , 703 , 1712 , 1825 , 946 , 1272 , 1613 , 1704 , 1382 , 1372 , 871 , 1730 , 1810 , 1590 , 1268 , 1239 , 966 , 1213 , 1700 , 1422 , 940 , 750 , 945 , 1370 , 1427 , 508 , 755 , 1381 , 1804 , 1729 , 1266 , 1809 , 1589 , 1504 , 859 , 1698 , 1327 , 1688 , 965 , 1488 , 911 , 938 , 888 , 1421 , 749 , 1611 , 506 , 1369 , 1648 , 1480 , 1802 , 479 , 1684 , 830 , 1260 , 1265 , 1211 , 1640 , 884 , 1697 , 1582 , 1366 , 937 , 1476 , 1379 , 1456 , 1587 , 1340 , 505 , 924 , 1682 , 1258 , 1636 , 963 , 1801 , 1231 , 934 , 829 , 639 , 882 , 855 , 1448 , 1150 , 1419 , 502 , 747 , 1581 , 1474 , 1798 , 734 , 1365 , 1676 , 1607 , 1338 , 922 , 1624 , 876 , 1634 , 1257 , 1311 , 1207 , 1681 , 1444 , 1244 , 933 , 881 , 1149 , 501 , 1254 , 827 , 992 , 1358 , 1473 , 1797 , 733 , 1674 , 447 , 1432 , 1620 , 1337 , 760 , 1579 , 921 , 1566 , 1392 , 1363 , 874 , 1633 , 1415 , 976 , 1442 , 743 , 1242 , 1920 , 847 , 494 , 1334 , 918 , 1253 , 931 , 1357 , 1856 , 1147 , 1673 /*, 1384 , 499 , 1428 , 1618 , 756 , 1565 , 1592 , 1795 , 731 , 873 , 860 , 968 , 1670 ,702 , 1441 , 1199 , 1241 , 1824 , 823 , 493 , 1612 , 1333 , 870 , 917 , 1575 , 1238 , 1212 , 1251 , 1355 , 944 , 1617 , 1426 , 383 , 754 , 1380 , 1728 , 1808 , 1588 , 858 , 1563 , 1669 , 701 , 1143 , 1326 , 964 , 910 , 1420 , 869 , 748 , 1610 , 727 , 1368 , 491 , 1331 , 1237 , 478 , 915 , 1264 , 1210 , 815 , 1696 , 1425 , 1183 , 753 , 936 , 1378 , 857 , 1586 , 1325 , 504 , 909 , 1351 , 962 , 1800 , 1230 , 1667 , 699 , 1609 , 828 , 638 , 854 , 1418 , 477 , 746 , 867 , 1559 , 1580 , 1364 , 1209 , 1135 , 1606 , 1235 , 719 , 255 , 487 , 1377 , 1680 , 1585 , 1256 , 1310 /*, 1206 , 932 , 961 , 880 , 1229 , 1148 , 1323 , 799 , 637 , 853 , 907 , 500 , 826 , 1417 , 745 , 1472 , 1796 , 732 , 446 , 695 , 1336 , 475 , 1578 , 1605 , 920 , 1362 , 1632 , 1414 , 1551 , 742 , 1309 , 1205 , 846 , 1119 , 1252 , 930 , 1356 , 1146 , 1672 , 825 , 1227 , 498 , 1319 , 635 , 445 , 851 , 903 , 1564 , 1794 , 730 , 872 , 1577 , 1361 , 1440 , 1198 , 1603 , 1240 , 471 , 687 , 1413 , 822 , 741 , 492 , 1332 , 916 , 1574 , 845 , 1307 , 1203 , 929 , 1145 , 1250 , 1354 , 497 , 1616 , 382 , 1087 , 1793 , 1223 , 729 , 443 , 631 , 1562 , 1668 , 700 , 1197 , 1142 , 821 , 1411 , 739 , 868 , 726 , 1573 , 463 , 843 , 490 , 671 , 1330 , 1236 , 1303 , 914 , 1249 , 1353 , 814 , 1424 , 1182 , 381 , 752 , 856 , 439 , 1561 , 1141 , 1324 , 1195 , 623 , 908 , 819 , 1350 , 1666 , 698 , 725 , 1608 , 1571 , 489 , 1329 , 839 , 913 , 476 , 866 , 1558 , 1208 , 1134 , 813 , 1181 , 1295 , 1234 , 379 , 718 , 254 , 486 , 1376 , 1139 , 1584 , 1349 , 431 , 1191 , 1665 , 697 , 960 , 723 , 1228 , 607 , 1322 , 798 , 636 , 852 , 906 , 865 , 1557 , 1416 , 1133 , 744 , 811 , 1233 , 694 , 1179 , 474 , 1604 , 717 , 253 , 485 , 375 , 1550 , 1347 , 1308 , 1204 , 1321 , 1118 , 797 , 415 , 905 , 1555 , 575 , 1131 , 693 , 824 , 1226 , 473 , 807 , 1318 , 634 , 444 , 850 , 1175 , 715 , 251 , 902 , 483 , 1576 , 1549 , 367 , 1360 , 1602 , 470 , 686 , 1412 , 1117 , 740 , 795 , 844 , 1306 , 1202 , 1225 , 1127 , 691 , 928 , 1317 , 633 , 849 , 901 , 1144 , 711 , 247 , 1167 , 496 , 1547 , 1086 , 1792 , 1222 , 1601 , 728 , 442 , 351 , 469 , 685 , 630 , 1115 , 1196 , 791 , 820 , 1305 , 1201 , 1410 , 738 , 1315 , 1572 , 899 , 462 , 842 , 670 , 1302 , 239 , 1085 , 1543 , 1221 , 441 , 1248 , 467 , 629 , 683 , 319 , 1352 , 1111 , 380 , 783 , 1409 , 737 , 438 , 1560 , 1140 , 1194 , 461 , 841 , 622 , 818 , 669 , 1301 , 223 , 724 , 1083 , 1219 , 1570 , 627 , 488 , 679 , 1328 , 1103 , 838 , 912 , 437 , 812 , 1180 , 1294 , 1193 , 459 , 621 , 378 , 817 , 667 , 191 , 1299 , 1079 , 1138 , 1569 , 1348 , 430 , 1190 , 837 , 1664 , 696 , 722 , 606 , 435 , 1293 , 455 , 864 , 1556 , 377 , 619 , 1132 , 663 , 127 , 1071 , 810 , 1232 , 1137 , 1178 , 716 , 429 , 252 , 1189 , 835 , 484 , 374 , 721 , 605 , 1291 , 1346 , 615 , 655 , 1055 , 1320 , 796 , 809 , 414 , 904 , 1177 , 427 , 1187 , 1554 , 373 , 574 , 1130 , 603 , 692 , 1287 , 1345 , 472 , 806 , 1174 , 714 , 250 , 482 , 413 , 1548 , 366 , 423 , 1553 , 371 , 573 , 1129 , 599 , 1116 , 794 , 805 , 1173 , 713 , 249 , 481 , 411 , 365 , 1224 , 1126 , 690 , 1316 , 632 , 848 , 571 , 900 , 591 , 710 , 246 , 1166 , 793 , 803 , 1171 , 1546 , 407 , 1600 , 350 , 468 , 684 , 363 , 1125 , 689 , 1114 , 567 , 790 , 709 , 245 , 1165 , 1304 , 1200 , 1545 , 399 , 349 , 1314 , 359 , 1123 , 898 , 1113 , 559 , 789 , 238 , 707 , 243 , 1163 , 1084 , 1542 , 1220 , 440 , 466 , 628 , 682 , 318 , 347 , 1313 , 1110 , 897 , 543 , 782 , 787 , 237 , 1159 , 1408 , 736 , 1541 , 465 , 681 , 317 , 343 , 460 , 840 , 1109 , 668 , 1300 , 781 , 222 , 235 , 1082 , 1218 , 1539 , 626 , 678 , 315 , 335 , 1102 , 1107 , 779 , 221 , 231 , 1081 , 436 , 1217 , 625 , 1192 , 311 , 677 , 458 , 620 , 816 , 1101 , 666 , 190 , 1298 , 775 , 219 , 1078 , 1568 , 303 , 675 , 457 , 836 , 1099 , 665 , 189 , 1297 , 215 , 1077 , 434 , 1292 , 287 , 454 , 376 , 618 , 1095 , 662 , 126 , 187 , 207 , 1070 , 1075 , 433 , 1136 , 453 , 428 , 1188 , 834 , 617 , 661 , 125 , 183 , 720 , 1069 , 604 , 1290 , 451 , 614 , 833 , 654 , 659 , 123 , 175 , 1054 , 1067 , 808 , 1289 , 1176 , 426 , 1186 , 613 , 653 , 119 , 159 , 372 , 1053 , 1063 , 602 , 1286 , 1344 , 611 , 425 , 1185 , 111 , 651 , 1051 , 601 , 412 , 1285 , 422 , 1552 , 95 , 647 , 370 , 1047 , 572 , 1128 , 598 , 1283 , 421 , 63 , 804 , 1172 , 712 , 248 , 369 , 1039 , 480 , 597 , 410 , 364 , 419 , 570 , 590 , 595 , 409 , 792 , 802 , 1170 , 569 , 589 , 406 , 362 , 1124 , 801 , 688 , 1169 , 566 , 587 , 405 , 708 , 361 , 244 , 1164 , 1544 , 565 , 583 , 398 , 403 , 348 , 358 , 1122 , 1112 , 558 , 563 , 397 , 788 , 357 , 706 , 242 , 1162 , 1121 , 557 , 395 , 346 , 355 , 705 , 1312 , 241 , 1161 , 896 , 542 , 555 , 391 , 786 , 345 , 236 , 1158 , 1540 , 541 , 551 , 464 , 680 , 316 , 785 , 342 , 1157 , 1108 , 539 , 780 , 341 , 234 , 1155 , 1538 , 535 , 314 , 334 , 339 , 233 , 1106 , 1537 , 527 , 778 , 313 , 333 , 220 , 230 , 1080 , 1105 , 1216 , 624 , 310 , 676 , 777 , 331 , 229 , 1100 , 774 , 309 , 327 , 218 , 227 , 302 , 773 , 307 , 674 , 456 , 217 , 1098 , 664 , 301 , 771 , 188 , 1296 , 673 , 214 , 1076 , 1097 , 286 , 299 , 213 , 1094 , 285 , 295 , 186 , 206 , 211 , 1093 , 1074 , 432 , 283 , 185 , 452 , 205 , 616 , 1091 , 1073 , 660 , 124 , 279 , 182 , 203 , 1068 , 271 , 181 , 199 , 450 , 832 , 658 , 122 , 174 , 179 , 449 , 1066 , 1288 , 657 , 121 , 173 , 1065 , 612 , 652 , 118 , 158 , 171 , 1052 , 1062 , 117 , 157 , 167 , 1061 , 610 , 424 , 1184 , 110 , 650 , 115 , 155 , 1050 , 1059 , 609 , 600 , 109 , 649 , 1284 , 151 , 1049 , 94 , 646 , 107 , 143 , 1046 , 93 , 103 , 645 , 1282 , 1045 , 420 , 62 , 91 , 643 , 1281 , 368 , 1038 , 1043 , 61 , 596 , 87 , 1037 , 418 , 59 , 79 , 1035 , 417 , 55 , 594 , 408 , 1031 , 47 , 593 , 568 , 31 , 588 , 800 , 1168 , 586 , 404 , 360 , 585 , 564 , 582 , 402 , 581 , 401 , 562 , 579 , 396 , 561 , 356 , 1120 , 556 , 394 , 354 , 704 , 240 , 1160 , 393 , 554 , 353 , 390 , 344 , 553 , 389 , 540 , 550 , 387 , 784 , 549 , 1156 , 538 , 547 , 537 , 340 , 1154 , 534 , 1153 , 533 , 338 , 232 , 1536 , 526 , 531 , 337 , 312 , 525 , 332 , 1104 , 523 , 776 , 519 , 330 , 228 , 329 , 308 , 326 , 226 , 325 , 225 , 772 , 306 , 323 , 216 , 305 , 300 , 770 , 672 , 1096 , 769 , 298 , 212 , 297 , 284 , 294 , 210 , 1092 , 293 , 209 , 282 , 291 , 184 , 204 , 1090 , 1072 , 281 , 1089 , 278 , 202 , 277 , 201 , 270 , 275 , 180 , 198 , 269 , 197 , 267 , 178 , 195 , 448 , 263 , 177 , 656 , 120 , 172 , 1064 , 170 , 169 , 116 , 156 , 166 , 1060 , 165 , 114 , 154 , 163 , 113 , 1058 , 608 , 153 , 108 , 648 , 1057 , 150 , 1048 , 149 , 106 , 142 , 147 , 105 , 141 , 92 , 102 , 644 , 139 , 1044 , 101 , 135 , 90 , 99 , 642 , 1280 , 89 , 1042 , 641 , 60 , 86 , 1041 , 85 , 1036 , 58 , 78 , 83 , 57 , 77 , 1034 , 416 , 54 , 75 , 1033 , 53 , 71 , 1030 , 46 , 51 , 1029 , 592 , 45 , 1027 , 30 , 43 , 29 , 39 , 27 , 23 , 15 , 584 , 580 , 400 , 578 , 577 , 560 , 392 , 352 , 552 , 388 , 386 , 385 , 548 , 546 , 545 , 536 , 1152 , 532 , 530 , 336 , 529 , 524 , 522 , 521 , 518 , 517 , 515 , 328 , 324 , 224 , 322 , 321 , 304 , 768 , 296 , 292 , 208 , 290 , 289 , 280 , 1088 , 276 , 200 , 274 , 273 , 268 , 196 , 266 , 194 , 265 , 193 , 262 , 176 , 261 , 259 , 168 , 164 , 162 , 161 , 112 , 152 , 1056 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 , 133 , 98 , 131 , 97 , 88 , 640 , 1040 , 84 , 82 , 81 , 56 , 76 , 74 , 1032 , 73 , 52 , 70 , 69 , 50 , 67 , 1028 , 49 , 44 , 1026 , 1025 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 576 , 384 , 544 , 528 , 520 , 516 , 514 , 513 , 320 , 288 , 272 , 264 , 192 , 260 , 258 , 257 , 160 , 144 , 136 , 132 , 130 , 129 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 1024 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 512 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@-2dB   
	//2047 , 2046 , 2045 , 2043 , 2039 , 2031 , 2015 , 1983 , 1919 , 1791 , 1535 , 2044 , 2042 , 2041 , 2038 , 2037 , 2035 , 2030 , 2029 , 2027 , 2023 , 2014 , 2013 , 2011 , 2007 , 1999 , 1982 , 1981 , 1979 , 1975 , 1967 , 1951 , 1918 , 1917 , 1915 , 1911 , 1903 , 1887 , 1023 , 1790 , 1789 , 1787 , 1783 , 1775 , 1855 , 1759 , 1534 , 1533 , 1531 , 1527 , 1727 , 1519 , 2040 , 2036 , 2034 , 2033 , 2028 , 2026 , 2025 , 2022 , 2021 , 2019 , 2012 , 2010 , 2009 , 2006 , 2005 , 2003 , 1998 , 1997 , 1995 , 1991 , 1980 , 1978 , 1977 , 1974 , 1973 , 1971 , 1966 , 1965 , 1963 , 1959 , 1503 , 1950 , 1949 , 1947 , 1943 , 1916 , 1914 , 1913 , 1910 , 1909 , 1907 , 1663 , 1902 , 1901 , 1899 , 1935 , 1895 , 1886 , 1885 , 1883 , 1879 , 1022 , 1021 , 1019 , 1015 , 1471 , 1871 , 1788 , 1786 , 1785 , 1782 , 1781 , 1779 , 1007 , 1774 , 1773 , 1771 , 1854 , 1853 , 1851 , 1767 , 1847 , 1758 , 1757 , 1755 , 1751 , 1839 , 991 , 1407 , 1532 , 1530 , 1529 , 1526 , 1525 , 1743 , 1726 , 1725 , 1523 , 1723 , 1518 , 1517 , 1515 , 1719 , 1511 , 1823 , 2032 , 2024 , 2020 , 2018 , 2017 , 2008 , 2004 , 2002 , 2001 , 1996 , 1994 , 1993 , 1990 , 1989 , 1976 , 1972 , 1970 , 1969 , 1987 , 1964 , 1962 , 1961 , 1958 , 1957 , 959 , 1502 , 1501 , 1955 , 1499 , 1948 , 1946 , 1945 , 1942 , 1941 , 1939 , 1711 , 1912 , 1908 , 1906 , 1905 , 1662 , 1661 , 1900 , 1898 , 1897 , 1934 , 1933 , 1894 , 1495 , 1893 , 1659 , 1931 , 1891 , 1279 , 1884 , 1882 , 1881 , 1878 , 1877 , 1020 , 1018 , 1017 , 1014 , 1013 , 1875 , 1470 , 1469 , 1655 , 1011 , 1927 , 1467 , 1870 , 1869 , 1784 , 1780 , 1778 , 1777 , 1006 , 1005 , 1772 , 1487 , 1770 , 1867 , 1769 , 1695 , 1852 , 1850 , 1849 , 1766 , 1765 , 1003 , 1846 , 1845 , 895 , 1763 , 1756 , 1754 , 1463 , 1753 , 1843 , 1647 , 1750 , 1749 , 1863 , 1838 , 1837 , 990 , 989 , 1747 , 999 , 1406 , 1405 , 1835 , 1528 , 987 , 1524 , 1742 , 1724 , 1522 , 1722 , 1741 , 1521 , 1721 , 1455 , 1516 , 1403 , 1514 , 1513 , 1718 , 1717 , 1739 , 1510 , 1822 , 1631 , 1509 , 2016 , 2000 , 1992 , 1821 , 1988 , 1831 , 1968 , 1986 , 1960 , 983 , 1956 , 958 , 1500 , 1954 , 1985 , 1715 , 767 , 1498 , 957 , 1944 , 1940 , 1953 , 1399 , 1497 , 1938 , 1507 , 1710 , 1819 , 1937 , 1904 , 1709 , 1660 , 1896 , 1735 , 1932 , 1494 , 1892 , 1658 , 1930 , 1890 , 955 , 1493 , 1278 , 1880 , 1657 , 1439 , 1929 , 1876 , 1889 , 1016 , 1277 , 1012 , 1874 , 1468 , 1707 , 975 , 1654 , 1010 , 1926 , 1873 , 1466 , 1815 , 1868 , 1653 , 1009 , 1391 , 1776 , 1925 , 1491 , 1004 , 1486 , 1866 , 1599 , 1768 , 1694 , 1465 , 1848 , 1764 , 1002 , 1275 , 1844 , 951 , 1485 , 1865 , 1693 , 894 , 1762 , 1462 , 1752 , 1842 , 1001 , 1646 , 1703 , 1651 , 1748 , 511 , 1923 , 893 , 1862 , 1761 , 1461 , 1841 , 1836 , 988 , 1746 , 998 , 1645 , 1483 , 1807 , 1691 , 1404 , 1861 , 1834 , 1271 , 986 , 1740 , 1745 , 997 , 1520 , 1720 , 1454 , 1402 , 891 , 943 , 1512 , 1375 , 1459 , 1833 , 1716 , 985 , 1738 , 1643 , 1630 , 1508 , 1859 , 1453 , 1820 , 1830 , 1401 , 982 , 1479 , 995 , 1687 , 1984 , 1714 , 766 , 956 , 1737 , 1952 , 1398 , 1496 , 1506 , 1629 , 1818 , 1263 , 1936 , 887 , 1708 , 1829 , 1734 , 981 , 954 , 1713 , 1639 , 1492 , 765 , 1451 , 1656 , 1438 , 1397 , 1928 , 927 , 1888 , 1505 , 1276 , 1343 , 1817 , 1706 , 974 , 1733 , 1627 , 1872 , 1814 , 1827 , 1679 , 953 , 1652 , 1008 , 1390 , 979 , 1924 , 1490 , 1437 , 1598 , 1464 , 763 , 1274 , 1395 , 879 , 950 , 1484 , 1705 , 973 , 1864 , 1692 , 1247 , 1447 , 1813 , 1000 , 1731 , 1702 , 1389 , 1650 , 1489 , 510 , 1597 , 1922 , 892 , 1760 , 1460 , 1840 , 1435 , 1623 , 1273 , 949 , 1644 , 1482 , 1806 , 1690 , 1860 , 971 , 759 , 1270 , 1701 , 1649 , 1811 , 1744 , 996 , 509 , 1921 , 1387 , 890 , 863 , 942 , 1595 , 1374 , 1458 , 1832 , 1215 , 984 , 1642 , 1481 , 1805 , 1689 , 947 , 1269 , 1858 , 1452 , 1431 , 1400 , 1615 , 1478 , 994 , 1686 , 889 , 1699 , 941 , 1736 , 1373 , 1457 , 507 , 967 , 751 , 1628 , 1641 , 1262 , 886 , 1828 , 1383 , 1857 , 1591 , 980 , 1803 , 1712 , 1638 , 1267 , 764 , 1450 , 831 , 1477 , 993 , 1685 , 1396 , 926 , 1504 , 1423 , 1151 , 1342 , 1816 , 939 , 1371 , 1261 , 1732 , 885 , 1626 , 503 , 1826 , 1678 , 952 , 735 , 978 , 1637 , 1449 , 1436 , 762 , 1799 , 1583 , 1475 , 925 , 1683 , 1394 , 878 , 1704 , 972 , 1341 , 1246 , 1446 , 1812 , 1625 , 935 , 1730 , 1259 , 1367 , 883 , 1388 , 1825 , 1677 , 1488 , 1596 , 977 , 495 , 1434 , 1635 , 761 , 1622 , 1272 , 948 , 1393 , 877 , 923 , 703 , 1245 , 1445 , 1339 , 970 , 1567 , 758 , 1729 , 1700 , 1648 , 1810 , 508 , 1675 , 1920 , 1255 , 1386 , 1359 , 862 , 1594 , 1621 , 1433 , 479 , 1214 , 875 , 1480 , 1804 , 1688 , 946 , 919 , 1243 , 1443 , 969 , 1268 , 757 , 1335 , 1430 , 1614 , 639 , 1809 , 1385 , 888/* , 1671 , 1698 , 861 , 1593 , 940 , 1372 , 1456 , 506 , 966 , 1619 , 750 , 1213 , 1640 , 945 , 871 , 1382 , 1429 , 1239 , 755 , 1856 , 1590 , 447 , 1613 , 1802 , 911 , 1327 , 1266 , 830 , 1476 , 992 , 1684 , 1697 , 859 , 505 , 965 , 749 , 1422 , 1150 , 1211 , 938 , 1370 , 1260 , 1381 , 884 , 1589 , 502 , 1427 , 1801 , 1611 , 1231 , 1265 , 829 , 734 , 383 , 1636 , 1448 , 1311 , 963 , 1798 , 1421 , 747 , 855 , 1582 , 1474 , 1149 , 924 , 1682 , 937 , 1369 , 1207 , 1340 , 1379 , 501 , 1587 , 1624 , 1607 , 733 , 934 , 1258 , 827 , 1366 , 882 , 1824 , 1676 , 255 , 1797 , 976 , 494 , 1581 , 1419 , 1473 , 1681 , 743 , 1147 , 847 , 1634 , 760 , 1199 , 1392 , 876 , 499 , 922 , 702 , 1244 , 1444 , 933 , 1338 , 1257 , 1365 , 731 , 881 , 1566 , 823 , 1728 , 493 , 1795 , 1579 , 1415 , 1633 , 1674 , 1254 , 1143 , 1358 , 921 , 1620 , 1183 , 1432 , 701 , 1337 , 931 , 478 , 1565 , 1363 , 727 , 874 , 815 , 918 , 1242 , 1442 , 491 , 968 , 756 , 1575 , 1334 , 1673 , 1253 , 1135 , 1357 , 638 , 1808 , 699 , 1384 , 1670 , 477 , 860 , 1563 , 1592 , 873 , 719 , 1618 , 917 , 1241 , 799 , 1441 , 487 , 1212 , 1333 , 1251 , 944 , 1119 , 1355 , 870 , 637 , 1428 , 1238 , 754 , 695 , 446 , 1612 , 1669 , 910 , 475 , 1559 , 1617 , 1326 , 915 , 1696 , 858 , 1331 , 504 , 964 , 1087 , 748 , 1351 , 869 , 635 , 1237 , 753 , 1210 , 687 , 445 , 909 , 1667 , 1551 , 471 , 1380 , 1325 , 1588 , 1426 , 1800 , 1610 , 857 , 1230 , 1264 , 828 , 867 , 631 , 382 , 1235 , 671 , 1209 , 443 , 907 , 1310 , 463 ,/* 962 , 1420 , 746 , 1323 , 854 , 1148 , 1425 , 936 , 1609 , 1368 , 1206 , 1229 , 623 , 1378 , 500 , 381 , 1586 , 439 , 903 , 1309 , 961 , 1319 , 745 , 853 , 1606 , 732 , 826 , 1205 , 1227 , 254 , 607 , 1796 , 1377 , 1580 , 379 , 1418 , 1472 , 431 , 1585 , 1680 , 742 , 1146 , 1307 , 846 , 851 , 1605 , 825 , 1198 , 498 , 1223 , 1203 , 575 , 253 , 375 , 415 , 1417 , 932 , 1256 , 1364 , 741 , 730 , 1145 , 1303 , 880 , 845 , 822 , 1603 , 492 , 1197 , 1794 , 497 , 1578 , 251 , 1414 , 1632 , 367 , 1142 , 1295 , 739 , 729 , 843 , 821 , 920 , 1182 , 700 , 1195 , 1793 , 1577 , 1336 , 247 , 1413 , 930 , 351 , 1564 , 1362 , 1141 , 726 , 839 , 814 , 819 , 490 , 1181 , 1191 , 1574 , 1672 , 1252 , 239 , 319 , 1411 , 929 , 1361 , 1134 , 1139 , 725 , 1356 , 813 , 489 , 698 , 1179 , 1573 , 223 , 476 , 1562 , 872 , 1133 , 718 , 723 , 916 , 1240 , 798 , 1440 , 811 , 486 , 697 , 1175 , 1332 , 1571 , 1250 , 191 , 1561 , 1118 , 1354 , 1131 , 717 , 636 , 797 , 807 , 485 , 694 , 1167 , 127 , 1668 , 1249 , 474 , 1558 , 1117 , 1127 , 1353 , 715 , 1616 , 914 , 795 , 483 , 693 , 1330 , 473 , 1086 , 1557 , 1115 , 1350 , 868 , 711 , 634 , 791 , 913 , 1236 , 752 , 686 , 444 , 691 , 908 , 1329 , 1666 , 1550 , 470 , 1085 , 1555 , 1111 , 1349 , 1324 , 633 , 783 , 685 , 856 , 1665 , 1549 , 469 , 1083 , 1103 , 1347 , 866 , 630 , 1234 , 670 , 1208 , 683 , 442 , 906 , 462 , 1547 , 467 , 1079 , 1322 , 865 , 629 , 1233 , 1424 , 669 , 679 , 441 , 905 , 1608 , 461 , 1071 , 1543 , 1228 , 1321 , 622 , 627 , 380 , 667 , 438 , 902 , 1308 , 1055 , 459 , 960 , 1318 , 744 , 621 , 852 , 663 , 437 , 901 , 455 , 1204 , 1226 , 606 , 1317 , 619 , 1376 , 378 , 655 , 430 , 435 , 1584 , 899 , 1306 , 1225 , 605 , 1315 , 615 , 850 , 1604 , 377 , 824 , 429 , 1305 , 1222 , 1202 , 574 , 603 , 252 , 849 , 374 , 414 , 1416 , 427 , 740 , 1144 , 1302 , 844 , 1221 , 1201 , 573 , 599 , 1602 , 373 , 413 , 423 , 1196 , 1301 , 496 , 1219 , 571 , 591 , 250 , 1601 , 366 , 371 , 411 , 1294 , 738 , 1299 , 728 , 842 , 567 , 249 , 820 , 365 , 407 , 1194 , 1293 , 737 , 1792 , 841 , 559 , 1576 , 246 , 1412 , 350 , 363 , 399 , 1140 , 1193 , 1291 , 543 , 838 , 245 , 818 , 349 , 359 , 1180 , 1190 , 1287 , 837 , 238 , 243 , 318 , 817 , 1410 , 347 , 928 , 1360 , 1138 , 724 , 1189 , 835 , 237 , 812 , 317 , 1409 , 343 , 488 , 1137 , 1178 , 1187 , 1572 , 222 , 235 , 315 , 335 , 1132 , 722 , 1177 , 221 , 231 , 311 , 810 , 696 , 1174 , 721 , 1570 , 190 , 219 , 303 , 809 , 1560 , 1130 , 716 , 1173 , 1569 , 189 , 215 , 796 , 287 , 806 , 484 , 1129 , 1166 , 1171 , 126 , 187 , 1248 , 207 , 805 , 1116 , 1126 , 1352 , 1165 , 714 , 125 , 183 , 794 , 803 , 482 , 1125 , 692 , 1163 , 713 , 123 , 175 , 793 , 472 , 1556 , 1114 , 1123 , 481 , 1159 , 710 , 119 , 159 , 790 , 912 , 1113 , 690 , 709 , 111 , 1328 , 789 , 1084 , 1554 , 1110 , 1348 , 689 , 707 , 95 , 632 , 782 , 787 , 1553 , 1109 , 684 , 63 , 1664 , 781 , 1548 , 468 , 1082 , 1102 , 1107 , 1346 , 779 , 1081 , 1101 , 682 , 1345 , 775 , 1546 , 466 , 1078 , 1099 , 681 , 864 , 628 , 1545 , 465 , 1077 , 1232 , 1095 , 668 , 678 , 440 , 904 , 460 , 1070 , 1542 , 1075 , 677 , 1320 , 626 , 1069 , 1541 , 666 , 675 , 625 , 1054 , 1067 , 458 , 1539 , 665 , 620 , 1053 , 1063 , 457 , 662 , 436 , 900 , 1051 , 454 , 661 , 1316 , 618 , 1047 , 453 , 654 , 659 , 434 , 617 , 898 , 1039 , 451 , 653 , 1224 , 433 , 604 , 1314 , 614 , 897 , 376 , 651 , 428 , 1313 , 613 , 1304 , 647 , 602 , 611 , 848 , 426 , 601 , 1220 , 1200 , 572 , 425 , 598 , 372 , 412 , 422 , 597 , 1300 , 1218 , 421 , 570 , 590 , 595 , 1600 , 370 , 410 , 1217 , 419 , 569 , 589 , 1298 , 369 , 409 , 566 , 587 , 248 , 1297 , 364 , 406 , 565 , 583 , 1292 , 736 , 405 , 840 , 558 , 563 , 362 , 398 , 403 , 557 , 1192 , 1290 , 361 , 397 , 542 , 555 , 244 , 1289 , 348 , 358 , 395 , 541 , 551 , 1286 , 357 , 391 , 539 , 836 , 242 , 816 , 1285 , 346 , 355 , 535 , 241 , 1188 , 1283 , 345 , 527 , 834 , 236 , 316 , 1408 , 342 , 833 , 1136 , 1186 , 341 , 234 , 314 , 1185 , 334 , 339 , 233 , 1176 , 313 , 333 , 220 , 230 , 310 , 331 , 229 , 720 , 309 , 327 , 218 , 227 , 302 , 307 , 808 , 217 , 301 , 1172 , 1568 , 188 , 214 , 286 , 299 , 213 , 1128 , 285 , 295 , 1170 , 186 , 206 , 211 , 283 , 1169 , 804 , 185 , 205 , 279 , 1164 , 124 , 182 , 203 , 271 , 802 , 181 , 199 , 1124 , 1162 , 712 , 801 , 122 , 174 , 179 , 1161 , 792 , 121 , 173 , 1122 , 480 , 1158 , 118 , 158 , 171 , 1121 , 1157 , 117 , 157 , 167 , 1112 , 1155 , 708 , 110 , 115 , 155 , 788 , 109 , 151 , 688 , 706 , 94 , 107 , 143 , 786 , 705 , 93 , 103 , 1552 , 1108 , 785 , 62 , 91 , 780 , 61 , 87 , 1106 , 59 , 79 , 1105 , 778 , 55 , 1080 , 1100 , 777 , 1344 , 47 , 774 , 31 , 1098 , 773 , 680 , 1097 , 771 , 1544 , 464 , 1076 , 1094 , 1093 , 1074 , 1091 , 676 , 1073 , 1068 , 1540 , 674 , 624 , 673 , 1066 , 1538 , 664 , 1065 , 1537 , 1052 , 1062 , 456 , 1061 , 1050 , 1059 , 660 , 1049 , 1046 , 452 , 658 , 1045 , 616 , 657 , 1038 , 1043 , 450 , 652 , 1037 , 432 , 449 , 896 , 1035 , 650 , 1031 , 1312 , 612 , 649 , 646 , 610 , 645 , 609 , 643 , 600 , 424 , 596 , 420 , 594 , 593 , 1216 , 418 , 568 , 588 , 417 , 368 , 408 , 586 , 585 , 1296 , 564 , 582 , 581 , 404 , 562 , 579 , 561 , 402 , 556 , 401 , 360 , 396 , 554 , 553 , 1288 , 394 , 540 , 550 , 393 , 549 , 356 , 390 , 538 , 547 , 389 , 537 , 1284 , 354 , 387 , 534 , 353 , 533 , 240 , 1282 , 344 , 526 , 531 , 1281 , 525 , 523 , 832 , 519 , 340 , 1184 , 338 , 337 , 232 , 312 , 332 , 330 , 329 , 228 , 308 , 326 , 325 , 226 , 306 , 323 , 225 , 305 , 216 , 300 , 298 , 297 , 212 , 284 , 294 , 293 , 210 , 282 , 291 , 1168 , 209 , 281 , 184 , 204 , 278 , 277 , 202 , 270 , 275 , 201 , 269 , 180 , 198 , 267 , 800 , 197 , 263 , 178 , 195 , 1160 , 177 , 120 , 172 , 170 , 1120 , 1156 , 169 , 116 , 156 , 166 , 1154 , 165 , 114 , 1153 , 154 , 163 , 113 , 153 , 108 , 150 , 149 , 106 , 142 , 147 , 105 , 141 , 704 , 92 , 102 , 139 , 101 , 135 , 784 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 1104 , 54 , 75 , 53 , 71 , 776 , 46 , 51 , 45 , 30 , 43 , 29 , 39 , 772 , 27 , 23 , 1096 , 770 , 15 , 769 , 1092 , 1090 , 1089 , 1072 , 672 , 1064 , 1536 , 1060 , 1058 , 1057 , 1048 , 1044 , 656 , 1042 , 1041 , 1036 , 448 , 1034 , 1033 , 1030 , 1029 , 648 , 1027 , 644 , 608 , 642 , 641 , 592 , 416 , 584 , 580 , 578 , 577 , 560 , 400 , 552 , 392 , 548 , 546 , 545 , 388 , 536 , 386 , 385 , 352 , 532 , 530 , 1280 , 529 , 524 , 522 , 521 , 518 , 517 , 515 , 336 , 328 , 324 , 322 , 321 , 224 , 304 , 296 , 292 , 290 , 289 , 208 , 280 , 276 , 274 , 273 , 200 , 268 , 266 , 265 , 196 , 262 , 261 , 194 , 259 , 193 , 176 , 168 , 164 , 1152 , 162 , 161 , 112 , 152 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 , 133 , 98 , 131 , 97 , 88 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 768 , 11 , 7 , 1088 , 1056 , 1040 , 1032 , 1028 , 1026 , 1025 , 640 , 576 , 544 , 384 , 528 , 520 , 516 , 514 , 513 , 320 , 288 , 272 , 264 , 260 , 258 , 257 , 192 , 160 , 144 , 136 , 132 , 130 , 129 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 1024 , 512 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , */// @0dB 
	//1023 , 1022 , 1021 , 1019 , 1015 , 1007 , 1020 , 1018 , 991 , 1017 , 1014 , 1013 , 959 , 1006 , 1011 , 1005 , 895 , 1003 , 990 , 1016 , 989 , 767 , 999 , 1012 , 958 , 987 , 1010 , 957 , 511 , 983 , 1004 , 1009 , 894 , 955 , 975 , 1002 , 893 , 951 , 1001 , 988 , 891 , 766 , 998 , 943 , 765 , 986 , 997 , 887 , 927 , 985 , 763 , 956 , 510 , 995 , 879 , 982 , 509 , 759 , 1008 , 863 , 954 , 981 , 974 , 507 , 751 , 953 , 831 , 979 , 892 , 950 , 973 , 503 , 735 , 1000 , 949 , 971 , 890 , 495 , 703 , 942 , 947 , 967 , 889 , 479 , 639 , 941 , 764 , 996 , 886 , 447 , 926 , 939 , 885 , 383 , 984 , 762 , 994 , 925 , 878 , 935 , 883 , 255 , 761 , 993 , 923 , 877 , 508 , 758 , 919 , 862 , 875 , 980 , 757 , 911 , 861 , 871 , 506 , 750 , 755 , 952 , 830 , 978 , 859 , 505 , 749 , 829 , 977 , 855 , 972 , 502 , 734 , 747 , 827 , 847 , 501 , 733 , 743 , 948 , 823 , 970 , 494 , 499 , 702 , 731 , 815 , 969 , 493 , 701 , 727 , 946 , 799 , 966 , 888 , 478 , 491 , 638 , 699 , 719 , 945 , 965 , 940 , 477 , 487 , 637 , 695 , 963 , 446 , 475 , 635 , 687 , 445 , 471 , 938 , 631 , 671 , 884 , 382 , 443 , 463 , 937 , 623 , 924 , 381 , 439 , 934 , 607 , 882 , 254 , 379 , 431 , 575 , 933 , 881 , 760 , 992 , 253 , 922 , 375 , 415 , 931 , 876 , 251 , 921 , 367 , 247 , 918 , 351 , 874 , 239 , 319 , 917 , 873 , 756 , 223 , 910 , 915 , 860 , 870 , 191 , 909 , 869 , 754 , 127 , 907 , 858 , 867 , 753 , 504 , 903 , 857 , 748 , 828 , 976 , 854 , 853 , 746 , 826 , 846 , 851 , 745 , 500 , 825 , 845 , 732 , 742 , 822 , 843 , 741 , 498 , 821 , 839 , 730 , 739 , 497 , 814 , 819 , 729 , 968 , 492 , 813 , 700 , 726 , 798 , 811 , 725 , 490 , 797 , 807 , 698 , 718 , 723 , 944 , 489 , 795 , 697 , 717 , 964 , 476 , 486 , 636 , 791 , 694 , 715 , 485 , 783 , 693 , 711 , 962 , 474 , 483 , 634 , 686 , 691 , 961 , 473 , 633 , 685 , 444 , 470 , 630 , 670 , 683 , 469 , 629 , 669 , 679 , 442 , 462 , 467 , 936 , 622 , 627 , 667 , 441 , 461 , 621 , 663 , 380 , 438 , 459 , 606 , 619 , 655 , 437 , 455 , 605 , 615 , 378 , 430 , 435 , 574 , 603 , 932 , 377 , 880 , 429 , 573 , 599 , 252 , 374 , 414 , 427 , 571 , 591 , 930 , 373 , 413 , 423 , 567 , 929 , 250 , 920 , 366 , 371 , 411 , 559 , 249 , 365 , 407 , 543 , 246 , 350 , 363 , 399 , 245 , 349 , 359 , 238 , 243 , 318 , 347 , 916 , 872 , 237 , 317 , 343 , 222 , 235 , 315 , 335 , 914 , 221 , 231 , 311 , 913 , 190 , 219 , 303 , 908 , 189 , 868 , 215 , 287 , 126 , 187 , 207 , 906 , 125 , 183 , 866 , 752 , 905 , 123 , 175 , 865 , 902 , 119 , 159 , 856 , 901 , 111 , 899 , 95 , 63 , 852 , 850 , 744 , 849 , 824 , 844 , 842 , 740 , 841 , 820 , 838 , 738 , 837 , 737 , 496 , 818 , 835 , 728 , 817 , 812 , 810 , 724 , 809 , 796 , 806 , 722 , 805 , 721 , 488 , 794 , 803 , 696 , 716 , 793 , 790 , 714 , 789 , 713 , 484 , 782 , 787 , 692 , 710 , 781 , 709 , 482 , 779 , 690 , 707 , 481 , 775 , 960 , 689 , 472 , 632 , 684 , 682 , 681 , 468 , 628 , 668 , 678 , 677 , 466 , 626 , 666 , 675 , 465 , 625 , 665 , 440 , 460 , 620 , 662 , 661 , 458 , 618 , 654 , 659 , 457 , 617 , 653 , 436 , 454 , 604 , 614 , 651 , 453 , 613 , 647 , 434 , 451 , 602 , 611 , 433 , 601 , 376 , 428 , 572 , 598 , 597 , 426 , 570 , 590 , 595 , 425 , 569 , 589 , 372 , 412 , 422 , 566 , 587 , 928 , 421 , 565 , 583 , 370 , 410 , 419 , 558 , 563 , 369 , 409 , 557 , 248 , 364 , 406 , 542 , 555 , 405 , 541 , 551 , 362 , 398 , 403 , 539 , 361 , 397 , 535 , 244 , 348 , 358 , 395 , 527 , 357 , 391 , 242 , 346 , 355 , 241 , 345 , 236 , 316 , 342 , 341 , 234 , 314 , 334 , 339 , 233 , 313 , 333 , 220 , 230 , 310 , 331 , 912 , 229 , 309 , 327 , 218 , 227 , 302 , 307 , 217 , 301 , 188 , 214 , 286 , 299 , 213 , 285 , 295 , 186 , 206 , 211 , 283 , 185 , 205 , 279 , 124 , 182 , 203 , 271 , 904 , 181 , 199 , 122 , 174 , 179 , 864 , 121 , 173 , 118 , 158 , 171 , 900 , 117 , 157 , 167 , 110 , 115 , 155 , 898 , 109 , 151 , 94 , 107 , 897 , 143 , 93 , 103 , 62 , 91 , 61 , 87 , 59 , 79 , 55 , 47 , 31 , 848 , 840 , 836 , 736 , 834 , 833 , 816 , 808 , 804 , 720 , 802 , 801 , 792 , 788 , 712 , 786 , 785 , 780 , 708 , 778 , 777 , 706 , 480 , 774 , 705 , 688 , 773 , 771 , 680 , 676 , 674 , 464 , 673 , 624 , 664 , 660 , 658 , 657 , 456 , 616 , 652 , 650 , 649 , 452 , 612 , 646 , 645 , 450 , 610 , 643 , 449 , 609 , 432 , 600 , 596 , 594 , 593 , 424 , 568 , 588 , 586 , 585 , 420 , 564 , 582 , 581 , 418 , 562 , 579 , 417 , 561 , 368 , 408 , 556 , 554 , 553 , 404 , 540 , 550 , 549 , 402 , 538 , 547 , 401 , 537 , 360 , 396 , 534 , 533 , 394 , 526 , 531 , 393 , 525 , 356 , 390 , 523 , 389 , 519 , 354 , 387 , 353 , 240 , 344 , 340 , 338 , 337 , 232 , 312 , 332 , 330 , 329 , 228 , 308 , 326 , 325 , 226 , 306 , 323 , 225 , 305 , 216 , 300 , 298 , 297 , 212 , 284 , 294 , 293 , 210 , 282 , 291 , 209 , 281 , 184 , 204 , 278 , 277 , 202 , 270 , 275 , 201 , 269 , 180 , 198 , 267 , 197 , 263 , 178 , 195 , 177 , 120 , 172 , 170 , 169 , 116 , 156 , 166 , 165 , 114 , 154 , 163 , 113 , 153 , 108 , 150 , 149 , 106 , 896 , 142 , 147 , 105 , 141 , 92 , 102 , 139 , 101 , 135 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 71 , 46 , 51 , 45 , 30 , 43 , 29 , 39 , 27 , 23 , 15 , 832 , 800 , 784 , 776 , 704 , 772 , 770 , 769 , 672 , 656 , 648 , 644 , 642 , 448 , 641 , 608 , 592 , 584 , 580 , 578 , 416 , 577 , 560 , 552 , 548 , 546 , 400 , 545 , 536 , 532 , 530 , 529 , 392 , 524 , 522 , 521 , 388 , 518 , 517 , 386 , 515 , 385 , 352 , 336 , 328 , 324 , 322 , 321 , 224 , 304 , 296 , 292 , 290 , 289 , 208 , 280 , 276 , 274 , 273 , 200 , 268 , 266 , 265 , 196 , 262 , 261 , 194 , 259 , 193 , 176 , 168 , 164 , 162 , 161 , 112 , 152 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 , 133 , 98 , 131 , 97 , 88 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 768 , 640 , 576 , 544 , 528 , 520 , 516 , 514 , 384 , 513 , 320 , 288 , 272 , 264 , 260 , 258 , 257 , 192 , 160 , 512 , 256 , 144 , 136 , 132 , 130 , 129 , 128 , 96 , 80 , 72 , 68 , 66 , 65 , 64 , 48 , 40 , 36 , 34 , 33 , 32 , 24 , 20 , 18 , 17 , 16 , 12 , 10 , 9 , 8 , 6 , 5 , 4 , 3 , 2 , 1 , 0 , //@5dB
	//1023 , 1022 , 1021 , 1019 , 1015 , 1007 , 991 , 959 , 1020 , 1018 , 1017 , 1014 , 1013 , 1011 , 1006 , 1005 , 1003 , 895 , 990 , 989 , 999 , 987 , 983 , 958 , 767 , 957 , 955 , 1016 , 1012 , 975 , 1010 , 1009 , 1004 , 1002 , 951 , 894 , 1001 , 511 , 893 , 988 , 998 , 986 , 997 , 891 , 943 , 985 , 982 , 995 , 766 , 956 , 887 , 981 , 954 , 765 , 927 , 974 , 953 , 1008 , 979 , 763 , 879 , 950 , 973 , 1000 , 510 , 892 , 949 , 971 , 759 , 996 , 509 , 890 , 863 , 942 , 984 , 947 , 994 , 889 , 941 , 967 , 507 , 751 , 886 , 980 , 764 , 831 , 993 , 926 , 939 , 885 , 503 , 952 , 735 , 978 , 762 , 925 , 878 , 972 , 935 , 883 , 977 , 495 , 761 , 948 , 877 , 923 , 703 , 970 , 758 , 508 , 862 , 479 , 875 , 946 , 919 , 969 , 757 , 639 , 888 , 861 , 940 , 966 , 506 , 750 , 945 , 871 , 755 , 447 , 911 , 830 , 992 , 859 , 965 , 505 , 749 , 938 , 884 , 502 , 829 , 734 , 383 , 963 , 747 , 855 , 924 , 937 , 501 , 934 , 733 , 827 , 882 , 976 , 255 , 494 , 743 , 760 , 847 , 876 , 499 , 922 , 702 , 933 , 881 , 731 , 823 , 493 , 921 , 701 , 931 , 478 , 727 , 874 , 815 , 918 , 491 , 968 , 756 , 638 , 699 , 860 , 477 , 873 , 719 , 917 , 799 , 487 , 944 , 870 , 637 , 754 , 695 , 446 , 910 , 475 , 915 , 858 , 964 , 504 , 748 , 869 , 635 , 753 , 445 , 687 , 909 , 471 , 857 , 828 , 867 , 631 , 382 , 671 , 443 , 907 , 962 , 463 , 746 , 854 , 936 , 623 , 500 , 381 , 439 , 903 , 961 , 745 , 853 , 732 , 826 , 254 , 607 , 379 , 431 , 742 , 846 , 851 , 825 , 498 , 575 , 253 , 932 , 375 , 415 , 741 , 880 , 730 , 845 , 822 , 492 , 497 , 251 , 367 , 739 , 729 , 843 , 821 , 920 , 700 , 930 , 247 , 351 , 726 , 839 , 814 , 819 , 490 , 239 , 929 , 319 , 725 , 813 , 489 , 698 , 223 , 476 , 872 , 718 , 723 , 916 , 798 , 811 , 486 , 697 , 191 , 717 , 636 , 797 , 807 , 485 , 694 , 127 , 474 , 715 , 914 , 795 , 483 , 693 , 473 , 868 , 711 , 634 , 913 , 752 , 791 , 444 , 686 , 691 , 908 , 470 , 633 , 783 , 685 , 856 , 469 , 866 , 630 , 670 , 683 , 442 , 906 , 462 , 467 , 865 , 629 , 669 , 679 , 441 , 905 , 461 , 622 , 627 , 380 , 667 , 438 , 902 , 459 , 960 , 744 , 852 , 621 , 663 , 437 , 901 , 455 , 606 , 619 , 378 , 655 , 430 , 435 , 899 , 605 , 615 , 850 , 824 , 377 , 429 , 574 , 603 , 252 , 849 , 374 , 414 , 427 , 740 , 844 , /*573 , 599 , 373 , 413 , 423 , 496 , 571 , 591 , 250 , 366 , 371 , 411 , 738 , 728 , 842 , 567 , 820 , 249,/* 365 , 407 , 737 , 841 , 559 , 246 , 350 , 363 , 399 , 838 , 543 , 818 , 245 , 349 , 359 , 837 , 238 , 928 , 817 , 243 , 318 , 347 , 724 , 835 , 812 , 237 , 317 , 343 , 488 , 222/* , 235 , 315 , 335 , 722 , 221 , 231 , 810 , 311 , 696 , 721 , 190 , 219 , 809 , 303 , 716 , 189 , 796 , 215 , 806 , 287 , 484 , 126 , 187 , 207 , 805 , 714 , 125 , 183 , 794 , 803 , 482 , 692 , 713 , 123 , 175 , 793 , 472 , 481 , 710 , 119 , 159 , 912 , 790 , 690 , 709 , 111 , 789 , 689 , 707 , 95 , 632 , 782 , 787 , 684 , 63 , 781 , 468 , 779 , 682 , 775 , 466 , 681 , 864 , 628 , 465 , 668 , 678 , 440 , 904 , 460 , 677 , 626 , 666 , 675 , 625 , 458 , 665 , 620 , 457 , 662 , 436 , 900 , 454 , 661 , 618 , 453 , 654 , 659 , 434 , 617 , 898 , 451 , 653 , 433 , 604 , 614 , 897 , 376 , 651 , 428 , 613 , 647 , 602 , 611 , 848 , 426 , 601 , 572 , 425 , 598 , 372 , 412 , 422 , 597 , 421 , 570 , 590 , 595 , 370 , 410 , 419 , 569 , 589 , 369 , 409 , 566 , 587 , 248 , 364 , 406 , 565 , 583 , 736 , 840 , 405 , 558 , 563 , 362 , 398 , 403 , 557 , 361 , 397 , 542 , 555 , 244 , 348 , 358 , 395 , 541 , 551 , 357 , 836 , 391 , 539 , 816 , 242 , 346 , 355 , 535 , 241 , 345 , 527 , 834 , 236 , 316 , 342 , 833 , 341 , 234 , 314 , 334 , 339 , 233 , 313 , 333 , 220 , 230 , 310 , 331 , 720 , 229 , 309 , 327 , 218 , 227 , 808 , 302 , 307 , 217 , 301 , 188 , 214 , 286 , 299 , 213 , 285 , 295 , 186 , 206 , 211 , 283 , 804 , 185 , 205 , 279 , 124 , 182 , 203 , 271 , 802 , 181 , 199 , 712 , 801 , 122 , 174 , 179 , 792 , 121 , 173 , 480 , 118 , 158 , 171 , 117 , 157 , 167 , 708 , 110 , 115 , 155 , 788 , 109 , 151 , 688 , 706 , 94 , 107 , 143 , 786 , 705 , 93 , 103 , 785 , 62 , 91 , 780 , 61 , 87 , 59 , 79 , 778 , 55 , 777 , 47 , 774 , 31 , 773 , 680 , 771 , 464 , 676 , 674 , 624 , 673 , 664 , 456 , 660 , 452 , 658 , 616 , 657 , 450 , 652 , 432 , 449 , 896 , 650 , 612 , 649 , 646 , 610 , 645 , 609 , 643 , 600 , 424 , 596 , 420 , 594 , 593 , 418 , 568 , 588 , 417 , 368 , 408 , 586 , 585 , 564 , 582 , 581 , 404 , 562 , 579 , 561 , 402 , 556 , 401 , 360 , 396 , 554 , 553 , 394 , 540 , 550 , 393 , 549 , 356 , 390 , 538 , 547 , 389 , 537 , 354 , 387 , 534 , 353 , 533 , 240 , 344 , 526 , 531 , 525 , 523 , 832 , 519 , 340 , 338 , 337 , 232 , 312 , 332 , 330 , 329 , 228 , 308 , 326 , 325 , 226 , 306 , 323 , 225 , 305 , 216 , 300 , 298 , 297 , 212 , 284 , 294 , 293 , 210 , 282 , 291 , 209 , 281 , 184 , 204 , 278 , 277 , 202 , 270 , 275 , 201 , 269 , 180 , 198 , 267 , 800 , 197 , 263 , 178 , 195 , 177 , 120 , 172 , 170 , 169 , 116 , 156 , 166 , 165 , 114 , 154 , 163 , 113 , 153 , 108 , 150 , 149 , 106 , 142 , 147 , 105 , 704 , 141 , 92 , 102 , 139 , 101 , 784 , 135 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 71 , 776 , 46 , 51 , 45 , 30 , 43 , 29 , 39 , 772 , 27 , 23 , 770 , 15 , 769 , 672 , 656 , 448 , 648 , 644 , 608 , 642 , 641 , 592 , 416 , 584 , 580 , 578 , 577 , 560 , 400 , 552 , 392 , 548 , 546 , 545 , 388 , 536 , 386 , 385 , 352 , 532 , 530 , 529 , 524 , 522 , 521 , 518 , 517 , 515 , 336 , 328 , 324 , 322 , 321 , 224 , 304 , 296 , 292 , 290 , 289 , 208 , 280 , 276 , 274 , 273 , 200 , 268 , 266 , 265 , 196 , 262 , 261 , 194 , 259 , 193 , 176 , 168 , 164 , 162 , 161 , 112 , 152 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 , 133 , 98 , 131 , 97 , 88 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 768 , 11 , 7 , 640 , 576 , 544 , 384 , 528 , 520 , 516 , 514 , 513 , 320 , 288 , 272 , 264 , 260 , 258 , 257 , 192 , 160 , 144 , 136 , 132 , 130 , 129 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 512 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@3dB
	//1023 , 1022 , 1021 , 1019 , 1015 , 1007 , 991 , 959 , 895 , 1020 , 1018 , 1017 , 1014 , 1013 , 1011 , 1006 , 1005 , 1003 , 999 , 990 , 989 , 987 , 983 , 767 , 958 , 957 , 955 , 975 , 951 , 894 , 893 , 511 , 943 , 891 , 1016 , 1012 , 1010 , 1009 , 1004 , 1002 , 1001 , 998 , 997 , 988 , 887 , 986 , 985 , 995 , 982 , 927 , 981 , 766 , 765 , 956 , 979 , 954 , 879 , 974 , 953 , 763 , 973 , 950 , 949 , 971 , 759 , 892 , 510 , 947 , 863 , 942 , 890 , 1008 , 509 , 1000 , 941 , 967 , 889 , 996 , 886 , 751 , 507 , 984 , 994 , 939 , 926 , 980 , 885 , 831 , 993 , 764 , 978 , 925 , 878 , 503 , 952 , 762 , 935 , 883 , 972 , 735 , 977 , 877 , 948 , 923 , 761 , 970 , 758 , 495 , 946 , 862 , 875 , 969 , 703 , 508 , 919 , 757 , 940 , 966 , 888 , 945 , 861 , 750 , 506 , 871 , 479 , 965 , 938 , 755 , 884 , 911 , 830 , 859 , 992 , 749 , 639 , 505 , 924 , 937 , 963 , 502 , 829 , 934 , 882 , 734 , 447 , 747 , 976 , 855 , 876 , 922 , 501 , 760 , 933 , 881 , 733 , 827 , 494 , 743 , 921 , 874 , 847 , 383 , 499 , 968 , 702 , 918 , 931 , 756 , 731 , 823 , 493 , 944 , 860 , 873 , 701 , 917 , 870 , 478 , 255 , 964 , 727 , 754 , 491 , 815 , 910 , 858 , 748 , 638 , 699 , 504 , 869 , 915 , 477 , 753 , 936 , 719 , 962 , 487 , 909 , 857 , 799 , 828 , 637 , 446 , 867 , 695 , 475 , 746 , 854 , 961 , 500 , 907 , 932 , 635 , 880 , 445 , 732 , 745 , 826 , 853 , 471 , 687 , 742 , 903 , 920 , 846 , 382 , 443 , 498 , 631 , 825 , 851 , 930 , 463 , 671 , 730 , 822 , 741 , 845 , 492 , 381 , 497 , 872 , 439 , 623 , 929 , 729 , 700 , 916 , 821 , 739 , 843 , 254 , 379 , 726 , 431 , 490 , 607 , 814 , 819 , 253 , 839 , 375 , 698 , 725 , 868 , 914 , 476 , 415 , 489 , 575 , 813 , 752 , 718 , 486 , 908 , 251 , 697 , 367 , 723 , 856 , 913 , 798 , 811 , 636 , 866 , 694 , 474 , 717 , 485 , 247 , 351 , 797 , 960 , 807 , 906 , 865 , 693 , 473 , 715 , 483 , 239 , 444 , 634 , 319 , 795 , 744 , 852 , 470 , 905 , 686 , 691 , 711 , 223 , 633 , 791 , 902 , 469 , 685 , 442 , 630 , 191 , 824 , 850 , 783 , 901 , 462 , 670 , 467 , 683 , 441 , 127 , 629 , 740 , 844 , 849 , 380 , 899 , 461 , 496 , 669 , 679 , 438 , 622 ,/* 627 , 928 , 728 , 459 , 667 , 437 , 820 , 738 , 621 , 842 , 378 , 455 , 663 , 430 , 435 , 737 , 606 , 619 , 841 , 377 , 655 , 818 , 429 , 605/* , 615 , 252 , 838 , 374 , 724 , 414 , 817 , 427 , 488 , 574 , 812 , 603 , 837 , 373 , 413 , 423 , 573 , 599 , 250 , 835 , 696 , 366 , 371 , 722 , 912 , 411 , 571 , 810 , 591 , 249 , 365 , 721 , 407 , 716 , 567 , 484 , 809 , 246 , 350 , 363 , 796 , 399 , 559 , 806 , 245 , 349 , 359 , 864 , 692 , 472 , 714 , 543 , 805 , 482 , 238 , 243 , 318 , 347 , 794 , 713 , 803 , 481 , 237 , 904 , 317 , 343 , 690 , 793 , 710 , 222 , 235 , 315 , 632 , 335 , 790 , 689 , 709 , 221 , 231 , 311 , 468 , 684 , 789 , 707 , 190 , 219 , 303 , 782 , 787 , 900 , 189 , 215 , 287 , 466 , 682 , 781 , 440 , 126 , 628 , 187 , 207 , 465 , 681 , 848 , 779 , 898 , 125 , 183 , 460 , 668 , 678 , 775 , 897 , 123 , 626 , 175 , 677 , 119 , 159 , 625 , 458 , 666 , 675 , 436 , 111 , 620 , 457 , 665 , 95 , 454 ,/* 662 , 434 , 736 , 63 , 618 , 453 , 840 , 661 , 433 , 376 , 617 , 451 , 654 , 659 , 428 , 604 , 614 , 653 , 613 , 651 , 816 , 426 , 602 , 611 , 647 , 836 , 425 , 372 , 601 , 412 , 422 , 572 , 598 , 834 , 421 , 370 , 597 , 833 , 410 , 419 , 369 , 570 , 590 , 595 , 248 , 409 , 364 , 720 , 569 , 589 , 406 , 566 , 587 , 808 , 405 , 362 , 565 , 583 , 398 , 403 , 361 , 558 , 563 , 244 , 397 , 348 , 358 , 557 , 395 , 357 , 542 , 555 , 804 , 242 , 391 , 346 , 355 , 541 , 551 , 241 , 345 , 712 , 539 , 802 , 480 , 236 , 316 , 342 , 535 , 801 , 792 , 341 , 527 , 234 , 314 , 334 , 339 , 688 , 233 , 313 , 333 , 708 , 220 , 230 , 310 , 331 , 229 , 788 , 309 , 327 , 706 , 218 , 227 , 302 , 307 , 705 , 217 , 786 , 301 , 188 , 214 , 785 , 286 , 299 , 213 , 780 , 285 , 295 , 186 , 206 , 211 , 283 , 464 , 680 , 185 , 205 , 778 , 279 , 124 , 182 , 203 , 271 , 777 , 181 , 199 , 774 , 896 , 122 , 174 , 179 , 773 , 676 , 121 , 173 , 771 , 118 , 158 , 624 , 171 , 674 , 117 , 157 , 167 , 673 , 110 , 115 , 155 , 456 , 664 , 109 , 151 , 94 , 107 , 143 , 93 , 103 , 62 , 91 , 452 , 660 , 61 , 87 , 432 , 59 , 79 , 616 , 450 , 55 , 658 , 449 , 47 , 657 , 31 , 652 , 612 , 650 , 649 , 610 , 646 , 609 , 424 , 645 , 600 , 643 , 420 , 596 , 832 , 418 , 368 , 594 , 417 , 593 , 408 , 568 , 588 , 586 , 585 , 404 , 564 , 582 , 581 , 402 , 360 , 562 , 579 , 401 , 561 , 396 , 556 , 394 , 356 , 554 , 393 , 553 , 390 , 354 , 540 , 550 , 389 , 353 , 549 , 240 , 387 , 344 , 538 , 547 , 537 , 534 , 800 , 533 , 340 , 526 , 531 , 525 , 338 , 523 , 337 , 519 , 232 , 312 , 332 , 330 , 329 , 228 , 308 , 326 , 325 , 226 , 306 , 323 , 704 , 225 , 305 , 216 , 300 , 784 , 298 , 297 , 212 , 284 , 294 , 293 , 210 , 282 , 291 , 209 , 281 , 184 , 204 , 278 , 277 , 202 , 270 , 776 , 275 , 201 , 269 , 180 , 198 , 267 , 197 , 263 , 178 , 195 , 772 , 177 , 120 , 172 , 770 , 769 , 170 , 169 , 116 , 156 , 166 , 165 , 672 , 114 , 154 , 163 , 113 , 153 , 108 , 150 , 149 , 106 , 142 , 147 , 105 , 141 , 92 , 102 , 139 , 101 , 135 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 71 , 448 , 46 , 51 , 656 , 45 , 30 , 43 , 29 , 39 , 27 , 23 , 15 , 648 , 608 , 644 , 642 , 641 , 416 , 592 , 584 , 580 , 578 , 400 , 577 , 560 , 392 , 552 , 388 , 352 , 548 , 386 , 385 , 546 , 545 , 536 , 532 , 530 , 529 , 524 , 522 , 521 , 336 , 518 , 517 , 515 , 328 , 324 , 322 , 321 , 224 , 304 , 296 , 292 , 290 , 289 , 208 , 280 , 276 , 274 , 273 , 200 , 268 , 266 , 265 , 196 , 262 , 261 , 194 , 259 , 193 , 176 , 768 , 168 , 164 , 162 , 161 , 112 , 152 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 , 133 , 98 , 131 , 97 , 88 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 640 , 576 , 384 , 544 , 528 , 520 , 516 , 514 , 513 , 320 , 288 , 272 , 264 , 260 , 258 , 257 , 192 , 160 , 144 , 136 , 132 , 130 , 129 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 512 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 ,*///@2dB
	//1023 , 1022 , 1021 , 1019 , 1015 , 1007 , 991 , 959 , 895 , 767 , 1020 , 1018 , 1017 , 1014 , 1013 , 1011 , 1006 , 1005 , 1003 , 999 , 990 , 989 , 987 , 983 , 975 , 958 , 957 , 955 , 951 , 943 , 511 , 894 , 893 , 891 , 887 , 927 , 879 , 766 , 765 , 763 , 863 , 759 , 1016 , 1012 , 1010 , 1009 , 1004 , 1002 , 1001 , 998 , 997 , 995 , 988 , 986 , 985 , 982 , 981 , 979 , 751 , 974 , 973 , 971 , 956 , 954 , 953 , 831 , 950 , 949 , 967 , 947 , 942 , 941 , 939 , 510 , 509 , 507 , 735 , 935 , 892 , 890 , 889 , 503 , 886 , 885 , 926 , 925 , 883 , 923 , 878 , 877 , 875 , 919 , 495 , 703 , 764 , 762 , 871 , 862 , 761 , 861 , 758 , 757 , 859 , 755 , 911 , 1008 , 1000 , 996 , 994 , 984 , 993 , 980 , 978 , 479 , 750 , 977 , 749 , 972 , 970 , 969 , 855 , 952 , 830 , 948 , 966 , 747 , 946 , 829 , 965 , 945 , 639 , 940 , 938 , 508 , 506 , 937 , 734 , 827 , 505 , 963 , 733 , 934 , 888 , 502 , 743 , 933 , 884 , 847 , 924 , 882 , 501 , 922 , 447 , 881 , 731 , 876 , 921 , 823 , 874 , 931 , 918 , 494 , 873 , 499 , 702 , 917 , 493 , 870 , 760 , 860 , 727 , 701 , 756 , 858 , 869 , 815 , 754 , 910 , 915 , 491 , 992 , 857 , 383 , 478 , 976 , 699 , 748 , 753 , 909 , 968 , 854 , 867 , 477 , 746 , 828 , 719 , 964 , 944 , 638 , 853 , 487 , 936 , 907 , 826 , 504 , 695 , 962 , 745 , 799 , 732 , 637 , 475 , 742 , 932 , 846 , 500 , 851 , 825 , 255 , 961 , 446 , 880 , 730 , 920 , 822 , 741 , 903 , 845 , 930 , 635 , 872 , 498 , 445 , 471 , 687 , 729 , 916 , 492 , 821 , 929 , 726 , 700 , 739 , 497 , 843 , 868 , 814 , 631 , 443 , 914 , 490 , 856 , 819 , 382 , 725 , 698 , 463 , 752 , 671 , 908 , 866 , 813 , 913 , 839 , 476 , 489 , 718 , 381 , 697 , 439 , 852 , 486 , 623 , 723 , 906 , 865 , 694 , 744 , 798 , 717 , 811 , 636 , 474 , 485 , 379 , 850 , 824 , 905 , 254 , 960 , 693 , 431 , 797 , 607 , 740 , 902 , 844 , 473 , 634 , 715 , 807 , 444 , 849 , 470 , 686 , 728 , 253 , 483 , 375 , 820 , 691 , 928 , 795 , 901 , 633 , 415 , 738 , 496 , 842 , 711 , 575 , 469 , 685 , 630 , 442 , 251 , 367 , 818 , 724 , 899 , 791 , 737 , 462 , 841 , 670 , 812 , 467 , 629 , /*683 , 441 , 912 , 838 , 488 , 247 , 817 , 380 , 696 , 438 , 461 , 351 , 622 , 722 , 669 , 783 , 864 , 837 , 627 , 679 , 810 , 716 , 239 , 437 , 459 , 621 , 721 , 484 , 378 , 667 , 319 , 904 , 692 , 835 , 430 , 796 , 809 , 606 , 472 , 435 , 714 , 619 , 377 , 223 , 806 , 455 , 663 , 848 , 429 , 252 , 482 , 374 , 605 , 690 , 794 , 713 , 900 , 805 , 615 , 632 , 414 , 191 , 655 , 481 , 710 , 373 , 427 , 574 , 468 , 684 , 603 , 689 , 250 , 793 , 803 , 366 , 413 , 127 , 898 , 790 , 709 , 736 , 840 , 371 , 573 , 423 , 599 , 249 , 466 , 628 , 682 , 365 , 440 , 411 , 246 , 897 , 789 , 707 , 816 , 571 , 460 , 591 , 350 , 668 , 465 , 782 , 681 , 363 , 407 , 245 , 787 , 836 , 626 , 567 , 678 , 349 , 238 , 781 , 436 , 359 , 458 , 620 , 399 , 720 , 243 , 666 , 625 , 559 , 677 , 318 , 347 , 834 , 237 , 779 , 808 , 457 ,/* 665 , 543 , 675 , 434 , 317 , 618 , 376 , 343 , 222 , 454 , 833 , 775 , 235 , 662 , 428 , 433 , 315 , 617 , 335 , 604 , 221 , 453 , 231 , 661 , 712 , 804 , 614 , 311 , 190 , 219 , 451 , 654 , 480 , 659 , 372 , 426 , 602 , 613 , 303 , 688 , 189 , 215 , 792 , 653 , 425 , 802 , 412 , 611 , 601 , 287 , 126 , 187 , 207 , 708 , 370 , 572 , 651 , 422 , 801 , 598 , 248 , 125 , 183 , 647 , 369 , 364 , 421 , 410 , 597 , 896 , 788 , 123 , 706 , 175 , 570 , 419 , 409 , 590 , 595 , 119 , 159 , 705 , 464 , 680 , 569 , 362 , 406 , 244 , 589 , 786 , 111 , 566 , 361 , 405 , 348 , 587 , 785 , 95 , 780 , 565 , 358 , 398 , 403 , 242 , 583 , 63 , 624 , 558 , 563 , 676 , 357 , 397 , 241 , 346 , 236 , 778 , 557 , 355 , 395 , 456 , 345 , 664 , 542 , 777 , 555 , 674 , 316 , 391 , 342 , 832 , 774 , 234 , 541 , 551 , 673 , 341 , 773 , 233 , 539 , 432 , 314 , 616 , 334 , 339 , 220 , 452 , 230 , 535 , 771 , 660 , 313 , 333 , 527 , 229 , 310 , 331 , 218 , 450 , 227 , 658 , 309 , 327 , 217 , 449 , 612 , 302 , 657 , 307 , 188 , 214 , 652 , 301 , 424 , 213 , 610 , 600 , 286 , 299 , 186 , 206 , 211 , 650 , 609 , 285 , 295 , 800 , 185 , 205 , 649 , 283 , 124 , 182 , 203 , 646 , 368 , 420 , 279 , 181 , 199 , 596 , 645 , 271 , 122 , 174 , 179 , 643 , 418 , 121 , 408 , 173 , 594 , 417 , 118 , 158 , 704 , 171 , 568 , 593 , 117 , 157 , 167 , 588 , 110 , 115 , 155 , 360 , 109 , 151 , 404 , 586 , 784 , 94 , 107 , 143 , 564 , 585 , 93 , 103 , 402 , 582 , 62 , 91 , 401 , 562 , 581 , 356 , 61 , 87 , 396 , 561 , 240 , 579 , 59 , 79 , 556 , 354 , 55 , 394 , 344 , 353 , 47 , 393 , 776 , 554 , 31 , 390 , 553 , 389 , 540 , 550 , 672 , 387 , 549 , 340 , 772 , 232 , 538 , 547 , 537 , 338 , 534 , 770 , 337 , 312 , 533 , 769 , 332 , 526 , 531 , 228 , 525 , 330 , 523 , 226 , 329 , 308 , 519 , 225 , 326 , 216 , 448 , 325 , 656 , 306 , 323 , 305 , 300 , 212 , 298 , 210 , 297 , 608 , 209 , 284 , 294 , 184 , 204 , 293 , 648 , 282 , 291 , 202 , 281 , 201 , 278 , 180 , 198 , 277 , 644 , 197 , 270 , 275 , 178 , 195 , 269 , 642 , 177 , 267 , 120 , 641 , 172 , 263 , 416 , 170 , 592 , 169 , 116 , 156 , 166 , 165 , 114 , 154 , 163 , 113 , 153 , 108 , 150 , 149 , 106 , 142 , 147 , 105 , 141 , 584 , 92 , 102 , 139 , 101 , 135 , 90 , 99 , 400 , 89 , 580 , 60 , 86 , 85 , 560 , 578 , 58 , 78 , 83 , 577 , 57 , 77 , 54 , 75 , 53 , 71 , 352 , 46 , 51 , 392 , 45 , 30 , 43 , 29 , 39 , 552 , 27 , 388 , 23 , 15 , 386 , 548 , 385 , 546 , 545 , 536 , 336 , 532 , 768 , 530 , 529 , 524 , 522 , 328 , 521 , 518 , 224 , 517 , 515 , 324 , 322 , 304 , 321 , 296 , 208 , 292 , 290 , 289 , 280 , 200 , 276 , 196 , 274 , 273 , 194 , 268 , 193 , 176 , 266 , 265 , 640 , 262 , 261 , 259 , 168 , 164 , 162 , 161 , 112 , 152 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 , 133 , 98 , 131 , 97 , 88 , 84 , 82 , 576 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 384 , 544 , 528 , 520 , 516 , 514 , 513 , 320 , 288 , 272 , 192 , 264 , 260 , 258 , 257 , 160 , 144 , 136 , 132 , 130 , 129 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 512 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@0dB
	//1023 , 1022 , 1021 , 1019 , 1015 , 1007 , 991 , 959 , 895 , 767 , 511 , 1020 , 1018 , 1017 , 1014 , 1013 , 1011 , 1006 , 1005 , 1003 , 999 , 990 , 989 , 987 , 983 , 975 , 958 , 957 , 955 , 951 , 943 , 927 , 894 , 893 , 891 , 887 , 879 , 863 , 766 , 765 , 763 , 759 , 751 , 831 , 735 , 510 , 509 , 507 , 1016 , 1012 , 1010 , 1009 , 1004 , 1002 , 1001 , 998 , 997 , 995 , 988 , 986 , 985 , 503 , 982 , 981 , 979 , 974 , 973 , 971 , 967 , 956 , 954 , 953 , 950 , 949 , 947 , 942 , 941 , 939 , 935 , 703 , 495 , 926 , 925 , 923 , 919 , 892 , 890 , 889 , 886 , 885 , 883 , 878 , 877 , 875 , 871 , 911 , 479 , 862 , 861 , 859 , 855 , 639 , 764 , 762 , 761 , 758 , 757 , 755 , 750 , 749 , 747 , 830 , 829 , 847 , 447 , 827 , 743 , 823 , 734 , 733 , 731 , 727 , 508 , 506 , 505 , 1008 , 1000 , 996 , 994 , 993 , 984 , 502 , 980 , 978 , 501 , 977 , 972 , 970 , 969 , 966 , 965 , 499 , 815 , 952 , 948 , 946 , 945 , 940 , 938 , 937 , 963 , 934 , 702 , 933 , 701 , 494 , 493 , 931 , 699 , 924 , 491 , 922 , 921 , 918 , 917 , 383 , 719 , 888 , 884 , 882 , 915 , 881 , 876 , 874 , 873 , 695 , 870 , 869 , 487 , 910 , 478 , 909 , 477 , 799 , 860 , 858 , 867 , 857 , 907 , 475 , 854 , 638 , 853 , 637 , 760 , 756 , 754 , 753 , 748 , 851 , 746 , 687 , 635 , 745 , 903 , 471 , 828 , 846 , 446 , 826 , 742 , 845 , 255 , 445 , 825 , 741 , 822 , 732 , 730 , 821 , 843 , 443 , 729 , 631 , 739 , 726 , 504 , 992 , 500 , 976 , 819 , 463 , 968 , 725 , 964 , 498 , 814 , 671 , 944 , 936 , 962 , 497 , 932 , 700 , 813 , 492 , 961 , 930 , 839 , 698 , 439 , 490 , 920 , 916 , 382 , 718 , 623 , 723 , 929 , 697 , 489 , 914 , 880 , 381 , 811 , 872 , 717 , 694 , 868 , 486 , 913 , 908 , 476 , 798 , 866 , 856 , 693 , 906 , 474 , 485 , 797 , 431 , 379 , 852 , 715 , 636 , 865 , 905 , 473 , 752 , 807 , 850 , 686 , 634 , 744 , 607 , 691 , 902 , 470 , 483 , 795 , 844 , 254 , 444 , 824 , 740 , 849 , 685 , 633 , 375 , 711 , 901 , 820 , 469 , 842 , 442 , 728 , 630 , 738 , 253 , 415 , 818 , 791 , 683 , 462 , 724 , 841 , 441 , 670 , 629 , 737 , 899 , 467 , 575 , 496 , 812 , 251 , 960 , 838 , 438 , 367 , 817 , 461 , 622 , 722 , 928 , 669 , 696 , 488 , 380 , 627 , 810 , 716 , 679 , 837 , 437 , 912 , 783 , 621 , 721 , 247 , 459 , 692 , 667 , 484 , 796 , 809 , 430 , 378 , 714 , 351 , 864 , 904 , 835 , 472 , 435 , 806 , 619 , 606 , 690 , 482 , 429 , 377 , 794 , 713 , 455 , 239 , 848 , 663 , 684 , 632 , 374 , 805 , 710 , 900 , 468 , 605 , 689 , 252 , 481 , 319 , 793 , 414 , 615 , 427 , 373 , 709 , 790 , 682 , 840 , 440 , 803 , 628 , 736 , 898 , 655 , 466 , 223 , 574 , 250 , 603 , 413 , 366 , 816 ,/* 460 , 789 , 681 , 668 , 423 , 371 , 707 , 897 , 465 , 626 , 573 , 678 , 836 , 249 , 436 , 365 , 599 , 411 , 782 , 620 , 720 , 191 , 246 , 458 , 787 , 666 , 625 , 808 , 350 , 677 , 571 , 834 , 434 , 781 , 363 , 407 , 245 , 457 , 591 , 618 , 665 , 127 , 428 , 376 , 712 , 349 , 454 , 675 , 238 , 567 ,/* 662 , 833 , 433 , 804 , 779 , 359 , 243 , 617 , 399 , 604 , 688 , 480 , 318 , 792 , 453 , 347 , 237 , 614 , 661 , 426 , 372 , 559 , 708 , 775 , 802 , 654 , 317 , 222 , 451 , 602 , 412 , 613 , 235 , 425 , 343 , 659 , 788 , 680 , 543 , 801 , 422 , 653 , 221 , 370 , 706 , 315 , 896 , 464 , 601, 572 , 611 , 248 , 231 , 335 , 364 , 598 , 410 , 421 , 190 , 369 , 651 , 705 , 219 , 786 , 311 , 624 , 676 , 570 , 597 , 409 , 189 , 419 , 780 , 362 , 647 , 215 , 785 , 406 , 244 , 303 , 456 , 590 , 664 , 569 , 595 , 126 , 187 , 348 , 674 , 361 , 566 , 405 , 207 , 832 , 589 , 432 , 287 , 778 , 358 , 242 , 616 , 125 , 398 , 183 , 673 , 565 , 403 , 452 , 587 , 777 , 346 , 236 , 660 , 357 , 241 , 123 , 558 , 397 , 175 , 774 , 563 , 583 , 345 , 316 , 355 , 119 , 450 , 395 , 557 , 612 , 159 , 234 , 773 , 424 , 342 , 658 , 542 , 800 , 111 , 449 , 391 , 555 , 652 , 771 , 233 , 220 , 341 , 657 , 314 , 600 , 610 , 541 , 95 , 230 , 551 , 334 , 339 , 313 , 420 , 609 , 368 , 650 , 704 , 63 , 539 , 218 , 229 , 333 , 310 , 649 , 535 , 217 , 596 , 227 , 408 , 331 , 309 , 188 , 418 , 646 , 214 , 784 , 527 , 327 , 302 , 307 , 417 , 568 , 645 , 213 , 594 , 301 , 186 , 360 , 643 , 404 , 206 , 211 , 593 , 588 , 286 , 299 , 185 , 124 , 205 , 285 , 295 , 182 , 672 , 564 , 402 , 203 , 586 , 776 , 283 , 181 , 356 , 240 , 401 , 199 , 122 , 585 , 396 , 279 , 174 , 179 , 562 , 121 , 582 , 344 , 271 , 173 , 354 , 561 , 118 , 581 , 394 , 556 , 158 , 171 , 772 , 353 , 117 , 579 , 393 , 157 , 167 , 110 , 115 , 448 , 390 , 554 , 155 , 770 , 232 , 340 , 656 , 109 , 389 , 553 , 151 , 769 , 540 , 94 , 107 , 387 , 550 , 143 , 338 , 312 , 93 , 103 , 549 , 337 , 608 , 62 , 538 , 91 , 228 , 547 , 332 , 61 , 537 , 87 , 648 , 59 , 79 , 534 , 216 , 226 , 330 , 308 , 55 , 533 , 225 , 329 , 47 , 526 , 531 , 326 , 306 , 31 , 525 , 416 , 325 , 305 , 644 , 212 , 523 , 323 , 300 , 519 , 642 , 210 , 592 , 298 , 641 , 184 , 209 , 297 , 204 , 284 , 294 , 293 , 202 , 282 , 291 , 180 , 201 , 281 , 400 , 198 , 584 , 278 , 178 , 197 , 277 , 177 , 195 , 120 , 270 , 275 , 172 , 269 , 560 , 580 , 267 , 170 , 352 , 263 , 169 , 116 , 578 , 392 , 156 , 166 , 577 , 165 , 114 , 154 , 163 , 113 , 153 , 108 , 388 , 552 , 150 , 768 , 149 , 106 , 386 , 142 , 147 , 105 , 385 , 141 , 92 , 102 , 548 , 139 , 101 , 336 , 135 , 90 , 99 , 546 , 89 , 545 , 60 , 536 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 532 , 224 , 328 , 53 , 71 , 46 , 51 , 530 , 45 , 529 , 30 , 43 , 524 , 29 , 39 , 324 , 304 , 27 , 522 , 23 , 322 , 521 , 15 , 321 , 518 , 517 , 515 , 640 , 208 , 296 , 292 , 290 , 200 , 289 , 280 , 196 , 276 , 176 , 194 , 193 , 274 , 273 , 268 , 266 , 265 , 262 , 168 , 261 , 259 , 576 , 164 , 162 , 112 , 161 , 152 , 148 , 146 , 104 , 145 , 384 , 140 , 138 , 100 , 137 , 134 , 98 , 133 , 97 , 131 , 88 , 544 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 528 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 520 , 21 , 14 , 19 , 320 , 13 , 11 , 516 , 7 , 514 , 513 , 288 , 192 , 272 , 264 , 260 , 258 , 257 , 160 , 144 , 136 , 132 , 96 , 130 , 129 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 512 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@-1.5dB
	//1023 , 1022 , 1021 , 1019 , 1015 , 1007 , 991 , 959 , 895 , 767 , 511 , 1020 , 1018 , 1017 , 1014 , 1013 , 1011 , 1006 , 1005 , 1003 , 999 , 990 , 989 , 987 , 983 , 975 , 958 , 957 , 955 , 951 , 943 , 927 , 894 , 893 , 891 , 887 , 879 , 863 , 831 , 766 , 765 , 763 , 759 , 751 , 735 , 703 , 510 , 509 , 507 , 503 , 495 , 479 , 1016 , 1012 , 1010 , 1009 , 1004 , 1002 , 1001 , 998 , 997 , 995 , 988 , 986 , 985 , 982 , 981 , 979 , 974 , 973 , 971 , 967 , 956 , 954 , 953 , 950 , 949 , 947 , 942 , 941 , 939 , 935 , 639 , 926 , 925 , 923 , 919 , 911 , 892 , 890 , 889 , 886 , 885 , 883 , 878 , 877 , 875 , 871 , 447 , 862 , 861 , 859 , 855 , 847 , 830 , 829 , 827 , 764 , 762 , 761 , 758 , 757 , 755 , 750 , 749 , 747 , 823 , 743 , 734 , 733 , 731 , 727 , 815 , 383 , 719 , 702 , 701 , 699 , 508 , 506 , 505 , 502 , 501 , 499 , 494 , 493 , 695 , 491 , 799 , 487 , 478 , 477 , 475 , 687 , 471 , 1008 , 1000 , 996 , 994 , 993 , 984 , 980 , 978 , 977 , 972 , 970 , 969 , 966 , 965 , 952 , 948 , 946 , 945 , 963 , 940 , 938 , 937 , 934 , 933 , 638 , 637 , 931 , 635 , 255 , 924 , 922 , 921 , 918 , 917 , 915 , 910 , 909 , 888 , 884 , 882 , 881 , 631 , 876 , 874 , 873 , 907 , 870 , 869 , 867 , 463 , 446 , 860 , 445 , 858 , 857 , 854 , 853 , 443 , 671 , 851 , 903 , 846 , 845 , 843 , 439 , 623 , 828 , 826 , 825 , 760 , 756 , 754 , 753 , 748 , 746 , 745 , 822 , 821 , 742 , 741 , 819 , 839 , 739 , 732 , 730 , 729 , 726 , 725 , 814 , 813 , 723 , 382 , 381 , 811 , 431 , 379 , 718 , 717 , 700 , 698 , 715 , 607 , 697 , 504 , 500 , 498 , 497 , 492 , 694 , 490 , 807 , 693 , 798 , 489 , 797 , 486 , 375 , 485 , 691 , 795 , 476 , 474 , 483 , 415 , 473 , 711 , 686 , 685 , 470 , 469 , 992 , 976 , 968 , 964 , 944 , 962 , 936 , 932 , 636 , 930 , 961 , 634 , 254 , 920 , 683 , 929 , 916 , 633 , 253 , 791 , 914 , 367 , 467 , 913 , 908 , 880 , 630 , 872 , 575 , 906 , 868 , 866 , 629 , 462 , 905 , 444 , 856 , 251 , 865 , 852 , 442 , 461 , 670 , 850 , 441 , 902 , 669 , 849 , 844 , 901 , 679 , 842 , 627 , 438 , 841 , 459 , 622 , 437 , 824 , 752 , 744 , 783 , 820 , 667 , 740 , 621 , 818 , 838 , 247 , 738 , 899 , 728 , 817 , 351 , 837 , 724 , 812 , 737 , 722 , 435 , 380 , 810 , 430 , 378 , 619 , 721 , 809 , 716 , 455 , 429 , 835 , 377 , 714 , 663 , 606 , 696 , 496 , 806 , 692 , 488 , 796 , 239 , 713 , 605 , 374 , 484 , 690 , 794 , 805 , 427 , 615 , 482 , 414 , 472 , 710 , 373 , 684 , 689 , 793 , 319 , 468 , 481 , 603 , 413 , 960 , 709 , 682 , 928 , 655 , 632 , 252 , 790 , 803 , 366 , 466 , 912 , 574 , 371 , 423 , 628 , 904 , 250 , 864 , 460 , 681 , 223 , 789 , 440 , 668 , 848 , 365 , 465 , 707 , 411 , 900 , 678 , 573 , 626 , 599 , 840 , 458 , 436 , 249 , 782 , 666 , 620 , 246 , 787 , 898 , 677 , 625 , 363 , 816 , 350 , 836 , 736 , 457 , 571 , 434 , 781 , 407 , 665 , 618 , 720 , 245 , 808 , 897 , 454 , 428 , 191 , 591 , 349 , 834 , 376 , 662 , 675 , 433 , 238 , 617 , 712 , 359 , 604 , 779 , 453 , 567 , 804 , 426 , 243 , 661 , 614 , 833 , 347 , 372 , 399 , 688 , 792 , 318 , 237 , 480 , 602 , 412 , 708 , 425 , 654 , 127 , 802 , 451 , 613 , 775 , 659 , 370 , 422 , 317 , 559 , 680 , 222 , 788/* , 601 , 235 , 343 , 364 , 464 , 653 , 706 , 801 , 410 , 572 , 598 , 611 , 369 , 248 , 421 , 221 , 315 , 786 , 676 , 705 , 624 , 409 , 362 , 543 , 651 , 597 , 231 , 456 , 570 , 335 , 780 , 406 , 664 , 419 , 244 , 896 , 190 , 785 , 590 , 219 , 348 , 311 , 361 , 674 , 432 , 569 , 595 , 647 , 405 , 616 , 358 , 778 ,/* 452 , 566 , 242 , 189 , 660 , 589 , 832 , 346 , 673 , 215 , 398 , 303 , 236 , 357 , 403 , 777 , 565 , 424 , 126 , 241 , 187 , 450 , 612 , 587 , 774 , 345 , 658 , 397 , 316 , 207 , 558 , 355 , 287 , 600 , 234 , 342 , 563 , 652 , 125 , 449 , 800 , 183 , 773 , 583 , 657 , 610 , 395 , 368 , 420 , 557 , 220 , 233 , 314 , 341 , 123 , 704 , 408 , 542 , 771 , 650 , 596 , 609 , 230 , 175 , 334 , 391 , 555 , 418 , 313 , 339 , 119 , 784 , 541 , 218 , 649 , 229 , 310 , 360 , 333 , 159 , 568 , 551 , 594 , 646 , 417 , 404 , 217 , 539 , 309 , 188 , 111 , 227 , 588 , 331 , 672 , 214 , 593 , 645 , 302 , 356 , 402 , 776 , 307 , 535 , 95 , 564 , 327 , 240 , 186 , 213 , 643 , 586 , 301 , 344 , 396 , 401 , 206 , 527 , 63 , 354 , 185 , 286 , 211 , 585 , 562 , 299 , 124 , 448 , 182 , 772 , 205 , 582 , 656 , 394 , 353 , 285 , 556 , 561 , 295 , 232 , 340 , 181 , 203 , 581 , 122 , 393 , 283 , 770 , 608 , 174 , 390 , 554 , 179 , 199 , 579 , 121 , 312 , 279 , 338 , 769 , 173 , 118 , 389 , 553 , 540 , 648 , 228 , 332 , 271 , 337 , 158 , 550 , 171 , 416 , 117 , 387 , 216 , 157 , 538 , 167 , 308 , 549 , 110 , 226 , 115 , 330 , 592 , 644 , 155 , 537 , 547 , 109 , 225 , 329 , 306 , 151 , 534 , 94 , 107 , 326 , 212 , 642 , 143 , 305 , 300 , 533 , 93 , 103 , 325 , 400 , 641 , 526 , 531 , 62 , 91 , 323 , 184 , 210 , 584 , 298 , 525 , 61 , 87 , 209 , 204 , 297 , 523 , 59 , 79 , 352 , 284 , 560 , 294 , 519 , 55 , 180 , 202 , 580 , 293 , 47 , 392 , 282 , 201 , 291 , 31 , 281 , 178 , 198 , 578 , 120 , 278 , 177 , 768 , 197 , 577 , 172 , 388 , 277 , 552 , 195 , 270 , 275 , 336 , 170 , 269 , 116 , 386 , 169 , 267 , 385 , 156 , 166 , 548 , 263 , 114 , 165 , 113 , 154 , 536 , 163 , 546 , 108 , 224 , 153 , 328 , 545 , 150 , 106 , 149 , 105 , 142 , 304 , 147 , 532 , 92 , 102 , 324 , 141 , 101 , 640 , 139 , 530 , 90 , 99 , 322 , 135 , 529 , 89 , 321 , 524 , 60 , 86 , 208 , 85 , 296 , 522 , 58 , 78 , 83 , 521 , 57 , 77 , 518 , 54 , 75 , 517 , 53 , 71 , 292 , 515 , 46 , 51 , 45 , 200 , 290 , 30 , 43 , 280 , 289 , 29 , 39 , 27 , 176 , 23 , 196 , 576 , 15 , 276 , 194 , 193 , 274 , 273 , 268 , 168 , 266 , 384 , 265 , 262 , 261 , 164 , 259 , 112 , 162 , 161 , 152 , 544 , 148 , 104 , 146 , 145 , 140 , 100 , 138 , 98 , 137 , 97 , 134 , 528 , 88 , 133 , 320 , 131 , 84 , 82 , 81 , 520 , 56 , 76 , 74 , 73 , 516 , 52 , 70 , 69 , 514 , 50 , 67 , 513 , 49 , 44 , 42 , 41 , 288 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 192 , 272 , 264 , 260 , 258 , 257 , 160 , 144 , 136 , 96 , 132 , 130 , 129 , 80 , 72 , 68 , 66 , 512 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@-3dB
	//1023 , 1022 , 1021 , 1019 , 1015 , 1007 , 991 , 959 , 895 , 767 , 511 , 1020 , 1018 , 1017 , 1014 , 1013 , 1011 , 1006 , 1005 , 1003 , 999 , 990 , 989 , 987 , 983 , 975 , 958 , 957 , 955 , 951 , 943 , 927 , 894 , 893 , 891 , 887 , 879 , 863 , 831 , 766 , 765 , 763 , 759 , 751 , 735 , 703 , 510 , 509 , 507 , 503 , 495 , 479 , 639 , 447 , 1016 , 1012 , 1010 , 1009 , 1004 , 1002 , 1001 , 998 , 997 , 995 , 988 , 986 , 985 , 982 , 981 , 979 , 974 , 973 , 971 , 967 , 956 , 954 , 953 , 950 , 949 , 947 , 942 , 941 , 939 , 935 , 926 , 925 , 923 , 919 , 383 , 911 , 892 , 890 , 889 , 886 , 885 , 883 , 878 , 877 , 875 , 871 , 862 , 861 , 859 , 855 , 847 , 830 , 829 , 827 , 823 , 815 , 764 , 762 , 761 , 758 , 757 , 755 , 750 , 749 , 747 , 743 , 734 , 733 , 731 , 727 , 719 , 799 , 702 , 701 , 699 , 695 , 687 , 255 , 508 , 506 , 505 , 502 , 501 , 499 , 494 , 493 , 491 , 487 , 478 , 477 , 475 , 471 , 671 , 638 , 637 , 635 , 631 , 463 , 623 , 446 , 445 , 443 , 439 , 431 , 607 , 1008 , 1000 , 996 , 994 , 993 , 984 , 980 , 978 , 977 , 972 , 970 , 969 , 966 , 965 , 963 , 952 , 948 , 946 , 945 , 940 , 938 , 937 , 934 , 933 , 931 , 924 , 922 , 921 , 918 , 917 , 915 , 382 , 381 , 379 , 910 , 909 , 907 , 888 , 884 , 882 , 881 , 876 , 874 , 873 , 375 , 870 , 869 , 867 , 860 , 858 , 857 , 854 , 853 , 903 , 851 , 846 , 845 , 843 , 415 , 839 , 828 , 826 , 825 , 822 , 821 , 819 , 814 , 813 , 811 , 367 , 807 , 760 , 756 , 754 , 753 , 748 , 746 , 745 , 742 , 741 , 739 , 575 , 732 , 730 , 729 , 726 , 725 , 723 , 718 , 717 , 715 , 798 , 797 , 795 , 711 , 700 , 698 , 697 , 791 , 694 , 693 , 691 , 686 , 685 , 254 , 683 , 253 , 251 , 351 , 679 , 504 , 500 , 498 , 497 , 492 , 247 , 490 , 489 , 486 , 485 , 483 , 476 , 474 , 473 , 470 , 469 , 467 , 783 , 670 , 669 , 636 , 634 , 633 , 667 , 630 , 629 , 462 , 461 , 627 , 459 , 622 , 621 , 619 , 663 , 444 , 442 , 441 , 239 , 438 , 437 , 435 , 455 , 430 , 429 , 615 , 319 , 427 , 606 , 605 , 603 , 655 , 992 , 976 , 968 , 964 , 962 , 944 , 936 , 932 , 930 , 961 , 929 , 920 , 916 , 914 , 380 , 378 , 913 , 377 , 908 , 906 , 905 , 423 , 880 , 872 , 374 , 868 , 866 , 373 , 865 , 856 , 852 , 902 , 850 , 901 , 849 , 223 , 844 , 842 , 414 , 841 , 413 , 371 , 838 , 899 , 837 , 824 , 820 , 818 , 599 , 817 , 411 , 812 , 810 , 809 , 835 , 366 , 365 , 806 , 805 , 752 , 744 , 740 , 738 , 574 , 728 , 737 , 724 , 722 , 573 , 721 , 363 , 716 , 714 , 796 , 803 , 794 , 713 , 793 , 571 , 710 , 696 , 407 , 790 , 692 , 709 , 690 , 789 , 689 , 684 , 591 , 682 , 252 , 250 , 681 , 249 , 707 , 787 , 350 , 359 , 191 , 678 , 349 , 496 , 246 , 488 , 484 , 677 , 482 , 245 , 472 , 481 , 468 , 567 , 466 , 782 , 668 , 632 , 465 , 781 , 666 , 628 , 347 , 460 , 626 , 665 , 399 , 675 , 458 , 243 , 625 , 457 , 620 , 618 , 662 , 440 , 238 , 436 , 779 , 617 , 661 , 434 , 454 , 237 , 433 , 453 , 428 , 614 , 318 , 426 , 343 , 559 , 613 , 317 , 604 , 425 , 659 , 235 , 602 , 654 , 127 , 960 , 928 , 451 , 912 , 775 , 376 , 904 , 422 , 601 , 653 , 372 , 864 , 900 , 848 , 222 , 840 , 412 , 370 , 611 , 421 , 898 , 836 , 315 , 598 , 816 , 410 , 221 , 808 , 834 , 369 , 364 , 897 , 335 , 804 , 231 , 597 , 409 , 651 , 736 , 572 , 720 , 362 , 833 , 802 , 712 , 543 , 792 , 419 , 570 , 406 , 219 , 708 , 788 , 361 , 688 , 311 , 590 , 801 , 680 , 248 , 595 , 569 , 706 , 786 , 358 , 405 , 190 , 348 , 676 , 244 , 647 , 589 , 480 , 566 , 464 , 780 , 705 , 785 , 346 , 357 , 189 , 664 , 398 , 674 , 215 , 242 , 624 , 456 , 403 , 565 , 778 , 303 , 616 , 660 , 236 , 345 , 587 , 432 , 452 , 397 , 673 , 241 , 355 , 187 , 342 , 558 , 612 , 316 , 777 , 424 , 658 , 234 , 563 , 126 , 450 , 774 , 207 , 600 , 652 , 395 , 341 , 557 , 583 , 610 , 420 , 314 , 657 , 233 , 220 , 287 , 125 , 368 , 449 , 773 , 183 , 896 , 334 , 230 , 596 , 408 , 650 , 832 , 609 , 339 , 542 , 313 , 555 , 418 , 391 , 218 , 360 , 310 , 123 , 800 , 333 , 229 , 771 , 649 , 594 , 568 , 404 , 541 , 417 , 175 , 646 , 588 , 217 , 309 , 704 , 551 , 784 , 356 , 188 , 593 , 214 , 331 , 227 , 402 , 564 , 119 , 645 , 302 , 539 , 344 , 586 , 396 , 672 , 240 , 307 , 159 , 213 , 354 , 186 , 401 , 776 , 327 , 562 , 301 , 643 , 585 , 206 , 111 , 535 , 394 , 340 , 556 , 353 , 185 , 582 , 211 , 656 , 232 , 286 , 124 , 561 , 448 , 772 , 182 , 299 , 205 , 393 , 608 , 338 , 312 , 554 , 527 , 581 , 95 , 390 , 285 , 181 , 122 , 332 , 228 , 770 , 648 , 203 , 295 , 540 , 337 , 553 , 416 , 174 , 579 , 389 , 216 , 283 , 308 , 63 , 179 , 121 , 550 , 769 , 592 , 330 , 199 , 226 , 173 , 118 , 644 , 538 , 387 , 279 , 549 , 306 , 329 , 158 , 212 , 225 , 400 , 117 , 171 , 537 , 326 , 300 , 642 , 547 , 584 , 271 , 110 , 305 , 534 , 157 , 352 , 184 , 115 , 210 , 167 , 325 , 641 , 560 , 109 , 533 , 298 , 155 , 204 , 392 , 209 , 526 , 323 , 580 , 94 , 284 , 107 , 297 , 531 , 180 , 151 , 202 , 294 , 525 , 93 , 336 , 552/* , 578 , 103 , 388 , 143 , 282 , 201 , 293 , 62 , 178 , 120 , 523 , 91 , 768 , 577 , 198 , 281 , 172 , 61 , 291 , 177 , 386 , 519 , 87 , 278 , 548 , 197 , 328 , 224 , 59 , 385 , 116 , 79 , 170 , 277 , 536 , 195 , 546 , 55 , 270 , 304 , 156 , 169 , 275 , 114 , 545 , 166 , 47 , 324 , 269 , 640 , 108 , 532 , 113 , 154 , 31 , 165 , 267 , 208 , 322 , 153 , 163 , 263 , 106 , 296 , 530 , 150 , 321 , 105 , 524 , 529 , 92 , 149 , 102 , 142 , 147 , 200 , 292 , 101 , 522 , 90 , 141 , 576 , 99 , 521 , 89 , 280 , 139 , 60 , 290 , 176 , 518 , 86 , 135 , 196 , 289 , 517 , 85 , 58 , 384 , 515 , 78 , 83 , 276 , 194 , 57 , 77 , 193 , 54 , 75 , 168 , 274 , 53 , 71 , 273 , 544 , 46 , 51 , 268 , 45 , 112 , 30 , 164 , 43 , 266 , 29 , 39 , 265 , 152 , 27 , 162 , 262 , 23 , 161 , 261 , 320 , 15 , 259 , 104 , 528 , 148 , 146 , 145 , 100 , 140 , 98 , 520 , 88 , 138 , 97 , 137 , 134 , 288 , 133 , 516 , 84 , 131 , 514 , 82 , 513 , 81 , 56 , 76 , 192 , 74 , 73 , 52 , 70 , 272 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 264 , 37 , 26 , 35 , 25 , 22 , 160 , 260 , 21 , 14 , 19 , 258 , 13 , 257 , 11 , 7 , 144 , 96 , 136 , 132 , 130 , 129 , 512 , 80 , 72 , 68 , 66 , 48 , 65 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 256 , 10 , 9 , 6 , 5 , 3 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@-5dB
	
	//511 , 510 , 509 , 507 , 503 , 495 , 479 , 447 , 383 , 508 , 506 , 505 , 502 , 501 , 499 , 494 , 493 , 491 , 487 , 478 , 477 , 475 , 471 , 255 , 446 , 445 , 443 , 463 , 439 , 382 , 381 , 431 , 379 , 504 , 500 , 498 , 497 , 492 , 490 , 489 , 375 , 486 , 485 , 476 , 415 , 474 , 483 , 473 , 470 , 469 , 254 , 253 , 367 , 467 , 444 , 251 , 442 , 462 , 441 , 461 , 438 , 437 , 459 , 247 , 351 , 435 , 380 , 430 , 378 , 429 , 377 , 455 , 496 , 239 , 488 , 374 , 484 , 427 , 373 , 414 , 482 , 472 , 319 , 468 , 413 , 252 , 481 , 366 , 371 , 466 , 423 , 250 , 223 , 440 , 365 , 460 , 411 , 465 , 436 , 249 , 458 , 246 , 363 , 350 , 434 , 407 , 457 , 245 , 428 , 191 , 349 , 376 , 454 , 433 , 238 , 359 , 426 , 243 , 453 , 347 , 372 , 399 , 318 , 237 , 425 , 412 , 127 , 480 , 370 , 451 , 422 , 317 , 222 , 235 , 343 , 364 , 410 , 464 , 369 , 248 , 421 , 315 , 221 , 409 , 362 , 231 , 335 , 406 , 456 , 419 , 244 , 190 , 348 , 219 , 311 , 361 , 432 , 405 , 358 , 242 , 189 , 452 , 346 , 215 , 398 , 303 , 236 , 357 , 403 , 424 , 126 , 241 , 187 , 345 , 397 , 450 , 316 , 207 , 355 , 287 , 234 , 342 , 125 , 183 , 449 , 395 , 368 , 420 , 233 , 314 , 341 , 220 , 123 , 408 , 230 , 175 , 334 , 391 , 418 , 313 , 339 , 119 , 218 , 229 , 310 , 360 , 333 , 159 , 417 , 404 , 217 , 309 , 188 , 111 , 227 , 331 , 214 , 302 , 356 , 402 , 307 , 95 , 327 , 240 , 186 , 213 , 301 , 344 , 396 , 401 , 206 , 63 , 354 , 185 , 286 , 211 , 299 , 124 , 182 , 205 , 448 , 394 , 353 , 285 , 295 , 232 , 340 , 181 , 203 , 122 , 393 , 283 , 174 , 390 , 179 , 199 , 121 , 312 , 279 , 338 , 173 , 118 , 389 , 228 , 332 , 271 , 337 , 158 , 171 , 416 , 387 , 117 , 216 , 157 , 308 , 167 , 110 , 226 , 115 , 330 , 155 , 109 , 225 , 329 , 306 , 151 , 94 , 107 , 326 , 212 , 305 , 300 , 143 , 93 , 103 , 325 , 400 , 62 , 91 , 323 , 184 , 210 , 298 , 61 , 87 , 209 , 204 , 297 , 59 , 79 , 352 , 284 , 294 , 55 , 180 , 202 , 293 , 392 , 47 , 282 , 201 , 291 , 31 , 281 , 178 , 198 , 120 , 278 , 177 , 197 , 172 , 388 , 277 , 195 , 270 , 275 , 336 , 170 , 386 , 116 , 269 , 169 , 267 , 385 , 156 , 166 , 263 , 114 , 165 , 113 , 154 , 163 , 108 , 224 , 328 , 153 , 150 , 106 , 149 , 105 , 304 , 142 , 147 , 92 , 102 , 324 , 141 , 101 , 139 , 90 , 99 , 322 , 135 , 89 , 321 , 60 , 86 , 208 , 85 , 296 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 71 , 292 , 46 , 51 , 45 , 200 , 290 , 30 , 43 , 280 , 289 , 29 , 39 , 27 , 176 , 23 , 196 , 15 , 276 , 194 , 193 , 274 , 273 , 268 , 168 , 266 , 384 , 265 , 262 , 261 , 164 , 259 , 112 , 162 , 161 , 152 , 148 , 104 , 146 , 145 , 140 , 100 , 138 , 98 ,/* 137 , 97 , 134 , 88 , 133 , 320 , 131 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 288 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 192 , 272 , 264 , 260 , 258 , 257 , 160 , 144 , 136 , 96 , 132 , 130 , 129 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@0dB
	//511 , 510 , 509 , 507 , 503 , 495 , 479 , 447 , 383 , 255 , 508 , 506 , 505 , 502 , 501 , 499 , 494 , 493 , 491 , 487 , 478 , 477 , 475 , 471 , 463 , 446 , 445 , 443 , 439 , 431 , 382 , 381 , 379 , 375 , 415 , 367 , 254 , 253 , 251 , 351 , 504 , 500 , 498 , 247 , 497 , 492 , 490 , 489 , 486 , 485 , 483 , 476 , 474 , 473 , 470 , 469 , 467 , 462 , 461 , 459 , 444 , 442 , 441 , 239 , 438 , 437 , 435 , 455 , 430 , 429 , 319 , 427 , 380 , 378 , 377 , 423 , 374 , 373 , 223 , 414 , 413 , 371 , 411 , 366 , 365 , 363 , 407 , 252 , 250 , 249 , 350 , 359 , 191 , 349 , 246 , 496 , 488 , 484 , 245 , 482 , 472 , 481 , 468 , 466 , 465 , 347 , 460 , 399 , 458 , 243 , 457 , 440 , 238 , 436 , 434 , 454 , 237 , 433 , 453 , 428 , 318 , 426 , 343 , 317 , 425 , 235 , 127 , 451 , 376 , 422 , 372 , 222 , 412 , 370 , 421 , 315 , 410 , 221 , 369 , 364 , 335 , 231 , 409 , 362 , 419 , 406 , 219 , 361 , 311 , 248 , 358 , 405 , 190 , 348 , 244 , 480 , 464 , 346 , 357 , 189 , 398 , 215 , 242 , 456 , 403 , 303 , 236 , 345 , 432 , 452 , 397 , 241 , 355 , 187 , 342 , 316 , 424 , 234 , 126 , 450 , 207 , 395 , 341 , 420 /*, 314 , 233 , 220 , 287 , 125 , 368 , 449 , 183 , 334 , 230 , 408 , 313 , 339 , 418 , 391 , 218 , 360 , 310 , 123 , 333 , 229 , 404 , 417 , 175 , 217 , 309 , 356 , 188 ,/* 214 , 331 , 227 , 402 , 119 , 302 , 344 , 396 , 240 , 307 , 159 , 213 , 354 , 186 , 401 , 327 , 301 , 206 , 111 , 394 , 340 , 353 , 185 , 211 , 232 , 286 , 124 , 448 , 182 , 299 , 205 , 393 , 312 , 338 , 95 , 390 , 285 , 181 , 122 , 332 , 228 , 203 , 295 , 337 , 416 , 174 , 216 , 389 , 283 , 308 , 179 , 63 , 121 , 330 , 199 , 226 , 173 , 118 , 387 , 279 , 306 , 329 , 158 , 212 , 225 , 400 , 117 , 171 , 326 , 300 , 271 , 305 , 110 , 157 , 352 , 184 , 115 , 210 , 167 , 325 , 109 , 298 , 155 , 204 , 392 , 209 , 323 , 94 , 284 , 107 , 297 , 180 , 151 , 202 , 294 , 93 , 336 , 103 , 388 , 143 , 282 , 201 , 293 , 178 , 62 , 120 , 91 , 198 , 281 , 172 , 177 , 61 , 291 , 386 , 87 , 278 , 197 , 328 , 224 , 59 , 385 , 116 , 170 , 79 , 277 , 195 , 55 , 270 , 304 , 156 , 169 , 275 , 114 , 166 , 47 , 324 , 269 , 108 , 113 , 154 , 165 , 31 , 267 , 208 , 322 , 153 , 163 , 263 , 106 , 296 , 150 , 321 , 105 , 92 , 149 , 102 , 142 , 147 , 200 , 292 , 101 , 90 , 141 , 99 , 89 , 280 , 139 , 176 , 60 , 290 , 86 , 135 , 196 , 289 , 85 , 58 , 384 , 78 , 83 , 276 , 194 , 57 , 77 , 193 , 54 , 75 , 168 , 274 , 53 , 71 , 273 , 46 , 51 , 268 , 45 , 112 , 164 , 30 , 43 , 266 , 29 , 39 , 265 , 152 , 27 , 162 , 262 , 23 , 161 , 261 , 320 , 15 , 259 , 104 , 148 , 146 , 145 , 100 , 140 , 98 , 88 , 138 , 97 , 137 , 134 , 288 , 133 , 84 , 131 , 82 , 81 , 56 , 76 , 192 , 74 , 73 , 52 , 70 , 272 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 264 , 37 , 26 , 35 , 25 , 22 , 160 , 260 , 21 , 14 , 19 , 258 , 13 , 257 , 11 , 7 , 144 , 96 , 136 , 132 , 130 , 129 , 80 , 72 , 68 , 66 , 48 , 65 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 256 , 10 , 9 , 6 , 5 , 3 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@-2dB
	//511 , 510 , 509 , 507 , 503 , 495 , 479 , 447 , 508 , 506 , 505 , 502 , 501 , 499 , 494 , 493 , 491 , 383 , 478 , 477 , 487 , 475 , 446 , 255 , 471 , 445 , 504 , 500 , 498 , 443 , 492 , 497 , 463 , 490 , 382 , 489 , 439 , 476 , 381 , 486 , 474 , 485 , 379 , 473 , 431 , 254 , 470 , 483 , 444 , 375 , 253 , 469 , 442 , 415 , 496 , 462 , 251 , 467 , 441 , 367 , 488 , 438 , 461 , 380 , 247 , 437 , 484 , 351 , 459 , 378 , 472 , 430 , 435 , 239 , 482 , 377 , 455 , 429 , 374 , 319 , 252 , 468 , 481 , 414 , 223 , 373 , 427 , 250 , 466 , 440 , 366 , 413 , 371 , 460 , 423 , 191 , 249 , 465 , 365 , 411 , 246 , 436 , 350 , 458 , 127 , 363 , 245 , 407 , 349 , 434 , 457 , 238 , 376 , 359 , 243 , 454 , 428 , 399 , 318 , 433 , 347 , 237 , 480 , 453 , 222 , 317 , 372 , 426 , 235 , 343 , 451 , 221 , 315 , 412 , 425 , 231 , 335 , 370 , 422 , 190 , 248 , 219 , 311 , 464 , 364 , 410 , 369 , 421 , 189 , 215 , 303 , 409 , 126 , 419 , 187 , 362 , 207 , 244 , 287 , 406 , 125 , 348 , 183 , 361 , 456 , 405 , 358 , 242 , 123 , 175 , 398 , 403 , 432 , 346 , 236 , 357 , 241 , 119 , 159 , 397 , 452 , 345 , 355 , 111, 316 , 395 , 234 , 342 , 95 , 391 , 450 , 233 , 341 , 220 , 63 , 314 , 424 , 449 , 230 , 334 , 339 , 313 , 229 , 333 , 218 , 310 , 227 , 331 , 217 , 368 , 309 , 420 , 188 , 327 , 214 , 302 , 307 , 408 , 213 , 301 , 418 , 186 , 206 , 211 , 286 , 299 , 417 , 185 , 205 , 285 , 295 , 124 , 182 , 360 , 203 , 283 , 404 , 181 , 199 , 279 , 122 , 174 , 179 , 271 , 402 , 121 , 173 , 356 , 401 , 240 , 118 , 158 , 171 , 396 , 117 , 157 , 167 , 344 , 354 , 110 , 115 , 155 , 394 , 353 , 109 , 151 , 393 , 94 , 107 , 143 , 390 , 93 , 103 , 232 , 340 , 389 , 62 , 91 , 387 , 448 , 61 , 87 , 338 , 59 , 79 , 312 , 337 , 55 , 228 , 332 , 47 , 31 , 226 , 330 , 216 , 308 , 225 , 329 , 326 , 306 , 325 , 305 , 323 , 212 , 300 , 210 , 298 , 416 , 209 , 184 , 297 , 204 , 284 , 294 , 293 , 202 , 282 , 291 , 201 , 180 , 281 , 198 , 278 , 197 , 178 , 277 , 195 , 177 , 270 , 275 , 120 , 172 , 269 , 267 , 400 , 170 , 263 , 169 , 116 , 156 , 166 , 165 , 114 , 154 , 163 , 113 , 153 , 352 , 108 , 150 , 392 , 149 , 106 , 142 , 147 , 105 , 141 , 92 , 102 , 139 , 101 , 135 , 388 , 90 , 99 , 89 , 386 , 60 , 86 , 385 , 85 , 58 , 78 , 83 , 336 , 57 , 77 , 54 , 75 , 53 , 71 , 46 , 51 , 45 , 30 , 43 , 29 , 39 , 27 , 224 , 23 , 328 , 15 , 324 , 304 , 322 , 321 , 208 , 296 , 292 , 290 , 289 , 200 , 280 , 196 , 276 , 194 , 176 , 274 , 193 , 273 , 268 , 266 , 265 , 262 , 261 , 168 , 259 , 164 , 162 , 161 , 112 , 152 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 ,/* 133 , 98 , 131 , 97 , 88 , 384 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 320 , 288 , 192 , 272 , 264 , 260 , 258 , 257 , 160 , 144 , 136 , 132 , 130 , 129 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@2dB
	//511 , 510 , 509 , 507 , 503 , 495 , 508 , 506 , 479 , 505 , 502 , 447 , 501 , 494 , 383 , 499 , 493 , 478 , 255 , 491 , 504 , 477 , 487 , 446 , 475 , 500 , 445 , 471 , 382 , 498 , 443 , 463 , 492 , 381 , 497 , 439 , 254 , 379 , 431 , 490 , 253 , 375 , 476 , 415 , 489 , 486 , 251 , 367 , 474 , 485 , 247 , 351 , 473 , 483 , 444 , 239 , 319 , 470 , 223 , 469 , 442 , 191 , 462 , 467 , 441 , 127 , 380 , 496 , 461 , 438 , 459 , 437 , 378 , 455 , 430 , 435 , 377 , 429 , 252 , 374 , 414 , 488 , 427 , 373 , 413 , 423 , 250 , 366 , 371 , 411 , 249 , 365 , 407 , 484 , 246 , 350 , 363 , 399 , 245 , 349 , 359 , 472 , 482 , 238 , 243 , 318 , 347 , 481 , 237 , 317 , 343 , 222 , 235 , 315 , 335 , 221 , 231 , 468 , 311 , 190 , 219 , 303 , 189 , 215 , 287 , 466 , 440 , 126 , 187 , 207 , 465 , 125 , 460 , 183 , 123 , 175 , 119 , 159 , 458 , 436 , 111 , 457 , 95 , 454 , 434 , 63 , 453 , 433 , 376 , 451 , 428 , 426 , 425 , 372 , 412 , 422 , 421 , 370 , 410 , 419 , 369 , 248 , 409 , 364 , 406 , 405 , 362 , 398 , 403 , 361 , 244 , 397 , 348 , 358 , 395 , 357 , 242 , 391 , 346 , 355 , 241 , 345 , 480 , 236 ,/* 316 , 342 , 341 , 234 , 314 , 334 , 339 , 233 , 313 , 333 , 220 , 230 , 310 , 331 , 229 , 309 , 327 , 218 , 227 , 302 , 307 , 217 , 301 , 188 , 214 , 286 , 299 , 213 , 285 , 295 , 186 , 206 , 211 , 283 , 464 , 185 , 205 , 279 , 124 , 182 , 203 , 271 , 181 , 199 , 122 , 174 , 179 , 121 , 173 , 118 , 158 , 171 , 117 , 157 , 167 , 110 , 115 , 155 , 456 , 109 , 151 , 94 , 107 , 143 , 93 , 103 , 62 , 91 , 452 , 61 , 87 , 432 , 59 , 79 , 450 , 55 , 449 , 47 , 31 , 424 , 420 , 418 , 368 , 417 , 408 , 404 , 402 , 360 , 401 , 396 , 394 , 356 , 393 , 390 , 354 , 389 , 353 , 240 , 387 , 344 , 340 , 338 , 337 , 232 , 312 , 332 , 330 , 329 , 228 , 308 , 326 , 325 , 226 , 306 , 323 , 225 , 305 , 216 , 300 , 298 , 297 , 212 , 284 , 294 , 293 , 210 , 282 , 291 , 209 , 281 , 184 , 204 , 278 , 277 , 202 , 270 , 275 , 201 , 269 , 180 , 198 , 267 , 197 , 263 , 178 , 195 , 177 , 120 , 172 , 170 , 169 , 116 , 156 , 166 , 165 , 114 , 154 , 163 , 113 , 153 , 108 , 150 , 149 , 106 , 142 , 147 , 105 , 141 , 92 , 102 , 139 , 101 , 135 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 448 , 71 , 46 , 51 , 45 , 30 , 43 , 29 , 39 , 27 , 23 , 15 , 416 , 400 , 392 , 388 , 352 , 386 , 385 , 336 , 328 , 324 , 322 , 321 , 224 , 304 , 296 , 292 , 290 , 289 , 208 , 280 , 276 , 274 , 273 , 200 , 268 , 266 , 265 , 196 , 262 , 261 , 194 , 259 , 193 , 176 , 168 , 164 , 162 , 161 , 112 , 152 , 148 , 146 , 145 , 104 , 140 , 138 , 137 , 100 , 134 , 133 , 98 , 131 , 97 , 88 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 384 , 320 , 288 , 272 , 264 , 260 , 258 , 257 , 192 , 160 , 144 , 136 , 132 , 130 , 129 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 256 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 ,  //@5dB
	
	//255 , 254 , 253 , 251 , 247 , 239 , 252 , 250 , 223 , 249 , 246 , 245 , 191 , 238 , 243 , 237 , 127 , 222 , 235 , 248 , 221 , 231 , 244 , 190 , 219 , 242 , 189 , 215 , 236 , 241 , 126 , 187 , 207 , 234 , 125 , 183 , 233 , 220 , 123 , 230 , 175 , 119 , 218 , 229 , 159 , 217 , 188 , 111 , 227 , 214 , 95 , 240 , 186 , 213 , 206 , 63 ,/* 185 , 211 , 124 , 182 , 205 , 232 , 181 , 203 , 122 , 174 , 179 , 199 , 121 , 173 , 118 , 228 , 158 , 171 , 117 , 216 , 157 , 167 , 110 , 226 , 115 , 155 , 109 , 225 , 151 , 94 , 107 , 212 , 143 , 93 , 103 , 62 , 91 , 184 , 210 , 61 ,/* 87 , 209 , 204 , 59 , 79 , 55 , 180 , 202 , 47 , 201 , 31 , 178 , 198 , 120 , 177 , 197 , 172 , 195 , 170 , 116 , 169 , 156 , 166 , 114 , 165 , 113 , 154 , 163 , 108 , 224 , 153 , 150 , 106 , 149 , 105 , 142 , 147 , 92 , 102 , 141 , 101 , 139 , 90 , 99 , 135 , 89 , 60 , 86 , 208 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 71 , 46 , 51 , 45 , 200 , 30 , 43 , 29 , 39 , 27 , 176 , 23 , 196 , 15 , 194 , 193 , 168 , 164 , 112 , 162 , 161 , 152 , 148 , 104 , 146 , 145 , 140 , 100 , 138 , 98 , 137 , 97 , 134 , 88 , 133 , 131 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 192 , 160 , 144 , 136 , 96 , 132 , 130 , 129 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@3dB
	//255 , 254 , 253 , 251 , 247 , 239 , 223 , 252 , 250 , 249 , 246 , 245 , 191 , 238 , 243 , 237 , 127 , 235 , 222 , 221 , 248 , 231 , 244 , 219 , 190 , 242 , 189 , 236 , 215 , 241 , 187 , 126 , 234 , 207 , 125 , 233 , 220 , 183 , 230 , 123 , 218 , 175 , 229 , 217 , 119 , 188 , 227 , 214 , 159 , 240 , 111 , 186 , 213 , 206 , 185 , 211 , 95 , 124 , 232 , 182 , 205 , 63 , 181 , 122 , 203 , 174 , 228 , 179 , 121 , 199 , 216 , 173 , 118 , 226 , 158 , 117 , 171 , 225 , 110 , 157 , 212 , 115 , 167 , 109 , 155 , 184 , 210 , 94 , 107 , 151 , 204 , 209, 93 , 103 , 143 , 62/* , 91 , 180 , 202 , 61 , 87 , 201 , 178 , 120 , 59 , 79 , 198 , 172 , 177 , 55 , 197 , 47 , 195 , 116 , 170 , 31 , 224 , 169 , 156 , 114 , 166 , 113 , 165 , 108 , 154 , 163 , 153 , 106 ,/* 150 , 105 , 149 , 208 , 92 , 102 , 142 , 147 , 101 , 141 , 90 , 99 , 139 , 89 , 135 , 60 , 86 , 200 , 85 , 58 , 78 , 83 , 57 , 77 , 176 , 54 , 75 , 196 , 53 , 71 , 46 , 51 , 194 , 45 , 193 , 30 , 43 , 29 , 39 , 168 , 27 , 23 , 15 , 112 , 164 , 162 , 152 , 161 , 104 , 148 , 146 , 145 , 100 , 140 , 98 , 138 , 97 , 88 , 137 , 134 , 133 , 131 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 192 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 160 , 144 , 96 , 136 , 132 , 130 , 129 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , */// @2dB
	//255 , 254 , 253 , 251 , 247 , 239 , 223 , 191 , 252 , 250 , 249 , 246 , 245 , 243 , 238 , 237 , 235 , 127 , 222 , 221 , 231 , 219 , 190 , 215 , 189 , 248 , 244 , 187 , 242 , 207 , 241 , 236 , 234 , 126 , 183 , 233 , 125 , 220 , 230 , 218 , 229 , 123 , 175 , 217 , 214 , 188 , 227 , 119 , 213 , 186 , 159 , 206 , 240 , 185 , 211 , 111 , 182 , 205 , 232 , 124 , 181 , 203 , 228 , 122 , 95 , 174 , 216 , 179 , 121 , 226 , 173 , 199 , 118 , 212 , 63 , 225 , 158 , 117 , 171 , 184 , 210 , 157 , 110 , 204 , 115 , 167 , 209 , 109 , 155 , 180 , 202 , 94 , 107 , 151 , 178 , 201 , 120 , 93 , 172 , 198 , 103 , 177 , 143 , 62 , 91 , 224 , 197 , 116 , 170 , 61 , 87 , 195 , 156 , 169 , 59 , 114 , 166 , 79 , 208 , 108 , 154 , 55 , 113 , 165 , 153 , 47 , 163 , 106 , 150 , 200 , 31 , 92 , 105 , 149 , 102 , 176 , 142 , 147 , 90 , 101 , 196 , 141 , 89 , 99 , 60 , 139 , 86 , 194 , 135 , 168 , 85 , 193 , 58 , 78 , 83 , 57 , 77 , 54 , 112 , 164 , 75 , 53 , 152 , 71 , 46 , 51 , 162 , 45 , 161 , 30 , 43 , 104 , 148 , 29 , 39 , 27 , 146 , 23 , 100 , 145 , 15 , 140 , 88 , 98 , 138 , 97/* , 137 , 134 , 84 , 133 , 192 , 131 , 82 , 56 , 81 , 76 , 74 , 52 , 73 , 70 , 50 , 69 , 49 , 67 , 44 , 160 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 144 , 14 , 19 , 13 , 11 , 7 , 96 , 136 , 132 , 130 , 129 , 80 , 72 , 68 , 48 , 66 , 65 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 128 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@0dB

	//127 , 126 , 125 , 123 , 119 , 111 , 124 , 122 , 95 , 121 , 118 , 117 , 63 , 110 , 115 , 109 , 107 , 94 , 120 , 93 , 103 , 116 , 62 , 91 , 61 , 114 , 87 , 108 , 113 , 59 , 79 , 106 , 55 , 105 , 92 , 102 , 47 , 90 , 101 , 31 , 89 , 60 , 99 , 86 , 112 , 58  , 85 , 78/* , 57 , 83 , 54  , 77 , 104 , 53 , 75 , 46 , 51 , 71 , 45 , 100 , 30 , 43 , 88 , 29 , 98 , 39 , 27 , 97 , 23 , 84 , 15 , 56 , 82 , 81 , 76 , 52 , 74 , 73 , 50 , 70 , 49 , 69 , 44 , 67 , 42 , 41 , 28 , 38 , 37 , 26 , 96 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 80 , 72 , 48 , 68 , 66 , 65 , 40 , 36 /*, 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0*///@1.5dB
	//127 , 126 , 125 , 123 , 119 , 111 , 124 , 122 , 95 , 121 , 118 , 63 , 117 , 110 , 115 , 109 , 94 , 107 , 120 , 93 , 103 , 62 , 116 , 91 , 61 , 87 , 114 , 108 , 59 , 113 , 79 , 55 , 106 , 92 , 105 , 47 , 102 , 31 , 90 , 101 , 89 , 60 , 99 , 86 , 58 , 85 , 112 , 78/* , 57 , 83 , 54 , 77 , 53 , 75 , 104 , 46 , 51 , 71 , 45 , 30 , 43 , 100 , 29 , 39 , 88 , 27 , 98 , 23 , 97 , 15 , 84 , 56 , 82 , 81 , 76 , 52 , 74 , 73 , 50 , 70 , 49 , 69 , 44 , 67 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 96 , 21 , 14 , 19 , 13 , 11 , 7 , 80 , 72 , 48 , 68 , 66 , 65 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@2dB
	//127 , 126 , 125 , 123 , 119 , 111 , 124 , 122 , 95 , 121 , 118 , 63 , 117 , 110 , 115 , 109 , 94 , 107 , 120 , 93 , 103 , 62 , 91 , 116 , 61 , 87 , 59 , 114 , 79 , 108 , 55 , 113 , 47 , 106 , 31 , 92 , 105 , 102 , 90 , 101 , 89 , 99 , 60 , 86 , 85 , 58 , 78 , 83/* , 57 , 77 , 54 , 112 , 75 , 53 , 71 , 46 , 51 , 45 , 30 , 43 , 104 , 29 , 39 , 27 , 23 , 100 , 15 , 88 , 98 , 97 , 84 , 82 , 56 , 81 , 76 , 74 , 52 , 73 , 70 , 50 , 69 , 49 , 67 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 96 , 80 , 72 , 68 , 48 , 66 , 65 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@3dB
	//127 , 126 , 125 , 123 , 119 , 111 , 124 , 95 , 122 , 63 , 121 , 118 , 117 , 110 , 115 , 109 , 94 , 107 , 93 , 103 , 62 , 120 , 91 , 61 , 87 , 59 , 116 , 79 , 55 , 114 , 47 , 108 , 113 , 31 , 106 , 105 , 92 , 102 , 101 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 /*, 57 , 77 , 54 , 75 , 53 , 71 , 46 , 51 , 45 , 112 , 30 , 43 , 29 , 39 , 27 , 23 , 15 , 104 , 100 , 98 , 88 , 97 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 96 , 80 , 72 , 68 , 66 , 48 , 65 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@4dB
	//127 , 126 , 125 , 123 , 119 , 111 , 95 , 124 , 63 , 122 , 121 , 118 , 117 , 110 , 115 , 109 , 94 , 107 , 93 , 103 , 62 , 91 , 61 , 87 , 120 , 59 , 79 , 55 , /*47 , 116 , 31 , 114 , 113/* , 108 , 106 , 105 , 92 , 102 , 101 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 71 , 46 , 51 , 45 , 30 , 43 , 29 , 39 , 27 , 23 , 15 , 112 , 104 , 100 , 98 , 97 , 88 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 96 , 80 , 72 , 68 , 66 , 65 , 48 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@5dB


	//63 , 62 , 61 , 59 , 55 , 47 , 60 , 31 , 58 , 57 , 54 , 53 , 46 , 51 , 45 , 30 , 43 , 29 , 56 , 39 , 27 , 52 , 23 , 15/* , 50 , 44 , 49 , 42 , 28 , 41 , 38 , 26 , 37 , 25 , 35 , 22 , 21 , 14 , 19 , 13 , 11 , 48 , 7 , 40 , 36 , 24 , 34 , 33 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , //@2dB
	//63 , 62 , 61 , 59 , 55 , 47 , 60 , 31 , 58 , 57 , 54 , 53 , 46 , 51 , 45 , 30 , 43 , 29 , 39 , 56 , 27 , 23 , 52 , 15 , 50 , 44 , 49 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 48 , 40 , 36 , 34 , 24 , 33 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@2.5dB
	63 , 62 , 61 , 59 , 55 , 47 , 60 , 31 , 58 , 57 , 54 , 53 , 46 , 51/* , 45 , 30 , 43 , 29 , 39 , 27 , 56 , 23 , 15 , 52 , 50 , 44 , 49 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 48 , 40 , 36 , 34 , 24 , 33 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , *///@3dB


	//127 , 126 , 125 , 123 , 119 , 111 , 124 , 95 , 122 , 63 , 121 , 118 , 117 , 110 , 115 , 109 , 94 , 107 , 93 , 103 , 62 , 120 , 91 , 61 , 87 , 59 , 116 , 79 , 55 , 114 , 47 , 108 , 113 , 31 , 106 , 105 , 92 , 102 , 101 , 90 , 99 , 89 , 60 , 86 , 85 , 58 , 78 , 83 , 57 , 77 , 54 , 75 , 53 , 71 , 46 , 51 , 45 , 112 , 30 , 43 , 29 , 39 , 27 , 23 , 15 , 104 , 100 , 98 , 88 , 97 , 84 , 82 , 81 , 56 , 76 , 74 , 73 , 52 , 70 , 69 , 50 , 67 , 49 , 44 , 42 , 41 , 28 , 38 , 37 , 26 , 35 , 25 , 22 , 21 , 14 , 19 , 13 , 11 , 7 , 96 , 80 , 72 , 68 , 66 , 48 , 65 , 40 , 36 , 34 , 33 , 24 , 20 , 18 , 17 , 12 , 10 , 9 , 6 , 5 , 3 , 64 , 32 , 16 , 8 , 4 , 2 , 1 , 0 , //@4dB
	//kirill Positions:
	//0 , 1 , 2 , 3 , 
	//0 , 1 , 2 , 4 , 3 , 5 , 6 , 7 , 
	//0 , 1 , 2 , 4 , 8 , 3 , 5 , 6 , 9 , 10 , 12 , 7 , 11 , 13 , 14 , 15 ,
	//0 , 1 , 2 , 4 , 8 , 16 , 3 , 5 , 6 , 9 , 10 , 17 , 12 , 18 , 20 , 7 , 24 , 11 , 13 , 19 , 14 , 21 , 22 , 25 , 26 , 28 , 15 , 23 , 27 , 29 , 30 , 31 ,
	///*0 , 1 , 2 , 4 , 8 , 16 , 3 , 32 , 5 , 6 , 9 , 10 , 17 , 12 , 18 , 33 , 20 , 34 , 7 , 24 , 36 , 11 , 40 , 13 , 19 , 14 , 48 , 21 , 35 , 22 , 25 , 37 , 26 , 38 , 41 , 28 , 42 , 15 , 49 , 44 ,*/ 50 , 23 , 52 , 27 , 39 , 56 , 29 , 43 , 30 , 45 , 51 , 46 , 53 , 54 , 57 , 58 , 31 , 60 , 47 , 55 , 59 , 61 , 62 , 63 , 
	//0 , 1 , 2 , 4 , 8 , 16 , 3 , 32 , 5 , 6 , 9 , 64 , 10 , 17 , 12 , 18 , 33 , 20 , 34 , 7 , 24 , 36 , 65 , 11 , 66 , 40 , 13 , 19 , 68 , 14 , 48 , 21 , 72 , 35 , 22 , 25 , 37 , 80 , 26 , 38 , 67 , 41 , 28 , 96 , 69 , 42 , 15 , 49 , 70 , 44 , 73 , 50 , 23 , 74 , 52 , 81 , 27 , 76 , 39 , 82 , 56 , 29 , 97 , 84 , 43 , 30 , 98 , 71 , 45 , 88 , 51 , 100 , 46 , 75 , 53 , 104 , 77 , 54 , 83 , 57 ,*/ 78 , 112 , 85 , 58 , 31 , 99 , 86 , 60 , 89 , 101 , 47 , 90 , 102 , 105 , 92 , 55 , 106 , 79 , 113 , 59 , 108 , 114 , 87 , 61 , 116 , 62 , 91 , 103 , 120 , 93 , 107 , 94 , 109 , 115 , 110 , 117 , 63 , 118 , 121 , 122 , 95 , 124 , 111 , 119 , 123 , 125 , 126 , 127 ,
	//0 , 1 , 2 , 4 , 8 , 16 , 3 , 32 , 5 , 6 , 9 , 64 , 10 , 17 , 12 , 18 , 128 , 33 , 20 , 34 , 7 , 24 , 36 , 65 , 11 , 66 , 40 , 13 , 19 , 68 , 14 , 129 , 48 , 21 , 72 , 130 , 35 , 22 , 25 , 132 , 37 , 80 , 26 , 38 , 67 , 136 , 41 , 28 , 96 , 69 , 42 , 15 , 144 , 49 , 70 , 44 , 73 , 131 , 50 , 23 , 74 , 160 , 133 , 52 , 81 , 27 , 76 , 134 , 39 , 82 , 137 , 56 , 29 , 192 , 97 , 138 , 84 , 43 , 30 , 145 , 98 , 71 , 140 , 45 , 88 , 146 , 51 , 100 , 46 , 75 , 161 , 148 , 53 , 104 , 77 , 162 , 135 , 54 , 83 , 152 , 57 , 78 , 164 , 193 , 112 , 139 , 85 , 58 , 31 , 194 , 99 , 168 , 86 , 141 , 60 , 89 , 147 , 196 , 101 , 142 , 47 , 90 , 176 , 149 , 102 , 200 , 105 , 92 , 163 , 150 , 55 , 153 , 106 , 79 , 165 , 208 , 113 , 154 , 59 , 108 , 166 , 195 , 114 , 169 , 87 , 156 , 61 , 224 , 197 , 170 , 116 , 143 , 62 , 91 , 177 , 198 , 103 , 172 , 201 , 120 , 93 , 178 , 151 , 202 , 107 , 94 , 180 , 209 , 155 , 204 , 109 , 167 , 210 , 115 , 184 , 157 , 110 , 225 , 212 , 171 , 117 , 158 , 63 , 226 , 199 , 118 , 173 , 216 , 121 , 179 , 228 , 174 , 203 , 122 , 95 , 181 , 232 , 205 , 124 , 182 , 211 , 185 , 206 , 111 , 240 , 213 , 186 , 159 , 227 , 214 , 119 , 188 , 217 , 229 , 175 , 218 , 123 , 230 , 233 , 220 , 125 , 183 , 234 , 207 , 126 , 241 , 187 , 236 , 242 , 215 , 189 , 244 , 190 , 219 , 231 , 248 , 221 , 235 , 222 , 127 , 237 , 243 , 238 , 245 , 191 , 246 , 249 , 250 , 223 , 252 , 239 , 247 , 251 , 253 , 254 , 255 , 
	//0 , 1 , 2 , 4 , 8 , 16 , 3 , 32 , 5 , 6 , 9 , 64 , 10 , 17 , 12 , 18 , 128 , 33 , 20 , 34 , 7 , 24 , 36 , 65 , 11 , 256 , 66 , 40 , 13 , 19 , 68 , 14 , 129 , 48 , 21 , 72 , 130 , 35 , 22 , 25 , 132 , 37 , 80 , 26 , 38 , 257 , 67 , 136 , 41 , 28 , 258 , 96 , 69 , 42 , 15 , 144 , 49 , 260 , 70 , 44 , 73 , 131 , 50 , 23 , 264 , 74 , 160 , 133 , 52 , 81 , 27 , 76 , 134 , 39 , 272 , 82 , 137 , 56 , 29 , 259 , 192 , 97 , 138 , 84 , 43 , 30 , 145 , 288 , 98 , 261 , 71 , 140 , 45 , 88 , 146 , 51 , 262 , 100 , 46 , 265 , 75 , 161 , 148 , 53 , 320 , 266 , 104 , 77 , 162 , 135 , 54 , 273 , 83 , 152 , 57 , 268 , 78 , 164 , 274 , 193 , 112 , 139 , 85 , 58 , 31 , 384 , 289 , 194 , 99 , 276 , 168 , 86 , 141 , 60 , 89 , 147 , 290 , 263 , 196 , 101 , 142 , 47 , 280 , 90 , 176 , 149 , 292 , 102 , 321 , 267 , 200 , 105 , 92 , 163 , 150 , 55 , 322 , 153 , 296 , 106 , 269 , 79 , 165 , 275 , 208 , 113 , 154 , 324 , 59 , 270 , 108 , 166 , 385 , 304 , 195 , 114 , 277 , 169 , 87 , 156 , 61 , 328 , 386 , 291 , 224 , 278 , 197 , 170 , 116 , 143 , 62 , 281 , 91 , 177 , 388 , 293 , 198 , 103 , 336 , 172 , 282 , 201 , 120 , 93 , 178 , 151 , 294 , 323 , 392 , 297 , 202 , 107 , 284 , 94 , 180 , 209 , 352 , 155 , 325 , 298 , 271 , 204 , 109 , 167 , 400 , 305 , 210 , 115 , 184 , 326 , 157 , 300 , 110 , 329 , 387 , 306 , 225 , 279 , 212 , 171 , 117 , 158 , 63 , 330 , 416 , 226 , 389 , 308 , 199 , 118 , 337 , 173 , 283 , 216 , 121 , 332 , 179 , 390 , 295 , 228 , 338 , 174 , 393 , 312 , 203 , 122 , 285 , 95 , 181 , 448 , 353 , 394 , 340 , 299 , 232 , 286 , 205 , 124 , 182 , 401 , 211 , 354 , 185 , 327 , 396 , 301 , 206 , 111 , 344 , 402 , 307 , 240 , 213 , 186 , 356 , 159 , 302 , 331 , 417 , 227 , 404 , 309 , 214 , 119 , 188 , 217 , 360 , 333 , 418 , 391 , 310 , 229 , 339 , 175 , 408 , 313 , 218 , 123 , 334 , 420 , 230 , 449 , 368 , 395 , 341 , 314 , 233 , 287 , 220 , 125 , 183 , 450 , 355 , 424 , 342 , 234 , 397 , 316 , 207 , 126 , 345 , 403 , 241 , 452 , 187 , 357 , 398 , 303 , 236 , 346 , 432 , 242 , 405 , 215 , 358 , 189 , 456 , 361 , 348 , 419 , 406 , 311 , 244 , 190 , 409 , 219 , 362 , 335 , 421 , 231 , 464 , 369 , 410 , 315 , 248 , 221 , 364 , 422 , 451 , 370 , 425 , 343 , 235 , 412 , 317 , 222 , 127 , 480 , 453 , 426 , 372 , 399 , 318 , 237 , 347 , 433 , 243 , 454 , 359 , 428 , 238 , 457 , 376 , 349 , 434 , 407 , 245 , 191 , 458 , 363 , 350 , 436 , 246 , 465 , 411 , 249 , 460 , 365 , 423 , 466 , 371 , 440 , 250 , 413 , 223 , 366 , 481 , 468 , 427 , 373 , 414 , 319 , 252 , 482 , 455 , 374 , 429 , 239 , 472 , 377 , 435 , 484 , 430 , 459 , 378 , 351 , 437 , 247 , 488 , 461 , 380 , 438 , 467 , 441 , 251 , 462 , 367 , 496 , 469 , 442 , 415 , 253 , 483 , 470 , 375 , 444 , 254 , 473 , 485 , 431 , 474 , 379 , 486 , 489 , 476 , 381 , 439 , 490 , 463 , 382 , 497 , 443 , 492 , 498 , 471 , 445 , 255 , 500 , 446 , 475 , 487 , 504 , 477 , 491 , 478 , 383 , 493 , 499 , 494 , 501 , 447 , 502 , 505 , 506 , 479 , 508 , 495 , 503 , 507 , 509 , 510 , 511 , 
	//0 , 1 , 2 , 4 , 8 , 16 , 3 , 32 , 5 , 6 , 9 , 64 , 10 , 17 , 12 , 18 , 128 , 33 , 20 , 34 , 7 , 24 , 36 , 65 , 11 , 256 , 66 , 40 , 13 , 19 , 68 , 14 , 129 , 48 , 21 , 72 , 130 , 35 , 22 , 25 , 512 , 132 , 37 , 80 , 26 , 38 , 257 , 67 , 136 , 41 , 28 , 258 , 96 , 69 , 42 , 15 , 144 , 49 , 260 , 70 , 44 , 73 , 131 , 50 , 23 , 264 , 74 , 160 , 513 , 133 , 52 , 81 , 27 , 76 , 514 , 134 , 39 , 272 , 82 , 137 , 56 , 29 , 516 , 259 , 192 , 97 , 138 , 84 , 43 , 30 , 145 , 288 , 98 , 261 , 71 , 520 , 140 , 45 , 88 , 146 , 51 , 262 , 100 , 46 , 265 , 75 , 161 , 528 , 148 , 53 , 320 , 266 , 104 , 77 , 162 , 515 , 135 , 54 , 273 , 83 , 152 , 57 , 268 , 78 , 544 , 164 , 517 , 274 , 193 , 112 , 139 , 85 , 58 , 31 , 518 , 384 , 289 , 194 , 99 , 276 , 168 , 86 , 521 , 141 , 60 , 89 , 147 , 290 , 576 , 263 , 196 , 101 , 522 , 142 , 47 , 280 , 90 , 176 , 529 , 149 , 292 , 102 , 321 , 524 , 267 , 200 , 105 , 92 , 163 , 530 , 150 , 55 , 322 , 153 , 296 , 106 , 269 , 79 , 640 , 545 , 165 , 532 , 275 , 208 , 113 , 154 , 324 , 59 , 270 , 108 , 546 , 166 , 519 , 385 , 304 , 195 , 114 , 277 , 169 , 87 , 536 , 156 , 61 , 328 , 548 , 386 , 291 , 224 , 577 , 278 , 197 , 170 , 116 , 523 , 143 , 62 , 281 , 91 , 177 , 768 , 578 , 388 , 293 , 198 , 103 , 552 , 336 , 172 , 525 , 282 , 201 , 120 , 93 , 178 , 531 , 151 , 294 , 580 , 323 , 526 , 392 , 297 , 202 , 107 , 284 , 94 , 641 , 560 , 180 , 533 , 209 , 352 , 155 , 325 , 298 , 584 , 271 , 204 , 109 , 642 , 547 , 167 , 534 , 400 , 305 , 210 , 115 , 184 , 326 , 537 , 157 , 300 , 110 , 329 , 644 , 549 , 387 , 306 , 225 , 592 , 279 , 212 , 171 , 117 , 538 , 158 , 63 , 330 , 550 , 416 , 769 , 226 , 579 , 389 , 308 , 199 , 648 , 118 , 553 , 337 , 173 , 540 , 283 , 216 , 121 , 332 , 179 , 770 , 608 , 390 , 295 , 228 , 581 , 554 , 338 , 174 , 527 , 393 , 312 , 203 , 122 , 285 , 95 , 656 , 561 , 181 , 772 , 582 , 448 , 353 , 556 , 394 , 340 , 299 , 232 , 585 , 286 , 205 , 124 , 643 , 562 , 182 , 535 , 401 , 211 , 354 , 185 , 327 , 776 , 586 , 396 , 301 , 206 , 111 , 672 , 344 , 645 , 564 , 402 , 307 , 240 , 593 , 213 , 186 , 356 , 539 , 159 , 302 , 588 , 331 , 646 , 551 , 417 , 784 , 227 , 594 , 404 , 309 , 214 , 649 , 119 , 568 , 188 , 541 , 217 , 360 , 333 , 418 , 771 , 704 , 609 , 391 , 310 , 229 , 650 , 596 , 555 , 339 , 175 , 542 , 408 , 313 , 218 , 123 , 334 , 657 , 800 , 610 , 420 , 773 , 230 , 583 , 449 , 368 , 652 , 557 , 395 , 341 , 314 , 233 , 600 , 287 , 220 , 125 , 658 , 563 , 183 , 774 , 612 , 450 , 355 , 558 , 424 , 342 , 777 , 234 , 587 , 397 , 316 , 207 , 126 , 673 , 345 , 660 , 565 , 403 , 241 , 832 , 452 , 187 , 357 , 778 , 616 , 398 , 303 , 236 , 589 , 674 , 346 , 647 , 566 , 432 , 785 , 242 , 595 , 405 , 215 , 664 , 358 , 569 , 189 , 780 , 590 , 456 , 361 , 676 , 348 , 419 , 786 , 705 , 624 , 406 , 311 , 244 , 651 , 597 , 570 , 190 , 543 , 409 , 219 , 362 , 335 , 896 , 801 , 706 , 611 , 421 , 788 , 231 , 680 , 598 , 464 , 369 , 653 , 572 , 410 , 315 , 248 , 601 , 221 , 364 , 659 , 802 , 422 , 775 , 708 , 613 , 451 , 370 , 654 , 559 , 425 , 343 , 792 , 235 , 602 , 412 , 317 , 222 , 127 , 688 , 661 , 804 , 614 , 480 , 833 , 453 , 426 , 372 , 779 , 712 , 617 , 399 , 318 , 237 , 604 , 675 , 347 , 662 , 567 , 433 , 243 , 834 , 454 , 665 , 359 , 808 , 618 , 428 , 781 , 238 , 591 , 457 , 376 , 677 , 349 , 434 , 787 , 720 , 625 , 407 , 245 , 666 , 836 , 571 , 191 , 782 , 620 , 458 , 363 , 678 , 350 , 897 , 816 , 707 , 626 , 436 , 789 , 246 , 681 , 599 , 465 , 668 , 573 , 411 , 249 , 840 , 460 , 365 , 898 , 803 , 736 , 423 , 790 , 709 , 682 , 628 , 466 , 371 , 655 , 574 , 440 , 793 , 250 , 603 , 413 , 223 , 366 , 689 , 900 , 805 , 710 , 615 , 481 , 848 , 684 , 468 , 427 , 373 , 794 , 713 , 632 , 414 , 319 , 252 , 605 , 690 , 663 , 806 , 482 , 835 , 455 , 904 , 374 , 809 , 714 , 619 , 429 , 796 , 239 , 606 , 472 , 377 , 692 , 435 , 721 , 864 , 484 , 667 , 837 , 810 , 430 , 783 , 716 , 621 , 459 , 378 , 679 , 351 , 912 , 817 , 722 , 627 , 437 , 247 , 696 , 838 , 669 , 812 , 622 , 488 , 841 , 461 , 380 , 899 , 818 , 737 , 438 , 791 , 724 , 683 , 629 , 467 , 670 , 575 , 441 , 251 , 842 , 462 , 367 , 928 , 738 , 901 , 820 , 711 , 630 , 496 , 849 , 685 , 469 , 442 , 795 , 728 , 633 , 415 , 253 , 844 , 691 , 902 , 807 , 740 , 483 , 850 , 686 , 470 , 905 , 375 , 824 , 715 , 634 , 444 , 797 , 254 , 607 , 473 , 693 , 960 , 865 , 485 , 906 , 852 , 811 , 744 , 431 , 798 , 717 , 636 , 474 , 379 , 694 , 913 , 723 , 866 , 486 , 697 , 839 , 908 , 813 , 718 , 623 , 489 , 856 , 476 , 381 , 914 , 819 , 752 , 439 , 725 , 698 , 868 , 671 , 814 , 490 , 843 , 463 , 382 , 929 , 739 , 916 , 821 , 726 , 631 , 497 , 700 , 443 , 729 , 872 , 492 , 845 , 930 , 903 , 822 , 741 , 498 , 851 , 687 , 471 , 920 , 825 , 730 , 635 , 445 , 255 , 846 , 932 , 742 , 961 , 880 , 500 , 907 , 853 , 826 , 745 , 446 , 799 , 732 , 637 , 475 , 695 , 962 , 867 , 487 , 936 , 854 , 746 , 909 , 828 , 719 , 638 , 504 , 857 , 477 , 915 , 753 , 964 , 699 , 869 , 910 , 815 , 748 , 491 , 858 , 478 , 383 , 944 , 754 , 917 , 727 , 870 , 701 , 968 , 873 , 493 , 860 , 931 , 918 , 823 , 756 , 499 , 702 , 921 , 731 , 874 , 494 , 847 , 933 , 743 , 976 , 881 , 501 , 922 , 827 , 760 , 447 , 733 , 876 , 934 , 963 , 882 , 502 , 937 , 855 , 747 , 924 , 829 , 734 , 639 , 505 , 992 , 965 , 938 , 884 , 911 , 830 , 749 , 506 , 859 , 479 , 945 , 755 , 966 , 871 , 940 , 750 , 969 , 888 , 508 , 861 , 946 , 919 , 757 , 703 , 970 , 875 , 495 , 862 , 948 , 758 , 977 , 923 , 761 , 972 , 877 , 935 , 978 , 883 , 503 , 952 , 762 , 925 , 735 , 878 , 993 , 980 , 939 , 885 , 926 , 831 , 764 , 507 , 994 , 967 , 886 , 941 , 751 , 984 , 889 , 509 , 947 , 996 , 942 , 971 , 890 , 510 , 863 , 949 , 759 , 1000 , 973 , 892 , 950 , 979 , 953 , 763 , 974 , 879 , 1008 , 981 , 954 , 927 , 765 , 995 , 982 , 887 , 956 , 766 , 985 , 997 , 943 , 986 , 891 , 511 , 998 , 1001 , 988 , 893 , 951 , 1002 , 975 , 894 , 1009 , 955 , 1004 , 1010 , 983 , 957 , 767 , 1012 , 958 , 987 , 999 , 1016 , 989 , 1003 , 990 , 895 , 1005 , 1011 , 1006 , 1013 , 959 , 1014 , 1017 , 1018 , 991 , 1020 , 1007 , 1015 , 1019 , 1021 , 1022 , 1023 , 
	//0 , 1 , 2 , 4 , 8 , 16 , 3 , 32 , 5 , 6 , 9 , 64 , 10 , 17 , 12 , 18 , 128 , 33 , 20 , 34 , 7 , 24 , 36 , 65 , 11 , 256 , 66 , 40 , 13 , 19 , 68 , 14 , 129 , 48 , 21 , 72 , 130 , 35 , 22 , 25 , 512 , 132 , 37 , 80 , 26 , 38 , 257 , 67 , 136 , 41 , 28 , 258 , 96 , 69 , 42 , 15 , 144 , 49 , 260 , 70 , 44 , 73 , 131 , 50 , 23 , 1024 , 264 , 74 , 160 , 513 , 133 , 52 , 81 , 27 , 76 , 514 , 134 , 39 , 272 , 82 , 137 , 56 , 29 , 516 , 259 , 192 , 97 , 138 , 84 , 43 , 30 , 145 , 288 , 98 , 261 , 71 , 520 , 140 , 45 , 88 , 146 , 51 , 262 , 100 , 1025 , 46 , 265 , 75 , 161 , 528 , 148 , 53 , 320 , 1026 , 266 , 104 , 77 , 162 , 515 , 135 , 54 , 273 , 83 , 152 , 57 , 1028 , 268 , 78 , 544 , 164 , 517 , 274 , 193 , 112 , 139 , 85 , 58 , 31 , 1032 , 518 , 384 , 289 , 194 , 99 , 276 , 168 , 86 , 521 , 141 , 60 , 89 , 147 , 290 , 576 , 263 , 196 , 101 , 522 , 142 , 1040 , 47 , 280 , 90 , 176 , 529 , 149 , 292 , 102 , 321 , 1027 , 524 , 267 , 200 , 105 , 92 , 163 , 530 , 150 , 55 , 322 , 1056 , 153 , 296 , 1029 , 106 , 269 , 79 , 640 , 545 , 165 , 532 , 275 , 208 , 113 , 154 , 324 , 59 , 1030 , 270 , 108 , 546 , 1033 , 166 , 519 , 385 , 304 , 195 , 114 , 277 , 169 , 87 , 536 , 156 , 61 , 1088 , 328 , 1034 , 548 , 386 , 291 , 224 , 577 , 278 , 197 , 170 , 116 , 523 , 143 , 1041 , 62 , 281 , 91 , 177 , 1036 , 768 , 578 , 388 , 293 , 198 , 103 , 552 , 336 , 172 , 1042 , 525 , 282 , 201 , 120 , 93 , 178 , 531 , 151 , 294 , 580 , 323 , 1152 , 1057 , 526 , 392 , 297 , 202 , 1044 , 107 , 284 , 94 , 641 , 560 , 180 , 533 , 209 , 352 , 1058 , 155 , 325 , 298 , 1031 , 584 , 271 , 204 , 109 , 642 , 547 , 1048 , 167 , 534 , 400 , 305 , 210 , 115 , 184 , 326 , 537 , 1060 , 157 , 300 , 1089 , 110 , 329 , 1035 , 644 , 549 , 387 , 306 , 225 , 592 , 279 , 212 , 171 , 117 , 538 , 158 , 1280 , 63 , 1090 , 330 , 1064 , 550 , 416 , 1037 , 769 , 226 , 579 , 389 , 308 , 199 , 648 , 118 , 553 , 337 , 173 , 1043 , 540 , 283 , 216 , 121 , 1092 , 332 , 179 , 1038 , 770 , 608 , 390 , 295 , 228 , 581 , 554 , 338 , 1153 , 174 , 1072 , 527 , 393 , 312 , 203 , 1045 , 122 , 285 , 95 , 656 , 561 , 181 , 1096 , 772 , 582 , 448 , 353 , 1154 , 1059 , 556 , 394 , 340 , 299 , 232 , 1046 , 585 , 286 , 205 , 124 , 643 , 562 , 1049 , 182 , 535 , 401 , 211 , 354 , 1536 , 185 , 327 , 1156 , 776 , 1061 , 586 , 396 , 301 , 206 , 1104 , 111 , 672 , 344 , 1050 , 645 , 564 , 402 , 307 , 240 , 593 , 213 , 186 , 356 , 539 , 1062 , 159 , 1281 , 302 , 1091 , 588 , 331 , 1160 , 1065 , 646 , 551 , 417 , 1052 , 784 , 227 , 594 , 404 , 309 , 214 , 649 , 119 , 568 , 188 , 1282 , 541 , 1120 , 217 , 360 , 1093 , 1066 , 333 , 418 , 1039 , 771 , 704 , 609 , 391 , 310 , 229 , 650 , 596 , 555 , 339 , 1168 , 175 , 1073 , 542 , 408 , 313 , 218 , 1284 , 123 , 1094 , 334 , 657 , 1068 , 800 , 610 , 420 , 1097 , 773 , 230 , 583 , 449 , 368 , 1155 , 652 , 1074 , 557 , 395 , 341 , 314 , 233 , 1047 , 600 , 287 , 220 , 125 , 658 , 563 , 1288 , 183 , 1098 , 774 , 612 , 450 , 355 , 1184 , 1537 , 558 , 424 , 342 , 1157 , 777 , 234 , 1076 , 587 , 397 , 316 , 207 , 1105 , 126 , 673 , 345 , 1051 , 660 , 565 , 403 , 241 , 1100 , 832 , 1538 , 452 , 187 , 357 , 1158 , 778 , 1063 , 616 , 398 , 1296 , 303 , 236 , 1106 , 589 , 674 , 346 , 1161 , 1080 , 647 , 566 , 432 , 1053 , 785 , 242 , 595 , 405 , 215 , 664 , 358 , 569 , 1540 , 189 , 1283 , 1216 , 780 , 1121 , 590 , 456 , 361 , 1162 , 1108 , 1067 , 676 , 348 , 419 , 1054 , 786 , 705 , 624 , 406 , 311 , 244 , 651 , 597 , 570 , 1169 , 190 , 1312 , 543 , 409 , 1122 , 219 , 1285 , 362 , 1095 , 1544 , 335 , 1164 , 896 , 1069 , 801 , 706 , 611 , 421 , 1112 , 788 , 231 , 680 , 598 , 464 , 369 , 1170 , 653 , 1075 , 572 , 410 , 315 , 248 , 1286 , 601 , 1124 , 221 , 364 , 659 , 1070 , 802 , 1289 , 422 , 1099 , 775 , 708 , 613 , 451 , 370 , 1185 , 654 , 1552 , 559 , 425 , 343 , 1172 , 792 , 235 , 1077 , 602 , 412 , 317 , 222 , 1344 , 127 , 688 , 1290 , 661 , 1128 , 804 , 614 , 480 , 1101 , 833 , 1186 , 1539 , 453 , 426 , 372 , 1159 , 779 , 712 , 1078 , 617 , 399 , 1297 , 318 , 237 , 1107 , 604 , 675 , 347 , 1176 , 1081 , 662 , 567 , 433 , 1292 , 243 , 1102 , 834 , 1568 , 454 , 665 , 359 , 1188 , 808 , 1541 , 618 , 428 , 1298 , 1217 , 781 , 238 , 1136 , 591 , 457 , 376 , 1163 , 1109 , 1082 , 677 , 349 , 434 , 1055 , 787 , 720 , 625 , 407 , 245 , 666 , 836 , 571 , 1542 , 1408 , 191 , 1313 , 1218 , 782 , 1123 , 620 , 458 , 1300 , 363 , 1192 , 1110 , 1545 , 678 , 350 , 1165 , 897 , 1084 , 816 , 707 , 626 , 436 , 1113 , 789 , 246 , 681 , 599 , 465 , 1171 , 668 , 1314 , 573 , 411 , 1600 , 249 , 1287 , 1220 , 840 , 1125 , 1546 , 460 , 365 , 1166 , 898 , 1071 , 803 , 736 , 1304 , 423 , 1114 , 790 , 709 , 682 , 628 , 466 , 371 , 1200 , 655 , 1553 , 574 , 440 , 1173 , 793 , 250 , 1316 , 603 , 413 , 1126 , 223 , 1345 , 366 , 689 , 1548 , 1291 , 1224 , 900 , 1129 , 805 , 710 , 615 , 481 , 1116 , 848 , 1187 , 684 , 1554 , 468 , 427 , 373 , 1174 , 794 , 713 , 1079 , 632 , 414 , 319 , 252 , 1346 , 605 , 690 , 1177 , 1320 , 663 , 1130 , 806 , 1293 , 482 , 1103 , 835 , 1664 , 1569 , 455 , 904 , 374 , 1189 , 809 , 714 , 1556 , 619 , 429 , 1299 , 1232 , 796 , 239 , 1137 , 606 , 472 , 377 , 1178 , 1348 , 1083 , 692 , 435 , 1294 , 721 , 1132 , 864 , 1570 , 484 , 667 , 837 , 1190 , 810 , 1543 , 1409 , 430 , 1328 , 1219 , 783 , 716 , 1138 , 621 , 459 , 1301 , 378 , 1193 , 1111 , 1560 , 679 , 351 , 1180 , 912 , 1085 , 817 , 722 , 627 , 437 , 1352 , 247 , 696 , 838 , 1572 , 1410 , 669 , 1315 , 1248 , 812 , 1601 , 622 , 488 , 1302 , 1221 , 841 , 1194 , 1140 , 1547 , 461 , 380 , 1167 , 899 , 1086 , 818 , 737 , 1305 , 438 , 1115 , 791 , 724 , 683 , 629 , 467 , 1201 , 670 , 1792 , 575 , 441 , 1602 , 1412 , 251 , 1317 , 1222 , 842 , 1127 , 1576 , 462 , 1360 , 367 , 1196 , 928 , 1549 , 738 , 1306 , 1225 , 901 , 1144 , 820 , 711 , 630 , 496 , 1117 , 849 , 1202 , 685 , 1555 , 469 , 442 , 1175 , 795 , 728 , 1318 , 633 , 415 , 1604 , 253 , 1347 , 844 , 691 , 1550 , 1416 , 1321 , 1226 , 902 , 1131 , 807 , 740 , 1308 , 483 , 1118 , 850 , 1665 , 686 , 1584 , 470 , 905 , 375 , 1204 , 824 , 715 , 1557 , 634 , 444 , 1233 , 797 , 254 , 1376 , 607 , 473 , 1179 , 1349 , 1322 , 693 , 1608 , 1295 , 1228 , 960 , 1133 , 865 , 1666 , 1571 , 485 , 906 , 852 , 1191 , 811 , 744 , 1558 , 1424 , 431 , 1329 , 1234 , 798 , 717 , 1139 , 636 , 474 , 379 , 1208 , 1350 , 1561 , 694 , 1181 , 913 , 1324 , 723 , 1134 , 866 , 1353 , 486 , 697 , 839 , 1668 , 1573 , 1411 , 908 , 1330 , 1249 , 813 , 718 , 1616 , 623 , 489 , 1303 , 1236 , 856 , 1195 , 1141 , 1562 , 476 , 381 , 1182 , 914 , 1087 , 819 , 752 , 439 , 1354 , 725 , 698 , 868 , 1574 , 1440 , 671 , 1793 , 1250 , 814 , 1603 , 1413 , 490 , 1332 , 1223 , 843 , 1672 , 1142 , 1577 , 463 , 1361 , 382 , 1197 , 929 , 1564 , 739 , 1307 , 1240 , 916 , 1145 , 821 , 726 , 631 , 497 , 1356 , 1203 , 700 , 1794 , 443 , 1632 , 1414 , 729 , 1319 , 1252 , 872 , 1605 , 1578 , 492 , 1362 , 845 , 1198 , 930 , 1551 , 1417 , 1336 , 1227 , 903 , 1146 , 822 , 741 , 1309 , 498 , 1119 , 851 , 1680 , 687 , 1585 , 471 , 920 , 1205 , 825 , 730 , 1796 , 635 , 445 , 1606 , 1472 , 255 , 1377 , 846 , 1580 , 1418 , 1364 , 1323 , 1256 , 932 , 1609 , 742 , 1310 , 1229 , 961 , 1148 , 880 , 1667 , 1586 , 500 , 907 , 853 , 1206 , 826 , 745 , 1559 , 1425 , 446 , 1235 , 799 , 732 , 1378 , 637 , 475 , 1209 , 1351 , 1800 , 695 , 1610 , 1420 , 1325 , 1230 , 962 , 1135 , 867 , 1696 , 1368 , 487 , 936 , 854 , 1669 , 746 , 1588 , 1426 , 909 , 1331 , 1264 , 828 , 719 , 1617 , 638 , 504 , 1237 , 857 , 1210 , 1380 , 1563 , 477 , 1183 , 915 , 1326 , 753 , 1612 , 1355 , 964 , 699 , 869 , 1670 , 1575 , 1441 , 910 , 1808 , 1251 , 815 , 748 , 1618 , 1428 , 491 , 1333 , 1238 , 858 , 1673 , 1143 , 1592 , 478 , 383 , 1212 , 944 , 1565 , 754 , 1241 , 917 , 1384 , 727 , 870 , 1357 , 1442 , 701 , 1795 , 1728 , 1633 , 1415 , 968 , 1334 , 1253 , 873 , 1674 , 1620 , 1579 , 493 , 1363 , 860 , 1199 , 931 , 1566 , 1432 , 1337 , 1242 , 918 , 1147 , 823 , 756 , 499 , 1358 , 1681 , 702 , 1824 , 921 , 1634 , 1444 , 731 , 1797 , 1254 , 874 , 1607 , 1473 , 494 , 1392 , 847 , 1676 , 1581 , 1419 , 1365 , 1338 , 1257 , 933 , 1624 , 743 , 1311 , 1244 , 976 , 1149 , 881 , 1682 , 1587 , 501 , 922 , 1207 , 827 , 760 , 1798 , 447 , 1636 , 1474 , 733 , 1379 , 876 , 1582 , 1448 , 1366 , 1801 , 1258 , 934 , 1611 , 1421 , 1340 , 1231 , 963 , 1150 , 882 , 1697 , 1369 , 502 , 937 , 855 , 1684 , 747 , 1589 , 1427 , 924 , 1265 , 829 , 734 , 1856 , 639 , 505 , 1476 , 1211 , 1381 , 1802 , 1640 , 1422 , 1327 , 1260 , 992 , 1613 , 1698 , 1370 , 965 , 938 , 884 , 1671 , 1590 , 1456 , 911 , 1809 , 1266 , 830 , 749 , 1619 , 1429 , 506 , 1239 , 859 , 1688 , 1382 , 1593 , 479 , 1213 , 945 , 1804 , 755 , 1614 , 1480 , 1385 , 966 , 871 , 1700 , 1372 , 1443 , 940 , 1810 , 1729 , 750 , 1648 , 1430 , 969 , 1335 , 1268 , 888 , 1675 , 1621 , 1594 , 508 , 861 , 1214 , 946 , 1567 , 1433 , 1243 , 919 , 1386 , 757 , 1359 , 1920 , 703 , 1825 , 1730 , 1635 , 1445 , 970 , 1812 , 1255 , 875 , 1704 , 1622 , 1488 , 495 , 1393 , 862 , 1677 , 1596 , 1434 , 1339 , 1272 , 948 , 1625 , 758 , 1245 , 977 , 1388 , 1683 , 1826 , 923 , 1446 , 761 , 1799 , 1732 , 1637 , 1475 , 972 , 1394 , 877 , 1678 , 1583 , 1449 , 1367 , 1816 , 1259 , 935 , 1626 , 1436 , 1341 , 1246 , 978 , 1151 , 883 , 1712 , 503 , 952 , 1685 , 762 , 1828 , 925 , 1638 , 1504 , 735 , 1857 , 878 , 1477 , 1450 , 1396 , 1803 , 1736 , 1641 , 1423 , 1342 , 1261 , 993 , 1628 , 1699 , 1371 , 980 , 939 , 885 , 1686 , 1591 , 1457 , 926 , 1267 , 831 , 764 , 1858 , 507 , 1478 , 1689 , 1383 , 1832 , 1642 , 1452 , 1805 , 1262 , 994 , 1615 , 1481 , 1400 , 967 , 886 , 1701 , 1373 , 1458 , 941 , 1811 , 1744 , 751 , 1649 , 1431 , 984 , 1269 , 889 , 1690 , 1860 , 1595 , 509 , 1215 , 947 , 1806 , 1644 , 1482 , 1387 , 996 , 1702 , 1374 , 1921 , 942 , 1840 , 1731 , 1650 , 1460 , 971 , 1813 , 1270 , 890 , 1705 , 1623 , 1489 , 510 , 863 , 1692 , 1597 , 1435 , 1273 , 949 , 1864 , 759 , 1484 , 1389 , 1922 , 1827 , 1760 , 1447 , 1000 , 1814 , 1733 , 1706 , 1652 , 1490 , 973 , 1395 , 892 , 1679 , 1598 , 1464 , 1817 , 1274 , 950 , 1627 , 1437 , 1247 , 979 , 1390 , 1713 , 953 , 1924 , 763 , 1829 , 1734 , 1639 , 1505 , 974 , 1872 , 879 , 1708 , 1492 , 1451 , 1397 , 1818 , 1737 , 1656 , 1438 , 1343 , 1276 , 1008 , 1629 , 1714 , 981 , 954 , 1687 , 1830 , 927 , 1506 , 765 , 1859 , 1479 , 1928 , 1398 , 1833 , 1738 , 1643 , 1453 , 1820 , 1263 , 995 , 1630 , 1496 , 1401 , 982 , 887 , 1716 , 1459 , 956 , 1745 , 766 , 1888 , 985 , 1508 , 1691 , 1861 , 1834 , 1454 , 1807 , 1740 , 1645 , 1483 , 1402 , 997 , 1703 , 1375 , 1936 , 943 , 1841 , 1746 , 1651 , 1461 , 986 , 1271 , 891 , 1720 , 1862 , 511 , 1693 , 1836 , 1646 , 1512 , 1865 , 998 , 1485 , 1404 , 1923 , 1842 , 1761 , 1462 , 1001 , 1815 , 1748 , 1707 , 1653 , 1491 , 988 , 893 , 1694 , 1599 , 1465 , 1275 , 951 , 1866 , 1486 , 1391 , 1952 , 1762 , 1925 , 1002 , 1844 , 1735 , 1654 , 1520 , 975 , 1873 , 894 , 1709 , 1493 , 1466 , 1819 , 1752 , 1657 , 1439 , 1277 , 1009 , 1868 , 1715 , 955 , 1926 , 1831 , 1764 , 1507 , 1004 , 1874 , 1710 , 1494 , 1929 , 1399 , 1848 , 1739 , 1658 , 1468 , 1821 , 1278 , 1010 , 1631 , 1497 , 983 , 1717 , 957 , 1984 , 767 , 1889 , 1509 , 1930 , 1876 , 1835 , 1768 , 1455 , 1822 , 1741 , 1660 , 1498 , 1403 , 1012 , 1718 , 1937 , 958 , 1747 , 1890 , 987 , 1510 , 1721 , 1863 , 1932 , 1837 , 1742 , 1647 , 1513 , 1880 , 999 , 1500 , 1405 , 1938 , 1843 , 1776 , 1463 , 1016 , 1749 , 1722 , 1892 , 989 , 1695 , 1838 , 1514 , 1867 , 1487 , 1406 , 1953 , 1763 , 1940 , 1003 , 1845 , 1750 , 1655 , 1521 , 990 , 895 , 1724 , 1467 , 1753 , 1896 , 1516 , 1869 , 1954 , 1927 , 1846 , 1765 , 1522 , 1005 , 1875 , 1711 , 1495 , 1944 , 1849 , 1754 , 1659 , 1469 , 1279 , 1011 , 1870 , 1956 , 1766 , 1985 , 1006 , 1904 , 1524 , 1931 , 1877 , 1850 , 1769 , 1470 , 1823 , 1756 , 1661 , 1499 , 1013 , 1719 , 959 , 1986 , 1891 , 1511 , 1960 , 1878 , 1770 , 1933 , 1852 , 1743 , 1662 , 1528 , 1881 , 1014 , 1501 , 1939 , 1777 , 1017 , 1988 , 1723 , 1893 , 1934 , 1839 , 1772 , 1515 , 1882 , 1502 , 1407 , 1968 , 1778 , 1941 , 1018 , 1751 , 1894 , 991 , 1725 , 1992 , 1897 , 1517 , 1884 , 1955 , 1942 , 1847 , 1780 , 1523 , 1020 , 1726 , 1945 , 1755 , 1898 , 1518 , 1871 , 1957 , 1767 , 2000 , 1007 , 1905 , 1525 , 1946 , 1851 , 1784 , 1471 , 1757 , 1900 , 1958 , 1987 , 1906 , 1526 , 1961 , 1879 , 1771 , 1948 , 1853 , 1758 , 1663 , 1529 , 1015 , 2016 , 1989 , 1962 , 1908 , 1935 , 1854 , 1773 , 1530 , 1883 , 1503 , 1969 , 1779 , 1019 , 1990 , 1895 , 1964 , 1774 , 1993 , 1912 , 1532 , 1885 , 1970 , 1943 , 1781 , 1021 , 1727 , 1994 , 1899 , 1519 , 1886 , 1972 , 1782 , 2001 , 1022 , 1947 , 1785 , 1996 , 1901 , 1959 , 2002 , 1907 , 1527 , 1976 , 1786 , 1949 , 1759 , 1902 , 2017 , 2004 , 1963 , 1909 , 1950 , 1855 , 1788 , 1531 , 2018 , 1991 , 1910 , 1965 , 1775 , 2008 , 1913 , 1533 , 1971 , 2020 , 1966 , 1995 , 1914 , 1534 , 1887 , 1973 , 1783 , 1023 , 2024 , 1997 , 1916 , 1974 , 2003 , 1977 , 1787 , 1998 , 1903 , 2032 , 2005 , 1978 , 1951 , 1789 , 2019 , 2006 , 1911 , 1980 , 1790 , 2009 , 2021 , 1967 , 2010 , 1915 , 1535 , 2022 , 2025 , 2012 , 1917 , 1975 , 2026 , 1999 , 1918 , 2033 , 1979 , 2028 , 2034 , 2007 , 1981 , 1791 , 2036 , 1982 , 2011 , 2023 , 2040 , 2013 , 2027 , 2014 , 1919 , 2029 , 2035 , 2030 , 2037 , 1983 , 2038 , 2041 , 2042 , 2015 , 2044 , 2031 , 2039 , 2043 , 2045 , 2046 , 2047 , 

	
	};
	vector <int> output(k1,0);
	for(int i=0;i<FrozenBits.size();i++)
	{
		output[FrozenBits[i]]=1;
	}
	return output;
}
vector <int> Polar_information_bit_creator(int k1,vector <int> Frozen_bits_polar, double condition)
{
	vector <int> output(k1,1);
	for (int i=0;i<k1;i++)
	{
		if(Frozen_bits_polar[i]==0)
		{
			double a=(rand()+0.001)/(RAND_MAX+0.001);
			if(a<condition)
			{
				output[i]=-1;
			}
		}
	}
	return output;
}

int precoding_info_bits(int k1,vector <int> Frozen_bits_polar)
{
	int output=k1;
	for(int i=0;i<k1;i++)
	{
		output-=Frozen_bits_polar[i];
	}
	return output;
}
//Polar Codes
/////////////////////////////////////////////////////

