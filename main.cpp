#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>
using namespace std;


const int N = 20;

void gauss_jordan(double a[][N+1], int n) {
    for (int i = 0; i < n; i++) {
        // find pivot row and swap
        int pivot = i;
        for (int j = i+1; j < n; j++) {
            if (abs(a[j][i]) > abs(a[pivot][i])) {
                pivot = j;
            }
        }
        if (pivot != i) {
            for (int j = 0; j <= n; j++) {
                swap(a[i][j], a[pivot][j]);
            }
        }

        // reduce pivot row
        for (int j = i+1; j <= n; j++) {
            a[i][j] /= a[i][i];
        }
        a[i][i] = 1;

        // eliminate other rows
        for (int j = 0; j < n; j++) {
            if (j != i) {
                for (int k = i+1; k <= n; k++) {
                    a[j][k] -= a[j][i] * a[i][k];
                }
                a[j][i] = 0;
            }
        }
    }
}

//////////////////...functions of T Saturation...//////////////////
//from 0.01 bar to 0.05 bar
double TsatL(double p){
	double T=16.1*log(p)+80.732;
	return T;
}
//from 0.05 bar to 0.1 bar
double TsatL2(double p){
	double T=21.18*log(p)+94.601;
	return T;
}
//from 0.5 bar to 8 bar
double TsatM (double p){
	double T=-0.0041*pow(p,6)+0.1182*pow(p,5)-1.3806*pow(p,4)+8.3997*pow(p,3)-29.223*pow(p,2)+67.089*p+54.426;
	return T;
}
//from 8 bar to 90 bar
double TsatH (double p){
	double T=pow(10,-7)*pow(p,5)-3*pow(10,-5)*pow(p,4)+0.0032*pow(p,3)-0.1866*pow(p,2)+7.0347*p+125.2;
	return T;
}
double Tsat(double P){
	double T=0;
	if (P<=0.05)
	T=TsatL(P);
	if (P>0.05 && P<=0.4)
	T=TsatL2(P);
	if (P>0.4 && P<=8)
	T=TsatM(P);
	if (P>8 && P<=90)
	T=TsatH(P);
	return T;
	}
//////////////////...functions of H Saturation...//////////////////
//from 0.01 bar to 0.05 bar
double HsatVL(double p){
	double H=29.237*log(p)+2647.7;
	return H;
}
//from 0.05 bar to 0.4 bar
double HsatVL1(double p){
	double H=37.018*log(p)+2669.1;
	return H;
}
//from 0.5 bar to 14 bar
double HsatVM (double p){
	double H=-0.0005*pow(p,6)+0.0244*pow(p,5)-0.4709*pow(p,4)+4.6477*pow(p,3)-25.421*pow(p,2)+83.309*p+2611.4;
	return H;
}
//from 14 bar to 110 bar
double HsatVH (double p){
	double H=3*pow(10,-8)*pow(p,5)-pow(10,-5)*pow(p,4)+0.0018*pow(p,3)-0.1427*pow(p,2)+4.947*p+2743.6;
	return H;
}
double HsatV(double P){
	double H=0;
	if (P<=0.05)
	H=HsatVL(P);
	if (P>0.05 && P<=0.4)
	H=HsatVL1(P);
	if (P>0.4 && P<=8)
	H=HsatVM(P);
	if (P>8 && P<90)
	H=HsatVH(P);
	return H;
	
}
//////////////////...functions of P Saturation in Specific H of water...//////////////////
//from 0.01 bar to 0.1 bar
double PsatLHL(double h){
	double P=pow(10,-8)*pow(h,3)-3*pow(10,-7)*h*h+0.0002*h+0.0048;
	return P;
}
//from 0.1 bar to 1 bar
double PsatLHM(double h){
	double P=4*pow(10,-8)*pow(h,3)-2*pow(10,-5)*h*h+0.0042*h-0.3188;
	return P;
}
//from 1 bar to 110 bar
double PsatLHH(double h){
	double P=-5*pow(10,-11)*pow(h,4)+2*pow(10,-7)*pow(h,3)-0.0002*h*h+0.1069*h-18.051;
	return P;
}
double PsatLH(double h){
	double P=0;
	if (h<=200)
	P=PsatLHL(h);
	if (h>200 && h<=420)
	P=PsatLHM(h);
	if (h>420 && h<=1500)
	P=PsatLHH(h);
	return P;
	
}

//////////////////...functions of H Liquid Saturation...//////////////////
//from 0.01 bar to 0.1 bar
double HsatLL(double p){
	double H=268667*pow(p,3)-61204*p*p+5528.7*p-17.116;
	return H;
}
//from 0.1 bar to 1 bar
double HsatLM1 (double p){
	double H=-485.88*pow(p,4)+1350.1*pow(p,3)-1458.2*pow(p,2)+890.32*p+121.18;
	return H;
}
//from 1 bar to 4 bar
double HsatLM2 (double p){
	double H=2.3588*pow(p,3)-27.925*pow(p,2)+151.65*p+294.19;
	return H;
}
//from 4 bar to 10 bar
double HsatLM3 (double p){
	double H=0.1189*pow(p,3)-3.9556*pow(p,2)+63.051*p+408.78;
	return H;
}
//from 10 bar to 110 bar
double HsatLH (double p){
	double H=-8*pow(10,-6)*pow(p,4)+0.0025*pow(p,3)-0.2988*pow(p,2)+20.987*p+586.89;
	return H;
}
double HsatL(double P){
	double H=0;
	if (P<=0.1)
	H=HsatLL(P);
	if (P>0.1 && P<=1)
	H=HsatLM1(P);
	if (P>1 && P<=4)
	H=HsatLM2(P);
	if (P>4 && P<=10)
	H=HsatLM3(P);
	if (P>10 && P<110)
	H=HsatLH(P);
	return H;
	
}
//////////////////...functions of S Liquid Saturation...//////////////////
//from 0.01 bar to 0.05 bar
double SsatLL(double p){
	double S=0.2301*log(p)+1.1627;
	return S;
}
//from 0.05 bar to 0.4 bar
double SsatLL1(double p){
	double S=0.2682*log(p)+1.2671;
	return S;
}
//from 0.5 bar to 4 bar
double SsatLM1 (double p){
	double S=-0.0061*pow(p,4)+0.0698*pow(p,3)-0.3124*pow(p,2)+0.7703*p+0.7791;
	return S;
}
//from 4 bar to 10 bar
double SsatLM2 (double p){
	double S=-0.0035*pow(p,2)+0.1079*p+1.4082;
	return S;
}
//from 10 bar to 110 bar
double SsatLH (double p){
	double S=-2*pow(10,-8)*pow(p,4)+6*pow(10,-6)*pow(p,3)-0.0007*pow(p,2)+0.0458*p+1.7617;
	return S;
}
double SsatL(double P){
	double S=0;
	if (P<=0.05)
	S=SsatLL(P);
	if (P>0.05 && P<=0.4)
	S=SsatLL1(P);
	if (P>0.4 && P<=4)
	S=SsatLM1(P);
	if (P>4 && P<=10)
	S=SsatLM2(P);
	if (P>10 && P<110)
	S=SsatLH(P);
	return S;
	
}
//////////////////...functions of S Vapor Saturation...//////////////////
//from 0.01 bar to 0.4 bar
double SsatVL(double p){
	double S=-0.354*log(p)+7.3385;
	return S;
}
//from 0.5 bar to 10 bar
double SsatVM (double p){
	double S=-0.336*log(p)+7.3599;
	return S;
}

//from 10 bar to 110 bar
double SsatVH (double p){
	double S=-1*pow(10,-6)*pow(p,3)+0.0003*pow(p,2)-0.0277*p+6.8057;
	return S;
}
double SsatV(double P){
	double S=0;
	if (P<=0.4)
	S=SsatVL(P);
	if (P>0.4 && P<=10)
	S=SsatVM(P);
	if (P>10 && P<110)
	S=SsatVH(P);
	return S;
	
}
//////////////////...functions of Cp...//////////////////
	double Cp(double T,double P){
	
    double Tstar=540; // Kelvin
    double Pstar=10; // bar
    double R=0.461526; // kJ/kg.K
    double Ti=T+273.15;
 

    double Tha=Tstar/Ti;
    double Pi=P/Pstar;

	////////...Region2...///////////
	
    ////////...Coefficients of Gammao...///////////
    int Jo[10];
    double no[10];
    
    Jo[0]=0;Jo[1]=0;Jo[2]=1;Jo[3]=-5;Jo[4]=-4;Jo[5]=-3;Jo[6]=-2;Jo[7]=-1;Jo[8]=2;Jo[9]=3;
    
    no[1]=-0.96927686500217*10;
    no[2]=0.10086655968018 *100;
    no[3]=-0.56087911283020/100;
    no[4]=0.71452738081455/10;
    no[5]=-0.40710498223928;
    no[6]=0.14240819171444*10;
    no[7]=-0.43839511319450*10;
    no[8]=-0.28408632460772;
    no[9]=0.21268463753307/10;

    
    ////////...Calculation...///////////
    double Gammaott=0;
    double Cp=0;
    double Gammaint=0;
    for (int i=1;i<10;i++){
    	Gammaint=no[i]*Jo[i]*(Jo[i]-1)*pow(Tha,Jo[i]-2);
    	Gammaott=Gammaott+Gammaint;
	}
 	////////...Coefficients of Gammar...///////////
 	int I[44];
    int J[44];
    double n[44];
    I[0]=0;I[1]=1;I[2]=1;I[3]=1;I[4]=1;I[5]=1;I[6]=2;I[7]=2;I[8]=2;
    I[9]=2;I[10]=2;I[11]=3;I[12]=3;I[13]=3;I[14]=3;
    I[15]=3;I[16]=4;I[17]=4;I[18]=4;I[19]=5;
    I[20]=6;I[21]=6;I[22]=6;
    I[23]=7;I[24]=7;I[25]=7;
    I[26]=8;I[27]=8;I[28]=9;
    I[29]=10;I[30]=10;I[31]=10;
    I[32]=16;I[33]=16;I[34]=18;
    I[35]=20;I[36]=20;I[37]=20;
    I[38]=21;I[39]=22;I[40]=23;
    I[41]=24;I[42]=24;I[43]=24;
    
    J[0]=0;J[1]=0;J[2]=1;J[3]=2;J[4]=3;J[5]=6;J[6]=1;J[7]=2;J[8]=4;
    J[9]=7;J[10]=36;J[11]=0;J[12]=1;J[13]=3;J[14]=6;J[15]=35;J[16]=1;
    J[17]=2;J[18]=3;J[19]=7;J[20]=3;J[21]=16;J[22]=35;J[23]=0;J[24]=11;
    J[25]=25;J[26]=8;J[27]=36;J[28]=13;J[29]=4;J[30]=10;J[31]=14;J[32]=29;J[33]=50;J[34]=57;
    J[35]=20;J[36]=35;J[37]=48;J[38]=21;J[39]=53;J[40]=39;J[41]=26;J[42]=40;J[43]=58;
    
    n[1]=-0.17731742473213/100;
    n[2]=-0.17834862292358/10;
    n[3]=-0.45996013696365/10;
    n[4]=-0.57581259083432/10;
    n[5]=-0.50325278727930/10;
    n[6]=-0.33032641670203*pow(10,-4);
    n[7]=-0.18948987516315/1000;
    n[8]=-0.39392777243355/100;
    n[9]=-0.43797295650573/10;
    n[10]=-0.26674547914087/10000;
    n[11]=0.20481737692309*pow(10,-7);
    n[12]=0.43870667284435*pow(10,-6);
    n[13]=-0.32277677238570*pow(10,-4);
    n[14]=-0.15033924542148/100;
    n[15]=-0.40668253562649/10;
    n[16]=-0.78847309559367*pow(10,-9);
    n[17]=0.12790717852285*pow(10,-7);
    n[18]=0.48225372718507 *pow(10,-6);
    n[19]=0.22922076337661 *pow(10,-5);
    n[20]=-0.16714766451061 *pow(10,-10);
    n[21]=-0.21171472321355 *pow(10,-2);
    n[22]=-0.23895741934104 *pow(10,2);
    n[23]=-0.59059564324270*pow(10,-17);
    n[24]=-0.12621808899101 *pow(10,-5);
    n[25]=-0.38946842435739 *pow(10,1);
    n[26]=0.11256211360459 *pow(10,-10);
    n[27]=-0.82311340897998 *pow(10,1);
    n[28]=0.19809712802088 *pow(10,-7);
    n[29]=0.10406965210174*pow(10,-18);
    n[30]=-0.10234747095929*pow(10,-12);
    n[31]=-0.10018179379511*pow(10,-8);
    n[32]=-0.80882908646985 *pow(10,-10);
    n[33]=0.10693031879409;
    n[34]=-0.33662250574171;
    n[35]=0.89185845355421*pow(10,-24);
    n[36]=0.30629316876232*pow(10,-12);
    n[37]=-0.42002467698208*pow(10,-7);
    n[38]=-0.59056029685639 *pow(10,-25);
    n[39]=0.37826947613457 *pow(10,-5);
    n[40]=-0.12768608934681 *pow(10,-14);
    n[41]=0.73087610595061*pow(10,-28);
    n[42]=0.55414715350778*pow(10,-16);
    n[43]=-0.94369707241210 *pow(10,-6);
	double Gammartt=0;
 	double Gammaint2=0;
    for (int i=1;i<44;i++){
    	Gammaint2=n[i]*pow(Pi,I[i])*J[i]*(J[i]-1)*pow((Tha-0.5),(J[i]-2));
    	Gammartt=Gammartt+Gammaint2;
	}
	Cp=-R*Tha*Tha*(Gammaott+Gammartt);
 	


    return Cp;
}
////////...Entropy of Superheated Steam Function...///////////

double Ssup(double T,double P){
	
    double Tstar=540; // Kelvin
    double Pstar=10; // bar
    double R=0.461526; // kJ/kg.K
    double Ti=T+273.15;
 

    double Tha=Tstar/Ti;
    double Pi=P/Pstar;

	////////...Region2...///////////
	
    ////////...Coefficients of Gammao...///////////
    int Jo[10];
    double no[10];
    
    Jo[0]=0;Jo[1]=0;Jo[2]=1;Jo[3]=-5;Jo[4]=-4;Jo[5]=-3;Jo[6]=-2;Jo[7]=-1;Jo[8]=2;Jo[9]=3;
    
    no[1]=-0.96927686500217*10;
    no[2]=0.10086655968018 *100;
    no[3]=-0.56087911283020/100;
    no[4]=0.71452738081455/10;
    no[5]=-0.40710498223928;
    no[6]=0.14240819171444*10;
    no[7]=-0.43839511319450*10;
    no[8]=-0.28408632460772;
    no[9]=0.21268463753307/10;

    
 	////////...Coefficients of Gammar...///////////
 	int I[44];
    int J[44];
    double n[44];
    I[0]=0;I[1]=1;I[2]=1;I[3]=1;I[4]=1;I[5]=1;I[6]=2;I[7]=2;I[8]=2;
    I[9]=2;I[10]=2;I[11]=3;I[12]=3;I[13]=3;I[14]=3;
    I[15]=3;I[16]=4;I[17]=4;I[18]=4;I[19]=5;
    I[20]=6;I[21]=6;I[22]=6;
    I[23]=7;I[24]=7;I[25]=7;
    I[26]=8;I[27]=8;I[28]=9;
    I[29]=10;I[30]=10;I[31]=10;
    I[32]=16;I[33]=16;I[34]=18;
    I[35]=20;I[36]=20;I[37]=20;
    I[38]=21;I[39]=22;I[40]=23;
    I[41]=24;I[42]=24;I[43]=24;
    
    J[0]=0;J[1]=0;J[2]=1;J[3]=2;J[4]=3;J[5]=6;J[6]=1;J[7]=2;J[8]=4;
    J[9]=7;J[10]=36;J[11]=0;J[12]=1;J[13]=3;J[14]=6;J[15]=35;J[16]=1;
    J[17]=2;J[18]=3;J[19]=7;J[20]=3;J[21]=16;J[22]=35;J[23]=0;J[24]=11;
    J[25]=25;J[26]=8;J[27]=36;J[28]=13;J[29]=4;J[30]=10;J[31]=14;J[32]=29;J[33]=50;J[34]=57;
    J[35]=20;J[36]=35;J[37]=48;J[38]=21;J[39]=53;J[40]=39;J[41]=26;J[42]=40;J[43]=58;
    
    n[1]=-0.17731742473213/100;
    n[2]=-0.17834862292358/10;
    n[3]=-0.45996013696365/10;
    n[4]=-0.57581259083432/10;
    n[5]=-0.50325278727930/10;
    n[6]=-0.33032641670203*pow(10,-4);
    n[7]=-0.18948987516315/1000;
    n[8]=-0.39392777243355/100;
    n[9]=-0.43797295650573/10;
    n[10]=-0.26674547914087/10000;
    n[11]=0.20481737692309*pow(10,-7);
    n[12]=0.43870667284435*pow(10,-6);
    n[13]=-0.32277677238570*pow(10,-4);
    n[14]=-0.15033924542148/100;
    n[15]=-0.40668253562649/10;
    n[16]=-0.78847309559367*pow(10,-9);
    n[17]=0.12790717852285*pow(10,-7);
    n[18]=0.48225372718507 *pow(10,-6);
    n[19]=0.22922076337661 *pow(10,-5);
    n[20]=-0.16714766451061 *pow(10,-10);
    n[21]=-0.21171472321355 *pow(10,-2);
    n[22]=-0.23895741934104 *pow(10,2);
    n[23]=-0.59059564324270*pow(10,-17);
    n[24]=-0.12621808899101 *pow(10,-5);
    n[25]=-0.38946842435739 *pow(10,1);
    n[26]=0.11256211360459 *pow(10,-10);
    n[27]=-0.82311340897998 *pow(10,1);
    n[28]=0.19809712802088 *pow(10,-7);
    n[29]=0.10406965210174*pow(10,-18);
    n[30]=-0.10234747095929*pow(10,-12);
    n[31]=-0.10018179379511*pow(10,-8);
    n[32]=-0.80882908646985 *pow(10,-10);
    n[33]=0.10693031879409;
    n[34]=-0.33662250574171;
    n[35]=0.89185845355421*pow(10,-24);
    n[36]=0.30629316876232*pow(10,-12);
    n[37]=-0.42002467698208*pow(10,-7);
    n[38]=-0.59056029685639 *pow(10,-25);
    n[39]=0.37826947613457 *pow(10,-5);
    n[40]=-0.12768608934681 *pow(10,-14);
    n[41]=0.73087610595061*pow(10,-28);
    n[42]=0.55414715350778*pow(10,-16);
    n[43]=-0.94369707241210 *pow(10,-6);
    
    ////////...Calculation...///////////
    double Gammaot=0;
    double Gammaotint=0;
    for (int i=1;i<10;i++){
    	Gammaotint=no[i]*Jo[i]*pow(Tha,Jo[i]-1);
    	Gammaot=Gammaot+Gammaotint;
	}
	double Gammao=0;
    double Gammaoint=0;
    for (int i=1;i<10;i++){
    	Gammaoint=no[i]*pow(Tha,Jo[i]);
    	Gammao=Gammao+Gammaoint;
	}
    Gammao=log(Pi)+Gammao;
    
	double Gammart=0;
 	double Gammartint=0;
    for (int i=1;i<44;i++){
    	Gammartint=n[i]*pow(Pi,I[i])*J[i]*pow((Tha-0.5),(J[i]-1));
    	Gammart=Gammart+Gammartint;
	}
	double Gammar=0;
 	double Gammarint=0;
    for (int i=1;i<44;i++){
    	Gammarint=n[i]*pow(Pi,I[i])*pow((Tha-0.5),(J[i]));
    	Gammar=Gammar+Gammarint;
	}
 	double Ssuper=R*Tha*(Gammaot+Gammart)-R*(Gammao+Gammar);


    return Ssuper;
}

////////...Enthalpy of Superheated Steam Function...///////////
	double Hsup(double T,double P){
	
    double Tstar=540; // Kelvin
    double Pstar=10; // bar
    double R=0.461526; // kJ/kg.K
    double Ti=T+273.15;
 

    double Tha=Tstar/Ti;
    double Pi=P/Pstar;

	////////...Region2...///////////
	
    ////////...Coefficients of Gammao...///////////
    int Jo[10];
    double no[10];
    
    Jo[0]=0;Jo[1]=0;Jo[2]=1;Jo[3]=-5;Jo[4]=-4;Jo[5]=-3;Jo[6]=-2;Jo[7]=-1;Jo[8]=2;Jo[9]=3;
    
    no[1]=-0.96927686500217*10;
    no[2]=0.10086655968018 *100;
    no[3]=-0.56087911283020/100;
    no[4]=0.71452738081455/10;
    no[5]=-0.40710498223928;
    no[6]=0.14240819171444*10;
    no[7]=-0.43839511319450*10;
    no[8]=-0.28408632460772;
    no[9]=0.21268463753307/10;


 	////////...Coefficients of Gammar...///////////
 	int I[44];
    int J[44];
    double n[44];
    I[0]=0;I[1]=1;I[2]=1;I[3]=1;I[4]=1;I[5]=1;I[6]=2;I[7]=2;I[8]=2;
    I[9]=2;I[10]=2;I[11]=3;I[12]=3;I[13]=3;I[14]=3;
    I[15]=3;I[16]=4;I[17]=4;I[18]=4;I[19]=5;
    I[20]=6;I[21]=6;I[22]=6;
    I[23]=7;I[24]=7;I[25]=7;
    I[26]=8;I[27]=8;I[28]=9;
    I[29]=10;I[30]=10;I[31]=10;
    I[32]=16;I[33]=16;I[34]=18;
    I[35]=20;I[36]=20;I[37]=20;
    I[38]=21;I[39]=22;I[40]=23;
    I[41]=24;I[42]=24;I[43]=24;
    
    J[0]=0;J[1]=0;J[2]=1;J[3]=2;J[4]=3;J[5]=6;J[6]=1;J[7]=2;J[8]=4;
    J[9]=7;J[10]=36;J[11]=0;J[12]=1;J[13]=3;J[14]=6;J[15]=35;J[16]=1;
    J[17]=2;J[18]=3;J[19]=7;J[20]=3;J[21]=16;J[22]=35;J[23]=0;J[24]=11;
    J[25]=25;J[26]=8;J[27]=36;J[28]=13;J[29]=4;J[30]=10;J[31]=14;J[32]=29;J[33]=50;J[34]=57;
    J[35]=20;J[36]=35;J[37]=48;J[38]=21;J[39]=53;J[40]=39;J[41]=26;J[42]=40;J[43]=58;
    
    n[1]=-0.17731742473213/100;
    n[2]=-0.17834862292358/10;
    n[3]=-0.45996013696365/10;
    n[4]=-0.57581259083432/10;
    n[5]=-0.50325278727930/10;
    n[6]=-0.33032641670203*pow(10,-4);
    n[7]=-0.18948987516315/1000;
    n[8]=-0.39392777243355/100;
    n[9]=-0.43797295650573/10;
    n[10]=-0.26674547914087/10000;
    n[11]=0.20481737692309*pow(10,-7);
    n[12]=0.43870667284435*pow(10,-6);
    n[13]=-0.32277677238570*pow(10,-4);
    n[14]=-0.15033924542148/100;
    n[15]=-0.40668253562649/10;
    n[16]=-0.78847309559367*pow(10,-9);
    n[17]=0.12790717852285*pow(10,-7);
    n[18]=0.48225372718507 *pow(10,-6);
    n[19]=0.22922076337661 *pow(10,-5);
    n[20]=-0.16714766451061 *pow(10,-10);
    n[21]=-0.21171472321355 *pow(10,-2);
    n[22]=-0.23895741934104 *pow(10,2);
    n[23]=-0.59059564324270*pow(10,-17);
    n[24]=-0.12621808899101 *pow(10,-5);
    n[25]=-0.38946842435739 *pow(10,1);
    n[26]=0.11256211360459 *pow(10,-10);
    n[27]=-0.82311340897998 *pow(10,1);
    n[28]=0.19809712802088 *pow(10,-7);
    n[29]=0.10406965210174*pow(10,-18);
    n[30]=-0.10234747095929*pow(10,-12);
    n[31]=-0.10018179379511*pow(10,-8);
    n[32]=-0.80882908646985 *pow(10,-10);
    n[33]=0.10693031879409;
    n[34]=-0.33662250574171;
    n[35]=0.89185845355421*pow(10,-24);
    n[36]=0.30629316876232*pow(10,-12);
    n[37]=-0.42002467698208*pow(10,-7);
    n[38]=-0.59056029685639 *pow(10,-25);
    n[39]=0.37826947613457 *pow(10,-5);
    n[40]=-0.12768608934681 *pow(10,-14);
    n[41]=0.73087610595061*pow(10,-28);
    n[42]=0.55414715350778*pow(10,-16);
    n[43]=-0.94369707241210 *pow(10,-6);

 	double Gammaot=0;
 	double Gammart=0;
 	double Gammaint3=0;
 	for (int i=1;i<10;i++){
    	Gammaint3=no[i]*Jo[i]*pow(Tha,(Jo[i]-1));
    	Gammaot=Gammaot+Gammaint3;
	}
	double Gammaint4=0;
	for (int i=1;i<44;i++){
    	Gammaint4=n[i]*pow(Pi,I[i])*J[i]*pow((Tha-0.5),(J[i]-1));
    	Gammart=Gammart+Gammaint4;
	}
	 double Hsuper=R*Ti*Tha*(Gammart+Gammaot);


    return Hsuper;
}
int main (){
	double Etha=0.865; // Polytropic Efficiency of Turbine Stages
	double EthaGen=0.98; // Generator Efficiency
	////////////...Inlet Parameters...///////////
	double Pin=88.26; // bar
	double Tin=535; // C
	double Hin=3475; // kj/kg
	double Qin=0, Qnom=206.67;
	cout<<"Please Enter Steam Flow (t/h)(Min:154, Max:235, THA:206.67):";
	cin>>Qin;
	double PerL=Qin/Qnom*100; // ton/hour
		
	////////////////////////////////////////////////...Gland Steam Parameters...//////////////////////////////////////////
   double ma,mb,mc,md,me,mf, mm;
   double ha,hb,hc,hd,he,hff, hm;
   ma=1.17; ha=3475;
   mb=0.13; hb=3475;
   mc=0.0179*Qin-0.156; hc=3415;
   md=0.0071*Qin+0.2087; hd=3415;
   me=0.002*Qin-0.0974; he=3415;
   mf=0.15; hff=3226;
   mm=0.44; hm=3300;
   
	//cout<<"----------\nInlet Parameters\nMin="<<Qin<<"\nPin="<<Pin<<"\nHin="<<Hin<<"\nTin="<<Tin<<endl;
	//cout<<"----------\nGland Steam Parameters:\nMa="<<ma<<"\nMb="<<mb<<"\nMc="<<mc<<"\nMd="<<md<<"\nMe="<<me<<"\nMf="<<mf<<"\nMm="<<mm<<endl;

///.............................................................Boiler....................................................................///

	/////////////////////////...fuel composition...///////////////////////////
	//Natural Gas Composition (Volume Fraction):
	double H2N=0;
	double CH4N=84;
	double C2H6N=5.35;
	double C3H8N=1.65;
	double C4H10N=0.3;
	double C5H12N=0.1;
	double CON=0.28;
	double CO2N=8;
	double O2N=0;
	double N2N=0.32;
	double H2ON=0;
	double NGHHV=40455; //kJ/kg
	double NGHHVV=35564; //kJ/Nm3
	double SumN=H2N+CH4N+C2H6N+C3H8N+C4H10N+C5H12N+CON+CO2N+O2N+N2N+H2ON;
	double NGEAR=5;
//	cout<<SumN;
	//Coke Oven Gas Composition (Volume Fraction):
	double H2c=56;
	double CH4c=18.8;
	double C2H6c=1.8;
	double C3H8c=0;
	double C4H10c=0;
	double C5H12c=0;
	double COc=6;
	double CO2c=8.6;
	double O2c=0.3;
	double N2c=8.5;
	double H2Oc=0;
	double COGHHV=26134; //kJ/kg
	double COGHHVV=14700; //kJ/Nm3
	double Sumc=H2c+CH4c+C2H6c+C3H8c+C4H10c+C5H12c+COc+CO2c+O2c+N2c+H2Oc;
	double COGEAR=7;
//	cout<<Sumc;
	//Blast Furnace Gas Composition (Volume Fraction):
	double H2b=6.5;
	double CH4b=0.3;
	double C2H6b=0;
	double C3H8b=0;
	double C4H10b=0;
	double C5H12b=0;
	double COb=24;
	double CO2b=17.6;
	double O2b=0.2;
	double N2b=51.4;
	double H2Ob=0;
	double BFGHHV=2961.85; //kJ/kg
	double BFGHHVV=3845; //kJ/Nm3
	double Sumb=H2b+CH4b+C2H6b+C3H8b+C4H10b+C5H12b+COb+CO2b+O2b+N2b+H2Ob;
	double BFEAR=10;
//	cout<<Sumb;

	////////////////////////...Molar Mass...////////////////////////
	double MN=14;
	double MO=16;
	double MC=12;
	double MH=1;
	double Ma=28.84;
	double MMf=0;

	double MMCH4=MC+4*MH;
	double MMC2H6=2*MC+4*MH;
	double MMC3H8=3*MC+8*MH;
	double MMC4H10=4*MC+10*MH;
	double MMC5H12=5*MC+12*MH;
	double MMH2=2*MH;
	double MMCO=MC+MO;
	double MMCO2=MC+2*MO;
	double MMN2=2*MN;
	double MMH2O=2*MH+MO;
	
	double MMNG=(H2N*MMH2+CH4N*MMCH4+C2H6N*MMC2H6+C3H8N*MMC3H8+C4H10N*MMC4H10+C5H12N*MMC5H12+CON*MMCO+CO2N*MMCO2+N2N*MMN2+H2ON*MMH2O)/100;
	//cout<<"Molar Mass of NG="<<MMNG<<" gr/mole";
	double MMCOG=(H2c*MMH2+CH4c*MMCH4+C2H6c*MMC2H6+C3H8c*MMC3H8+C4H10c*MMC4H10+C5H12c*MMC5H12+COc*MMCO+CO2c*MMCO2+N2c*MMN2+H2Oc*MMH2O)/100;
	//cout<<"Molar Mass of COG="<<MMCOG<<" gr/mole";
	double MMBFG=(H2b*MMH2+CH4b*MMCH4+C2H6b*MMC2H6+C3H8b*MMC3H8+C4H10b*MMC4H10+C5H12b*MMC5H12+COb*MMCO+CO2b*MMCO2+N2b*MMN2+H2Ob*MMH2O)/100;
	//cout<<"Molar Mass of BFGG="<<MMBFG<<" gr/mole";
	
	//Combined Fuel:
	double Tg=0;
	double Ta=0;
	double Phi1, Phi;
	Phi1=0;

	cout<<"Please Enter Ambient Air Temperature (C):";
	cin>>Ta;
	cout<<"Please Enter Relative Humidity (0-100):";
	cin>>Phi1;
	if (Phi1>100 || Phi1<0){
		cout<<"Please Choose a Number Between 0-100:";
		cin>>Phi1;
	}

	cout<<"\nPlease Enter Flue Gas Temperature (C):";
	cin>>Tg;
	Phi=Phi1/100;
	double NGSp=0, COGSp=0, BFGSp=0;
	cout<<"Please Enter share of N.G, C.O.G and B.F.G from Generated Power (Sum must be 100):";
	cout<<"\nNatural Gas Share of Boiler's Power is (0-100):";
	cin>>NGSp;
	cout<<"Cock Oven Gas Share of Boiler's Power is (0-100):";
	cin>>COGSp;
	cout<<"Blast Furnace Gas Share of Boiler's Power is (0-100):";
	cin>>BFGSp;

	//NGSp=50; COGSp=25; BFGSp=25;	
	double BTP=100000; //Boiler's Total Power (kwh)
	double BTPJ=BTP*3600; // Boiler's Total Power (kJ)
	double NGv=NGSp/100*BTPJ/NGHHVV; //Nm3
	double COGv=COGSp/100*BTPJ/COGHHVV; //Nm3
	double BFGv=BFGSp/100*BTPJ/BFGHHVV; //Nm3
	double totalv=NGv+COGv+BFGv; //Nm3
	double NGS1=0, COGS1=0, BFGS1=0;
	NGS1=NGv/totalv;
	COGS1=COGv/totalv;
	BFGS1=BFGv/totalv;
	double FuelMole=totalv/22.414;
	
	double NGS=NGS1, COGS=COGS1, BFGS=BFGS1;
	int NGSi=NGS*100, BFGSi=BFGS*100, COGSi=COGS*100;
	cout<<"\nVolume Fraction of NG="<<NGSi<<" Nm3/100 Nm3 of Fuel"<<"\nVolume Fraction of COG="<<COGSi<<" Nm3/100 Nm3 of Fuel"<<"\nVolume Fraction of BFG="<<BFGSi<<" Nm3/100 Nm3 of Fuel";
	double EAR=NGS*NGEAR+COGS*COGEAR+BFGS*BFEAR; //Minimum Required Excess Air
	double H2=NGS*H2N+COGS*H2c+H2b*BFGS; //Mole of H2 in Fuel
	double CH4=NGS*CH4N+COGS*CH4c+CH4b*BFGS;
	double C2H6=NGS*C2H6N+COGS*C2H6c+C2H6b*BFGS;
	double C3H8=NGS*C3H8N+COGS*C3H8c+C3H8b*BFGS;
	double C4H10=NGS*C4H10N+COGS*C4H10c+C4H10b*BFGS;
	double C5H12=NGS*C5H12N+COGS*C5H12c+C5H12b*BFGS;
	double CO=NGS*CON+COGS*COc+COb*BFGS;
	double CO2=NGS*CO2N+COGS*CO2c+CO2b*BFGS;
	double O2=NGS*O2N+COGS*O2c+O2b*BFGS;
	double N2=NGS*N2N+COGS*N2c+N2b*BFGS;
	double H2O=NGS*H2ON+COGS*H2Oc+H2Ob*BFGS;
	double FuelHHV=NGS*NGHHVV+COGS*COGHHVV+BFGS*BFGHHVV; //kJ/Nm3
	double Sum=H2+CH4+C2H6+C3H8+C4H10+C5H12+CO+CO2+O2+N2+H2O;
	//cout<<Sum;
	
	MMf=(H2*2*MH+CH4*MC+CH4*4*MH+C2H6*2*MC+C2H6*6*MH+C3H8*3*MC+C3H8*8*MH+C4H10*4*MC+C4H10*10*MH+C5H12*5*MC+C5H12*12*MH+CO*MC+CO*MO+CO2*MC+CO2*2*MO+O2*2*MO+N2*2*MN+H2O*2*MH+H2O*MO)/100;
	
	////////////////////////...Equation...////////////////////////
	double x=0; // CO2 Coefficient
	double y=0; // H2O Coefficient
	double z=0; // N2 Coefficient
	double z1=0; // N Coefficient in Lefthand Side
	double t1=0; // H2O Coefficient in Lefthand Side
	double w=0; // O2 Coefficient
	double v=0; // CO Coefficient
	double w1=0; // O Coefficient in lefthand Side
	double ath=0; //Stoichiometric Air Coefficient
	double EA1=0; // Excess Air Percent

	
	//Stoichiometric Equation:
	w1=CO+2*CO2+2*O2;
	z1=N2;
	t1=H2O;
	x=CH4+2*C2H6+3*C3H8+4*C4H10+5*C5H12+CO+CO2;
	y=(2*H2+4*CH4+6*C2H6+8*C3H8+10*C4H10+12*C5H12)/2+t1;
	ath=(2*x+y-w1)/2;
	z=N2+ath*3.76;
	cout<<"\n\nStoichiometric Equation:\n"<<H2<<" H2+"<<CH4<<" CH4+"<<C2H6<<" C2H6+"<<C3H8<<" C3H8+"<<C4H10<<" C4H10+"<<C5H12<<" C5H12+"<<CO<<" CO+\n"<<CO2<<" CO2+"<<O2<<" O2+"<<N2<<" N2+"<<t1<<" H2O+"<<ath<<"(O2+3.76 N2) --->\n"<<x<<" CO2+"<<y<<" H2O+"<<z<<" N2";
	
	//Excess Air Situation:
	cout<<"\n\n\nMinimum Excess Air Required is:"<<EAR<<" %";
	cout<<"\nPlease Enter Excess Air Percentage (0-100):";
	cin>>EA1;
	double EA=EA1-EAR;
	
	if (EA>=0){
	z=3.76*(ath*(1+EA/100))+z1;
	w=EA/100*ath;
	cout<<"\n\nCombustion Equation:\n"<<H2<<" H2+"<<CH4<<" CH4+"<<C2H6<<" C2H6+"<<C3H8<<" C3H8+"<<C4H10<<" C4H10+"<<C5H12<<" C5H12+"<<CO<<" CO+\n"<<CO2<<" CO2+"<<O2<<" O2+"<<N2<<" N2+"<<t1<<" H2O+"<<(1+EA/100)<<" X "<<ath<<"(O2+3.76 N2) --->\n"<<x<<" CO2+"<<y<<" H2O+"<<z<<" N2+"<<w<<" O2";
	}
	
	if (EA<0){
	y=(2*H2+4*CH4+6*C2H6+8*C3H8+10*C4H10+12*C5H12)/2;
	x=2*ath*(1+EA/100)+w1-(CH4+2*C2H6+3*C3H8+4*C4H10+5*C5H12+CO+CO2)-y;
	v=CH4+2*C2H6+3*C3H8+4*C4H10+5*C5H12+CO+CO2-x;
	z=3.76*(ath*(1+EA/100))+z1;
	if (x<0)
	{x=0;
	v=CH4+2*C2H6+3*C3H8+4*C4H10+5*C5H12+CO+CO2;}
	//cout<<"y="<<y<<"\n x="<<x<<"\n v="<<v<<"\n z="<<z<<"\n w="<<w;
	cout<<"\n\nCombustion Equation:\n"<<H2<<" H2+"<<CH4<<" CH4+"<<C2H6<<" C2H6+"<<C3H8<<" C3H8+"<<C4H10<<" C4H10+"<<C5H12<<" C5H12+"<<CO<<" CO+\n"<<CO2<<" CO2+"<<O2<<" O2+"<<N2<<" N2+"<<t1<<" H2O+"<<(1+EA/100)<<" X "<<ath<<"(O2+3.76 N2) --->\n"<<x<<" CO2+"<<v<<" CO+"<<y<<" H2O+"<<z<<" N2";
	}

	////////////////////////...Density...////////////////////////
	double RhoCH4=MMCH4/22.414;
	double RhoC2H6=MMC2H6/22.414;
	double RhoC3H8=MMC3H8/22.414;
	double RhoC4H10=MMC4H10/22.414;
	double RhoC5H12=MMC5H12/22.414;
	double RhoCO=MMCO/22.414;
	double RhoCO2=MMCO2/22.414;
	double RhoH2=MMH2/22.414;
	double RhoN2=MMN2/22.414;
	double RhoH2O=MMH2O/22.414;
	
	////////////////////////...Thermal Properties...////////////////////////
	double HHVV=NGS*NGHHVV+BFGS*BFGHHVV+COGS*COGHHVV; //kJ/kg
	double HHV=HHVV/MMf*22.414;
	double Cpfg=1.068;
	cout<<"\n\n\nMoles of Fuel= "<<FuelMole;
	cout<<"\nFuel LHV= "<<HHV;
	cout<<"\nMolar Mass of Fuel="<<MMf;
	cout<<"\nDensity of Fuel="<<MMf/22.414;
	////////////////////////...Losses...////////////////////////
	cout<<"\n\n......Boiler Losses......\n\n";
	///Stack Losses:
	cout<<"Stack Losses:";
	//Dry Gas Loss
	double Wg=(w*MO*2+z*MN*2+x*(MC+2*MO)+v*(MC+MO))/1000;
	double DGL=(Wg*(Tg-Ta)*Cpfg*100)/(HHV*MMf*100/1000);
	double dlge=DGL*(HHV*MMf*100/1000)/4.179;
	//cout<<"\n  Dry Gas Loss:"<<dlge<<" kJ";
	cout<<"\n  Dry Gas Loss:"<<DGL<<" %";
	//Moisture Loss
	double Hs=2646;
	double Ha=88;
	double ML=(9*H2*(2*MH)+H2O*(2*MH+MO))/(MMf*100)*(Hs-Ha)/HHV*100;
	cout<<"\n  Moisture Loss:"<<ML<<" %";
	//Humidity Loss
	double Pg=2;
	double Pt=84;
	double wh=0.622*Phi*Pg/(Pt-Phi*Pg);
	double HL=wh*(1+EA/100)*ath*Ma/1000*4.76*(Hs-Ha)/HHV/(MMf/1000*100)*100;
	cout<<"\n  Humidity Loss:"<<HL<<" %";
	//Radiation Loss
	double RL=(pow(10,-7)*Qin*Qin-6*pow(10,-5)*Qin+0.012)*100;
	cout<<"\n\n  Radiation Loss:"<<RL<<" %";
	//Uncountable Loss
	double UL=0.25;
	cout<<"\n  Uncountable Loss:"<<UL<<" %";
	//Manufacturer's Margine
	double MM=0.25;
	cout<<"\n  Manufacturer's Margine:"<<MM<<" %";
	//Incomplete Combustion Loss
	double ICL=0;
	double hfd=0;
	double COF=v/(x+v);
	//cout<<"\n      COF="<<COF;
	double CarbonWP=0;
	double FuelW=CH4*RhoCH4+C2H6*RhoC2H6+C3H8*RhoC3H8+C4H10*RhoC4H10+C5H12*RhoC5H12+CO*RhoCO+CO2*RhoCO2;
	double CH4WP=CH4*RhoCH4/FuelW;
	double C2H6WP=C2H6*RhoC2H6/FuelW;
	double C3H8WP=C3H8*RhoC3H8/FuelW;
	double C4H10WP=C4H10*RhoC4H10/FuelW;
	double C5H12WP=C5H12*RhoC5H12/FuelW;
	double COWP=CO*RhoCO/FuelW;
	double CO2WP=CO2*RhoCO2/FuelW;
	//double SumWP=CH4WP+C2H6WP+C3H8WP+C4H10WP+C5H12WP+COWP+CO2WP;
	//cout<<"\n      Ch4wp="<<CH4WP;
	CarbonWP=MC/MMCH4*CH4WP+2*MC/MMC2H6*C2H6WP+3*MC/MMC3H8*C3H8WP+4*MC/MMC4H10*C4H10WP+5*MC/MMC5H12*C5H12WP+MC/MMCO*COWP+MC/MMCO2*CO2WP;
	//cout<<"\n      CarbonWP="<<CarbonWP;
	hfd=10160*2.326;
	ICL=COF*CarbonWP*hfd/HHV/(MMf/1000*100)*100;
	cout<<"\n  Incomplete Combustion Loss:"<<ICL<<" %";
	/// Sigma Loss
	double SL=DGL+ML+HL+RL+UL+MM+ICL;
	cout<<"\n\nSum of Boiler Losses="<<SL<<" %";
	cout<<"\n\nEfficiency of Boiler="<<100-SL<<" %";
	double Ethab=100-SL;
	//cout<<"\n\nHeat Generated= "<<HHV*MMf*FuelMole/3600<<" kWh";
	
///......................................................Turbine..............................................................///
	//////...Stage 6...///
	double EthaPoly=Etha;
	if (PerL<100)
	EthaPoly=0.0014*Qin+0.55;
	double ER6=0.0007*Qin+0.6622;
	double Ps6=Pin*pow(ER6,6);
	double Cp6=Cp(Tin,Ps6);
	double Cv6=Cp6-8.3145/18.0153;
	double k6=Cp6/Cv6;
	double Ts6=(Tin+273.15)*pow((Ps6/Pin),(k6-1)/k6*EthaPoly)-273.15;
	double Hs6=Hsup(Ts6,Ps6);
	double Ss6=Ssup(Ts6,Ps6);
	double x6=0;
	if (Ss6>SsatV(Ps6))
	 x6=1;
	if (Ss6<SsatV(Ps6)){
	x6=(Ss6-SsatL(Ps6))/(SsatV(Ps6)-SsatL(Ps6));
	}
	//cout<<"----------\nStage6\nPs6="<<Ps6<<"\nTs6="<<Ts6<<"\nHs6="<<Hs6<<"\nSs6="<<Ss6<<"\nX6="<<x6;
	//////...Stage 8...///
	Etha=0.855;
	double ER8=0.8;
	double Ps8=Ps6*pow(ER8,2);
	double Cp8=Cp(Ts6,Ps8);
	double Cv8=Cp8-8.3145/18.0153;
	double k8=Cp8/Cv8;
	double Ts8=(Ts6+273.15)*pow((Ps8/Ps6),(k8-1)/k8*Etha)-273.15;
	double Hs8=Hsup(Ts8,Ps8);
	double Ss8=1.0*Ss6;
	double x8=0;
	if (Ss8>SsatV(Ps8))
	 x8=1;
	if (Ss8<SsatV(Ps8)){
	x8=(Ss8-SsatL(Ps8))/(SsatV(Ps8)-SsatL(Ps8));
	}
	//cout<<"\n----------\nStage8\nPs8="<<Ps8<<"\nTs8="<<Ts8<<"\nHs8="<<Hs8<<"\nSs8="<<Ss8<<"\nX8="<<x8;
	//////...Stage 10...///
	Etha=0.86;
	double ER10=0.775;
	double Ps10=Ps8*pow(ER10,2);
	double Cp10=Cp(Ts8,Ps10);
	double Cv10=Cp10-8.3145/18.0153;
	double k10=Cp10/Cv10;
	double Ts10=(Ts8+273.15)*pow((Ps10/Ps8),(k10-1)/k10*Etha)-273.15;
	double Hs10=Hsup(Ts10,Ps10);
	double Ss10=1.0*Ss8;
	double x10=0;
	if (Ss10>SsatV(Ps10))
	 x10=1;
	if (Ss10<SsatV(Ps10)){
	x10=(Ss10-SsatL(Ps10))/(SsatV(Ps10)-SsatL(Ps10));
	}
	//cout<<"\n----------\nStage10\nPs10="<<Ps10<<"\nTs10="<<Ts10<<"\nHs10="<<Hs10<<"\nSs10="<<Ss10<<"\nX10="<<x10;
	//////...Stage 14...///
	Etha=0.91;
	double ER14=0.735;
	double Ps14=Ps10*pow(ER14,4);
	double Cp14=Cp(Ts10,Ps14);
	double Cv14=Cp14-8.3145/18.0153;
	double k14=Cp14/Cv14;
	double Ts14=(Ts10+273.15)*pow((Ps14/Ps10),(k14-1)/k14*Etha)-273.15;
	double Hs14=Hsup(Ts14,Ps14);
	double Ss14=1.0*Ss10;
	double x14=0;
	if (Ss14>SsatV(Ps14))
	 x14=1;
	if (Ss14<SsatV(Ps14)){
	x14=(Ss14-SsatL(Ps14))/(SsatV(Ps14)-SsatL(Ps14));
}
	//cout<<"\n----------\nStage14\nPs14="<<Ps14<<"\nTs14="<<Ts14<<"\nHs14="<<Hs14<<"\nSs14="<<Ss14<<"\nX14="<<x14;
	//////...Stage 16...///
	Etha=0.90;
	double ER16=0.65;
	double Ps16=Ps14*pow(ER16,2);
	double Cp16=Cp(Ts14,Ps16);
	double Cv16=Cp16-8.3145/18.0153;
	double k16=Cp16/Cv16;
	double Ts16=(Ts14+273.15)*pow((Ps16/Ps14),(k16-1)/k16*Etha)-273.15;
	if (Ts16<Tsat(Ps16))
	Ts16=Tsat(Ps16);
	double Hs16=Hsup(Ts16,Ps16);
	double Ss16=1.0*Ss14;
	double x16=0;
	if (Ss16>SsatV(Ps16))
	 x16=1;
	if (Ss16<SsatV(Ps16)){
	x16=(Ss16-SsatL(Ps16))/(SsatV(Ps16)-SsatL(Ps16));
}
	if(x16<1)
	Hs16=x16*(HsatV(Ps16)-HsatL(Ps16))+HsatL(Ps16);
	//cout<<"\n----------\nStage16\nPs16="<<Ps16<<"\nTs16="<<Ts16<<"\nHs16="<<Hs16<<"\nSs16="<<Ss16<<"\nX16="<<x16;
	//////...Stage 18...///
	double ER18=0.56;
	double Ps18=Ps16*pow(ER18,2);
	double Cp18=Cp(Ts16,Ps18);
	double Cv18=Cp18-8.3145/18.0153;
	double k18=Cp18/Cv18;
	double Ts18=(Ts16+273.15)*pow((Ps18/Ps16),(k18-1)/k18*Etha)-273.15;
	double Hs18=Hsup(Ts18,Ps18);
	if (Ts18<Tsat(Ps18)){
	Ts18=Tsat(Ps18);
}
	double Ss18=1.0*Ss16;
	double x18=0;
	if (Ss18>SsatV(Ps18))
	 x18=1;
	if (Ss18<SsatV(Ps18)){
	x18=(Ss18-SsatL(Ps18))/(SsatV(Ps18)-SsatL(Ps18));
}
	if(x18<1)
	Hs18=x18*(HsatV(Ps18)-HsatL(Ps18))+HsatL(Ps18);

	//cout<<"\n----------\nStage18\nPs18="<<Ps18<<"\nTs18="<<Ts18<<"\nHs18="<<Hs18<<"\nSs18="<<Ss18<<"\nX18="<<x18;

	//////...Stage 20...///
	double ER20=-0.0005*Qin+0.620;
	double Ps20=Ps18*pow(ER20,2);
	double Cp20=Cp(Ts18,Ps20);
	double Cv20=Cp20-8.3145/18.0153;
	double k20=Cp20/Cv20;
	double Ts20=(Ts18+273.15)*pow((Ps20/Ps18),(k20-1)/k20*Etha)-273.15;
	double Hs20=Hsup(Ts20,Ps20);
	if (Ts20<Tsat(Ps20)){
	Ts20=Tsat(Ps20);
}
	
	double Ss20, x20;
	Ss20=1.01*Ss18;
	x20=(Ss20-SsatL(Ps20))/(SsatV(Ps20)-SsatL(Ps20));
	if(x20<1)
	Hs20=x20*(HsatV(Ps20)-HsatL(Ps20))+HsatL(Ps20);
	//cout<<"\n----------\nStage20\nPs20="<<Ps20<<"\nTs20="<<Ts20<<"\nHs20="<<Hs20<<"\nSs20="<<Ss20<<"\nX20="<<x20;

  
  ////////////////////////////////////////////////...Turbine Parameters...//////////////////////////////////////////
	double m1=Qin;
	double m13=0,m14=0,m15=0,m16=0,m17=0,m18=0;

	////////.... Calculating Extractions Flow Rate....//////////
	double h[26];
	for(int i=0; i<26; i++){
		h[i]=0;
	}
	h[13]=Hs6;h[14]=Hs8;h[15]=Hs10;h[16]=Hs14;h[17]=Hs16;h[18]=Hs18;h[3]=Hs20;
	h[2]=Hin;
	
	
	double m3=m1-m13-m14-m15-m16-m17-m18;
//	cout<<"\n----------\nFlow Rate of 6th Ex.="<<m18<<"\nFlow Rate of 5th Ex.="<<m17<<"\nFlow Rate of 4th Ex.="<<m16<<"\nFlow Rate of 3rd Ex.="<<m15<<"\nFlow Rate of 2nd Ex.="<<m14<<"\nFlow Rate of 1st Ex.="<<m13;

  ////////////////////////////////////////////////...Pumps...//////////////////////////////////////////
   double Powercp=190; //kW
   double Ethacp=0.7;
   
   double Powerfp=1600; //kW
   double Ethafp=0.7;

	////////////////////////////////////////////////...System of Equations Inputs ...//////////////////////////////////////////

	double hf[16];
	hf[1]=Hin;
	hf[2]=Hs6;
	hf[3]=Hs8;
	hf[4]=Hs10;
	hf[5]=Hs14;
	hf[6]=Hs16;
	hf[7]=Hs18;
	hf[8]=Hs20;
	hf[9]=HsatL(Ps20);
	double hde=681;
	hf[16]=HsatL(0.92*Ps6);
	hf[17]=HsatL(0.92*Ps8);
	hf[18]=HsatL(0.92*Ps14);
	hf[19]=HsatL(0.92*Ps16);
	hf[20]=HsatL(0.92*Ps18);
	
	double dh1h=0, dh2h=0, dh1l=0, dh2l=0, dh3l=0;
	dh3l=0.2174*Qin+44.698;
	dh2l=167.57*Qin-6133.2;
	if (dh2l<0)
	dh2l=0;
	dh1l=117.09*Qin-3084;
	if (dh1l<0)
	dh1l=0;
	dh2h=1.0088*Qin-66.015;
	if (dh2h<0)
	dh2h=0;
	dh1h=0.1909*Qin+59.613;

	double b[N][N+1];
	
	//Energy Conservation of Turbine
	b[0][0]=3.6;
	b[0][1]=1;
	b[0][2]=1;
	b[0][3]=1;
	b[0][4]=1;
	b[0][5]=1;
	b[0][6]=1;
	b[0][7]=1;
	b[0][8]=0;
	b[0][9]=0;
	b[0][10]=0;
	b[0][11]=0;
	b[0][12]=0;
	b[0][13]=0;
	b[0][14]=0;
	b[0][15]=0;
	b[0][16]=0;
	b[0][17]=0;
	b[0][18]=0;
	b[0][19]=0;
	b[0][20]=m1*Hin-1.03*(ma+mb+mc+md+me+mf+mm)*Hin;
	
	//Energy Conservation of Condensser
	b[1][0]=0;
	b[1][1]=0;
	b[1][2]=0;
	b[1][3]=0;
	b[1][4]=0;
	b[1][5]=0;
	b[1][6]=0;
	b[1][7]=-1/hf[8]*hf[9];
	b[1][8]=1;
	b[1][9]=0;
	b[1][10]=0;
	b[1][11]=0;
	b[1][12]=0;
	b[1][13]=0;
	b[1][14]=0;	
	b[1][15]=0;
	b[1][16]=0;
	b[1][17]=0;
	b[1][18]=0;
	b[1][19]=0;
	b[1][20]=0;
		
	//Energy Conservation of Condensate Pump
	b[2][0]=0;
	b[2][1]=0;
	b[2][2]=0;
	b[2][3]=0;
	b[2][4]=0;
	b[2][5]=0;
	b[2][6]=0;
	b[2][7]=0;
	b[2][8]=1;
	b[2][9]=-1;
	b[2][10]=0;
	b[2][11]=0;
	b[2][12]=0;
	b[2][13]=0;
	b[2][14]=0;
	b[2][15]=0;
	b[2][16]=0;
	b[2][17]=0;
	b[2][18]=0;
	b[2][19]=0;
	b[2][20]=-3.6*Ethacp*Powercp;
		
	//Energy Conservation of 3rd Lp Heater
	b[3][0]=0;
	b[3][1]=0;
	b[3][2]=0;
	b[3][3]=0;
	b[3][4]=0;
	b[3][5]=0;
	b[3][6]=1;
	b[3][7]=0;
	b[3][8]=0;
	b[3][9]=1;
	b[3][10]=-1;
	b[3][11]=0;
	b[3][12]=0;
	b[3][13]=0;
	b[3][14]=0;
	b[3][15]=0;
	b[3][16]=0;
	b[3][17]=0;
	b[3][18]=1;
	b[3][19]=-1;
	b[3][20]=-me*he-0.31*3226-0.44*3475;
	
	//Energy Conservation of 2nd Lp Heater
	b[4][0]=0;
	b[4][1]=0;
	b[4][2]=0;
	b[4][3]=0;
	b[4][4]=0;
	b[4][5]=1;
	b[4][6]=0;
	b[4][7]=0;
	b[4][8]=0;
	b[4][9]=0;
	b[4][10]=1;
	b[4][11]=-1;
	b[4][12]=0;
	b[4][13]=0;
	b[4][14]=0;
	b[4][15]=0;
	b[4][16]=0;
	b[4][17]=1;
	b[4][18]=-1;
	b[4][19]=1;
	b[4][20]=0;
	
	//Energy Conservation of 1st Lp Heater
	b[5][0]=0;
	b[5][1]=0;
	b[5][2]=0;
	b[5][3]=0;
	b[5][4]=1;
	b[5][5]=0;
	b[5][6]=0;
	b[5][7]=0;
	b[5][8]=0;
	b[5][9]=0;
	b[5][10]=0;
	b[5][11]=1;
	b[5][12]=-1;
	b[5][13]=0;
	b[5][14]=0;
	b[5][15]=0;
	b[5][16]=0;
	b[5][17]=-1;
	b[5][18]=0;
	b[5][19]=0;
	b[5][20]=-md*hd;
	
	//Energy Conservation of Deaerator
	b[6][0]=0;
	b[6][1]=0;
	b[6][2]=0;
	b[6][3]=1;
	b[6][4]=0;
	b[6][5]=0;
	b[6][6]=0;
	b[6][7]=0;
	b[6][8]=0;
	b[6][9]=0;
	b[6][10]=0;
	b[6][11]=0;
	b[6][12]=1;
	b[6][13]=0;
	b[6][14]=0;
	b[6][15]=0;
	b[6][16]=1;
	b[6][17]=0;
	b[6][18]=0;
	b[6][19]=0;
	b[6][20]=Qin*hde-3.6*Ethafp*Powerfp-ha*ma;
	
	//Energy Conservation of 2nd HP Heater
	b[7][0]=0;
	b[7][1]=0;
	b[7][2]=-1;
	b[7][3]=0;
	b[7][4]=0;
	b[7][5]=0;
	b[7][6]=0;
	b[7][7]=0;
	b[7][8]=0;
	b[7][9]=0;
	b[7][10]=0;
	b[7][11]=0;
	b[7][12]=0;
	b[7][13]=1;
	b[7][14]=0;
	b[7][15]=-1;
	b[7][16]=1;
	b[7][17]=0;
	b[7][18]=0;
	b[7][19]=0;
	b[7][20]=Qin*hde+mc*hc;
	
	//Energy Conservation of 1st Hp Heater
	b[8][0]=0;
	b[8][1]=1;
	b[8][2]=0;
	b[8][3]=0;
	b[8][4]=0;
	b[8][5]=0;
	b[8][6]=0;
	b[8][7]=0;
	b[8][8]=0;
	b[8][9]=0;
	b[8][10]=0;
	b[8][11]=0;
	b[8][12]=0;
	b[8][13]=1;
	b[8][14]=-1;
	b[8][15]=-1;
	b[8][16]=0;
	b[8][17]=0;
	b[8][18]=0;
	b[8][19]=0;
	b[8][20]=0;
		
	//Mass Conservation of Turbine
	b[9][0]=0;
	b[9][1]=1/Hs6;
	b[9][2]=1/Hs8;
	b[9][3]=1/Hs10;
	b[9][4]=1/Hs14;
	b[9][5]=1/Hs16;
	b[9][6]=1/Hs18;
	b[9][7]=1/Hs20;
	b[9][8]=0;
	b[9][9]=0;
	b[9][10]=0;
	b[9][11]=0;
	b[9][12]=0;
	b[9][13]=0;
	b[9][14]=0;
	b[9][15]=0;
	b[9][16]=0;
	b[9][17]=0;
	b[9][18]=0;
	b[9][19]=0;
	b[9][20]=Qin-ma-mb-mc-md-me-mf;
	
	//Enthalpy Increase in 3rd Lp Heater
	b[10][0]=0;
	b[10][1]=0;
	b[10][2]=0;
	b[10][3]=0;
	b[10][4]=0;
	b[10][5]=0;
	b[10][6]=0;
	b[10][7]=-dh3l/hf[8];
	b[10][8]=0;
	b[10][9]=-1;
	b[10][10]=1;
	b[10][11]=0;
	b[10][12]=0;
	b[10][13]=0;
	b[10][14]=0;
	b[10][15]=0;
	b[10][16]=0;
	b[10][17]=0;
	b[10][18]=0;
	b[10][19]=0;
	b[10][20]=(0.75+0.28)*dh3l;
	
	//Enthalpy Increase in 2nd Lp Heater
	b[11][0]=0;
	b[11][1]=0;
	b[11][2]=0;
	b[11][3]=0;
	b[11][4]=0;
	b[11][5]=0;
	b[11][6]=0;
	b[11][7]=0;
	b[11][8]=0;
	b[11][9]=0;
	b[11][10]=-1;
	b[11][11]=1;
	b[11][12]=0;
	b[11][13]=0;
	b[11][14]=0;
	b[11][15]=0;
	b[11][16]=0;
	b[11][17]=0;
	b[11][18]=0;
	b[11][19]=0;
	b[11][20]=dh2l;
	
	//Enthalpy Increase in 1st Lp Heater
	b[12][0]=0;
	b[12][1]=0;
	b[12][2]=0;
	b[12][3]=0;
	b[12][4]=0;
	b[12][5]=0;
	b[12][6]=0;
	b[12][7]=0;
	b[12][8]=0;
	b[12][9]=0;
	b[12][10]=0;
	b[12][11]=-1;
	b[12][12]=1;
	b[12][13]=0;
	b[12][14]=0;
	b[12][15]=0;
	b[12][16]=0;
	b[12][17]=0;
	b[12][18]=0;
	b[12][19]=0;
	b[12][20]=dh1l;

	//Enthalpy Increase in 2nd HP Heater
	b[13][0]=0;
	b[13][1]=0;
	b[13][2]=0;
	b[13][3]=0;
	b[13][4]=0;
	b[13][5]=0;
	b[13][6]=0;
	b[13][7]=0;
	b[13][8]=0;
	b[13][9]=0;
	b[13][10]=0;
	b[13][11]=0;
	b[13][12]=0;
	b[13][13]=1;
	b[13][14]=0;
	b[13][15]=0;
	b[13][16]=0;
	b[13][17]=0;
	b[13][18]=0;
	b[13][19]=0;
	b[13][20]=dh2h*Qin+681*Qin;

	//Enthalpy Increase in 1st HP Heater
	b[14][0]=0;
	b[14][1]=0;
	b[14][2]=0;
	b[14][3]=0;
	b[14][4]=0;
	b[14][5]=0;
	b[14][6]=0;
	b[14][7]=0;
	b[14][8]=0;
	b[14][9]=0;
	b[14][10]=0;
	b[14][11]=0;
	b[14][12]=0;
	b[14][13]=-1;
	b[14][14]=1;
	b[14][15]=0;
	b[14][16]=0;
	b[14][17]=0;
	b[14][18]=0;
	b[14][19]=0;
	b[14][20]=dh1h*Qin;
	
	//Mass Conservation in 1st HP Heater
	b[15][0]=0;
	b[15][1]=hf[16]/hf[2];
	b[15][2]=0;
	b[15][3]=0;
	b[15][4]=0;
	b[15][5]=0;
	b[15][6]=0;
	b[15][7]=0;
	b[15][8]=0;
	b[15][9]=0;
	b[15][10]=0;
	b[15][11]=0;
	b[15][12]=0;
	b[15][13]=0;
	b[15][14]=0;
	b[15][15]=-1;
	b[15][16]=0;
	b[15][17]=0;
	b[15][18]=0;
	b[15][19]=0;
	b[15][20]=0;
	
	//Mass Conservation in 2nd HP Heater
	b[16][0]=0;
	b[16][1]=0;
	b[16][2]=hf[17]/hf[3];
	b[16][3]=0;
	b[16][4]=0;
	b[16][5]=0;
	b[16][6]=0;
	b[16][7]=0;
	b[16][8]=0;
	b[16][9]=0;
	b[16][10]=0;
	b[16][11]=0;
	b[16][12]=0;
	b[16][13]=0;
	b[16][14]=0;
	b[16][15]=hf[17]/hf[16];
	b[16][16]=-1;
	b[16][17]=0;
	b[16][18]=0;
	b[16][19]=0;
	b[16][20]=0;
	
	//Mass Conservation in 1st LP Heater
	b[17][0]=0;
	b[17][1]=0;
	b[17][2]=0;
	b[17][3]=0;
	b[17][4]=hf[18]/hf[5];
	b[17][5]=0;
	b[17][6]=0;
	b[17][7]=0;
	b[17][8]=0;
	b[17][9]=0;
	b[17][10]=0;
	b[17][11]=0;
	b[17][12]=0;
	b[17][13]=0;
	b[17][14]=0;
	b[17][15]=0;
	b[17][16]=0;
	b[17][17]=-1;
	b[17][18]=0;
	b[17][19]=0;
	b[17][20]=0;
	
	//Mass Conservation in 2nd LP Heater
	b[18][0]=0;
	b[18][1]=0;
	b[18][2]=0;
	b[18][3]=0;
	b[18][4]=0;
	b[18][5]=hf[19]/hf[6];
	b[18][6]=0;
	b[18][7]=0;
	b[18][8]=0;
	b[18][9]=0;
	b[18][10]=0;
	b[18][11]=0;
	b[18][12]=0;
	b[18][13]=0;
	b[18][14]=0;
	b[18][15]=0;
	b[18][16]=0;
	b[18][17]=hf[19]/hf[18];
	b[18][18]=-1;
	b[18][19]=0;
	b[18][20]=0;
	
	//Mass Conservation in 3rd LP Heater
	b[19][0]=0;
	b[19][1]=0;
	b[19][2]=0;
	b[19][3]=0;
	b[19][4]=0;
	b[19][5]=0;
	b[19][6]=hf[20]/hf[7];
	b[19][7]=0;
	b[19][8]=0;
	b[19][9]=0;
	b[19][10]=0;
	b[19][11]=0;
	b[19][12]=0;
	b[19][13]=0;
	b[19][14]=0;
	b[19][15]=0;
	b[19][16]=0;
	b[19][17]=0;
	b[19][18]=hf[20]/hf[19];
	b[19][19]=-1;
	b[19][20]=0;

	    double a[N][N+1] = {
        {b[0][0], b[0][1], b[0][2], b[0][3], b[0][4], b[0][5], b[0][6], b[0][7], b[0][8], b[0][9], b[0][10], b[0][11], b[0][12], b[0][13], b[0][14], b[0][15], b[0][16], b[0][17], b[0][18], b[0][19], b[0][20]},
        {b[1][0], b[1][1], b[1][2], b[1][3], b[1][4], b[1][5], b[1][6], b[1][7], b[1][8], b[1][9], b[1][10], b[1][11], b[1][12], b[1][13], b[1][14], b[1][15], b[1][16], b[1][17], b[1][18], b[1][19], b[1][20]},
        {b[2][0], b[2][1], b[2][2], b[2][3], b[2][4], b[2][5], b[2][6], b[2][7], b[2][8], b[2][9], b[2][10], b[2][11], b[2][12], b[2][13], b[2][14], b[2][15], b[2][16], b[2][17], b[2][18], b[2][19], b[2][20]},
        {b[3][0], b[3][1], b[3][2], b[3][3], b[3][4], b[3][5], b[3][6], b[3][7], b[3][8], b[3][9], b[3][10], b[3][11], b[3][12], b[3][13], b[3][14], b[3][15], b[3][16], b[3][17], b[3][18], b[3][19], b[3][20]},
        {b[4][0], b[4][1], b[4][2], b[4][3], b[4][4], b[4][5], b[4][6], b[4][7], b[4][8], b[4][9], b[4][10], b[4][11], b[4][12], b[4][13], b[4][14], b[4][15], b[4][16], b[4][17], b[4][18], b[4][19], b[4][20]},
        {b[5][0], b[5][1], b[5][2], b[5][3], b[5][4], b[5][5], b[5][6], b[5][7], b[5][8], b[5][9], b[5][10], b[5][11], b[5][12], b[5][13], b[5][14], b[5][15], b[5][16], b[5][17], b[5][18], b[5][19], b[5][20]},
        {b[6][0], b[6][1], b[6][2], b[6][3], b[6][4], b[6][5], b[6][6], b[6][7], b[6][8], b[6][9], b[6][10], b[6][11], b[6][12], b[6][13], b[6][14], b[6][15], b[6][16], b[6][17], b[6][18], b[6][19], b[6][20]},
        {b[7][0], b[7][1], b[7][2], b[7][3], b[7][4], b[7][5], b[7][6], b[7][7], b[7][8], b[7][9], b[7][10], b[7][11], b[7][12], b[7][13], b[7][14], b[7][15], b[7][16], b[7][17], b[7][18], b[7][19], b[7][20]},
        {b[8][0], b[8][1], b[8][2], b[8][3], b[8][4], b[8][5], b[8][6], b[8][7], b[8][8], b[8][9], b[8][10], b[8][11], b[8][12], b[8][13], b[8][14], b[8][15], b[8][16], b[8][17], b[8][18], b[8][19], b[8][20]},
        {b[9][0], b[9][1], b[9][2], b[9][3], b[9][4], b[9][5], b[9][6], b[9][7], b[9][8], b[9][9], b[9][10], b[9][11], b[9][12], b[9][13], b[9][14], b[9][15], b[9][16], b[9][17], b[9][18], b[9][19], b[9][20]},
        {b[10][0], b[10][1], b[10][2], b[10][3], b[10][4], b[10][5], b[10][6], b[10][7], b[10][8], b[10][9], b[10][10], b[10][11], b[10][12], b[10][13], b[10][14], b[10][15], b[10][16], b[10][17], b[10][18], b[10][19], b[10][20]},
        {b[11][0], b[11][1], b[11][2], b[11][3], b[11][4], b[11][5], b[11][6], b[11][7], b[11][8], b[11][9], b[11][10], b[11][11], b[11][12], b[11][13], b[11][14], b[11][15], b[11][16], b[11][17], b[11][18], b[11][19], b[11][20]},
        {b[12][0], b[12][1], b[12][2], b[12][3], b[12][4], b[12][5], b[12][6], b[12][7], b[12][8], b[12][9], b[12][10], b[12][11], b[12][12], b[12][13], b[12][14], b[12][15], b[12][16], b[12][17], b[12][18], b[12][19], b[12][20]},
        {b[13][0], b[13][1], b[13][2], b[13][3], b[13][4], b[13][5], b[13][6], b[13][7], b[13][8], b[13][9], b[13][10], b[13][11], b[13][12], b[13][13], b[13][14], b[13][15], b[13][16], b[13][17], b[13][18], b[13][19], b[13][20]},
        {b[14][0], b[14][1], b[14][2], b[14][3], b[14][4], b[14][5], b[14][6], b[14][7], b[14][8], b[14][9], b[14][10], b[14][11], b[14][12], b[14][13], b[14][14], b[14][15], b[14][16], b[14][17], b[14][18], b[14][19], b[14][20]},
        {b[15][0], b[15][1], b[15][2], b[15][3], b[15][4], b[15][5], b[15][6], b[15][7], b[15][8], b[15][9], b[15][10], b[15][11], b[15][12], b[15][13], b[15][14], b[15][15], b[15][16], b[15][17], b[15][18], b[15][19], b[15][20]},
        {b[16][0], b[16][1], b[16][2], b[16][3], b[16][4], b[16][5], b[16][6], b[16][7], b[16][8], b[16][9], b[16][10], b[16][11], b[16][12], b[16][13], b[16][14], b[16][15], b[16][16], b[16][17], b[16][18], b[16][19], b[16][20]},
        {b[17][0], b[17][1], b[17][2], b[17][3], b[17][4], b[17][5], b[17][6], b[17][7], b[17][8], b[17][9], b[17][10], b[17][11], b[17][12], b[17][13], b[17][14], b[17][15], b[17][16], b[17][17], b[17][18], b[17][19], b[17][20]},
        {b[18][0], b[18][1], b[18][2], b[18][3], b[18][4], b[18][5], b[18][6], b[18][7], b[18][8], b[18][9], b[18][10], b[18][11], b[18][12], b[18][13], b[18][14], b[18][15], b[18][16], b[18][17], b[18][18], b[18][19], b[18][20]},
        {b[19][0], b[19][1], b[19][2], b[19][3], b[19][4], b[19][5], b[19][6], b[19][7], b[19][8], b[19][9], b[19][10], b[19][11], b[19][12], b[19][13], b[19][14], b[19][15], b[19][16], b[19][17], b[19][18], b[19][19], b[19][20]},
      
        
    
	};

    gauss_jordan(a, N);
    
	double MakeupPercent=0;
	cout<<"\nPlease Enter Percentage of Make-up Water (0-100):";
	cin>>MakeupPercent;
	double PowerGen=EthaGen*a[0][N];
	double Flow1stHpEx=a[1][N]/hf[2];
	double Flow2ndHpEx=a[2][N]/hf[3];
	double FlowDeaEx=a[3][N]/hf[4];
	double Flow1stLpEx=a[4][N]/hf[5];
	double Flow2ndLpEx=a[5][N]/hf[6];
	double Flow3rdLpEx=a[6][N]/hf[7];
	double FlowCond=a[7][N]/hf[8];
	double FlowCondOut=(1.0+MakeupPercent/100)*FlowCond+mm+0.31+0.28;
	double EnthalpyCondOut=a[8][N]/FlowCondOut; 
	//double EnthalpyCondOut=hf[9];    
	//double FlowCondOut=a[8][N]/hf[9];
	double EnthalpyCP=a[9][N]/FlowCondOut;
	double Enthalpy3rdLP=a[10][N]/FlowCondOut;
	double FlowInt=FlowCondOut+Flow1stLpEx+Flow2ndLpEx+Flow3rdLpEx+me+md;
	double Enthalpy2ndLP=a[11][N]/FlowInt;
	double Enthalpy1stLP=a[12][N]/FlowInt;
	double Enthalpy2ndHP=a[13][N]/Qin;
	double Enthalpy1stHP=a[14][N]/Qin;
	
/*
    
	cout << "\n----------\nPower Generated = " << PowerGen <<"  kW"<< endl;
	
	cout << "\nFlow of 1st Hp FWH Ext = " << Flow1stHpEx <<"  t/h"<< endl;
	
	cout << "Flow of 2nd HP FWH Ext = " << Flow2ndHpEx <<"  t/h"<< endl;
	
	cout << "Flow of Deaerator Ext = " << FlowDeaEx <<"  t/h"<< endl;
	
	cout << "Flow of 1st LP FWH Ext = " << Flow1stLpEx <<"  t/h"<< endl;
	
	cout << "Flow of 2nd LP FWH Ext = " << Flow2ndLpEx <<"  t/h"<< endl;
	
	cout << "Flow of 3rd LP FWH Ext = " << Flow3rdLpEx <<"  t/h"<< endl;
	
	cout << "\nCondensser Flow = " << FlowCond <<"  t/h"<< endl;
	
	cout << "Condensser Flow + Extractions + Gland Steam = " << Flow1stHpEx+Flow2ndHpEx+FlowDeaEx+Flow1stLpEx+Flow2ndLpEx+Flow3rdLpEx+FlowCond+ma+mb+mc+md+me+mf+mm <<"  t/h"<< endl;
	
	cout << "Condensser Outlet Enthalpy = " << EnthalpyCondOut <<"  kj/kg"<< endl;

	cout << "Condensser Outlet Flow = " << FlowCondOut <<"  t/h"<< endl;

	cout << "Enthalpy at Condensate Pump Outlet = " << EnthalpyCP <<"  Kj/Kg"<< endl;
	
	cout << "\nEnthalpy at 3rd LP F.W.H Outlet = " << Enthalpy3rdLP <<"  Kj/Kg"<< endl;

	cout << "Enthalpy at 2nd LP F.W.H Outlet = " << Enthalpy2ndLP <<"  Kj/Kg"<< endl;
	
	cout << "Enthalpy at 1st LP F.W.H Outlet = " << Enthalpy1stLP <<"  Kj/Kg"<< endl;
	
	cout << "\nEnthalpy at Deaerator's Outlet = " << hde <<"  Kj/Kg"<< endl;
	
	cout << "\nEnthalpy at 2nd HP F.W.H Outlet = " << Enthalpy2ndHP <<"  Kj/Kg"<< endl;
	
	cout << "Enthalpy at 1st HP F.W.H Outlet = " << Enthalpy1stHP <<"  Kj/Kg"<< endl;
	
	cout << "\nCondenser Heat Load = " << FlowCond*(hf[8]-hf[9])/3.6/1000<<"  MW"<< endl;
	
	*/
	  ////////////////////////////////////////////////...Condenser...//////////////////////////////////////////	
  	double mcw=9650; //ton/hour
  	double Cpcw=4.1797; //kJ/kg.K
  	double Tcwout=42.62;
  	double Tcwin=32.84;
  	double DeltaTcw=(Tcwout-Tcwin); //C
  	
  	double m4=0; //Condenser Outlet
  	double Ethac=(-0.373*PerL*PerL*PerL+111.17*PerL*PerL-3470.7*PerL+442226)/1000000;
  
	double Qout=Ethac*mcw*Cpcw*DeltaTcw;
	//cout << "Heat Disipated By the Cooling System = " << Qout/3.6/1000/Ethac<<"  MW"<< endl;
	

	
	////////////////////////////////////////////////...Boiler Inputs ...//////////////////////////////////////////
	
	double QBinf=Qin*(Hin-Enthalpy1stHP);
	double QBina=QBinf/Ethab/3600;
	
	double EffCycle=PowerGen/1000/QBina;


	
	//cout <<"\n-----------------------\n";
	//cout << "Boiler Heat Input = " << QBina<<"  MW"<< endl;
	//cout << "Power Generated = " << PowerGen/1000 <<"  MW"<< endl;
	//cout << "Cycle Efficiency = " << EffCycle*100<<" %"<< endl;
	//cout << "Natural Gas Consumption = " << Qngt<<" Nm3/hr"<< endl;
	
	double Qcond=FlowCond*(hf[8]-hf[9]);
	double deltaTct=Qcond/9650/4.179*1.2;
	//cout << "Delta T Cooling Water = " << deltaTct<<" C"<< endl;

  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////...Second Iteration.../////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
 
	//cout<<"\n----------Second Iteration----------\n"<<endl;
	cout<<"Calculating...\n";
	//cout<<"Enter Cooling Water inlet Temperature (C):";
  	//cin>>Tcwin;
  	double Twet=8*pow(10,-5)*Ta*Ta*Ta-0.0084*Ta*Ta+0.7427*Ta-1.2072+4;
  	double CTAP=8.34; //Cooling Tower Approach
	Tcwin=Twet+CTAP;
  	double Sin=Ssup(Tin,Pin);
	cout<<"----------\nInlet Parameters\nMin="<<Qin<<"\nPin="<<Pin<<"\nHin="<<Hin<<"\nTin="<<Tin<<"\nSin="<<Ssup(Tin,Pin)<<endl;
	cout<<"----------\nGland Steam Parameters:\nMa="<<ma<<"\nMb="<<mb<<"\nMc="<<mc<<"\nMd="<<md<<"\nMe="<<me<<"\nMf="<<mf<<"\nMm="<<mm<<endl;
	  ////////////////////////////////////////////////...Condenser...//////////////////////////////////////////	
  	mcw=9650; //ton/hour
  	Cpcw=4.1797; //kJ/kg.K
  	
  	Ethac=(-0.373*PerL*PerL*PerL+111.17*PerL*PerL-3470.7*PerL+442226)/1000000;
	
	// Condenser Parameters
	double A=3800; //m2
	double U0=4.4;
	double Cf=0.85;
	double Cw=4*pow(10,-6)*Tcwin*Tcwin*Tcwin-0.0007*Tcwin*Tcwin+0.0355*Tcwin+0.5033;
	double U=Cw*Cf*U0;
	double NTU=A*U/Cpcw/mcw*3.6;
	double Coeff=mcw/3.6*Cpcw*(1-exp(-NTU));
    //Condenser Temperature
    double Coeff1=Coeff/FlowCond*3.6;
    double as=0.0036;
    double bs=2.0266+Coeff1;
    double cs=-2495.2-Coeff1*Tcwin;
    double deltas=bs*bs-4*as*cs;
    double Tcs=(-1*bs+sqrt(deltas))/2/as;
    
    cout << "\n----------\nCooling Parameters:\nDry Bulb Temperature = " << Ta<<"  (C)"<< endl;
    cout << "Wet Bulb Temperature = " << Twet<<"  (C)"<< endl;
    cout << "Cooling Water Temperature = " << Tcwin<<"  (C)"<< endl;
	cout << "Condenser Temperature = " << Tcs<<"  (C)"<< endl;
	
	//////...Stage 6...///
	EthaPoly=0.875;
	if (PerL<100)
	EthaPoly=0.0014*Qin+0.55;
	 ER6=0.0007*Qin+0.6622;
	 Ps6=Pin*pow(ER6,6);
	 Cp6=Cp(Ts6,Ps6);
	 Cv6=Cp6-8.3145/18.0153;
	 k6=Cp6/Cv6;
	 Ts6=(Tin+273.15)*pow((Ps6/Pin),(k6-1)/k6*EthaPoly)-273.15;
	 Hs6=Hsup(Ts6,Ps6);
	 Ss6=Ssup(Ts6,Ps6);
	 x6=0;
	if (Ss6>SsatV(Ps6))
	 x6=1;
	if (Ss6<SsatV(Ps6)){
	x6=(Ss6-SsatL(Ps6))/(SsatV(Ps6)-SsatL(Ps6));
	}
	cout<<"----------\nStage6\nPs6="<<Ps6<<"\nTs6="<<Ts6<<"\nHs6="<<Hs6<<"\nSs6="<<Ss6<<"\nX6="<<x6;
	//////...Stage 8...///
	Etha=0.855;
	 ER8=0.8;
	 Ps8=Ps6*pow(ER8,2);
	 Cp8=Cp(Ts8,Ps8);
	 Cv8=Cp8-8.3145/18.0153;
	 k8=Cp8/Cv8;
	 Ts8=(Ts6+273.15)*pow((Ps8/Ps6),(k8-1)/k8*Etha)-273.15;
	 Hs8=Hsup(Ts8,Ps8);
	 Ss8=1.01*Ss6;
	 x8=0;
	if (Ss8>SsatV(Ps8))
	 x8=1;
	if (Ss8<SsatV(Ps8)){
	x8=(Ss8-SsatL(Ps8))/(SsatV(Ps8)-SsatL(Ps8));
	}
	cout<<"\n----------\nStage8\nPs8="<<Ps8<<"\nTs8="<<Ts8<<"\nHs8="<<Hs8<<"\nSs8="<<Ss8<<"\nX8="<<x8;
	//////...Stage 10...///
	Etha=0.875;
	 ER10=0.775;
	 Ps10=Ps8*pow(ER10,2);
	 Cp10=Cp(Ts10,Ps10);
	 Cv10=Cp10-8.3145/18.0153;
	 k10=Cp10/Cv10;
	 Ts10=(Ts8+273.15)*pow((Ps10/Ps8),(k10-1)/k10*Etha)-273.15;
	 Hs10=Hsup(Ts10,Ps10);
	 Ss10=1.01*Ss8;
	 x10=0;
	if (Ss10>SsatV(Ps10))
	 x10=1;
	if (Ss10<SsatV(Ps10)){
	x10=(Ss10-SsatL(Ps10))/(SsatV(Ps10)-SsatL(Ps10));
	}
	cout<<"\n----------\nStage10\nPs10="<<Ps10<<"\nTs10="<<Ts10<<"\nHs10="<<Hs10<<"\nSs10="<<Ss10<<"\nX10="<<x10;
	//////...Stage 14...///
	Etha=1;
	 ER14=0.735;
	 Ps14=Ps10*pow(ER14,4);
	 Cp14=Cp(Ts14,Ps14);
	 Cv14=Cp14-8.3145/18.0153;
	 k14=Cp14/Cv14;
	 Ts14=(Ts10+273.15)*pow((Ps14/Ps10),(k14-1)/k14*Etha)-273.15;
	 Hs14=Hsup(Ts14,Ps14);
	 Ss14=1.01*Ss10;
	 x14=0;
	if (Ss14>SsatV(Ps14))
	 x14=1;
	if (Ss14<SsatV(Ps14)){
	x14=(Ss14-SsatL(Ps14))/(SsatV(Ps14)-SsatL(Ps14));
}
	cout<<"\n----------\nStage14\nPs14="<<Ps14<<"\nTs14="<<Ts14<<"\nHs14="<<Hs14<<"\nSs14="<<Ss14<<"\nX14="<<x14;
	//////...Stage 16...///
	Etha=1;
	 ER16=0.65;
	 Ps16=Ps14*pow(ER16,2);
	 Cp16=Cp(Ts16,Ps16);
	 Cv16=Cp16-8.3145/18.0153;
	 k16=Cp16/Cv16;
	 Ts16=(Ts14+273.15)*pow((Ps16/Ps14),(k16-1)/k16*Etha)-273.15;
	if (Ts16<Tsat(Ps16))
	Ts16=Tsat(Ps16);
	 Hs16=Hsup(Ts16,Ps16);
	 Ss16=1.01*Ss14;
	 x16=0;
	if (Ss16>SsatV(Ps16))
	 x16=1;
	if (Ss16<SsatV(Ps16)){
	x16=(Ss16-SsatL(Ps16))/(SsatV(Ps16)-SsatL(Ps16));
}
	if(x16<1)
	Hs16=x16*(HsatV(Ps16)-HsatL(Ps16))+HsatL(Ps16);
	cout<<"\n----------\nStage16\nPs16="<<Ps16<<"\nTs16="<<Ts16<<"\nHs16="<<Hs16<<"\nSs16="<<Ss16<<"\nX16="<<x16;
	//////...Stage 18...///
	 ER18=0.56;
	 Ps18=Ps16*pow(ER18,2);
	 Cp18=Cp(Ts18,Ps18);
	 Cv18=Cp18-8.3145/18.0153;
	 k18=Cp18/Cv18;
	 Ts18=(Ts16+273.15)*pow((Ps18/Ps16),(k18-1)/k18*Etha)-273.15;
	 Hs18=Hsup(Ts18,Ps18);
	if (Ts18<Tsat(Ps18)){
	Ts18=Tsat(Ps18);
}
	 Ss18=1.01*Ss16;
	 x18=0;
	if (Ss18>SsatV(Ps18))
	 x18=1;
	if (Ss18<SsatV(Ps18)){
	x18=(Ss18-SsatL(Ps18))/(SsatV(Ps18)-SsatL(Ps18));
}
	if(x18<1)
	Hs18=x18*(HsatV(Ps18)-HsatL(Ps18))+HsatL(Ps18);

	cout<<"\n----------\nStage18\nPs18="<<Ps18<<"\nTs18="<<Ts18<<"\nHs18="<<Hs18<<"\nSs18="<<Ss18<<"\nX18="<<x18;

	//////...Stage 20...///
	 ER20=0.195979039369815+0.0129231897692486*Tcs-0.00125398302006083*Qin;
	 Ps20=Ps18*pow(ER20,2);
	 Cp20=Cp(Ts20,Ps20);
	 Cv20=Cp20-8.3145/18.0153;
	 k20=Cp20/Cv20;
	 Ts20=(Ts18+273.15)*pow((Ps20/Ps18),(k20-1)/k20*Etha)-273.15;
	 Hs20=Hsup(Ts20,Ps20);
	if (Ts20<Tsat(Ps20)){
	Ts20=Tsat(Ps20);
}
	
	 Ss20, x20;
	Ss20=1.01*Ss18;
	x20=(Ss20-SsatL(Ps20))/(SsatV(Ps20)-SsatL(Ps20));
	if(x20<1)
	Hs20=x20*(HsatV(Ps20)-HsatL(Ps20))+HsatL(Ps20);
	cout<<"\n----------\nStage20\nPs20="<<Ps20<<"\nTs20="<<Ts20<<"\nHs20="<<Hs20<<"\nSs20="<<Ss20<<"\nX20="<<x20;
	double xp20, Hsp20, Ssp20;
	Ssp20=Sin;
	xp20=(Ssp20-SsatL(Ps20))/(SsatV(Ps20)-SsatL(Ps20));
	Hsp20=xp20*(HsatV(Ps20)-HsatL(Ps20))+HsatL(Ps20);
	
	//Turbine Isentropic Efficiency
	double TIE=(Hin-Hs20)/(Hin-Hsp20)*100;
	double TIEs=(Hin-Hs20)/((Hin-Hs20)+(273.15+Ts20)*(Ss20-Sin))*100;
	cout<<"\n----------\nTurbine Isentropic Efficiency\nHs="<<Hsp20<<"\nTurbine Isentropic Efficiency="<<TIE<<"\nTurbine Isentropic Efficiency (Entropy Formula)="<<TIEs;

  ////////////////////////////////////////////////...Turbine Parameters...//////////////////////////////////////////
	 m1=Qin;
	 m13=0,m14=0,m15=0,m16=0,m17=0,m18=0;

	////////.... Calculating Extractions Flow Rate....//////////

	h[13]=Hs6;h[14]=Hs8;h[15]=Hs10;h[16]=Hs14;h[17]=Hs16;h[18]=Hs18;h[3]=Hs20;
	h[2]=Hin;
	
	
	 m3=m1-m13-m14-m15-m16-m17-m18;
//	cout<<"\n----------\nFlow Rate of 6th Ex.="<<m18<<"\nFlow Rate of 5th Ex.="<<m17<<"\nFlow Rate of 4th Ex.="<<m16<<"\nFlow Rate of 3rd Ex.="<<m15<<"\nFlow Rate of 2nd Ex.="<<m14<<"\nFlow Rate of 1st Ex.="<<m13;
	
  ////////////////////////////////////////////////...Pumps...//////////////////////////////////////////
    Powercp=190; //kW
    Ethacp=0.7;
   
    Powerfp=1600; //kW
    Ethafp=0.7;

	////////////////////////////////////////////////...System of Equations Inputs ...//////////////////////////////////////////

	hf[1]=Hin;
	hf[2]=Hs6;
	hf[3]=Hs8;
	hf[4]=Hs10;
	hf[5]=Hs14;
	hf[6]=Hs16;
	hf[7]=Hs18;
	hf[8]=Hs20;
	hf[9]=HsatL(Ps20);
	 hde=681;
	hf[16]=HsatL(0.92*Ps6);
	hf[17]=HsatL(0.92*Ps8);
	hf[18]=HsatL(0.92*Ps14);
	hf[19]=HsatL(0.92*Ps16);
	hf[20]=HsatL(0.92*Ps18);
	

	
	dh1h=0, dh2h=0, dh1l=0, dh2l=0, dh3l=0;
	dh3l=0.2174*Qin+44.698;
	dh2l=167.57*Qin-6133.2;
	if (dh2l<0)
	dh2l=0;
	dh1l=117.09*Qin-3084;
	if (dh1l<0)
	dh1l=0;
	dh2h=1.0088*Qin-66.015;
	if (dh2h<0)
	dh2h=0;
	dh1h=0.1909*Qin+59.613;
	

	 b[N][N+1];
	 for (int i=0; i<N+1; i++)
	 for (int j=0; j<N+2; j++)
	 b[i][j]=0;
	
	//Energy Conservation of Turbine
	b[0][0]=3.6;
	b[0][1]=1;
	b[0][2]=1;
	b[0][3]=1;
	b[0][4]=1;
	b[0][5]=1;
	b[0][6]=1;
	b[0][7]=1;
	b[0][8]=0;
	b[0][9]=0;
	b[0][10]=0;
	b[0][11]=0;
	b[0][12]=0;
	b[0][13]=0;
	b[0][14]=0;
	b[0][15]=0;
	b[0][16]=0;
	b[0][17]=0;
	b[0][18]=0;
	b[0][19]=0;
	b[0][20]=m1*Hin-1.03*(ma+mb+mc+md+me+mf+mm)*Hin;
	
	//Energy Conservation of Condensser
	b[1][0]=0;
	b[1][1]=0;
	b[1][2]=0;
	b[1][3]=0;
	b[1][4]=0;
	b[1][5]=0;
	b[1][6]=0;
	b[1][7]=-1/hf[8]*hf[9];
	b[1][8]=1;
	b[1][9]=0;
	b[1][10]=0;
	b[1][11]=0;
	b[1][12]=0;
	b[1][13]=0;
	b[1][14]=0;	
	b[1][15]=0;
	b[1][16]=0;
	b[1][17]=0;
	b[1][18]=0;
	b[1][19]=0;
	b[1][20]=0;
		
	//Energy Conservation of Condensate Pump
	b[2][0]=0;
	b[2][1]=0;
	b[2][2]=0;
	b[2][3]=0;
	b[2][4]=0;
	b[2][5]=0;
	b[2][6]=0;
	b[2][7]=0;
	b[2][8]=1;
	b[2][9]=-1;
	b[2][10]=0;
	b[2][11]=0;
	b[2][12]=0;
	b[2][13]=0;
	b[2][14]=0;
	b[2][15]=0;
	b[2][16]=0;
	b[2][17]=0;
	b[2][18]=0;
	b[2][19]=0;
	b[2][20]=-3.6*Ethacp*Powercp;
		
	//Energy Conservation of 3rd Lp Heater
	b[3][0]=0;
	b[3][1]=0;
	b[3][2]=0;
	b[3][3]=0;
	b[3][4]=0;
	b[3][5]=0;
	b[3][6]=1;
	b[3][7]=0;
	b[3][8]=0;
	b[3][9]=1;
	b[3][10]=-1;
	b[3][11]=0;
	b[3][12]=0;
	b[3][13]=0;
	b[3][14]=0;
	b[3][15]=0;
	b[3][16]=0;
	b[3][17]=0;
	b[3][18]=1;
	b[3][19]=-1;
	b[3][20]=-me*he-0.31*3226-0.44*3475;
	
	//Energy Conservation of 2nd Lp Heater
	b[4][0]=0;
	b[4][1]=0;
	b[4][2]=0;
	b[4][3]=0;
	b[4][4]=0;
	b[4][5]=1;
	b[4][6]=0;
	b[4][7]=0;
	b[4][8]=0;
	b[4][9]=0;
	b[4][10]=1;
	b[4][11]=-1;
	b[4][12]=0;
	b[4][13]=0;
	b[4][14]=0;
	b[4][15]=0;
	b[4][16]=0;
	b[4][17]=1;
	b[4][18]=-1;
	b[4][19]=1;
	b[4][20]=0;
	
	//Energy Conservation of 1st Lp Heater
	b[5][0]=0;
	b[5][1]=0;
	b[5][2]=0;
	b[5][3]=0;
	b[5][4]=1;
	b[5][5]=0;
	b[5][6]=0;
	b[5][7]=0;
	b[5][8]=0;
	b[5][9]=0;
	b[5][10]=0;
	b[5][11]=1;
	b[5][12]=-1;
	b[5][13]=0;
	b[5][14]=0;
	b[5][15]=0;
	b[5][16]=0;
	b[5][17]=-1;
	b[5][18]=0;
	b[5][19]=0;
	b[5][20]=-md*hd;
	
	//Energy Conservation of Deaerator
	b[6][0]=0;
	b[6][1]=0;
	b[6][2]=0;
	b[6][3]=1;
	b[6][4]=0;
	b[6][5]=0;
	b[6][6]=0;
	b[6][7]=0;
	b[6][8]=0;
	b[6][9]=0;
	b[6][10]=0;
	b[6][11]=0;
	b[6][12]=1;
	b[6][13]=0;
	b[6][14]=0;
	b[6][15]=0;
	b[6][16]=1;
	b[6][17]=0;
	b[6][18]=0;
	b[6][19]=0;
	b[6][20]=Qin*hde-3.6*Ethafp*Powerfp-ha*ma;
	
	//Energy Conservation of 2nd HP Heater
	b[7][0]=0;
	b[7][1]=0;
	b[7][2]=-1;
	b[7][3]=0;
	b[7][4]=0;
	b[7][5]=0;
	b[7][6]=0;
	b[7][7]=0;
	b[7][8]=0;
	b[7][9]=0;
	b[7][10]=0;
	b[7][11]=0;
	b[7][12]=0;
	b[7][13]=1;
	b[7][14]=0;
	b[7][15]=-1;
	b[7][16]=1;
	b[7][17]=0;
	b[7][18]=0;
	b[7][19]=0;
	b[7][20]=Qin*hde+mc*hc;
	
	//Energy Conservation of 1st Hp Heater
	b[8][0]=0;
	b[8][1]=1;
	b[8][2]=0;
	b[8][3]=0;
	b[8][4]=0;
	b[8][5]=0;
	b[8][6]=0;
	b[8][7]=0;
	b[8][8]=0;
	b[8][9]=0;
	b[8][10]=0;
	b[8][11]=0;
	b[8][12]=0;
	b[8][13]=1;
	b[8][14]=-1;
	b[8][15]=-1;
	b[8][16]=0;
	b[8][17]=0;
	b[8][18]=0;
	b[8][19]=0;
	b[8][20]=0;
		
	//Mass Conservation of Turbine
	b[9][0]=0;
	b[9][1]=1/Hs6;
	b[9][2]=1/Hs8;
	b[9][3]=1/Hs10;
	b[9][4]=1/Hs14;
	b[9][5]=1/Hs16;
	b[9][6]=1/Hs18;
	b[9][7]=1/Hs20;
	b[9][8]=0;
	b[9][9]=0;
	b[9][10]=0;
	b[9][11]=0;
	b[9][12]=0;
	b[9][13]=0;
	b[9][14]=0;
	b[9][15]=0;
	b[9][16]=0;
	b[9][17]=0;
	b[9][18]=0;
	b[9][19]=0;
	b[9][20]=Qin-ma-mb-mc-md-me-mf;
	
	//Enthalpy Increase in 3rd Lp Heater
	b[10][0]=0;
	b[10][1]=0;
	b[10][2]=0;
	b[10][3]=0;
	b[10][4]=0;
	b[10][5]=0;
	b[10][6]=0;
	b[10][7]=-dh3l/hf[8];
	b[10][8]=0;
	b[10][9]=-1;
	b[10][10]=1;
	b[10][11]=0;
	b[10][12]=0;
	b[10][13]=0;
	b[10][14]=0;
	b[10][15]=0;
	b[10][16]=0;
	b[10][17]=0;
	b[10][18]=0;
	b[10][19]=0;
	b[10][20]=(0.75+0.28)*dh3l;
	
	//Enthalpy Increase in 2nd Lp Heater
	b[11][0]=0;
	b[11][1]=0;
	b[11][2]=0;
	b[11][3]=0;
	b[11][4]=0;
	b[11][5]=0;
	b[11][6]=0;
	b[11][7]=0;
	b[11][8]=0;
	b[11][9]=0;
	b[11][10]=-1;
	b[11][11]=1;
	b[11][12]=0;
	b[11][13]=0;
	b[11][14]=0;
	b[11][15]=0;
	b[11][16]=0;
	b[11][17]=0;
	b[11][18]=0;
	b[11][19]=0;
	b[11][20]=dh2l;
	
	//Enthalpy Increase in 1st Lp Heater
	b[12][0]=0;
	b[12][1]=0;
	b[12][2]=0;
	b[12][3]=0;
	b[12][4]=0;
	b[12][5]=0;
	b[12][6]=0;
	b[12][7]=0;
	b[12][8]=0;
	b[12][9]=0;
	b[12][10]=0;
	b[12][11]=-1;
	b[12][12]=1;
	b[12][13]=0;
	b[12][14]=0;
	b[12][15]=0;
	b[12][16]=0;
	b[12][17]=0;
	b[12][18]=0;
	b[12][19]=0;
	b[12][20]=dh1l;

	//Enthalpy Increase in 2nd HP Heater
	b[13][0]=0;
	b[13][1]=0;
	b[13][2]=0;
	b[13][3]=0;
	b[13][4]=0;
	b[13][5]=0;
	b[13][6]=0;
	b[13][7]=0;
	b[13][8]=0;
	b[13][9]=0;
	b[13][10]=0;
	b[13][11]=0;
	b[13][12]=0;
	b[13][13]=1;
	b[13][14]=0;
	b[13][15]=0;
	b[13][16]=0;
	b[13][17]=0;
	b[13][18]=0;
	b[13][19]=0;
	b[13][20]=dh2h*Qin+681*Qin;

	//Enthalpy Increase in 1st HP Heater
	b[14][0]=0;
	b[14][1]=0;
	b[14][2]=0;
	b[14][3]=0;
	b[14][4]=0;
	b[14][5]=0;
	b[14][6]=0;
	b[14][7]=0;
	b[14][8]=0;
	b[14][9]=0;
	b[14][10]=0;
	b[14][11]=0;
	b[14][12]=0;
	b[14][13]=-1;
	b[14][14]=1;
	b[14][15]=0;
	b[14][16]=0;
	b[14][17]=0;
	b[14][18]=0;
	b[14][19]=0;
	b[14][20]=dh1h*Qin;
	
	//Mass Conservation in 1st HP Heater
	b[15][0]=0;
	b[15][1]=hf[16]/hf[2];
	b[15][2]=0;
	b[15][3]=0;
	b[15][4]=0;
	b[15][5]=0;
	b[15][6]=0;
	b[15][7]=0;
	b[15][8]=0;
	b[15][9]=0;
	b[15][10]=0;
	b[15][11]=0;
	b[15][12]=0;
	b[15][13]=0;
	b[15][14]=0;
	b[15][15]=-1;
	b[15][16]=0;
	b[15][17]=0;
	b[15][18]=0;
	b[15][19]=0;
	b[15][20]=0;
	
	//Mass Conservation in 2nd HP Heater
	b[16][0]=0;
	b[16][1]=0;
	b[16][2]=hf[17]/hf[3];
	b[16][3]=0;
	b[16][4]=0;
	b[16][5]=0;
	b[16][6]=0;
	b[16][7]=0;
	b[16][8]=0;
	b[16][9]=0;
	b[16][10]=0;
	b[16][11]=0;
	b[16][12]=0;
	b[16][13]=0;
	b[16][14]=0;
	b[16][15]=hf[17]/hf[16];
	b[16][16]=-1;
	b[16][17]=0;
	b[16][18]=0;
	b[16][19]=0;
	b[16][20]=0;
	
	//Mass Conservation in 1st LP Heater
	b[17][0]=0;
	b[17][1]=0;
	b[17][2]=0;
	b[17][3]=0;
	b[17][4]=hf[18]/hf[5];
	b[17][5]=0;
	b[17][6]=0;
	b[17][7]=0;
	b[17][8]=0;
	b[17][9]=0;
	b[17][10]=0;
	b[17][11]=0;
	b[17][12]=0;
	b[17][13]=0;
	b[17][14]=0;
	b[17][15]=0;
	b[17][16]=0;
	b[17][17]=-1;
	b[17][18]=0;
	b[17][19]=0;
	b[17][20]=0;
	
	//Mass Conservation in 2nd LP Heater
	b[18][0]=0;
	b[18][1]=0;
	b[18][2]=0;
	b[18][3]=0;
	b[18][4]=0;
	b[18][5]=hf[19]/hf[6];
	b[18][6]=0;
	b[18][7]=0;
	b[18][8]=0;
	b[18][9]=0;
	b[18][10]=0;
	b[18][11]=0;
	b[18][12]=0;
	b[18][13]=0;
	b[18][14]=0;
	b[18][15]=0;
	b[18][16]=0;
	b[18][17]=hf[19]/hf[18];
	b[18][18]=-1;
	b[18][19]=0;
	b[18][20]=0;
	
	//Mass Conservation in 3rd LP Heater
	b[19][0]=0;
	b[19][1]=0;
	b[19][2]=0;
	b[19][3]=0;
	b[19][4]=0;
	b[19][5]=0;
	b[19][6]=hf[20]/hf[7];
	b[19][7]=0;
	b[19][8]=0;
	b[19][9]=0;
	b[19][10]=0;
	b[19][11]=0;
	b[19][12]=0;
	b[19][13]=0;
	b[19][14]=0;
	b[19][15]=0;
	b[19][16]=0;
	b[19][17]=0;
	b[19][18]=hf[20]/hf[19];
	b[19][19]=-1;
	b[19][20]=0;

	    double a1[N][N+1] = {
        {b[0][0], b[0][1], b[0][2], b[0][3], b[0][4], b[0][5], b[0][6], b[0][7], b[0][8], b[0][9], b[0][10], b[0][11], b[0][12], b[0][13], b[0][14], b[0][15], b[0][16], b[0][17], b[0][18], b[0][19], b[0][20]},
        {b[1][0], b[1][1], b[1][2], b[1][3], b[1][4], b[1][5], b[1][6], b[1][7], b[1][8], b[1][9], b[1][10], b[1][11], b[1][12], b[1][13], b[1][14], b[1][15], b[1][16], b[1][17], b[1][18], b[1][19], b[1][20]},
        {b[2][0], b[2][1], b[2][2], b[2][3], b[2][4], b[2][5], b[2][6], b[2][7], b[2][8], b[2][9], b[2][10], b[2][11], b[2][12], b[2][13], b[2][14], b[2][15], b[2][16], b[2][17], b[2][18], b[2][19], b[2][20]},
        {b[3][0], b[3][1], b[3][2], b[3][3], b[3][4], b[3][5], b[3][6], b[3][7], b[3][8], b[3][9], b[3][10], b[3][11], b[3][12], b[3][13], b[3][14], b[3][15], b[3][16], b[3][17], b[3][18], b[3][19], b[3][20]},
        {b[4][0], b[4][1], b[4][2], b[4][3], b[4][4], b[4][5], b[4][6], b[4][7], b[4][8], b[4][9], b[4][10], b[4][11], b[4][12], b[4][13], b[4][14], b[4][15], b[4][16], b[4][17], b[4][18], b[4][19], b[4][20]},
        {b[5][0], b[5][1], b[5][2], b[5][3], b[5][4], b[5][5], b[5][6], b[5][7], b[5][8], b[5][9], b[5][10], b[5][11], b[5][12], b[5][13], b[5][14], b[5][15], b[5][16], b[5][17], b[5][18], b[5][19], b[5][20]},
        {b[6][0], b[6][1], b[6][2], b[6][3], b[6][4], b[6][5], b[6][6], b[6][7], b[6][8], b[6][9], b[6][10], b[6][11], b[6][12], b[6][13], b[6][14], b[6][15], b[6][16], b[6][17], b[6][18], b[6][19], b[6][20]},
        {b[7][0], b[7][1], b[7][2], b[7][3], b[7][4], b[7][5], b[7][6], b[7][7], b[7][8], b[7][9], b[7][10], b[7][11], b[7][12], b[7][13], b[7][14], b[7][15], b[7][16], b[7][17], b[7][18], b[7][19], b[7][20]},
        {b[8][0], b[8][1], b[8][2], b[8][3], b[8][4], b[8][5], b[8][6], b[8][7], b[8][8], b[8][9], b[8][10], b[8][11], b[8][12], b[8][13], b[8][14], b[8][15], b[8][16], b[8][17], b[8][18], b[8][19], b[8][20]},
        {b[9][0], b[9][1], b[9][2], b[9][3], b[9][4], b[9][5], b[9][6], b[9][7], b[9][8], b[9][9], b[9][10], b[9][11], b[9][12], b[9][13], b[9][14], b[9][15], b[9][16], b[9][17], b[9][18], b[9][19], b[9][20]},
        {b[10][0], b[10][1], b[10][2], b[10][3], b[10][4], b[10][5], b[10][6], b[10][7], b[10][8], b[10][9], b[10][10], b[10][11], b[10][12], b[10][13], b[10][14], b[10][15], b[10][16], b[10][17], b[10][18], b[10][19], b[10][20]},
        {b[11][0], b[11][1], b[11][2], b[11][3], b[11][4], b[11][5], b[11][6], b[11][7], b[11][8], b[11][9], b[11][10], b[11][11], b[11][12], b[11][13], b[11][14], b[11][15], b[11][16], b[11][17], b[11][18], b[11][19], b[11][20]},
        {b[12][0], b[12][1], b[12][2], b[12][3], b[12][4], b[12][5], b[12][6], b[12][7], b[12][8], b[12][9], b[12][10], b[12][11], b[12][12], b[12][13], b[12][14], b[12][15], b[12][16], b[12][17], b[12][18], b[12][19], b[12][20]},
        {b[13][0], b[13][1], b[13][2], b[13][3], b[13][4], b[13][5], b[13][6], b[13][7], b[13][8], b[13][9], b[13][10], b[13][11], b[13][12], b[13][13], b[13][14], b[13][15], b[13][16], b[13][17], b[13][18], b[13][19], b[13][20]},
        {b[14][0], b[14][1], b[14][2], b[14][3], b[14][4], b[14][5], b[14][6], b[14][7], b[14][8], b[14][9], b[14][10], b[14][11], b[14][12], b[14][13], b[14][14], b[14][15], b[14][16], b[14][17], b[14][18], b[14][19], b[14][20]},
        {b[15][0], b[15][1], b[15][2], b[15][3], b[15][4], b[15][5], b[15][6], b[15][7], b[15][8], b[15][9], b[15][10], b[15][11], b[15][12], b[15][13], b[15][14], b[15][15], b[15][16], b[15][17], b[15][18], b[15][19], b[15][20]},
        {b[16][0], b[16][1], b[16][2], b[16][3], b[16][4], b[16][5], b[16][6], b[16][7], b[16][8], b[16][9], b[16][10], b[16][11], b[16][12], b[16][13], b[16][14], b[16][15], b[16][16], b[16][17], b[16][18], b[16][19], b[16][20]},
        {b[17][0], b[17][1], b[17][2], b[17][3], b[17][4], b[17][5], b[17][6], b[17][7], b[17][8], b[17][9], b[17][10], b[17][11], b[17][12], b[17][13], b[17][14], b[17][15], b[17][16], b[17][17], b[17][18], b[17][19], b[17][20]},
        {b[18][0], b[18][1], b[18][2], b[18][3], b[18][4], b[18][5], b[18][6], b[18][7], b[18][8], b[18][9], b[18][10], b[18][11], b[18][12], b[18][13], b[18][14], b[18][15], b[18][16], b[18][17], b[18][18], b[18][19], b[18][20]},
        {b[19][0], b[19][1], b[19][2], b[19][3], b[19][4], b[19][5], b[19][6], b[19][7], b[19][8], b[19][9], b[19][10], b[19][11], b[19][12], b[19][13], b[19][14], b[19][15], b[19][16], b[19][17], b[19][18], b[19][19], b[19][20]},
      
        
    
	};

    gauss_jordan(a1, N);

	 PowerGen=EthaGen*a1[0][N];
	 Flow1stHpEx=a1[1][N]/hf[2];
	 Flow2ndHpEx=a1[2][N]/hf[3];
	 FlowDeaEx=a1[3][N]/hf[4];
	 Flow1stLpEx=a1[4][N]/hf[5];
	 Flow2ndLpEx=a1[5][N]/hf[6];
	 Flow3rdLpEx=a1[6][N]/hf[7];
	 FlowCond=a1[7][N]/hf[8];
	 FlowCondOut=(1.0+MakeupPercent/100)*FlowCond+mm+0.31+0.28;
	 EnthalpyCondOut=a1[8][N]/FlowCondOut; 
	// EnthalpyCondOut=hf[9];    
	// FlowCondOut=a1[8][N]/hf[9];
	 EnthalpyCP=a1[9][N]/FlowCondOut;
	 Enthalpy3rdLP=a1[10][N]/FlowCondOut;
	 FlowInt=FlowCondOut+Flow1stLpEx+Flow2ndLpEx+Flow3rdLpEx+me+md;
	 Enthalpy2ndLP=a1[11][N]/FlowInt;
	 Enthalpy1stLP=a1[12][N]/FlowInt;
	 Enthalpy2ndHP=a1[13][N]/Qin;
	 Enthalpy1stHP=a1[14][N]/Qin;
	

    
	cout << "\n----------\nPower Generated = " << PowerGen <<"  kW"<< endl;
	
	cout << "\nFlow of 1st Hp FWH Ext = " << Flow1stHpEx <<"  t/h"<< endl;
	
	cout << "Flow of 2nd HP FWH Ext = " << Flow2ndHpEx <<"  t/h"<< endl;
	
	cout << "Flow of Deaerator Ext = " << FlowDeaEx <<"  t/h"<< endl;
	
	cout << "Flow of 1st LP FWH Ext = " << Flow1stLpEx <<"  t/h"<< endl;
	
	cout << "Flow of 2nd LP FWH Ext = " << Flow2ndLpEx <<"  t/h"<< endl;
	
	cout << "Flow of 3rd LP FWH Ext = " << Flow3rdLpEx <<"  t/h"<< endl;
	
	cout << "\nCondensser Flow = " << FlowCond <<"  t/h"<< endl;
	
	cout << "Condensser Flow + Extractions + Gland Steam = " << Flow1stHpEx+Flow2ndHpEx+FlowDeaEx+Flow1stLpEx+Flow2ndLpEx+Flow3rdLpEx+FlowCond+ma+mb+mc+md+me+mf+mm <<"  t/h"<< endl;
	
	cout << "Condensser Outlet Enthalpy = " << EnthalpyCondOut <<"  kj/kg"<< endl;

	cout << "Condensser Outlet Flow = " << FlowCondOut <<"  t/h"<< endl;

	cout << "Enthalpy at Condensate Pump Outlet = " << EnthalpyCP <<"  Kj/Kg"<< endl;
	
	cout << "\nEnthalpy at 3rd LP F.W.H Outlet = " << Enthalpy3rdLP <<"  Kj/Kg"<< endl;

	cout << "Enthalpy at 2nd LP F.W.H Outlet = " << Enthalpy2ndLP <<"  Kj/Kg"<< endl;
	
	cout << "Enthalpy at 1st LP F.W.H Outlet = " << Enthalpy1stLP <<"  Kj/Kg"<< endl;
	
	cout << "\nEnthalpy at Deaerator's Outlet = " << hde <<"  Kj/Kg"<< endl;
	
	cout << "\nEnthalpy at 2nd HP F.W.H Outlet = " << Enthalpy2ndHP <<"  Kj/Kg"<< endl;
	
	cout << "Enthalpy at 1st HP F.W.H Outlet = " << Enthalpy1stHP <<"  Kj/Kg"<< endl;
	
	cout << "\nCondenser Heat Load = " << (FlowCond*hf[8]-FlowCondOut*hf[9])/3.6/1000<<"  MW"<< endl;
	
	cout << "Heat Disipated By the Cooling System = " << Qout/3.6/1000/Ethac<<"  MW"<< endl;
	
	
	////////////////////////////////////////////////...Boiler Inputs ...//////////////////////////////////////////
	
	
	 QBinf=Qin*(Hin-Enthalpy1stHP);
	 QBina=QBinf/Ethab/3600*100;
	
	 EffCycle=PowerGen/1000/QBina;
	
	
	cout <<"\n-----------------------\n";
	cout << "Boiler Heat Input = " << QBina<<"  MW"<< endl;
	cout << "Power Generated = " << PowerGen/1000 <<"  MW"<< endl;
	cout << "Cycle Efficiency = " << EffCycle*100<<" %"<< endl;
	
	 Qcond=FlowCond*(hf[8]-hf[9]);
	 deltaTct=Qcond/9650/4.179*1.2;
	cout << "Delta T Cooling Water = " << deltaTct<<" C"<< endl;
	cout << "Condenser Approach = " << Ts20-(Tcwin+deltaTct)<<" C"<< endl;
	////////////////////////////////////////////////...Printing ...//////////////////////////////////////////
	
	ofstream myfile;
    myfile.open ("Results.txt");
	myfile <<"Pressure\n"<<Ps6;
	myfile <<"\n"<<Ps8;
	myfile <<"\n"<<Ps10;
	myfile <<"\n"<<Ps14;
	myfile <<"\n"<<Ps16;
	myfile <<"\n"<<Ps18;
	myfile <<"\n"<<Ps20;
	myfile <<"\n\nTemperature\n"<<Ts6<<"\n"<<Ts8<<"\n"<<Ts10<<"\n"<<Ts14<<"\n"<<Ts16<<"\n"<<Ts18<<"\n"<<Ts20;
	myfile <<"\n\nEnthalpy\n"<<Hs6<<"\n"<<Hs8<<"\n"<<Hs10<<"\n"<<Hs14<<"\n"<<Hs16<<"\n"<<Hs18<<"\n"<<Hs20;
	//myfile <<"\n\nQuality\n"<<x6<<"\n"<<x8<<"\n"<<x10<<"\n"<<x14<<"\n"<<x16<<"\n"<<x18<<"\n"<<x20;
	myfile <<"\n\nFlow of Extractions\n"<<Flow1stHpEx<<"\n"<<Flow2ndHpEx<<"\n"<<FlowDeaEx<<"\n"<<Flow1stLpEx<<"\n"<<Flow2ndLpEx<<"\n"<<Flow3rdLpEx<<"\n"<<FlowCond;
	myfile <<"\n\nPower\n"<<PowerGen;
	myfile <<"\n\nEnthalpy of Condensate\n"<< EnthalpyCondOut <<"\n"<< EnthalpyCP <<"\n"<< Enthalpy3rdLP <<"\n"<<Enthalpy2ndLP<<"\n"<<Enthalpy1stLP<<"\n"<<hde<<"\n"<<Enthalpy2ndHP<<"\n"<<Enthalpy1stHP;
	myfile <<"\n\nCondenser Heat Load\n"<<(FlowCond*hf[8]-FlowCondOut*hf[9])/3.6/1000<<"\n"<<Qout/3600/Ethac;
	myfile <<"\n\nBoiler Heat Input\n"<<QBina<<"\n\nPower Generated\n"<<PowerGen/1000<<"\n\nCycle Efficiency\n"<<EffCycle*100;
	myfile << "\nDelta T Cooling Water\n " <<deltaTct<< endl;

  myfile.close();
	
	ofstream myfile1;
    myfile1.open ("EES Output.txt");
	myfile1 <<"T[1]="<<Tsat(Ps20);
	myfile1 <<"\nS[1]=entropy(Water, T="<<Tsat(Ps20)<<", x=0)";
	myfile1 <<"\nT[2]="<<Enthalpy3rdLP<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[2]=1.05*S[1]";
	myfile1 <<"\nT[3]="<<Enthalpy2ndLP<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[3]=1.05*S[2]";
	myfile1 <<"\nT[4]="<<Enthalpy1stLP<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[4]=1.05*S[3]";
	myfile1 <<"\nT[5]="<<hde<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[5]=entropy(Water,T=T[5],x=0)";
	myfile1 <<"\nT[6]="<<Enthalpy2ndHP<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[6]=1.05*S[5]";
	myfile1 <<"\nT[7]="<<Enthalpy1stHP<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[7]=1.05*S[6]";
	myfile1 <<"\nT[8]=t_sat(Water,P="<<Pin*100<<")";
	myfile1 <<"\nS[8]=entropy(Water,T=T[8],x=0)";
	myfile1 <<"\nT[9]=t_sat(Water,P="<<Pin*100<<")";
	myfile1 <<"\nS[9]=entropy(Water,T=T[9],x=1)";
	myfile1 <<"\nT[10]="<<Tin;
	myfile1 <<"\nS[10]="<<Ssup(Tin,Pin);
	myfile1 <<"\nT[11]="<<Ts6;
	myfile1 <<"\nS[11]="<<Ss6;
	
	myfile1 <<"\nT[12]=t_sat(Water,P="<<Ps6*100<<")";
	myfile1 <<"\nS[12]=entropy(Water,T=T[12],x=1)";
	myfile1 <<"\nT[13]=t_sat(Water,P="<<Ps6*100<<")";
	myfile1 <<"\nS[13]=entropy(Water,T=T[13],x=0)";
	myfile1 <<"\nT[14]=t_sat(Water,P="<<Ps8*100<<")";
	myfile1 <<"\nS[14]=entropy(Water,T=T[14],x=0)";
	myfile1 <<"\nT[15]=t_sat(Water,P="<<Ps6*100<<")";
	myfile1 <<"\nS[15]=entropy(Water,T=T[13],x=0)";
	myfile1 <<"\nT[16]=t_sat(Water,P="<<Ps6*100<<")";
	myfile1 <<"\nS[16]=entropy(Water,T=T[12],x=1)";
	myfile1 <<"\nT[17]="<<Ts6;
	myfile1 <<"\nS[17]="<<Ss6;
	myfile1 <<"\nT[18]="<<Ts8;
	myfile1 <<"\nS[18]="<<Ss8;
	
	myfile1 <<"\nT[19]=t_sat(Water,P="<<Ps8*100<<")";
	myfile1 <<"\nS[19]=entropy(Water,T=T[19],x=1)";
	myfile1 <<"\nT[20]=t_sat(Water,P="<<Ps8*100<<")";
	myfile1 <<"\nS[20]=entropy(Water,T=T[20],x=0)";
	myfile1 <<"\nT[21]=t_sat(Water,P="<<Ps10*100<<")";
	myfile1 <<"\nS[21]=entropy(Water,T=T[21],x=0)";
	myfile1 <<"\nT[22]=t_sat(Water,P="<<Ps8*100<<")";
	myfile1 <<"\nS[22]=entropy(Water,T=T[22],x=0)";
	myfile1 <<"\nT[23]=t_sat(Water,P="<<Ps8*100<<")";
	myfile1 <<"\nS[23]=entropy(Water,T=T[23],x=1)";
	myfile1 <<"\nT[24]="<<Ts8;
	myfile1 <<"\nS[24]="<<Ss8;
	myfile1 <<"\nT[25]="<<Ts10;
	myfile1 <<"\nS[25]="<<Ss10;
	
	myfile1 <<"\nT[26]=t_sat(Water,P="<<Ps10*100<<")";
	myfile1 <<"\nS[26]=entropy(Water,T=T[26],x=1)";
	myfile1 <<"\nT[27]=t_sat(Water,P="<<Ps10*100<<")";
	myfile1 <<"\nS[27]=entropy(Water,T=T[27],x=0)";
	myfile1 <<"\nT[28]="<<hde<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[28]=S[5]";
	myfile1 <<"\nT[29]=t_sat(Water,P="<<Ps10*100<<")";
	myfile1 <<"\nS[29]=entropy(Water,T=T[29],x=0)";
	myfile1 <<"\nT[30]=t_sat(Water,P="<<Ps10*100<<")";
	myfile1 <<"\nS[30]=entropy(Water,T=T[30],x=1)";
	myfile1 <<"\nT[31]="<<Ts10;
	myfile1 <<"\nS[31]="<<Ss10;
	myfile1 <<"\nT[32]="<<Ts14;
	myfile1 <<"\nS[32]="<<Ss14;
	
	myfile1 <<"\nT[33]=t_sat(Water,P="<<Ps14*100<<")";
	myfile1 <<"\nS[33]=entropy(Water,T=T[33],x=1)";
	myfile1 <<"\nT[34]=t_sat(Water,P="<<Ps14*100<<")";
	myfile1 <<"\nS[34]=entropy(Water,T=T[34],x=0)";
	myfile1 <<"\nT[35]=t_sat(Water,P="<<Ps16*100<<")";
	myfile1 <<"\nS[35]=entropy(Water,T=T[42],x=0)";
	myfile1 <<"\nT[36]=t_sat(Water,P="<<Ps14*100<<")";
	myfile1 <<"\nS[36]=entropy(Water,T=T[36],x=0)";
	myfile1 <<"\nT[37]=t_sat(Water,P="<<Ps14*100<<")";
	myfile1 <<"\nS[37]=entropy(Water,T=T[37],x=1)";
	myfile1 <<"\nT[38]="<<Ts14;
	myfile1 <<"\nS[38]="<<Ss14;
	myfile1 <<"\nT[39]="<<Ts16;
	myfile1 <<"\nS[39]="<<Ss16;
	
	myfile1 <<"\nT[40]=t_sat(Water,P="<<Ps16*100<<")";
	myfile1 <<"\nS[40]=entropy(Water,T=T[40],x=0)";
	myfile1 <<"\nT[41]=t_sat(Water,P="<<Ps18*100<<")";
	myfile1 <<"\nS[41]=entropy(Water,T=T[47],x=0)";
	myfile1 <<"\nT[42]=t_sat(Water,P="<<Ps16*100<<")";
	myfile1 <<"\nS[42]=entropy(Water,T=T[42],x=0)";
	myfile1 <<"\nT[43]="<<Ts16;
	myfile1 <<"\nS[43]="<<Ss16;
	myfile1 <<"\nT[44]="<<Ts18;
	myfile1 <<"\nS[44]="<<Ss18;
	
	myfile1 <<"\nT[45]=t_sat(Water,P="<<Ps18*100<<")";
	myfile1 <<"\nS[45]=entropy(Water,T=T[45],x=0)";
	myfile1 <<"\nT[46]="<<Enthalpy2ndLP<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[46]=S[3]";
	myfile1 <<"\nT[47]=t_sat(Water,P="<<Ps18*100<<")";
	myfile1 <<"\nS[47]=entropy(Water,T=T[47],x=0)";
	myfile1 <<"\nT[48]="<<Ts18;
	myfile1 <<"\nS[48]="<<Ss18;
	myfile1 <<"\nT[49]="<<Ts20;
	myfile1 <<"\nS[49]="<<Ss20;
	
	myfile1 <<"\nT[50]=t_sat(Water,P="<<Ps20*100<<")";
	myfile1 <<"\nS[50]=entropy(Water,T=T[50],x=0)";
	myfile1 <<"\nT[51]="<<EnthalpyCondOut<<"/specheat(Water, T=70, P=100)";
	myfile1 <<"\nS[51]=S[1]";


  myfile1.close();
		

	return 0;
}

