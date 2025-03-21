#ifndef SLIDINGMINPDPLUGIN_H
#define SLIDINGMINPDPLUGIN_H


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include "PluginProxy.h"

#define TOLERANCE 1.0e-10

/* Period parameters */  
//#define N 624
//#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

class SlidingMinPDPlugin : public Plugin {
	      public:
           void input(std::string);;
           void run();;
           void output(std::string);;

   private:
int z_rndu=137;
unsigned w_rndu=13757;
unsigned long mt[624]; /* the array for the state vector  */
int mti=625; /* mti==N+1 means mt[N] is not initialized */
        FILE    *fp1;
        char    *cc, inputFile[300], outputFile1[300], outputFile2[300];
        char    temp[LINELIMIT], temp2[MAXLENGHT], timechar[3];
        int             i, k, n,j,win_Count, fr_Dim, align, opt, res_no =0;
        double  **dist, ***dist_frags;
        double  GOP, GEP, match, mismatch,  mutrate;
        unsigned int    maxDim;
        Mintaxa *minresults[NUMSEQS];
        int             wSize=0, stSize,numWins;
        short int bootreps=100;
int closestRel=0;
double PCCThreshold=0;
int model=TN93;
int seed=-3;
double alpha=0.5;
int fullBootscan=false;
int distPenalty=0;
int BootThreshold=96;
int BootTieBreaker=1;
int crossOpt=1;
int Bootstrap=0;
int recOn = 1;
int reportDistances=0;
double gapPenalty = 1;/* 1 is default: means gap columns are ignored are not counted*/
int printBoot=0;
char codonFile[200];
int clustBoot=0;
int *codonList;

Fasta	*seqs[NUMSEQS]; /* Declare an array of Fasta structures */ 
BKPNode*BootscanCandidates(FILE *fp4,  Fasta **seqs, int  *Frag_Seq, int step, int windowSize , short int bootReps,
						 int numSeqs, int s, int numWins, int min_seq, double min_seq_dist, int align, int maxDim, int op, int *avgBS);;
//void ErrorMessageExit(char *errMsg,int opt);
void freeBKPNode(BKPNode *bkpLL);
int BitCriteria(int Num, int Bit);/* Is bit set in Num */
int TwoCrossovers(int Num, int MaxBit, int *Begin, int *End);/* Is bit set in Num */
double SumofSquaredValues(double ***values, int **values2, int points, int data, int s,double *m_x);
char **  MatrixInit(int s, int dim1, int dim2);
double JC69distance(int gaps, int matches, int al_len, double alpha);
double K2P80distance(int p1, int p2, int gaps, int matches, int al_len, double alpha);
double TN93distance(int pT, int pC, int pG, int pA, int p1, int p2, int gaps,int matches, int al_len, double alpha, int s, int t);
double PairwiseDistance(char **seqsChars,int t, int s, int windowSize, int model);
int DistOnly(Fasta **seqs, double **dist, double ***dist_frags, int maxDim, int n,int win_Count, int model, int windowSize, int step, int numWins);
int SaveMinResults(Mintaxa **minresults, int *resno, int a, int c, double dist, double div);
void OptionalRecPrintout(FILE *fp4,  RecRes *recRes,Fasta **seqs,int *Frag_Seq,FDisNode (*recSolutions)[3][2], int i, int l, int r, int p, int q, int f, int index, int last, int BKP, int win_Count,int align,int fr_size);
double FillArrayInOrder(int *arrayInOrder, double ***dist_frags, double **bootVals, int  *Frag_Seq, int s, int k, int noSeqs, int bootOpt);
BKPNode * ArrayToBKPNode( int  *F, int win_Count, int fr_size,int align, int bootOpt, int step);
BKPNode * GetRecSolutionsII( double ***dist_frags, double **bootVals, Fasta **seqs, int  *Frag_Seq, int noSeqs, int s, int win_Count, int min_seq, int fr_size, double min_seq_dist, int align, int bootOpt, int step, int *avgB); 
double FindBestinStretch( double ***dist_frags, double **bootVals, int  *Frag_Seq,   int bootOpt,int noSeqs, int s, int win_Count, int candidate,int lb,int ub,double largest);
BKPNode * GetRecSolutionsI(double ***dist_frags, double **bootVals, Fasta **seqs, int  *Frag_Seq,  int noSeqs, int s, int win_Count, int min_seq, int fr_size, double min_seq_dist, int align, int bootOpt, int step,int *avgBS);
int PickSeqofLargerAvgDis(double ***dist_frags, int k_idx, int f_idx, int min_idx,  int s, int win_Count, int k, int f,double min_d);
int PickSeqMinDist(double ***dist_frags, int  *Frag_Seq, int min_idx,int noSeqs, int k_idx, int f_idx, int s, int win_Count, int k, int f,double min_d);
void PrintFragments(FILE *fp3, int  *Frag_Seq, Fasta **seqs, double **dist, double ***dist_frags, int win_Count, int i, int s, int align );
int PickCandidates( Fasta **seqs, int  *Frag_Seq, double ***dist_frags, int win_Count,int min_idx, int s, int startAS, double min_d);
BKPNode * Check4Recombination(FILE *fp3, FILE *fp4, int  *Frag_Seq,  Fasta **seqs, double **dist, double ***dist_frags, int startAS, int s, int win_Count,  int min_idx,  int align, int maxDim, int opt, int step, int windowSize , short int bootReps, double min_d);
double GetBootstrapValues(Fasta **seqs, double ***bootVals,  int s, int a, int bootreps);
double GetClusteredBootstrapValues(Mintaxa *mintaxa, int max_times, Fasta **seqs, double ***bootVals,  int s, int min_idx, int bootreps);
int CheckClustering(Mintaxa *mintaxa,int max_times,Fasta **seqs,int s,int min_idx);
void StoreBootstrapValues(int *boots, Mintaxa **minresults, int resno, int start);
int OutputResults(int win_Count, int s, int startAS, double iBase, double min_d, int min_idx, int *times_max, double *div, BKPNode *bkpLL, double ***bootVals, Mintaxa *mintaxa,Mintaxa **minresults, Fasta **seqs, double **dist, int maxDim, int n, FILE *fp1, double mutrate, int fr_size, int *resno, int align, int avgBS,int bootreps);
int GetMinDist(Mintaxa **minresults, Fasta **seqs, double **dist, double ***dist_frags, int maxDim, int n, char *File1, char *File2, double mutrate, int win_Count, int *resno, int align, int opt, int step, int winSize);
void BuildPartialNJTrees(Fasta **seqs, Mintaxa **minresults, int resno, double **dist, char *File1, int n);
int PrepareWeights(int *wmod,int bootreps, int windowSize);
int BootRip(short int bootreps,  int model, double alpha, int noSeqs, int windowSize, char *seqs,  double *dist,  int *wmod);
int GetBootNJdistances(int bootreps,int noSeqs, double *dstMat,Fasta **seqs, int distPenalty);
int DoBootscan(double ***dist_frags,int **plotVal,int numWins, Fasta **seqsBS,  int maxDim, int noSeqs,int stepSize, int windowSize,int closestRel, int opt, short int bootreps, int fullBootscan);
 int BootscanShrinkPool( FILE *fp4, Fasta **seqs, int  *Frag_Seq, int step, int windowSize , short int bootReps,int numSeqs, int s, int numWins, int maxDim, int bootThreshold, int spikeLen);
int CopyToSeqsBootscan(Fasta **seqsBS, Fasta **seqs,int numSeqs, int seqNo, int *maxBSDim);
int GetParentsFromAnswer(Fasta **seqs,char *temp,int *parents, int n);
int Parents_Missed(int *Frag_Seq,int numSeqs,int *parents,int p);
int GetAvgBootstrapDist(Fasta **seqs,double **dist, double ***bootVals, int n, int bootreps);
int Bootscan(char *File1, char *File2, double ***dist_frags, int numWins, Mintaxa **minresults,Fasta **seqs, double **dist, int maxDim, int n,int step, int windowSize, int closestRel, int *resno, short int bootreps, int align);
void ReadUntil(FILE *fv, char stopChar, char *what);
void PrintUsage();
void ReadParams( char *inputFile, char *outputFile, int *opt, int *align,  int *wSize, int *stSize, char* ifile);
void CheckTimes(Fasta **seqs, int n);

void heapify_mintaxa ( Mintaxa **list , int newnode );
void heapsort_mintaxa ( Mintaxa **list, int last );

void heapify_Fasta ( Fasta **list , int newnode );
void heapsort_Fasta ( Fasta **list, int last );
void PrintTree(FILE *fv, TNode *node, Fasta **seqs,int *boots, int printboot);
void AddTipOrNode(int *w,int *p,int *max_w, int chidx, double brlen,TNode **NodeStorage, int Join);
void NJTree(TNode **NodeStorage, double **DistMatrix, int UBound, int OutgroupIn0, int *seqIds, int *max_w, int *w);
int MatricizeTopology(TNode *node, int *topoMatrix, int n, int **dist, double avgBrLen,int *max);
double GetAvgBrLen(TNode *node, int *nodes);
void SetSeed (unsigned long seed);
double rndu (void);

double rndgamma (double s);
double zeroin(double ax, double bx, double (*f)(double x));
int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);

};

#endif
