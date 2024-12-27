#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sim_bp.h"
#include <math.h>
#include <vector>

/*  argc holds the number of command line arguments
    argv[] holds the commands themselves

    Example:-
    sim bimodal 6 gcc_trace.txt
    argc = 4
    argv[0] = "sim"
    argv[1] = "bimodal"
    argv[2] = "6"
    ... and so on
*/
//int BimodalPredictionTable[1000000];
std::vector<int> BimodalPredictionTable(10000000,2);
std::vector<int> GsharePredictionTable(10000000,2);
std::vector<int> ChooserTable(10000000, 1);


int mispredictions=0;
int hyb_binomal_mispred=0, hyb_gshare_mispred=0;
int GHR=0;
void branch_predict(char prediction, unsigned long index, char outcome, bool ishybrid)
{ if(prediction=='b')
  {
    if (outcome=='t')
        { //printf("H");
             if (BimodalPredictionTable[index]<2)
            {
                
                mispredictions++;
            }

            if(BimodalPredictionTable[index]<3)
                BimodalPredictionTable[index]++;

        }

         if (outcome=='n')
        { 
             if (BimodalPredictionTable[index]>=2)
            {
                
                mispredictions++;
            }
                if(BimodalPredictionTable[index]>0)
                BimodalPredictionTable[index]--;
        }
  }   

  if(prediction=='g')
  {
    if (outcome=='t')
        { 
             if (GsharePredictionTable[index]<2)
            {
                
                mispredictions++;
            }

            if(GsharePredictionTable[index]<3)
                GsharePredictionTable[index]++;

        }

         if (outcome=='n')
        { 
             if (GsharePredictionTable[index]>=2)
            {
                
                mispredictions++;
            }
                if(GsharePredictionTable[index]>0)
                 GsharePredictionTable[index]--;
        }
  }   

}


void bimodal_branch_predict(unsigned long addr, char outcome, int m, bool ishybrid)
{   
    addr=addr>>2;
   
    unsigned long mask = (1 << m) - 1;
    
    unsigned long index = addr & mask;
    

    branch_predict('b',index, outcome, ishybrid);
     
}


void updateGHR(char outcome, int n)
{
    GHR=GHR>>1;
    unsigned long temp_MSB = 1<<(n-1);
    if (outcome == 't')
    {  
        GHR = temp_MSB|GHR;
        
    }
}


void g_share_branch_predictor(unsigned long addr, char outcome, int m, int n, bool ishybrid)
{
    addr = addr>> 2; 
    unsigned long M_mask = (1 << m) - 1;
    unsigned long N_mask = (1 << n) - 1;
    unsigned long M_bits = addr & M_mask;
    unsigned long PC_bits = (M_bits >> (m - n)) & N_mask;
    unsigned long index = ((PC_bits ^ (GHR & N_mask)) << (m - n)) | (M_bits & ((1 << (m - n)) - 1));
    
    branch_predict('g', index, outcome, ishybrid);
    updateGHR(outcome, n);

    
}

void hybrid_branch_predictor(unsigned long addr, char outcome, int k, int m2, int m1, int n) 
{ //printf("h");
    
    
    unsigned long Kaddr = addr >> 2; 
    unsigned long chooser_index = (Kaddr) & ((1 << k) - 1);

    unsigned long mask = (1 << m2) - 1;
    
    unsigned long Binomalindex = Kaddr & mask;

    unsigned long M_mask = (1 << m1) - 1;
    unsigned long N_mask = (1 << n) - 1;
    unsigned long M_bits = Kaddr & M_mask;
    unsigned long PC_bits = (M_bits >> (m1 - n)) & N_mask;
    unsigned long Gshareindex = ((PC_bits ^ (GHR & N_mask)) << (m1 - n)) | (M_bits & ((1 << (m1 - n)) - 1));
    
    unsigned long temp_bi = BimodalPredictionTable[Binomalindex];
    unsigned long temp_gs = GsharePredictionTable[Gshareindex];

    if (ChooserTable[chooser_index] >= 2)
    {
        
        branch_predict('g', Gshareindex, outcome, false);
        updateGHR(outcome, n);

    } 
    else
    {
       branch_predict('b', Binomalindex, outcome, false);
       updateGHR(outcome, n);
    }

    if(((temp_bi>=2 && outcome == 'n') || (temp_bi<2 && outcome == 't')) && ((temp_gs>=2 && outcome == 't') || (temp_gs<2 && outcome == 'n')))

{ //printf("Test\n");
    if(ChooserTable[chooser_index]<3)
    ChooserTable[chooser_index]++;

}
 else if(((temp_bi>=2 && outcome == 't') || (temp_bi<2 && outcome == 'n')) && ((temp_gs>=2 && outcome == 'n') || (temp_gs<2 && outcome == 't')))

{
    if(ChooserTable[chooser_index]>0)
    ChooserTable[chooser_index]--;

}  
    
}

    


int main (int argc, char* argv[])
{  
    
    int total_table_vals;
    FILE *FP;               // File handler
    char *trace_file;       // Variable that holds trace file name;
    bp_params params;       // look at sim_bp.h header file for the the definition of struct bp_params
    char outcome;           // Variable holds branch outcome
    unsigned long int addr; // Variable holds the address read from input file
    bool bimodal= false,gshare = false,hybrid = false;
    
    if (!(argc == 4 || argc == 5 || argc == 7))
    {
        printf("Error: Wrong number of inputs:%d\n", argc-1);
        exit(EXIT_FAILURE);
    }
    
    params.bp_name  = argv[1];
    
    // strtoul() converts char* to unsigned long. It is included in <stdlib.h>
    if(strcmp(params.bp_name, "bimodal") == 0)              // Bimodal
    {
        if(argc != 4)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M2       = strtoul(argv[2], NULL, 10);
        trace_file      = argv[3];
        printf("COMMAND\n%s %s %lu %s\n", argv[0], params.bp_name, params.M2, trace_file);
        BimodalPredictionTable.resize(pow(2,params.M2) , 2);
        bimodal = true;
    }
    else if(strcmp(params.bp_name, "gshare") == 0)          // Gshare
    {
        if(argc != 5)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.M1       = strtoul(argv[2], NULL, 10);
        params.N        = strtoul(argv[3], NULL, 10);
        trace_file      = argv[4];
        printf("COMMAND\n%s %s %lu %lu %s\n", argv[0], params.bp_name, params.M1, params.N, trace_file);
        GsharePredictionTable.resize(pow(2,params.M1) , 2);
        gshare = true;

    }
    else if(strcmp(params.bp_name, "hybrid") == 0)          // Hybrid
    {
        if(argc != 7)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc-1);
            exit(EXIT_FAILURE);
        }
        params.K        = strtoul(argv[2], NULL, 10);
        params.M1       = strtoul(argv[3], NULL, 10);
        params.N        = strtoul(argv[4], NULL, 10);
        params.M2       = strtoul(argv[5], NULL, 10);
        trace_file      = argv[6];
        printf("COMMAND\n%s %s %lu %lu %lu %lu %s\n", argv[0], params.bp_name, params.K, params.M1, params.N, params.M2, trace_file);
        BimodalPredictionTable.resize(pow(2,params.M2) , 2);
        GsharePredictionTable.resize(pow(2,params.M1) , 2);
        hybrid = true;
        gshare = true;
        bimodal = true;
    }
    else
    {
        printf("Error: Wrong branch predictor name:%s\n", params.bp_name);
        exit(EXIT_FAILURE);
    }
    
    // Open trace_file in read mode
    FP = fopen(trace_file, "r");
    if(FP == NULL)
    {
        // Throw error and exit if fopen() failed
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }
    
    char str[2];

int totalcount=0;

    while(fscanf(FP, "%lx %s", &addr, str) != EOF)
    {
        totalcount++;
        outcome = str[0];
       // if (outcome == 't')
       //     printf("%lx %s\n", addr, "t");           // Print and test if file is read correctly
      //  else if (outcome == 'n')
      //      printf("%lx %s\n", addr, "n");          // Print and test if file is read correctly
        /*************************************
            Add branch predictor code here
        **************************************/
         if(strcmp(params.bp_name, "bimodal") == 0)     
       { 
         bimodal_branch_predict(addr, outcome, params.M2, false);
       }
       else if(strcmp(params.bp_name, "gshare") == 0)     
       {
         g_share_branch_predictor(addr, outcome, params.M1, params.N, false);
       } 
       else if(strcmp(params.bp_name, "hybrid") == 0) 
       {
        hybrid_branch_predictor(addr,outcome, params.K ,params.M2, params.M1, params.N);
       }
   
    }
    printf("OUTPUT\n");
    printf("number of predictions:    %d\nnumber of mispredictions: %d\n", totalcount, mispredictions);
    printf("misprediction rate:       %0.2f%%", (float(mispredictions*100)/float(totalcount)));

if(hybrid)
{
printf("\nFINAL CHOOSER CONTENTS");

 int klimit=pow(2,params.K);
    for(int i=0; i<klimit; i++)
     {
     printf("\n%d      %d", i,ChooserTable[i]);
     }

}
if(gshare)
{
printf("\nFINAL GSHARE CONTENTS");

 int m1limit=pow(2,params.M1);
    for(int i=0; i<m1limit; i++)
     {
     printf("\n%d      %d", i,GsharePredictionTable[i]);
     }

}
if(bimodal)
{
printf("\nFINAL BIMODAL CONTENTS");

 int m2limit=pow(2,params.M2);
    for(int i=0; i<m2limit; i++)
     {
     printf("\n%d      %d", i,BimodalPredictionTable[i]);
     }

}
     return 0; 
 }

