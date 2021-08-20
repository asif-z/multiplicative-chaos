# include <cstdlib>
#include <iostream>
#include <math.h> 
#include <cmath>
#include <stdio.h>
#include <iomanip> 
#include <algorithm>
#include <cstring>
#include <chrono>
#include <sstream>
#include <string>
#include <fstream>
#include <random>
#include <omp.h>
#include <mpi.h>




using namespace std;
using namespace std::chrono;

const int N = 20000;
/*If a partition of n has n or more number of part k, then (X(k))^(n)/(n)! would be smaller than the machine level precision for a sufficiently large n.
So, we can assign 0 to such contributions. In our case with double data type, we can choose such n to be 30. */
const int approx_const = 30;   
const int samples = 10000000;

/*Considering the limitation of memory, each node of the computation processes a sample size of 200 at one time. */
const int chunk = 200;   
 

const string output_sample_mean_fileName = "sample_mean_A_N_upto_20k_with_10mil_pm1_samples.txt"; 

/*The following function is the implementation of our algorithm to generate A_N.*/
/* X_k_N is an array of individual sample sequence. That is, X(k) for 1<=k<=N. */
void generate_A_N(int N, double  *fact, double *X_k_N, double * A_k_n){

 
   A_k_n[0] =  1.0;
   

   double A_k_n_temp[N+1];
   
  
   A_k_n_temp[0] = 1.0;
   
   /* In k^th iteration, A_k_n_temp contains the contribution to A(n) by the partitions with largest part k. While transitioning 
   from k^th to (k+1)^th iteration, A_k_n is updated so that it constains contribution to A(n) by the partitions with largest part <=k.
   After the final iteration (i.e with k=N), A_k_n actually becomes the sequence of A(n) for 1 <=n <= N. 
*/


   double X_1 = X_k_N[1];
   


   for (int i = 1; i <= approx_const; i++){  
       A_k_n[i] = pow(X_1, (double) i)*fact[i]; 
       
   }

    

   for(int k = 2; k <= N/2; k++){
       double X_k = X_k_N[k];
       double X_k_N_m[N/k];

       for(int m = 0; m <= std:: min(N/k, approx_const); m++){
           X_k_N_m[m] = pow(X_k, (double) m) * fact[m];
       }

        
    

       for(int n = k; n < 2*k; n++){
           
           A_k_n_temp[n] = X_k * A_k_n[n-k];
         }

        for (int n = 2*k; n <= std:: min(N, approx_const*k); n++){
             double temp_sum = 0;
               for(int m = 1; m <= n/k ; m++){
                   temp_sum += A_k_n[n-m*k]*X_k_N_m[m];
               }
               A_k_n_temp[n] = temp_sum;

        }
              
           
    

       for(int i = k; i <= std:: min(N, approx_const*k); i++){
           A_k_n[i] += A_k_n_temp[i];
       }

        

   } 

   for (int k = N/2+1; k <= N; k++){
       double X_k = X_k_N[k];
       
       for(int n = k; n <= N; n++){
           A_k_n[n] += X_k * A_k_n[n-k];
       }
   }


}

/* The following function finally computes the sample_mean and write it in an output text file.*/
void write_sample_mean(double *sample_mean_sum){
    std::ofstream myFile;
    myFile.open(output_sample_mean_fileName, std::ios::app);
    myFile.precision(15);
 
    for(int i = 0; i <= N; i++){
        myFile << sample_mean_sum[i]/samples;
        if(i!= N){
        myFile << ",";
        }

    }

    myFile << "\n";
    
}


/* This function sums up A(n) and store in an array, which is finally divided by the sample size to compute the sample_mean.*/

void compute_sample_mean_sum(double A_N[][N+1], double * sample_mean_sum_each_node){

    for(int i = 0; i < chunk; i++){
        for(int j = 0; j <= N; j++){
            double temp = A_N[i][j];
            sample_mean_sum_each_node[j] += std::abs(temp);
        }
        
    }

}

/*The following function distributes the computation among all the cores within each node to generate A(n) parallely.*/

 void generate_data_parallely(double A_N_distribution[][N+1], double *fact){
        
            #pragma omp parallel for shared(A_N_distribution, fact)
            for(int i = 0; i < chunk; i++){
                double A_N[N+1]  = {0.0};
                generate_A_N(N, fact, A_N_distribution[i], A_N);
                memcpy(A_N_distribution[i], A_N, sizeof(A_N));
            }
} 



void generate_distribution(double * root, double * fact, double * sample_mean_sum_each_node, int nodes, int status){
    

    /* We use mt19937 random generator to generate normally distributed X(j). This generator requires a seed, which is also 
    generated randomly using random_device.*/
    std::random_device generator_seed;
    std::mt19937 generator(generator_seed());
    std::normal_distribution<double> distribution(0.0, 1.0);

    double sample_sets[chunk][N+1];

    

    for(int i =0; i < N+1; i++){
        sample_mean_sum_each_node[i] = 0.0;
    }

    /*The total sample size is divided into multiple nodes and each node operate on one chunk at a time because of memory limitations.
    The variable c counts the total number of iterations in each node.  */
   
    for(int c = 0; c < (samples/nodes/chunk); c++){

        /*The following loop populate each chunk of sample array with randomly generated standard normal samples. Since the contribution
     of part j is always a function of X(j)/sqrt(j), we store value normalized by squared-root.*/
        for(int i = 0; i < chunk ; i++){
            sample_sets[i][0] = 1;
            
             
            for(int j = 1; j<= N; j++){

             
                int rv = -1;

                if(distribution(generator)){
                    rv = 1;
                }

                sample_array[i][j] = rv * root[j];
        
            }
            
        }

         
        
        generate_data_parallely(sample_sets, fact);
        

        compute_sample_mean_sum(sample_sets, sample_mean_sum_each_node);

        /*The following if-clause is just to track the completed computations for a big job.*/
        if(status == 0){
            
            cout << "number of completed iterations: " << (c+1)*chunk*nodes << endl;
        }
        
    }

}





int main(int argc, char *argv[]){

    time_t t_1;
    time_t t_2;

    
    /*Rank is a numeric label provided to all operating nodes in the computation. Rank=0 is the root node. Size is the total
    number of nodes.*/

    int rank, size;

    double fact[approx_const+1];

    fact[0] = 1.0;

    /*The following loop populates reciprocals of factorials using in build gamma function in C++ library.*/
    for(int i = 1; i <= approx_const; i++){

        fact[i] = 1./tgammal((long double)i+1);
        

    }

    double root[N+1];

    for(int i = 1; i <= N; i++){
        root[i] = 1./sqrt((long double) i);
        
    } 


    

    double sample_mean_sum[N+1] = {0.0};



    double *global_sample_mean_sum_array = NULL;

    double local_sample_mean_sum_array[N+1];
     
    


    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    
    
    if(rank == 0){
        t_1 = std::time(0);
        cout << "number of MPI processes: " << size << endl;
        cout << "number of maximum possible threads: "<< omp_get_max_threads()<< endl;
        global_sample_mean_sum_array = new double[(N+1)*size];
        
        
    }

    /*MPI_Scatter distributes the computation accross all nodes. */


    MPI_Scatter(global_sample_mean_sum_array, N+1, MPI_DOUBLE, local_sample_mean_sum_array, N+1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    /*The following function operates in all communicating nodes. As seen above, this function calls all othe functions in this program.*/
    generate_distribution(root, fact, local_sample_mean_sum_array, size, rank);

    MPI_Gather(local_sample_mean_sum_array, N+1, MPI_DOUBLE, global_sample_mean_sum_array, N+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*MPI_Gather gathers data from all communicating nodes.*/
    MPI_Finalize(); 

    if(rank == 0){

        for(int i =0 ; i < size; i++){
            for(int j =0; j< N+1; j++){
                sample_mean_sum[j] += global_sample_mean_sum_array[i*(N+1)+j];
            }

        } 

        delete [] global_sample_mean_sum_array;
        write_sample_mean(sample_mean_sum);

        t_2 = std::time(0);
        cout << "total time of completion: " << t_2-t_1 << endl;
     

    }

   
   return 0;

}
