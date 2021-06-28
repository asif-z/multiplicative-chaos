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
#include <complex>
#include <cmath>



using namespace std;
using namespace std::chrono;

/* const int N = 50000;
const int approx_const = 30;
const int samples = 50;
const int chunk = 5;  */
 
const int N = 5000;
/*If a partition has a part k greater than or equal to n times, then (X(j))^(n)/(n)! would be smaller than the machine level precision.
So, we can assign 0 to such contributions. In our case with double data type, such n is 30. */
const int approx_const = 30;
const int samples = 100;
/*Because of memory limit, each node of the computation only process a sample size of 200 at one time. */
const int chunk = 20; 

 


const string output_expectation_fileName = "expectation_A_n_100mil_complex_normal_samples_upto_5k.txt";


void generate_A_n(int N, double  *fact, complex <double> *X_k_n, complex <double> * A_k_n){

 
   A_k_n[0] =  1.0 ;
  
   

   complex <double> A_k_n_temp[N+1];
  
   A_k_n_temp[0] = 1.0;
   
   
    

   complex <double> x_1 = X_k_n[1];
   


   for (int i = 1; i <= approx_const; i++){  
       A_k_n[i] = pow(x_1,  i)*fact[i]; 
       
   }

    

   for(int k = 2; k <= N/2; k++){
       complex <double> X_k = X_k_n[k];
       complex <double> X_k_n_m[N/k];

       for(int m = 0; m <= std:: min(N/k, approx_const); m++){
           X_k_n_m[m] = pow(X_k, m) * fact[m];
       }

        
    

       for(int n = k; n < 2*k; n++){
           
           A_k_n_temp[n] = X_k * A_k_n[n-k];
        }
            
        for (int n = 2*k; n <= std:: min(N, approx_const*k); n++){
             complex <double> temp_sum = 0.0;
             
               for(int m = 1; m <= n/k ; m++){
                   temp_sum += A_k_n[n-m*k]*X_k_n_m[m];
               }
               A_k_n_temp[n] = temp_sum;

        }
              
           

        for(int i = k; i <= std:: min(N, approx_const*k); i++){
           A_k_n[i] += A_k_n_temp[i];
       } 
 

   } 
    
   for (int k = N/2+1; k <= N; k++){
      complex <double> X_k = X_k_n[k];
       
       for(int n = k; n <= N; n++){
           A_k_n[n] += X_k * A_k_n[n-k];
       }
   }

    
}

/* The following function finally computes the expectation and write it in an output text file.*/
void write_expectation(double *expectation_sum){
    std::ofstream myFile;
    myFile.open(output_expectation_fileName, std::ios::app);
    myFile.precision(15);
 
    for(int i = 0; i <= N; i++){
        myFile << expectation_sum[i]/samples;
        if(i!= N){
        myFile << ",";
        }

    }

    myFile << "\n";
    
}


/* This function sums up A(n) and store in an array, which is finally divided by the sample size to compute the expectation.*/
void compute_moments(complex <double> A_n_distribution[][N+1], double * expectation_sum_each_node){

    for(int i = 0; i < chunk; i++){
        for(int j = 0; j <= N; j++){
            expectation_sum_each_node[j] += std::abs(A_n_distribution[i][j]);
        }
    
        
    }

}


/*The following function distributes the computation among all the cores within each node to generate A(n) parallely.*/

 void generate_data_parallely(complex <double> A_n_distribution[][N+1], double *fact){
        
            #pragma omp parallel for shared(A_n_distribution, fact)
            for(int i = 0; i < chunk; i++){
                complex <double> A_n[N+1]  = {0.0};
                generate_A_n(N, fact, A_n_distribution[i], A_n);
                memcpy(A_n_distribution[i], A_n, sizeof(A_n));
            }
} 

void generate_distribution(double * root, double * fact, double * expectation_sum_each_node, int nodes, int status){
    

    /* We use mt19937 random generator to generate normally distributed X(j). This generator requires a seed, which is also 
    generated randomly using random_device.*/

    std::random_device generator_seed;
    std::mt19937 generator(generator_seed());
    std::normal_distribution<double> distribution(0.0, root[2]);

   complex <double> sample_array[chunk][N+1];

    

    for(int i =0; i < N+1; i++){
        expectation_sum_each_node[i] = 0.0;
    }

    double real = 0.0;
    double imaginary = 0.0;
     /*The total sample size is divided into multiple nodes and each node operate on one chunk at a time because of memory limitations.
    The variable c counts the total number of iterations in each node.  */
   
    for(int c = 0; c < (samples/nodes/chunk); c++){

          /*The following loop populate each chunk of sample array with randomly generated standard normal samples. Since the contribution
     of part j is always a function of X(j)/sqrt(j), we store value normalized by squared-root.*/

        for(int i = 0; i < chunk ; i++){
            sample_array[i][0] = 1.0;
            
            for(int j = 1; j<= N; j++){

                real = distribution(generator);
                imaginary = distribution(generator);  

    

                
                
                sample_array[i][j] = (complex<double>(real, imaginary))*root[j];
            
        
            }
            
        }

        
         
        
        generate_data_parallely(sample_array, fact);

    
        
        compute_moments(sample_array, expectation_sum_each_node);

        
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


    

    double expectation_sum[N+1] = {0.0};

    double *global_expectation_sum_array = NULL;

    double local_expectation_sum_array[N+1];
     
    


    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    
    
    if(rank == 0){
        t_1 = std::time(0);
        cout << "number of MPI processes: " << size << endl;
        cout << "number of maximum possible threads: "<< omp_get_max_threads()<< endl;
        global_expectation_sum_array = new double[(N+1)*size];
        
        
    }

        /*MPI_Scatter distributes the computation accross all nodes. */

    MPI_Scatter(global_expectation_sum_array, N+1, MPI_DOUBLE, local_expectation_sum_array, N+1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

/*The following function operates in all communicating nodes. As seen above, this function calls all othe functions in this program.*/  
    generate_distribution(root, fact, local_expectation_sum_array, size, rank);

 /*MPI_Gather gathers data from all communicating nodes.*/

    MPI_Gather(local_expectation_sum_array, N+1, MPI_DOUBLE, global_expectation_sum_array, N+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Finalize(); 

    if(rank == 0){

        for(int i =0 ; i < size; i++){
            for(int j =0; j< N+1; j++){
                expectation_sum[j] += global_expectation_sum_array[i*(N+1)+j];
            }

        } 

        delete [] global_expectation_sum_array;
        write_expectation(expectation_sum);

        t_2 = std::time(0);
        cout << "total time of completion: " << t_2-t_1 << endl;
     

    }

   
   return 0;

}