#include <iostream> //cout, endl
#include <fstream> // file read
#include <limits> // double max
#include <stdlib.h> // Random
#include <assert.h> // Assert
#include <iomanip>
#include <vector>
#include <iomanip>

#include "basic.h"
#include "kmeansDP.h"

// #define DEBUG 1
using std::setw;
using std::cout; using std::endl;
using std::ifstream;
using std::left;
using std::string;
using std::vector;
using std::setprecision;

// Prototypes
void readData(const char* path, points *&data, double *&domain);
void performExp(points *data, double* domain, int k, double epsilon, int seed , double* result, int maxH, double T, double ratio, int n_sample);

int main(int argc, char* args[]){
  const int print_width = 9;
	const int iterExp = 100;
	cout << "iterExp = " << iterExp << endl;
  string experiment_mode;
  // const string experiment_mode = "n_sample";
	// Load data
	points *data;
	double* domain;
	if(argc<2)	readData("test", data, domain);
	else	readData(args[1], data, domain);
	if (argc < 3) {
	  cout << "invalid arguments !!!!!!!!" << endl;
    std::exit(1);
  }
  else{
    experiment_mode = args[2];
  }
  cout << "experiment_mode = " << experiment_mode << endl;
	
	// Check the domain of the data
	for(int i = 0; i< data->dim; ++i)
		cout<<i<<"-th dim:["<<domain[2*i]<<", "<<domain[2*i+1]<<"]"<<endl;
	cout<<endl;

	int sizeExp = 100;
	double *result = new double[sizeExp];
	
	// Experiment parameter
	int seed = 2;
	int n_cluster_min = 5;
	int n_cluster_max = 6;
	int n_cluster_step = 1;

	// Tree
	double T_min = data->size/1000.;   // min-threshold
	double T_max = T_min+1;

  double epsilon_max, epsilon_min, ratio_min, ratio_max;
  int min_n_sample = 30, max_n_sample = 30;
	// ratio for tree & noise
	if (experiment_mode == "ratio"){
    if (false) {
      if (argc < 4) {
        cout << "invalid arguments !!!!!!!!" << endl;
        std::exit(1);
      }
      epsilon_min = std::stod(args[3]); // default: 0.05
      epsilon_max = std::stod(args[3]);
    }
    epsilon_min = 0.05;
    epsilon_max = 2;
    ratio_min = 0.05;
    ratio_max = 1.0;
  }
  else if (experiment_mode == "eps"){
    epsilon_min = 0.05;
    epsilon_max = 2;
    ratio_min = 0.15;
    ratio_max = 0.16;
  }
  else if (experiment_mode == "n_sample"){
    epsilon_min = 0.05;
    epsilon_max = 0.06;
    ratio_min = 0.15;
    ratio_max = 0.16;
    min_n_sample = 1;
    max_n_sample = 30;
  }
  else{
    cout << "invalid experiment_mode !!!!!!!!!!!!! " << endl;
  }

	//int minH = 4*log(data->size)/(log(t)*data->dim);
	int minH = log(data->size)/(data->dim)-1;  // depth
	//int minH = 3;
	int maxH = minH+1;


	// Perform Experiment	
	cout<<"kmeans vs NaiveDP vs EUGDP vs Proposed (uni, rel)"<<endl;
	vector<string> column_names = {"k", "ratio", "eps", "maxH", "T", "n_sample", "uni", "rel", "vs_Eu", "vs_uni", "toNoNoise"};
	for (int i = 0; i<column_names.size();i++){
	  cout << setw(print_width) << left << column_names[i];
  }
  cout << endl;
  for (int n_sample = min_n_sample; n_sample <= max_n_sample; n_sample ++){
    for(double ratio = ratio_min; ratio<ratio_max; ratio+=0.05){
      for(int t =T_min; t < T_max; t*= 2){
        for (int h = minH; h<maxH; ++h){
          for (double epsilon = epsilon_min; epsilon <= epsilon_max; epsilon *= 2){
            for(int k= n_cluster_min; k < n_cluster_max; k+= n_cluster_step){

              double per1, per2, per3, per4, percent_of_uni_quad_to_the_no_noise;
              double avg1 = 0.0, avg2 = 0.0, avg3=0., avg4 = 0., avg5 = 0.;
              for(int i = 0; i< iterExp;++i){
                // Initialize the result
                for(int temp = 0;temp<sizeExp;++temp)	result[temp]=0.;
                
                performExp(data, domain, k, epsilon, seed+i, result, h, t, ratio, n_sample);
                double uniform_quadtree = result[8];
                if (false){
                  cout << "no_privacy: " << result[0] << "   " <<
                      "existing_work: " << result[2] << "   " <<
                      "uni: " << result[3] << "   " <<
                      "rel: " << result[4] << "\n"
                      << "uniQuad: " << uniform_quadtree << endl;
                }
                // relative increments. The worse one must be preceded.
                per1 = ((result[2] - result[3]) / result[0] * 100.); //uni
                per2 = ((result[2] - result[4]) / result[0] * 100.); //rel
                per3 = ((result[2] - uniform_quadtree) / result[0] * 100.); //rel
                per4 = ((result[3] - uniform_quadtree) / result[0] * 100.); //rel
                percent_of_uni_quad_to_the_no_noise = ((uniform_quadtree - result[0]) / result[0] * 100.); //rel
                avg1 += per1;
                avg2 += per2;
                avg3 += per3;
                avg4 += per4;
                avg5 += percent_of_uni_quad_to_the_no_noise;
                // cout<<per1<<", "<<per2<<" ";
              }
              avg1 /= iterExp;
              avg2 /= iterExp;
              avg3 /= iterExp;
              avg4 /= iterExp;
              avg5 /= iterExp;

              // Output the result
              cout << setw(print_width) << left << std::fixed << setprecision(3) << k;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << ratio;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << epsilon;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << h;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << t;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << n_sample;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << avg1;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << avg2;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << avg3;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << avg4;
              cout << setw(print_width) << left << std::fixed << setprecision(3) << avg5;
              cout << endl;
              
              // cout<<"k = "<<k<<" ratio = "<<ratio<<" eps = "<<epsilon<<" maxH = "<< h <<" T = "<<setw(3)<<t<<"  ";
              /*
              for(int i = 0;i<sizeExp-3;++i)
                cout<< setw(12) << result[i]/iterExp <<" ";
              cout<<" EUG "<<setw(5)<<result[5]/iterExp <<" K_ "<<setw(5)<< result[6]/iterExp;
              double per1 = ((result[2] -result[3])/result[0] * 100.); //uni
              int per2 = int((result[2] -result[4])/result[0] * 100.); //rel
              */
              // cout<<"   uni_total: "<<setw(3)<<avg1<<"  "<<"rel_total: "<<setw(3)<<avg2<<"  "<< "uniform_quad: "<<setw(3) << avg3 << "   vs_uni:" << setw(3) << avg4 << endl;

            }
          }
        }
      }
    }
  }

	// de-allocation
	delete[] result;
	delete[] domain;
	delete data;

	cout<<"    Done!   "<<endl;
	return 0;
}

void performExp(points *data, double* domain, int k, double epsilon, int seed, double* result, int maxH, double T, double ratio, int n_sample){

	points *initial = new points(k, data->dim);
	initCenter(initial, domain, data->dim, seed);
#ifdef  DEBUG
	cout<<"Initial Center    ";
	printPoints(initial);
	cout<<endl;
#endif
	// For k-means result
	points *center1 = new points(k, data->dim);

	// Perform k-means clusterings
	double SSE1 = kmeans(data, initial, center1);
#ifdef  DEBUG
	cout<<"Kmeans 1    ";
	printPoints(center1);	
	cout<<endl;
#endif

	double SSE2 = kmeansNaiveDP(data, initial, center1, domain, epsilon);
#ifdef  DEBUG
	cout<<"Kmeans 2    ";
	printPoints(center1);
	cout<<endl;
#endif

	int n_bucket = kmeansEUGDP(data, initial, center1, domain, epsilon);
	double SSE3 = calcSSE(data, center1);
#ifdef  DEBUG
	cout<<"Kmeans 3    ";
	printPoints(center1);	
	cout<<endl;
#endif

	int n_bucket1 = kmeansDP(data, initial, center1, domain, epsilon, seed, maxH, T, 1, ratio);
	double SSE4 = calcSSE(data, center1);
	#ifdef  DEBUG
	cout<<"Kmeans 4    ";
	printPoints(center1);
	cout<<endl;
#endif
	int n_bucket2 = kmeansDP(data, initial, center1, domain, epsilon, seed, maxH, T, 2, ratio);
	double SSE5 = calcSSE(data, center1);
#ifdef  DEBUG
	cout<<"Kmeans 5    ";
	printPoints(center1);
	cout<<endl;
#endif

	// run algorithm that uniform points are in each leaf of the quadtree
	int n_bucket3 = uniform_kmeansDP(data, initial, center1, domain, epsilon, seed, maxH, T, 1, ratio, n_sample);
	double SSE6 = calcSSE(data, center1);
#ifdef  DEBUG
	cout<<"uniform    ";
	printPoints(center1);
	cout<<endl;
  cout << endl<<endl;
#endif
	// Evaluation
	//cout<<"SSE @ k = "<<k<<" "<<SSE1<<", "<<SSE2<<", "<<SSE3<<", "<<SSE4<<endl;		
	result[0] += SSE1 / data->size;
	result[1] += SSE2 / data->size;
	result[2] += SSE3 / data->size;   // kmeansEUGDP
	result[3] += SSE4 / data->size;   // uni-quad
	result[4] += SSE5 / data->size;   // rel-quad
	result[5] += n_bucket;
	result[6] += n_bucket1;
	result[7] += n_bucket2;
	result[8] += SSE6 / data->size;   // the result of the uniform_quadtree
	result[9] += n_bucket3;
	delete initial;
	delete center1;
}

void readData(const char* path, points *&data, double *&domain){
	cout<<"Reading the data! :  "<<path<<endl;
	ifstream inputFile(path);
	if(inputFile.is_open()){
		int size, dim, count = 0;
		inputFile>>size;
		inputFile>>dim;
		data = new points(size, dim);
		domain = new double[dim*2];
		for(int i = 0; i<dim;++i){
			domain[2*i]=DBL_MAX;
			domain[2*i+1]=DBL_MIN;
		}
		// Data load
		while(inputFile>>data->d[count]){
			domain[2*(count%dim)] = MIN(data->d[count], domain[2*(count%dim)]);
			domain[2*(count%dim)+1] = MAX(data->d[count], domain[2*(count%dim)+1]);
			++count;
		}
		double area = 1.0;
		for(int i = 0; i<dim;++i)
			area *= (domain[2*i+1]-domain[2*i]);
		cout<<"Size, Dim = "<<size<<", "<<dim<<endl;
		cout<<"Density = "<<size/area<<endl;
		cout<<endl;
		inputFile.close();
	}
	else
		cout<<"Unable to open file: "<<path<<endl;
	return;
}
