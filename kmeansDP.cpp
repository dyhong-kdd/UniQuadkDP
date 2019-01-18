#include "kmeansDP.h"
#include "util_gen.h"
#include "dpqt.h"
#include <iostream>
#include <iomanip>

//#define EUGDP_DEBUG 1

using std::cout;
using std::endl;
using std::setw;

// SSE: the sum of squared errors
double calcSSE(const points *data, const points *center){
	double SSE = 0.0;
	double temp_dist, temp_min;
	for(int i = 0; i < data->size; ++i){
		// Find the closest center
		temp_min = DBL_MAX;
		for(int j = 0; j < center->size; ++j){
			temp_dist = dist(data->p[i], center->p[j], data->dim);
			if(temp_min > temp_dist)
				temp_min = temp_dist;
		}
		SSE += temp_min; // Assign Error
	}
	return SSE;
}

// For initialize the center
void initCenter(points *center, double *domain, int dim, int seed){
  static std::mt19937 gen(seed);

  double len = 0;
  for (int i = 0; i< dim; ++i){
    len += pow(domain[2*i+1]-domain[2*i],2);
  }
  len = pow(len, 0.5);
  //cout<<"len = "<<len<<endl;

  bool pass = true;
  while(true){
    for(int iter = 0; iter<100; ++iter){
    // sample once
      for(int i = 0; i< center->size; ++i){
          for(int j = 0; j<dim; ++j){
            std::uniform_real_distribution<double> uniDist(domain[2*j],domain[2*j+1]);
            double uni = uniDist(gen);
            center->p[i][j]= uni;
          }
        }
        double temp_dist = 0.;
        pass = true;
        for(int i = 0; i< center->size-1; ++i){
          for(int j = i+1; j< center->size; ++j){
            temp_dist = center->distPoint(i,j);
            //cout<<"dist"<<temp_dist<<endl;
            if(temp_dist < len)
              pass = false;
          }
        }
        if(pass==true){
        //cout<<"Pass"<<endl;
        return;
      }		
    }
    len *= 0.7;
    //cout<<"len = "<<len<<endl;
  }
}

void kmeansIter(const points *data, points *center, double *domain, double epsilon){
	points *sum = new points(center->size, center->dim);
	sum->fillZero();
	
	int idx = -1;
	double tempMin, tempDist;

	double *count = new double[center->size];
	for(int i = 0; i< center->size;++i)
		count[i] = 0;

	for(int i =0;i< data->size; ++i ){
		//Find the closest center (idx)
		tempMin = DBL_MAX;
		for(int j = 0; j<sum->size; ++j){
			tempDist = dist(data->p[i], center->p[j], data->dim);
			if(tempMin > tempDist){
				idx = j;
				tempMin = tempDist;
			}
		}
		// idx: the center index of the minimum distance with i-th point
		// Normally, weigths are just 1 unless it is Quad-Tree node or Histogram bucket.
		count[idx] += data->weights[i];
		sum->add(data, i, idx);
	}

	// Insert Noise
	// For 1. NaiveDP
	if(epsilon > 0.){
		for(int i = 0; i<sum->size; ++i){
			count[i] += laplace(0, 1./(0.5*epsilon));
			for(int j = 0; j<sum->dim; ++j){
				sum->p[i][j] += laplace(0, MAX(ABS(domain[2*j]),ABS(domain[2*j+1]))/(0.5*epsilon));
			}
		}
	}

	// Calcualte new center
	for(int i = 0; i< center->size;++i){
		if (count[i] <= 1) continue;
		for(int j = 0; j< center->dim;++j)
			center->p[i][j] = sum->p[i][j]/count[i];
	}

	delete sum;
	delete[] count;
}

// 0. No DP. Original K-means
double kmeans(const points *data, const points *initial, points *outCenter, int maxIter){
	points *center = new points(initial->size, initial->dim);
	center->copy(initial);
	double tempSSE = 0, minSSE = DBL_MAX;
	for (int i = 0; i< maxIter; ++i){
		// 1 iteration
		kmeansIter(data, center);
		tempSSE = calcSSE(data,center);
		if (tempSSE == minSSE)	break;
		minSSE = tempSSE;
	}
	outCenter->copy(center);

	delete center;
	return minSSE;
}

// 1. Naive DP
double kmeansNaiveDP(const points *data, const points *initial, points *outCenter, double *domain, double epsilon, int maxIter){
	points *center = new points(initial->size, initial->dim);
	center->copy(initial);

	for (int i = 0; i< maxIter; ++i)
		kmeansIter(data, center, domain, epsilon/maxIter);

	double SSE = calcSSE(data,center);
	outCenter->copy(center);
	delete center;
	return SSE;
}

// 2. EUG DP
// Devide the domain axis with M equal-width buckets
double kmeansEUGDP(const points *data, const points *center, points *outCenter, double *domain, double epsilon){

	// Histogram initialize
	double ratio = 0.00;
	int size = data->size;
	//int size = data->size + laplace(0, 1/(ratio*epsilon));
	
	int M =int(pow(size*epsilon/10., 2.0/(2.+data->dim)))+1;
	int n_bucket = int(pow(M, data->dim));
	double *histo = new double[n_bucket];

	for(int i = 0; i<n_bucket;++i) 
		histo[i]=0;

	// Make Histogram
	int idx = 0;
	for(int i =0;i< data->size; ++i ){
		idx = data->includeIn(i, M, domain);
		histo[idx]++;
	}
	
	// Debug
#ifdef EUGDP_DEBUG
	int sum_valid = 0;
	cout<<"EUGDP -> M = "<< M <<" #bucket = "<<n_bucket<<endl;
	cout<<"print buckets"<<endl;
	cout<<"---------------"<<endl;
	for(int i = 0;i < n_bucket; ++i ){
		if (i % M ==0 && i !=0 )	cout<<endl;
		sum_valid += histo[i];
		cout<<setw(5)<< histo[i]<<" ";
	}
	cout<<endl<<"data_size = "<< data->size<<" sum = "<<sum_valid<<endl;
	cout<<"---------------"<<endl;
#endif

	// Insert Noise in each bucket
	for(int i = 0; i< n_bucket; ++i)
		histo[i] += laplace(0, 1./((1-ratio)*epsilon)); // sensitivity is 1 in counting
	
	// Generate bucket points
	//cout<<"EUGDP:: generate bucket"<<endl;
	points *bucket = new points(n_bucket, data->dim);
	memcpy(bucket->weights, histo, sizeof(double) * n_bucket);
	for(int i = 0; i<n_bucket;++i){
		int idx = i;
		for(int j = 0; j<data->dim;++j){
			bucket->p[i][j] = domain[2*j] + (idx%M+0.5) * ((domain[2*j+1]-domain[2*j])/M);
			idx /= M;
		}
	}

	// Perform K-means over Histogram
	//cout<<"EUGDP:: kmeans"<<endl;
	double SSE_bucket = kmeans(bucket, center, outCenter);

	delete[] histo;
	delete bucket;
	
	return n_bucket;
}

// 3. Proposed Alg.
// Using Quad-Tree
double kmeansDP(const points *data, const points *center, points *outCenter, double *domain, double epsilon, int seed, int maxH, double T, int n_Type, double ratio){
	// Build QuadTree
	Noise noiseType;
	if(n_Type ==1)
		noiseType = Noise::uni; // uni
	else
		noiseType = Noise::rel; // geo, uni
	
	double** bound = new double*[2];
	bound[0] = new double[data->dim];
	bound[1] = new double[data->dim];
	for(int i = 0; i < data->dim;++i){
		bound[0][i] = domain[2*i];
		bound[1][i] = domain[2*i+1];
	}

	//double ratio = 1.0;

	DPQTTree tree(data->size, data->p, bound, T, data->dim, maxH, noiseType, ratio*epsilon, seed);

	// Find all leaves
	vector<DPQTNode*>* leaves = new vector<DPQTNode*>();
	tree.root->getLeaves(leaves);
	
	// Allocate Points to save leaf node
	int n_bucket = leaves->size();
	points *bucket = new points(n_bucket, data->dim);

	for(int i = 0; i< n_bucket;++i){
		double **bound =(*leaves)[i]->bound;
		for (int j = 0; j<data->dim;++j)
			bucket->p[i][j] = (bound[0][j]+bound[1][j])/2;
		bucket->weights[i] = (*leaves)[i]->count + laplace(0, 1./((1.-ratio)*epsilon) );
		// bucket->weights[i] = (*leaves)[i]->noisyCount;
	}

	if (false){

	}

	// Perform K-means
	double SSE_bucket = kmeans(bucket, center, outCenter);

	delete leaves;	
	delete[] bound[0];
	delete[] bound[1];
	delete[] bound;
	delete bucket;
	return n_bucket;
}

double uniform_kmeansDP(const points *data, const points *center, points *outCenter, double *domain, double epsilon, int seed, int maxH, double T, int n_Type, double ratio, int n_sample){
	// Build QuadTree
	const int div = 50, def_n_uniform_points = n_sample;

	Noise noiseType;
	if(n_Type ==1)
		noiseType = Noise::uni; // uni
	else
		noiseType = Noise::rel; // geo, uni

	double** bound = new double*[2];
	bound[0] = new double[data->dim];
	bound[1] = new double[data->dim];
	for(int i = 0; i < data->dim;++i){
		bound[0][i] = domain[2*i];
		bound[1][i] = domain[2*i+1];
	}

	//double ratio = 1.0;

	DPQTTree tree(data->size, data->p, bound, T, data->dim, maxH, noiseType, ratio*epsilon, seed);

	// Find all leaves
	vector<DPQTNode*>* leaves = new vector<DPQTNode*>();
	tree.root->getLeaves(leaves);

	// Allocate Points to save leaf node
	int n_bucket = leaves->size();
	points *bucket = new points(n_bucket, data->dim);

	int noisy_total_count = 0;
	for(int i = 0; i< n_bucket;++i){
		double **bound =(*leaves)[i]->bound;
		for (int j = 0; j<data->dim;++j)
			bucket->p[i][j] = (bound[0][j]+bound[1][j])/2;
		bucket->weights[i] = (*leaves)[i]->count + laplace(0, 1./((1.-ratio)*epsilon), seed);
    // noisy_total_count += std::abs(std::floor(bucket->weights[i] / div)) + 1;
    noisy_total_count += def_n_uniform_points ;
		// bucket->weights[i] = (*leaves)[i]->noisyCount;
	}

	points *uniform_points = new points(noisy_total_count, data->dim);
	int uniform_points_count = 0;
	for(int i = 0; i< n_bucket;++i){
		// double n_uniform_points = std::abs(std::floor(bucket->weights[i] / div)) + 1;
    double **bound =(*leaves)[i]->bound;
		double n_uniform_points = def_n_uniform_points;

		// cout<<n_uniform_points<<endl;

		for (int j=0; j<(int)n_uniform_points; j++){
			// insert uniform noise for each dimension
			for (int k=0; k<data->dim;k++){
				double bound_range = std::abs(bound[0][k] - bound[1][k]);
				double uniform_noise = generate_uniform_noise(seed) * bound_range;
				uniform_points->p[uniform_points_count][k] = (bound[0][k]+bound[1][k])/2 + uniform_noise;
			}
			uniform_points->weights[uniform_points_count] = bucket->weights[i] / n_uniform_points;

			// cout<<uniform_points->weights[uniform_points_count]<<endl;

			uniform_points_count ++ ;
			if (uniform_points_count > noisy_total_count){
				cout<<"exceed!!"<<endl;
			}
		}
	}
  // cout << "cap: " << noisy_total_count << " cnt: " << uniform_points_count << endl;
	// printPoints(uniform_points);

	// Perform K-means
	// cout << "my ======" << endl;
	double SSE_bucket = kmeans(uniform_points, center, outCenter);
	// cout << "end my ======" << endl;

	delete leaves;
	delete[] bound[0];
	delete[] bound[1];
	delete[] bound;
	delete bucket;
	return n_bucket;
}

