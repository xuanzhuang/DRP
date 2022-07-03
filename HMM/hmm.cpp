#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>  //for high_resolution_clock

using namespace std;

/* 
    HMM main function overview
    (1) read input data
    (2) HMM forward equations (alphas)
    (3) HMM backward equations (betas)
    (4) calculate the probabilities and output results
*/

/* constant variable definition */
const int n_founders = 8; //number of founders
const int n_states = 36; //number of states
const double trans_scale = 1;
const double bb_alpha = 0.5; //used for beta binomial
const double bb_beta = 0.5; //used for beta binomial
const double err = 0.005; // Genotyping error rate
const double max_prob_threshold = 0.95;

struct two_states {
    int a;
    int b;
};

const string founder[n_founders] = {"F1","F2","F3","F4","F5","F6","F7","F8"};

const string founder_state[n_states] = {"F1F1","F1F2","F1F3","F1F4","F1F5","F1F6","F1F7","F1F8","F2F2","F2F3","F2F4","F2F5","F2F6","F2F7","F2F8","F3F3","F3F4","F3F5","F3F6","F3F7","F3F8","F4F4","F4F5","F4F6","F4F7","F4F8","F5F5","F5F6","F5F7","F5F8","F6F6","F6F7","F6F8","F7F7","F7F8","F8F8"};
const string additive_column[n_founders] = {"FF1","FF2","FF3","FF4","FF5","FF6","FF7","FF8"};

const two_states founder_two_states[36] = {
    {0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7},
    {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7},
    {2, 2}, {2, 3}, {2, 4}, {2, 5}, {2, 6}, {2, 7},
    {3, 3}, {3, 4}, {3, 5}, {3, 6}, {3, 7},
    {4, 4}, {4, 5}, {4, 6}, {4, 7},
    {5, 5}, {5, 6}, {5, 7},
    {6, 6}, {6, 7},
    {7, 7}
};

// data structure for each record, used by the forward and backward functions
struct record {
    string arm;
    long pos;
    int K;
    int N;
    double pb[n_founders];
    double rf;
};

/* function definition */
void emission(int K, int N, double (&pb)[n_founders], double (&emission_p)[n_states]);
void transition(double rf, double (&transition_p)[n_states][n_states]);
double elnsum(double x, double y);
void print_transition_p(double (&transition_p)[n_states][n_states]);
void print_emission_p(double (&emission_p)[n_states]);


/* The first argument is the sample file name */
int main(int argc, char** argv)
{
    if(argc!=4)
    {
        cout << "hmm takes three inputs. Refer to snakemake for details." << endl;
        return 0;
    }
    
    // Record the start time
    auto start = chrono::high_resolution_clock::now();
    
    ifstream infile(argv[1]);
    string outfile_name = argv[2];
    string outfile_plot_name = argv[3];
    
    long file_size = 0;
    
    long i = 0; //index of sample file record

    double emission_p[n_states] = {};
    
    double transition_p[n_states][n_states];
    
    double init_p[36] = {1.0/64, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32,
                         1.0/64, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32,
                         1.0/64, 1.0/32, 1.0/32, 1.0/32, 1.0/32, 1.0/32,
                         1.0/64, 1.0/32, 1.0/32, 1.0/32, 1.0/32,
                         1.0/64, 1.0/32, 1.0/32, 1.0/32,
                         1.0/64, 1.0/32, 1.0/32,
                         1.0/64, 1.0/32,
                         1.0/64};
    
    // Step 1: read the file line-by-line into a vector, which will be used by the forward and backward functions
    vector<record> data;

    if (infile.is_open())
    {
        while(!infile.eof())
        {
            data.push_back(record());

            if(infile >> data[i].arm >> data[i].pos >> data[i].K >> data[i].N >> data[i].pb[0] >> data[i].pb[1] >> data[i].pb[2] >> data[i].pb[3] >> data[i].pb[4] >> data[i].pb[5] >> data[i].pb[6] >> data[i].pb[7] >> data[i].rf)
            {
            
                //cout << data[i].arm << " " << data[i].pos << " " << data[i].K << " " << data[i].N << " " << data[i].pb[0] << " " << data[i].pb[1] << " " << data[i].pb[2] << " " << data[i].pb[3] << " " << data[i].pb[4] << " " << data[i].pb[5] << " " << data[i].pb[6] << " " << data[i].pb[7] << " " << data[i].rf << endl;
                i++;
            }
        }
        
        file_size = i;
        cout << "Read in " << file_size << " records." << endl;
        infile.close();
    }
    else
    {
        cout << "Unable to open file " << argv[1] << endl;
        return 0;
    }
    
    if(file_size == 0) return 0;
    
    vector< vector<double> > alpha(file_size, vector<double>(n_states, 0));
    
    vector< vector<double> > beta(file_size, vector<double>(n_states, 0));
        

    // Step 2: HMM forward equations (alphas)
    for (i=0; i<file_size; i++)
    {
        // for the first position
        if (i == 0)
        {
            emission(data[i].K, data[i].N, data[i].pb, emission_p);
            
            for (int j=0; j<n_states; j++)
            {
                alpha[i][j] = init_p[j] + emission_p[j];
            }
        }
        else
        {
            // Get the transition probabilities for this position and the previous position
            transition(data[i-1].rf, transition_p);
                        
            // Get the emission probabilities for this position
            emission(data[i].K, data[i].N, data[i].pb, emission_p);
            
            for (int x=0; x<n_states; x++)
            {
                double sum = 0;
                
                // add up values for all the founders
                for (int y=0; y<n_states; y++)
                {
                    // note that the alpha is already eln from above
                    double a1_trans = alpha[i-1][y] + transition_p[x][y];
                    sum = elnsum(sum, a1_trans);
                }
                
                alpha[i][x] = sum + emission_p[x];

            }
        }
    }
    
    // Step 3: HMM backward equations (betas), count backwards
    for (i=file_size-1; i>=0; i--)
    {
        // for the last position
        if (i == file_size-1)
        {
            for (int j=0; j<n_states; j++)
            {
                beta[i][j] = 0; //log(1)
            }
        }
        else
        {
            // Get transition for position and position ahead
            transition(data[i].rf, transition_p);
                
            // Get emission for position ahead
            emission(data[i+1].K, data[i+1].N, data[i+1].pb, emission_p);
                
            for (int x=0; x<n_states; x++)
            {
                double sum =0;
                
                // add up values for all the founders
                for (int y=0; y<n_states; y++)
                {
                    double b_tm1 = beta[i+1][y] + transition_p[x][y] + emission_p[y];
                    sum = elnsum(sum, b_tm1);
                }
                    
                beta[i][x] = sum;

            }
        }
    }
    
    //Step 4: Compute the probability of state for each position; Get Y for each position and generate outputs
    
    ofstream out_plot_file (outfile_plot_name);
    if (!out_plot_file.is_open())
    {
        cout << "Error: Unable to open file " << outfile_plot_name << endl;
    }
    
    ofstream out_file (outfile_name);
    if (!out_file.is_open())
    {
        cout << "Error: Unable to open file " << outfile_name << endl;
    }
	
    //print output file header
    out_file << "index,chr,pos,";
    for (int x=0; x<n_states; x++)
    {
        out_file << founder_state[x] << ",";
    }
    for (int x=0; x<n_founders; x++)
    {
        out_file << founder[x] << ",";
    }
    out_file << endl;
    
    // print output file body
    for (i=0; i<file_size; i++)
    {
        out_file << i << "," << data[i].arm << "," << data[i].pos << "," ;
		
        double addsum[n_founders] = {0};
        
        double max_y_value = 0;
        int max_state = -1;
		
        for (int x=0; x<n_states; x++)
        {
            double sum =0;
                
            for (int y=0; y<n_states; y++)
            {
                double value = alpha[i][y] + beta[i][y];
                sum = elnsum(sum, value);
            }
            
            // Get Y out of log space
            double Y_value = exp( alpha[i][x] + beta[i][x] - sum );
            
            if(Y_value > max_y_value)
            {
                max_y_value = Y_value;
                max_state = x;
            }
            
            for (int j=0; j<n_founders; j++)
            {
                int found1 = founder_two_states[x].a;
                int found2 = founder_two_states[x].b;
                
                if ((found1 == j) && (found2 == j))
                {
                    addsum[j] = addsum[j] + Y_value;
                    
                }
                else if ((found1 == j) || (found2 == j))
                {
                    addsum[j] = addsum[j] + 0.5*Y_value;
                }
            }
            
            out_file << Y_value << ",";
        }
        
        if (max_y_value >= max_prob_threshold)
        {
            out_plot_file << data[i].arm << "," << data[i].pos << "," << founder_state[max_state] << "," << founder_two_states[max_state].a << "," << founder_two_states[max_state].b << endl;
        }
        
        //output 8 additives
        for (int j=0; j<n_founders; j++)
        {
            out_file << addsum[j] << ",";
        }
        
        out_file << endl;
    }
	
    out_file.close();
    out_plot_file.close();
    
    // Record end time
    auto finish = chrono::high_resolution_clock::now();
    
    chrono::duration<double> elapsed = finish - start;
    cout << "Elapsed time: " << elapsed.count() << " s\n";
    
	return 0;
}

void emission(int K, int N, double (&pb)[n_founders], double (&emission_p)[n_states])
{
    double p1 = 0;
    double p2 = 0;
    double prob_sum = 0;
    double binom_coef = 0;
    
    if (N<=0)
    {
        cout << "Warning: N <= 0 !" << endl; // error checking; should not happen since we remove N=0 in the inputs
        
        for (int x=0; x<n_states; x++)
        {
            emission_p[x] = 0; //? correct?
        }
        return;
    }
    
    for (int x=0; x<n_states; x++)
    {
        p1 = pb[founder_two_states[x].a];
        p2 = pb[founder_two_states[x].b];
        
        //add up all the possible genotypes
        double v1 = exp( lgamma(bb_alpha+K) + lgamma(N+bb_beta-K) - lgamma(bb_alpha+N+bb_beta) ); //calculate beta function using gamma function 
        double v2 = exp( lgamma(bb_alpha) + lgamma(bb_beta) - lgamma(bb_alpha+bb_beta) );
        prob_sum = log( (p1 * p2 * pow(1-err, K) * pow(err, N-K)) + ((1-p1) * p2 * v1/v2) + (p1 * (1-p2) * (v1/v2)) + ((1-p1) * (1-p2) * pow(err, K) * pow(1-err, N-K)) );
        
        // get binomial coefficient
        binom_coef = lgamma(N+1) - lgamma(K+1) - lgamma(N-K+1);
        
        // calculate probability
        emission_p[x] = binom_coef + prob_sum;
    }
    
    return;
}

void transition(double rf, double (&transition_p)[n_states][n_states])
{

    //computer 1 - e^(-r)
    double trans_prob = 1 - exp(-rf*trans_scale);
    
    for (int x=0; x<n_states; x++)
    {
        for (int y=0; y<n_states; y++)
        {
            if (x == y) //founder state at position 1 = founder state at position 2
            {
                // Probability of staying in state
                if(founder_two_states[x].a == founder_two_states[x].b) //case ii -> ii
                {
                    transition_p[x][y] = log( (1 - trans_prob + trans_prob/n_founders) * (1 - trans_prob + trans_prob/n_founders) );
                }
                else //case ij -> ij
                {
                    transition_p[x][y] = log( (1 - trans_prob + trans_prob/n_founders) * (1 - trans_prob + trans_prob/n_founders) + (trans_prob/n_founders) * (trans_prob/n_founders) );
                }
            }
            else //founder state at position 1 != founder state at position 2
            {
                //If original is homozygous:
                if (founder_two_states[x].a == founder_two_states[x].b) //position 1 = ii
                {
                    // and if new is homozygous
                    if (founder_two_states[y].a == founder_two_states[y].b) //position 2 = jj (and position 1 = ii), i.e., ii -> jj
                    {
                        // Probability of moving to another homozygous state
                        transition_p[x][y] = log( trans_prob/n_founders * trans_prob/n_founders );
                    }
                    else // position 2 != jj (and position 1 = ii)
                    {
                        // Original homozygous, new is heterozygous with one copy same as original
                        if ( (founder_two_states[x].a == founder_two_states[y].a) || (founder_two_states[x].a == founder_two_states[y].b) ) // ii -> ij or ii -> ji
                        {
                            transition_p[x][y] = log( 2 * (1 - trans_prob + trans_prob/n_founders) * (trans_prob/n_founders) );
                        }
                        else // ii -> jk
                        {
                            // All other moves from homozygous state involve 2 events, prob near 0
                            transition_p[x][y] = log( 2 * trans_prob/n_founders * trans_prob/n_founders );                     
                        }
                    }
                }
                else //position 1 !=ii
                {
                    // If original heterozygous, new is one of homo states (only btw two represented in het)
                    if ( ((founder_two_states[y].a == founder_two_states[x].a)
                          && (founder_two_states[y].b == founder_two_states[x].a))
                         || ((founder_two_states[y].a == founder_two_states[x].b)
                          && (founder_two_states[y].b == founder_two_states[x].b)) ) //ij -> ii or ij->jj
                    {
                        transition_p[x][y] = log( (1-trans_prob + trans_prob/n_founders) * (trans_prob/n_founders) );
                    }
                    else
                    {
                        // Original heterozygous, new is het with one founder same as original het
                        if ((founder_two_states[y].a == founder_two_states[x].a) ||
                            (founder_two_states[y].a == founder_two_states[x].b) ||
                            (founder_two_states[y].b == founder_two_states[x].a) ||
                            (founder_two_states[y].b == founder_two_states[x].b)) //ij -> ik or jk or ki or kj
                        {
                            transition_p[x][y] = log( (1-trans_prob + trans_prob/n_founders) * (trans_prob/n_founders) + (trans_prob/n_founders) * (trans_prob/n_founders) );
                        }
                        else
                        {
                            if (founder_two_states[y].a == founder_two_states[y].b) // ij -> kk
                            {
                                transition_p[x][y] = log( trans_prob/n_founders * trans_prob/n_founders );
                            }
                            else //ij -> kl
                            {
                                transition_p[x][y] = log( 2 * trans_prob/n_founders * trans_prob/n_founders );
                            }
                        }
                    }
                }
            }
            
        } // end for (int y=0; y<n_states; y++)

    } // end (int x=0; x<n_states; x++)
    
    return;
}

double elnsum(double x, double y)
{
    // x, y are ln based probabilities

    if(x==0)
    {
        return y;
    }
    if(y==0)
    {
        return x;
    }
    if (x > y)
    {
        return x + log(1 + exp(y-x));
    }
    else
    {
        return y + log(1 + exp(x-y));
    }
}

void print_transition_p(double (&transition_p)[n_states][n_states])
{
    cout << "transition_p" << endl;
    
    for (int x=0; x<n_states; x++)
    {
        for (int y=0; y<n_states; y++)
        {
            cout << transition_p[x][y] << " ";
        }
        
        cout << endl;
    }
    cout << endl;
    
    return;
}

void print_emission_p(double (&emission_p)[n_states])
{
    cout << "emission_p" << endl;
    
    for (int x=0; x < n_states; x++)
    {
        cout << emission_p[x] << " ";
    }

    cout << endl;
    cout << endl;
    return;
}


