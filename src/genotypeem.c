#include <Rmath.h>

/**
 * calculates expectation values for the possible genotypes
 * @param data_vec
 *          a vector of data that we use to calculate expectation
 * @param data_len
 *          the size of the data vector
 * @param geno_count
 *          determines the number of genotypes whose probability we
 *          are estimating
 * @param taus
 *          TODO fill in doc
 * @param mus
            TODO fill in doc
 * @param sigmas
            TODO fill in doc
 * @param expectation_matrix
 *          this matrix is destructively updated with the expectation that
 *          a data item belongs to a genotype. These probabilities can be
 *          accessed like expectation_matrix[data_index + geno_index * data_len]
 */
void calc_expectation(
    double *data_vec,
    int data_len,
    int geno_count,
    double *taus,
    double *mus,
    double *sigmas,
    double *expectation_matrix)
{
    double sqrt_sigmas[geno_count];
    for(int geno_index = 0; geno_index < geno_count; geno_index++)
    {
        sqrt_sigmas[geno_index] = sqrt(sigmas[geno_index]);
    }
    
    for(int data_index = 0; data_index < data_len; data_index++)
    {
        double datum = data_vec[data_index];
        double curr_sum = 0.0;
        for(int geno_index = 0; geno_index < geno_count; geno_index++)
        {
            double tau = taus[geno_index];
            
            if(tau > 0.0)
            {
                double mu = mus[geno_index];
                double sqrt_sigma = sqrt_sigmas[geno_index];
                
                double curr_val = tau * dnorm(datum, mu, sqrt_sigma, 0);
                curr_sum += curr_val;
                expectation_matrix[data_index + geno_index * data_len] = curr_val;
            }
            else
            {
                expectation_matrix[data_index + geno_index * data_len] = 0.0;
            }
        }
        
        // make sure all of the genos sum up to 1.0 so that they can be
        // interpreted as probabilities
        for(int geno_index = 0; geno_index < geno_count; geno_index++)
        {
            expectation_matrix[data_index + geno_index * data_len] /= curr_sum;
        }
    }
}

void calc_expectation_two_genos_from_r(
    double *data_vec,
    int *data_len,
    double *tau_1_and_2,
    double *mu1,
    double *mu2,
    double *sigma1,
    double *sigma2,
    double *return_expectation_matrix)
{
    double mus[2] = {*mu1, *mu2};
    double sigmas[2] = {*sigma1, *sigma2};
    
    calc_expectation(
        data_vec,
        *data_len,
        2,
        tau_1_and_2,
        mus,
        sigmas,
        return_expectation_matrix);
}

/**
 * Calculates the cross product of a vector against itself
 * @param vec   the vector
 * @param len   the vector's length
 */
double self_crossprod(double vec[], int len)
{
    double result = 0.0;
    for(int i = 0; i < len; i++)
    {
        result += vec[i] * vec[i];
    }
    
    return result;
}
