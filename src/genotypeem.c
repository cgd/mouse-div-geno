#include <Rmath.h>
#include <stdlib.h>

#define EM_THRESHOLD 0.001
#define EM_ITER_MAX 51

/**
 * calculates expectation values for the possible genotypes
 * @param data_vec
 *          a vector of data that we use to calculate expectation. The
 *          vector length should be equal to sample_count
 * @param sample_count
 *          the number of samples
 * @param geno_count
 *          determines the number of genotypes whose probability we
 *          are estimating
 * @param taus
 *          TODO fill in doc
 * @param mus
 *          TODO fill in doc
 * @param sigmas
 *          TODO fill in doc
 * @param expectation_matrix
 *          this matrix is destructively updated with the expectation that
 *          a data item belongs to a genotype. These probabilities can be
 *          accessed like expectation_matrix[sample_index + geno_index * sample_count]
 */
void calc_expectation(
    double *data_vec,
    int sample_count,
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
    
    for(int data_index = 0; data_index < sample_count; data_index++)
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
                expectation_matrix[data_index + geno_index * sample_count] = curr_val;
            }
            else
            {
                expectation_matrix[data_index + geno_index * sample_count] = 0.0;
            }
        }
        
        // make sure all of the genos sum up to 1.0 so that they can be
        // interpreted as probabilities
        for(int geno_index = 0; geno_index < geno_count; geno_index++)
        {
            expectation_matrix[data_index + geno_index * sample_count] /= curr_sum;
        }
    }
}

void calc_expectation_two_genos_from_r(
    double *data_vec,
    int *sample_count,
    double *taus,
    double *mu1,
    double *mu2,
    double *sigma1,
    double *sigma2,
    double *expectation_matrix)
{
    double mus[2] = {*mu1, *mu2};
    double sigmas[2] = {*sigma1, *sigma2};
    
    calc_expectation(
        data_vec,
        *sample_count,
        2,
        taus,
        mus,
        sigmas,
        expectation_matrix);
}

void calc_expectation_three_genos_from_r(
    double *data_vec,
    int *sample_count,
    double *taus,
    double *mu1,
    double *mu2,
    double *mu3,
    double *sigma1,
    double *sigma2,
    double *sigma3,
    double *expectation_matrix)
{
    double mus[3] = {*mu1, *mu2, *mu3};
    double sigmas[3] = {*sigma1, *sigma2, *sigma3};
    
    calc_expectation(
        data_vec,
        *sample_count,
        3,
        taus,
        mus,
        sigmas,
        expectation_matrix);
}

double weighted_mean(double *data, double *weights, int len)
{
    double sum_of_weighted_data = 0.0;
    double sum_of_weights = 0.0;
    for(int i = 0; i < len; i++)
    {
        sum_of_weighted_data += data[i] * weights[i];
        sum_of_weights += weights[i];
    }
    
    return sum_of_weighted_data / sum_of_weights;
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

/**
 * calculates the weighted variance of the given data
 * @param data  the data that we are calculating the variance for
 */
double weighted_variance(double *data, double *weights, int len)
{
    
    double sum_of_weights = 0.0;
    for(int i = 0; i < len; i++)
    {
        sum_of_weights += weights[i];
    }
    
    double scaled_weights[len];
    double sum_sq_scaled_weights = 0.0;
    double weighted_sum = 0.0;
    for(int i = 0; i < len; i++)
    {
        double curr_scaled_weight = weights[i] / sum_of_weights;
        scaled_weights[i] = curr_scaled_weight;
        
        weighted_sum += curr_scaled_weight * data[i];
        sum_sq_scaled_weights += curr_scaled_weight * curr_scaled_weight;
    }
    
    double xs[len];
    for(int i = 0; i < len; i++)
    {
        xs[i] = sqrt(scaled_weights[i]) * (data[i] - weighted_sum);
    }
    
    return self_crossprod(xs, len) / (1.0 - sum_sq_scaled_weights);
}

void maximize_expectation(
    double *expectation_matrix,
    double *data_vec,
    int sample_count,
    int geno_count,
    double *taus,
    double *mus,
    double *sigmas)
{
    for(int geno_index = 0; geno_index < geno_count; geno_index++)
    {
        taus[geno_index] = 0.0;
        for(int data_index = 0; data_index < sample_count; data_index++)
        {
            taus[geno_index] += expectation_matrix[data_index + geno_index * sample_count];
        }
        taus[geno_index] /= sample_count;
        
        mus[geno_index] = weighted_mean(
            data_vec,
            expectation_matrix + geno_index * sample_count,
            sample_count);
        
        sigmas[geno_index] = weighted_variance(
            data_vec,
            expectation_matrix + geno_index * sample_count,
            sample_count);
    }
}

/**
 * Simply copys from the from array to the to array
 */
void copy_doubles(double *from, double *to, int length)
{
    for(int i = 0; i < length; i++)
    {
        to[i] = from[i];
    }
}

/**
 * This function runs the EM algorithm on the given input. The contents of
 * expectation matrix and taus will be updated to contain the most likely
 * values as determined by the EM. The outher pointer values may also be
 * modified but thier contents are not defined
 * @param init_sigma
 *          initial value to use for sigma's TODO: more detail here
 * @param mus
 *          initial (hint) mean liklihood per genotype (the length of this array
 *          must be equal to geno_count). TODO: more detail here
 * @param expectation_matrix
 *          this is a (sample_count x geno_count) matrix with probabilities
 *          that a sample has a particular genotype. This means that the
 *          genotype rows should sum to one. The initial values passed in for
 *          this matrix do not matter.
 * @param data_vec
 *          the data per-sample (vector size should match sample_count)
 * @param taus
 *          length matches geno_count TODO: more detail here
 * @param sample_count
 *          the number of samples that we are genotyping in this EM call
 * @param geno_count
 *          the number of possible genotypes
 */
void run_em(
    double init_sigma,
    double *mus,
    double *expectation_matrix,
    double *data_vec,
    double *taus,
    int sample_count,
    int geno_count)
{
    // some initialization
    double init_tau = 1.0 / (double)geno_count;
    double sigmas[geno_count];
    for(int i = 0; i < geno_count; i++)
    {
        taus[i] = init_tau;
        sigmas[i] = init_sigma;
    }
    
    for(int i = 0, ok = 1; i < EM_ITER_MAX && ok; i++)
    {
        // copy prev values so that we can check the EM stopping condition
        double old_mus[geno_count];
        copy_doubles(mus, old_mus, geno_count);
        
        // the E step: updates the expectation_matrix
        calc_expectation(
            data_vec,
            sample_count,
            geno_count,
            taus,
            mus,
            sigmas,
            expectation_matrix);
        
        for(int geno_index = 0; geno_index < geno_count; geno_index++)
        {
            for(int sample_index = 0; sample_index < sample_count; sample_index++)
            {
                int curr_index = sample_index + geno_index * sample_count;
                if(isnan(expectation_matrix[curr_index]))
                {
                    expectation_matrix[curr_index] = init_tau;
                }
            }
        }
        
        // the M step: updates taus, mus and sigmas
        maximize_expectation(
            expectation_matrix,
            data_vec,
            sample_count,
            geno_count,
            taus,
            mus,
            sigmas);
        
        // determine whether or not we've hit the EM algo's stopping condition
        ok = 0;
        for(int geno_index = 0; geno_index < geno_count; geno_index++)
        {
            double tau = taus[geno_index];
            double mu = mus[geno_index];
            double old_mu = old_mus[geno_index];
            double sigma = sigmas[geno_index];
            
            if(isnan(sigma) || isinf(sigma))
            {
                ok = 0;
                break;
            }
            
            if(tau > 0 && fabs(mu - old_mu) >= EM_THRESHOLD)
            {
                ok = 1;
            }
        }
        
        if(ok)
        {
            double avg_sigma = 0.0;
            for(int geno_index = 0; geno_index < geno_count; geno_index++)
            {
                avg_sigma += sigmas[geno_index];
            }
            avg_sigma /= (double)geno_count;
            
            for(int geno_index = 0; geno_index < geno_count; geno_index++)
            {
                sigmas[geno_index] = (sigmas[geno_index] + avg_sigma) / 2.0;
            }
        }
    }
}

void run_em_from_r(
    double *init_sigma,
    double *mus,
    double *expectation_matrix,
    double *data_vec,
    double *taus,
    int *sample_count,
    int *geno_count)
{
    run_em(
        *init_sigma,
        mus,
        expectation_matrix,
        data_vec,
        taus,
        *sample_count,
        *geno_count);
}

void maximize_expectation_two_genos_from_r(
    double *expectation_matrix,
    double *data_vec,
    int *sample_count,
    double *taus,
    double *mus,
    double *sigmas)
{
    maximize_expectation(
        expectation_matrix,
        data_vec,
        *sample_count,
        2,
        taus,
        mus,
        sigmas);
}

void maximize_expectation_three_genos_from_r(
    double *expectation_matrix,
    double *data_vec,
    int *sample_count,
    double *taus,
    double *mus,
    double *sigmas)
{
    maximize_expectation(
        expectation_matrix,
        data_vec,
        *sample_count,
        3,
        taus,
        mus,
        sigmas);
}
