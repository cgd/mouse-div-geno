#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define UNKNOWN_GENO_CODE -1

/**
 * calculates all pairwise euclidian distance between vector pairs and returns
 * the result as a lower triangle. The length of the returned matrix is
 * (vectors_length - 1) with the sub-arrays starting with length = 1 and
 * increasing by 1 to (vectors_length - 1). Note that the diagonal (which would
 * be all 0's) is not included
 */
double** euclidian_distances(
    int vectors_length,
    double vec_data1[],
    double vec_data2[])
{
    double **dist_triangle = calloc(sizeof(double*), vectors_length - 1);
    for(int row = 1; row < vectors_length; row++)
    {
        double *curr_row = calloc(sizeof(double), row);
        dist_triangle[row - 1] = curr_row;
        for(int col = 0; col < row; col++)
        {
            double diff1 = vec_data1[row] - vec_data1[col];
            double diff2 = vec_data2[row] - vec_data2[col];
            
            curr_row[col] = sqrt(diff1 * diff1 + diff2 * diff2);
        }
    }
    
    return dist_triangle;
}

/**
 * Frees the arrays allocated by euclidian_distances
 */
void free_euclidian_distances(int vectors_length, double **dist_triangle)
{
    for(int row = 1; row < vectors_length; row++)
    {
        free(dist_triangle[row - 1]);
    }
    free(dist_triangle);
}

/**
 * this function allows you to treat lower_triangle as a symetric matrix rather
 * than a triangle (so you can swap i1 and i2 and should get the same result).
 * It is not leagal to ask for a number along the diagonal (where i1 == i2).
 */
double symmetric_get(double **lower_triangle, int i1, int i2)
{
    // it is not valid to ask for something along the diagonal
    assert(i1 != i2 && i1 >= 0 && i2 >= 0);
    
    int index1;
    int index2;
    if(i1 > i2)
    {
        index1 = i1 - 1;
        index2 = i2;
    }
    else
    {
        index1 = i2 - 1;
        index2 = i1;
    }
    
    return lower_triangle[index1][index2];
}

/**
 * This function finds the missing genotypes and assigns them by using a nearest
 * neighbor approach.
 */
void vdist(int vectors_length, double vec_data1[], double vec_data2[], int genotypes[])
{
    // which genotypes are unknown?
    int assigned_geno_count = 0;
    int assigned_geno_indices[vectors_length];
    int unassigned_geno_indices[vectors_length];
    for(int i = 0; i < vectors_length; i++)
    {
        if(genotypes[i] == UNKNOWN_GENO_CODE)
        {
            unassigned_geno_indices[i - assigned_geno_count] = i;
        }
        else
        {
            assigned_geno_indices[assigned_geno_count] = i;
            assigned_geno_count++;
        }
    }
    
    // if all genotypes are unknown then there is nothing that we can do
    if(assigned_geno_count == 0 || vectors_length - assigned_geno_count == 0)
    {
        return;
    }
    
    // calculate all pairwise distances
    double **dist_triangle = euclidian_distances(vectors_length, vec_data1, vec_data2);
    
    // for each iteration we will find the unassigned genotype which is nearest
    // to an assigned genotype. We then set the unassigned genotype to match
    // the assigned genotype
    for( ; assigned_geno_count < vectors_length; assigned_geno_count++)
    {
        int unassigned_geno_count = vectors_length - assigned_geno_count;
        
        int min_assigned_index = assigned_geno_indices[0];
        int min_unassigned_index_index = 0;
        int min_unassigned_index = unassigned_geno_indices[min_unassigned_index_index];
        double min_dist = symmetric_get(dist_triangle, min_assigned_index, min_unassigned_index);
        
        for(int unassigned_index_index = 0; unassigned_index_index < unassigned_geno_count; unassigned_index_index++)
        {
            int unassigned_index = unassigned_geno_indices[unassigned_index_index];
            for(int assigned_index_index = 0; assigned_index_index < assigned_geno_count; assigned_index_index++)
            {
                int assigned_index = assigned_geno_indices[assigned_index_index];
                double curr_dist = symmetric_get(dist_triangle, assigned_index, unassigned_index);
                if(curr_dist < min_dist)
                {
                    min_assigned_index = assigned_index;
                    min_unassigned_index_index = unassigned_index_index;
                    min_unassigned_index = unassigned_index;
                    min_dist = curr_dist;
                }
            }
        }
        
        // now that we have the min pair, assign the genotype and update
        // the assigned/unassigned indices
        genotypes[min_unassigned_index] = genotypes[min_assigned_index];
        assigned_geno_indices[assigned_geno_count] = min_unassigned_index;
        for(int i = min_unassigned_index_index; i < (unassigned_geno_count - 1); i++)
        {
            unassigned_geno_indices[i] = unassigned_geno_indices[i + 1];
        }
    }
    
    free_euclidian_distances(vectors_length, dist_triangle);
}

void vdist_from_r(int *vectors_length, double vec_data1[], double vec_data2[], int genotypes[])
{
    vdist(*vectors_length, vec_data1, vec_data2, genotypes);
}
