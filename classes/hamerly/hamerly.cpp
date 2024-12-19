#include "hamerly.h"
#include "../point/point.h"

#include <string>
#include <cmath>
#include <omp.h>

hamerly::hamerly( int _d, int _N, int _K ) : lloyd( _d, _N, _K )
{
    C = 0;

    maxdist = -1;
    secmaxdist = -1;
    maxdist_index = -1;

    cdistances = new float[K];
    criticals = NULL;

    critical_in_threads = NULL;
    nC = NULL;
}

hamerly::~hamerly()
{
    delete[] criticals;

    delete[] points_per_class;
    delete[] cdistances;

    if ( critical_in_threads != NULL )
        delete[] critical_in_threads;

    if ( nC != NULL )
        delete[] nC;
}


void hamerly::first_assignation( int nThreads )
{
    for ( int i = 0; i < K; i++ )
    {
        points_per_class[i] = 0;
        for ( int j = 0; j < d; j++ )
            average_per_class[i][j] = 0;
    }

    #pragma omp parallel num_threads( nThreads )
    {
        int id = omp_get_thread_num();
        int p_per_thread = N / nThreads;
        int reminder = N % nThreads;
        int offset = reminder;

        if ( id < reminder )
        {
            p_per_thread += 1;
            offset = 0;
        }

        for ( int i = 0; i < p_per_thread; i++ )
        {
            Point * p = &points[id * p_per_thread + i + offset];
            float mindist = p -> distance( centroids[0] );
            float secmindist = p -> distance( centroids[1] );
            int cent = 0;

            if (mindist > secmindist)
            {
                float temp = secmindist;
                mindist = secmaxdist;
                secmaxdist = temp;
                cent = 1;
            }

            for ( int j = 2; j < K; j++ )
            {
                // Calculate the distance from the point
                float actdist = p -> distance( centroids[j] );

                // Check if it is lower than the lower distance
                if ( actdist < mindist )
                {
                    mindist = actdist;
                    cent = j;
                }
                else if ( actdist < secmindist )
                    secmindist = mindist;
            }

            // assign the point to the closest centroid
            p -> set_centroid( &centroids[ cent ], cent );
            p -> ub = mindist;
            p -> lb = secmindist;

            // update average of all points in that class
            points_per_class[cent] += 1;
            
            for ( int j = 0; j < d; j++ )
                average_per_class[cent][j] += p -> get_coord( j );
        }
    }
}

void hamerly::refresh_assignation( int nThreads )
{
    #pragma omp parallel num_threads( nThreads )
    {
        int id = omp_get_thread_num();
        int p_per_thread = C / nThreads;
        int reminder = C % nThreads;
        int offset = reminder;

        if ( id < reminder )
        {
            p_per_thread += 1;
            offset = 0;
        }

        for ( int i = 0; i < p_per_thread; i++ )
        {
            Point * p = criticals[id * p_per_thread + i + offset];
            float mindist = p -> distance( centroids[0] );
            float secmindist = p -> distance( centroids[1] );
            int cent = 0;

            if (mindist > secmindist)
            {
                float temp = secmindist;
                mindist = secmaxdist;
                secmaxdist = temp;
                cent = 1;
            }


            for ( int j = 2; j < K; j++ )
            {
                // Calculate the distance from the point
                float actdist = p -> distance( centroids[j] );

                // Check if it is lower than the lower distance
                if ( actdist < mindist )
                {
                    secmindist = mindist;
                    mindist = actdist;
                    cent = j;
                }
                else if ( secmindist == -1 )
                    secmindist = actdist;
            }

            // update average of all points in that class
            if ( cent != p -> centroid_index )
            {
                // the point is changing its class so the number of point per each class must be updated
                points_per_class[ p -> centroid_index ] -= 1; 
                points_per_class[ cent ] += 1;

                // since the point is in another class it no longer affect the average of the previous class
                for ( int j = 0; j < d; j++ )
                {
                    float coord = p -> get_coord(j);
                    average_per_class[ p -> centroid_index ][j] -= coord;
                    average_per_class[ cent ][j] += coord;
                }
            }

            // assign the point to the closest centroid
            p -> set_centroid( &centroids[ cent ], cent );
            p -> ub = mindist;
            p -> lb = secmindist;
        }
    }
}

void hamerly::assign_points( int nThreads )
{
    if ( C == 0 )
        first_assignation( nThreads );
    else
        refresh_assignation( nThreads );
}

void hamerly::update_centroids( int nThreads )
{
    #pragma omp parallel num_threads( nThreads )
    {
        int id = omp_get_thread_num();
        int c_per_thread = K / nThreads;
        int reminder = K % nThreads;
        int offset = reminder;

        // account for the remainder
        if ( id < reminder )
        {
            c_per_thread += 1;
            offset = 0;
        }

        for ( int i = 0; i < c_per_thread; i++ )
        {
            int cent_index = id * c_per_thread + i + offset;
            Point * c = &centroids[ cent_index ];
            int sumofsquare = 0;
            
            for ( int j = 0; j < d; j++ )
            {
                float new_coord = average_per_class[cent_index][j] / points_per_class[cent_index];
                float old_coord = c -> get_coord( j );
                sumofsquare += pow( new_coord - old_coord, 2 );

                c -> set_coord( j, new_coord ); 
            }

            cdistances[cent_index] = sumofsquare / points_per_class[cent_index];
        }
    }

    // find the maximum distance and the second highest one
    maxdist = cdistances[0];
    secmaxdist = cdistances[1];
    maxdist_index = 0;

    if (maxdist < secmaxdist)
    {
        float temp = maxdist;
        maxdist = secmaxdist;
        secmaxdist = temp;
        maxdist_index = 1;
    }

    for ( int i = 2; i < K; i++ )
    {
        float act_dist = cdistances[i];

        if ( act_dist > maxdist )
        {
            maxdist = act_dist;
            maxdist_index = i;
        }
        else if ( act_dist > secmaxdist )
            secmaxdist = act_dist;
    }

    update_criticals( nThreads );
}

void hamerly::update_criticals( int nThread )
{
    if ( critical_in_threads != NULL )
    {
        delete[] critical_in_threads;
        delete[] nC;
    }

    critical_in_threads = new list<Point*>[nThread];
    nC = new int[nThread];

    for ( int i = 0; i < nThread; i++ )
        nC[i] = 0;

    #pragma omp parallel num_threads( nThread )
    {
        int id = omp_get_thread_num();

        // calulculate points per each thread
        int p_per_thread = N / nThread;
        int remainder = N % nThread;
        int offset = remainder;

        // account for the remainder
        if ( id < remainder )
        {
            p_per_thread += 1;
            offset = 0;
        }

        for ( int i = 0; i < p_per_thread; i++ )
        {
            Point * p = &points[id * p_per_thread + i + offset];

            // update ub
            p -> ub = p -> ub + ( p -> get_centroid() ) -> dist;

            // update lb
            if ( p -> centroid_index == maxdist_index )
                p -> lb -= secmaxdist;
            else
                p -> lb -= maxdist;

            // if critical
            if ( p -> ub > p -> lb )
            {
                // add critical to the list
                critical_in_threads[id].push_back( p );

                // increase the number of 
                nC[id] += 1;
            }
        }
    }

    // calculate the total number of criticals found by each thread
    C = 0;
    for ( int i = 0; i < nThread; i++ )
        C += nC[i];

    // initiate the array criticals 
    if ( criticals != NULL )
        delete[] criticals;

    criticals = new Point*[C];

    #pragma omp parallel num_threads( nThread )
    {
        int id = omp_get_thread_num();

        // calculate the index to start from
        int index = 0;
        for ( int i = 0; i < id; i++ )
            index += nC[i];

        // insert the value inside the list
        for ( int i = 0; i < nC[id]; i++ )
        {
            criticals[index + i] = critical_in_threads[id].back();
            critical_in_threads[id].pop_back();
        }
    }
}

void hamerly::print_criticals()
{
    for ( int i = 0; i < C; i++ )
        cout << *criticals[i] << endl;
}

int hamerly::get_criticalsNumber()
{
    return C;
}
