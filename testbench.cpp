#include <iostream>
#include "classes/lloyd/lloyd.h"
#include "classes/hamerly/hamerly.h"
#include "classes/point/point.h"
#include "classes/dataset/dataset.h"
#include <chrono>

using namespace std;

#define N 50000 // Number of points
#define K 3 // Number of centroids
#define D 3 // Data dimensionality
#define T 3 // number of threads
#define ITE 100 // number of iterations

int main()
{   
    dataset data(D, N, K);
    data.rnd_dataset( 0.0f, 50.0f );
    data.export_dataset( "dataset/points.txt", "dataset/centroids.txt" );
    data.load_dataset( "dataset/points.txt", "dataset/centroids.txt" );

    lloyd mod1( D, N, K );
    hamerly mod2( D, N, K, 3 );


    for ( int i = 0; i < N; i++ )
    {
        Point temp( D, data.get_point(i) );
        mod1.insert_point( i, &temp );
        mod2.insert_point( i, &temp );
    }

    for ( int i = 0; i < K; i++ )
    {
        Point temp( D, data.get_centroid(i) );
        mod1.insert_centroid( i, &temp );
        mod2.insert_centroid( i, &temp );
    } 

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    for ( int i = 0; i < ITE; i++ )
    {
        // mod1.assign_points();
        // mod1.update_centroids();

        mod2.assign_points();
        mod2.update_centroids();
        cout << "i:" << i << endl;

        // if ( mod1 == mod2 )
        //     cout << "1 == 2" << endl;
        // else
        //     cout << "1 != 2" << endl;
    }

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    
    cout << "time: " << time_span.count() << " seconds" << endl;
}