#include <iostream>
#include <chrono>
#include "../classes/dataset/dataset.h"
#include "../classes/hamerly/hamerly.h"
#include "../classes/point/point.h"
#include "../classes/settings.h"

int main()
{
    // load the dataset
    dataset data( D, N, K );
    data.rnd_dataset( 0.0, 50.0 );

    hamerly mod( D, N, K, T );

    // insert the points inside the model
    for ( int i = 0; i < N; i++ )
    {
        Point temp( D, data.get_point( i ) );
        mod.insert_point( i, &temp );
    }

    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

    // insert the centroids inside the model
    for ( int i = 0; i < K; i++ )
    {
        Point temp( D, data.get_centroid( i ) );
        mod.insert_centroid( i, &temp );  
    }

    for ( int i = 0; i < ITE; i++ )
    {
        mod.assign_points();
        mod.update_centroids();
    }

    chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = chrono::duration_cast<chrono::duration<double>>(end - start);
    double time = elapsed.count();

    cout << time << endl;
}

// g++ -o openMP.exe openMP.cpp ../classes/dataset/dataset.cpp ../classes/hamerly/hamerly.cpp ../classes/lloyd/lloyd.cpp ../classes/point/point.cpp -fopenmp

// ./openMP.exe