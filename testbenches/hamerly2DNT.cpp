#include <iostream>
#include "../classes/hamerly/hamerly.h"
#include "../classes/point/point.h"
#include "../classes/dataset/dataset.h"
#include <chrono>
#include <fstream>

using namespace std;

#define N 50000 // Number of points
#define K 3 // Number of centroids
#define D 2 // Data dimensionality
#define ITE 100 // number of iterations
#define STEPS 500 // increase in the points 

int main( int argc, char* argv[] )
{
    if ( argc != 2 )
        throw runtime_error( "invalid arguments" );

    int threads = stoi( argv[1] );

    string centroids_file = "../datasets/centroids2D.txt";
    string points_file = "../datasets/points2D.txt";
    string output_file = "../testbench_out/hamerly2D" + to_string(threads) + "T.txt";

    dataset data( D, N, K );
    data.load_dataset( points_file, centroids_file ); 

    int iterations = N / STEPS;

    ofstream file;
    file.open( output_file );

    if ( !file.is_open() )
        throw runtime_error( "couldn't open file " + output_file );

    int nPoints = STEPS;
    cout << "ite: " << iterations << endl;
    for ( int i = 0; i < iterations; i++ )
    {
        hamerly model( D, nPoints, K, threads );

        for ( int j = 0; j < nPoints; j++ )
        {
            Point temp( D, data.get_point(j) );
            model.insert_point( j, &temp );
        }

        for ( int j = 0; j < K; j++ )
        {
            Point temp( D, data.get_centroid(j) );
            model.insert_centroid( j, &temp );
        } 

        cout << "running test with : " << nPoints << " points" << endl;
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        for ( int j = 0; j < ITE; j++ )
        {
            model.assign_points();
            model.update_centroids();
        }

        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);

        file << "points: " << nPoints <<  "; time: " << time_span.count() << " seconds" << endl;

        nPoints += STEPS;
    }

    file.close();
}