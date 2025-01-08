#include <iostream>
#include <chrono>
#include <fstream>
#include "classes/dataset/dataset.h"
#include "classes/hamerly/hamerly.h"
#include "classes/lloyd/lloyd.h"
#include "classes/point/point.h"

using namespace std;

#define T 9 // max number of threads
#define ITE 100 // number of iterations per cycle

#define N 1143 // number of points
#define D 13 // data dimensionality
#define K 10 // number of centroids
#define STEP 1  // step size to increase the number of point per test

int main()
{
    // load dataset
    dataset data( D, N, K );
    //data.load_dataset( "dataset/points.txt", "dataset/centroids.txt" );
    data.rnd_dataset( 0.0, 50.0 );

    int tests = D / STEP;
    int ndim = STEP;

    ofstream lloydfile;
    lloydfile.open( "testres/lloyd.txt" );

    if ( !lloydfile.is_open() )
        throw runtime_error( "Couldn't open testres/lloyd.txt" );

    for ( int i = 0; i < tests; i++ )
    {
        cout << "running with " << ndim << " dimension" << endl;
        // execute lloyd
        lloyd modlloyd( ndim, N, K );

        for ( int j = 0; j < N; j++ )
        {
            Point temp( ndim, data.get_point( j ) );
            modlloyd.insert_point( j, &temp );
        }

        for ( int j = 0; j < K; j++ )
        {
            Point temp( ndim, data.get_centroid( j ) );
            modlloyd.insert_centroid( j, &temp );           
        }

        chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
        for ( int j = 0; j < ITE; j++ )
        {
            modlloyd.assign_points();
            modlloyd.update_centroids();
        }

        chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = chrono::duration_cast<chrono::duration<double>>(end - start);
        double lloydtime = elapsed.count();

        lloydfile << "dimnesion: " << ndim << "; time: " << lloydtime << " seconds" << endl;

        // execute the hamerly algorithm with threads from 1 to T
        double prevtime = lloydtime;
        for ( int j = 0; j < T; j++ )
        {
            hamerly modham( ndim, N, K, j + 1 );

            for ( int k = 0; k < N; k++ )
            {
                Point temp( ndim, data.get_point( k ) );
                modham.insert_point( k, &temp );
            }

            for ( int k = 0; k < K; k++ )
            {
                Point temp( ndim, data.get_centroid( k ) );
                modham.insert_centroid( k, &temp );           
            }

            chrono::high_resolution_clock::time_point startham = chrono::high_resolution_clock::now();
            for ( int k = 0; k < ITE; k++ )
            {
                modham.assign_points();
                modham.update_centroids();
            }
            
            chrono::high_resolution_clock::time_point endham = chrono::high_resolution_clock::now();
            auto elapsedham = chrono::duration_cast<chrono::duration<double>>(endham - startham);
            double hamtime = elapsedham.count();

            ofstream hamfile;
            hamfile.open( "testres/hamerly" + to_string( j + 1 ) + "T.txt", ios::app );

            if ( !hamfile.is_open() )
                throw runtime_error( "Couldn't open testres/hamerly" + to_string( j + 1 ) + "T.txt" );

            hamfile << "points: " << ndim << "; time: " << hamtime << " seconds" << endl;
            hamfile.close();

            if( i == tests - 1 )
            {
                double speedup = 1 - (hamtime / prevtime);

                ofstream speedupfile;
                speedupfile.open( "testres/speedup.txt", ios::app );

                if ( !speedupfile.is_open() )
                    throw runtime_error( "Couldn't open testres/speedup.txt" );

                speedupfile << "lloyd vs hamerly " << (j + 1) << ": " << speedup << " %" << endl;
                speedupfile.close();

                //prevtime = hamtime;
            }
        }

        ndim += STEP;
    }

    lloydfile.close();
}

// g++ -o testbench.exe testbench.cpp ./classes/dataset/dataset.cpp ./classes/hamerly/hamerly.cpp ./classes/lloyd/lloyd.cpp ./classes/point/point.cpp -fopenmp