#include <iostream>
#include <chrono>
#include "../classes/dataset/dataset.h"
#include "../classes/point/point.h"
#include "../classes/lloyd/lloyd.h"
#include "../classes/MPIHamerly/MPIHamerly.h"
#include "../classes/settings.h"

using namespace std;

int main( int argc, char ** argv )
{
    // load dataset
    dataset data( D, N, K );
    data.rnd_dataset( 0.0, 50.0 );

    MPIHamerly mod( &argc, &argv );

    // insert points into the model
    for ( int i = 0; i < N; i++ )
    {
        mod.init_point( i, data.get_point(i) );
    }

    // insert centroids into the model
    for ( int i = 0; i < K; i++ )
    {
        mod.init_centroid( i, data.get_centroid(i) );
    }

    mod.start_time();
    for ( int i = 0; i < ITE; i++ )
    {
        mod.assign_points();
        mod.update_centroids();
    }
    mod.elapsed();

}

// g++ -o hammpi.exe hammpi.cpp ../classes/dataset/dataset.cpp ../classes/lloyd/lloyd.cpp ../classes/MPIHamerly/MPIHamerly.cpp ../classes/point/point.cpp -I "Path/to/MPI/Include" -L "Path/to/MPI/Lib" -lmsmpi

// mpiexec.exe -n 3 .\hammpi.exe
