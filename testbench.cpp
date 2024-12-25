#include <iostream>
#include "classes/lloyd/lloyd.h"
#include "classes/hamerly/hamerly.h"
#include "classes/point/point.h"

using namespace std;

#define N 100 // Number of points
#define K 3 // Number of centroids
#define D 2 // Data dimensionality
#define T 3 // number of threads
#define ITE 10 // number of iterations

int main()
{   
    lloyd mod1( D, N, K );
    mod1.generate_points( 0, 50 );
    mod1.generate_centroids( 0, 50 );
  
    hamerly mod2( D, N, K, 3 );
    for ( int i = 0; i < N; i++ )
        mod2.insert_point( i, mod1.get_point( i ) );

    for ( int i = 0; i < K; i++ )
        mod2.insert_centroid( i, mod1.get_centroid( i ) );   

    for ( int i = 0; i < ITE; i++ )
    {
        mod1.assign_points();
        mod2.assign_points();

        mod1.update_centroids();
        mod2.update_centroids();

        if ( mod1 == mod2 )
            cout << "1 == 2" << endl;
        else
            cout << "1 != 2" << endl;

        mod1.export_data( "./testbench_inout/out1.txt" );
        mod2.export_data( "./testbench_inout/out2.txt" );
    }
    
    cout << "done" << endl << endl;
}

// lloyd mod1;
// mod1.load_data( "./testbench_inout/input.txt" );

// hamerly mod2;
// mod2.load_data( "./testbench_inout/input.txt" );

// for ( int i = 0; i < ITE; i++ )
// {
//     mod1.assign_points();
//     mod1.update_centroids();

//     mod2.assign_points(1);
//     mod2.update_centroids(1);
//     mod2.print_criticals();

//     if ( mod1 == mod2 )
//         cout << "1 == 2" << endl;
//     else
//         cout << "1 != 2" << endl;
// }