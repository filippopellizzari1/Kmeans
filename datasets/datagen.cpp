#include <iostream>
#include "../classes/dataset/dataset.h"

using namespace std;

#define N 100000 // Number of points
#define K 3 // Number of centroids
#define D 3 // Data dimensionality

int main()
{
    dataset data(D, N, K);
    data.rnd_dataset( 0.0f, 50.0f );
    data.export_dataset( "points3D.txt", "centroids3D.txt" );
}

// ../classes/dataset/dataset.cpp