#include "dataset.h"

#include <random>
#include <fstream>
#include <iostream>
using namespace std;

#define EPS 0.0001

double dataset::random_double( double min, double max )
{
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 eng(rd()); // Seed the generator
    std::uniform_real_distribution<> distr( min, max ); // range definition

    return distr(eng);
}

dataset::dataset( int _d, int _N, int _K )
{
    d = _d;
    N = _N;
    K = _K;

    points = new double*[N];
    for ( int i = 0; i < N; i++ )
        points[i] = new double[d];

    centroids = new double*[K];
    for ( int i = 0; i < K; i++ )
        centroids[i] = new double[d];
}

dataset::~dataset()
{
    for ( int i = 0; i < N; i++ )
        delete[] points[i];
    delete[] points;

        for ( int i = 0; i < K; i++ )
        delete[] centroids[i];
    delete[] centroids;
}

void dataset::rnd_dataset( double min, double max )
{
    for ( int i = 0; i < N; i++ )
        for ( int j = 0; j < d; j++ )
            points[i][j] = random_double( min, max );

    for ( int i = 0; i < K; i++ )
        for ( int j = 0; j < d; j++ )
            centroids[i][j] = random_double( min, max );
}

void dataset::export_dataset( string points_file, string centroids_file )
{
    ofstream file;
    file.open( points_file );

    if ( !file.is_open() )
        throw runtime_error( "Couldn't open " + points_file );

    for ( int i = 0; i < N; i++ )
    {
        for ( int j = 0; j < d; j++ )
        {
            file << points[i][j];

            if ( j < d - 1 )
                file << ",";
        }

        if ( i < N - 1 )
            file << endl;
    }

    file.close();
    file.open( centroids_file );
    
    if ( !file.is_open() )
        throw runtime_error( "Couldn't open " + centroids_file );

    for ( int i = 0; i < K; i++ )
    {
        for ( int j = 0; j < d; j++ )
        {
            file << centroids[i][j];

            if ( j < d - 1 )
                file << ",";
        }

        if ( i < K - 1)
            file << endl;
    }

    file.close();
}

bool dataset::double_eq( double f1, double f2 )
{
    bool upper_cond = f1 < f2 + EPS;
    bool lower_cond = f1 > f2 - EPS;

    return upper_cond && lower_cond;
}


bool operator == ( const dataset &d1, const dataset &d2 )  
{
    if ( d1.N != d2.N || d1.K != d2.K || d1.d != d2.d )
        return false;
    
    bool equal = true;

    for ( int i = 0; i < d1.N && equal; i++ )
        for ( int j = 0; j < d1.d && equal; j++ )
            if ( !dataset::double_eq( d1.points[i][j], d2.points[i][j] )  )
                equal = false;

    for ( int i = 0; i < d1.K && equal; i++ )
        for ( int j = 0; j < d1.d && equal; j++ )
            if ( !dataset::double_eq(d1.centroids[i][j], d2.centroids[i][j]) )
                equal = false;

    return equal;
}

void dataset::load_dataset( string points_file, string centroids_file )
{
    ifstream file;
    file.open( points_file );

    if ( !file.is_open() )
        throw runtime_error( "Couldn't open " + points_file );

    string line;

    int line_count = 0;
    while ( getline( file, line ) )
    {
        if ( line_count >= N )
            throw runtime_error( "Found more points than the declared number, in " + points_file );

        string coord_s;
        double coord;
        int last_comma = 0;
        for ( int i = 0; i < d; i++ )
        {
            int new_comma = line.find( ',', last_comma );
            
            coord_s = line.substr( last_comma, new_comma );
            coord = stof( coord_s );

            points[line_count][i] = coord;
            last_comma = new_comma + 1;
        }
        line_count++;
    }

    if ( line_count < N )
        throw runtime_error( "Found less points than the declared number, in " + points_file );

    file.close();
    file.open( centroids_file );

    if ( !file.is_open() )
        throw runtime_error( "Couldn't open " + centroids_file );

    line_count = 0;

    while ( getline( file, line ) )
    {
        if ( line_count >= K )
            throw runtime_error( "Found more points than the declared number, in " + centroids_file );

        string coord_s;
        double coord;
        int last_comma = 0;
        for ( int i = 0; i < d; i++ )
        {
            int new_comma = line.find( ',', last_comma );
            
            coord_s = line.substr( last_comma, new_comma );
            coord = stof( coord_s );

            centroids[line_count][i] = coord;
            last_comma = new_comma + 1;
        }
        line_count++;
    }

    if ( line_count < K )
        throw runtime_error( "Found less points than the declared number, in " + centroids_file );

    file.close();

}

double * dataset::get_point( int i )
{
    return points[i];
}

double * dataset::get_centroid( int i )
{
    return centroids[i];
}
