#include <iostream>
#include <cmath>
#include "point.h"

#define EPS 0.0001

Point::Point()
{
    d = 0;
    dist = 0;
    ub = 0;
    lb = 0;
    coords = NULL;
    centroid_index = -1;
}

Point::Point( int _d, double * _coords )
{
    dist = 0;
    ub = 0;
    lb = 0;
    d = _d;
    centroid_index = -1;

    coords = new double[d];
    for ( int i = 0; i < _d; i++ )
        coords[i] = _coords[i];
}

// Copy constructor
Point::Point( const Point &p )
{
    d = p.get_dim();
    ub = p.ub;
    lb = p.lb;
    dist = p.dist;
    centroid_index = p.centroid_index;

    coords = new double[d];

    for ( int i = 0; i < d; i++ )
        coords[i] =  p.get_coord( i );
}

Point::Point( int _d, double * _coords, double _dist, double _ub, double _lb, int _centroid_index)
{
    d = _d;
    dist = _dist;
    ub = _ub;
    lb = _lb;
    centroid_index = _centroid_index;
    
    coords = new double[d];
    for ( int i = 0; i < d; i++ )
        coords[i] = _coords[i];
}

Point::~Point()
{
    delete[] coords;
}

double Point::get_coord( int i ) const
{
    return coords[i];
}

void Point::set_coord( int i, double val )
{
    coords[i] = val;
}

double Point::distance( const Point &p )
{
    if ( d != p.get_dim() )
        throw invalid_argument( "Cannot calculate distance between points with different dimensionality" ); 

    double sumofsquare = 0;

    for ( int i = 0; i < d; i++ )
        sumofsquare += pow( coords[i] - p.get_coord( i ), 2 );

    return sqrt(sumofsquare);
}

int Point::get_dim() const
{
    return d;
}

void Point::set_centroid( Point * _centroid, int _centroid_index )
{
    centroid_index = _centroid_index;
}

ostream& operator <<( ostream &os, const Point &p )
{
    os << "Coords: ";

    for ( int i = 0; i < p.get_dim(); i++ )
    {
        os << p.get_coord( i );

        if ( i != p.get_dim() - 1)
            os << ", ";
        else
            os << "; ";
    }

    os << "C: " <<  p.centroid_index << "; D: " << p.get_dim() << "; LB: " << p.lb << "; UB: " << p.ub << "; DIST: " << p.dist;

    return os;
}

bool double_eq( double f1, double f2 )
{
    // double f1_copy = f1;
    // double f2_copy = f2;
    // while ( f1_copy > 0 || f2_copy > 0 )
    // {
    //     int digit1 = f1_copy;
    //     int digit2 = f2_copy;

    //     cout << digit1 << "," << digit2 << endl;

    //     f1_copy = f1_copy - digit1;
    //     f2_copy = f2_copy - digit2;
        
    //     f1_copy *= 10;
    //     f2_copy *= 10;
    // }

    // cout << endl;

    bool upper_cond = f1 < f2 + EPS;
    bool lower_cond = f1 > f2 - EPS;

    return upper_cond && lower_cond;
}

bool operator == ( const Point &p1, const Point &p2 )
{
    if ( p1.get_dim() != p2.get_dim() )
        return false;

    if ( p1.centroid_index != p2.centroid_index )
        return false;

    if ( p1.dist != p2.dist )
        return false;

    bool equal = true;

    for ( int i = 0; i < p1.get_dim() && equal; i++ )
        if ( !double_eq( p1.get_coord(i), p2.get_coord(i) ) )
            equal = false;

    return equal;
}

bool operator != ( const Point &p1, const Point &p2 )
{
    return !( p1 == p2 );
}

Point operator + ( const Point &p1, const Point &p2 )
{
    if ( p1.get_dim() != p2.get_dim() )
        throw invalid_argument( "Cannot sum points with different dimensionality" );

    double * coords = new double[p1.get_dim()];
    for ( int i = 0; i < p1.get_dim(); i++ )
        coords[i] = p1.get_coord( i ) + p2.get_coord( i );
    
    Point temp( p1.get_dim(), coords );
    delete[] coords;

    return temp;
}

Point& Point::operator= ( const Point &p )
{
    d = p.get_dim();
    ub = p.ub;
    lb = p.lb;
    dist = p.dist;
    centroid_index = p.centroid_index;

    if (coords != NULL)
        delete[] coords;

    coords = new double[d];

    for ( int i = 0; i < d; i++ )
        coords[i] =  p.get_coord( i );  

    return *this;
}