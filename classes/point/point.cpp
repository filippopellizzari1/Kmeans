#include <iostream>
#include <cmath>
#include "point.h"

Point::Point( int _d, float * _coords )
{
    dist = 0;
    ub = 0;
    lb = 0;
    centroid = NULL;
    d = _d;

    for ( int i = 0; i < _d; i++ )
        coords.push_back( _coords[i] );
}

// Copy constructor
Point::Point( const Point &p )
{
    d = p.get_dim();
    centroid = p.get_centroid();
    ub = p.ub;
    lb = p.lb;
    dist = p.dist;

    for ( int i = 0; i < d; i++ )
        coords.push_back( p.get_coord( i ) );
}

float Point::get_coord( int i ) const
{
    return coords[i];
}

float Point::distance( const Point &p1, const Point &p2 )
{
    if ( p1.get_dim() != p2.get_dim() )
        throw invalid_argument( "Cannot calculate distance between points with different dimensionality" ); 

    float sumofsquare = 0;

    for ( int i = 0; i < p1.get_dim(); i++ )
        sumofsquare += pow( p1.get_dim() + p2.get_dim(), 2 );

    return sqrt( sumofsquare );
}

int Point::get_dim() const
{
    return d;
}

void Point::set_centroid( Point * _centroid )
{
    centroid = _centroid;
}

Point * Point::get_centroid() const
{
    return centroid;
}

ostream& operator <<( ostream &os, const Point &p )
{
    os << "Point( ";

    for ( int i = 0; i < p.get_dim(); i++ )
    {
        os << p.get_coord( i );

        if ( i != p.get_dim() - 1)
            os << ", ";
        else
            os << " ";
    }

    os << ")";

    return os;
}

bool operator == ( const Point &p1, const Point &p2 )
{
    if ( p1.get_dim() != p2.get_dim() )
        return false;

    bool equal = true;

    for ( int i = 0; i < p1.get_dim() && equal; i++ )
        if ( p1.get_coord( i ) != p2.get_coord( i ) )
            equal = false;

    return equal;
}

Point operator + ( const Point &p1, const Point &p2 )
{
    if ( p1.get_dim() != p2.get_dim() )
        throw invalid_argument( "Cannot sum points with different dimensionality" );

    float * coords = new float[p1.get_dim()];
    for ( int i = 0; i < p1.get_dim(); i++ )
        coords[i] = p1.get_coord( i ) + p2.get_coord( i );
    
    Point temp( p1.get_dim(), coords );
    delete[] coords;

    return temp;
}
