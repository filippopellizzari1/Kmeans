#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>
using namespace std;

class Point
{
    private:
        int d;
        double * coords;

    public:
        double dist, ub, lb;
        int centroid_index;
        Point();
        Point( int _d, double * _coords );
        Point( const Point &p );
        Point( int _d, double * _coords, double _dist, double _ub, double _lb, int _centroid_index );
        ~Point();


        double get_coord( int i ) const;
        void set_coord( int i, double val );
        int get_dim() const;
        void set_centroid( Point * _centroid, int _centroid_index ); 
        Point * get_centroid() const; 
        double distance( const Point &p );

        friend bool operator == ( const Point &p1, const Point &p2 );
        friend bool operator != ( const Point &p1, const Point &p2 );
        friend Point operator + ( const Point &p1, const Point &p2 );
        Point& operator =(const Point &p );
};

ostream& operator <<( ostream &os, const Point &p ); 

#endif
