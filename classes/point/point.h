#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>
using namespace std;

class Point
{
    private:
        int d;
        Point * centroid;
        float * coords;

    public:
        float dist, ub, lb;
        int centroid_index;
        Point();
        Point( int _d, float * _coords );
        Point( const Point &p );
        ~Point();


        float get_coord( int i ) const;
        void set_coord( int i, float val );
        int get_dim() const;
        void set_centroid( Point * _centroid, int _centroid_index ); 
        Point * get_centroid() const; 
        float distance( const Point &p );

        friend bool operator == ( const Point &p1, const Point &p2 );
        friend Point operator + ( const Point &p1, const Point &p2 );
        Point& operator =(const Point &p );
};

ostream& operator <<( ostream &os, const Point &p ); 

#endif
