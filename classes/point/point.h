#ifndef __POINT_H__
#define __POINT_H__

#include <iostream>
#include <vector>
using namespace std;

class Point
{
    private:
        int d;
        float dist;
        vector<float> coords;
        Point * centroid;

    public:
        Point( int _d, float * _coords );
        Point( const Point &p );

        float get_coord( int i ) const;
        int get_dim() const;
        void set_centroid( Point * _centroid ); 
        Point * get_centroid(); 
        void set_disance( float _dist );
        float get_distance();
        float distance( const Point &p1, const Point &p2 );

        friend bool operator == ( const Point &p1, const Point &p2 );
        friend Point operator + ( const Point &p1, const Point &p2 );
};

ostream& operator <<( ostream &os, const Point &p ); 

#endif
