//
//  Point.hpp
//  assignment2
//
//  Created by yuyaolong on 15/10/22.
//  Copyright © 2015年 yuyaolong. All rights reserved.
//

#ifndef __POINT_H__
#define __POINT_H__

#include "Vector.h"

class Point {
public:
    double mass;
    Vector3d xposition;
    Vector3d velocity;
    Vector3d force;
    
    Vector2d textureCoordinate;
    
    Point(double m, Vector3d x, Vector3d v, Vector3d f, Vector2d t);
    
private:
    Point();
    
    
};

#endif 
