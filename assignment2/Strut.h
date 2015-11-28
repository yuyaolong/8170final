//
//  Strut.hpp
//  assignment2
//
//  Created by yuyaolong on 15/10/22.
//  Copyright © 2015年 yuyaolong. All rights reserved.
//

#ifndef __STRUT_H__
#define __STRUT_H__

class Strut {
public:
    double k;
    double d;
    double l;
    int vertexIndices[2];
    int faceIndices[2];
    
    Strut(double k, double d, double l, int vertexIndices[2], int faceIndices[2]);
private:
    Strut();
};

#endif 