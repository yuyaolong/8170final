//
//  Strut.cpp
//  assignment2
//
//  Created by yuyaolong on 15/10/22.
//  Copyright © 2015年 yuyaolong. All rights reserved.
//

#include "Strut.h"

Strut::Strut(double k, double d, double l, int vertexIndices[2], int faceIndices[2]):
                    k(k),d(d),l(l)
{
    for (int i=0; i<2; i++) {
        this->vertexIndices[i] = vertexIndices[i];
        this->faceIndices[i] = faceIndices[i];
    }
    
}