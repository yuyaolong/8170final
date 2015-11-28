//
//  Face.cpp
//  assignment2
//
//  Created by yuyaolong on 15/10/22.
//  Copyright © 2015年 yuyaolong. All rights reserved.
//

#include "Face.h"
Face::Face(int s[3], double v[3], int p[3]){
    
    for (int i=0; i<3; i++) {
        this->strutIndices[i] = s[i];
        this->vertexAngles[i] = v[i];
        this->pointIndices[i] = p[i];
    }
    
    
}