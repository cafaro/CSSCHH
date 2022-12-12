//
//  datagen.h
//  datagen
//
//  Created by marco on 15/01/14.
//  Copyright (c) 2014 marco. All rights reserved.
//

#ifndef datagen_datagen_h
#define datagen_datagen_h

struct DataParams {
    long         n;
    long         k;
    long         seedPrng;
    uint32_t     seedZipf;
    uint32_t     seed;
    double       skew;
    double       hza;
};


#endif
