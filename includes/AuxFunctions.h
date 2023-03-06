#include <cstdio>
#include "TTree.h"
#include "TGraph.h"
#include <vector>

#define PBARSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBARW 60

inline void PrintProgress(float percentage) {
    int val = (int)(percentage*100);
    int lpad =(int)(percentage*PBARW);
    int rpad = PBARW-lpad;
    printf("\r%3d%% [%.*s%*s]",val,lpad,PBARSTR,rpad,"");
    fflush(stdout);
}
