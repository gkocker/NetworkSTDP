//fast exponential function, from G.C. Cawley, Neural Comput. (2000),
//based on the original algorithm from N.N. Schraudolph, Neural Comput. (1999).
// using namespace std;

// #include <iostream>

static union 
{
	  double d;
	    struct {
#ifdef LITTLE_ENDIAN
	    int j,i;
#else 
		    int i,j;
#endif
	  } n;
} _eco;

#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)
