#ifndef _CDate
#define _CDate

#include "standard.h"

class CDate
{
public:
    
    int Year;
    int Month;
    int Day;
    
    double DecimalYear;

	CDate();

    CDate(double sdate);

	~CDate();
};

#endif