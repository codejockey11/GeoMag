#include "CDate.h"

CDate::CDate()
{
	memset(this, 0x00, sizeof(CDate));
}

CDate::CDate(double sdate)
{
	memset(this, 0x00, sizeof(CDate));

	DecimalYear = sdate;
}

CDate::~CDate()
{

}