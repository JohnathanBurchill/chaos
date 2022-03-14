#ifndef CDF_ATTRS_H
#define CDF_ATTRS_H

#include <cdf.h>


CDFstatus addgEntry(CDFid id, long attrNum, long entryNum, const char *entry);


typedef struct varAttr {
    char * name;
    char * type;
    char * units;
    char * desc;
    double validMin;
    double validMax;
    char * format;
} varAttr;

CDFstatus addVariableAttributes(CDFid id, varAttr attr);

void addAttributes(CDFid id, const char *calVersion, const char satellite, const char *version, double minTime, double maxTime);


#endif // CDF_ATTRS_H
