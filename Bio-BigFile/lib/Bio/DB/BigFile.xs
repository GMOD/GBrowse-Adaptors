#include <unistd.h>
#include <math.h>

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "udc.h"
#include "localmem.h"
#include "bigWig.h"

/* Let Perl redefine these */
#undef TRUE
#undef FALSE
#undef warn

#ifdef PERL_CAPI
#define WIN32IO_IS_STDIO
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef FCGI
 #include <fcgi_stdio.h>
#else
 #ifdef USE_SFIO
  #include <config.h>
 #else
  #include <stdio.h>
 #endif
 #include <perlio.h>
#endif

typedef struct bbiFile     *Bio__DB__BBIFile;
typedef struct bbiInterval *Bio__DB__bbiInterval;
typedef struct bbiIntervalList {
    struct lm          *lm;
    struct bbiInterval *head;
  } *Bio__DB__bbiIntervalList;
typedef struct bbiSummaryList {
  int                      size;
  struct bbiSummaryElement *summary;
} *Bio__DB__BigWigExtendedSummary;


MODULE = Bio::DB::BigFile PACKAGE = Bio::DB::BigFile PREFIX=bw_

Bio::DB::BBIFile
bw_bigWigFileOpen(packname="Bio::DB::BigFile",filename)
   char* packname
   char* filename
  PROTOTYPE: $$
  CODE:
  RETVAL = bigWigFileOpen(filename);
  OUTPUT:
  RETVAL

MODULE = Bio::DB::BigFile PACKAGE = Bio::DB::BBIFile PREFIX=bbi_
int
bbi_bigWigIntervalDump(bwf,chrom,start,end,maxCount=0,out=stdout)
    Bio::DB::BBIFile bwf
    char            *chrom
    unsigned int     start
    unsigned int     end
    int              maxCount
    FILE            *out
    CODE:
      RETVAL = bigWigIntervalDump(bwf,chrom,start,end,maxCount,out);
    OUTPUT:
      RETVAL

void
bbi_close(bbi)
    Bio::DB::BBIFile bbi
    CODE:
      bigWigFileClose(&bbi);

void
bbi_DESTROY(bbi)
    Bio::DB::BBIFile bbi
    CODE:
      bigWigFileClose(&bbi);

Bio::DB::bbiIntervalList
bbi_bigWigIntervalQuery(bwf,chrom,start,end)
    Bio::DB::BBIFile bwf
    char            *chrom
    unsigned int     start
    unsigned int     end
    PREINIT:
    struct bbiIntervalList *list;
    CODE:
    list = Newxz(list,1,struct bbiIntervalList);
    list->lm = lmInit(0);
    list->head = bigWigIntervalQuery(bwf,chrom,start,end,list->lm);
    RETVAL = list;
    OUTPUT:
      RETVAL

SV*
bbi_bigWigSummaryArray(bwf,chrom,start,end,summaryType=0,size)
   Bio::DB::BBIFile bwf
   char            *chrom
   unsigned int     start
   unsigned int     end
   unsigned int     summaryType
   unsigned int     size
  PREINIT:
    int     i;
    boolean result;
    double  *values;
    AV      *avref;
  CODE:
    values = Newxz(values,size,double);
    result = bigWigSummaryArray(bwf,chrom,start,end,summaryType,size,values);
    if (result != TRUE) {
      Safefree(values);
      XSRETURN_EMPTY;
    } else {
      avref = (AV*) sv_2mortal((SV*)newAV());
      for (i=0;i<size;i++)
	av_push(avref, newSVnv(values[i]));
      Safefree(values);
      RETVAL = (SV*) newRV((SV*)avref);
    }
  OUTPUT:
     RETVAL

Bio::DB::BigWigExtendedSummary
bbi_bigWigSummaryArrayExtended(bwf,chrom,start,end,size)
   Bio::DB::BBIFile bwf
   char            *chrom
   unsigned int     start
   unsigned int     end
   unsigned int     size
  PREINIT:
    int     i;
    boolean result;
    struct bbiSummaryElement          *summary;
    Bio__DB__BigWigExtendedSummary     summaryList;
    SV     *p;
  CODE:
   summary = Newxz(summary,size,struct bbiSummaryElement);
   result  = bigWigSummaryArrayExtended(bwf,chrom,start,end,size,summary);
   if (result != TRUE) {
     Safefree(summary);
     XSRETURN_EMPTY;
   } else {
      summaryList = Newxz(summaryList,1,struct bbiSummaryList);
      summaryList->size = size;
      summaryList->summary = summary;
      p = newSV(sizeof(summaryList));
      sv_setref_pv(p,"Bio::DB::BigWigExtendedSummary",(void*) summaryList);
      RETVAL = summaryList;
   }
  OUTPUT:  
    RETVAL


MODULE = Bio::DB::BigFile PACKAGE = Bio::DB::BigWigExtendedSummary PREFIX=bwes_

int
bwes_size(el);
      Bio::DB::BigWigExtendedSummary el
      CODE:
        RETVAL = el->size;
      OUTPUT:
        RETVAL

unsigned long
bwes_validCount(el,i)
      Bio::DB::BigWigExtendedSummary el
      int i
      CODE:
        if (i>el->size-1)
	  croak("Attempt to read past end of ExtendedSummary results %d > %d",i,el->size-1);
        RETVAL = el->summary[i].validCount;
      OUTPUT:
        RETVAL
    
double
bwes_minVal(el,i)
      Bio::DB::BigWigExtendedSummary el
      int i
      CODE:
        if (i>el->size-1)
	  croak("Attempt to read past end of ExtendedSummary results %d > %d",i,el->size-1);
        RETVAL = el->summary[i].minVal;
      OUTPUT:
        RETVAL
    
double
bwes_maxVal(el,i)
      Bio::DB::BigWigExtendedSummary el
      int i
      CODE:
        if (i>el->size-1)
	  croak("Attempt to read past end of ExtendedSummary results %d > %d",i,el->size-1);
        RETVAL = el->summary[i].maxVal;
      OUTPUT:
        RETVAL

double
bwes_sumData(el,i)
      Bio::DB::BigWigExtendedSummary el
      int i
      CODE:
        if (i>el->size-1)
	  croak("Attempt to read past end of ExtendedSummary results %d > %d",i,el->size-1);
        RETVAL = el->summary[i].sumData;
      OUTPUT:
        RETVAL
    
double
bwes_sumSquares(el,i)
      Bio::DB::BigWigExtendedSummary el
      int i
      CODE:
        if (i>el->size-1)
	  croak("Attempt to read past end of ExtendedSummary results %d > %d",i,el->size-1);
        RETVAL = el->summary[i].sumSquares;
      OUTPUT:
        RETVAL

void
bwes_DESTROY(el)
      Bio::DB::BigWigExtendedSummary el
      CODE:
      if (el->summary != NULL) {
	Safefree(el->summary);
      }
      el->summary = NULL;

MODULE = Bio::DB::BigFile PACKAGE = Bio::DB::bbiIntervalList PREFIX=bbil_

Bio::DB::bbiInterval
bbil_head(list)
  Bio::DB::bbiIntervalList list
  CODE:
     RETVAL = list->head;
  OUTPUT:
     RETVAL

void
bbil_DESTROY(list)
  Bio::DB::bbiIntervalList list
  CODE:
	if (list->lm != NULL) lmCleanup(&list->lm);

MODULE = Bio::DB::BigFile PACKAGE = Bio::DB::bbiInterval PREFIX=bbii_

Bio::DB::bbiInterval
bbii_next(interval)
      Bio::DB::bbiInterval interval
  CODE:
     RETVAL = interval->next;
  OUTPUT:
     RETVAL

unsigned int
bbii_start(interval)
      Bio::DB::bbiInterval interval
  CODE:
     RETVAL = interval->start;
  OUTPUT:
     RETVAL

unsigned int
bbii_end(interval)
      Bio::DB::bbiInterval interval
  CODE:
     RETVAL = interval->end;
  OUTPUT:
     RETVAL

double
bbii_value(interval)
      Bio::DB::bbiInterval interval
  CODE:
     RETVAL = interval->val;
  OUTPUT:
     RETVAL

