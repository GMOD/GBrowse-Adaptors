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

#include <unistd.h>
#include <math.h>
#include "bam.h"
#include "khash.h"
#include "faidx.h"

/* stolen from bam_aux.c */
#define MAX_REGION 1<<29

typedef faidx_t*      Bio__DB__Sam__Fai;
typedef bamFile       Bio__DB__Bam;
typedef bam_header_t* Bio__DB__Bam__Header;
typedef bam1_t*       Bio__DB__Bam__Alignment;
typedef bam_index_t*  Bio__DB__Bam__Index;
typedef struct {
  SV* callback;
  SV* data;
} fetch_callback_data;
typedef fetch_callback_data *fetch_callback_dataptr;
typedef struct {
  int    start;
  int    end;
  double width;
  int*   bin;
} coverage_graph;
typedef coverage_graph *coverage_graph_ptr;

void XS_pack_charPtrPtr( SV * arg, char ** array, int count) {
  int i;
  AV * avref;
  avref = (AV*)sv_2mortal((SV*)newAV());
  for (i=0; i<count; i++) {
    av_push(avref, newSVpv(array[i], strlen(array[i])));
  }
  SvSetSV( arg, newRV((SV*)avref));
}

int bam_fetch_fun (const bam1_t *b, void *data) {
  dSP;
  int count;

  fetch_callback_dataptr fcp;
  SV* callback;
  SV* callbackdata;
  SV* alignment_obj;
  int returnval;
  bam1_t *b2;

  fcp          = (fetch_callback_dataptr) data;
  callback     = fcp->callback;
  callbackdata = fcp->data;

  /* turn the bam1_t into an appropriate object */
  b2 = bam_dup1(b);
  // b2->hash = b->hash;

  alignment_obj = sv_setref_pv(newSV(sizeof(bam1_t)),"Bio::DB::Bam::Alignment",(void*) b2);

  /* set up subroutine stack for the call */
  ENTER;
  SAVETMPS;
  PUSHMARK(SP);
  XPUSHs(sv_2mortal(alignment_obj));
  XPUSHs(callbackdata);
  PUTBACK;

  /* execute the call */
  count = call_sv(callback,G_ARRAY);

  if (count == 0) {
    returnval = 0;
  } else if (count == 1) {
    returnval = POPi;
  } else
    croak("bam_fetch() callback must return empty or a scalar integer");

  SPAGAIN;

  PUTBACK;
  FREETMPS;
  LEAVE;

  return returnval;
}

int add_pileup_line (const bam1_t *b, void *data) {
  bam_plbuf_t *pileup = (bam_plbuf_t*) data;
  bam_plbuf_push(b,pileup);
  return 0;
}

int coverage_from_pileup_fun (uint32_t tid, 
			      uint32_t pos, 
			      int n, 
			      const bam_pileup1_t *pl, 
			      void *data) {
  coverage_graph_ptr  cgp;
  int                 bin;

  cgp = (coverage_graph_ptr) data;
  if (pos >= cgp->start && pos <= cgp->end) {
    bin = (pos-cgp->start)/cgp->width;
    cgp->bin[bin] += n;
  }

  return 0;
}

/* copied from bam_aux.c because "we need it" */
uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2])
{
       uint8_t *s;
       int y = tag[0]<<8 | tag[1];
       s = bam1_aux(b);
       while (s < b->data + b->data_len) {
               int type, x = (int)s[0]<<8 | s[1];
               s += 2;
               if (x == y) return s;
               type = toupper(*s); ++s;
               if (type == 'C') ++s;
               else if (type == 'S') s += 2;
               else if (type == 'I' || type == 'F') s += 4;
               else if (type == 'D') s += 8;
               else if (type == 'Z' || type == 'H') { while (*s) putchar(*s++); ++s; }
       }
       return 0;
}

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Sam::Fai PREFIX=fai_

Bio::DB::Sam::Fai
fai_open(packname="Bio::DB::Sam::Fai", filename)
  char * packname
  char * filename
 PROTOTYPE: $$
 CODE:
    RETVAL = fai_load(filename);
 OUTPUT:
    RETVAL

void
fai_destroy(fai)
  Bio::DB::Sam::Fai fai
  PROTOTYPE: $
  CODE:
    fai_destroy(fai);

SV*
fai_fetch(fai,reg)
  Bio::DB::Sam::Fai    fai
    const char *reg
  PROTOTYPE: $$$
  PREINIT:
    char     *seq;
    int       len;
  CODE:
    seq = fai_fetch(fai,reg,&len);
    if (seq == NULL)
       XSRETURN_EMPTY;
    RETVAL = newSVpv(seq,len);
    free((void*)seq);
  OUTPUT:
    RETVAL


MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam PREFIX=bam_

Bio::DB::Bam
bam_open(packname="Bio::DB::Bam", filename)
      char * packname
      char * filename
      PROTOTYPE: $$
      CODE:
        RETVAL = bam_open(filename,"r");
      OUTPUT:
      RETVAL

void
bam_DESTROY(bam)
   Bio::DB::Bam bam
PROTOTYPE: $
CODE:
   bam_close(bam);

int
bam_index_build(packname="Bio::DB::Bam", filename)
   char *      packname
   const char * filename
  CODE:
     RETVAL = bam_index_build(filename);
  OUTPUT:
     RETVAL

Bio::DB::Bam::Index
bam_index_open(packname="Bio::DB::Bam", filename)
      char * packname
      char * filename
    PROTOTYPE: $$
    CODE:
    RETVAL = bam_index_load(filename);
    OUTPUT:
    RETVAL

Bio::DB::Bam::Header
bam_header(bam)
    Bio::DB::Bam bam
    PROTOTYPE: $
    PREINIT:
      bam_header_t *bh;
    CODE:
      bgzf_seek(bam,0,0);
      bh = bam_header_read(bam);
      RETVAL = bh;
    OUTPUT:
      RETVAL

char*
bam_tell(bam)
    Bio::DB::Bam bam
PROTOTYPE: $
CODE:
    int64_t t = bam_tell(bam);
    char    string[128];
    sprintf(string,"%llu",t);
    RETVAL = string;
OUTPUT:
    RETVAL

void
bam_seek(bam,pos,dir)
    Bio::DB::Bam bam
    int pos
    int dir
PROTOTYPE: $$$
CODE:
    bam_seek(bam,pos,dir);

Bio::DB::Bam::Alignment
bam_read1(bam)
    Bio::DB::Bam bam
  PROTOTYPE: $
  PREINIT:
    bam1_t *b;
  CODE:
    b = bam_init1();
    if (bam_read1(bam,b) >= 0) {
      RETVAL = b;
    }
    else
       XSRETURN_EMPTY;
  OUTPUT:
    RETVAL      

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam::Alignment PREFIX=bam_

void
bam_DESTROY(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    bam_destroy1(b);

int
bam_tid(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.tid;
OUTPUT:
    RETVAL

int
bam_pos(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.pos;
OUTPUT:
    RETVAL

int
bam_calend(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
   RETVAL=bam_calend(&b->core,bam1_cigar(b));
OUTPUT:
   RETVAL    

int
bam_cigar2qlen(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
   RETVAL=bam_cigar2qlen(&b->core,bam1_cigar(b));
OUTPUT:
   RETVAL    

int
bam_qual(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.qual;
OUTPUT:
    RETVAL

int
bam_flag(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.flag;
OUTPUT:
    RETVAL

int
bam_n_cigar(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.n_cigar;
OUTPUT:
    RETVAL

int
bam_l_qseq(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.l_qseq;
OUTPUT:
    RETVAL

SV*
bam_qseq(b)
Bio::DB::Bam::Alignment b
PROTOTYPE: $
PREINIT:
    char* seq;
    int   i;
CODE:
    seq = Newxz(seq,b->core.l_qseq+1,char);
    for (i=0;i<b->core.l_qseq;i++) {
      seq[i]=bam_nt16_rev_table[bam1_seqi(bam1_seq(b),i)];
    }
    RETVAL = newSVpv(seq,b->core.l_qseq);
    Safefree(seq);
OUTPUT:
    RETVAL

int
bam_mtid(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.mtid;
OUTPUT:
    RETVAL

int
bam_mpos(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.mpos;
OUTPUT:
    RETVAL

int
bam_isize(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->core.isize;
OUTPUT:
    RETVAL

int
bam_l_aux(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->l_aux;
OUTPUT:
    RETVAL

SV*
bam_aux_get(b,tag)
   Bio::DB::Bam::Alignment b
   char*               tag
PROTOTYPE: $$
PREINIT:
   int           type;
   uint8_t       *s;
CODE:
   s    = bam_aux_get_core(b,tag);
   if (s==0)
      XSRETURN_EMPTY;
   type = *s++;
   if (type == 'c')
     RETVAL = newSViv((int32_t)*(int8_t*)s);
   else if (type=='C')
     RETVAL = newSViv((int32_t)*(uint8_t*)s);
   else if (type=='s' || type=='S')
     RETVAL = newSViv((int32_t)*(int16_t*)s);
   else if (type=='i' || type=='I')
     RETVAL = newSViv(*(int32_t*)s);
   else if (type=='f')
     RETVAL = newSVnv(*(float*)s);
   else if (type=='Z' || type=='H')
     RETVAL = newSVpv((char*)s,0);
   else
     XSRETURN_EMPTY;
OUTPUT:
   RETVAL

void
bam_aux_keys(b)
Bio::DB::Bam::Alignment b
PROTOTYPE: $
PREINIT:
   uint8_t *s;
   uint8_t type, key[2];
PPCODE:
   {
     s = bam1_aux(b);  /* s is a khash macro */
     while (s < b->data + b->data_len) {
       XPUSHs(sv_2mortal(newSVpv(s,2)));
       s   += 2; 
       type = *s++;
       if (type == 'A') { printf("A:%c", *s); ++s; }
       else if (type == 'C') { ++s; }
       else if (type == 'c') { ++s; }
       else if (type == 'S') { s += 2; }
       else if (type == 's') { s += 2; }
       else if (type == 'I') { s += 4; }
       else if (type == 'i') { s += 4; }
       else if (type == 'f') { s += 4; }
       else if (type == 'Z' || type == 'H') { while (*s) ++s; }
     }
   }

SV*
bam_data(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=newSVpv(b->data,b->data_len);
OUTPUT:
    RETVAL

int
bam_data_len(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->data_len;
OUTPUT:
    RETVAL

int
bam_m_data(b)
    Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=b->m_data;
OUTPUT:
    RETVAL

SV*
bam_qname(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=newSVpv(bam1_qname(b),0);
OUTPUT:
    RETVAL

int
bam_paired(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FPAIRED) != 0;
OUTPUT:
  RETVAL

int
bam_proper_pair(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FPROPER_PAIR) != 0;
OUTPUT:
  RETVAL

int
bam_unmapped(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FUNMAP) != 0;
OUTPUT:
  RETVAL

int
bam_munmapped(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
    RETVAL=(b->core.flag&BAM_FMUNMAP) != 0;
OUTPUT:
  RETVAL

int
bam_reversed(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
  RETVAL=bam1_strand(b);
OUTPUT:
  RETVAL

int
bam_mreversed(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
CODE:
  RETVAL=bam1_mstrand(b);
OUTPUT:
  RETVAL

SV*
bam_cigar(b)
  Bio::DB::Bam::Alignment b
PROTOTYPE: $
PREINIT:
    int        i;
    uint32_t  *c;
    AV        *avref;
CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    c     = bam1_cigar(b);
    for (i=0;i<b->core.n_cigar;i++)
      av_push(avref, newSViv(c[i]));
    RETVAL = (SV*) newRV((SV*)avref); 
OUTPUT:
  RETVAL

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam::Header PREFIX=bam_

int
bam_n_targets(bamh)
  Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  CODE:
    RETVAL = bamh->n_targets;
  OUTPUT:
    RETVAL

char**
bam_target_name(bamh)
  Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  PREINIT:
    int count_charPtrPtr;
  CODE:
    count_charPtrPtr=bamh->n_targets;
    RETVAL = bamh->target_name;
  OUTPUT:
    RETVAL

SV*
bam_target_len(bamh)
    Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  PREINIT:
    int i;
    AV * avref;
  CODE:
    avref = (AV*) sv_2mortal((SV*)newAV());
    for (i=0;i<bamh->n_targets;i++)
       av_push(avref, newSViv(bamh->target_len[i]));
    RETVAL = (SV*) newRV((SV*)avref); 
  OUTPUT:
    RETVAL

SV*
bam_text(bamh)
  Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  CODE:
    /* in case text is not null terminated, we copy it */
    RETVAL = newSVpv(bamh->text,bamh->l_text);
  OUTPUT:
    RETVAL

void
bam_parse_region(bamh,region)
    Bio::DB::Bam::Header bamh
    char*            region
    PROTOTYPE: $
    PREINIT:
       int seqid,start,end;
    PPCODE:
    {
      bam_parse_region(bamh,
		       region,
		       &seqid,
		       &start,
		       &end);
      if (seqid < 0)
	XSRETURN_EMPTY;
      else {
	EXTEND(sp,3);
	PUSHs(sv_2mortal(newSViv(seqid)));
	PUSHs(sv_2mortal(newSViv(start)));
	PUSHs(sv_2mortal(newSViv(end)));
      }
    }

void
bam_DESTROY(bamh)
  Bio::DB::Bam::Header bamh
  PROTOTYPE: $
  CODE:
    bam_header_destroy(bamh);

MODULE = Bio::DB::Sam PACKAGE = Bio::DB::Bam::Index PREFIX=bam_

int
bam_fetch(bai,bfp,ref,start,end,callback,callbackdata)
  Bio::DB::Bam::Index bai
  Bio::DB::Bam        bfp
  int   ref
  int   start
  int   end
  CV*   callback
  SV*   callbackdata
PREINIT:
  fetch_callback_data fcd;
CODE:
  {
    fcd.callback = (SV*) callback;
    fcd.data     = callbackdata;
    RETVAL = bam_fetch(bfp,bai,ref,start,end,&fcd,bam_fetch_fun);
  }
OUTPUT:
    RETVAL

SV*
bam_coverage(bai,bfp,ref,start,end,bins)
    Bio::DB::Bam::Index bai
    Bio::DB::Bam        bfp
    int             ref
    int             start
    int             end
    int             bins
PREINIT:
    coverage_graph  cg;
    bam_plbuf_t    *pileup;
    AV*             array;
    SV*             cov;
    int             i;
    bam_header_t   *bh;
CODE:
  {
      if (end >= MAX_REGION) {
          bgzf_seek(bfp,0,0);
          bh  = bam_header_read(bfp);
          end = bh->target_len[ref];
          bam_header_destroy(bh);
      }
      if (bins > (end-start))
         bins = end-start;

      /* coverage graph used to communicate to our callback
	  the region we are sampling */
      cg.start = start;
      cg.end   = end;
      cg.width = ((double)(end-start))/bins;
      cg.bin   = calloc(bins+1,sizeof(int));
      // Newxz(cg.bin,bins,int); /* gives compile warnings? */

      /* accumulate coverage into the coverage graph */
      pileup   = bam_plbuf_init(coverage_from_pileup_fun,(void*)&cg);
      bam_fetch(bfp,bai,ref,start,end,(void*)pileup,add_pileup_line);
      bam_plbuf_push(NULL,pileup);

      // is crashing?
      bam_plbuf_destroy(pileup);

      /* now normalize to coverage/bp and convert into an array */
      array = newAV();
      av_extend(array,bins);
      for  (i=0;i<bins;i++)
         av_store(array,i,newSVnv(((float)cg.bin[i])/cg.width));
      Safefree(cg.bin);
      //Safefree(cg.bin);

      RETVAL = (SV*) newRV((SV*)array);
  }
OUTPUT:
    RETVAL
