/*! \file IRIS_decoder.c
\brief Program to decode IRIS RAW data and metadata of one IRIS subtask to
intermediate data format (see ODIM_struct.h)

This code is modified from IRIS sigmet/src/utils/examples/change_raw.C by 
Harri Hohti of FMI and can be used only with IRIS libraries and headers.<BR> 
The original copyright information is in the source code.<BR>

<B>The program accepts four options:</B><BR>
<B>-?</B> : usage <BR>
<B>-v</B> : verbose output <BR>
<B>-d</B> : prints all metadata information <BR>
<B>-s</B> : only prints available quantities (IRIS and ODIM names). <B>Generates no output file</B>. <BR>

 After option(s) the next argument is the full input file path (IRIS RAW file) and the last
 argument is the full output file path.<BR>
 <B>Example:</B> ./IRIS_decoder -d $IRIS_PRODUCT_RAW/VAN101231235505.RAW1234 output_file.dat <BR>

 In operational use the program is typically installed in IRIS output pipe script, 
 which is connected to IRIS output device where RAW data is sent from IRIS output menu.

 The documentation of IRIS data types and structures are freely downloadable from
 Vaisala/Sigmet ftp site: ftp://ftp.sigmet.com/outgoing/manuals/program/3data.pdf 
 */


/*   ORIGINAL COPYRIGHT INFORMATION OF THE change_raw.C
 -------------------------------------------------------------------------
 *   COPYRIGHT (c) 2000, 2001, 2002, 2003, 2004, 2005, 2008, 2009 
 *             VAISALA INC., WESTFORD MASSACHUSETTS, U.S.A.
 * 
 * THIS SOFTWARE IS FURNISHED UNDER A LICENSE AND MAY BE USED AND  COPIED
 * ONLY  IN  ACCORDANCE WITH  THE  TERMS  OF  SUCH  LICENSE  AND WITH THE
 * INCLUSION OF THE ABOVE COPYRIGHT NOTICE.  THIS SOFTWARE  OR  ANY OTHER
 * COPIES  THEREOF MAY NOT BE PROVIDED OR OTHERWISE MADE AVAILABLE TO ANY
 * OTHER PERSON.  NO TITLE TO AND OWNERSHIP OF  THE  SOFTWARE  IS  HEREBY
 * TRANSFERED.<BR>
 _________________________________________________________________________*/


#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include "sig_data_types.h"
#include "sigtypes.h"
#include "dsp.h"
#include "headers.h"
#include "iris_task.h"
#include "ingest.h"
#include "product.h"
#include "setup.h"
#include "user_lib.h"
#include "dsp_lib.h"
#include "ODIM_struct.h"

#define SIGMET_SETUP_H 1
#define PRODPTR( IREC, IOFF ) \
  ((UINT1 *)(prod_c + (IREC * TAPE_RECORD_LEN) + IOFF))

struct full_time_spec FullYmds_Time;
struct ymds_time* const pYmds_Time = &FullYmds_Time.Ymds;

/*!\var MetaData *meta 
\brief pointer to MetaData structure */ 
MetaData *meta; 
int argF; /**<\brief Argument index pointing to input file */
SINT4  irec_c ;         /**<\brief Record number (history: 6144 byte tape records) */
SINT4  ioff_c ;         /**<\brief Offset within record */
UINT1 *prod_c ;         /**<\brief Pointer to RAW product */
UINT2  pol_code;
int scans; /**<\brief Number of scans in the input subtask RAW file */
int quantities; /**<\brief Number of quantities saved in the input subtask RAW file */
int SCAN_QUANTITIES=FALSE; /**<\brief Set TRUE if option -s given: 
Only prints available quantities, no output file generated */
int VERB=FALSE; /**<\brief Set TRUE if option -v given: Increases verbosity.  */
int DUMPALL=FALSE; /**<\brief Set TRUE if option -d given: Dumps all information and decodes to output file */
size_t fres;
double NyqV;
double NyqW;
double UnambV;
double antgain=-1,antgainH,antgainV;
double radomeloss=0,radomelossH,radomelossV;

FILE *DATAF; /**<\brief File pointer for temporary output data */
FILE *METAF; /**<\brief Output metadata and data are written to this file (in that order) */
UINT4 totsize; /**<\brief Total size of MetaData structure */
UINT1 POL_H,POL_V,POL_HV; /**<\brief POL_ variables are booleans indicating polarization used */
UINT1 IS_XHDR = 0;
UINT1 singlePRF;
int64_t binmethod_avg;
double RXlossH, RXlossV, TXlossH, TXlossV,radconstHV,nomTXpower;

/*!\fn void get_raw_bytes( SINT2 *buf_a, SINT4 icnt_a )
\brief Routine to extract bytes from RAW record, and static variables used by the routine
 */
void get_raw_bytes( SINT2 *buf_a, SINT4 icnt_a );
/** \brief Processes the RAW product */
static void product_raw(struct raw_product *pPRaw);
void usage( void );
/** \brief Sets name and ODIM code (see ODIM_struct.h) of quantities in MetaData structure */
void ProcessDatatype(SINT4 quantity, SINT4 *datatypes);
/** \brief Gives date (YYYYMMDD) and time (hhmmss) strings from ymds_time structure. <BR>
<B> All times in time related functions are in UTC! </B> */
time_t give_date_time(char *date, char *time, struct ymds_time ymds);
/** \brief Gives UNIX seconds from date and time strings. */
time_t sec_from_date_time(char *date, char *time);
/** \brief Prints data \#iQ (quantity) attributes from dataset (scan) \#iS */ 
void DumpDataAttributes(int iS, int iQ);
/** \brief Prints dataset (scan) \#iS attributes */ 
void DumpDatasetAttributes(int iS);
/** \brief Prints common attributes */ 
void DumpCommonAttributes(void);
/** \brief Prints all attributes (previous dumps combined) */ 
void DumpAllAttributes(void);

/* ================================================== */
/** Exit status will be "1" for any kind of error, "0" for successful return.
*/
int main( int argc, char *argv[] )
{
  MESSAGE istatus ; SINT4 iSize, iChan ; 
  struct raw_product *pRaw;

  setbuf(stdout,NULL);
  argF=1;
  {
    int i;

    if(argc==1) { usage(); exit(1); }
    for(i=1;i<argc;i++) if(argv[i][0]=='-')
    { 
      if(argv[i][1]=='h') { usage(); exit(1); }
      if(argv[i][1]=='s') SCAN_QUANTITIES = TRUE;
      if(argv[i][1]=='v') VERB = TRUE;
      if(argv[i][1]=='d') { DUMPALL = TRUE; VERB = TRUE; }
      argF++;
    }
  }

  istatus = imapopen( argv[argF], FALSE, (void**)(void*)&pRaw, &iSize, &iChan ) ;
  if( istatus != SS_NORMAL ) 
  {
    fprintf( stderr,  "Could not open '%s' for Read/Write.\n", argv[argF] ) ;
    return(1) ;
  }

  if(!SCAN_QUANTITIES) METAF=fopen(argv[argF+1],"w");

  meta=calloc(1,sizeof(MetaData));
  totsize=sizeof(MetaData);
  product_raw(pRaw);
  free(meta);
  istatus = imapclose( pRaw, iSize, iChan ) ;
  if( istatus != SS_NORMAL ) { return(2); }

  exit( EXIT_SUCCESS ) ;
}


/* ============================================================== */
/**
 * Extract information about a RAW product.  Entered with a pointer
 * to the beginning (first 6144-byte record) of the full product.
 */
static void product_raw(struct raw_product *pRaw)
{

  struct ingest_header *inghdr;
  struct product_hdr *prodhdr;
  SINT4 type_i,iAz, datatypes[64], 
        datatype, scan, scanlo, scanhi ;
  /* struct data_convert Convert; */
  char cdate[10]={0}, ctime[10]={0};
  UINT1 *ray_times;
  time_t csecs;

  DATAF=tmpfile();

  /* The first two "records" of the product consist of a product
   * header structure and an ingest header structure.  Each is padded
   * out to one record.
   */

  prod_c = (UINT1 *) pRaw ;
  inghdr = &(pRaw->Record[1].IHeader);

  ray_times=calloc(inghdr->icf.irtotl,1);

  /* Make sure that this file really looks like a raw product.  Check
   * the ID's in the product and ingest headers.
   */
  prodhdr=&pRaw->Record[0].PHeader;    

  /*  if((pRaw->Record[0].PHeader.hdr.id    != ST_PRODUCT_HDR) || */
  if((prodhdr->hdr.id    != ST_PRODUCT_HDR) ||
     (inghdr->hdr.id     != ST_INGEST_HDR ) ||
     (inghdr->tcf.hdr.id != ST_TASK_CONF  ) ) 
     {
       fprintf( stderr,  "ERROR: File headers contain invalid ID's\n" ) ; exit(1) ;
     }

  /* Count up the number of parameters that were recorded.  This is
   * needed in order to skip through the scans properly.
   */
  for( quantities=type_i=0 ; type_i < 128 ; type_i++ )
  {
      if( lDspMaskTest( &inghdr->tcf.dsp.DataMask, type_i) )
      { 
         if(type_i == DB_XHDR) 
         { 
            IS_XHDR = 1; 
            if(SCAN_QUANTITIES || VERB) printf("\nExtended header (XHDR) found, skipping.\n"); 
            continue; 
         }
         datatypes[quantities] = type_i; 
         quantities++ ; 
     }
  }


  NyqV = fNyquistVelocity( prodhdr->end.iprf, prodhdr->end.itrig,
                           prodhdr->end.ilambda, prodhdr->end.ipolar);
  NyqW = fNyquistWidth( prodhdr->end.iprf, prodhdr->end.ilambda,
                         prodhdr->end.ipolar );

  /*  if(VERB) printf("Nyquist velocity %.2f, Nyquist width %.2f\n",NyqV,NyqW); */

  /* Get scan limits for the data in this particular RAW product.  If
   * it is a single scan product, then use just the one that is
   * selected in the PSI structure.  Otherwise convert the entire set.
   */

  /*  if( pRaw->Record[0].PHeader.pcf.psi.raw.iflags & RAW_FLG_SWEEP ) */ 

  if( prodhdr->pcf.psi.raw.iflags & RAW_FLG_SWEEP ) 
  {
    /*    scanlo  = pRaw->Record[0].PHeader.pcf.psi.raw.isweep ; */
    scanlo  = prodhdr->pcf.psi.raw.isweep ;
    scanhi  = scanlo;
    scans=1;
  } else 
  {
    scanlo = 1 ;
    scanhi = inghdr->tcf.scan.isweeps ;
    scans=scanhi;
  }
  meta->how.scan_count=scans; /* V23 */
  ProcessDatatype(-1,datatypes);

  if(SCAN_QUANTITIES || VERB) 
  {
    /* Listing of available quantities */
    printf("\nQuantities measured\n");
    printf("IRIS     ODIM\n");
    printf("--------------------\n");
    for(type_i=0;type_i<quantities;type_i++)
    { 
      printf("%-8s %-8s\n",sdata_name6(datatypes[type_i]),meta->dataset[0].data[type_i].what.quantity);
    }
    printf("====================\n\n");
    if(!VERB) return;
  }


  /* Initializes root group attributes for ODIM HDF5 conversion */
  {
    long volinter,secs;
    double melt;
    char *envp=NULL, envstr[255], test_env[20];

    meta->scans=scans;
    /* /what attributes */
    volinter=atoi(getenv("ODIM_VOLUME_INTERVAL"))*60;
    secs=inghdr->icf.VolumeYmds.isec;
    /* volume nominal time is rounded down to nearest interval minute */
    inghdr->icf.VolumeYmds.isec=(secs/volinter)*volinter; 
    give_date_time(cdate,ctime,inghdr->icf.VolumeYmds);
    sprintf(meta->what.date,"%s",cdate);
    sprintf(meta->what.time,"%s",ctime);

    /* /where attributes */
    sprintf(meta->where.sitecode,"%.3s",inghdr->icf.sSitename);

    sprintf(test_env,"ODIM_%s_source",meta->where.sitecode);
    if(getenv(test_env)==NULL)
    {
      printf("\nThe IRIS RAW file comes from previously unknown radar having site name defined as \"%s\".\nSo there is no mandatory environment variable %s defined for it.\n",inghdr->icf.sSitename,test_env); 
      printf("Please add the ODIM_%s_* and IRIS_%s_* environment variables to your conversion environment.\nSee test.sh provided with the software.\n\n",meta->where.sitecode,meta->where.sitecode);
      
      exit(111);
    }

    meta->where.lon=fDegFromBin4(inghdr->icf.ilon);
    meta->where.lat=fDegFromBin4(inghdr->icf.ilat);
    meta->where.height=(double)inghdr->icf.ialtitude/100.0; /* cm -> m */
    meta->where.towerheight=(double)inghdr->icf.irad_hgt;

    /* /how attributes */
    sprintf(meta->how.sw_version,"%s",inghdr->icf.siris_version);

    /* get the antenna gains and losses for site (not in IRIS files) */
    {
       char *sitecode=meta->where.sitecode;
       char *agp=NULL, *agHp=NULL, *agVp=NULL; 
       char *rlop=NULL, *rloHp=NULL, *rloVp=NULL; 

           sprintf(envstr,"IRIS_%s_antgain",sitecode);
           agp = getenv(envstr);
           if(agp) antgain=atof(agp);

           sprintf(envstr,"IRIS_%s_antgainH",sitecode);
           agHp = getenv(envstr);
           if(agHp) antgainH=atof(agHp);

           sprintf(envstr,"IRIS_%s_antgainV",sitecode);
           agVp = getenv(envstr);
           if(agVp) antgainV=atof(agVp);

           if(agp) antgainH = antgainV = antgain;
           if(agHp && !agVp) antgainV = antgainH;
           if(!agHp && agVp) antgainH = antgainV;

           sprintf(envstr,"ODIM_%s_radomeloss",sitecode);
           rlop = getenv(envstr);
           if(rlop) radomeloss=atof(rlop);

           sprintf(envstr,"ODIM_%s_radomelossH",sitecode);
           rloHp = getenv(envstr);
           if(rloHp) radomelossH=atof(rloHp);

           sprintf(envstr,"ODIM_%s_radomelossV",sitecode);
           rloVp = getenv(envstr);
           if(rloVp) radomelossV=atof(rloVp);

           if(rlop) radomelossH = radomelossV = radomeloss;
           if(rloHp && !rloVp) radomelossV = radomelossH;
           if(!rloHp && rloVp) radomelossH = radomelossV;
    }
    meta->how.antgain  = antgain;
    meta->how.antgainH = antgainH;
    meta->how.antgainV = antgainV;

    meta->how.radomeloss  = radomeloss;
    meta->how.radomelossH = radomelossH;
    meta->how.radomelossV = radomelossV;

    sprintf(envstr,"ODIM_%s_poltype",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) sprintf(meta->how.poltype,"%s",envp);
    else
    { 
       envp=getenv("ODIM_poltype");
       if(envp) sprintf(meta->how.poltype,"%s",envp);
       else meta->how.poltype[0]=0;
    }

    sprintf(envstr,"ODIM_%s_TXtype",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) sprintf(meta->how.TXtype,"%s",envp);
    else
    { 
       envp=getenv("ODIM_TXtype");
       if(envp) sprintf(meta->how.TXtype,"%s",envp);
       else meta->how.TXtype[0]=0;
    }

    sprintf(envstr,"ODIM_%s_system",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) sprintf(meta->how.system,"%s",envp);
    else
    { 
       envp=getenv("ODIM_system");
       if(envp) sprintf(meta->how.system,"%s",envp);
       else meta->how.TXtype[0]=0;
    }

    sprintf(envstr,"ODIM_%s_OUR",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) meta->how.OUR = atof(envp);
    else meta->how.OUR=-1;

    meta->how.comment[0]=0;    
    sprintf(envstr,"ODIM_%s_how_comment",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) sprintf(meta->how.comment,"%s",envp);
    envp=getenv("ODIM_how_comment"); 
    if(envp) sprintf(meta->how.comment,"%s",envp);
    

    /* V23 IRIS specific 1.4 dB, correction to continuous calib signal vs pulsed */
    envp=getenv("FiniteBandwithLoss");
    if(envp) meta->how.CWloss = atof(envp); else meta->how.CWloss = 1.4; 
    
    sprintf(envstr,"IRIS_%s_RXlossH",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) RXlossH = atof(envp); else RXlossH = 0.0;
    meta->how.RXlossH = RXlossH - meta->how.CWloss;

    sprintf(envstr,"IRIS_%s_RXlossV",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) RXlossV = atof(envp); else RXlossV = 0.0; 
    meta->how.RXlossV = RXlossV - meta->how.CWloss;

    sprintf(envstr,"IRIS_%s_TXlossH",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) meta->how.TXlossH = atof(envp); else meta->how.TXlossH = 0.0; 

    sprintf(envstr,"IRIS_%s_TXlossV",meta->where.sitecode);
    envp=getenv(envstr);
    if(envp) meta->how.TXlossV = atof(envp); else meta->how.TXlossV = 0.0; 

    meta->how.beamwH=fDegFromBin4(inghdr->tcf.misc.iHorzBeamWidth);
    meta->how.beamwV=fDegFromBin4(inghdr->tcf.misc.iVertBeamWidth);

    meta->how.wavelength=(double)inghdr->tcf.misc.ilambda/100.0; /* 1/100 cm -> cm */
 
    if (IC_MELTING_UNKNOWN == inghdr->icf.iMeltingHeight) melt=DBL_MAX;
    else melt = 0.001*(double)(SINT2)(0x8000 ^ inghdr->icf.iMeltingHeight ); /* m -> km */
    meta->how.freeze=melt;

    nomTXpower=(double)inghdr->tcf.misc.ixmt_pwr/1000.0; /* W -> kW */
    meta->how.peakpwr = nomTXpower;
    meta->how.RAC=(double)inghdr->tcf.dsp.igas_atten/100000.0;
    meta->how.gasattn = meta->how.RAC;
  /* meta->how.dynrange ? */
  }

  /* Initialize record pointer to the first data record (third one in)
   * and loop for each scan that was recorded.  We will scan through
   * the entire product (to check for errors), and will perhaps print
   * some subset of it.
   */

  if(VERB && !DUMPALL) DumpCommonAttributes();

  irec_c = 2 ;                  /* Record number */
  ioff_c = 0 ;                  /* Offset within record */

  for( scan = scanlo ; scan <= scanhi ; scan++ ) 
  {
    struct ingest_data_header inghdrs[64] ; 
    SINT2 iQ,iS,tQ,azgates;
    char sTimeBuf[TIMENAME_SIZE];
    UINT1 *scandata[64];
    int databytes[64],first_ray=0,rotsgn=1;
    /* double first_az, last_az; */
    double azdiff;
    int CHANGE_QUANTITY_RESOLUTION;
    long N,scansize[64],min_raysecs,max_raysecs,raysecs;
    struct tm Sdd;

    min_raysecs=100000;
    max_raysecs=0;
    POL_H=0;
    POL_V=0;
    POL_HV=0;

    if(scans==1) iS=0; else iS=scan-1;
    meta->dataset[iS].quantities=quantities;
    meta->dataset[iS].how.scan_index=scan; /* V23 */

    /* Extract the INGEST data file headers for each of the parameters
     * that were recorded.  The headers appear sequentially in the
     * first record of each scan.
     */
    for( iQ=tQ=0 ; iQ < quantities+IS_XHDR ; iQ++,tQ++ ) 
    {
      get_raw_bytes( (SINT2*)&inghdrs[tQ], INGEST_DATA_HEADER_SIZE ) ;
      if(iQ==0 && IS_XHDR) tQ--;
    }

    if(VERB)
    {
      printf("================================================\n"); 
      printf( "Scan %2.2d began at: %s\n", scan, shhmmssddmonyyyy_r( &inghdrs[0].time, sTimeBuf ));
    }
    
    {
      UINT4 iLowPRF;
      
      switch( prodhdr->end.ipolar )
      {
         case POL_VERT_FIX:
           sprintf(meta->dataset[iS].how.polarization,"V");
           sprintf(meta->dataset[iS].how.polmode,"single-V");
           POL_V=1; 
           meta->dataset[iS].how.POL_V=1;
         break;
         case POL_ALTERNATING:
           sprintf(meta->dataset[iS].how.polarization,"H|V");
           sprintf(meta->dataset[iS].how.polmode,"switched-dual");
           POL_HV=1; 
           meta->dataset[iS].how.POL_HV=1;
         break;
         case POL_SIMULTANEOUS:
           sprintf(meta->dataset[iS].how.polarization,"H,V");
           sprintf(meta->dataset[iS].how.polmode,"simultaneous-dual");
           POL_HV=1;  
           meta->dataset[iS].how.POL_HV=1;
         break;
         default:
           sprintf(meta->dataset[iS].how.polarization,"H");
           sprintf(meta->dataset[iS].how.polmode,"single-H");
           POL_H=1; 
           meta->dataset[iS].how.POL_H=1;
         break;
      }

      /* /dataset/where attributes */
      meta->dataset[iS].where.bin_elangle = inghdrs[0].iangle;
      meta->dataset[iS].where.elangle = fDegFromBin2(inghdrs[0].iangle);
      meta->dataset[iS].where.rstart = 1.0e-5*(double)inghdr->tcf.rng.ibin_first; /* cm -> km */
      meta->dataset[iS].where.rscale=0.01*(double)inghdr->tcf.rng.ibin_out_step;  /* cm -> m */
      meta->dataset[iS].where.nrays=inghdr->icf.irtotl;
      /*      printf("ISNDX %d\n",inghdrs[0].isndx); */

      /* /dataset/how attributes */
      {
        short i;
        char c;

        for(i=0;i<12;i++) if((c=inghdr->tcf.end.stname[i])!=' ')  meta->dataset[iS].how.task[i]=c;
      }

      /* printf("flags %d %d\n",(unsigned short)inghdr->tcf.cal.iflags,(unsigned short)inghdr->tcf.cal.iflags2); */

      meta->dataset[iS].how.RXbandwidth = (double)inghdr->tcf.cal.iReceiverBandwidth/1000.0;
      meta->dataset[iS].how.nomTXpower = nomTXpower; 

      if(POL_H | POL_HV) 
      {
         meta->dataset[iS].how.MDSH = (double)inghdr->tcf.cal.iI0Horiz/100.0; /* 1/100 dBm -> dBm */
         meta->dataset[iS].how.radconstH=(double)inghdr->tcf.cal.iRadarConstantHoriz/100.0; /* 1/100 dB -> dB */
      }

      if(POL_V | POL_HV) 
      {
         meta->dataset[iS].how.MDSV = (double)inghdr->tcf.cal.iI0Vert/100.0; /* 1/100 dBm -> dBm */
         radconstHV=(double)inghdr->tcf.cal.iRadarConstantVert/100.0; /* 1/100 dB -> dB */
         meta->dataset[iS].how.radconstHV = radconstHV;
         meta->dataset[iS].how.radconstH = radconstHV + RXlossH;
      }

      meta->dataset[iS].how.pulsewidth = 0.01*(double)inghdr->tcf.dsp.ipw; /* 1/100 us -> us */
      meta->dataset[iS].how.NEZH=(double)inghdr->GParm.iz_calib/16.0;

      if(POL_HV) 
          meta->dataset[iS].how.HVratio=(double)inghdr->GParm.inse_hv_ratio/100.0;

      iLowPRF = NINT( fPrfLowFromHighCase( prodhdr->end.iprf, prodhdr->end.itrig ));

      switch( prodhdr->end.itrig )
      {
          default:
             meta->dataset[iS].how.lowprf=prodhdr->end.iprf;
             meta->dataset[iS].how.highprf=prodhdr->end.iprf;
             singlePRF=TRUE;
          break;
          case PRF_2_3: case PRF_3_4: case PRF_4_5:
             meta->dataset[iS].how.highprf=prodhdr->end.iprf;
             meta->dataset[iS].how.lowprf=iLowPRF;
             singlePRF=FALSE;
          break;
      }
 
      meta->dataset[iS].how.radhoriz = 1.0e-5*(double)inghdr->tcf.rng.ibin_last; /* cm -> km */
      meta->dataset[iS].how.UnambVel = UnambV = 0.0025 * meta->how.wavelength * prodhdr->end.iprf;
      meta->dataset[iS].how.NI = NyqV;
      meta->dataset[iS].how.NyqWidth = NyqW;


      meta->dataset[iS].how.SQI=(double)inghdr->tcf.cal.isqi_thr/256.0;
      meta->dataset[iS].how.CSR=-(double)inghdr->tcf.cal.iccr_thr/16.0;
      meta->dataset[iS].how.LOG=(double)inghdr->tcf.cal.izns_thr/16.0;
      meta->dataset[iS].how.SNRT=(double)inghdr->tcf.cal.isig_thr/16.0;
      meta->dataset[iS].how.PMI=(double)inghdr->tcf.cal.ipmi_thr/256.0;
      meta->dataset[iS].how.Cflags[0]=inghdr->tcf.cal.iflags;
      meta->dataset[iS].how.Cflags[1]=inghdr->tcf.cal.iflags2;
      meta->dataset[iS].how.LDR_bias=(double)inghdr->tcf.cal.ildr_bias/100.0;
      meta->dataset[iS].how.ZDR_bias=(double)inghdr->tcf.cal.izdr_bias/16.0;
      meta->dataset[iS].how.NEZV = meta->dataset[iS].how.NEZH + meta->dataset[iS].how.ZDR_bias;
 
      meta->dataset[iS].how.Vsamples=(double)inghdr->tcf.dsp.isamp;
      binmethod_avg = inghdr->tcf.rng.ibin_in_num/inghdr->tcf.rng.ibin_out_num;
      if(singlePRF) sprintf(meta->dataset[iS].how.azmethod,"NEAREST");
      else sprintf(meta->dataset[iS].how.azmethod,"AVERAGE");
      sprintf(meta->dataset[iS].how.elmethod,"NEAREST");
      meta->dataset[iS].how.binmethod_avg = binmethod_avg;
      if(binmethod_avg == 1) sprintf(meta->dataset[iS].how.binmethod,"NEAREST");
      else sprintf(meta->dataset[iS].how.binmethod,"AVERAGE");

      /* Average power is peakpwr*pulsewidth*(highprf+lowprf)/2, so the scale coeff. is 0.5*0.001 kW to get Watts */
      meta->dataset[iS].how.avgpwr=0.0005*meta->how.peakpwr * meta->dataset[iS].how.pulsewidth *
                                  (meta->dataset[iS].how.lowprf + meta->dataset[iS].how.highprf);

      sprintf(meta->dataset[iS].how.Dclutter,"%d",inghdr->tcf.dsp.idop_filter_first);
      DspStringFromPMode(meta->dataset[iS].how.ProcMode, inghdr->tcf.dsp.imajor_mode, NULL );
      DspStringFromPhaseMod( meta->dataset[iS].how.XMTphase , inghdr->tcf.dsp.iXmtPhaseSequence);
    }

    for( iQ=0 ; iQ < quantities ; iQ++ )
    {
         char thr[100];

         datatype=datatypes[iQ];
         if(iS==0)
         {
           if(datatype == DB_KDP   || 
              datatype == DB_RHOHV || 
              datatype == DB_SQI   || 
              datatype == DB_CCOR8 || 
              datatype == DB_PMI8)
           { 
              databytes[iQ]=2;
              ProcessDatatype(datatype,NULL);
           }  
           else
           { 
              databytes[iQ]=inghdrs[iQ].ibits_bin/8;
           }
         }

         switch(datatype)
         {

            default:
               sprintf(thr," ");
            break;

	    case DB_DBT: case DB_DBT2: case DB_DBTE8: case DB_DBTE16: case DB_DBTV8: case DB_DBTV16: 
              sprintf(thr,"%s",dsp_string_from_tcf(inghdr->tcf.cal.iuz_tcf,inghdr->tcf.cal.iuz_tcfMask));
            break;

            case DB_DBZ: case DB_DBZ2: case DB_DBZC: case DB_DBZC2: case DB_DBZE8: case DB_DBZE16: case DB_DBZV8: case DB_DBZV16:
              sprintf(thr,"%s",dsp_string_from_tcf(inghdr->tcf.cal.icz_tcf,inghdr->tcf.cal.icz_tcfMask));
            break;

            case DB_VEL: case DB_VEL2: case DB_VELC: case DB_VELC2:
              sprintf(thr,"%s",dsp_string_from_tcf(inghdr->tcf.cal.ivl_tcf,inghdr->tcf.cal.ivl_tcfMask));
            break;

            case DB_WIDTH: case DB_WIDTH2:
              sprintf(thr,"%s",dsp_string_from_tcf(inghdr->tcf.cal.iwd_tcf,inghdr->tcf.cal.iwd_tcfMask));
            break;

            case DB_ZDR: case DB_ZDR2:
              sprintf(thr,"%s",dsp_string_from_tcf(inghdr->tcf.cal.izdr_tcf,inghdr->tcf.cal.izdr_tcfMask));
            break;


	    /* Set polarization mode according to LDR data recorded if any */
            case DB_LDRH: case DB_LDRH2:
              sprintf(meta->dataset[iS].how.polmode,"LDR-H");
            break;

            case DB_LDRV: case DB_LDRV2:
              sprintf(meta->dataset[iS].how.polmode,"LDR-V");
            break;
         }

         /* Set filtering thresholds for quantities */
	 /*         printf("TCF mask %s\n",thr); */
         if(strstr(thr,"LOG")) 
            meta->dataset[iS].data[iQ].how.LOG = (double)inghdr->tcf.cal.izns_thr/16.0;
         else meta->dataset[iS].data[iQ].how.LOG = 0.0;

         if(strstr(thr,"CSR"))
            meta->dataset[iS].data[iQ].how.CSR = -(double)inghdr->tcf.cal.iccr_thr/16.0;
         else meta->dataset[iS].data[iQ].how.CSR = 1000.0;

         if(strstr(thr,"SQI")) 
            meta->dataset[iS].data[iQ].how.SQI = (double)inghdr->tcf.cal.isqi_thr/256.0;
         else meta->dataset[iS].data[iQ].how.SQI = 0.0;

         if(strstr(thr,"SIG"))
            meta->dataset[iS].data[iQ].how.SNRT = (double)inghdr->tcf.cal.isig_thr/16.0;
	 else meta->dataset[iS].data[iQ].how.SNRT = 0.0;

         if(strstr(thr,"PMI"))
            meta->dataset[iS].data[iQ].how.PMI=(double)inghdr->tcf.cal.ipmi_thr/256.0;
   	 else meta->dataset[iS].data[iQ].how.PMI=0.0;

         /* printf("Bytes of datatype %d=%d: %d\n",iQ,datatypes[iQ],databytes[iQ]); */
         meta->dataset[iS].data[iQ].what.bytes=databytes[iQ]; 

    }


    /* Read the data from each of the azimuth angles, and for each of the
     * quantities recorded.
     */
    azgates=inghdr->icf.irtotl;
    for( iAz=0 ; iAz < azgates ; iAz++ ) 
    {
       for( iQ=0 ; iQ < quantities ; iQ++ ) 
       {
         struct data_ray ray;
         int uncomp=1;
         SINT4 inlen, ioutlen; /* irecHold = irec_c, ioffHold = ioff_c ; */

         if(iQ==0 && IS_XHDR) uncomp=2; /* uncompress twice if XHDR present to skip it */
         do {
               uncompress_cowords( get_raw_bytes,
                                   pRaw->Record[0].PHeader.hdr.ibytes - ((irec_c * TAPE_RECORD_LEN) + ioff_c),
                                   &inlen, (SINT2 *)&ray, sizeof(ray)/2, &ioutlen ) ;
               uncomp--;
         } while(uncomp);

         if(!iAz)
         {
           /* scan size in bins: bins in ray * azimuth gates */ 
            scansize[iQ]=ray.hdr.ibincount * inghdr->icf.irtotl;
            scandata[iQ]=calloc(scansize[iQ],databytes[iQ]);
            /* fill the data array with 'undetect' */
            memset(scandata[iQ],255,scansize[iQ]*databytes[iQ]);
            /*            printf("bincount az:%d qu:%d=%d %d\n",iAz,iQ,datatypes[iQ],ray.hdr.ibincount); */

            meta->dataset[iS].where.nbins=ray.hdr.ibincount;
            if(iQ==0)
            {
               azdiff=fDegFromBin2(ray.hdr.iaz_end)-fPDegFromBin2(ray.hdr.iaz_start);
               if(azdiff<-azgates/2) azdiff+=azgates;
               if(azdiff>azgates/2) azdiff-=azgates;
               if(azdiff>0) rotsgn=1; else rotsgn=-1;
            }
         }
         N=iAz*ray.hdr.ibincount*databytes[iQ];
         if(!iQ)
         {
             raysecs = ray.hdr.itime;
             if(raysecs>=max_raysecs) 
             {
                max_raysecs=raysecs; 
                first_ray=iAz+rotsgn;
             }
	     /* V23: Start and stop azimuths of each ray */ 
             meta->dataset[iS].how.startazA[iAz] = fPDegFromBin2(ray.hdr.iaz_start); /* V23 */
             meta->dataset[iS].how.stopazA[iAz] = fPDegFromBin2(ray.hdr.iaz_end); /* V23 */
	     /*
             meta->dataset[iS].how.startelA[iAz]=fElDegFromBin2(ray.hdr.iel_start);
             meta->dataset[iS].how.stopelA[iAz]=fElDegFromBin2(ray.hdr.iel_end);
             */
         }

        /* If there is a ray here, then extract data. Otherwise
         * entire ray is missing and ray kept filled with 'undetect' value.
         */

         if( ioutlen > 0 )
         {

            CHANGE_QUANTITY_RESOLUTION = FALSE;
            datatype = datatypes[iQ];
            /*      printf("DATA %s AZ %d\n",sdata_name6(datatype),iAz); */
            switch( datatype ) 
            {

                /* ADD: run thru bins and convert to 16-bit */
                case DB_KDP:               /* KDP (1 byte) */
                {
                  long bI;
                  UINT2 newBinVal;
                  UINT1 *pBinVal = ray.data.iData1;
  
                  CHANGE_QUANTITY_RESOLUTION = TRUE;
                  newBinVal=1;
                  for(bI=0 ; bI < ray.hdr.ibincount; bI++,pBinVal++)
                  {
                    memcpy(&scandata[iQ][N+bI*2],&newBinVal,2); /* ADD CORRECT RELATION HERE */
                  }
                } break;

	        case DB_RHOHV: case DB_SQI: case DB_CCOR8: case DB_PMI8:  /* RhoHV etc (1 byte) */
                {
                  long bI;
                  UINT2 newBinVal;
                  UINT1 *pBinVal = ray.data.iData1;
  
                  CHANGE_QUANTITY_RESOLUTION = TRUE;
                  for(bI=0 ; bI < ray.hdr.ibincount; bI++,pBinVal++)
                  {
                    while(1)
                    {
                      if(*pBinVal==255) { newBinVal=65535; break; }
                      if(*pBinVal==0) { newBinVal=0; break; }
                      newBinVal=(UINT2)(1.0000001+65533.0*sqrt(((double)*pBinVal-1.0)/253.0));
                      break;
                    }
                    memcpy(&scandata[iQ][N+bI*2],&newBinVal,2);
                  }
                } break;
            }

            if(!CHANGE_QUANTITY_RESOLUTION)
            {
              /*               N=iAz*ray.hdr.ibincount*databytes[iQ]; */
               if(databytes[iQ]==1)
               {
                   memcpy(&scandata[iQ][N],ray.data.iData1,ray.hdr.ibincount*databytes[iQ]);
               }

               if(databytes[iQ]==2)
               {
                   memcpy(&scandata[iQ][N],ray.data.iData2,ray.hdr.ibincount*databytes[iQ]);
               }
            }
         } else 
         { 
           if(iQ==0) printf("RAY %d MISSING, SCAN %d\n",iAz,scan); 
         }
         totsize+=ray.hdr.ibincount*databytes[iQ];
       }  
    }

    if(first_ray==azgates) first_ray=0;
    if(first_ray<0) first_ray=azgates-1;

    meta->dataset[iS].how.antspeed=(double)rotsgn*fDegFromBin2(inghdr->tcf.scan.iscan_speed);
    meta->dataset[iS].how.rpm=meta->dataset[iS].how.antspeed/6.0;
    meta->dataset[iS].how.angres=(double)inghdr->tcf.scan.ires1000/1000.0;
    /*    printf("CALC a1gate %ld\n",(long)(first_az/meta->dataset[iS].how.angres)); */
    meta->dataset[iS].where.a1gate=first_ray;

    csecs=give_date_time(cdate,ctime,inghdrs[0].time);
    meta->dataset[iS].how.startepochs=csecs;
    sprintf(meta->dataset[iS].what.startdate,"%s",cdate);
    sprintf(meta->dataset[iS].what.starttime,"%s",ctime);
    csecs+=max_raysecs;
    meta->dataset[iS].how.endepochs=csecs;
    gmtime_r(&csecs,&Sdd);
    strftime(cdate,10,"%Y%m%d",&Sdd);
    strftime(ctime,10,"%H%M%S",&Sdd);

    sprintf(meta->dataset[iS].what.enddate,"%s",cdate);
    sprintf(meta->dataset[iS].what.endtime,"%s",ctime);

    if(VERB && !DUMPALL) DumpDatasetAttributes(iS);
    
    for( iQ=0 ; iQ < quantities ; iQ++ )
    { 
      if(VERB && !DUMPALL) DumpDataAttributes(iS,iQ);
       fwrite(&scandata[iQ][0],scansize[iQ],databytes[iQ],DATAF);
       free(scandata[iQ]);
    } 
    /* Done with this scan.  Discard the remainder of this block, if
     * any.
     */
    if( ioff_c ) { ioff_c = 0 ; irec_c++ ; }
  }

  {
     off_t filesize;
     unsigned char *totdata;

     filesize=ftell(DATAF);
     totdata=malloc(filesize);
     rewind(DATAF);
     fres=fread(totdata,filesize,1,DATAF);
     fwrite(meta,sizeof(MetaData),1,METAF);
     fwrite(totdata,filesize,1,METAF);
     /*     printf("%lu %lu\n",totsize,ftell(METAF)); */
     free(totdata);
  }
  fclose(DATAF);
  fclose(METAF);

  if(DUMPALL) DumpAllAttributes();

  return;
}



/* ================================================== */
/** Co-Routine to read the next run of bytes from the raw product file,
 * skipping and checking the record headers as we go.
 */
void get_raw_bytes( SINT2 *buf_a, SINT4 icnt_a )
{
  SINT4 icnt, iremain = icnt_a ; UINT1 *pbuf = (UINT1 *)buf_a ;

  while( iremain > 0 ) {
    /* If the offset is zero, then get the next record header and
     * compare the record numbers.  Exit if there is a mismatch.
     */
    if( ioff_c == 0 ) {
      struct raw_prod_bhdr bhdr ;
      memcpy( (void *)&bhdr, PRODPTR(irec_c, ioff_c), RAW_PROD_BHDR_SIZE ) ;
      if( bhdr.irec != irec_c ) {
        fprintf( stderr,  "Block header mismatch (%d) at block %d\n",
                 bhdr.irec, irec_c ) ; exit(1) ;
      }
      ioff_c += RAW_PROD_BHDR_SIZE ;
    }

    /* Extract as many bytes out of this record as we can.  Bump record
     * number and zero the offset if we read it all.
     */
    icnt = TAPE_RECORD_LEN - ioff_c ;
    if( icnt > iremain ) icnt = iremain ;

    memcpy( pbuf, PRODPTR(irec_c, ioff_c), icnt ) ;
    iremain -= icnt ; ioff_c += icnt ; pbuf += icnt ;

    if( ioff_c == TAPE_RECORD_LEN ) { ioff_c = 0 ; irec_c++ ; }
  }
}


void usage( void )
{
  printf( "\nCommand line options:\n" ) ;
  printf( " -h : Usage\n" ) ;
  printf( " -v : Verbose output\n" ) ;
  printf( " -d : Prints all meta data information\n" ) ;
  printf( " -s : Only Prints available quantities\n\n" ) ;
  printf( " After options the arguments are: input RAW file path, output file path\n");
  printf( " Example: ./IRIS_decoder -d $IRIS_PRODUCT_RAW/VAN101231235505.RAW1234 output_file.dat\n\n");
}

/** Assigns name and ODIM code to IRIS data type. If quantity index is given and datatypes pointer
    is NULL, only that datatype is processed. Otherwise all data types pointed with *datatypes 
    are processed. 
    The /datasetN/dataM/what of all scans 1-N are set to correct values M in MetaData structure.
*/ 
void ProcessDatatype(SINT4 quantity, SINT4 *datatypes)
{
   int iQ,iS,QN,datatype;
   DataWhat data_what;

   memset(&data_what,0,sizeof(DataWhat));  
   if(datatypes==NULL) QN=1; else QN=quantities;
      
   for(iQ=0 ; iQ < QN ; iQ++)
   {
           if(datatypes!=NULL) datatype = datatypes[iQ];
           else datatype=quantity;
 
           switch( datatype ) 
            {

                case DB_DBT:               /* Horiz total power (1 byte) */
                  data_what.QuantIdx = OQ_TH;
                  sprintf(data_what.quantity,"TH");
                break;

                case DB_DBZ:               /* Clutter Corrected reflectivity (1 byte) */
                  data_what.QuantIdx = OQ_DBZH;
                  sprintf(data_what.quantity,"DBZH");
                break;

                case DB_VEL:               /* Radial velocity (H) (1 byte) */
                  data_what.QuantIdx=OQ_VRADH;
                  sprintf(data_what.quantity,"VRADH");
                break;

                case DB_WIDTH:             /* Width (H) (1 byte) */
                  data_what.QuantIdx=OQ_WRADH;
                  sprintf(data_what.quantity,"WRADH");
                break;

                case DB_ZDR:               /* Differential reflectivity (1 byte) */
                  data_what.QuantIdx=OQ_ZDR;
                  sprintf(data_what.quantity,"ZDR");
                break;

                case DB_ZDRC:               /* Corrected differential reflectivity (1 byte) */
                  data_what.QuantIdx=OQ_ZDRC;
                  sprintf(data_what.quantity,"ZDRC");
                break;

                case DB_ZDRC2:               /* Corrected differential reflectivity (2 byte) */
                  data_what.QuantIdx=OQ_ZDRC2;
                  sprintf(data_what.quantity,"ZDRC");
                break;

                case DB_DBZC:              /* Fully corrected reflectivity (1 byte) */
                  data_what.QuantIdx=OQ_DBZHC;
                  sprintf(data_what.quantity,"DBZHC");
                break;

                case DB_DBT2:              /* Horiz uncorrected reflectivity (2 byte) */
                  data_what.QuantIdx=OQ_TH2;
                  sprintf(data_what.quantity,"TH");
                break;

                case DB_DBZ2:              /* Corrected reflectivity (2 byte) */
                  data_what.QuantIdx=OQ_DBZH2;
                  sprintf(data_what.quantity,"DBZH");
                break;

                case DB_VEL2:              /* Velocity (2 byte) */
                  data_what.QuantIdx=OQ_VRADH2;
                  sprintf(data_what.quantity,"VRADH");
                break;

                case DB_WIDTH2:            /* Width (2 byte) */
                  data_what.QuantIdx=OQ_WRADH2;
                  sprintf(data_what.quantity,"WRADH");
                break;

                case DB_ZDR2:              /* Differential reflectivity (2 byte) */
                  data_what.QuantIdx=OQ_ZDR2;
                  sprintf(data_what.quantity,"ZDR");
                break;

                case DB_KDP: case DB_KDP2:    /* Kdp (specific differential phase)(1 byte) */
                  data_what.QuantIdx=OQ_KDP2;
                  sprintf(data_what.quantity,"KDP");
                break;

                case DB_PHIDP:             /* PHIdp (differential phase)(1 byte) */
                  data_what.QuantIdx=OQ_PHIDP;
                  sprintf(data_what.quantity,"PHIDP");
                break;

                case DB_VELC:              /* Corrected Velocity (1 byte) */
                  data_what.QuantIdx=OQ_VRADDH;
                  sprintf(data_what.quantity,"VRADDH");
                break;

                case DB_SQI:               /* SQI (1 byte) converted to 2-byte */
                  data_what.QuantIdx=OQ_SQIH;
                  sprintf(data_what.quantity,"SQIH");
                break;

                case DB_RHOHV: case DB_RHOHV2:            /* RhoHV(0) (1 byte) */
                  data_what.QuantIdx=OQ_RHOHV2;
                  sprintf(data_what.quantity,"RHOHV");
                break;

                case DB_DBZC2:             /* Fully corrected reflectivity (2 byte) */
                  data_what.QuantIdx=OQ_DBZHC2;
                  sprintf(data_what.quantity,"DBZHC");
                break;

                case DB_VELC2:             /* Corrected Velocity (2 byte) */
                  data_what.QuantIdx=OQ_VRADDH2;
                  sprintf(data_what.quantity,"VRADDH");
                break;

                case DB_SQI2:              /* SQI (2 byte) */
                  data_what.QuantIdx=OQ_SQIH2;
                  sprintf(data_what.quantity,"SQIH");
                break;

                case DB_PHIDP2:            /* PHIdp (differential phase)(2 byte) */
                  data_what.QuantIdx=OQ_PHIDP2;
                  sprintf(data_what.quantity,"PHIDP");
                break;

                case DB_LDRH:              /* LDR H to V (1 byte) */
                  data_what.QuantIdx=OQ_LDR;
                  sprintf(data_what.quantity,"LDR");
                break;

                case DB_LDRH2:             /* LDR H to V (2 byte) */
                  data_what.QuantIdx=OQ_LDR2;
                  sprintf(data_what.quantity,"LDR");
                break;

                case DB_LDRV:              /* LDR V to H (1 byte) */
                  data_what.QuantIdx=OQ_LDRV;
                  sprintf(data_what.quantity,"LDRV");
                break;

                case DB_LDRV2:              /* LDR V to H (2 byte) */
                  data_what.QuantIdx=OQ_LDRV2;
                  sprintf(data_what.quantity,"LDRV");
                break;

                case DB_HCLASS:             /* HydroClass (1 byte) */
                  data_what.QuantIdx=OQ_HCLASS;
                  sprintf(data_what.quantity,"HCLASS");
                break;

                case DB_HCLASS2:             /* HydroClass (2 byte) */
                  data_what.QuantIdx=OQ_HCLASS2;
                  sprintf(data_what.quantity,"HCLASS");
                break;

                /* Enhanced total power (1 byte) */
                case DB_DBTV8:               
                  data_what.QuantIdx = OQ_TV;
                  sprintf(data_what.quantity,"TV");
		break; 

                /* Enhanced total power (2 byte) */
                case DB_DBTV16:              
                  data_what.QuantIdx=OQ_TV2;
                  sprintf(data_what.quantity,"TV");
                break;

                /* Enhanced total power (1 byte) */
                case DB_DBTE8:               
                  data_what.QuantIdx = OQ_TX;
                  sprintf(data_what.quantity,"TX");
		break; 

                /* Enhanced total power (2 byte) */
                case DB_DBTE16:              
                  data_what.QuantIdx=OQ_TX2;
                  sprintf(data_what.quantity,"TX");
                break;
		
                /* Enhanced reflectivity (1 byte) */
                case DB_DBZE8:              
                  data_what.QuantIdx=OQ_DBZX;
                  sprintf(data_what.quantity,"DBZX");
                break;

                /* Enhanced reflectivity (2 byte) */
                case DB_DBZE16:              
                  data_what.QuantIdx=OQ_DBZX2;
                  sprintf(data_what.quantity,"DBZX");
                break;

                /* Signal to noise ratio (1 byte) */
                case DB_SNR8:              
                  data_what.QuantIdx=OQ_SNR;
                  sprintf(data_what.quantity,"SNR");
                break;

                /* Signal to noise ratio (2 byte) */
                case DB_SNR16:              
                  data_what.QuantIdx=OQ_SNR2;
                  sprintf(data_what.quantity,"SNR");
                break;

                /* V-channel reflectivity (1 byte) */
                case DB_DBZV8:              
                  data_what.QuantIdx=OQ_DBZV;
                  sprintf(data_what.quantity,"DBZV");
                break;

                /* V-channel reflectivity (2 byte) */
                case DB_DBZV16:              
                  data_what.QuantIdx=OQ_DBZV2;
                  sprintf(data_what.quantity,"DBZV");
                break;

		/*------------------- V8.13.7 ---------------------*/

                /* PMI (Polarimetric meteo index) (1 byte) converted to 2-byte */
                case DB_PMI8:              
                  data_what.QuantIdx=OQ_PMI2;
                  sprintf(data_what.quantity,"PMI");
                break;

                /* PMI (Polarimetric meteo index) (2 byte) */
                case DB_PMI16:              
                  data_what.QuantIdx=OQ_PMI2;
                  sprintf(data_what.quantity,"PMI");


                break;
                /* The log receiver signal-to-noise ratio (1 byte) */
                case DB_LOG8:              
                  data_what.QuantIdx=OQ_LOG;
                  sprintf(data_what.quantity,"LOG");
                break;

                /* The log receiver signal-to-noise ratio (2 byte) */
                case DB_LOG16:              
                  data_what.QuantIdx=OQ_LOG2;
                  sprintf(data_what.quantity,"DBZV");
                break;

		/* Doppler channel clutter signal power (-CSR) (1 byte) */
                case DB_CSP8:              
                  data_what.QuantIdx=OQ_CSP;
                  sprintf(data_what.quantity,"CSP");
                break;

		/* Doppler channel clutter signal power (-CSR) (2 byte) */
                case DB_CSP16:              
                  data_what.QuantIdx=OQ_CSP2;
                  sprintf(data_what.quantity,"CSP");
                break;


                /* Cross correlation, uncorrected CCOR (1 byte) converted from 1-byte */
                case DB_CCOR8:              
                  data_what.QuantIdx=OQ_CCOR2;
                  sprintf(data_what.quantity,"CCOR");
                break;

                /* Cross correlation, uncorrected CCOR (2 byte) */
                case DB_CCOR16:              
                  data_what.QuantIdx=OQ_CCOR2;
                  sprintf(data_what.quantity,"CCOR");
                break;


                /* Attenuation of Zh (1 byte) */
                case DB_AH8:              
                  data_what.QuantIdx=OQ_PIA;
                  sprintf(data_what.quantity,"PIA");
                break;

                /* Attenuation of Zh (2 byte) */
                case DB_AH16:              
                  data_what.QuantIdx=OQ_PIA2;
                  sprintf(data_what.quantity,"PIA");
                break;


                /* Attenuation of Zv (1 byte) */
                case DB_AV8:              
                  data_what.QuantIdx=OQ_ATTV;
                  sprintf(data_what.quantity,"ATTV");
                break;

                /* Attenuation of Zv (2 byte) */
                case DB_AV16:              
                  data_what.QuantIdx=OQ_ATTV2;
                  sprintf(data_what.quantity,"ATTV");
                break;


                /* Attenuation of Zdr (1 byte) */
                case DB_AZDR8:              
                  data_what.QuantIdx=OQ_ATTZDR;
                  sprintf(data_what.quantity,"ATTZDR");
                break;

                /* Attenuation of Zdr (2 byte) */
                case DB_AZDR16:              
                  data_what.QuantIdx=OQ_ATTZDR2;
                  sprintf(data_what.quantity,"ATTZDR");
                break;
            }

           /* printf("%d %d %s %s\n",iQ,datatype,data_what.quantity,sdata_name6(datatype)); */
       for(iS=0;iS<scans;iS++)
       {
          if(datatypes == NULL)  meta->dataset[iS].data[quantity].what = data_what;
          else meta->dataset[iS].data[iQ].what = data_what;
       }
   }
}


time_t give_date_time(char *date, char *time, struct ymds_time ymds)
{
  int HH,MM,SS;

  HH=ymds.isec/3600;
  MM=(ymds.isec%3600)/60;
  SS=ymds.isec%60;

  sprintf(time,"%02d%02d%02d",HH,MM,SS);
  sprintf(date,"%4d%02d%02d",ymds.iyear,ymds.imon,ymds.iday);
  return(sec_from_date_time(date,time));
}

time_t sec_from_date_time(char *date, char *time)
{
   struct tm Sd,*Sdd;
   int y,m;
   time_t secs;
 
   sscanf(date,"%4d%2d%2d",&y,&m,&Sd.tm_mday);
   sscanf(time,"%2d%2d%2d",&Sd.tm_hour,&Sd.tm_min,&Sd.tm_sec);
   Sd.tm_year=y-1900;
   Sd.tm_mon=m-1;
   Sd.tm_isdst=0;
   Sdd=&Sd;
   secs=mktime(Sdd);
   return(secs);
}

void DumpCommonAttributes()
{
    printf("\nCommon attributes\n");
    printf("-----------------\n\n");

    printf("what.date          : %s\n",meta->what.date);
    printf("what.time          : %s\n",meta->what.time);
    printf("where.sitecode     : %s\n\n",meta->where.sitecode);

    printf("where.lon          : %.6f deg\n",meta->where.lon);
    printf("where.lat          : %.6f deg\n",meta->where.lat);
    printf("where.height       : %.1f m\n",meta->where.height);
    printf("where.towerheight  : %.1f m\n\n",meta->where.towerheight);

    printf("how.sw_version     : %s\n",meta->how.sw_version);
    printf("how.beamwidth      : %.3f deg\n",meta->how.beamwidth);
    printf("how.wavelength     : %.3f cm\n",meta->how.wavelength);
    if(meta->how.freeze < DBL_MAX)
    printf("how.freeze         : %.1f  km\n",meta->how.freeze);
    printf("how.peakpwr        : %.3f  kW\n",meta->how.peakpwr);
    printf("how.RAC            : %f dB/km\n",meta->how.RAC);
}

void DumpDatasetAttributes(int iS)
{
    SetWhat setwhat=meta->dataset[iS].what;
    SetWhere setwhere=meta->dataset[iS].where;
    How sethow=meta->dataset[iS].how;

    printf("\nScan #%d attributes\n",iS+1);
    printf("------------------\n\n");

    printf("dataset[%d].what.startdate       : %s\n",iS,setwhat.startdate);
    printf("dataset[%d].what.starttime       : %s\n",iS,setwhat.starttime);
    printf("dataset[%d].what.enddate         : %s\n",iS,setwhat.enddate);
    printf("dataset[%d].what.endtime         : %s\n\n",iS,setwhat.endtime);

    printf("dataset[%d].where.bin_elangle    : %ld\n",iS,(long)setwhere.bin_elangle);
    printf("dataset[%d].where.elangle        : %.2f deg\n",iS,setwhere.elangle);
    printf("dataset[%d].where.rstart         : %.3f km\n",iS,setwhere.rstart);
    printf("dataset[%d].where.rscale         : %.2f m\n",iS,setwhere.rscale);
    printf("dataset[%d].where.nrays          : %ld\n",iS,(long)setwhere.nrays);
    printf("dataset[%d].where.nbins          : %ld\n",iS,(long)setwhere.nbins);
    printf("dataset[%d].where.a1gate         : %ld\n\n",iS,(long)setwhere.a1gate);

    printf("dataset[%d].how.task             : %s\n",iS,sethow.task);
    printf("dataset[%d].how.binmethod_avg    : %ld\n",iS,(long)sethow.binmethod_avg);
    printf("dataset[%d].how.radhoriz         : %.2f km\n",iS,sethow.radhoriz);
    if(POL_H | POL_HV)
    printf("dataset[%d].how.MDSH (cal I0)    : %.3f dBm\n",iS,sethow.MDSH);
    if(POL_V | POL_HV)
    printf("dataset[%d].how.MDSV (cal I0)    : %.3f dBm\n",iS,sethow.MDSV);
    if(POL_H | POL_HV)
    printf("dataset[%d].how.radconstH        : %.3f dB\n",iS,sethow.radconstH);
    if(POL_V | POL_HV)
    printf("dataset[%d].how.radconstHV       : %.3f dB\n",iS,sethow.radconstHV);
    if(POL_HV)
    printf("dataset[%d].how.HVratio          : %.3f dBZ\n",iS,sethow.HVratio);
    printf("dataset[%d].how.NEZ              : %.3f dBZ\n",iS,sethow.NEZ);
    printf("dataset[%d].how.pulsewidth       : %.3f us\n",iS,sethow.pulsewidth);
    printf("dataset[%d].how.lowprf           : %.0f Hz\n",iS,sethow.lowprf);
    printf("dataset[%d].how.highprf          : %.0f Hz\n",iS,sethow.highprf);
    printf("dataset[%d].how.avgpwr           : %.1f W\n",iS,sethow.avgpwr);
    printf("dataset[%d].how.UnambVel         : %.3f m/s\n",iS,sethow.UnambVel);
    printf("dataset[%d].how.NI               : %.3f m/s\n",iS,sethow.NI);
    printf("dataset[%d].how.rpm              : %.3f RPM, %.3f deg/s\n",iS,sethow.rpm,sethow.rpm*6.0);
    printf("dataset[%d].how.angres           : %.3f deg\n",iS,sethow.angres);
    printf("dataset[%d].how.polarization     : %s\n",iS,sethow.polarization);
    printf("dataset[%d].how.Vsamples         : %ld\n",iS,(long)sethow.Vsamples);
    printf("dataset[%d].how.Dclutter         : filter #%s\n",iS,sethow.Dclutter);
    printf("dataset[%d].how.ProcMode         : %s\n",iS,sethow.ProcMode);
    printf("dataset[%d].how.XMTphase         : %s\n",iS,sethow.XMTphase);
    printf("dataset[%d].how.SQI              : %.2f\n",iS,sethow.SQI);
    printf("dataset[%d].how.CSR              : %.2f dB\n",iS,sethow.CSR);
    printf("dataset[%d].how.LOG              : %.2f dB\n",iS,sethow.LOG);
    printf("dataset[%d].how.SNRT              : %.2f dB\n",iS,sethow.SNRT);
    printf("dataset[%d].how.PMI              : %.2f\n\n",iS,sethow.PMI);
}

void DumpDataAttributes(int iS, int iQ)
{
    DataWhat datawhat=meta->dataset[iS].data[iQ].what;
    DataHow datahow=meta->dataset[iS].data[iQ].how;

    printf("\nQuantity index %s attributes\n",datawhat.quantity);
    printf("------------------------------\n\n");
  
    printf("dataset[%d].data[%d].what.quantity   : %s\n",iS,iQ,datawhat.quantity);
    printf("dataset[%d].data[%d].what.QuantIdx   : %d\n",iS,iQ,datawhat.QuantIdx);
    printf("dataset[%d].data[%d].what.bytes      : %d\n",iS,iQ,datawhat.bytes);

    if(datahow.SQI > 0.0)
    printf("dataset[%d].data[%d].how.SQI         : %.2f\n",iS,iQ,datahow.SQI);
    if(datahow.CSR < 999.0)
    printf("dataset[%d].data[%d].how.CSR         : %.2f dB\n",iS,iQ,datahow.CSR);
    if(datahow.LOG > 0.0)
    printf("dataset[%d].data[%d].how.LOG         : %.2f dB\n",iS,iQ,datahow.LOG);
    if(datahow.SNRT > 0.0)
    printf("dataset[%d].data[%d].how.SNRT         : %.2f dB\n",iS,iQ,datahow.SNRT);
    if(datahow.PMI > 0.0)
    printf("dataset[%d].data[%d].how.PMI         : %.2f\n",iS,iQ,datahow.PMI);
    printf("\n");
}

void DumpAllAttributes()
{
   int iS,iQ;

   DumpCommonAttributes();
   for(iS=0;iS<scans;iS++)
   { 
        DumpDatasetAttributes(iS);
        for(iQ=0;iQ<quantities;iQ++) DumpDataAttributes(iS,iQ);
   }
}   

