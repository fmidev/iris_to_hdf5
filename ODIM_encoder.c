/*! \file ODIM_encoder.c
\brief This code converts the intermediate radar data/metadata file generated
by IRIS_decoder.c to ODIM HDF5 format (see OPERA Working Document WD_2008_03).

The code uses hlhdf library provided by SMHI to convert metadata organized as
C-structures (see ODIM_struct.h) to HDF5. The hlhdf library is under LGPL licensing.

All program arguments after options are the pathnames to the intermediate
file(s). Program accepts one option -v for verbose output.
The program generates a HDF5 file named according to ODIM specification to output
directory given as ODIM_OUTPUT_DIR environment variable, or user specified filename
as ODIM_OUTPUT_FILE variable. In both cases the ODIM style filename is written to
file specified as ODIM_OUTPUT_DIR/ODIM_NAME_FILE. If user don't set the ODIM_NAME_FILE
environment variable, the default filename is ODIM_filename.txt. All other
input needed are also given as environment variables (see test.sh).
*/

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "ODIM_struct.h"

# define uchar unsigned char
# define FALSE 0
# define TRUE  1

hid_t H5out;
char groupattr[2000],A1,A2,outname[200],*origcenter,timestamp[100];
QuantCfg QCF[OQ_QUANTS_TOTAL];
uchar POL_H,POL_V,POL_HV,ALL_QUANTS=0;
short wanted_scanquants[MAX_SCANS][MAX_QUANTS]={{0}};
short VERB=FALSE,QUIET=FALSE;
size_t fres;
char boolstr[2][6]={"False","True"};
char flagname[2][16][50]={{{0}}};
char *envp;

/*-----------------------------------------------------------------------------------------*/
/** \brief Adds any HDF5 scalar numeric attribute named <I>*attr</I> to group named <I>*group</I> and sets it to value <I>val</I>, with wanted type */
int  add_attr_numeric_to_group(hid_t group, char *attr, void *val, hid_t type);
/** \brief Adds dataset to group */
hid_t add_dataset_to_group(hid_t group, char *name, int compress_level,int bytes, int rank, hsize_t *dims, void *data);
/** \brief Gives quantity code of ODIM quantity name <I>*Qstr</I> */
short getQuantityCode(char *Qstr);
/** \brief sets parameters (gain, offset, nodata, undetect) of all quantities */
void SetQuantityParams();
/** \brief Parses wanted quantities for encoding from ODIM_</I>SITE</I>_quantities environment
variable.<BR> <I>SITE</I> is the IRIS radar site code.<BR> The codes of quantities found are
put to <I>wanted_scanquants[scan index][quantity index]</I> table. */ 
void get_wanted_quantities(char *Qstr);

int main(int argc, char** argv)
{

  uchar *in_scandata=NULL,*outdata=NULL; 
  int S,Q,fI,iS,iQ,tS;
  int64_t vol_scan_number,scans,scans_total=0;
  /*  unsigned char strattr[1000]; */
  char *outdir=NULL,*outfile=NULL,*odimname=NULL,*compress_str=NULL;
  int compresslevel;
  char ODIM_namestr[200];
  char def_outdir[2]=".";
  char datagroup[200],
       setgroup[200];

  short argF,last_Q=0,radnum=0;
  /*  unsigned char quantflags[32]={0}; */
  FILE *METAF;
  char envname[1000]={0},sitecode[4]={0};
  MetaData *meta;

  RootWhat in_what;
  RootWhere in_where;
  How in_how;

  hid_t H5out=-1,G_root_what=-1,G_root_how=-1,G_root_where=-1;
  hsize_t scandims[2];

  /*-----------------------------------------------------------------------------------------*/

  setbuf(stdout,NULL);
  SetQuantityParams();
  origcenter=getenv("ODIM_ORIGCENTER");
  outdir=getenv("ODIM_OUTPUT_DIR");
  if(outdir==NULL) outdir=def_outdir;
  outfile=getenv("ODIM_OUTPUT_FILE");
  odimname=getenv("ODIM_NAME_FILE");
  compress_str=getenv("ODIM_COMPRESSION_LEVEL");
  if(compress_str==NULL) compresslevel=6; else compresslevel=atoi(compress_str);
  if(compresslevel < 0 || compresslevel > 9) compresslevel=6;

  argF=1;
  {
    int i;

    for(i=1;i<argc;i++) if(argv[i][0]=='-')
    { 
       if(argv[i][1]=='v') VERB = 1;
       if(argv[i][1]=='q') QUIET = 1;
       argF++;
    }
  }

  /* set the names of IRIS flag attributes */
  sprintf(flagname[0][0],"f_speckle_Z");
  sprintf(flagname[0][2],"f_speckle_V");
  sprintf(flagname[0][4],"f_rangenorm");
  sprintf(flagname[0][5],"f_beginpulse");
  sprintf(flagname[0][6],"f_endpulse");
  sprintf(flagname[0][7],"f_varpulses_dPRF");
  sprintf(flagname[0][8],"f_3lag_PP02");
  sprintf(flagname[0][9],"f_shipVcorr");
  sprintf(flagname[0][10],"f_unfold_Vc");
  sprintf(flagname[0][11],"f_fallsp_Vc");
  sprintf(flagname[0][12],"f_beamblock_Zc");
  sprintf(flagname[0][13],"f_Z_atten_Zc");
  sprintf(flagname[0][14],"f_targdetect_Zc");
  sprintf(flagname[0][15],"f_stormrel_Vc");
  sprintf(flagname[1][0],"f_dp_atten_Zc+ZDRc");
  sprintf(flagname[1][1],"f_dp_atten_Z+ZDR");

 
  meta=calloc(1,sizeof(MetaData));  

  vol_scan_number=0; /* index of scan really written to h5-file (origin 1 = ODIM scan_index) */
  /* looping thru data from subtasks */
  for(fI=argF; fI < argc; fI++)
  {

    /* argv[1] is the volume HDF5 file, all others are IRIS metadata/data files */
     METAF=fopen(argv[fI],"r");
     fres=fread(&meta[0],1,sizeof(MetaData),METAF);
     scans=meta->scans;
     if(VERB) printf("\n=========================================================================================\n");
     if(VERB) printf("\nFile %s, having %ld scans \n",argv[fI],(long)scans);

     in_what=meta->what;
     in_where=meta->where;
     in_how=meta->how;

     sprintf(sitecode,"%.3s",in_where.sitecode);
     /*    if(VERB) printf("SITECODE %s\n",sitecode); */

     if(!scans_total) /* common metadata for whole volume is taken from the first scan or subvolume */
     {
       {
         /* reading wanted quantities */
         char *Wstr=NULL;

         sprintf(envname,"ODIM_%s_quantities",sitecode);
         Wstr=getenv(envname);
         get_wanted_quantities(Wstr);
       }
       /*       for(S=0;S<wanted_quants;S++)if(VERB) printf("%s\n",wanted_quantarr[S]); */

        sprintf(timestamp,"%s%s",meta->dataset[0].what.startdate,meta->dataset[0].what.starttime);
        sprintf(envname,"ODIM_%s_source",sitecode);
        sprintf(in_what.source,"%s",getenv(envname));
        radnum=atoi(strstr(in_what.source,"RAD:")+6);

       /* !!! outname pitää olla tmpname, koska lopullista tiedostonimeä ei voi
          luoda, ennen kuin kaikki PPI:t on käyty läpi, eli A1 ja A2 on generoitu !!! 
          Siispä loppuun rename !!!
       */

       if(outfile != NULL) sprintf(outname,"%s/%s",outdir,outfile);
       else
       {
          /* construct the preliminary Odyssey filename */
          sprintf(ODIM_namestr,"T_PAxx%02d_C_%s_%s.h5",radnum,origcenter,timestamp);
          sprintf(outname,"%s/%s",outdir,ODIM_namestr);
       }

       H5out=H5Fcreate(outname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
       H5LTset_attribute_string(H5out,"/","Conventions",getenv("ODIM_Conventions"));
       G_root_what=H5Gcreate2(H5out,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
       G_root_where=H5Gcreate2(H5out,"where",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
       G_root_how=H5Gcreate2(H5out,"how",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

      /* ROOT GROUP /what   */
       H5LTset_attribute_string(H5out,"what","version",getenv("ODIM_what_version"));
       H5LTset_attribute_string(H5out,"what","date",in_what.date);  
       H5LTset_attribute_string(H5out,"what","time",in_what.time);
       H5LTset_attribute_string(H5out,"what","source",in_what.source);

      /* ROOT GROUP /where   */
       add_attr_numeric_to_group(G_root_where,"lon",&in_where.lon,H5T_NATIVE_DOUBLE);
       add_attr_numeric_to_group(G_root_where,"lat",&in_where.lat,H5T_NATIVE_DOUBLE);
       add_attr_numeric_to_group(G_root_where,"height",&in_where.height,H5T_NATIVE_DOUBLE);

      /* ROOT GROUP /how   */
        H5LTset_attribute_string(H5out,"how","system",in_how.system);
        H5LTset_attribute_string(H5out,"how","software","IRIS");
        H5LTset_attribute_string(H5out,"how","sw_version",in_how.sw_version);
        H5LTset_attribute_string(H5out,"how","TXtype",in_how.TXtype);
        add_attr_numeric_to_group(G_root_how,"beamwH",&in_how.beamwH,H5T_NATIVE_DOUBLE); /* V23 */
        add_attr_numeric_to_group(G_root_how,"beamwV",&in_how.beamwV,H5T_NATIVE_DOUBLE); /* V23 */
        add_attr_numeric_to_group(G_root_how,"wavelength",&in_how.wavelength,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_root_how,"antgainH",&in_how.antgainH,H5T_NATIVE_DOUBLE); /* V23 */
        add_attr_numeric_to_group(G_root_how,"antgainV",&in_how.antgainV,H5T_NATIVE_DOUBLE); /* V23 */
        /* towerheight moved from /root/where */
        add_attr_numeric_to_group(G_root_how,"towerheight",&in_where.towerheight,H5T_NATIVE_DOUBLE);
        H5LTset_attribute_string(H5out,"how","poltype",in_how.poltype); /* V23 */
        H5LTset_attribute_string(H5out,"how","simulated",getenv("ODIM_how_simulated"));

        /* common quality attributes */
        add_attr_numeric_to_group(G_root_how,"freeze",&in_how.freeze,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_root_how,"RXlossH",&in_how.RXlossH,H5T_NATIVE_DOUBLE); /* V23 */
        add_attr_numeric_to_group(G_root_how,"RXlossV",&in_how.RXlossV,H5T_NATIVE_DOUBLE); /* V23 */
        add_attr_numeric_to_group(G_root_how,"TXlossH",&in_how.TXlossH,H5T_NATIVE_DOUBLE); /* V23 */
        add_attr_numeric_to_group(G_root_how,"TXlossV",&in_how.TXlossV,H5T_NATIVE_DOUBLE); /* V23 */
        add_attr_numeric_to_group(G_root_how,"FiniteBandwithLoss",&in_how.CWloss,H5T_NATIVE_DOUBLE); /* IRIS specific */
        if(in_how.radomelossH > 0.0)
           add_attr_numeric_to_group(G_root_how,"radomelossH",&in_how.radomelossH,H5T_NATIVE_DOUBLE); /* V23 */
        if(in_how.radomelossV > 0.0)
           add_attr_numeric_to_group(G_root_how,"radomelossV",&in_how.radomelossV,H5T_NATIVE_DOUBLE); /* V23 */
        if(in_how.pointaccEL > 0.0)
           add_attr_numeric_to_group(G_root_how,"pointaccEL",&in_how.pointaccEL,H5T_NATIVE_DOUBLE);
        if(in_how.pointaccAZ > 0.0)
           add_attr_numeric_to_group(G_root_how,"pointaccAZ",&in_how.pointaccAZ,H5T_NATIVE_DOUBLE);
        if(in_how.malfunc[0])
           H5LTset_attribute_string(H5out,"how","malfunc",in_how.malfunc);
        if(in_how.radar_msg[0])
           H5LTset_attribute_string(H5out,"how","radar_msg",in_how.radar_msg);
        if(in_how.dynrange > 0.0)
           add_attr_numeric_to_group(G_root_how,"dynrange",&in_how.dynrange,H5T_NATIVE_DOUBLE);
        if(in_how.OUR > 0.0) 
           add_attr_numeric_to_group(G_root_how,"OUR",&in_how.OUR,H5T_NATIVE_DOUBLE);
        if(in_how.comment[0]) 
           H5LTset_attribute_string(H5out,"how","comment",in_how.comment);
        add_attr_numeric_to_group(G_root_how,"RAC",&in_how.RAC,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_root_how,"gasattn",&in_how.RAC,H5T_NATIVE_DOUBLE);

	/* add_attr_numeric_to_group(G_root_how,"RXbandwidth",&in_how.RXbandwidth,H5T_NATIVE_DOUBLE); per scan ? V23 */
        /* add_attr_numeric_to_group(G_root_how,"avgpwr",&in_how.avgpwr,H5T_NATIVE_DOUBLE);  per scan ? */
        /* add_attr_numeric_to_group(G_root_how,"antgain",&in_how.antgain,H5T_NATIVE_DOUBLE); obsolete V23 */
	/* add_attr_numeric_to_group(G_root_how,"radconstH",&in_how.radconstH,H5T_NATIVE_DOUBLE); per scan */
        /* add_attr_numeric_to_group(G_root_how,"radconstV",&in_how.radconstV,H5T_NATIVE_DOUBLE);  V23 */
        /* add_attr_numeric_to_group(G_root_how,"radconstHV",&in_how.radconstHV,H5T_NATIVE_DOUBLE);  obsolete V23 */
        /* add_attr_numeric_to_group(G_root_how,"TXpower",&in_how.TXpower,H5T_NATIVE_DOUBLE); array V23 */
        /* add_attr_numeric_to_group(G_root_how,"NI",&in_how.NI,H5T_NATIVE_DOUBLE); per scan */
        /* add_attr_numeric_to_group(G_root_how,"Vsamples",&in_how.Vsamples,H5T_NATIVE_LLONG); per scan */
        /* add_attr_numeric_to_group(G_root_how,"radhoriz",&in_how.radhoriz,H5T_NATIVE_DOUBLE); per scan */

     }

  /*------------------ looping thru scans ----------------------------------------------*/

     for(S=1;S<=scans;S++)
     {
        SetWhat in_setwhat;    
        SetWhere in_setwhere;    
        How in_sethow; 
        long nbins,nrays,quantities;
        ulong insize,outsize;
        short DPOL,eQ=0,acc_quants,wanted_quants_in_scan;
        short outbytes,AQ,WQ,aq,wq; 
        double relangle;
        hid_t G_dataset,G_dataset_what,G_dataset_where,G_dataset_how;

        iS=S-1; /* scan index of read data array */
        tS=S+scans_total; /* dataset index */ 

        in_setwhat=meta->dataset[iS].what;
        in_setwhere=meta->dataset[iS].where;
        in_sethow=meta->dataset[iS].how;

        wanted_quants_in_scan = wanted_scanquants[tS][0];

        if(!wanted_quants_in_scan) continue;
        acc_quants=0;
        for(wq=0;wq<MAX_QUANTS;wq++)
        {
          WQ=wanted_scanquants[tS][wq];
          if(!WQ) break;
          if(WQ<0) { acc_quants++; break; }
          for(aq=0;aq<MAX_QUANTS;aq++)
          {
             AQ=meta->dataset[iS].data[aq].what.QuantIdx;
             /*if(VERB) printf("AQ %d WQ %d\n",AQ,WQ); */
             if((WQ == AQ) || (WQ == AQ+TWOB) || (AQ == WQ+TWOB)) { acc_quants++; break; }
             if(!AQ) break;
          }
        } 
        if(!acc_quants) continue;

        vol_scan_number++;
	if(VERB) printf("\n\nENCODING SCAN #%d\n=====================================================\n",(int)vol_scan_number);

        POL_H=in_sethow.POL_H;
        POL_V=in_sethow.POL_V;
        POL_HV=in_sethow.POL_HV;

        if(strchr(in_sethow.polarization,'H') && 
           strchr(in_sethow.polarization,'V')) DPOL=1;
  
        /* DATASETs */
        sprintf(setgroup,"/dataset%d",(int)vol_scan_number);
        G_dataset=H5Gcreate2(H5out,setgroup,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        G_dataset_what=H5Gcreate2(G_dataset,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        G_dataset_where=H5Gcreate2(G_dataset,"where",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        G_dataset_how=H5Gcreate2(G_dataset,"how",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        
        /* /datasetS/what attributes */
        H5LTset_attribute_string(G_dataset,"what","product","SCAN");
        H5LTset_attribute_string(G_dataset,"what","startdate",in_setwhat.startdate);  
        H5LTset_attribute_string(G_dataset,"what","starttime",in_setwhat.starttime);  
        H5LTset_attribute_string(G_dataset,"what","enddate",in_setwhat.enddate);  
        H5LTset_attribute_string(G_dataset,"what","endtime",in_setwhat.endtime);  

        /* /datasetS/where attributes */
        relangle=((double)(int)((in_setwhere.elangle+0.01)*10.0))/10.0;
        add_attr_numeric_to_group(G_dataset_where,"elangle",&relangle,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_where,"nbins",&in_setwhere.nbins,H5T_NATIVE_LLONG);
        add_attr_numeric_to_group(G_dataset_where,"rstart",&in_setwhere.rstart,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_where,"rscale",&in_setwhere.rscale,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_where,"nrays",&in_setwhere.nrays,H5T_NATIVE_LLONG);
        add_attr_numeric_to_group(G_dataset_where,"a1gate",&in_setwhere.a1gate,H5T_NATIVE_LLONG);

        /* /datasetS/how attributes */
        add_attr_numeric_to_group(G_dataset_how,"scan_index",&vol_scan_number,H5T_NATIVE_LLONG);  
        add_attr_numeric_to_group(G_dataset_how,"startepochs",&in_sethow.startepochs,H5T_NATIVE_DOUBLE);  
        add_attr_numeric_to_group(G_dataset_how,"endepochs",&in_sethow.endepochs,H5T_NATIVE_DOUBLE);  
        H5LTset_attribute_string(G_dataset,"how","task",in_sethow.task);
        add_attr_numeric_to_group(G_dataset_how,"radhoriz",&in_sethow.radhoriz,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_how,"rpm",&in_sethow.rpm,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_how,"antspeed",&in_sethow.antspeed,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_how,"pulsewidth",&in_sethow.pulsewidth,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_how,"lowprf",&in_sethow.lowprf,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_how,"highprf",&in_sethow.highprf,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_how,"nomTXpower",&in_sethow.nomTXpower,H5T_NATIVE_DOUBLE);  /* = peakpwr */
        add_attr_numeric_to_group(G_dataset_how,"peakpwr",&in_sethow.nomTXpower,H5T_NATIVE_DOUBLE);
        add_attr_numeric_to_group(G_dataset_how,"avgpwr",&in_sethow.avgpwr,H5T_NATIVE_DOUBLE);
        H5LTset_attribute_string(G_dataset,"how","polarization",in_sethow.polarization); /* "H,V" used when polmode is "simultaneous-dual", H|V if switched */
        H5LTset_attribute_string(G_dataset,"how","Dclutter",in_sethow.Dclutter);
        H5LTset_attribute_string(G_dataset,"how","ProcMode",in_sethow.ProcMode);
        H5LTset_attribute_string(G_dataset,"how","XMTphase",in_sethow.XMTphase);
	add_attr_numeric_to_group(G_dataset_how,"RXbandwidth",&in_sethow.RXbandwidth,H5T_NATIVE_DOUBLE); /* per scan ? V23 */
        /* add_attr_numeric_to_group(G_dataset_how,"MDS",&in_sethow.MDSH,H5T_NATIVE_DOUBLE); */
        add_attr_numeric_to_group(G_dataset_how,"Vsamples",&in_sethow.Vsamples,H5T_NATIVE_LLONG);  
        add_attr_numeric_to_group(G_dataset_how,"NI",&in_sethow.NI,H5T_NATIVE_DOUBLE);  
        /*      add_attr_numeric_to_group(G_dataset_how,"NEW:NyqWidth",&in_sethow.NEW_NyqWidth,H5T_NATIVE_DOUBLE); */  
        if(POL_H | POL_HV) add_attr_numeric_to_group(G_dataset_how,"radconstH",&in_sethow.radconstH,H5T_NATIVE_DOUBLE);  
        add_attr_numeric_to_group(G_dataset_how,"binmethod_avg",&in_sethow.binmethod_avg,H5T_NATIVE_LLONG);  
        H5LTset_attribute_string(G_dataset,"how","binmethod",in_sethow.binmethod);  

        H5LTset_attribute_double(G_dataset,"how","startazA",in_sethow.startazA,in_setwhere.nrays);  
        H5LTset_attribute_double(G_dataset,"how","stopazA",in_sethow.stopazA,in_setwhere.nrays);

        H5LTset_attribute_string(G_dataset,"how","polmode",in_sethow.polmode); /* V23 */
        if(POL_H | POL_HV) /* MDSH/V not in spec yet */
        { 
          add_attr_numeric_to_group(G_dataset_how,"NEZH",&in_sethow.NEZH,H5T_NATIVE_DOUBLE);  
	  add_attr_numeric_to_group(G_dataset_how,"MDSH",&in_sethow.MDSH,H5T_NATIVE_DOUBLE);
        }

        if(POL_V | POL_HV)
        { 
          add_attr_numeric_to_group(G_dataset_how,"NEZV",&in_sethow.NEZV,H5T_NATIVE_DOUBLE);  
	  add_attr_numeric_to_group(G_dataset_how,"MDSV",&in_sethow.MDSV,H5T_NATIVE_DOUBLE);
        }

        if(POL_HV)
        {
           add_attr_numeric_to_group(G_dataset_how,"HVratio",&in_sethow.HVratio,H5T_NATIVE_DOUBLE);  
           add_attr_numeric_to_group(G_dataset_how,"radconstHV",&in_sethow.radconstHV,H5T_NATIVE_DOUBLE);  
           add_attr_numeric_to_group(G_dataset_how,"zdrcal",&in_sethow.ZDR_bias,H5T_NATIVE_DOUBLE);
           add_attr_numeric_to_group(G_dataset_how,"LDR_bias",&in_sethow.LDR_bias,H5T_NATIVE_DOUBLE);
        }

        /* Set the "scan_optimized" quantity according to single/dual PRF */
        { 
	  char optquant[20];

          if(in_sethow.lowprf != in_sethow.highprf) sprintf(optquant,"VRADH");
	     else sprintf(optquant,"DBZH");
          H5LTset_attribute_string(G_dataset,"how","scan_optimized",optquant);
	}

        /* set IRIS flags */
        {
  	   uint16_t flagI,fpos,flagbit;

           for(flagI=0;flagI<2;flagI++)
           {
             for(fpos=0;fpos<16;fpos++)
             {
               if(flagname[flagI][fpos][0] != 0)
               {
                 flagbit=((1<<fpos) & in_sethow.Cflags[flagI]) ? 1:0;
                 H5LTset_attribute_string(G_dataset,"how",flagname[flagI][fpos],boolstr[flagbit]);
               }
             }
	   }
	}
    


        H5Gclose(G_dataset_what);
        H5Gclose(G_dataset_where);
        H5Gclose(G_dataset_how);

        nrays=in_setwhere.nrays;
        nbins=in_setwhere.nbins;
        scandims[0]=(hsize_t)nrays;
        scandims[1]=(hsize_t)nbins;

        quantities=meta->dataset[iS].quantities;
        for(Q=1;Q<=quantities;Q++) /* looping thru quantities of a scan */
        {
           DataWhat in_datawhat;
           DataHow in_datahow;
           int binbytes,iW;
           short Encode=0;
           ulong iN,oN;
           ushort W;
           unsigned char B;
           double c_offset,c_gain,eps=1.0e-6,fB;
           double wanted_gain,wanted_offset,avail_gain,avail_offset;
           short avail_Q,wanted_Q,wanted_quants;

           iQ=Q-1;
           Encode=0;
           in_datawhat=meta->dataset[iS].data[iQ].what;
           in_datahow=meta->dataset[iS].data[iQ].how;
           binbytes=in_datawhat.bytes;
           /* if(VERB) printf("%s %d\n",QCF[in_datawhat.QuantIdx].in_quantity,binbytes); */
           insize=nrays*nbins*binbytes;
           in_scandata=malloc(insize);
           fres=fread(in_scandata,insize,1,METAF);
           /* compare in_datawhat.quantity and wanted quantities */

           if(ALL_QUANTS) wanted_quants=1; else wanted_quants=MAX_QUANTS;
           for(iW=0;iW<wanted_quants;iW++)
           {

             avail_Q=in_datawhat.QuantIdx;
             /*     if(VERB) printf("Available %s\n",
                     QCF[avail_Q].in_quantity);  */
             if(ALL_QUANTS) wanted_Q=avail_Q; else wanted_Q=wanted_scanquants[tS][iW]; 
             if(!wanted_Q) break;
             if(wanted_Q<0) { wanted_Q=avail_Q; iW=wanted_quants; }
             /*if(VERB) printf("Q PAIR:  Available %s, wanted %d:%s\n",
                QCF[avail_Q].in_quantity,iW,QCF[wanted_Q].in_quantity); */

             if((wanted_Q==OQ_SQIH && avail_Q==OQ_SQIH2) ||
                (wanted_Q==OQ_PMI && avail_Q==OQ_PMI2) ||
                (wanted_Q==OQ_CCOR && avail_Q==OQ_CCOR2) ||
                (wanted_Q==OQ_RHOHV && avail_Q==OQ_RHOHV2))
                { 
                   if(VERB) printf("\nNOTICE: linear conversion to 8-bit %s not possible,\nusing forced 16-bit output\n",
                            QCF[wanted_Q].quantity);
                    wanted_Q+=TWOB;
                }

             if(avail_Q == wanted_Q) 
             { 
                Encode=1;
                outbytes=binbytes;
                break; 
             }
             if(wanted_Q == avail_Q-TWOB)
             { 
                Encode=8;
                outbytes=1;
                break;  
             }

             if(wanted_Q == avail_Q+TWOB) 
             { 
                Encode=16;
                outbytes=2;
                break;  
             }
           }
           outsize=nbins*nrays*outbytes;

           if(Encode)
           {
                avail_gain=QCF[avail_Q].gain;
                wanted_gain=QCF[wanted_Q].gain;
                avail_offset=QCF[avail_Q].offset;
                wanted_offset=QCF[wanted_Q].offset;
                /*    if(VERB) printf("%s: wG = %f, wF = %f\n",QCF[wanted_Q].in_quantity,wanted_gain,wanted_offset); */
           } else { if(VERB) printf("SKIPPING %s\n---------------------\n",QCF[avail_Q].in_quantity); continue; }


           /* If conversion between 8/16 bit data is requested, the new gain and offset are calculated */
           if(Encode>1)
           { 
                outdata=malloc(outsize); 

                if(wanted_Q == OQ_VRADH2)
                {
                  avail_gain *= in_sethow.NI;
                  avail_offset *= in_sethow.NI;
                  /* if(VERB) printf("Nyq = %f, wG = %f, wF = %f\n",in_sethow.NEW_NyqVel,wanted_gain,wanted_offset);
                   */
                }
                if(wanted_Q == OQ_VRADH)
                {
                  avail_gain /= in_sethow.NI;
                  avail_offset /= in_sethow.NI;
                  /* if(VERB) printf("Nyq = %f, wG = %f, wF = %f\n",in_sethow.NEW_NyqVel,wanted_gain,wanted_offset);
                   */
                }

                if(wanted_Q == OQ_WRADH2)
                {
                  wanted_gain *= in_sethow.NyqWidth;
                  wanted_offset *= in_sethow.NyqWidth;
                }
                if(wanted_Q == OQ_WRADH)
                {
                  wanted_gain /= in_sethow.NyqWidth;
                  wanted_offset /= in_sethow.NyqWidth;
                }

                c_gain = avail_gain/wanted_gain;
                c_offset = eps + (avail_offset-wanted_offset)/wanted_gain;

                if(wanted_Q == OQ_VRADH)
                {
                  wanted_gain *= in_sethow.NI;
                  wanted_offset *= in_sethow.NI;
                }
                if(wanted_Q == OQ_WRADH)
                {
                  wanted_gain *= in_sethow.NyqWidth;
                  wanted_offset *= in_sethow.NyqWidth;
                }

               if(VERB) printf("\nWRITING %d-bit %s, conversion from %s\n",Encode,QCF[wanted_Q].in_quantity,QCF[avail_Q].in_quantity); 
               /* printf("%f %f\n",c_gain,c_offset); */
           }

           if(Encode==1)
           { 
              outdata=in_scandata;
              if(binbytes==1)
	      {
	         if(avail_Q==OQ_VRADH)
	         {
		    wanted_gain *= in_sethow.NI;
		    wanted_offset *= in_sethow.NI;
	         }

	         if(avail_Q==OQ_WRADH)
	         {
		    wanted_gain *= in_sethow.NyqWidth;
		    wanted_offset *= in_sethow.NyqWidth;
	         }
	      }
              /*  memcpy(&outdata[0],&in_scandata[0],outsize); */ 
             if(VERB) printf("\nWRITING %s\n",QCF[wanted_Q].in_quantity);
           } 

           /* If conversion between 8/16 bit data is requested, new output quantity values are calculated */
           if(Encode==8)
           {
             oN=0;
             for(iN=0;iN<insize;iN+=binbytes)
             {
                memcpy(&W,in_scandata+iN,2);
                while(1)
                {                
                  if(W == (ushort)QCF[avail_Q].nodata)   { B = (uchar)QCF[wanted_Q].nodata;   break; }
                  if(W == (ushort)QCF[avail_Q].undetect) { B = (uchar)QCF[wanted_Q].undetect; break; }
                  fB=c_gain*(double)W+c_offset;
                  if(fB<0) { B=0; break; }
                  if(fB>254) { B=254; break; }
                  B=(uchar)fB;
                  break;
                }
                outdata[oN]=B;
                oN++;
             }
           }     

           if(Encode==16)
           {
	     /*printf("%d %d\n",(ushort)QCF[wanted_Q].nodata,(ushort)QCF[wanted_Q].undetect);
	       printf("%f %f\n",c_gain,c_offset); */
             oN=0;
             for(iN=0;iN<insize;iN+=binbytes)
             {
                B=in_scandata[iN];
                while(1)
                {                
                  if(B == (uchar)QCF[avail_Q].undetect) { W = (ushort)QCF[wanted_Q].undetect; break; }
                  if(B == (uchar)QCF[avail_Q].nodata)   { W = (ushort)QCF[wanted_Q].nodata;   break; }
                  W=(ushort)(c_gain*(double)B+c_offset);
                  break;
                }
                memcpy(outdata+oN,&W,2);
                oN+=2;
             }
           }     

           if(Encode)
           { 
	       hid_t G_data,D_data,G_datawhat,G_datahow; /* ,G_datawhere; */

               QuantCfg out_datawhat = QCF[wanted_Q];
               double wanted_nodata,wanted_undetect;

               eQ++;
               last_Q=wanted_Q;
               sprintf(datagroup,"%s/data%d",setgroup,eQ);
               if(VERB) printf("DATAGROUP %s\n",datagroup);
               if(VERB) printf("___________________________________\n");

               G_data=H5Gcreate2(H5out,datagroup,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
               D_data=add_dataset_to_group(G_data,"data",compresslevel,outbytes,2,scandims,outdata);

               /* /datasetS/dataQ/data attributes */   
               if(outbytes==1)
               {
		  H5LTset_attribute_string(G_data,"data","CLASS","IMAGE");  
		  H5LTset_attribute_string(G_data,"data","IMAGE_VERSION","1.2");  
               }
               /* /datasetS/dataQ/what attributes */
               G_datawhat=H5Gcreate2(G_data,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
               G_datahow=H5Gcreate2(G_data,"how",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 
               wanted_nodata=(double)out_datawhat.nodata;
               wanted_undetect=(double)out_datawhat.undetect;

               add_attr_numeric_to_group(G_datawhat,"gain",&wanted_gain,H5T_NATIVE_DOUBLE);
               add_attr_numeric_to_group(G_datawhat,"nodata",&wanted_nodata,H5T_NATIVE_DOUBLE);
               add_attr_numeric_to_group(G_datawhat,"offset",&wanted_offset,H5T_NATIVE_DOUBLE);
               H5LTset_attribute_string(G_data,"what","quantity",out_datawhat.quantity);  
               add_attr_numeric_to_group(G_datawhat,"undetect",&wanted_undetect,H5T_NATIVE_DOUBLE);  

               if(in_datahow.SQI>0) add_attr_numeric_to_group(G_datahow,"SQI",&in_datahow.SQI,H5T_NATIVE_DOUBLE);
               if(in_datahow.CSR<1000.0)  add_attr_numeric_to_group(G_datahow,"CSR",&in_datahow.CSR,H5T_NATIVE_DOUBLE);
               if(in_datahow.LOG>0)  add_attr_numeric_to_group(G_datahow,"LOG",&in_datahow.LOG,H5T_NATIVE_DOUBLE);  
               if(in_datahow.SNRT>0)  add_attr_numeric_to_group(G_datahow,"SNRT",&in_datahow.SNRT,H5T_NATIVE_DOUBLE);  
               if(in_datahow.PMI>0)  add_attr_numeric_to_group(G_datahow,"PMI",&in_datahow.PMI,H5T_NATIVE_DOUBLE);  


               H5Dclose(D_data);
               H5Gclose(G_datawhat);
               H5Gclose(G_datahow);
               H5Gclose(G_data);

               if(Encode>1) free(outdata);
          }
          free(in_scandata);
        }

        if(eQ>1) A1='Z';
          else
          {
            switch (last_Q)
            {
               case OQ_DBZH: case OQ_DBZH2:
                 A1='G';
               break;

               case OQ_VRADH: case OQ_VRADH2: case OQ_VRADDH: case OQ_VRADDH2:
                 A1='H';
               break;

               case OQ_WRADH: case OQ_WRADH2:
                 A1='I';
               break;

               default:
                 A1='X';
               break;
            }
          }
     }
     scans_total+=scans;
     if(VERB) printf("%d scans total done\n",(int)scans_total);
  }

  if(!vol_scan_number) goto fail;

  add_attr_numeric_to_group(G_root_how,"scan_count",&scans_total,H5T_NATIVE_LLONG); /* scans total V23 */

  if(vol_scan_number>1)
  { 
     H5LTset_attribute_string(H5out,"what","object","PVOL");
     A2='Z';
  }
  else
  {  
     H5LTset_attribute_string(H5out,"what","object","SCAN");
     A2='A';
  }

   

  if(VERB)  printf("\n======================== HDF5 CREATED ==============================\n");

  /* construct the final Odyssey filename */
  if(outfile==NULL)
  {
          FILE *ODIM_NAME;
          char odimpath[500];
  
          sprintf(ODIM_namestr,"T_PA%c%c%02d_C_%s_%s.h5",A1,A2,radnum,origcenter,timestamp);
          if(odimname) 
	  {
             sprintf(odimpath,"%s",odimname);
             ODIM_NAME=fopen(odimpath,"w");
             fprintf(ODIM_NAME,"%s/%s",outdir,ODIM_namestr);
             fclose(ODIM_NAME);
	  }
  }

  H5Gclose(G_root_what);
  H5Gclose(G_root_where);
  H5Gclose(G_root_how);
  H5Fclose(H5out);

  /* rename the h5 file if ODIM name convention is used */
  if(outfile==NULL) 
  {
     char finalname[300]={0};

     sprintf(finalname,"%s/%s",outdir,ODIM_namestr);
     rename(outname,finalname);
     if(!QUIET) printf("%s\n",finalname);
  } 
  else if(!QUIET) printf("%s\n",outname);

  return(0);
 fail:
  if(VERB) printf("\n!!!!!!!!!!!!!!!!!  NO SUITABLE DATA FOR ENCODING !!!!!!!!!!!!!!!!!\n\n");
  return(1);
}

/*=============================================================================*/

int  add_attr_numeric_to_group(hid_t group, char *attr, void *val, hid_t type)
{
   hid_t Attr,space;
   int ret;

   space  = H5Screate(H5S_SCALAR);
   Attr = H5Acreate2(group, attr, type, space, H5P_DEFAULT, H5P_DEFAULT);
   ret = H5Awrite(Attr, type, val);
   H5Sclose(space);
   H5Aclose(Attr);
   return(ret);
}


hid_t add_dataset_to_group(hid_t group, char *name, int compress_level,int bytes, int rank, hsize_t *dims, void *data)
{
     hid_t dataspace,plist,dset,dtype=0;

     if(bytes==1) dtype=H5Tcopy(H5T_NATIVE_UCHAR);
     if(bytes==2) dtype=H5Tcopy(H5T_NATIVE_USHORT);

     dataspace=H5Screate_simple(rank, dims, NULL);
     plist = H5Pcreate(H5P_DATASET_CREATE);
     H5Pset_chunk(plist, rank, dims);
     H5Pset_deflate( plist, compress_level);
     dset = H5Dcreate2(group, name, dtype, dataspace,
            H5P_DEFAULT, plist, H5P_DEFAULT);
     H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
     H5Pclose(plist);
     H5Sclose(dataspace);

     return(dset);
  }

short getQuantityCode(char *Qstr)
{
   short Qi;

   for(Qi=0; Qi<OQ_QUANTS_TOTAL; Qi++)
   {
     if(strcmp(Qstr,QCF[Qi].in_quantity)==0) return(Qi);
   }
   return(0);
}

void SetQuantityParams()
{

  short i;
 
  for(i=0;i<OQ_QUANTS_TOTAL;i++)
  {
    memset( &QCF[i],0,sizeof(QuantCfg));
    QCF[i].in_quantity[0]=' ';
  }    

  sprintf(QCF[OQ_TH].in_quantity,"TH");
  sprintf(QCF[OQ_TH].quantity,"TH");
  QCF[OQ_TH].gain=0.5;
  QCF[OQ_TH].offset=-32;
  QCF[OQ_TH].nodata=255;
  QCF[OQ_TH].undetect=0;

  sprintf(QCF[OQ_TV].in_quantity,"TV");
  sprintf(QCF[OQ_TV].quantity,"TV");
  QCF[OQ_TV].gain=0.5;
  QCF[OQ_TV].offset=-32;
  QCF[OQ_TV].nodata=255;
  QCF[OQ_TV].undetect=0;

  sprintf(QCF[OQ_TX].in_quantity,"TX");
  sprintf(QCF[OQ_TX].quantity,"TX");
  QCF[OQ_TX].gain=0.5;
  QCF[OQ_TX].offset=-32;
  QCF[OQ_TX].nodata=255;
  QCF[OQ_TX].undetect=0;

  sprintf(QCF[OQ_DBZH].in_quantity,"DBZH");
  sprintf(QCF[OQ_DBZH].quantity,"DBZH");
  QCF[OQ_DBZH].gain=0.5;
  QCF[OQ_DBZH].offset=-32.0;
  QCF[OQ_DBZH].nodata=255;
  QCF[OQ_DBZH].undetect=0;

  sprintf(QCF[OQ_DBZV].in_quantity,"DBZV");
  sprintf(QCF[OQ_DBZV].quantity,"DBZV");
  QCF[OQ_DBZV].gain=0.5;
  QCF[OQ_DBZV].offset=-32.0;
  QCF[OQ_DBZV].nodata=255;
  QCF[OQ_DBZV].undetect=0;

  sprintf(QCF[OQ_DBZX].in_quantity,"DBZX");
  sprintf(QCF[OQ_DBZX].quantity,"DBZX");
  QCF[OQ_DBZX].gain=0.5;
  QCF[OQ_DBZX].offset=-32.0;
  QCF[OQ_DBZX].nodata=255;
  QCF[OQ_DBZX].undetect=0;

  sprintf(QCF[OQ_SNR].in_quantity,"SNR");
  sprintf(QCF[OQ_SNR].quantity,"SNR");
  QCF[OQ_SNR].gain=0.5;
  QCF[OQ_SNR].offset=-32.0;
  QCF[OQ_SNR].nodata=255;
  QCF[OQ_SNR].undetect=0;

  sprintf(QCF[OQ_LOG].in_quantity,"LOG");
  sprintf(QCF[OQ_LOG].quantity,"LOG");
  QCF[OQ_LOG].gain=0.5;
  QCF[OQ_LOG].offset=-32.0;
  QCF[OQ_LOG].nodata=255;
  QCF[OQ_LOG].undetect=0;

  sprintf(QCF[OQ_CSP].in_quantity,"CSP");
  sprintf(QCF[OQ_CSP].quantity,"CSP");
  QCF[OQ_CSP].gain=0.5;
  QCF[OQ_CSP].offset=-32.0;
  QCF[OQ_CSP].nodata=255;
  QCF[OQ_CSP].undetect=0;

  sprintf(QCF[OQ_ATTV].in_quantity,"ATTV");
  sprintf(QCF[OQ_ATTV].quantity,"ATTV");
  QCF[OQ_ATTV].gain=0.5;
  QCF[OQ_ATTV].offset=-32.0;
  QCF[OQ_ATTV].nodata=255;
  QCF[OQ_ATTV].undetect=0;

  sprintf(QCF[OQ_PIA].in_quantity,"PIA");
  sprintf(QCF[OQ_PIA].quantity,"PIA");
  QCF[OQ_PIA].gain=0.5;
  QCF[OQ_PIA].offset=-32.0;
  QCF[OQ_PIA].nodata=255;
  QCF[OQ_PIA].undetect=0;

  sprintf(QCF[OQ_ATTZDR].in_quantity,"ATTZDR");
  sprintf(QCF[OQ_ATTZDR].quantity,"ATTZDR");
  QCF[OQ_ATTZDR].gain=0.5;
  QCF[OQ_ATTZDR].offset=-32.0;
  QCF[OQ_ATTZDR].nodata=255;
  QCF[OQ_ATTZDR].undetect=0;

  sprintf(QCF[OQ_DBZHC].in_quantity,"DBZHC");
  sprintf(QCF[OQ_DBZHC].quantity,"DBZHC");
  QCF[OQ_DBZHC].gain=0.5;
  QCF[OQ_DBZHC].offset=-32.0;
  QCF[OQ_DBZHC].nodata=255;
  QCF[OQ_DBZHC].undetect=0;

  sprintf(QCF[OQ_VRADH].in_quantity,"VRADH");
  sprintf(QCF[OQ_VRADH].quantity,"VRADH");
  QCF[OQ_VRADH].gain=1.0/127.0;
  QCF[OQ_VRADH].offset=-128.0/127.0;
  QCF[OQ_VRADH].nodata=0;
  QCF[OQ_VRADH].undetect=0;

  sprintf(QCF[OQ_VRADH2].in_quantity,"VRADH2");
  sprintf(QCF[OQ_VRADH2].quantity,"VRADH");
  QCF[OQ_VRADH2].gain=0.01;
  QCF[OQ_VRADH2].offset=-327.68;
  QCF[OQ_VRADH2].nodata=65535;
  QCF[OQ_VRADH2].undetect=0;

  sprintf(QCF[OQ_VRADDH].in_quantity,"VRADDH");
  sprintf(QCF[OQ_VRADDH].quantity,"VRADDH");
  QCF[OQ_VRADDH].gain=75.0/127.0;
  QCF[OQ_VRADDH].offset=-(150.0/127.0+74.4);
  QCF[OQ_VRADDH].nodata=255;
  QCF[OQ_VRADDH].undetect=0;

  sprintf(QCF[OQ_VRADDH2].in_quantity,"VRADDH2");
  sprintf(QCF[OQ_VRADDH2].quantity,"VRADDH");
  QCF[OQ_VRADDH2].gain=0.01;
  QCF[OQ_VRADDH2].offset=-327.68;
  QCF[OQ_VRADDH2].nodata=65535;
  QCF[OQ_VRADDH2].undetect=0;

  sprintf(QCF[OQ_WRADH].in_quantity,"WRADH");
  sprintf(QCF[OQ_WRADH].quantity,"WRADH");
  QCF[OQ_WRADH].gain=1.0/256.0; 
  QCF[OQ_WRADH].offset=0;
  QCF[OQ_WRADH].nodata=255;
  QCF[OQ_WRADH].undetect=0;

  sprintf(QCF[OQ_WRADH2].in_quantity,"WRADH2");
  sprintf(QCF[OQ_WRADH2].quantity,"WRADH");
  QCF[OQ_WRADH2].gain=0.01;
  QCF[OQ_WRADH2].offset=0;
  QCF[OQ_WRADH2].nodata=65535;
  QCF[OQ_WRADH2].undetect=0;

  sprintf(QCF[OQ_ZDR].in_quantity,"ZDR");
  sprintf(QCF[OQ_ZDR].quantity,"ZDR");
  QCF[OQ_ZDR].gain=1.0/16.0;
  QCF[OQ_ZDR].offset=128.0/16.0;
  QCF[OQ_ZDR].nodata=255;
  QCF[OQ_ZDR].undetect=0;

  sprintf(QCF[OQ_ZDRC].in_quantity,"ZDRC");
  sprintf(QCF[OQ_ZDRC].quantity,"ZDRC");
  QCF[OQ_ZDRC].gain=1.0/16.0;
  QCF[OQ_ZDRC].offset=128.0/16.0;
  QCF[OQ_ZDRC].nodata=255;
  QCF[OQ_ZDRC].undetect=0;

  sprintf(QCF[OQ_TH2].in_quantity,"TH2");
  sprintf(QCF[OQ_TH2].quantity,"TH");
  QCF[OQ_TH2].gain=0.01;
  QCF[OQ_TH2].offset=-327.68;
  QCF[OQ_TH2].nodata=65535;
  QCF[OQ_TH2].undetect=0;

  sprintf(QCF[OQ_TV2].in_quantity,"TV2");
  sprintf(QCF[OQ_TV2].quantity,"TV");
  QCF[OQ_TV2].gain=0.01;
  QCF[OQ_TV2].offset=-327.68;
  QCF[OQ_TV2].nodata=65535;
  QCF[OQ_TV2].undetect=0;

  sprintf(QCF[OQ_TX2].in_quantity,"TX2");
  sprintf(QCF[OQ_TX2].quantity,"TX");
  QCF[OQ_TX2].gain=0.01;
  QCF[OQ_TX2].offset=-327.68;
  QCF[OQ_TX2].nodata=65535;
  QCF[OQ_TX2].undetect=0;

  sprintf(QCF[OQ_DBZH2].in_quantity,"DBZH2");
  sprintf(QCF[OQ_DBZH2].quantity,"DBZH");
  QCF[OQ_DBZH2].gain=0.01;
  QCF[OQ_DBZH2].offset=-327.68;
  QCF[OQ_DBZH2].nodata=65535;
  QCF[OQ_DBZH2].undetect=0;

  sprintf(QCF[OQ_DBZV2].in_quantity,"DBZV2");
  sprintf(QCF[OQ_DBZV2].quantity,"DBZV");
  QCF[OQ_DBZV2].gain=0.01;
  QCF[OQ_DBZV2].offset=-327.68;
  QCF[OQ_DBZV2].nodata=65535;
  QCF[OQ_DBZV2].undetect=0;

  sprintf(QCF[OQ_DBZX2].in_quantity,"DBZX2");
  sprintf(QCF[OQ_DBZX2].quantity,"DBZX");
  QCF[OQ_DBZX2].gain=0.01;
  QCF[OQ_DBZX2].offset=-327.68;
  QCF[OQ_DBZX2].nodata=65535;
  QCF[OQ_DBZX2].undetect=0;

  sprintf(QCF[OQ_SNR2].in_quantity,"SNR2");
  sprintf(QCF[OQ_SNR2].quantity,"SNR");
  QCF[OQ_SNR2].gain=0.01;
  QCF[OQ_SNR2].offset=-327.68;
  QCF[OQ_SNR2].nodata=65535;
  QCF[OQ_SNR2].undetect=0;

  sprintf(QCF[OQ_LOG2].in_quantity,"LOG2");
  sprintf(QCF[OQ_LOG2].quantity,"LOG");
  QCF[OQ_LOG2].gain=0.01;
  QCF[OQ_LOG2].offset=-327.68;
  QCF[OQ_LOG2].nodata=65535;
  QCF[OQ_LOG2].undetect=0;

  sprintf(QCF[OQ_CSP2].in_quantity,"CSP2");
  sprintf(QCF[OQ_CSP2].quantity,"CSP");
  QCF[OQ_CSP2].gain=0.01;
  QCF[OQ_CSP2].offset=-327.68;
  QCF[OQ_CSP2].nodata=65535;
  QCF[OQ_CSP2].undetect=0;

  sprintf(QCF[OQ_ATTV2].in_quantity,"ATTV2");
  sprintf(QCF[OQ_ATTV2].quantity,"ATTV");
  QCF[OQ_ATTV2].gain=0.01;
  QCF[OQ_ATTV2].offset=-327.68;
  QCF[OQ_ATTV2].nodata=65535;
  QCF[OQ_ATTV2].undetect=0;

  sprintf(QCF[OQ_PIA2].in_quantity,"PIA2");
  sprintf(QCF[OQ_PIA2].quantity,"PIA");
  QCF[OQ_PIA2].gain=0.01;
  QCF[OQ_PIA2].offset=-327.68;
  QCF[OQ_PIA2].nodata=65535;
  QCF[OQ_PIA2].undetect=0;

  sprintf(QCF[OQ_ATTZDR2].in_quantity,"ATTZDR2");
  sprintf(QCF[OQ_ATTZDR2].quantity,"ATTZDR");
  QCF[OQ_ATTZDR2].gain=0.01;
  QCF[OQ_ATTZDR2].offset=-327.68;
  QCF[OQ_ATTZDR2].nodata=65535;
  QCF[OQ_ATTZDR2].undetect=0;

  sprintf(QCF[OQ_DBZHC2].in_quantity,"DBZHC2");
  sprintf(QCF[OQ_DBZHC2].quantity,"DBZHC");
  QCF[OQ_DBZHC2].gain=0.01;
  QCF[OQ_DBZHC2].offset=-327.68;
  QCF[OQ_DBZHC2].nodata=65535;
  QCF[OQ_DBZHC2].undetect=0;

  sprintf(QCF[OQ_ZDR2].in_quantity,"ZDR2");
  sprintf(QCF[OQ_ZDR2].quantity,"ZDR");
  QCF[OQ_ZDR2].gain=0.01;
  QCF[OQ_ZDR2].offset=-327.68;
  QCF[OQ_ZDR2].nodata=65535;
  QCF[OQ_ZDR2].undetect=0;

  sprintf(QCF[OQ_ZDRC2].in_quantity,"ZDRC2");
  sprintf(QCF[OQ_ZDRC2].quantity,"ZDRC");
  QCF[OQ_ZDRC2].gain=0.01;
  QCF[OQ_ZDRC2].offset=-327.68;
  QCF[OQ_ZDRC2].nodata=65535;
  QCF[OQ_ZDRC2].undetect=0;

  sprintf(QCF[OQ_KDP2].in_quantity,"KDP2");
  sprintf(QCF[OQ_KDP2].quantity,"KDP");
  QCF[OQ_KDP2].gain=0.01;
  QCF[OQ_KDP2].offset=-327.68;
  QCF[OQ_KDP2].nodata=65535;
  QCF[OQ_KDP2].undetect=0;

  sprintf(QCF[OQ_PHIDP].in_quantity,"PHIDP");
  sprintf(QCF[OQ_PHIDP].quantity,"PHIDP");
  QCF[OQ_PHIDP].gain=180.0/254.0;
  QCF[OQ_PHIDP].offset=-QCF[OQ_PHIDP].gain;
  QCF[OQ_PHIDP].nodata=255;
  QCF[OQ_PHIDP].undetect=0;

  sprintf(QCF[OQ_SQIH2].in_quantity,"SQIH2");
  sprintf(QCF[OQ_SQIH2].quantity,"SQIH");
  QCF[OQ_SQIH2].gain=1.0/65533.0;
  QCF[OQ_SQIH2].offset=-QCF[OQ_SQIH2].gain;
  QCF[OQ_SQIH2].nodata=65535;
  QCF[OQ_SQIH2].undetect=0;

  sprintf(QCF[OQ_PMI2].in_quantity,"PMI2");
  sprintf(QCF[OQ_PMI2].quantity,"PMI");
  QCF[OQ_PMI2].gain=1.0/65533.0;
  QCF[OQ_PMI2].offset=-QCF[OQ_PMI2].gain;
  QCF[OQ_PMI2].nodata=65535;
  QCF[OQ_PMI2].undetect=0;

  sprintf(QCF[OQ_CCOR2].in_quantity,"CCOR2");
  sprintf(QCF[OQ_CCOR2].quantity,"CCOR");
  QCF[OQ_CCOR2].gain=1.0/65533.0;
  QCF[OQ_CCOR2].offset=-QCF[OQ_CCOR2].gain;
  QCF[OQ_CCOR2].nodata=65535;
  QCF[OQ_CCOR2].undetect=0;

  sprintf(QCF[OQ_RHOHV2].in_quantity,"RHOHV2");
  sprintf(QCF[OQ_RHOHV2].quantity,"RHOHV");
  QCF[OQ_RHOHV2].gain=1.0/65533.0;
  QCF[OQ_RHOHV2].offset=-QCF[OQ_RHOHV2].gain;
  QCF[OQ_RHOHV2].nodata=65535;
  QCF[OQ_RHOHV2].undetect=0;

  sprintf(QCF[OQ_RHOHV].in_quantity,"RHOHV");
  sprintf(QCF[OQ_RHOHV].quantity,"RHOHV");
  QCF[OQ_RHOHV].gain=1.0/65533.0;
  QCF[OQ_RHOHV].offset=-QCF[OQ_RHOHV].gain;
  QCF[OQ_RHOHV].nodata=65535;
  QCF[OQ_RHOHV].undetect=0;

  sprintf(QCF[OQ_PHIDP2].in_quantity,"PHIDP2");
  sprintf(QCF[OQ_PHIDP2].quantity,"PHIDP");
  QCF[OQ_PHIDP2].gain=360.0/65534.0;
  QCF[OQ_PHIDP2].offset=-QCF[OQ_PHIDP2].gain;
  QCF[OQ_PHIDP2].nodata=65535;
  QCF[OQ_PHIDP2].undetect=0;

  sprintf(QCF[OQ_LDR].in_quantity,"LDR");
  sprintf(QCF[OQ_LDR].quantity,"LDR");
  QCF[OQ_LDR].gain=0.2;
  QCF[OQ_LDR].offset=-QCF[OQ_LDR].gain-45.0;
  QCF[OQ_LDR].nodata=255;
  QCF[OQ_LDR].undetect=0;

  sprintf(QCF[OQ_LDR2].in_quantity,"LDR2");
  sprintf(QCF[OQ_LDR2].quantity,"LDR");
  QCF[OQ_LDR2].gain=0.01;
  QCF[OQ_LDR2].offset=-327.68;
  QCF[OQ_LDR2].nodata=65535;
  QCF[OQ_LDR2].undetect=0;

  sprintf(QCF[OQ_LDRV].in_quantity,"LDRV");
  sprintf(QCF[OQ_LDRV].quantity,"LDRV");
  QCF[OQ_LDRV].gain=0.2;
  QCF[OQ_LDRV].offset=-QCF[OQ_LDRV].gain-45.0;
  QCF[OQ_LDRV].nodata=255;
  QCF[OQ_LDRV].undetect=0;

  sprintf(QCF[OQ_LDRV2].in_quantity,"LDRV2");
  sprintf(QCF[OQ_LDRV2].quantity,"LDRV");
  QCF[OQ_LDRV2].gain=0.01;
  QCF[OQ_LDRV2].offset=-327.68;
  QCF[OQ_LDRV2].nodata=65535;
  QCF[OQ_LDRV2].undetect=0;

  sprintf(QCF[OQ_HCLASS].in_quantity,"HCLASS");
  sprintf(QCF[OQ_HCLASS].quantity,"HCLASS");
  QCF[OQ_HCLASS].gain=1.0;
  QCF[OQ_HCLASS].offset=0.0;
  QCF[OQ_HCLASS].nodata=255;
  QCF[OQ_HCLASS].undetect=0;

  sprintf(QCF[OQ_HCLASS2].in_quantity,"HCLASS2");
  sprintf(QCF[OQ_HCLASS2].quantity,"HCLASS");
  QCF[OQ_HCLASS2].gain=1.0;
  QCF[OQ_HCLASS2].offset=0.0;
  QCF[OQ_HCLASS2].nodata=65535;
  QCF[OQ_HCLASS2].undetect=0;


}


void get_wanted_quantities(char *Qstr)
{
         /* reading wanted quantities */
          char *votoken, *swtoken, *qutoken, *vostr, *swstr, *qustr;
          char *vo_saveptr, *sw_saveptr, *qu_saveptr; 
          char volim[2]=" ",swlim[2]=":",qulim[2]="," ;
          short Qcode,vo,sw,qu,swcount,swI,swlist[MAX_SCANS]={0};
          /*      char qulist[MAX_SCANS][MAX_QUANTS][20]={{{0}}}; */

          if(!Qstr || strcmp(Qstr,"ALL")==0 || strcmp(Qstr,"*")==0   || strcmp(Qstr,"*:*")==0) 
          {
            ALL_QUANTS=1;
            for(sw=1;sw<MAX_SCANS;sw++) wanted_scanquants[sw][0]=-1;
           if(VERB) printf("\nAll found quantities will be encoded\n\n");
            return;
          } else if(VERB) printf("\nSEARCHING: %s\n\n",Qstr);

          for(vo=1,vostr=&Qstr[0];;vo++,vostr=NULL)
          {
               swcount=0;
               votoken=strtok_r(vostr,volim,&vo_saveptr);
               if(votoken == NULL) break;
               for(sw=1,swstr=votoken;;sw++,swstr=NULL)
               {
                  swtoken=strtok_r(swstr,swlim,&sw_saveptr);
                  if(swtoken == NULL) break;
                  for(qu=1,qustr=swtoken;;qu++,qustr=NULL)
                  {
                    qutoken=strtok_r(qustr,qulim,&qu_saveptr);
                    if(qutoken==NULL) break;
                    if(sw==1) 
                    { 
                      if(qutoken[0]=='*')
                      { 
                         swcount=MAX_SCANS;
                         for(swI=0;swI<swcount;swI++) swlist[swI]=swI;
                      }
                      else
                      {
                        swlist[swcount]=atoi(qutoken);
                        swcount++;
                      }
                    } else for(swI=0;swI<swcount;swI++)
                             { 
                               if(qutoken[0]=='*') Qcode=-1; else
                                  Qcode=getQuantityCode(qutoken);
                               wanted_scanquants[swlist[swI]][qu-1]=Qcode;
                             }
                  }
               }
          }
          /*
          for(sw=1;sw<MAX_SCANS;sw++) for(qu=0;qu<MAX_QUANTS;qu++)
          {
            if(wanted_scanquants[sw][qu-1]) 
             if(VERB) printf("SCAN %d %d %d %s\n",sw,qu,wanted_scanquants[sw][qu-1],QCF[wanted_scanquants[sw][qu-1]].in_quantity);
          }
          */
 }

