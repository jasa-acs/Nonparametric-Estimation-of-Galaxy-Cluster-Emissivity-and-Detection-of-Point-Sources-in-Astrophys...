#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include "math.h"
/*#include "TFile.h"
#include "TH1F.h"*/

const double pi=3.14159265359;
using namespace std;

int findbin(double dist,double *bins){
    int bin=0;
    while (bins[bin]<dist) {
        bin++;
    }
    return bin;
}

int line_num(char *filename){
    FILE *f=fopen(filename, "r");;
    char c;
    int lines = 0;
    if(f == NULL) return 0;
    while((c = fgetc(f)) != EOF){
        if(c == '\n') lines++;
    }
    fclose(f);
    return lines;
}

int col_num(char *filename){
    ifstream fin(filename);
    string Line;
    stringstream iss; //declare the stream
    double num; //This will hold the data, it could be an array if you like.  It can be whatever type you need.
    int ncol=0;
    
    while (getline(fin, Line))  // Put each line into variable "Line"
    {
        iss << Line; // Put the line into a string stream;
        
        while (iss.good()) // while the stream is not empty
        {
            iss >> num; //get data from the stream. This will give you only up until the next whitespace.
            //Get numbers, strings, whatever you want.
            ncol++;
        }
    }
    return ncol;
}

double CalcMedian(int nval,double *array)
{
    double median;
    sort(array, array+nval);
    
    if (nval  % 2 == 0)
    {
        median = (array[nval / 2 - 1] + array[nval / 2]) / 2;
    }
    else
    {
        median = array[nval / 2];
    }
    
    return median;
}

void logbinning(double binsize,double maxrad,int nbin,double *bins,double *ebins,int &newnbin){
    do {
        newnbin=0;
        double binedge=binsize/60./2.;
        double db=0.0;
        int i=0;
        while (db<binsize/60.) {
            bins[i]=(i+0.5)*binsize/60.;
            ebins[i]=binsize/60./2.;
            double bb=nbin/log10(maxrad/binsize*60.)*(log10(binedge)-log10(binsize/60.));
            //double base=log10(binsize/60.)+log10(maxrad/binsize*60.)*(bb-0.5)/nbin;
            double base1=log10(binsize/60.)+log10(maxrad/binsize*60.)*(bb-1)/nbin;
            double base2=log10(binsize/60.)+log10(maxrad/binsize*60.)*bb/nbin;
            double be=1./2.*(pow(10.,base1)+pow(10.,base2));
            double ebe=pow(10.,base2)-pow(10.,base1);
            db=ebe;
            binedge=bins[i];
            i++;
            if (i>nbin) {
                break;
            }
        }
        binedge+=binsize/60./2.;
        double thisbin=binedge;
        int b2=1;
        while (thisbin<maxrad) {
            double base1=log10(binedge)+log10(maxrad/binsize*60.)*(b2-1)/nbin;
            double base2=log10(binedge)+log10(maxrad/binsize*60.)*b2/nbin;
            thisbin=1./2.*(pow(10.,base1)+pow(10.,base2));
            double ebe=1./2.*(pow(10.,base2)-pow(10.,base1));
            bins[i]=thisbin;
            ebins[i]=ebe;
            b2++;
            i++;
            if (i>nbin) {
                break;
            }
        }
        newnbin=i-1;
    }while (0);
}


void mk_sb_profile(double *img,double *exposure,double *profile,double *eprof,double *bins,double *ebins,int &nbin,long *axes,
                   double centroid_ra,double centroid_dec,double pixsize,double maxrad,double binsize,bool islogbin,double anglow,double anghigh){
    int nnn[nbin];
    for (int i=0;i<nbin;i++){
        eprof[i]=0.0;
        profile[i]=0.0;
        nnn[i]=0;
    }
    if (islogbin) {
        int newbin=0;
        logbinning(binsize,maxrad,nbin,bins,ebins,newbin);
        double *tb=new double[newbin+1];
        tb[0]=0.0;
        for (int i=1; i<newbin+1; i++) {
            tb[i]=bins[i-1]+ebins[i-1];
        }
        nbin=newbin;
        delete [] tb;
    }
    double angh=anghigh;
    double angl=anglow;
    if (anghigh<anglow){//We cross the zero
        anghigh+=2*pi-anglow;//angle with respect to anglow
    }
    else {
        anghigh-=anglow;
    }
    for (int i=0;i<axes[0];i++){
        for (int j=0;j<axes[1];j++){
            double posx=(i-centroid_ra)*pixsize*60;//arcmin
            double posy=(j-centroid_dec)*pixsize*60;
            double angle=atan(posy/posx);
            //Set all angles between 0 and 2pi
            if (posx<0.0){
                angle+=pi;
            }
            if ((posy<0.0)&&(posx>0.0)){
                angle+=2.*pi;
            }
            if (angh<angl && angle<angl){// we cross the zero
                angle+=2*pi-anglow;
            }
            else angle-=anglow; //Set the origin at anglow
            double dist=sqrt(posx*posx+posy*posy);
            if ((dist<maxrad)&&(exposure[j*axes[0]+i]>0.0)&&(angle>0.0)&&(angle<anghigh)){
                int bin;
                if (!islogbin) {
                    bin=(int)floor(dist/(binsize/60.)); //left-inclusive
                }
                else {
                    bin=findbin(dist,bins);
                }
                if (bin<nbin){
                    profile[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i];
                    eprof[bin]+=img[j*axes[0]+i]/exposure[j*axes[0]+i]/exposure[j*axes[0]+i];
                    nnn[bin]++;
                }
            }
        }
    }
    for (int i=0;i<nbin;i++){
        if (!islogbin) {
            bins[i]=(i+0.5)*binsize/60.;
            ebins[i]=binsize/60./2.;			
        }
        if (nnn[i]>0){
            profile[i]=profile[i]/nnn[i];
            eprof[i]=sqrt(eprof[i])/nnn[i];
            profile[i]/=pixsize*pixsize*60*60;
            eprof[i]/=pixsize*pixsize*60*60;
        }
        else {
            profile[i]=0.0;
            eprof[i]=0.0;
        }
    }
}

void invertp(double *inp,double *outp,int np){
    for (int i=0; i<np; i++) {
        int k=np-i-1;
        outp[k]=inp[i];
    }
}

double vij(double *binsh,int i,int j){
    double fact=0.0;
    if (j>i){
        double cd=binsh[i]*binsh[i]-binsh[j]*binsh[j];
        fact=4./3.*pi*pow(cd,3./2.);
    }
    else {
        fact=0.0;
    }
    return fact;
}

double* medsmooth(int nbin,double *profile,int width){
    if (nbin+1<width) {
        printf("Error: Size of smoothing window must be smaller than array size\n");
        return profile;
    }
    else if (width%2==0){
        printf("Error: Window width must be an odd integer\n");
        return profile;
    }
    else {
        double xx[width];
        double *smoothed=new double[nbin];
        int nm=width/2-1;
        for (int i=0; i<nbin; i++) {
            int tn=0;
            for (int j=i-nm; j<i+nm+1; j++) {
                if (j>-1 && j<nbin) {
                    xx[tn]=profile[j];
                    tn++;
                }
            }
            smoothed[i]=CalcMedian(tn,xx);
        }
        //Fix the ends
        double Y0=3.*profile[0]-2.*profile[1];
        xx[0]=Y0;
        xx[1]=profile[0];
        xx[2]=profile[1];
        smoothed[0]=CalcMedian(3,xx);
        Y0=3.*profile[nbin-1]-2.*profile[nbin-2];
        xx[0]=Y0;
        xx[1]=profile[nbin-1];
        xx[2]=profile[nbin-2];
        smoothed[nbin-1]=CalcMedian(3,xx);
        return smoothed;
    }
}

void deproject(double* profile, int nbin, double *bins, double *ebins, double *deprof){
    double b=1.0;///3.;
    double *binsh=new double[nbin+1];
    double *invpr=new double[nbin];
    double *psm=medsmooth(nbin,profile,5);
    invertp(psm,invpr,nbin);
    binsh[0]=0.1;
    for (int i=0;i<nbin;i++){
        double bh=bins[i]+ebins[i];
        binsh[i+1]=bh; //radius in pixel
    }
    double *radius=new double[nbin+1];
    invertp(binsh,radius,nbin+1);
    
    double *indpr=new double[nbin];
    //edge correction
    double Rm=binsh[nbin];
    double Rm1=binsh[nbin-1];
    double rin=Rm1;
    double rout=Rm;
    double term1=1.0;
    double term2=rout/rin*acos(rin/Rm);
    double term3=rout/Rm*sqrt(1-rin*rin/Rm/Rm);
    double term4=rout/rin-1.;
    double f=term1*(1.-2./pi*(term2-term3)/term4);
    double anb=pi*(radius[0]*radius[0]-radius[1]*radius[1]);
    double Nout=invpr[0]*anb/b/(vij(radius,0,1)-vij(radius,1,1)-vij(radius,0,0)+vij(radius,1,0));
    indpr[0]=Nout*(1.-f);
    for (int m=1; m<nbin; m++) {
        anb=pi*(radius[m]*radius[m]-radius[m+1]*radius[m+1]);
        double sum=0.0;
        for (int i=0; i<m; i++) {
            double vol=vij(radius,i,m+1)-vij(radius,i+1,m+1)-vij(radius,i,m)+vij(radius,i+1,m);
            sum+=indpr[i]*vol;
        }
        rin=radius[m+1];
        rout=radius[m];
        term1=(Rm1+Rm)*Rm*Rm1/((rin+rout)*rin*rout);
        term2=(rout/rin*acos(rin/Rm)-acos(rout/Rm));
        term3=rout/Rm*(sqrt(1-rin*rin/Rm/Rm)-sqrt(1-rout*rout/Rm/Rm));
        term4=rout/rin-1.;
        f=term1*(1.-2./pi*(term2-term3)/term4);
        if (m==nbin-1) {
            f=0.0;
        }
        indpr[m]=(invpr[m]*anb/b-sum)/(vij(radius,m,m+1)-vij(radius,m+1,m+1)-vij(radius,m,m)+vij(radius,m+1,m))-f*Nout;
    }
    double *dumm=medsmooth(nbin,indpr,5);
    invertp(dumm,deprof,nbin);
    
    delete [] dumm;
    delete [] psm;
    delete [] invpr;
    delete [] indpr;
    delete [] radius;
    delete [] binsh;
}

void read_file(char *filename,double *buffer,int npix){
    FILE *fl=fopen(filename,"r");
    string Line;
    stringstream iss; //declare the stream
    int pix=0;
    
    while (pix<npix)  // Put each line into variable "Line"
    {
        fscanf(fl,"%lf",&buffer[pix]);
        pix++;
    }
    fclose(fl);
}

void proffit(int argc, char **argv){
    do {
        int status=0;
        if (argc!=9){
            printf("Usage: \n");
            printf("deprojection imgfile expmap bkgmap cx cy anglow anghigh outfile\n");
            break;
        }
        double *img=NULL;
        bool isimg=false;
        double *exposure=NULL;
        bool isexp=false;
        double *backmap=NULL;
        double *profile=NULL;
        bool isprofile=false;
        double *bins=NULL;
        bool isbins=false;
        double *eprof=NULL;
        double *ebins=NULL;
        double *ecp=NULL;
        double *backprof=NULL;
        double *deprof=NULL;
        double *edeprof=NULL;
        bool islogbin=true;
        double binsize,maxrad,pixsize,centroid_ra,centroid_dec,rad2pix,bin2pix;
        *argv++;
        char *imgfile=*argv;
        *argv++;
        char *expmap=*argv;
        *argv++;
        char *bkgmap=*argv;
        *argv++;
        double cx=atof(*argv);
        *argv++;
        double cy=atof(*argv);
        *argv++;
        double anglow=atof(*argv)*pi/180.;
        *argv++;
        double anghigh=atof(*argv)*pi/180.;
        *argv++;
        char *outfile=*argv;
        centroid_ra=cx-1;
        centroid_dec=cy-1;
        binsize=4*60.;
        pixsize=1./60.;
        int nlin=line_num(imgfile);
        int ncol=col_num(imgfile);
        long npix=nlin*ncol;
        long axes[2];
        axes[0]=ncol;
        axes[1]=nlin;
        maxrad=nlin/2.*1.00;
        img=new double[npix];
        read_file(imgfile,img,npix);
        int nlt,nct;
        nlt=line_num(expmap);
        nct=col_num(expmap);
        if (nlt!=nlin || nct!=ncol) {
            printf("Error: Exposure map does not have the same dimension as image file\n");
            break;
        }
        nlt=line_num(bkgmap);
        nct=col_num(bkgmap);
        if (nlt!=nlin || nct!=ncol) {
            printf("Error: Background map does not have the same dimension as image file\n");
            break;
        }
        exposure=new double[npix];
        backmap=new double[npix];
        read_file(expmap,exposure,npix);
        read_file(bkgmap,backmap,npix);
        printf("Input file read successfuly\n");
        bin2pix=binsize/pixsize/3600.;
        rad2pix=maxrad/pixsize/60.;
        int nbin=(int)floor(rad2pix/bin2pix);
        profile=NULL;
        bins=NULL;
        eprof=NULL;
        ebins=NULL;
        profile=new double[nbin];
        bins=new double[nbin];
        eprof=new double[nbin];
        ebins=new double[nbin];
        mk_sb_profile(img,exposure,profile,eprof,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin,anglow,anghigh);
        nbin=(int)floor(rad2pix/bin2pix);
        backprof=NULL;
        bins=NULL;
        ebins=NULL;
        backprof=new double[nbin];
        bins=new double[nbin];
        ebins=new double[nbin];
        double *ebp=new double[nbin];
        double *eee=new double[nbin];
        mk_sb_profile(backmap,exposure,backprof,ebp,bins,ebins,nbin,axes,centroid_ra,centroid_dec,pixsize,maxrad,binsize,islogbin,anglow,anghigh);
        for (int i=0; i<nbin; i++) {
            profile[i]-=backprof[i];
        }
        delete [] ebp;
        delete [] eee;
        //deproject
        deprof=NULL;
        deprof=new double[nbin];
        edeprof=NULL;
        edeprof=new double[nbin];
        deproject(profile,nbin,bins,ebins,deprof);
        double *binsh=new double[nbin+1];
        binsh[0]=1e-2;
        for (int i=0;i<nbin;i++){
            binsh[i+1]=bins[i]+ebins[i];
        }
        char temp[200];
        /*sprintf(temp,"%s.root",outfile);
        TFile *ffr=new TFile(temp,"recreate");
        TH1F *hhp=new TH1F("profile","profile",nbin,binsh);
        TH1F *hhdep=new TH1F("hdep","hdep",nbin,binsh);*/
        sprintf(temp,"%s.txt",outfile);
        FILE *fde=fopen(temp,"w");
        fprintf(fde,"R [pixel]   SB [counts/s/pixel]  Error  Deprojected\n");
        for (int i=0;i<nbin;i++){
            /*hhp->SetBinContent(i+1,profile[i]);
            hhdep->SetBinContent(i+1,deprof[i]);*/
            fprintf(fde,"%g  %g  %g  %g\n",bins[i],profile[i],eprof[i],deprof[i]);
         }
        fclose(fde);
        /*hhp->Write();
        hhdep->Write();
        ffr->Close();*/
    }
    while (0);
}

int main(int argc, char **argv) {
    proffit(argc,argv);
}

