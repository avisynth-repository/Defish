/*
  "Defish", barrel and pincushion distortion correction filter for AviSynth.

  Copyright (c) 2010 David Horman
  Multithreaded version copyright (c) 2011 Andrey "Efenstor" Pivovarov

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  The original author can be contacted at: david(dot)horman(at)jerseymail(dot)co(dot)uk
  The multithreaded version author can be contacted at: efenstor(at)efenstor(dot)net

  Discussion at http://forum.doom9.org/showthread.php?t=152860
*/

#include <math.h>
#include <stdio.h>
#include <windows.h>
#include "avs_headers\avisynth.h"

#define MIN(X,Y) ((X) > (Y) ? (Y) : (X))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MINMAX(X,Y,Z) (MAX(X,MIN(Y,Z)))

class defish : public GenericVideoFilter {

public:
  defish(PClip _child, double _fov, double _scale, double _aspect, const char *_direction, const char *_scaling, double _a, double _b, double _c, bool _pin, int _threads, IScriptEnvironment* env);
  ~defish();

  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);

  static DWORD WINAPI ThreadProc(__in LPVOID params);

private:
  double *xmap;
  double *ymap;
  int threads;
  struct procParams {
    unsigned int *dstp;
    unsigned int *srcp;
    int width;
    int height;
    int dpitch;
    int spitch;
	int yFrom;
	int yTo;
	defish* pcl;
  };
};

defish::defish(PClip _child, double _fov, double _scale, double _aspect, const char *_direction, const char *_scaling, double _a, double _b, double _c, bool _pin, int _threads, IScriptEnvironment* env) : GenericVideoFilter(_child) {
  if (!vi.IsRGB32()) env->ThrowError("defish: input must be RGB32.");
  if (vi.width&1 || vi.height&1) env->ThrowError("defish: input must have even dimensions.");
  if (_pin && (_a+_b+_c>1)) env->ThrowError("When 'pin' is true, a+b+c must be less than 1.0.");
  _scale=MAX(0,_scale);

  int direction=3;
  if (strcmp(_direction,"x")==0) direction=1;
  else if (strcmp(_direction,"y")==0) direction=2;

  xmap=(double *)calloc((vi.width>>1)*(vi.height>>1),sizeof(double));
  ymap=(double *)calloc((vi.width>>1)*(vi.height>>1),sizeof(double));

  double fov2=(_fov/180.0)*3.14159265358979323846; // radians
  double z=(vi.width/2-0.5)/tan(fov2/2);

  double r_x_to_theta=(fov2/2)/(vi.width/2-0.5); // *scale

  double xp,yp,d,f;

  if (_fov!=0) {
    double fity;
    double autoscale=vi.width/(fov2*z);

    if (direction==1) {
      fity=1;
    } else {
      double d=(vi.height/2-0.5)*_aspect;
      if (_fov>0) fity=((atan(d/z)/(r_x_to_theta))/d);
      else if (_fov<0) fity=((atan(d/z)/(r_x_to_theta))/d);
      else fity=1.0;

      fity=fity/_aspect;
    }

    if (strcmp(_scaling,"fity")==0) {
      autoscale=fity;
    } else if (strcmp(_scaling,"fitx")==0) {
      autoscale=1;
    } else if (strcmp(_scaling,"fitxy")==0) {
      if (_fov>=0) autoscale=MAX(1,fity); else autoscale=MIN(1,fity);
    }

    int p=0;
    for (int y=0; y<(vi.height>>1); y++) {
      yp=(y+0.5)*_aspect/_scale;
      double yp2=yp*yp;
      for (int x=0; x<(vi.width>>1); x++) {
        xp=(x+0.5)/_scale;
        double xp2=xp*xp;
        d=sqrt(xp2+yp2);
        if (_fov>0) f=((atan(d/z)/(r_x_to_theta*autoscale))/d);
        else if (_fov<0) f=1.0/((atan(d/z)/(r_x_to_theta*autoscale))/d);
        else f=1.0;

        if (direction&1) xmap[p]=xp*f; else xmap[p]=xp;
        if (direction&2) ymap[p]=yp*f/_aspect; else ymap[p]=yp;
        p++;
      }
    }
  } else {
    double _d=1-(_a+_b+_c);
    double norm_dx=1.0/((vi.width>>1)-1.5);
    double norm_dy=1.0/(((vi.height>>1)-1.5)*_aspect);
    double norm_d=MAX(norm_dx,norm_dy);

    int p=0;
    for (int y=0; y<(vi.height>>1); y++) {
      yp=(y+0.5)*_aspect/_scale;
      double yp2=yp*yp;
      for (int x=0; x<(vi.width>>1); x++) {
        xp=(x+0.5)/_scale;
        double xp2=xp*xp;
        d=sqrt(xp2+yp2)*norm_d;
        double new_d=_a*d*d*d*d+_b*d*d*d+_c*d*d+_d*d;
        double f;
        if (_pin) f=d/new_d; else f=new_d/d;
        if (direction&1) xmap[p]=xp*f; else xmap[p]=xp;
        if (direction&2) ymap[p]=yp*f/_aspect; else ymap[p]=yp;
        p++;
      }
    }
  }

  if(!_threads)
  {
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	threads = si.dwNumberOfProcessors;
  } else threads = _threads;
}

defish::~defish() {
  free(xmap);
  free(ymap);
}

PVideoFrame __stdcall defish::GetFrame(int n, IScriptEnvironment* env) {
  int i;

  PVideoFrame src = child->GetFrame(n, env);
  PVideoFrame dst = env->NewVideoFrame(vi);

  procParams pb;
  pb.dstp=(unsigned int *)dst->GetWritePtr();
  pb.srcp=(unsigned int *)src->GetReadPtr();
  pb.width=dst->GetRowSize()/4;
  pb.dpitch=dst->GetPitch()/4;
  pb.spitch=src->GetPitch()/4;
  pb.height=dst->GetHeight();
  pb.pcl=this;
  int blockHeight = pb.height/2/threads;

  HANDLE *h = (HANDLE*)malloc(sizeof(HANDLE)*threads);
  procParams *p = (procParams*)malloc(sizeof(procParams)*threads);
  for(i=0; i<threads; i++) {
	memcpy(&p[i], &pb, sizeof(procParams));
	p[i].yFrom = blockHeight*i;
	if(i<threads-1) p[i].yTo = blockHeight*i+blockHeight-1;
	else p[i].yTo = (pb.height/2)-1;
    h[i] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)ThreadProc, &p[i], 0, NULL);
  }

  for(i=0; i<threads; i++) {
    WaitForSingleObject(h[i], INFINITE);
  }

  free(h);
  free(p);
  return dst;
}

DWORD WINAPI defish::ThreadProc(__in LPVOID param)
{
  double xweights[4];

  procParams *pp = (procParams*)param;
  defish *pcl = pp->pcl;
  unsigned int *dstp = pp->dstp;
  unsigned int *srcp = pp->srcp;
  int width = pp->width;
  int height = pp->height;
  int dpitch = pp->dpitch;
  int spitch = pp->spitch;
  int yFrom = pp->yFrom;
  int yTo = pp->yTo;

  int midpoint=(width>>1)+(height>>1)*dpitch;
  int u_offset=midpoint+dpitch*yFrom;
  int b_offset=midpoint-dpitch*yFrom-dpitch;

  int p=yFrom*(width>>1);
  for (int y=yFrom; y<=yTo; y++) {
    double oyc=y+0.5;
    for (int x=0; x<width>>1; x++) {
      int px,py;
      double oxc=x+0.5;

      double bxc=pcl->xmap[p];
      double byc=pcl->ymap[p];
      p++;

      px=(int)(0x1000000+bxc-1.5)-0x1000000;
      py=(int)(0x1000000+byc-1.5)-0x1000000;

      int qx=MIN(px+4,width>>1);
      int qy=MIN(py+4,height>>1);
      int pqx=qx-px;
      int pqy=qy-py;

      if (pqx>0 && pqy>0) {
        int tr_base=midpoint+py*spitch+px;
        int tl_base=midpoint+py*spitch-1-px;
        int br_base=midpoint-(1+py)*spitch+px;
        int bl_base=midpoint-(1+py)*spitch-1-px;
        double dx=px+0.5-bxc;
        double dy=py+0.5-byc;
        double yweight;

        for (int i=0; i<4; i++) {
          double d=fabs(dx);
          if (d<1) xweights[i]=((d-1.8)*d-0.2)*d+1.0;
          else if (d<2) xweights[i]=((-1.0/3.0*(d-1)+0.8)*(d-1)-7.0/15.0)*(d-1);
          else xweights[i]=0;
          dx+=1;
        }

        double tr_r=0,tr_g=0,tr_b=0;
        double tl_r=0,tl_g=0,tl_b=0;
        double br_r=0,br_g=0,br_b=0;
        double bl_r=0,bl_g=0,bl_b=0;
        for (int iy=0; iy<pqy; iy++) {
          double d=fabs(dy);
          dy+=1;
          if (d<1) yweight=((d-1.8)*d-0.2)*d+1.0;
          else if (d<2) yweight=((-1.0/3.0*(d-1)+0.8)*(d-1)-7.0/15.0)*(d-1);
          else continue;
          for (int ix=0; ix<pqx; ix++) {
            double weight=xweights[ix]*yweight;

            int pixel=srcp[tr_base+ix];
            int r=pixel&0xff;
            int g=(pixel>>8)&0xff;
            int b=(pixel>>16)&0xff;
            tr_r+=r*weight;
            tr_g+=g*weight;
            tr_b+=b*weight;

            pixel=srcp[tl_base-ix];
            r=pixel&0xff;
            g=(pixel>>8)&0xff;
            b=(pixel>>16)&0xff;
            tl_r+=r*weight;
            tl_g+=g*weight;
            tl_b+=b*weight;

            pixel=srcp[br_base+ix];
            r=pixel&0xff;
            g=(pixel>>8)&0xff;
            b=(pixel>>16)&0xff;
            br_r+=r*weight;
            br_g+=g*weight;
            br_b+=b*weight;

            pixel=srcp[bl_base-ix];
            r=pixel&0xff;
            g=(pixel>>8)&0xff;
            b=(pixel>>16)&0xff;
            bl_r+=r*weight;
            bl_g+=g*weight;
            bl_b+=b*weight;
          } // ix
          tr_base+=spitch;
          tl_base+=spitch;
          br_base-=spitch;
          bl_base-=spitch;
        } // iy

        int r=(int)MINMAX(0,tr_r+0.5,255);
        int g=(int)MINMAX(0,tr_g+0.5,255);
        int b=(int)MINMAX(0,tr_b+0.5,255);
        dstp[u_offset+x]=(b<<16)|(g<<8)|r;

        r=(int)MINMAX(0,tl_r+0.5,255);
        g=(int)MINMAX(0,tl_g+0.5,255);
        b=(int)MINMAX(0,tl_b+0.5,255);
        dstp[u_offset-1-x]=(b<<16)|(g<<8)|r;

        r=(int)MINMAX(0,br_r+0.5,255);
        g=(int)MINMAX(0,br_g+0.5,255);
        b=(int)MINMAX(0,br_b+0.5,255);
        dstp[b_offset+x]=(b<<16)|(g<<8)|r;

        r=(int)MINMAX(0,bl_r+0.5,255);
        g=(int)MINMAX(0,bl_g+0.5,255);
        b=(int)MINMAX(0,bl_b+0.5,255);
        dstp[b_offset-1-x]=(b<<16)|(g<<8)|r;
      } else {
        dstp[u_offset+x]=0;
        dstp[u_offset-1-x]=0;
        dstp[b_offset+x]=0;
        dstp[b_offset-1-x]=0;
      }
    } // x
    u_offset+=dpitch;
    b_offset-=dpitch;
  } // y

  return 0;
}

AVSValue __cdecl Create_defish(AVSValue args, void* user_data, IScriptEnvironment* env) {
  return new defish(
    args[0].AsClip(),
    args[1].AsFloat(0),
    args[2].AsFloat(1),
    args[3].AsFloat(1),
    args[4].AsString("xy"),
    args[5].AsString(""),
    args[6].AsFloat(0),
    args[7].AsFloat(0),
    args[8].AsFloat(0),
    args[9].AsBool(false),
    args[10].AsInt(0),
    env
  );
}

const AVS_Linkage *AVS_linkage = 0;
extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
  AVS_linkage = vectors;
  env->AddFunction("defish", "c[fov]f[scale]f[aspect]f[direction]s[scaling]s[a]f[b]f[c]f[pin]b[threads]i", Create_defish, 0);
  return "'Defish'";
}

