#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#include "nrcomplex.h"

fcomplexNR Cadd(fcomplexNR a, fcomplexNR b)
{
    fcomplexNR c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

fcomplexNR Csub(fcomplexNR a, fcomplexNR b)
{
    fcomplexNR c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}


fcomplexNR Cmul(fcomplexNR a, fcomplexNR b)
{
    fcomplexNR c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

fcomplexNR ComplexNR(float re, float im)
{
    fcomplexNR c;
    c.r=re;
    c.i=im;
    return c;
}

fcomplexNR Conjg(fcomplexNR z)
{
    fcomplexNR c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

fcomplexNR Cdiv(fcomplexNR a, fcomplexNR b)
{
    fcomplexNR c;
    float r,den;
    if (fabs(b.r) >= fabs(b.i)) {
        r=b.i/b.r;
        den=b.r+r*b.i;
        c.r=(a.r+r*a.i)/den;
        c.i=(a.i-r*a.r)/den;
    } else {
        r=b.r/b.i;
        den=b.i+r*b.r;
        c.r=(a.r*r+a.i)/den;
        c.i=(a.i*r-a.r)/den;
    }
    return c;
}

float Cabs(fcomplexNR z)
{
    float x,y,ans,temp;
    x=fabs(z.r);
    y=fabs(z.i);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    } else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

fcomplexNR Csqrt(fcomplexNR z)
{
    fcomplexNR c;
    float x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    } else {
        x=fabs(z.r);
        y=fabs(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        } else {
            r=x/y;
            w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        } else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

fcomplexNR RCmul(float x, fcomplexNR a)
{
    fcomplexNR c;
    c.r=x*a.r;
    c.i=x*a.i;
    return c;
}

#else /* ANSI */
/* traditional - K&R */

fcomplexNR Cadd(a,b)
fcomplexNR a,b;
{
    fcomplexNR c;
    c.r=a.r+b.r;
    c.i=a.i+b.i;
    return c;
}

fcomplexNR Csub(a,b)
fcomplexNR a,b;
{
    fcomplexNR c;
    c.r=a.r-b.r;
    c.i=a.i-b.i;
    return c;
}


fcomplexNR Cmul(a,b)
fcomplexNR a,b;
{
    fcomplexNR c;
    c.r=a.r*b.r-a.i*b.i;
    c.i=a.i*b.r+a.r*b.i;
    return c;
}

fcomplexNR ComplexNR(re,im)
float im,re;
{
    fcomplexNR c;
    c.r=re;
    c.i=im;
    return c;
}

fcomplexNR Conjg(z)
fcomplexNR z;
{
    fcomplexNR c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

fcomplexNR Cdiv(a,b)
fcomplexNR a,b;
{
    fcomplexNR c;
    float r,den;
    if (fabs(b.r) >= fabs(b.i)) {
        r=b.i/b.r;
        den=b.r+r*b.i;
        c.r=(a.r+r*a.i)/den;
        c.i=(a.i-r*a.r)/den;
    } else {
        r=b.r/b.i;
        den=b.i+r*b.r;
        c.r=(a.r*r+a.i)/den;
        c.i=(a.i*r-a.r)/den;
    }
    return c;
}

float Cabs(z)
fcomplexNR z;
{
    float x,y,ans,temp;
    x=fabs(z.r);
    y=fabs(z.i);
    if (x == 0.0)
        ans=y;
    else if (y == 0.0)
        ans=x;
    else if (x > y) {
        temp=y/x;
        ans=x*sqrt(1.0+temp*temp);
    } else {
        temp=x/y;
        ans=y*sqrt(1.0+temp*temp);
    }
    return ans;
}

fcomplexNR Csqrt(z)
fcomplexNR z;
{
    fcomplexNR c;
    float x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    } else {
        x=fabs(z.r);
        y=fabs(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        } else {
            r=x/y;
            w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        } else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

fcomplexNR RCmul(x,a)
fcomplexNR a;
float x;
{
    fcomplexNR c;
    c.r=x*a.r;
    c.i=x*a.i;
    return c;
}

#endif /* ANSI */
