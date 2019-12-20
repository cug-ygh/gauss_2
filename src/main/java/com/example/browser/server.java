package com.example.browser;

import static java.lang.Math.*;

public class server {
    //一些基本的椭球参数
    //分别是海福特椭球，北京54坐标系基准椭球，西安80坐标系基准椭球，WGS-84椭球
    //长半轴
    public static double[]a={6378388,6378245,6378140,6378137};
    //短半轴
    public static double[]b={6356911.9461279,6356863.018773,6356755.2881575,6356752.3142451};
    //第一偏心率
    private static double e;
    //第二偏心率
    private static double e2;
    //中央子午线的经度
    // public double L0=;


    //大地坐标求解高斯坐标，正算
    //L代表经度，B代表纬度,kind代表椭球类型(0代表海福特椭球，1代表北京54坐标系基准椭球，2代表西安80坐标系基准椭球，
    // 3代表WGS-84椭球),projN代表投影带度
    public static  double[] GeotoGauss(double longitude , double latitude, int kind,int projN){

        //第一偏心率
        e=sqrt(a[kind]*a[kind]+b[kind]*b[kind])/a[kind];
        //第二偏心率
        e2=sqrt(a[kind]*a[kind]+b[kind]*b[kind])/b[kind];
        //角度变弧度
        double L=longitude*PI/180;
        double B=latitude*PI/180;

        double t=tan(B);

        double eta=sqrt((e2*e2)*(cos(B)*cos(B)));

        double p=180/PI*3600;

        //基本常量计算
        double m0=a[kind]*(1-pow(e,2));
        double m2=3/2*pow(e,2)*m0;
        double m4=5/4*pow(e,2)*m2;
        double m6=7/6*pow(e,2)*m4;
        double m8=8/7*pow(e,2)*m6;

        double a0=m0+m2/2+3/8*m4+5/16*m6+35/128*m8;
        double a2=m2/2+m4/2+15/32*m6+7/16*m8;
        double a4=m4/8+3/16*m6+7/32*m8;
        double a6=m6/32+m8/16;
        double a8=m8/128;

        //子午线弧长
        double X=a0*B-sin(B)*cos(B)*(a2-a4+a6+(2*a4-16/3*a6)*pow(sin(B),2)+16/3*a6*pow(sin(B),4));
        //子午圈曲率半径
        double N=a[kind]*pow((1-pow(e,2)*pow(sin(B),2)),-1/2);
        //中央子午线经度
        double L0 = calculator(L,projN);
        double l=L-L0;
        double x=X+N/(2*pow(p,2))*sin(B)*cos(B)*pow(l,2)+N/(24*pow(p,4))*sin(B)*pow(cos(B),3)*(5-pow(t,2)+9*pow(eta,2)+4*pow(eta,4))
                *pow(l,4)+N/(720*pow(p,6))*sin(B)*pow(cos(B),5)*(61-58*pow(t,2)+pow(t,4))*pow(l,6);
        double y=N/p*cos(B)*l+N/(6*pow(p,3))*pow(cos(B),3)*(1-pow(t,2)+pow(eta,2))*pow(l,3)+N/(120*pow(p,5))*pow(cos(B),5)*(5-18*t*t
                +pow(t,4)+14*eta*eta-58*pow(eta,2)*pow(t,2))*pow(l,5);
        //System.out.println(x+","+y);
        return new double[]{x,y};
    }
    public static double calculator(double l,int projn){
        double res=0;
        if(projn==3){
            long n=round(l/3);
            res=n*3;
        }
        if(projn==6){
            //带号
            long n=round((l+3)/6);

            res=6*n-3;
        }


        return res;
    }
    public static double[] Gauss2Geo(double x, double y, int kind, int n, int num){

        //第一偏心率
        e=sqrt(a[kind]*a[kind]+b[kind]*b[kind])/a[kind];
        //第二偏心率
        e2=sqrt(a[kind]*a[kind]+b[kind]*b[kind])/b[kind];

        //基本常量计算
        double m0=a[kind]*(1-pow(e,2));
        double m2=3/2*pow(e,2)*m0;
        double m4=5/4*pow(e,2)*m2;
        double m6=7/6*pow(e,2)*m4;
        double m8=8/7*pow(e,2)*m6;

        double a0=m0+m2/2+3/8*m4+5/16*m6+35/128*m8;
        double a2=m2/2+m4/2+15/32*m6+7/16*m8;
        double a4=m4/8+3/16*m6+7/32*m8;
        double a6=m6/32+m8/16;
        double a8=m8/128;

        double bf=calculatorbf(x,a0,a2,a4,a6);

        double tf=tan(bf);
        double wf=sqrt(1-e*e*sin(bf)*sin(bf));
        double nf=a[kind]/wf;
        double mf=a[kind]*(1-pow(e,2))/pow(wf,3);
        double itaf=e2*cos(bf);
        //计算纬度
        double B=bf-tf*y*y/(2*mf*nf)+tf*(5+3*tf*tf+pow(itaf,2)-9*pow(itaf,2)*tf*tf)*pow(y,4)/(24*mf
        *pow(nf,3))-tf*(61+90*tf*tf+45*pow(tf,4))*pow(y,6)/(720*mf*pow(nf,5));
        //计算经度
        double l=y/(nf*cos(bf))-(1+2*tf*tf+itaf*itaf)*pow(y,3)/(6*pow(nf,3)*cos(bf))+(5+28*tf*tf+24*pow(tf,4)
        +6*itaf*itaf+8*itaf*itaf*tf*tf)*pow(y,5)/(120*pow(nf,5)*cos(bf));
        double L0=0;
        //6度带和三度带
        if(n==3){
            L0=num*3;
        }
        else if(n==6){
            L0=num*6-3;
        }
        double L=l+L0;
        //System.out.println(B+":"+L);
        return new double[]{L,B};
    }

    public static double calculatorbf(double x,double a0,double a1,double a2,double a3) {
        double bf=x/a0;
        double bf1=bf;
        double bf2=0;
        double bf3=0;
        boolean stop=false;
        do{
            stop=false;
            bf3=0.0-a1*sin(2*bf1)/2.0+a2*sin(4*bf1)/4.0+a3*sin(6*bf1)/6.0;
            bf2=(x-bf3)/a0;
            if(abs(bf2-bf1)>(PI*pow(10.0,-8)/(36*18))){
                stop=true;
                bf1=bf2;
            }

        }while(stop==true);
        bf=bf1;
        return bf;
    }

    public static void main(String []args){

//        double []a=GeotoGauss(132.51,30.2,1,3);
//        System.out.println(a[0]);
        Gauss2Geo(2,2,0,3,2);
    }
}
