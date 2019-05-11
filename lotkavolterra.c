#include<stdio.h> // Equazione Lortka-Volterra RK4
#include<stdlib.h>
#include<math.h>
#include<time.h>
void LotkaVolterraStandard(double a,double b,double c, double d,double dt,double x0,double y0);
void LotkaVolterraSpecial(double a, double b, double c,double d, double e1, double e2,double
dt,double x0,double y0);
void LotkaVolterraPesca(double a,double b,double c, double d,double e,double f,double dt,double
x0,double y0);
void Errore(double a,double b,double c, double d,double dt,double x0,double y0);
void FunzioneH(double a,double b,double c, double d,double dt,double x0,double y0);
void Inizializzazione();
int main(){
Inizializzazione();
}
void LotkaVolterraStandard(double a,double b,double c, double d,double dt,double x0,double y0){
double k1,k2,k3,k4,l1,l2,l3,l4,X,Y,t,H;
int i=0;
double J0[1][1],J1[1][1];
J0[0][0]=a;
J0[0][1]=0;
J0[1][0]=0;
J0[1][1]=-c;
J1[0][0]=0;
J1[0][1]=-(b*c)/d;
J1[1][0]=(a*d)/b;
J1[1][1]=0;
printf("I punti di equilibrio sono:\n1)(x,y)=(0,0)\n2)(x,y)=(%lf,%lf)\n\n",c/d,a/b);
printf("La matrice jacobiana calcolata nel primo punto di equilibrio (0,0) è: \n\n %lf%lf \n",J0[0][0],J0[0][1]);
printf("J(0,0)=\n");
printf(" %lf %lf \n\n",J0[1][0],J0[1][1]);
printf("La matrice jacobiana calcolata nel secondo punto di equilibrio (%.2lf,%.2lf) è: \n\n%lf %lf \n",c/d,a/b,J1[0][0],J1[0][1]);
printf("J(%.2lf,%.2lf)=\n",c/d,a/b);
printf(" %lf %lf \n",J1[1][0],J1[1][1]);
 FILE *fp;
 fp=fopen("LVstd.dat","w+");
fprintf(fp,"#x0 y0 t H\n");
 t=0;
for(i=0;i<=1000000;i++){
 H=-(c*log(x0)-d*x0-b*y0+a*log(y0));
k1=x0*(a-b*y0);
l1=y0*(-c+d*x0);
k2=(x0+0.5*dt*k1)*(a-b*(y0+0.5*dt*l1));
l2=(y0+0.5*dt*l1)*(-c+d*(x0+0.5*dt*k1));
k3=(x0+0.5*dt*k2)*(a-b*(y0+0.5*dt*l2));
l3=(y0+0.5*dt*l2)*(-c+d*(x0+0.5*dt*k2));
k4=(x0+dt*k3)*(a-b*(y0+dt*l3));
l4=(y0+dt*l3)*(-c+d*(x0+dt*k3));
X=x0+(1./6)*dt*(k1+2*k2+2*k3+k4);
Y=y0+(1./6)*dt*(l1+2*l2+2*l3+l4);
fprintf(fp,"%.10lf %.10lf %lf %.20lf\n",x0,y0,t,H);
x0=X;
y0=Y;
t=t+dt;
}
fclose(fp);
}
void LotkaVolterraSpecial(double a, double b, double c,double d, double e1, double e2,double
dt,double x0,double y0){
double m1,m2,m3,m4,n1,n2,n3,n4,X2,Y2,t2,H,X,Y;
int h=0;
double J0[1][1],J1[1][1];
X=(a*e2+b*c)/(b*d+e1*e2);
Y=(a*d-c*e1)/(b*d+e1*e2);
J0[0][0]=a;
J0[0][1]=0;
J0[1][0]=0;
J0[1][1]=-c;
J1[0][0]=-e1*X;
J1[0][1]=-b*X;
J1[1][0]=d*Y;
J1[1][1]=-e2*Y;
printf("I punti di equilibrio sono:\n1)(x,y)=(0,0)\n2)(x,y)=(%lf,%lf)\n3) (x,y)=(%lf,%lf)\n4)(x,y)=(%lf,%lf)\n\n",(a*e2+b*c)/(b*d+e1*e2),(a*d-c*e1)/(b*d+e1*e2),0.,-c/e2,a/e1,0.);
printf("La matrice jacobiana calcolata nel primo punto di equilibrio (0,0) è: \n\n%lf %lf \n",J0[0][0],J0[0][1]);
printf("J(0,0)=\n");
printf(" %lf %lf \n\n",J0[1][0],J0[1][1]);
printf("La matrice jacobiana calcolata nel secondo punto di equilibrio (%lf,%lf) è: \n\n %lf %lf \n",X,Y,J1[0][0],J1[0][1]);
printf("J(%lf,%lf)=\n",X,Y);
printf(" %lf %lf \n",J1[1][0],J1[1][1]);
 FILE *fq;
 fq=fopen("LVsp.dat","w+");
 t2=0;
fprintf(fq,"#x0 y0 t H\n");
for(h=0;h<=100000;h++){
H=c*log(x0)-d*x0-b*y0+a*log(y0);
m1=x0*(a-b*y0-e1*x0);
n1=y0*(-c+d*x0-e2*y0);
m2=(x0+0.5*dt*m1)*(a-b*(y0+0.5*dt*n1)-e1*(x0+0.5*dt*m1));
n2=(y0+0.5*dt*n1)*(-c+d*(x0+0.5*dt*m1)-e2*(y0+0.5*dt*n1));
//printf("%lf %lf\n",m2,n2);
m3=(x0+0.5*dt*m2)*(a-b*(y0+0.5*dt*n2)-e1*(x0+0.5*dt*m2));
n3=(y0+0.5*dt*n2)*(-c+d*(x0+0.5*dt*m2)-e2*(y0+0.5*dt*n2));
m4=(x0+dt*m3)*(a-b*(y0+dt*n3)-e1*(x0+dt*m3));
n4=(y0+dt*n3)*(-c+d*(x0+dt*m3)-e2*(y0+dt*n3));
X2=x0+(1./6)*dt*(m1+2*m2+2*m3+m4);
Y2=y0+(1./6)*dt*(n1+2*n2+2*n3+n4);
fprintf(fq,"%lf %lf %lf %lf\n",x0,y0,t2,H);
x0=X2;
y0=Y2;
t2=t2+dt;
}
fclose(fq);
}
void LotkaVolterraPesca(double a,double b,double c, double d,double e,double f,double dt,double
x0,double y0){
double m1,m2,m3,m4,n1,n2,n3,n4,X2,Y2,t2,H,X,Y;
int h=0;
double J0[1][1],J1[1][1];
J0[0][0]=a-e;
J0[0][1]=0;
J0[1][0]=0;
J0[1][1]=-c-f;
J1[0][0]=0;
J1[0][1]=(-b*(c+f))/d;
J1[1][0]=(d*(a-e))/b;
J1[1][1]=0;
printf("I punti di equilibrio sono:\n1)(x,y)=(0,0)\n2)(x,y)=(%lf,%lf)\n\n",(c+f)/d,(a-e)/b);
printf("La matrice jacobiana calcolata nel primo punto di equilibrio (0,0) è: \n\n%lf %lf \n",J0[0][0],J0[0][1]);
printf("J(0,0)=\n");
printf(" %lf %lf \n\n",J0[1][0],J0[1][1]);
printf("La matrice jacobiana calcolata nel secondo punto di equilibrio (%lf,%lf) è: \n\n%lf %lf \n",(c+f)/d,(a-e)/b,J1[0][0],J1[0][1]);
printf("J(%lf,%lf)=\n",(c+f)/d,(a-e)/b);
printf(" %lf %lf \n",J1[1][0],J1[1][1]);
 FILE *fu;
 fu=fopen("LVfish.dat","w+");
 t2=0;
fprintf(fu,"#x0 y0 t H\n");
for(h=0;h<=100000;h++){
H=-((c+f)*log(x0)-d*x0-b*y0+(a-e)*log(y0));
m1=x0*(a-b*y0-e);
n1=y0*(-c+d*x0-f);
m2=(x0+0.5*dt*m1)*(a-b*(y0+0.5*dt*n1)-e);
n2=(y0+0.5*dt*n1)*(-c+d*(x0+0.5*dt*m1)-f);
//printf("%lf %lf\n",m2,n2);
m3=(x0+0.5*dt*m2)*(a-b*(y0+0.5*dt*n2)-e);
n3=(y0+0.5*dt*n2)*(-c+d*(x0+0.5*dt*m2)-f);
m4=(x0+dt*m3)*(a-b*(y0+dt*n3)-e);
n4=(y0+dt*n3)*(-c+d*(x0+dt*m3)-f);
X2=x0+(1./6)*dt*(m1+2*m2+2*m3+m4);
Y2=y0+(1./6)*dt*(n1+2*n2+2*n3+n4);
fprintf(fu,"%lf %lf %lf %lf\n",x0,y0,t2,H);
x0=X2;
y0=Y2;
t2=t2+dt;
}
fclose(fu);
}
void Errore(double a,double b,double c, double d,double dt,double x0,double y0){
FILE *fs;
fs=fopen("LVerr.dat","w+");
 fprintf(fs,"#t dt H H0|(H-H0)/H0|\n");
double k1,k2,k3,k4,l1,l2,l3,l4,X,Y,t,H,H0,xi,yi;
H0=-(c*log(x0)-d*x0-b*y0+a*log(y0));
//printf("%lf\n",log(2.7182818284));
xi=x0;
yi=y0;
do{
x0=xi;
y0=yi;
t=0;
while(t<100){
//printf("%lf\n",t);
H=-(c*log(x0)-d*x0-b*y0+a*log(y0));
k1=x0*(a-b*y0);
l1=y0*(-c+d*x0);
k2=(x0+0.5*dt*k1)*(a-b*(y0+0.5*dt*l1));
l2=(y0+0.5*dt*l1)*(-c+d*(x0+0.5*dt*k1));
k3=(x0+0.5*dt*k2)*(a-b*(y0+0.5*dt*l2));
l3=(y0+0.5*dt*l2)*(-c+d*(x0+0.5*dt*k2));
k4=(x0+dt*k3)*(a-b*(y0+dt*l3));
l4=(y0+dt*l3)*(-c+d*(x0+dt*k3));
X=x0+(1./6)*dt*(k1+2*k2+2*k3+k4);
Y=y0+(1./6)*dt*(l1+2*l2+2*l3+l4);
x0=X;
y0=Y;
 t=t+dt;
if(t>2-0.5*dt && t<2+0.5*dt){
//printf(" %lf %.14lf %lf %lf%.20lf\n",t,dt,H0,H0,fabs(((H-H0)/H0)));
 fprintf(fs," %lf %.14lf %.10lf %.10lf%.20lf\n",t,dt,H,H0,fabs(((H-H0)/H0)));
}
}
dt=dt/2;
}while(dt>0.0001);
fclose(fs);
}
void FunzioneH(double a,double b,double c, double d,double dt,double x0,double y0){
double k1,k2,k3,k4,l1,l2,l3,l4,X,Y,t,H,H0,k,l;
double tmp1,tmp2;
tmp1=x0;
tmp2=y0;
int i=0;
 FILE *ft;
 ft=fopen("LVfun.dat","w+");
fprintf(ft,"#x0 y0 H0\n");
for(k=1;k<=tmp1;k++){
printf("%lf\n",k);
x0=k;
for(l=1;l<=tmp2;l++){
 H0=-(c*log(x0)-d*x0-b*l+a*log(l));
 t=0;
 y0=l;
for(i=0;i<10000;i++){
 H=c*log(x0)-d*x0-b*y0+a*log(y0);
k1=x0*(a-b*y0);
l1=y0*(-c+d*x0);
k2=(x0+0.5*dt*k1)*(a-b*(y0+0.5*dt*l1));
l2=(y0+0.5*dt*l1)*(-c+d*(x0+0.5*dt*k1));
k3=(x0+0.5*dt*k2)*(a-b*(y0+0.5*dt*l2));
l3=(y0+0.5*dt*l2)*(-c+d*(x0+0.5*dt*k2));
k4=(x0+dt*k3)*(a-b*(y0+dt*l3));
l4=(y0+dt*l3)*(-c+d*(x0+dt*k3));
X=x0+(1./6)*dt*(k1+2*k2+2*k3+k4);
Y=y0+(1./6)*dt*(l1+2*l2+2*l3+l4);
fprintf(ft,"%lf %lf %lf\n",x0,y0,H0);
x0=X;
y0=Y;
t=t+dt;
}
}
}
fclose(ft);
}
void Inizializzazione(){
double dt,x0,y0,e1,e2;
double a,b,c,d,f1,f2,e,f;
int q,w,i;
double seed;
seed=time(0);
srand(seed);
printf("\nQuale modello vuoi studiare?\n\n1)Modello classico...........................(premere1)\n2)Modello con competizione intraspecie.......(premere 2)\n3)Modello con superpredatore................(premere 3)\n\n");
scanf("%d",&q);
printf("\n");
printf("In che modo vuoi procedere?\n\n1)Inserimento numeri da default.......(premere1)\n2)Inseriemento numeri da tastiera.....(premere 2)\n3)Inserimento numerirandom...........(premere 3)\n\n");
scanf("%d",&w);
printf("\n");
if(w==1){
dt=0.001;
x0=8;
 y0=5;
if(q==1){
 a=4;
b=2;
c=3;
d=1;
LotkaVolterraStandard(a,b,c,d,dt,x0,y0);
 FunzioneH(a,b,c,d,0.001,x0,y0);
}
if(q==2){
a=6;
b=3;
c=4;
d=1;
f1=0.01;
f2=0.01;
LotkaVolterraSpecial(a,b,c,d,f1,f2,dt,x0,y0);
}
if(q==3){
a=4;
b=2;
c=3;
d=1;
e=2;
f=1;
LotkaVolterraPesca(a,b,c,d,e,f,dt,x0,y0);
FunzioneH(a-e,b,c+f,d,0.001,x0,y0);
}
Errore(a,b,c,d,0.03,x0,y0);
}
if(w==2){
 printf("dimmi il passo di integrazione (dt)\n");
 scanf("%lf",&dt);
 printf("dimmi il numero iniziale di prede (x0)\n");
 scanf("%lf",&x0);
 printf("dimmi il numero iniziale di predatori (y0)\n");
 scanf("%lf",&y0);
printf("dimmi a\n");
 scanf("%lf",&a);
 printf("dimmi b\n");
 scanf("%lf",&b);
 printf("dimmi c\n");
 scanf("%lf",&c);
 printf("dimmi d\n");
 scanf("%lf",&d);
if(q==1){
LotkaVolterraStandard(a,b,c,d,dt,x0,y0);
FunzioneH(a,b,c,d,dt,x0,y0);
}
if(q==2){
 printf("dimmi f1\n");
scanf("%lf",&f1);
 printf("dimmi f2\n");
 scanf("%lf",&f2);
 LotkaVolterraSpecial(a,b,c,d,f1,f2,dt,x0,y0);
}
if(q==3){
printf("dimmi e\n");
scanf("%lf",&e);
printf("dimmi f\n");
scanf("%lf",&f);
LotkaVolterraPesca(a,b,c,d,e,f,dt,x0,y0);
FunzioneH((a-e),b,(c+f),d,0.001,x0,y0);
}
Errore(a,b,c,d,0.03,x0,y0);
}
if(w==3){
a=20*rand();
b=20*rand();
c=20*rand();
d=20*rand();
x0=(int)(20*rand());
y0=(int)(20*rand());
dt=0.001;
printf("I numeri estratti sono:\n\na=%lf b=%lf c=%lf d=%lf x0=%lfy0=%lf\n\n",a,b,c,d,x0,y0);
if(q==1){
LotkaVolterraStandard(a,b,c,d,dt,x0,y0);
FunzioneH(a,b,c,d,dt,x0,y0);
}
if(q==2){
 f1=20*rand();
 f2=20*rand();
 printf("f1=%lf f2=%lf \n\n",f1,f2);
 LotkaVolterraSpecial(a,b,c,d,f1,f2,dt,x0,y0);
}
if(q==3){
e=20*rand();
f=20*rand();
LotkaVolterraPesca(a,b,c,d,e,f,dt,x0,y0);
FunzioneH(a-e,b,c+f,d,0.001,x0,y0);
}
Errore(a,b,c,d,0.03,x0,y0);
}
}
//FINE Equazione Lotka-Volterra RK4
