
/*                                            |--------------|---------------------------------|
                                              |         Name | Tanmay Kirtania                 |
                                              |           ID | 2016-1-60-071                   |
                                              |         Dept | Computer Science & Engineering  |
                                              |    Institute | East West University            |
                                              |--------------|---------------------------------|                            */

#include <bits/stdc++.h>
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include <windows.h>

#define im int
#define jay main
#define bye return
#define O void
#define u cout
#define v cin
#define z endl
#define w setw
#define f(i, n, test) for(int i=n; i<test; ++i)
#define fe(i, n, test) for(i=n; i<=test; ++i)
#define fd(i, n, test) for(i=n; i<test; ++i)
#define fde(i, n, test) for(i=n; i<=test; ++i)

using namespace std;

void CFLinear();
void CFPolynomial();
void CFTranscendental();
void Interpolation();
void Menu();

int n;

im jay()
{
    int ch;
    while(true)
    {
        Menu();
        v>>ch;
        switch(ch)
        {
        case 1:
            CFLinear();
            u<<endl<<" Get back to the menu"<<z;
            system("pause");
            system("cls");
            break;
        case 2:
            CFPolynomial();
            u<<z<<"Get back to the menu"<<z;
            system("pause");
            system("cls");
            break;
        case 3:
            CFTranscendental();
            u<<z<<" Get back to the menu"<<z;
            system("pause");
            system("cls");
            break;
        case 4:
            Interpolation();
            u<<z<<" Get back to the menu"<<z;
            system("pause");
            system("cls");
            break;
        case 5:
            exit(1);
            break;
        default:
            u<<z<<"Invalid Choice! Press 1 to 5"<<z<<z;
            u<<z<<" Get back to the menu"<<z;
            system("pause");
            system("cls");
        }
    }
    bye 0;
}

O Menu()
{
    u<<"\t\t1. Curve Fitting Linear"<<z;
    u<<"\t\t2. Curve Fitting Polynomial"<<z;
    u<<"\t\t3. Curve Fitting Transcendental"<<z;
    u<<"\t\t4. Newton's Interpolation"<<z;
    u<<"\t\t5. Exit"<<z;
    u<<z<<"Enter your choice: ";
}

O CFLinear()
{
    u.precision(4);
    u.setf(ios::fixed);

    u<<z<<"Program for curve fitting Linear"<<z;
    system("pause");
    system("cls");

    u<<"\nEnter the number of data pairs: ";
    v>>n;
    system("cls");
    u<<"Status:: Number of data pairs is entered :"<<n<<z;

    double x[n],y[n],a,b;
    u<<"\nEnter the values of x-axis: ";
    f(i,0,n)
        v>>x[i];
    system("cls");
    u<<"<>Status:: Values of X are entered"<<z;

    u<<"\nEnter the values of y-axis: ";
    f(i,0,n)
        v>>y[i];
    system("cls");
    u<<"<>Status:: Values of Y are entered"<<z;

    double sum_x=0,sum_x2=0,sum_y=0,sum_xy=0;
    f(i,0,n)
    {
        sum_x+=x[i];
        sum_y+=y[i];
        sum_x2+=pow(x[i],2);
        sum_xy=sum_xy+x[i]*y[i];
    }

    b=(n*sum_xy-sum_x*sum_y)/(n*sum_x2-sum_x*sum_x);
    a=(sum_x2*sum_y-sum_x*sum_xy)/(sum_x2*n-sum_x*sum_x);

    double y_fit[n];
    f(i,0,n)
        y_fit[i]=b*x[i]+a;

    u<<"-----------------------------------------------------------------\n";
    u<<"Iterations"<<w(5)<<"x"<<w(19)<<"y(observed)"<<w(19)<<"y(fitted)"<<z;
    u<<"-----------------------------------------------------------------\n";
    f(i,0,n)
        u<<i+1<<"."<<w(17)<<x[i]<<w(15)<<y[i]<<w(18)<<y_fit[i]<<z;
    u<<"\nThe linear fit line is of the form:\n"<<a<<" + "<<b<<"x"<<z;
}

O CFPolynomial()
{
    int i,j,k,N;
    u.precision(4);
    u.setf(ios::fixed);

    u<<z<<"Program for Curve fitting Polynomial"<<z;
    system("pause");
    system("cls");

    u<<"\nEnter the number of data pairs: ";
    v>>N;
    system("cls");
    u<<"Status::Number of input is entered: "<<N<<z;

    double x[N],y[N];
    u<<"Enter the values of X: ";
    fd(i,0,N)
        v>>x[i];
    system("cls");
    u<<"Status::Values of X is entered"<<z;

    u<<"Enter the values of Y: ";
    fd(i,0,N)
        v>>y[i];
    system("cls");
    u<<"Status::Values of Y is entered"<<z;

    u<<"\nEnter the degree of Polynomial to be used for the fit: ";
    v>>n;
    system("cls");
    u<<"Status::enterer highest order of the equation is: "<<n<<z;

    double X[2*n+1];
    fd(i,0,2*n+1)
    {
        X[i]=0;
        fd(j,0,N)
            X[i]=X[i]+pow(x[j],i);
    }

    double B[n+1][n+2],a[n+1];
    fde(i,0,n)
        fde(j,0,n)
            B[i][j]=X[i+j];

    double Y[n+1];
    fd(i,0,n+1)
    {
        Y[i]=0;
        fd(j,0,N)
            Y[i]=Y[i]+pow(x[j],i)*y[j];
    }
    fde(i,0,n)
        B[i][n+1]=Y[i];
    n++;
    system("cls");
    u<<"status::Augmented Matrix of the following linear system:\n";

    fd(i,0,n)
    {
        fde(j,0,n)
            u<<B[i][j]<<w(16);
        u<<"\n";
    }
    fd(i,0,n)
        fd(k,i+1,n)
            if (B[i][i]<B[k][i])
                fde(j,0,n)
                {
                    double temp=B[i][j];
                    B[i][j]=B[k][j];
                    B[k][j]=temp;
                }

    fd(i,0,n-1)
        for (k=i+1;k<n;k++)
            {
                double t=B[k][i]/B[i][i];
                fde(j,0,n)
                    B[k][j]=B[k][j]-t*B[i][j];
            }
    for (i=n-1;i>=0;i--)
    {
        a[i]=B[i][n];
        fd(j,0,n)
            if (j!=i)
                a[i]=a[i]-B[i][j]*a[j];
        a[i]=a[i]/B[i][i];
    }
    u<<"\nThe values of the coefficients are as follows:\n";
    fd(i,0,n)
        u<<"x^"<<i<<"="<<a[i]<<z;

    u<<"\nHence the fitted Polynomial is given by:\ny=";
    fd(i,0,n)
        u<<" + ("<<a[i]<<")"<<"x^"<<i;
    u<<"\n";
    u<<"";
}

O CFTranscendental()
{
    float sum_x=0,sum_y=0,sum_xy=0,sum_x2=0;
    float A,B,a,b;
    u<<z<<"Program for curve fitting Transcendental"<<z;
    system("pause");
    system("cls");

    printf("\n Enter the number of pairs:");
    v>>n;
    system("cls");
    u<<"Status::The number of data pairs entered is: "<<n<<z;

    float x[n],y[n],X[n],Y[n];
    printf("\n Enter the values of x:");
    f(i,0,n)
        v>>x[i];
    system("cls");
    u<<"Status::Values of X are entered"<<z;

    printf("\n Enter the values of y:");
    f(i,0,n)
        v>>y[i];
    system("cls");
    u<<"Status::Values of Y are entered"<<z;

    f(i,0,n)
    {
        X[i]=log1p(x[i]-1);
        Y[i]=log1p(y[i]-1);
    }

    f(i,0,n)
    {
        sum_x+=X[i];
        sum_x2+=X[i]*X[i];
        sum_y+=Y[i];
        sum_xy+=X[i]*Y[i];
    }

    B=((n*sum_xy)-(sum_x*sum_y))/((n*sum_x2)-(sum_x*sum_x));
    A=(sum_y/n)-B*(sum_x/n);
    a=exp(A);

    u<<z<<z<<"----------------------------------------------------------------------------------------------------";
    u<<z<<"\tx"<<w(12)<<"y"<<w(16)<<"ln(x)"<<w(19)<<"ln(y)"<<w(19)<<"(ln(x))^2"<<w(22)<<"ln(x).ln(y)"<<z;
    u<<"----------------------------------------------------------------------------------------------------"<<z;
    for(int i=0;i<n;++i)
        u<<"\t"<<x[i]<<w(12)<<y[i]<<w(16)<<X[i]<<w(19)<<Y[i]<<w(19)<<X[i]*X[i]<<w(22)<<X[i]*Y[i]<<z;
    u<<"----------------------------------------------------------------------------------------------------"<<z;
    u<<"Sum= "<<w(32)<<sum_x<<w(19)<<sum_y<<w(19)<<sum_x2<<w(22)<<sum_xy<<z;

    printf("\n\n The curve is Y= %4.3fx^%4.3f\n",a,B);
}

O Interpolation()
{
    int j;
    float a,sum=0,mul;

    u<<z<<"Program for Newton's Interpolation Method"<<z;
    system("pause");
    system("cls");

    u<<"Enter the number of points: ";
    v>>n;
    float x[n],f[n];

    u<<"\nEnter the values of x: ";
    for(int i=0;i<n;i++)
        v>>x[i];
    u<<"\nEnter the values of y corresponding f(x): ";
    for(int i=0;i<n;i++)
        v>>f[i];

    u<<"\nEnter x for calculation: ";
    v>>a;

    fd(j,0,n-1)
    {
        for(int i=n-1;i>j;--i)
            f[i]=(f[i]-f[i-1])/(x[i]-x[i-j-1]);
    }

    for(int i=n-1;i>=0;--i)
    {
        mul=1;
        fd(j,0,i)
            mul*=(a-x[j]);

        mul*=f[j];
        sum+=mul;
    }

    u<<z<<"The result is, y=f(x)= "<<sum<<z;
}
