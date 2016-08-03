#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#define ARRAYD(p) ((double *) (((PyArrayObject *)p)->data))


int entero(double x){
	double y=x+1.0;
	int i;
	i=(int)x;
	int s;
	s=(int)y;
	double ss,ii;
	ss=(double)s;
	ii=(double)i;
	if(ss-x<x-ii){
		return s;
	}
	else{
		return i;
	}

}

static PyObject *Cfunctions_Conv(PyObject *self,PyObject *args){

double rot;
int len;
double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

double *CA;
PyObject *PyA;

double *CB;
PyObject *PyB;

double *MX;
PyObject *PyMX;

double *MN;
PyObject *PyMN;

PyArg_ParseTuple(args,"OOOOOOOdi",&PyLam,&PyFlux,&PyNFlux,&PyA,&PyB,&PyMN,&PyMX,&rot,&len);
CLam = ARRAYD(PyLam);
CFlux = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);
CA = ARRAYD(PyA);
CB = ARRAYD(PyB);
MN = ARRAYD(PyMN);
MX = ARRAYD(PyMX);

PyObject *myList = PyList_New(len);
PyObject *ReferenceValue;

double PI= 3.14159265;
double lux = 299792.458;
double suma,norm,G,c1,c2,c3,x;
int i = 0, ini = 0, fin = 0, j, min, max;

/*while( i < len ){

	if( CLam[i] >= 3530.0 && ini == 0){
		ini = i;
	}	
	if( CLam[len-1-i] <= 8350.0 && fin == 0){
		fin = len-1-i;
	}
	if( ini != 0 && fin != 0 ){
		break;
	}
	i++;

}
*/
while( i < len ){
	if( CLam[i] >= 3530.0){
		ini = i;
		break;
		}
	i++;
	}
i = len - 1;
while(i >= 0){
	if( CLam[i] <= 8350.0){
		fin = i;
		break;
		}
	i--;
	}

i=0;
while( i < ini ){
	ReferenceValue = PyFloat_FromDouble(CFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	i++;
}

i=fin+1;
while( i < len ){
	ReferenceValue = PyFloat_FromDouble(CFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	i++;
}

i=ini;

while( i <= fin ){
	min = i;
	j=i;
	while( j-1 > 0 ){
		if( CLam[j-1] < MN[i] ){
			min = j;
			break;
		}
		j=j-1;
	}
	max = i;
	j=i;
	while( j+1 < len-1 ){
		if( CLam[j+1] > MX[i] ){
			max = j;
			break;
		}
		j++;
	}
	
	c1 = 2.*(1.-CA[i]-CB[i]);
	c2 = 0.5*PI*(CA[i]+2.*CB[i]);
	c3 = 4.*CB[i]/3.;

	suma = 0.0;
	norm = 0.0;
	j=min;
	while( j <= max ){
		x = ((CLam[i]-CLam[j])/CLam[i])*(lux/rot);
		G = c1*sqrt(1.-x*x)+c2*(1.-x*x)-c3*pow((1.-x*x),1.5);
		suma = suma + CFlux[j]*G;
		norm = norm + G;
		j++;	
	}
	CNFlux[i] = suma/norm;

	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);

	i++;

}
PyObject *MyResult = Py_BuildValue("O",myList);

Py_DECREF(myList);   
return MyResult;

}

static PyObject *Cfunctions_res(PyObject *self,PyObject *args){


int len;
double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

double *dev;
PyObject *Pydev;

PyArg_ParseTuple(args,"OOOOi",&PyLam,&PyFlux,&PyNFlux,&Pydev,&len);
CLam = ARRAYD(PyLam);
CFlux = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);
dev = ARRAYD(Pydev);
PyObject *myList = PyList_New(len);
PyObject *ReferenceValue;
double suma,norm,G;
int i = 0, ini = 0, fin = 0, j, min, max;
while( i < len ){

	if( CLam[i] >= 3700.0 && ini == 0){
		ini = i;
	}	
	if( CLam[len-1-i] <= 8000.0 && fin == 0){
		fin = len-1-i;
	}
	if( ini != 0 && fin != 0 ){
		break;
	}
	i++;

}


i=ini;

while( i <= fin ){

	j=i;
	while( j >= 0 ){
		if( CLam[j-1] < CLam[i] - 5*dev[i] ){
			min = j;
			break;
		}
		j=j-1;
	}
	j=i;
	while( j < len ){
		if( CLam[j+1] > CLam[i] + 5*dev[i] ){
			max = j;
			break;
		}
		j++;
	}
	
	
	suma = 0.0;
	norm = 0.0;
	j=min;
	while( j <= max ){
		G = exp(-0.5*(CLam[i]-CLam[j])*(CLam[i]-CLam[j])/(dev[i]*dev[i]));
		suma = suma + CFlux[j]*G;
		norm = norm + G;
		j++;	
	}
	
	CNFlux[i] = suma/norm;
	
	i++;

}

i=0;
while(i<len){

	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	i++;

}
PyObject *MyResult = Py_BuildValue("O",myList);

Py_DECREF(myList);   
return MyResult;

}

static PyObject *Cfunctions_macro(PyObject *self,PyObject *args){

double mac;
int len,lar;

double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

double *CA;
PyObject *PyA;

double *CB;
PyObject *PyB;

double *MX;
PyObject *PyMX;

double *MN;
PyObject *PyMN;

PyArg_ParseTuple(args,"OOOOOOOdii",&PyLam,&PyFlux,&PyNFlux,&PyA,&PyB,&PyMN,&PyMX,&mac,&len,&lar);
CLam = ARRAYD(PyLam);
CFlux = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);
CA = ARRAYD(PyA);
CB = ARRAYD(PyB);
MN = ARRAYD(PyMN);
MX = ARRAYD(PyMX);

PyObject *myList = PyList_New(len);
PyObject *ReferenceValue;

double suma,norm,y,x,wavemac,m,n;
int i = 0, j, min, max,k=0;



i=0;

while( i < len ){
	min = i;
	j=i;
	while( j-1 > 0 ){
		if( CLam[j-1] < MN[i] ){
			min = j;
			break;
		}
		else{
			min = j;
		}
		j=j-1;
	}
	max = i;
	j=i;
	while( j+1 < len-1 ){
		if( CLam[j+1] > MX[i] ){
			max = j;
			break;
		}
		else{
			max = j;		
		}
		j++;
	}

	wavemac = CLam[i] * mac / 299792.458;
	
	suma = 0.0;
	norm = 0.0;
	j=min;
	while( j <= max ){
		x = sqrt( pow( ( CLam[i] - CLam[j] ),2 ) ) / wavemac;
		k = 0;
		while(k < lar - 1){
			if(x >= CA[k] && x < CA[k+1]){
				m = ( CA[k] - CA[k+1] ) / ( CB[k] - CB[k+1] );
				n = CB[k] -m*CA[k];
				y = m*x + n;
				break;
			}
			else{
				y = 0.;
			}
			k++;
		}
		suma = suma + CFlux[j]*y;
		norm = norm + y;
		j++;	
	}
	CNFlux[i] = suma/norm;
	
	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);

	i++;

}

PyObject *MyResult = Py_BuildValue("O",myList);
Py_DECREF(myList);   
return MyResult;

}

static PyObject *Cfunctions_macrot(PyObject *self,PyObject *args){

double mac,rot,dx,dy;
int len;

double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

double *CA;
PyObject *PyA;

double *MX;
PyObject *PyMX;

double *MN;
PyObject *PyMN;

PyArg_ParseTuple(args,"OOOOOOddidd",&PyLam,&PyFlux,&PyNFlux,&PyA,&PyMN,&PyMX,&mac,&rot,&len,&dx,&dy);
CLam = ARRAYD(PyLam);
CFlux = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);
CA = ARRAYD(PyA);
MN = ARRAYD(PyMN);
MX = ARRAYD(PyMX);

PyObject *myList = PyList_New(len);
PyObject *ReferenceValue;

double norm,y,x,lux=299792.458,dlamR,I_nu,F1,F2,sumaj,normj,fmac,wtemp,cos,sin;
int i = 0, j, min, max;

fmac = mac/lux;

i=0;

while( i < len ){
	//printf("%f\t%f\n",mac,rot);
	if(i > 0 && i < len-1){
		//if(i%1000==0.){
		//	printf("%f\t%f\n",CLam[i],CNFlux[i]);
		//}
		//printf("%f\t%f\n",CLam[i],CNFlux[i]);
		min = i;
		j=i;
		while( j-1 > 0 ){
			if( CLam[j-1] < MN[i] ){
				min = j;
				break;
			}
			else{
				min = j;
			}
			j=j-1;
		}
		max = i;
		j=i;
		while( j+1 < len-1 ){
			if( CLam[j+1] > MX[i] ){
				max = j;
				break;
			}
			else{
				max = j;		
			}
			j++;
		}
		
		//printf("%d\t%d\n",min,max);
		//printf("%f\t%f\t%f\n",lux*(CLam[min]-CLam[i])/CLam[i],CLam[i],lux*(CLam[max]-CLam[i])/CLam[i]);
		x = -1.;
		norm = 0.;
		while( x <= 1. ){
			//printf("%f\t%f\n",CLam[i],CLam[i]*(1+x*rot/lux));
			wtemp = CLam[i]*(1+x*rot/lux);
			dlamR =x*rot/lux;
			y = -sqrt(1.-pow(x,2));
			while( y <= sqrt(1.-pow(x,2)) ){
				j = min;
				sumaj = 0.;
				normj = 0.;
				while(j <= max){
					cos = sqrt(1. - pow(x,2) - pow(y,2) );
					sin = sqrt(pow(x,2)+pow(y,2));
					I_nu = CFlux[j] * ( 1. - CA[j]+ CA[j] * cos );
					//I_nu = 1.;
					//printf("%f\t%f\t%f\n",x,y,I_nu);
					if(fmac>0){
						if(sqrt(pow(cos,2)) > 1e-20 && sqrt(pow(sin,2)) > 1e-20) {
							//printf("opt1\n");
							F1 = exp( -pow( ((CLam[i]-CLam[j])/CLam[i] - dlamR)/ (fmac * cos ),2 ) ) / cos;
							F2 = exp( -pow( ((CLam[i]-CLam[j])/CLam[i] - dlamR)/ (fmac * sin ),2 ) ) / sin;
							//printf("%f\n",I_nu*(F1+F2));
							
							sumaj = sumaj + I_nu * (F1+F2) * 0.5 * (CLam[j+1] - CLam[j-1]);
							normj = normj + (F1 + F2) * 0.5 * (CLam[j+1]-CLam[j-1]);
						}
						else if(sqrt(pow(sin,2)) > 1e-20){
							//printf("opt2\n");
							F2 = exp( -pow( ((CLam[i]-CLam[j])/CLam[i] - dlamR)/ (fmac * sin ),2 ) ) / sin;
							//printf("%f\n",I_nu*(F1+F2));
							//printf("%f\n",I_nu*(F1+F2));
							sumaj = sumaj + I_nu*F2*0.5*(CLam[j+1]-CLam[j-1]);
							normj = normj + F2*0.5*(CLam[j+1]-CLam[j-1]);
						}
						else if(sqrt(pow(cos,2)) > 1e-20){
							//printf("opt3\n");
							F1 = exp( -pow( ((CLam[i]-CLam[j])/CLam[i] - dlamR)/ (fmac * cos ),2 ) ) / cos;
							//printf("%f\n",I_nu*(F1+F2));
							//printf("%f\n",I_nu*(F1+F2));
							sumaj = sumaj + I_nu*F1*0.5*(CLam[j+1]-CLam[j-1]);
							normj = normj + F1*0.5*(CLam[j+1]-CLam[j-1]);
						}

					}
					else{
						if( pow(wtemp-CLam[j],2) <= pow(wtemp -  CLam[j+1],2) && pow(wtemp -  CLam[j],2) <= pow(wtemp -  CLam[j-1],2) ){
							//printf("%f\t%f\n",wtemp,CLam[j]);
							sumaj = sumaj + I_nu*0.5*(CLam[j+1]-CLam[j-1]);
							normj = normj + 0.5*(CLam[j+1]-CLam[j-1]);
						}

					}
					j++;
			
				}
				if(sumaj >0.0){
					CNFlux[i] = CNFlux[i] + (sumaj/normj)*dx*dy;
					norm = norm + dx*dy;
				}
				y = y + dy;
			}
			x = x + dx;
		}
		
		if(norm > 0){
			CNFlux[i] = CNFlux[i]/norm;
		}
		else{
			CNFlux[i] = CFlux[i];
		}
		
	}
	else{
		CNFlux[i] = CFlux[i];
	}
	//printf("%f\n",CNFlux[i]);
	//printf("%f\n",CNFlux[i]);
	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	
	i++;

}


PyObject *MyResult = Py_BuildValue("O",myList);
Py_DECREF(myList);   
return MyResult;

}

// the following only works in Python 2

static PyMethodDef CfunctionsMethods[] = {
{"Conv", Cfunctions_Conv, METH_VARARGS, "¿que hace esto?"},{"res", Cfunctions_res, METH_VARARGS, "¿que hace esto?"},{"macro", Cfunctions_macro, METH_VARARGS, "¿que hace esto?"},{"macrot", Cfunctions_macrot, METH_VARARGS, "¿que hace esto?"},
{NULL, NULL, 0, NULL}
};

void initCfunctions(void){
(void) Py_InitModule("Cfunctions", CfunctionsMethods);
}

// something like this should work in Python 3,
// but MR is not familiar enough with C extensions to get the syntax right

//static struct PyMethodDef Cfunctions = 
//{
//    PyMethodDef_HEAD_INIT,
//    "Cfunctions",
//    "",
//    -1,
//    CfunctionsMethods
//};
//
//PyMODINIT_FUNC PyInit_Cfunctions(void)
//{
//    return PyModule_Create(&Cfunctions);
//}
