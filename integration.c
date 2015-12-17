#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#define ARRAYD(p) ((double *) (((PyArrayObject *)p)->data))

struct my_g_params { double med; double sig; gsl_spline *spline; gsl_interp_accel *acc;};
struct Vint_params { double li; double mac; double rot; gsl_spline *spline; gsl_interp_accel *acc; gsl_spline *spline2; gsl_interp_accel *acc2; double Li; double Lf; double pres;};
struct conv_params { double x; double li; double mac; double rot; gsl_spline *spline; gsl_interp_accel *acc; gsl_spline *spline2; gsl_interp_accel *acc2; double Li; double Lf; double pres;};
struct mult_params { double y; double x; double li; double mac; double rot; gsl_spline *spline; gsl_interp_accel *acc; gsl_spline *spline2; gsl_interp_accel *acc2; double Li; double Lf;};
struct my_conv_params { double med; gsl_spline *spline; gsl_interp_accel *acc; gsl_spline *spline2; gsl_interp_accel *acc2;};

double my_g (double l, void * p) {
	//printf("%f\n",l);
	struct my_g_params * params = (struct my_g_params *)p;
	double med   = (params->med);
	double sig   = (params->sig);
	gsl_spline *spline = (params->spline);
	gsl_interp_accel *acc = (params->acc);
	double pi = 3.14159265359;
	double x = (l-med) / sig;
	double F_l;
	//printf ("%f\t%f\n", wav[0],flx[0]);

	//gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	//gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, len);
	//gsl_spline_init (spline, wav, flx, len);
	F_l = gsl_spline_eval (spline, l, acc);

	//printf ("%f\t%f\n", l,F_l);
	//gsl_spline_free (spline);
        //gsl_interp_accel_free (acc);
	return F_l * exp(-0.5*pow(x,2)) / sqrt(2.0*pi*sig*sig);
}

double my_conv (double l, void * p) {
	//printf("%f\n",l);
	struct my_conv_params * params = (struct my_conv_params *)p;
	double med   = (params->med);
	gsl_spline *spline = (params->spline);
	gsl_interp_accel *acc = (params->acc);
	gsl_spline *spline2 = (params->spline2);
	gsl_interp_accel *acc2 = (params->acc2);

	double x = (l - med)/med;
	double F_l,I_l;
	//printf("%f\n",x);
	F_l = gsl_spline_eval (spline, l, acc);
	I_l = gsl_spline_eval (spline2, x, acc2) / med;

	return F_l * I_l;
}

double the_integral(double li, double R, double *l, double *f, int len, gsl_spline *spline, gsl_interp_accel *acc){

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double sig,result,error,low_lim,up_lim;
	sig = li/(R*2.3548201);
	gsl_function F;
	struct my_g_params params = { li, sig, spline, acc };
	F.function = &my_g;
	F.params = &params;
	if(li-5*sig < l[0]){
		low_lim = l[0];
	}
	else{
		low_lim = li-5*sig;
	}
	
     	if(li+5*sig > l[len-1]){
		up_lim = l[len-1];
	}
	else{
		up_lim = li+5*sig;
	}	
	gsl_integration_qags(&F, low_lim, up_lim, 0, 1e-3, 1000, w, &result, &error); 
     
	//printf ("%f\n", f[0]);
     
	gsl_integration_workspace_free (w);
     
	return result;
}

double the_integral2(double li, double *d, double *inst, double *l, double *f, int len, double dinf, double dsup, gsl_spline *spline, gsl_interp_accel *acc, gsl_spline *spline2, gsl_interp_accel *acc2){
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result,error,low_lim,up_lim;
	gsl_function F;
	struct my_conv_params params = { li, spline, acc, spline2, acc2 };
	F.function = &my_conv;
	F.params = &params;
	
	if(li-dinf < l[0]){
		low_lim = l[0];
	}
	else{
		low_lim = li-dinf;
	}
	
     	if(li+dsup > l[len-1]){
		up_lim = l[len-1];
	}
	else{
		up_lim = li+dsup;
	}
	
	gsl_integration_qags(&F, low_lim, up_lim, 0, 1e-3, 1000, w, &result, &error); 
     
	//printf ("%f\n", f[0]);
     
	gsl_integration_workspace_free (w);
     
	return result;
}


double mult(double l, void * p) {
	//printf("MULT\n");
	
	double I, Fr, Ft,result,med,sig;
	double lux = 299792.458;
	double pi = 3.14159265359;
	double At = 0.5, Ar = 0.5, cost,sint,dlam,ldc;
	struct mult_params * params_mult = (struct mult_params *)p;
	double x     = (params_mult->x);
	double y     = (params_mult->y);
	double li    = (params_mult->li);
	double mac   = (params_mult->mac);
	double rot   = (params_mult->rot);
	double Li   = (params_mult->Li);
	double Lf   = (params_mult->Lf);

	gsl_spline *spline = (params_mult->spline);
	gsl_interp_accel *acc = (params_mult->acc);

	gsl_spline *spline2 = (params_mult->spline2);
	gsl_interp_accel *acc2 = (params_mult->acc2);

	
	cost = sqrt( 1. - pow(x,2) - pow(y,2) );
	sint = sqrt(pow(x,2) + pow(y,2));
	
	
	//printf("%f\t%f\t%f\n",x,y,l);
	if(l >= Li && l <= Lf){
		ldc = gsl_spline_eval(spline2, l, acc2);
		I = gsl_spline_eval(spline, l, acc) * ( 1. - ldc + ldc * cost );
	}
	else if(l < Li){
		ldc = gsl_spline_eval(spline2, Li, acc2);
		I = gsl_spline_eval(spline, Li, acc) * ( 1. - ldc + ldc * cost );
	}
	else{
		ldc = gsl_spline_eval(spline2, Lf, acc2);
		I = gsl_spline_eval(spline, Lf, acc) * ( 1. - ldc + ldc * cost );
	}
	
	
	dlam = (l - li)/ li - x*rot/lux;
	sig = 0.005;


	if(mac == 0.){
		med = li*(1 + x*rot/lux);
		result = I * exp(-0.5*pow( (l-med)/sig,2))/sqrt(2*pi*sig*sig);

	}
	else{
		Fr = (Ar/(sqrt(pi)*mac*cost)) * exp (- pow(dlam*lux/(mac*cost),2) );
		Ft = (At/(sqrt(pi)*mac*sint)) * exp (- pow(dlam*lux/(mac*sint),2) );
		result = I*(Fr + Ft);
		if(cost == 0. ){
			result = I*Ft;
		}
		else if (sint == 0.){
			result = I*Fr;
		}
		else if (cost == 0. && sint == 0. ){
			result = 0.;
		}
		else{
			result = I*(Fr + Ft);
		}
		result = lux * result / (li*(Ar+At));
	}
	//printf("%f\t%f\t%f\t%f\n",I,Fr,Ft,result);
	return result;

}

double conv(double y, void * p) {
	double lux = 299792.458;
	//printf("CONV\n");
	double result_conv, error_conv,ldc;
	struct conv_params * params_conv = (struct conv_params *)p;
	double x     = (params_conv->x);
	double mac   = (params_conv->mac);
	double rot   = (params_conv->rot);
	double li   = (params_conv->li);
	double Li   = (params_conv->Li);
	double Lf   = (params_conv->Lf);
	double pres   = (params_conv->pres);
	gsl_spline *spline = (params_conv->spline);
	gsl_interp_accel *acc = (params_conv->acc);
	gsl_spline *spline2 = (params_conv->spline2);
	gsl_interp_accel *acc2 = (params_conv->acc2);

	if(mac==0. && rot == 0.){
		ldc = gsl_spline_eval(spline2, li, acc2);
		result_conv = gsl_spline_eval(spline, li, acc) * ( 1. - ldc + ldc * sqrt(1-pow(x,2)-pow(y,2)) );
	}	
	else{
		gsl_integration_workspace * w_conv = gsl_integration_workspace_alloc (100);
		gsl_function Fmult;
		struct mult_params params_mult = { y, x, li, mac, rot, spline, acc, spline2, acc2, Li, Lf};
		Fmult.function = &mult;
		Fmult.params = &params_mult;
		//printf("%f\t%f\n",li*(1.-3.*mac/lux ),Li);

	
		gsl_integration_qags(&Fmult, li*(1.-(3.*mac+rot)/lux ), li*(1.+(3.*mac+rot)/lux ), 0, pres, 100, w_conv, &result_conv, &error_conv);
		gsl_integration_workspace_free (w_conv);
	}
	return result_conv;
	

}

double Vint (double x, void * p) {
	//printf("%f\n",x);
	//printf("VERTICAL\n");
	double result_Vint, error_Vint;
	struct Vint_params * params_Vint = (struct Vint_params *)p;
	double mac   = (params_Vint->mac);
	double rot   = (params_Vint->rot);
	double li   = (params_Vint->li);
	double Li   = (params_Vint->Li);
	double Lf   = (params_Vint->Lf);
	double pres   = (params_Vint->pres);
	gsl_spline *spline = (params_Vint->spline);
	gsl_interp_accel *acc = (params_Vint->acc);
	gsl_spline *spline2 = (params_Vint->spline2);
	gsl_interp_accel *acc2 = (params_Vint->acc2);

	gsl_integration_workspace * w_Vint = gsl_integration_workspace_alloc (100);
	gsl_function FVint;
	struct conv_params params_conv = { x, li, mac, rot, spline, acc, spline2, acc2, Li, Lf,pres };
	FVint.function = &conv;
	FVint.params = &params_conv;
	gsl_integration_qags(&FVint, (-1.)*sqrt(1. - pow(x,2)), sqrt(1. - pow(x,2)), 0, pres, 100, w_Vint, &result_Vint, &error_Vint);
	gsl_integration_workspace_free (w_Vint);
	return result_Vint;
}

double disk_integration(double li, double mac, double rot, double *l, double *f, int len, gsl_spline *spline, gsl_interp_accel *acc, gsl_spline *spline2, gsl_interp_accel *acc2, double Li, double Lf,double pres){
	double result, error;
	double pi = 3.14159265359;
	//printf("DISCO\n");
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100);

	
	gsl_function F;
	struct Vint_params params = { li, mac, rot, spline, acc, spline2, acc2, Li, Lf, pres };

	F.function = &Vint;
	F.params = &params;
	
	gsl_integration_qags(&F, -1., 1., 0, pres, 100, w, &result, &error); 
     
	//printf ("%f\n", f[0]);
     
	gsl_integration_workspace_free (w);
     
	return result/pi;
}
  
static PyObject *integration_InstConv(PyObject *self,PyObject *args){
int n_datos;
double R;

double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

PyArg_ParseTuple(args,"OOOid",&PyLam,&PyFlux,&PyNFlux,&n_datos,&R);
CLam   = ARRAYD(PyLam);
CFlux  = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);

PyObject *myList = PyList_New(n_datos);
PyObject *ReferenceValue;

gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n_datos);
gsl_spline_init (spline, CLam, CFlux, n_datos);



int i;
i = 0;
//printf("%f\n",R);
while(i < n_datos){
	
	//printf("%f\n",CLam[i]);
	if(CLam[i] - 3. * CLam[i]/(2.3*R) > CLam[0] && CLam[i] + 3. * CLam[i]/(2.3*R) < CLam[n_datos - 1]){ 
		CNFlux[i] = the_integral(CLam[i],R,CLam,CFlux,n_datos,spline,acc);
	}
	else{
		CNFlux[i] = CFlux[i];
	}
	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	i++;
}
gsl_spline_free (spline);
gsl_interp_accel_free (acc);

PyObject *MyResult = Py_BuildValue("O",myList);
Py_DECREF(myList);   
return MyResult;

}

static PyObject *integration_InstConvVarGau(PyObject *self,PyObject *args){
int n_datos;

double *CR;
PyObject *PyR;

double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

PyArg_ParseTuple(args,"OOOOi",&PyLam,&PyFlux,&PyNFlux,&PyR,&n_datos);
CLam   = ARRAYD(PyLam);
CFlux  = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);
CR = ARRAYD(PyR);

PyObject *myList = PyList_New(n_datos);
PyObject *ReferenceValue;

gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n_datos);
gsl_spline_init (spline, CLam, CFlux, n_datos);



int i;
i = 0;
//printf("%f\n",R);
while(i < n_datos){
	
	//printf("%f\n",CLam[i]);
	if(CLam[i] - 3. * CLam[i]/(2.3*CR[i]) > CLam[0] && CLam[i] + 3. * CLam[i]/(2.3*CR[i]) < CLam[n_datos - 1]){ 
		CNFlux[i] = the_integral(CLam[i],CR[i],CLam,CFlux,n_datos,spline,acc);
		
	}
	else{
		CNFlux[i] = CFlux[i];
	}
	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	i++;
}
gsl_spline_free (spline);
gsl_interp_accel_free (acc);

PyObject *MyResult = Py_BuildValue("O",myList);
Py_DECREF(myList);   
return MyResult;

}


static PyObject *integration_NotParamInstConv(PyObject *self,PyObject *args){

int n_datos, n_datos2;

double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

double *CDel;
PyObject *PyDel;

double *CInst;
PyObject *PyInst;

PyArg_ParseTuple(args,"OOOOOii",&PyLam,&PyFlux,&PyNFlux,&PyDel,&PyInst,&n_datos,&n_datos2);
CLam   = ARRAYD(PyLam);
CFlux  = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);
CDel = ARRAYD(PyDel);
CInst = ARRAYD(PyInst);

PyObject *myList = PyList_New(n_datos);
PyObject *ReferenceValue;


gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, n_datos2);
gsl_spline_init (spline2, CDel, CInst, n_datos2);

gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n_datos);
gsl_spline_init (spline, CLam, CFlux, n_datos);
//printf("Hi\n");


int i;
double dinf,dsup;
i = 0;
while(i < n_datos){
	//printf("%f\n",CLam[i]);
	dinf = CLam[i]*CDel[1]*(-1);
	dsup = CLam[i]*CDel[n_datos2 -2];
	//printf("%f\t%f\t%f\n",dinf,CLam[i],dsup);
	CNFlux[i] = the_integral2(CLam[i],CDel,CInst,CLam,CFlux,n_datos,dinf,dsup,spline,acc,spline2,acc2);
	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	i++;
}
gsl_spline_free (spline);
gsl_interp_accel_free (acc);
gsl_spline_free (spline2);
gsl_interp_accel_free (acc2);

PyObject *MyResult = Py_BuildValue("O",myList);
Py_DECREF(myList);   
return MyResult;

}

static PyObject *integration_MacRot(PyObject *self,PyObject *args){
int n_datos;
double mac,rot,pres;

double *CLam;
PyObject *PyLam;

double *CFlux;
PyObject *PyFlux;

double *CNFlux;
PyObject *PyNFlux;

double *CLdc;
PyObject *PyLdc;

PyArg_ParseTuple(args,"OOOOiddd",&PyLam,&PyFlux,&PyNFlux,&PyLdc,&n_datos,&mac,&rot,&pres);
CLam   = ARRAYD(PyLam);
CFlux  = ARRAYD(PyFlux);
CNFlux = ARRAYD(PyNFlux);
CLdc = ARRAYD(PyLdc);

PyObject *myList = PyList_New(n_datos);
PyObject *ReferenceValue;

gsl_interp_accel *acc = gsl_interp_accel_alloc ();
gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n_datos);
gsl_spline_init (spline, CLam, CFlux, n_datos);

gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, n_datos);
gsl_spline_init (spline2, CLam, CLdc, n_datos);
double ltest;
int i;
i = 0;

while(i < n_datos){
	
	//printf("%d\n",i);
	//printf("%f\n",CLam[i]);
	//ltest = gsl_spline_eval(spline2, 5500., acc2);
	//printf("HI %f\n",ltest);
	CNFlux[i] = disk_integration(CLam[i],mac,rot,CLam,CFlux,n_datos,spline,acc,spline2,acc2,CLam[0],CLam[n_datos-1],pres);
	ReferenceValue = PyFloat_FromDouble(CNFlux[i]);
	PyList_SET_ITEM(myList,i,ReferenceValue);
	printf("%d\t%f\t%f\t%f\n",i,CLam[i],CFlux[i],CNFlux[i]);
	i++;
}
gsl_spline_free (spline);
gsl_interp_accel_free (acc);
gsl_spline_free (spline2);
gsl_interp_accel_free (acc2);
PyObject *MyResult = Py_BuildValue("O",myList);
Py_DECREF(myList);   
return MyResult;

}

static PyMethodDef integrationMethods[] = {
{"InstConv", integration_InstConv, METH_VARARGS, "¿que hace esto?"},{"InstConvVarGau", integration_InstConvVarGau, METH_VARARGS, "¿que hace esto?"},{"NotParamInstConv", integration_NotParamInstConv, METH_VARARGS, "¿que hace esto?"},{"MacRot", integration_MacRot, METH_VARARGS, "¿que hace esto?"},
{NULL, NULL, 0, NULL}
};

void initintegration(void){
(void) Py_InitModule("integration", integrationMethods);
}
