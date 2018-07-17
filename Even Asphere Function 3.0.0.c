/************************************************
* Even Asphere Function.c						*
* Written by Richard Flynn						*
* version 2.9 on Jan 16, 2014					*
* Last Revised:  Jan 22, 2014					*
************************************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Clearing out obsolete comments left over from legacy versions
	 In some cases, replacing them with current information
	 No code was harmed in the modification of this file.
*/

/* Scaled version of Zemax's Even Asphere Surface Definition
  Based on MIL-HDBK-141, p. 5-10 to 5-11.
  Independent testing found this to be consistent with Zemax surface.
  In particular, both surfaces fail in the same way for extreme aspheres
    where high-angle, undeviated input rays intersect multiple points on
	the surface.
  On the other hand, this code runs significantly faster than NRL-written
    code which handles such cases.  Speed and compatibility with Zemax
	chosen over accuracy for extreme cases.
*/

/************************************************
Sag is given by:

z = c*r^2/ (1 + sqrt( 1 - (1+k)*c^2*r^2 ) 
       + a1*r^2 + a2*r^4 + a3*r^6 +...+ a8*r^16

where
	c = curvature
	k = conic constant
	a1..a8 = aspheric coefficients

all the a_n terms are the even asphere terms.  The a_n terms have units.

NOTE: the code in this module takes the Even Asphere coefficients as given,
but by this point "Decode FD.c" has already SCALED the values entered by the
user in Zemax's data columns.

In particular:
   a1  =  a1_Typed / 100
   a2  =  a2_Typed / 10000
     etc.

In other words, the coefficients here are 10^(-2n) times smaller than the
values typed by the user in OpticStudio. The scaling terms are provided so
that the numerical values of the surface coefficients are not so small.  It
is not uncommon for 10 mm diameter optics described by Zemax's Even Asphere
8th or 10th order surface terms to have  values near 10^-20.  While Zemax
appears to recognize its own Even Asphere variables and handle their magnitudes
appropriately for optimization, experience has shown that values this small in
User Defined surfaces become difficult to handle for Zemax's nonlinear optimizer
algorithms.  By rescaling, the coefficients become much bigger numbers, and much
easier for Zemax's optimizers to handle.
************************************************/

#include <stdio.h>
#include <math.h>
#include "usersurf.h"

// CONSTANT and DATA STRUCTURE definitions:  Contained in the following header.
#include "Asphere GRIN.h"	// Defines data structures & global constants.

// ******* FUNCTIONS TO BE EXPORTED ***********
double Get_sag(double r, SurfPtype *SurfPs);

// ******* LOCAL FUNCTIONS ***********
double Poly_Eval_Even(double x, int m, double Cs[]);
double Poly_Eval_Even_From_Zero(double x, int m, double Cs[]);


// ******* Function definitions ********

double Get_sag(double r, SurfPtype *SurfPs) {
	// Calculate the sag of an even asphere
	// z = c*r^2/ (1 + sqrt( 1 - (1+k)*c^2*r^2 ) + a1*r^2 + a2*r^4 + a3*r^6 +...+ a8*r^16

	int m = AUMAX;		// Order of even asphere polynomial
	double sag;			// Sag value, to accumulate
	double rsq;			// r^2 for speed
	double D, Den;		// Calculation vars
	
	// Calculate conic term
	rsq = r*r;
	D = 1.0 - (1.0 + SurfPs->k)*SurfPs->cv*SurfPs->cv*rsq;
	if (D < 0.0) D = 0.0;	// Round D up to 0.0, to prevent negative square root
	Den = 1.0 + sqrt(D);	// Denominator
	sag = SurfPs->cv*rsq / Den;

	// Calculate even asphere term
	sag = sag + Poly_Eval_Even(r, m, SurfPs->As);

	return sag;
}


int Ray_Trace_To_Even_Asphere(SurfPtype *SurfPs, USER_DATA *UD, FIXED_DATA *FD) {
	// Traces a ray from the original UD position to intersect the even asphere surface
	// then refracts the ray into the new material
	//   Returns FD->surf if ray lands outside lens semi-diameter
	//   Returns -FD->surf if total internal reflection
	//   Otherwise returns 0
	//   Does not account for multiple intersections with surface (Zemax Even Aspheres do not either)
	/* This is a two-step process, taken from derivation in document MIL-HDBK-141
		(I've modified the equations to include the (1+k) conic term and d*S^2 term)

	Sag equation: Z = c*S^2/(1+sqrt(1-(1+k)*c^2*S^2) + d*S^2 + e*S^4 + f*S^6 + g*S^8 + ...

	Step 1) Propagating from "tangent plane" to "tangent sphere".
	Following MIL-HDBK-141, p. 5-10 to 5-11.
	Note: Thickness & initial z position are zero, so some simplifications:
	("_" denotes subscript)

	Substitutions:
	Adivn = A/n_-1
	ncI = n_-1*cos(I)

	Eq (1): d_-1/n_-1 = 0
	Eq (2): Y_T = Y_-1
	Eq (3): X_T = Y_-1		// So... Eq(1-3) are unused.
	Eq (9): H = c(X_-1^2 + Y_-1^2)
	Eq(10): B = M_-1 - c(Y_-1*L_-1 + X_-1*K_-1)
	Eq (7): ncI = n_-1*sqrt[(B/n_-1)^2 - c*H]
	or:     ncI = sqrt[B^2 - c*H*n_-1^2]
	Eq (8): Adivn = H/(B + ncI)
	Eq (4): X = X_-1 + Adivn*K_-1
	Eq (5): Y = Y_-1 + Adivn*L_-1
	Eq (6): Z = Adivn*M_-1
	and we're done there, since we're still within same material there is no refraction.
	Note: K, L, M are the *optical* direction cosines, e.g. the direction cosine multiplied by 
	index of the current medium!  K/n = actual direction cosine, etc.
	K = K_-1
	L = L_-1
	M = M_-1

	Special case: No conic coefficients!
	Go ahead and refract into the spherical material.

	Substitutions:
	ncGam = n*cos(Gamma)

	Eq(11): ncGam = n*sqrt[(ncI/n)^2 - (n_-1/n)^2 + 1]
	or:     ncGam = sqrt[(ncI)^2 - (n_-1)^2 + n^2]
	Eq(12): Gam = ncGam - ncI
	Eq(13): K = K_-1 - X*c*Gam
	Eq(14): L = L_-1 - Y*c*Gam
	Eq(15): M = M_-1 - (Z*c - 1)*Gam


	Step 2) Propagating from "tangent sphere" to aspheric surface
	Following MIL-HDBK-141, p. 5-18 to 5-19.
	Note: Had to modify Eqs(18) & (17) to account for conic term.  Double-checked derivation 

	to find that this substitution is sufficient.

	Substitutions:
	dAp = delta_A_prime / n_-1
	Ssq = S^2

	Iterative Section:
	while (|dAp| > min_dAp) {
	Eq(16): Ssq = X^2 + Y^2
	Eq(18): W = sqrt(1 - (1+k)*c^2*Ssq)
	Eq(17): F = Z - sag(S, k, As) // Here it's just the sag equation
	Eq(17)alt: F = Z - [ c*Ssq/(1+sqrt(1-(1+k)*c^2*Ssq)) + d*S^2 + e*S^4 + f*S^6 + g*S^8 + ... ]
	Eq(22): E = c + W [ 2*d + 4e*S^2 + 6f*S^4 + 8g*S^6 + ... ]
	Eq(23): U = - X*E
	Eq(24): V = - Y*E
	Eq(25): dAp = -F*W / (K_-1*U + L_-1*V + M_-1*W)
	Eq(19): X_+1 = X + dAp*K_-1
	Eq(20): Y_+1 = Y + dAp*L_-1
	Eq(21): Z_+1 = Z + dAp*M_-1
	}
	repeats until dAp gets sufficiently small to desired level

	Substitutions:
	GncI = G*n_-1*cos(I)	// Note: This is the value *before* refraction
	GncIp = G*n*cos(I')	// Note: This is the value *after* refraction

	Non-Iterative Section:
	Eq(26): G = sqrt(U^2 + V^2 + W^2)
	Eq(27): GncI = K_-1*U + L_-1*V + M_-1*W
	Eq(28): GncIp = n*sqrt[(GncI/n)^2 - G^2*(n_-1/n)^2 + G^2]
	or    : GncIp = sqrt[GncI^2 - G^2*(n_-1)^2 + G^2*n^2]
	Eq(29): P = (GncIp - GncI)/G^2
	Eq(30): K = K_-1 + U*P
	Eq(31): L = L_-1 + V*P
	Eq(32): M = M_-1 + W*P

	where X,Y,Z are the offsets from input "tangent sphere" and K,L,M are direction cosines after refraction.
	In Zemax, want to use FD->n1 as n_-1 and FD->n2 as n!  Due to weird way of handling GRIN.
	*/

	// Iteration Limits
	const double min_dAp = 1.0e-8;

	// Step 1) Vars
	double rad;		// Radix of a square root
	double X, Y;	// Coordinates
	double Z = 0.0;	// Z begins at tangent plane = 0
	double H, B;	// Vars for step 1)
	double ncI;		// n_-1*cos(I).  I is angle of incidence.
	double Adivn;	// A/n_-1
	double K, L, M;	// Optical direction cosines (direction cosine * index)
	double np, n;	// Prior and current refractive index
	double c;		// Spherical surface curvature
	int IsAsphere;	// Flag for whether it's an asphere
	int IsZDir = 0;	// Flag for whether ray is parallel to Z axis
	int i;			// Counter
	double ncGam;	// n*cos(Gamma)
	double Gam;		// Gamma
	// Step 2) Vars
	int its;		// Number of iterations
	double Xo, Yo;	// Inital X, Y values
	double dAp;		// Residual error in path length, or Delta[A']/n_-1
	double Ssq;		// Radial cylindrical coordinate, squared
	double S;		// Radial cylindrical coordinate
	double W, U, V;	// Components of surface normal vector (not normalized)
	double F, E;	// Intermediate variables: F = residual sag, E = Approx. curvature of aspheric surface
	double GncI;	// G*n_-1*cos(I), related to angle of incidence, I
	double GncIp;	// G*n*cos(I'), related to angle of refraction, I'
	double G, P;	// G = length of surface normal vector, P = Parallel to surface normal, Gamma, divided by G
	double Cs[AUMAX];	// Coefficients for derivative


	// Get variable values from outside data
	X = UD->x;
	Y = UD->y;
	Xo = UD->x;
	Yo = UD->y;
	np = FD->n1;
	n = FD->n2;
	K = UD->l * np;
	L = UD->m * np;
	M = UD->n * np;
	c = FD->cv;

	// Check for rays parallel to Z axis first
	//	Speeds up some later steps
	if (UD->l == 0.0 && UD->m == 0.0) IsZDir = 1;

	// ******* Step 1) *******
	// Propagating from "tangent plane" to "tangent sphere"
	//
	// Check to see whether any aspheric coefficients are in play
	IsAsphere = 0;
	if (FD->k != 0.0) {	// If there's conic term, it's an asphere
		IsAsphere = 1;
	}
	else {				// Or, if any asphere terms are non-zero, it's an asphere
		for(i=0;i<AUMAX;i++)
			if (SurfPs->As[i] != 0.0) IsAsphere = 1;
	}
	// Propagate ray to spherical surface
	if (IsZDir && IsAsphere) {	// Simplify calculation if ray is parallel to Z, and surface is aspheric
								//	If surface is spherical, it's faster & simpler just to run the normal processing below
		// X and Y stay unchanged, and Z is just the sag
		Ssq = X*X + Y*Y;
		S = sqrt(Ssq);
		Z = Get_sag(S, SurfPs);		// Propagates right to the sag, whether aspheric or spherical
		// Note: If we ran this for a spherical surface, it's missing the calculation of ncI below, which breaks the refraction algorithm
		// Broken into two sections, for Steps 1) and 2), in case we can skip the aspheric part of Step 2)
	}
	else {			// Normal processing when ray is not in Z direction
		H = c*(X*X + Y*Y);
		B = M - c*(Y*L + X*K);
		rad = B*B - c*H*np*np;
		if (rad < 0.0) return(-1000);	// Return error if negative square root
		ncI = sqrt(rad);
		Adivn = H/(B + ncI);
		X = X + Adivn*K;
		Y = Y + Adivn*L;
		Z = Adivn*M;
	}
	// Done propagating
	// If it's not aspheric, go ahead and refract the exit ray and return final values
	if (IsAsphere == 0) {
		// If we're outside the semi-diameter, throw error of missed surface
		// if (sqrt(X*X + Y*Y) > FD->sdia) return (FD->surf);	// SKIPPING THIS, so vignetting handles missed surfaces: more like default Zemax behavior
		//
		// Refract the ray
		rad = ncI*ncI - np*np + n*n;
		if (rad < 0.0) return(-FD->surf);	// Indicates total internal reflection, return Zemax preferred -FD->surf
		//	**** Note: this incorrectly uses the vertex value of refractive index, rather than the value at X, Y, Z
		ncGam = sqrt(rad);
		Gam = ncGam - ncI;
		K = K - X*c*Gam;
		L = L - Y*c*Gam;
		M = M - (Z*c - 1)*Gam;
		// Record final values in UD
		/*	What Zemax wants updated in a real ray trace:
			(from case 5: of us_grin1.c)
			UD-> x,y,z
			UD-> l,m,n
			UD-> ln,mn,nn
			UD-> path
			If ray trace fails, return -FD->surf.	*/
		UD->x = X;
		UD->y = Y;
		UD->z = UD->z + Z;
		UD->l = K / n;
		UD->m = L / n;
		UD->n = M / n;
		UD->path = Adivn * np;
		// Surface normals taken from page 5-9 of document
		// But with opposite sign convention, copied from us_grin1.c, lines 151-153
		UD->ln = c*X;
		UD->mn = c*Y;
		UD->nn = (c*Z - 1);
		// Return (note: we've already checked for missed surfaces & TIR)
		//	Actually, we'll still trace missed surfaces, more like Zemax homogeneous lens behavior
		return 0;
	}
	// If it is aspheric, continue to Step 2)
	//
	// ******* Step 2) *******
	// Propagating from "tangent sphere" to aspheric surface
	//
	if (IsZDir) {	// Skip Iterative Section when ray is in Z direction
		// Calculate W, U and V to get surface normal, for refracting ray later
		// Ssq, S and Z were calculated earlier, in the first section
		W = sqrt(1 - (1+FD->k)*c*c*Ssq);
		// Setting coefficients of this equation: E = c + W [ 2d + 4e*S^2 + 6f*S^4 + 8g*S^6 + ... ]
		for (i=0;i<AUMAX;i++)	// Set up the coefficients
			Cs[i] = 2*(i+1) * SurfPs->As[i];
		E = c + W * Poly_Eval_Even_From_Zero(S, AUMAX, Cs);	// Done
		U = - X*E;
		V = - Y*E;
	}
	else {			// Normal processing when ray is not in Z direction
		// Iterative Section
		dAp = min_dAp*1000;		// Initialize to enter loop successfully
		its = 0;
		while ((fabs(dAp) > min_dAp) && its < 1000)
		{
			Ssq = X*X + Y*Y;
			S = sqrt(Ssq);
			W = sqrt(1 - (1+FD->k)*c*c*Ssq);
			F = Z - Get_sag(S, SurfPs);			// Here it's just the sag equation
			// Setting coefficients of this equation: E = c + W [ 2d + 4e*S^2 + 6f*S^4 + 8g*S^6 + ... ]
			for (i=0;i<AUMAX;i++)	// Set up the coefficients
				Cs[i] = 2*(i+1) * SurfPs->As[i];
			E = c + W * Poly_Eval_Even_From_Zero(S, AUMAX, Cs);	// Done
			U = - X*E;
			V = - Y*E;
			dAp = -F*W / (K*U + L*V + M*W);
			X = X + dAp*K;
			Y = Y + dAp*L;
			Z = Z + dAp*M;
			its++;
		}
	}
	// If we're outside the semi-diameter, throw error of missed surface
	// if (sqrt(X*X + Y*Y) > FD->sdia) return (FD->surf);	// SKIPPING THIS, so vignetting handles missed surfaces: more like default Zemax behavior
	//
	// Non-Iterative Section
	//
	// Refract the ray
	G = sqrt(U*U + V*V + W*W);
	GncI = K*U + L*V + M*W;
	rad = GncI*GncI - G*G*np*np + G*G*n*n;
	// Test for total internal reflection, return if found:
	//	note: sin(I) >= n2/n1, is equivalent to
	//        n1^2 - n2^2 >= (n1*cos(I))^2
	//        which turns out to be indicated by the sign of variable "rad" above
	if (rad < 0.0) return (-FD->surf);
	//	If there is total internal reflection, return negative surface value
	//	**** Note: this incorrectly uses the vertex value of refractive index, rather than the value at X, Y, Z
	GncIp = sqrt(rad);
	P = (GncIp - GncI)/G/G;
	K = K + U*P;
	L = L + V*P;
	M = M + W*P;
	// Record final values in UD
	/*	What Zemax wants updated in a real ray trace:
		(from case 5: of us_grin1.c)
		UD-> x,y,z
		UD-> l,m,n
		UD-> ln,mn,nn
		UD-> path
		If ray trace fails with TIR, return -FD->surf.
		If ray trace fails due to missed surface, return FD->surf	*/
	UD->x = X;
	UD->y = Y;
	UD->z = UD->z + Z;
	UD->l = K / n;
	UD->m = L / n;
	UD->n = M / n;
	UD->path = sqrt((X-Xo)*(X-Xo) + (Y-Yo)*(Y-Yo) + Z*Z);
	// Surface normals: Page 5-17 notes that U/G, V/G, W/G are direction cosines of the surface normal
	// Direction cosines are equal to the unit vector components
	// Zemax wants these to point "outward" thus the minus signs
	UD->ln = -U/G;
	UD->mn = -V/G;
	UD->nn = -W/G;
	// Return (note: we've already checked for missed surfaces & TIR)
	//	Actually, we'll still trace missed surfaces, more like Zemax homogeneous lens behavior
	return 0;
}


// Given polynomial coefficients Cs[0]..Cs[m-1], Poly_Eval_Even() computes:
//   Poly_Eval_Even = Cs[0] x^2 + Cs[1] x^4 + Cs[2] x^6 + ... + Cs[m-1] x^m*2
//
double Poly_Eval_Even(double x, int m, double Cs[]) {
	// m is number of terms
	double sum, xsq;
	int j;

	if (m==0) return Cs[m];
	xsq = x*x;
	sum = Cs[m-1]*xsq;
	for (j=m-2; j>=0; j--)
		sum = (sum + Cs[j])*xsq;
	return sum;
}

// Given polynomial coefficients Cs[0]..Cs[m-1], Poly_Eval_Even_From_Zero() computes:
//   Poly_Eval_Even_From_Zero = Cs[0] + Cs[1] x^2 + Cs[2] x^4 + ... + Cs[m-1] x^(m-1)*2
//
double Poly_Eval_Even_From_Zero(double x, int m, double Cs[]) {
	// m is number of terms
	double sum, xsq;
	int j;

	if (m==0) return Cs[m];
	xsq = x*x;
	sum = Cs[m-1];
	for (j=m-2; j>=0; j--)
		sum = sum*xsq + Cs[j];
	return sum;
}

