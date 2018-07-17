/************************************************
* Stack to Composition Functions.c				*
* Written by Guy Beadie	& Richard Flynn			*
* version 2.0 on Jun 30, 2008					*
* version 2.1 on Nov 23, 2010					*
* Last Revised:  Feb  5, 2014					*
************************************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Clearing out obsolete comments left over from legacy versions
	 In some cases, replacing them with current information
	 No code was harmed in the modification of this file.
*/

#include <math.h>
#include <stdio.h>

// CONSTANT and DATA STRUCTURE definitions:  Contained in the following header.
#include "Asphere GRIN.h"	// Defines data structures & global constants.

// LOCAL CONSTANTS
#define PI 3.1415926535897932
#define AMAX 11				// Max polynomial size in Poly_Deriv()

// ******* FUNCTIONS TO BE EXPORTED ***********
int Stack_n_and_dndzeta(double *n, double *dndzeta, double zeta, GrinPtype *GrinPs, MatPtype *MatPs, int JustIndex);

// ******* LOCAL FUNCTIONS ***********
int Flatline_GRIN(double *n, double *dndzeta, double zeta, double nPs[], int JustIndex);
int Polynomial_GRIN(double *n, double *dndzeta, double zeta, double nPs[], int JustIndex);
double Poly_Eval(double x, int m, double Cs[]);
double Poly_Deriv(double x, int m, double Cs[]);


// ******* NONLOCAL FUNCTIONS ***********
double Dispersion_Eqn(int Type, double L, double Cs[]);

/**************  DESCRIPTION *****************
These functions describe the variation in material composition as a function of normalized
stack height, zeta.  Zeta used to be defined as position in the stack z divided by the mold
thickness m_To, such that valid stack positions were bound between [0,1].  Though the
nomenclature has persisted, we no longer enter a separate value for m_To, now simply assumed
to be 1 unit of Zemax lens dimension (millimeter, inch, etc.)  Therefore, there is no "out
of bounds" value for zeta anymore.

The composite is presumed to be a mixture of material A and material B, parameterized by the
volume fraction of material A: etaA.  If the refractive index of the pure materials is given
by nA and nB, then the index of a composite is given by:

    n(etaA) = sqrt[ etaA * nA^2  +  (1-etaA) * nB^2 ]

Composition of the stack is parameterized by the index of refraction at a single, design
wavelength Lref.  The variation of index is described by a polynomial:

    n(zeta, Lref)  ==  nref  =  sum{ n[j] * zeta^j }

where j ranges from 0 to the maximum polynomial order, and n[j] are the index coefficients
entered by the user in Zemax (or OpticStudio) to define the GRIN distribution.

To find the index at any other wavelength L, the weight fraction etaA is found first via:

    etaA  =  [ nref^2 - nB(Lref)^2 ] / [ nA(Lref)^2 - nB(Lref)^2 ]

Then, the index n(L) is found from the composite equation given above:

    n(L)  =  sqrt[ etaA * nA(L)^2  +  (1-etaA) * nB(L)^2 ]

Note the careful distinction between the index values nA,nB at the design wavelength Lref
versus the index values at the wavelength L actively being raytraced.

To calculate the index gradient, dn/dz, "simple" application of the chain rule is used
to get the right quantity.

*****************************/

int Stack_n_and_dndzeta(double *n, double *dndzeta, double zeta, GrinPtype *GrinPs, MatPtype *MatPs, int JustIndex)
{
	/* PASSED VALUES
	*n			// Expects return value of refractive index
	*dndzeta	// Expects return value of d(n)/d(zeta)
	zeta		// Normalized position within stack
	*GrinPs		// Gradient index parameters, several of which are used below
	*MatPs		// Material parameters, several of which are used below
	JustIndex	// Flag: if == 0, calculate n & dndzeta.  Else, only calculate n.
	*/

	// LOCAL VARIABLES
	double Lref, L;				// Reference wavelength, wavelength
	double nref;				// Reference refractive index
	double ns;					// Refractive index calculated at current wavelength
	int MatTypeCurr;			// Current material type within loop
	double *MatCsCurr;			// Points to array of coefficients of the current material within loop
	double nrefA, nrefB;		// Index of two materials at reference wavelength
	double nA, nB;				// Index of two materials at current wavelength
	double etaA;				// Material fraction of nA vs nB
	double nref_sq_diff;		// Intermediate calculation of nB^2-nA^2
	double dnrefdzeta;			// d(nref)/d(zeta)
	double dndetaA;				// d(n)/d(etaA)
	double detaAdnref;			// d(etaA/d(nref)
	int radixshift;				// Row-wise shift into a 2D array represented by a 1D array
	int ret_val = 0;			// Error return variable
	double sign_nref;			// Sign of the current nref variable, to be transmitted to calculated n & dndzeta

	// *** Main Body of Function *** //

	// * Initialization Block * //

	// Read in wavelength & reference wavelength
	Lref	= MatPs->L0;
	L		= MatPs->L;

	// * Interpret GRIN Distribution Function Block * //

	// Based on chosen GRIN stack distribution,
	//	calculate nref & dnrefdzeta, the reference refractive index & its gradient

	switch( GrinPs->nFnType )	// Switch based on the Grin Index Distribution Function Type
	{							//	It's an enumerated type, see Asphere GRIN.h
	case Flatline:
		// Calculate nref & dnref/dzeta
		ret_val = Flatline_GRIN(&nref, &dnrefdzeta, zeta, GrinPs->nPs, JustIndex);
		break;
	case Polynomial:	// Polynomial distribution
		// Calculate nref & dnref/dzeta
		ret_val = Polynomial_GRIN(&nref, &dnrefdzeta, zeta, GrinPs->nPs, JustIndex);
		break;
	default:	// Failed to choose a valid distribution
		return -1;			// Return error -1
		break;
	}
	if (ret_val != 0) return ret_val;	// Check for errors
	
	// * Calculate n & dn/dzeta Block * //

	// Calculate reference index & index of material A
	radixshift	= 0;										// Row-wise shift into a 1D array used as a 2D array
	MatTypeCurr = (int) MatPs->DispPs[radixshift + 0];		// Decode material type
	MatCsCurr	= &MatPs->DispPs[radixshift + 1];			// Decode Cs from DispPs
	nrefA	= Dispersion_Eqn(MatTypeCurr, Lref, MatCsCurr);	// Find index at reference wavelength
	if (nrefA == -1.0) {	// Throw error if dispersion function gives a bad result
		return -1; }
	if (L != Lref) { 
		nA	= Dispersion_Eqn(MatTypeCurr, L, MatCsCurr);	// Find index at current wavelength, store value
		if (nA == -1.0) {	// Throw error if dispersion function gives a bad result
			return -1; }
	}
	else {		// Don't bother with dispersion equation second time if L == Lref... it's just nrefA
		nA  = nrefA;
	}

	// Calculate reference index & index of material B
	radixshift	= MATP_MEM_LEN;								// Row-wise shift into a 1D array used as a 2D array
	MatTypeCurr = (int) MatPs->DispPs[radixshift + 0];		// Decode material type
	MatCsCurr	= &MatPs->DispPs[radixshift + 1];			// Decode Cs from DispPs
	nrefB	= Dispersion_Eqn(MatTypeCurr, Lref, MatCsCurr);	// Find index at reference wavelength
	if (nrefB == -1.0) {	// Throw error if dispersion function gives a bad result
		return -1; }
	if (L != Lref) { 
		nB	= Dispersion_Eqn(MatTypeCurr, L, MatCsCurr);	// Find index at current wavelength, store value
		if (nB == -1.0) {	// Throw error if dispersion function gives a bad result
			return -1; }	
	}
	else {		// Don't bother with dispersion equation second time if L == Lref... it's just nrefB
		nB  = nrefB;
	}

	// Check for errors
	if (nrefA == nrefB)	// Throw error if materials have identical index values
	{					// Return an error to the user
		*n = 0;
		*dndzeta = 0;
		return -1;
	}

	// Calculate refractive index at wavelength L
	nref_sq_diff = nrefB*nrefB - nrefA*nrefA;			// nB^2 - nA^2, for re-use
	// Use two materials
	// Determine sign of nref, sign_nref
	if (nref < 0.0) sign_nref = -1.0;
	else sign_nref = 1.0;
	// Calculate index at wavelength L, based on weight fraction of material A (vs material B)
	etaA = ( nrefB*nrefB - nref*nref ) / nref_sq_diff;			// Calculate the weight fraction of material A
	ns = sign_nref*sqrt( etaA * nA*nA + (1 - etaA) * nB*nB );	// Calculate the index at current wavelength
																	//	sign_nref takes care of +-sqrt() cases, based on nref's sign
	// Calculate d(n)/d(zeta)
	if (JustIndex != 0)	// If JustIndex is set (not zero), set d(n)/d(zeta) to zero
	{
		*dndzeta = 0;
	}
	else	// JustIndex == 0, so calculate the slope of the refractive index w.r.t. zeta
	{
		// Calculate d(n)/d(zeta)
			// Follow the chain rule: d(n)/d(zeta) = d(n)/d(etaA)*d(etaA)/d(nref)*d(nref)/d(zeta)
			// Earlier we calculated d(nref)/d(zeta)
		dndetaA = ( nA*nA - nB*nB ) / ( 2.0 * ns);			// Calculate d(n)/d(etaA)
		detaAdnref = -2.0 * nref / nref_sq_diff;			// Calculate d(etaA)/d(nref)
		*dndzeta = dndetaA*detaAdnref*dnrefdzeta;			// Return d(n)/d(zeta), calculated by chain rule
	}

	// Check n_Bflag
	if (GrinPs->n_BFlag != 0) {		// When n_Bflag set
		if (ns > nB) {
			ns = nB;		// High boundary: If index exceeds boundary material, set equal to boundary
			*dndzeta = 0;	// Set zero gradient since now in homogeneous material
		}
		if (ns < nA) {
			ns = nA;		// Low boundary: If index exceeds boundary material, set equal to boundary
			*dndzeta = 0;	// Set zero gradient since now in homogeneous material
		}
	}

	// Calculate return values
	*n = ns;		// Pass ns as the n value

	return ret_val;
}


/****************************************************************************/
//			LOCAL FUNCTION DEFINITIONS:
/****************************************************************************/

// ************* GRIN distribution functions:

int Flatline_GRIN(double *n, double *dndzeta, double zeta, double nPs[], int JustIndex)
{
	/* PASSED VALUES
	*n			// Expects return value of refractive index
	*dndzeta	// Expects return value of d(n)/d(zeta)
	zeta		// Normalized position within stack
	nPs[]		// Gradient index function coefficients.  Should be passed GrinPs->nPs
	JustIndex	// Flag: if == 0, calculate n & dndzeta.  Else, only calculate n.
	*/

	// LOCAL CONSTANTS
	
	// Calculate refractive index, via polynomial
	*n = 1.5;	// Find reference index at current position

	// Gradient is zero, since value is constant!
	*dndzeta = 0;

	return 0;
}

int Polynomial_GRIN(double *n, double *dndzeta, double zeta, double nPs[], int JustIndex)
{
	/* PASSED VALUES
	*n			// Expects return value of refractive index
	*dndzeta	// Expects return value of d(n)/d(zeta)
	zeta		// Normalized position within stack
	nPs[]		// Gradient index function coefficients.  Should be passed GrinPs->nPs
	JustIndex	// Flag: if == 0, calculate n & dndzeta.  Else, only calculate n.
	*/

	// LOCAL CONSTANTS
	const int maxPoly = MAX_POLY;	// Max degree of polynomial to employ
									// BE SURE maxPoly <= IP_NUM, or arrays could overflow
	
	// Calculate refractive index, via polynomial
	*n = Poly_Eval(zeta,maxPoly,nPs);	// Find reference index at current position

	if (JustIndex != 0)	// If JustIndex is set (not zero), return from here and set d(n)/d(zeta) to zero
	{
		*dndzeta = 0;
		return 0;
	}
	// Else, JustIndex == 0, calculate d(n)/d(zeta)

	// Calculate gradient of refractive index along stack, via polynomial derivative
	*dndzeta = Poly_Deriv(zeta, maxPoly, nPs);	// Calculate d(nref)/d(zeta)

	return 0;
}

// ************* Support functions for GRIN distributions:


// Given polynomial coefficients Cs[0]..Cs[m], Poly_Eval() computes:
//   Poly_Eval = Cs[0] + Cs[1] x^1 + Cs[2] x^2 + ... + Cs[m] x^m
//

double Poly_Eval(double x, int m, double Cs[]) {
	// m is the _index_ value of the max term
	double sum;
	int j;

	if (m==0) return Cs[m];
	sum = Cs[m];
	for (j=m-1; j>=0; j--)
		sum = sum*x + Cs[j];
	return sum;
}


// Given polynomial coefficients Cs[0]..Cs[m], Poly_Deriv() computes
// the derivative of the original polynomial:
//   Poly_Eval = Cs[1] + 2*Cs[2] x^1 + 3*Cs[3] x^2 + ... + m*Cs[m] x^(m-1)
//
double Poly_Deriv(double x, int m, double Cs[]) {
	// m is the _index_ value of the max term of the original polynomial
	double sum, dCs[AMAX-1];
	int j;

	if (m==0) return 0.0;
	for (j=0; j<=m-1; j++)
		dCs[j] = (j+1)*Cs[j+1];
	sum = Poly_Eval(x, m-1, dCs);
	return sum;
}



