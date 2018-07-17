/************************************************
* Space to Stack Functions.c					*
* Written by Richard Flynn & Guy Beadie			*
* version 2.0 on Jun 30, 2008					*
* version 2.1 on Nov 23, 2010					*
* Last Revised:  May 12, 2014					*
************************************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Clearing out obsolete comments left over from legacy versions
	 In some cases, replacing them with current information
	 No code was harmed in the modification of this file.
*/

#include <math.h>
#include <stdio.h>
#include "usersurf.h"

// CONSTANT and DATA STRUCTURE definitions:  Contained in the following header.
#include "Asphere GRIN.h"	// Defines data structures & global constants.

// LOCAL CONSTANTS
#define PI 3.1415926535897932

// ******* FUNCTIONS TO BE EXPORTED ***********
int Get_n_and_Grad_n(double *n, double GradN[], double x, double y, double z, FIXED_DATA *FD, GrinPtype *GrinPs, int JustIndex);

// ******* LOCAL FUNCTIONS ***********
double Get_LayerZ_and_Gradient(double GradLZ[], double x, double y, double z, MoldPtype *MoldPs, int JustIndex);

// ******* NONLOCAL FUNCTIONS ***********
int Stack_n_and_dndzeta(double *n, double *dndzeta, double zeta, GrinPtype *GrinPs, MatPtype *MatPs, int JustIndex);
int Decode_MoldPs(MoldPtype *MoldPs, FIXED_DATA *FD);
int Decode_MatPs(MatPtype *MatPs, FIXED_DATA *FD);
int Decode_ThermPs(ThermPtype *ThermPs, FIXED_DATA *FD);
int Thermal_Scale_MoldPs(MoldPtype *MoldPs, ThermPtype *ThermPs);

/****************************************************************************/
//			EXPORTED FUNCTION DEFINITIONS:
/****************************************************************************/

//  Note: m_To has now been removed.  Zeta now refers to z position within GRIN sphere, relative to the front lens vertex.

	// OLD COMMENTS on calculating index (still illustrates the basic approach):

// The steps to take here are:
//	(1) Convert Zemax's (x,y,z) to LayerZ: the layer height in the GRIN stack
//  (2) From zeta = LayerZ/m_To, get the fractional weight of material A at this position, etaA
//	(3) Compute: n = sqrt(nA*nA*etaA + nB*nB*(1-etaA))
//
	// OLD COMMENTS on calculating gradient of index (still illustrates the basic approach):

  // Now to compute the gradient of the index function.  The real location (x,y,z) is converted to a stack
  // layer LayerZ(x,y,z).  This stack layer is converted to a normalized stack location zeta = LayerZ / m_To.
  // Finally, the index is computed as a function of zeta.  Therefore, to get the gradient of the index, we
  // compute:

  //	dn/dx  =  d(Index)/d(etaA) * d(etaA)/d(zeta) * d(zeta)/d(LayerZ) * d(LayerZ)/dx
  //	dn/dy  =  d(Index)/d(etaA) * d(etaA)/d(zeta) * d(zeta)/d(LayerZ) * d(LayerZ)/dy
  //	dn/dz  =  d(Index)/d(etaA) * d(etaA)/d(zeta) * d(zeta)/d(LayerZ) * d(LayerZ)/dz
  //
  // d(Index)/d(etaA) is computed from the expression:
  // 
  //	n = sqrt( etaA * nA^2  +  (1-etaA) * nB^2 )
  //
  // from which we get
  //
  //	d(Index) / d(etaA)  =  (nA^2 - nB^2) / (2 n)
  //
  // The second term is obtained by Get_dEtaA_dz(), the third term is simply 1/m_To, and the gradient of deltaZ
  // is given by Get_Gradient_LayerZ():

int Get_n_and_Grad_n(double *n, double GradN[], double x, double y, double z, FIXED_DATA *FD, GrinPtype *GrinPs, int JustIndex)
{
	/* Passed values
	*n			// Expects return value of refractive index
	GradN[]		// Expects return value of refractive index gradient components in x,y,z
	x, y, z		// x, y and z positions within lens, respectively
	*FD			// Fixed data structure, lots of info.  To be decomposed into local data structures
	JustIndex	// Flag: if == 0, calculate n & GradN[].  Else, only calculate n.
	*/

	// Local Variables
	int ret_val=0;		// Return value, init as 0
	int i;				// Loop counter
	double LayerZ;		// Position within stack
	double GradLZ[3];	// Gradient of stack curvature
	double dndLayerZ;	// d(n)/d(LayerZ)

	// Structures Variables derived from FD
		// Create memory space for the structures
		//	GrinPs are passed as a parameter to the function
	MoldPtype MoldPsVals;
	MatPtype MatPsVals;
		// Create pointers to those structures, so you can use the "->" operator below.  Otherwise, "." operator is necessary.
	MoldPtype *MoldPs = &MoldPsVals;
	MatPtype *MatPs = &MatPsVals;
		// Create Thermal Model memory structure
	ThermPtype ThermPsVals;
	ThermPtype *ThermPs = &ThermPsVals;

	// Read in mold parameters, MoldPs:
	ret_val = Decode_MoldPs(MoldPs,FD);
	// Read in glass parameters, MatPs:
	ret_val = Decode_MatPs(MatPs,FD) || ret_val;

	if (ret_val != 0)		// Throw an error if problems with decoding
		return ret_val;				// currently decoding never passes an error...

	// Apply thermal deformation to MoldPs
	Decode_ThermPs(ThermPs, FD);	// Read in the thermal parameters
	Thermal_Scale_MoldPs(MoldPs, ThermPs);

	// ** Main Body of Function ** //
	
	LayerZ = Get_LayerZ_and_Gradient(GradLZ, x, y, z, MoldPs, JustIndex);	// Find position within stack & gradient along stack

	ret_val = Stack_n_and_dndzeta(n, &dndLayerZ, LayerZ, GrinPs, MatPs, JustIndex);	// Calculates n & dn/dzeta

	// Calculate GradN if JustIndex is set to zero
	if ( (ret_val == 0)	&& (JustIndex == 0) )	// If no errors have occured and JustIndex=0, calculate GradN[]
	{ 
		// GradLZ has already been calculated above, when JustIndex == 0
		for (i=0; i<3; i++)
			GradN[i] = dndLayerZ * GradLZ[i];	// Return the index gradient in 3 dimensions
												// dn/da = d(n)/d(LayerZ) * d(LayerZ)/d(a), where a = (x,y,z) respectively
	}

	return ret_val;
}


/****************************************************************************/
//			LOCAL FUNCTION DEFINITIONS:
/****************************************************************************/

double Get_LayerZ_and_Gradient(double GradLZ[], double x, double y, double z, MoldPtype *MoldPs, int JustIndex)
{
	// Universal variables
	double rsq, zorigin;
	double LayerZ;
	double recip_RGrin;				// Calculate 1 / RGrin
	double m_cv;					// Local values for readability
	double zdiff;					// Temporary difference value to save a subtract
	double cv_sign;					// Sign of the mold curvature

	// Copy local values
	m_cv		= MoldPs->m_cv;

	// Calculations
	if (m_cv == 0.0)					// Flat GRIN stack
	{
		// Calculate Gradient of LayerZ
		if (JustIndex == 0)				// If we're calculating gradient
		{
			GradLZ[0] = 0.0;
			GradLZ[1] = 0.0;
			GradLZ[2] = 1.0;
		}		
		return z;		// This is the LayerZ value
	}

	rsq = x*x+y*y;	// Only need r^2 if GRIN stack isn't flat

	// Spheroids Only
	if (m_cv > 0.0)			// Handles sign flips for direction of curvature
		cv_sign = 1.0;		// Positive curvature
	else
		cv_sign = -1.0;		// Negative curvature
	zorigin = 1 / m_cv;		// Unpolished origin, Originally: 1.0 / (m_cv*(1.0+m_k))
	zdiff = z-zorigin;
	// Calculate LayerZ
	LayerZ = zorigin - cv_sign*sqrt(zdiff*zdiff + rsq);
	// Calculate Gradient of LayerZ
	if (JustIndex == 0)				// If we're calculating gradient
	{
		recip_RGrin = 1.0 / sqrt(zdiff*zdiff + rsq);	// Reciprocal of R_GRIN
		GradLZ[0] = -cv_sign * x * recip_RGrin;
		GradLZ[1] = -cv_sign * y * recip_RGrin;
		GradLZ[2] = -cv_sign * zdiff * recip_RGrin;
	}
	return LayerZ;
}

