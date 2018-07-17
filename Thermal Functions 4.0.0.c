/********************************
* Thermal Functions.c			*
* Written by Richard Flynn		*
* First Version: May 12, 2014	*
* Last Revised : May 12, 2014	*
********************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Added some explanation of the functions
	 No code was harmed in the modification of this file.
*/

/*
The three functions defined here transform the parameters which describe the GRIN distribution
in a manner consistent with uniform expansion (or contraction) of the lens element.

For example, the ashperic surface coefficients As[j], used to describe the lens sag via

   sag = (conic sag) + sum{ As[j] * r^j }

are modified to account for thermal expansion by a linear scale factor 's' according to:

   sag = (Zemax-handled conic sag) + sum{ AsMod[j] * r^j }

where

   s = 1 + TCE * dT

   AsMod[j] = As[j] * s^(1-j)

Consider the calculation of sag at the edge of a lens.  At elevated temperature, the edge
of the lens shifts to:

   Rmax(Thot) = Rmax(Troom) * [1 + TCE * (Thot-Troom)]

Zemax handles the expansion of the conic sag, or edge thickness.  The aspheric contribution
is similarly scaled by 's', but also needs to be referenced to its original location
from the un-expanded lens.  In other words:

  Aspheric sag at Rmax(Thot) = s * Aspheric sag at Rmax(Troom)

  sag(Rmax,Thot) = (Zemax-handled conic sag)  +  s * sum{ As[j] * Rmax(Troom)^j }
                 = (Zemax-handled conic sag)  +  sum{ s * As[j] * [Rmax(Thot)/s]^j }
				 = (Zemax-handled conic sag)  +  sum{ As[j] * s / s^j * Rmax(Thot)^j }
				 = (Zemax-handled conic sag)  +  sum{ AsMod[j] * Rmax(Thot)^j }

Similar arguments are followed to define the other functions in this module.
*/

#include <math.h>
#include "usersurf.h"

// CONSTANT and DATA STRUCTURE definitions:  Contained in the following header.
#include "Asphere GRIN.h"	// Defines data structures & global constants.

// ******* FUNCTIONS TO BE EXPORTED ***********
int Thermal_Scale_GrinPs(GrinPtype *GrinPs, ThermPtype *ThermPs);
int Thermal_Scale_MoldPs(MoldPtype *MoldPs, ThermPtype *ThermPs);
int Thermal_Scale_SurfPs(SurfPtype *SurfPs, ThermPtype *ThermPs);

// These functions scale GRIN lens parameters accoring to Thermal Coefficient of Expansion (TCE)
	// Note: The constants for array readout are defined in Asphere GRIN.h

int Thermal_Scale_GrinPs(GrinPtype *GrinPs, ThermPtype *ThermPs)
	// Scales index polynomial
{
	double s;		// Scaling factor for Thermal Expansion
	double s_inv;	// 1/s
	double fac;		// Multiplication factor
	int max;		// Max order
	int j;			// Loop variable

	if (ThermPs->T == ThermPs->T0) return 0;	// Skip function if there is no delta T

	// Get scaling factor
	//	Based on formula: L' = L(1 + TCE*dT) = L*s
	//	s = (1 + TCE*dT)
	s = ThermPs->s;

	// Re-scale index polynomial
	//	Scaling of index polynomial uses formula: b_i = a_i / s^i
	//	Basically, we're streching the distribution along the z axis
	//  We'll begin with i=1 -> fac = 1/s and divide by s for increasing coeff order
	// Find max order actually used
	//	Loop from max order, until reach a non-zero value
	max = IP_NUM-1;
	while ( (GrinPs->nPs[max] == 0.0) && (max > 0) )
		max--;		// Decrement max, to min value of 0
	// Now re-scale the non-zero nPs values
	if (max > 0)	// Skip process if max == 0 ... there are no non-zero coeffs with i >= 1
	{
		// Initialize re-scaling
		fac = 1/s;
		s_inv = fac;
		// Loop
		for (j = 1; j <= max; j++)
		{
			GrinPs->nPs[j] = GrinPs->nPs[j]*fac;
			fac = fac*s_inv;
		}
	}

	return 0;
}

int Thermal_Scale_MoldPs(MoldPtype *MoldPs, ThermPtype *ThermPs)
	// Scales mold curvature
{
	// Re-scale mold curvature
	//	We are scaling up mold radius by scaling factor s
	//	This is equivalent to dividing mold curvature, cv, by s
	MoldPs->m_cv = MoldPs->m_cv / ThermPs->s;

	return 0;
}

int Thermal_Scale_SurfPs(SurfPtype *SurfPs, ThermPtype *ThermPs)
	// Scales even asphere coefficients
{
	double s;		// Scaling factor for Thermal Expansion
	double ssq_inv;	// 1/s^2
	double fac;		// Multiplication factor
	int max;		// Max order
	int j;			// Loop variable

	if (ThermPs->T == ThermPs->T0) return 0;	// Skip function if there is no delta T

	// Get scaling factor
	//	Based on formula: L' = L(1 + TCE*dT) = L*s
	//	s = (1 + TCE*dT)
	s = ThermPs->s;

	// Re-scale even asphere coefficients
	//	Scaling of a surface shape polynomial uses formula: b_i = a_i / s^(i-1)
	//	Basically, we're stretching along z (1/s^i) and making bigger bump amplitude (*s)
	//	However, even aspheres are even powers only, so i = (j+1)*2
	//  Furthermore, we'll begin with i=2 -> fac = 1/s and divide by s^2 for increasing coeff order
	// Find max order actually used
	//	Loop from max order, until reach a non-zero value
	max = AUMAX-1;
	while ( (SurfPs->As[max] == 0.0) && (max > -1) )
		max--;		// Decrement max, to min value of -1
	// Now re-scale the non-zero As values
	if (max > -1)	// Skip process if max == -1 ... there are no non-zero coeffs
	{
		// Initialize re-scaling
		fac = 1/s;
		ssq_inv = 1/s/s;
		// Loop
		for (j = 0; j <= max; j++)
		{
			SurfPs->As[j] = SurfPs->As[j]*fac;
			fac = fac*ssq_inv;
		}
	}

	return 0;
}
