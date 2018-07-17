/********************************
* Decode FD.c					*
* Written by Richard Flynn		*
* First Version: Oct  8, 2010	*
* version 2.1 on Nov 23, 2010	*
* Last Revised : May 12, 2014	*
********************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Clearing out obsolete comments left over from legacy versions
	 In some cases, replacing them with current information
	 No code was harmed in the modification of this file.
*/

#include <math.h>
#include "usersurf.h"

// CONSTANT and DATA STRUCTURE definitions:  Contained in the following header.
#include "Asphere GRIN.h"	// Defines data structures & global constants.

// ******* FUNCTIONS TO BE EXPORTED ***********
int Decode_MoldPs(MoldPtype *MoldPs, FIXED_DATA *FD);
int Decode_MatPs(MatPtype *MatPs, FIXED_DATA *FD);
int Decode_GrinPs(GrinPtype *GrinPs, FIXED_DATA *FD);
int Decode_SurfPs(SurfPtype *SurfPs, FIXED_DATA *FD);
int Decode_ThermPs(ThermPtype *ThermPs, FIXED_DATA *FD);

// These functions decode the FD data structure from Zemax into local data structures
	// Note: The constants for array readout are defined in Asphere GRIN.h

int Decode_MoldPs(MoldPtype *MoldPs, FIXED_DATA *FD)
{
	int PosInFD;		// Read-in counter

	// Read in mold parameters, MoldPs:
	PosInFD = MOLDP_START;
	MoldPs->m_cv		= FD->xdata[PosInFD++];
	return 0;
}

int Decode_MatPs(MatPtype *MatPs, FIXED_DATA *FD)
{
	// Wavelength & reference wavelength are also handled by MatPs

	// Read in lens design wavelength, the wavelength at which the index function is accurate
	//	When design wavelength is set negative, we want to copy
	//	FD->pwavelength into the design wavelength
	if ( FD->param[2] < 0 )	// For negative values
		MatPs->L0			= FD->pwavelength;	// Read in pwavelength as lens design wavelength
	else					// For positive values
		MatPs->L0			= FD->param[2];		// Read in lens design wavelength

	// Read in the current ray trace wavelength
	MatPs->L				= FD->wavelength;				// Read in the current wavelength

	// Read in DispPs array location
	MatPs->DispPs				= &FD->xdata[MATP_START];	// Read in material dispersion parameters
																// Copies pointer to DispPs sub-array within
																// FD->xdata
	// DispPs represents a 2D array of size [NMats][MATP_MEM_LEN]
	/* Tip for addressing DispPs:
	DispPs is a 1D array masquerading as a 2D array, here's how to address it:
	Want:	DispPs[i][j]
	Use:	DispPs[i*MATP_MEM_LEN + j]
	*/

	return 0;
}

int Decode_GrinPs(GrinPtype *GrinPs, FIXED_DATA *FD)
{
	GrinPs->n_BFlag				= (int) FD->xdata[GRIN_FLAG_START];	// Read in the boundary flag for the GRIN index function

	GrinPs->nFnType				= (int) FD->xdata[GRIN_TYPE_START];	// Read in the type of GRIN index distribution used
																	// nFnType is of an enumerated type, but casting as (int) works
	GrinPs->nPs					= &(FD->xdata[GRIN_START]);			// Refer to nPs[] by reference

	return 0;
}

int Decode_SurfPs(SurfPtype *SurfPs, FIXED_DATA *FD)
{
	int i;			// Loop counter
	int pos;		// Tracks position within As[] for read-in

	// Read in most of SurfPs directly from same members in FD
	//
	SurfPs->cv			= FD->cv;
	SurfPs->k			= FD->k;
	SurfPs->thic		= FD->thic;
	SurfPs->surf		= FD->surf;
	SurfPs->n1			= FD->n1;
	// Read in the rest of SurfPs from FD->param and FD->xdata
	//
	// Loop for As[0..5]
	for (i=0; i<=5; i++)
		SurfPs->As[i]	= FD->param[3+i];	// As[0..5] are in FD->param[3..8]
	// Loop for As[6..AUMAX-1]
	pos = 6;		// Start reading at As[6]
	for (i = SURFP_START ; i < ( SURFP_START + SURFP_LEN ) ; i++)
	{
		SurfPs->As[pos] = FD->xdata[i];		// As[6..AUMAX-1] are in FD->xdata[SURFP_START..( SURFP_START + SURFP_LEN - 1 )]
		pos++;		// Increment pos
	}

	return 0;
}

int Decode_ThermPs(ThermPtype *ThermPs, FIXED_DATA *FD)
{
	// Read in global thermal parameters
	ThermPs->T0			= FD->xdata[TEMPP_START];
	ThermPs->T			= FD->xdata[TEMPP_START+1];
	ThermPs->TCE			= FD->xdata[TEMPP_START+2];

	// Calculate scaling factor
	//	Based on formula: L' = L(1 + TCE*dT) = L*s
	//	s = (1 + TCE*dT)
	ThermPs->s = 1 + ThermPs->TCE*1.0e-6*(ThermPs->T - ThermPs->T0);	// TCE as-entered has 1e-6 factor missing

	return 0;
}