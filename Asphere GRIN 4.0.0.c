/************************************************
* Asphere GRIN 4.0.0							*
* Written by Richard Flynn & Guy Beadie			*
* version 2.0 on Jun 30, 2008					*
* version 2.1 on Nov 23, 2010					*
* Last Revised:  May 12, 2014					*
************************************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Clearing out obsolete comments left over from legacy versions
	 In some cases, replacing them with current information
	 No code was harmed in the modification of this file.
	 .. ALMOST no code.  Default material variable name
	    in Asphere GRIN.h, called in case 7 here, changed
		from SAN17 to SAN, to reflect true source of dispersion
		equation - it's from Zemax's SAN definition, not PolymerPlus's
		SAN17 definition.
*/

/* v4.0.0 - RF
Adding thermal modeling to GRIN lens

The thermal model automatically accounts for the thermal changes to:
 - Even asphere coefficient terms
 - GRIN polynomial coefficient terms
 - Refractive index changes to GRIN materials

Aspheric surface coefficients merely need to be scaled by the overall thermal expansion
(or contraction) of the lens.

As a reminder, the gradient index in this code is modeled by the mixing of two pure
materials, A and B, such that the index n is found by

  n(x,y,z,L) = sqrt[ etaA * nA(L)^2  + (1-etaA) * nB(L)^2 ]

where etaA is the volume fraction of material A (versus material B) at position
(x,y,z) and nA(L),nB(L) are the indices of refraction of materials A and B at
wavelength L.

Temperature-dependent values are computed based on the assumptions that the whole lens
volume expands unfiormly, that GRIN material expanded to a new location has the same
relative material mixture as it did at its "room temperature" location, and that the
pure materials A and B have their temperature-dependent index values calculated as if
they were bulk materials at that temperature.

In other words, at a volume location (X,Y,Z) within a non-room-temperature lens, the
first step is project back to the room-temperature volume location (x,y,z).  From the
lens definition, defined at room temperature, one can identify the material volume
fraction etaA at this point.  One finds the temperature-dependent index at (X,Y,Z) by
computing:

  n(T) = sqrt[ etaA * nA(T)^2  +  (1-etaA) * nB(T)^2 ]

In practice, the code works a bit differently.  Rather than compute (x,y,z) from (X,Y,Z)
while leaving the GRIN polynomial coefficients constant, the polynomial GRIN coefficients
are scaled with temperature (under the assumption of uniform thermal expansion) and
(X,Y,Z) are used to compute etaA.

This allows legacy code in "Space to Stack Functions.c" and "Stack to Composition
Functions.c" to remain nearly unchanged upon adopting the thermal model.
*/

/* v2.0.10guy - tweaked by Guy Beadie 101117
Changes made to the default definition of the material types, to reflect
recent updates to the models introduced in Dispersion Funtions 2.0.10guy.c
*/

/*
Errors fixed preceded by // ERR_FIXED:
						 // <prev line(s) of code
						 <new line(s) of code
*/

/*
Designed to incorporate a Zemax-LIKE Even Asphere surface, a conically-molded GRIN stack, a material
distribution in the GRIN stack given by a 1D polynomial function, and two user-parameter-
defined materials that make up the composite GRIN.

The Even Asphere definition is different from Zemax in that we scale our coefficients:
   Zemax values Az[n]:   sag(r)  =  sum{ Az[n] * y^n }
   Our values Anrl[n]:   sag(r)  =  sum{ Anrl[n] * 10^(-n) * y^n }

In other words, our even asphere coefficients are 10^(n) times bigger than Zemax's coefficients
for polynomial terms of power 'n'.  This helps the nonlinear optimization algorithms, which
otherwise have trouble distinguishing between values as tiny as 10^-22 or 10^-24, which are
not uncommon for the higher-order sag terms in even asphere surfaces.

The code is spread out over 9 files:
	(1) Asphere GRIN.c - this file, containing the DLL definition
	(2) Even Asphere Function.c - contains surface-specific functions for generalized aspheric polish
	(3) Space to Stack Functions.c - translates Zemax coordinates (x,y,z) into stack positions &
									returns index values
	(4) Stack to Composition Functions.c - converts from stack position to local material composition
	(5) Dispersion Functions.c - contains material dispersion equations for pure-material indices
	(6) Decode FD.c - translates the FD data structure into smaller logical structures
	(7) Thermal Functions.c - calculates thermal perturbation of Asphere GRIN parameters
	(8) usersurf.h - header file with Zemax-specific structure definitions for user defined DLL surfaces
	(9) Asphere GRIN.h - header file which defines NRL-specific Asphere GRIN constants, structures, types
*/

#include <windows.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "usersurf.h"

// CONSTANT and DATA STRUCTURE definitions:  Contained in the following header.
#include "Asphere GRIN.h"	// Defines data structures & global constants.

/*
ERROR HANDLING (from Zemax manual):

The convention ZEMAX uses internally is to return a zero value if the DLL computed
a meaningful result, and no error occurred. Otherwise, the DLL should return -1.
The exception is when ray tracing, either paraxial or real. If the ray misses the
surface, it should return the surface number. If the ray total internally reflects
(TIR's) then the return value should be the negative of the surface number.
*/

// ******* NONLOCAL FUNCTIONS ***********
double Get_sag(double r, SurfPtype *SurfPs);
int Ray_Trace_To_Even_Asphere(SurfPtype *SurfPs, USER_DATA *UD, FIXED_DATA *FD);
int Get_n_and_Grad_n(double *Index, double GradN[], double x, double y, double z, FIXED_DATA *FD, GrinPtype *GrinPs, int JustIndex);
int Decode_SurfPs(SurfPtype *SurfPs, FIXED_DATA *FD);
int Decode_GrinPs(GrinPtype *GrinPs, FIXED_DATA *FD);
int Decode_MoldPs(MoldPtype *MoldPs, FIXED_DATA *FD);
int Decode_ThermPs(ThermPtype *ThermPs, FIXED_DATA *FD);
int Thermal_Scale_GrinPs(GrinPtype *GrinPs, ThermPtype *ThermPs);
int Thermal_Scale_MoldPs(MoldPtype *MoldPs, ThermPtype *ThermPs);
int Thermal_Scale_SurfPs(SurfPtype *SurfPs, ThermPtype *ThermPs);

// ******* LOCAL FUNCTIONS ***********
int InitializeMaterial(double *DispPpos, enum MaterialType Mat);
int __declspec(dllexport) APIENTRY UserDefinedSurface(USER_DATA *UD, FIXED_DATA *FD);


/*************************************************************************
**************************************************************************
**************************************************************************
**************************************************************************
**************************************************************************/

BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
{
   return TRUE;
}


int __declspec(dllexport) APIENTRY UserDefinedSurface(USER_DATA *UD, FIXED_DATA *FD)
	{
	int i, j, ret_val;
	int PosInFD;				// Tracks position within FD
	int breakflag;				// Flag for exiting a case statement, set within a loop
	int JustIndex = 0;			// Flag: only calculate index (not dn/dz) if == 1 (actually, !=0)
	int NMats = MAT_NUM;		// Number of materials being used.  Controls some loops.
	int tempi;					// Temporary int variable for intermediate calculations
	char temp[10], temp1len[2];
	char tempchar;
	double power;
	double r, sag;
	double GRIN_n, Grad_n[3];
	SurfPtype SurfPdata;			// Create an instance of SurfP data
	SurfPtype *SurfPs = &SurfPdata;	// Create pointer to that SurfP data
	GrinPtype GrinPdata;			// Create an instance of GrinP data
	GrinPtype *GrinPs = &GrinPdata;	// Create pointer to that GrinP data
	MoldPtype MoldPdata;			// Create an instance of MoldP data (to determine asymmetry)
	MoldPtype *MoldPs = &MoldPdata;	// Create pointer to that MoldP data (same purpose)
	ThermPtype ThermPdata;			// Create an instance of ThermP data
	ThermPtype *ThermPs = &ThermPdata;	// Create pointer to that ThermP data

	// Begin main body of the DLL:
	switch(FD->type)
   	{
      case 0:
      	/* ZEMAX is requesting general information about the surface */
         switch(FD->numb)
         	{
            case 0:
            	/* ZEMAX wants to know the name of the surface */
		         /* do not exceed 12 characters */
		         strcpy(UD->string,VERSION_NAME);
               break;
            case 1:
            	/* ZEMAX wants to know if this surface is rotationally symmetric */
					// is symmetric, so return any character in the string
				strcpy(UD->string,"1");
               break;
            case 2:
            	/* ZEMAX wants to know if this surface is a gradient index media */
               /* it IS, so return any character in the string */
            	 strcpy(UD->string,"1");
            	break;
            }
         break;
      case 1:
      	/* ZEMAX is requesting the names of the parameter columns */
         /* the value FD->numb will indicate which value ZEMAX wants. */
         /* Only "q" in parameter 1 is used for this surface type */
         /* returning a null string indicates that the parameter is unused. */
         switch(FD->numb)
         	{
            case 1:
            	/* All GRINs must use parameter 1 as Delta T !!!!!!!!!!!!!! */
            	strcpy(UD->string, "Delta T");
               break;
            case 2:
            	strcpy(UD->string, "Design WaveL");
               break;
            case 3:
            	strcpy(UD->string, "2nd Order Asph");
               break;
            case 4:
            	strcpy(UD->string, "4th Order Asph");
               break;
            case 5:
            	strcpy(UD->string, "6th Order Asph");
               break;
            case 6:
            	strcpy(UD->string, "8th Order Asph");
               break;
            case 7:
            	strcpy(UD->string, "10th Order Asph");
               break;
            case 8:
            	strcpy(UD->string, "12th Order Asph");
               break;
            default:
            	UD->string[0] = '\0';
            	break;
            }
      	break;
      case 2:
      	/* ZEMAX is requesting the names of the extra data columns */
         /* the value FD->numb will indicate which value ZEMAX wants. */
         /* returning a null string indicates that the extradata value is unused. */
		//
		/* Hybrid if/then + switch structure to allow for ranges based on constants */
		//
		// Check for nPs[0] to nPs[IP_NUM - 1]
		//
		// Read in & decode the GrinPs...need values to name columns on the fly
		Decode_GrinPs(GrinPs, FD);
		// Continue with nPs[0] to nPs[IP_NUM - 1]
		if ( ( GRIN_START <= FD->numb )  &&  ( FD->numb < (GRIN_START + IP_NUM) ) )
		{
			switch(GrinPs->nFnType)	// Switch based on the Grin Index Distribution Function Type
			{						//	It's an enumerated type, see Asphere GRIN.h
			case Flatline:
				// Name columns nPs[1], etc, by default
				tempi = FD->numb - (GRIN_START);	// Store position within nPs[] in tempi
													//	Start counting from 0, so no addition to GRIN_START
				_itoa(tempi,temp,10);				// Translate position value into a string
													// Name "n_Ps[1]", etc for polynomial coefficients
				strcpy(UD->string, "n/a n_Ps[");	// Use normal naming convention for all rows, prepend "n/a"
				strcat(UD->string,temp);
				strcat(UD->string,"]");
				// Reads, e.g. n_Ps[0]
				break;
			case Polynomial:
				// Name columns nPs[1], etc, by default
				tempi = FD->numb - (GRIN_START);	// Store position within nPs[] in tempi
													//	Start counting from 0, so no addition to GRIN_START
				_itoa(tempi,temp,10);				// Translate position value into a string
													// Name "n_Ps[1]", etc for polynomial coefficients
				if (tempi > MAX_POLY)		// Test the position value against the max polynomial coefficient
					strcpy(UD->string, "n/a n_Ps[");	// If too high, prepend "n/a" to name of column
				else
					strcpy(UD->string, "n_Ps[");		// If not too high use normal naming convention
				strcat(UD->string,temp);
				strcat(UD->string,"]");
				// Reads, e.g. n_Ps[0]
				break;
			default:
				// Name columns nPs[1], etc, by default
				_itoa(FD->numb - (GRIN_START),temp,10);		// Store position within nPs[]
															// Start counting from 0, so no addition to GRIN_START
				strcpy(UD->string, "n_Ps[");
				strcat(UD->string,temp);
				strcat(UD->string,"]");
				// Reads, e.g. n_Ps[0]
				break;
			}
			break;
		}
		//
		// Check for s_As[6] to s_As[AUMAX - 1]
		//
		if ( ( SURFP_START <= FD->numb )  &&  ( FD->numb < (SURFP_START + SURFP_LEN) ) )
		{
			tempi = FD->numb - (SURFP_START - 6);			// Store position within SurfPs->As[]
															// Start counting from 6, thus the -6
			_itoa(2*(tempi+1),temp,10);						// s_As[i] stores the 2*(i+1) order even asphere term
															// Thus we convert 2*(i+1) to string here
			strcpy(UD->string,temp);
			strcat(UD->string, "th Order Asph");
			// Reads, e.g. 14th Order Asph
			break;
		}
		//
		// Check for DispPs[0][0] to DispPs[NMats-1][MATP_MEM_LEN-1]
		//
		NMats = MAT_NUM;	// Note, we're now using exactly two materials

		breakflag = 0;
		for (i=0; i<NMats ; i++)	// "i" represents the material number
		{
			// Check for DispPs[i][0], which is the dispersion model type, n_Type
			// Column name is Mat[X] Disp Fn, rather than n_Type.  Easier meaning to user
			if ( (MATP_START + i*MATP_MEM_LEN) == FD->numb )
			{
				_itoa((i+1),temp,10);	// Store material number, i for Mat[i]
										// The value shown to user is i+1 (count from 1, not zero)
				strcpy(UD->string, "Mat[");
				strcat(UD->string,temp);
				strcat(UD->string,"] Disp Fn");
				// Reads, e.g. Mat[0] Disp Fn
				//	The variable itself will be MatPs->DispPs in the first array position for each material
				breakflag = 1;
				break;	// Leave the for loop if you find a value in this range
			}
			// Check for DispPs[i][j>0], the glass parameters, DispPs[][]
			if ( ( (MATP_START + i*MATP_MEM_LEN) < FD->numb )  &&  ( FD->numb < (MATP_START + GP_NUM + 1 + i*MATP_MEM_LEN) ) )	// nPs[1] to nPs[7]
			{
				_itoa((i+1),temp,10);	// Store material number, i for Mat[i]
										// The value shown to user is i+1 (count from 1, not zero)
				j = FD->numb - (MATP_START + i*MATP_MEM_LEN);	// Store DispPs[i][j] position j within material
				tempchar = (char)('A'+ j-1);	// Store an uppercase character corresponding to j
												// when j==1, temp1len[0]=='A'
				sprintf(temp1len,"%c",tempchar);
				strcpy(UD->string, "Mat[");
				strcat(UD->string,temp);
				strcat(UD->string,"] ");
				strcat(UD->string,temp1len);
				// Reads, e.g. Mat[0] C
				breakflag = 1;
				break;	// Leave the for loop if you find a value in this range
			}
			// Check for Thermal Material Parameters in DispPs[i][j]
			if ( ( (MATP_START + GP_NUM + 1 + i*MATP_MEM_LEN) <= FD->numb )  &&  ( FD->numb < (MATP_START + MATP_MEM_LEN + i*MATP_MEM_LEN) ) )	// Thermal Material Parameters
			{
				_itoa((i+1),temp,10);	// Store material number, i for Mat[i]
										// The value shown to user is i+1 (count from 1, not zero)
				j = FD->numb - (MATP_START + i*MATP_MEM_LEN);	// Store DispPs[i][j] position j within material
				j = j - (GP_NUM + 1);	// Change meaning of j: now a location within thermal material parameters from 0..8
				// Label based on individual values of j
				//	Add Mat[X] label when j=3..8
				if ( (j >= 3) && (j <= 8) )
				{
					strcpy(UD->string, "Mat[");
					strcat(UD->string,temp);
					strcat(UD->string,"] ");
				}
				//	Add [X] label when j=0 or 1
				if (j == 0 || j == 1 )
				{
					strcpy(UD->string, "[");
					strcat(UD->string,temp);
					strcat(UD->string,"] ");
				}
				switch(j)
					{
					case 0:
            			strcat(UD->string, "TCE x 1E-6");
						break;
					case 1:
            			strcat(UD->string, "Char. Temp");
						break;
					case 2:
            			strcpy(UD->string, "Thermal n Model");
						break;
					case 3:
            			strcat(UD->string, "D0");
						break;
					case 4:
            			strcat(UD->string, "D1");
						break;
					case 5:
            			strcat(UD->string, "D2");
						break;
					case 6:
            			strcat(UD->string, "E0");
						break;
					case 7:
            			strcat(UD->string, "E1");
						break;
					case 8:
            			strcat(UD->string, "Lam_tk");
						break;
					default:
						UD->string[0] = '\0';
            			break;
					}
				// Reads, e.g. Mat[1] TCE
				breakflag = 1;
				break;	// Leave the for loop if you find a value in this range
			}
		}
		if (breakflag) break;	// If the value was found within the loop, exit the case statement
		//
		// Check for values in fixed positions
		//
        switch(FD->numb)
			{
			case (MOLDP_START):
            	strcpy(UD->string, "m_cv");
				break;
			case (TEMPP_START):
            	strcpy(UD->string, "Initial Temp");
				break;
			case (TEMPP_START+1):
            	strcpy(UD->string, "Current Temp");
				break;
			case (TEMPP_START+2):
            	strcpy(UD->string, "TCE Lens x 1E-6");
				break;
			case GRIN_TYPE_START:
				strcpy(UD->string, "GRINProfileType");
				// Describes the variable GrinPs->nFnType. "GRINProfileType" is more descriptive to user.
				break;
			case GRIN_FLAG_START:
				strcpy(UD->string, "n_BFlag");
				break;
            default:
				UD->string[0] = '\0';
            	break;
            }
      	break;
      case 3:
      	/* ZEMAX wants to know the sag of the surface */

			// Read in surface parameters:
			Decode_SurfPs(SurfPs,FD);		// Read in the surface parameters
			Decode_ThermPs(ThermPs, FD);	// Read in the thermal parameters
			Thermal_Scale_SurfPs(SurfPs, ThermPs);	// Thermally scale the SurfPs

			UD->sag1 = 0.0;
			UD->sag2 = 0.0;
			r = sqrt(UD->x*UD->x + UD->y*UD->y);

			sag = Get_sag(r,SurfPs);
			UD->sag1 = sag;
			UD->sag2 = UD->sag1;

      	break;
      case 4:
      	/* ZEMAX wants a paraxial ray trace to this surface */
         /* x, y, z, and the optical path are unaffected, at least for this surface type */
         /* for paraxial ray tracing, the return z coordinate should always be zero. */
         /* paraxial surfaces are always planes with the following normals */
         /* we will ignore the aspheric terms, even the quadratic one, since it has a */
         /* meaning that is hard to interpret if q != 0.0 */

         UD->ln =  0.0;
         UD->mn =  0.0;
         UD->nn = -1.0;

		 // Notes on calculating refractive index (no code used)
 		 /***************************************************************************
		 * CORRECTION 3/4/2011: ZEMAX wants FD->n1 and FD->n2 for refraction.   	*
		 *	 Using any other values confuses the optimizer!							*
		 * For a GRIN surface, ZEMAX does __NOT__ want the real refracted angle.    *
		 *   For compatibility with GRIN propagation, we must return the refracted  *
		 *   angle as if the index at the surface is FD->n2.					    *
		 ***************************************************************************/

		 power = (FD->n2 - FD->n1)*FD->cv;
         if ((UD->n) != 0.0)
         	{
            (UD->l) = (UD->l)/(UD->n);
            (UD->m) = (UD->m)/(UD->n);

            (UD->l) = (FD->n1*(UD->l) - (UD->x)*power)/(FD->n2);
            (UD->m) = (FD->n1*(UD->m) - (UD->y)*power)/(FD->n2);

            /* normalize */
            (UD->n) = sqrt(1/(1 + (UD->l)*(UD->l) + (UD->m)*(UD->m) ) );
            /* de-paraxialize */
            (UD->l) = (UD->l)*(UD->n);
            (UD->m) = (UD->m)*(UD->n);
            }
         break;
 
	  case 5:
      	/* ZEMAX wants a real ray trace to this surface */
		 // Read in surface parameters:
		 Decode_SurfPs(SurfPs,FD);		// Read in the surface parameters
		 Decode_ThermPs(ThermPs, FD);	// Read in the thermal parameters
		 Thermal_Scale_SurfPs(SurfPs, ThermPs);	// Thermally scale the SurfPs

		 // Single function handles ray propagation, finding surface normal, and refracting ray
		 //  Returns FD->surf if ray lands outside lens semi-diameter
		 //  Returns -FD->surf if total internal reflection
		 //  Otherwise returns 0
		 ret_val = Ray_Trace_To_Even_Asphere(SurfPs, UD, FD);

		 // Notes on calculating refractive index (no code used)
 		 /***************************************************************************
		 * CORRECTION 3/4/2011: ZEMAX wants FD->n1 and FD->n2 for refraction.   	*
		 *	 Using any other values confuses the optimizer!							*
		 * For a GRIN surface, ZEMAX does __NOT__ want the real refracted angle.    *
		 *   For compatibility with GRIN propagation, we must return the refracted  *
		 *   angle as if the index at the surface is FD->n2.					    *
		 ***************************************************************************/

		 return ret_val;
         break;

      case 6:
      	/* ZEMAX wants the index, dn/dx, dn/dy, and dn/dz at the given x, y, z. */

         /* OLD ZEMAX DEFAULT COMMENT: This is only required for gradient index surfaces, so return dummy values */

		  // Calculate the index and its gradient
		 JustIndex = 0;	// Calculate index & its gradient below, Get_n_and_Grad_n()
 		 ret_val = Decode_GrinPs(GrinPs, FD);	// Decode GrinPs to pass to Get_n_and_Grad_n() function
		 if (ret_val != 0) return ret_val;	

		 Decode_ThermPs(ThermPs, FD);			// Read in the thermal parameters
		 Thermal_Scale_GrinPs(GrinPs, ThermPs);	// Thermally scale the GrinPs

		 ret_val = Get_n_and_Grad_n(&GRIN_n,Grad_n,UD->x,UD->y,UD->z,FD,GrinPs,JustIndex);
		 if (ret_val != 0) return ret_val;

         UD->index = GRIN_n;
         UD->dndx = Grad_n[0];
         UD->dndy = Grad_n[1];
         UD->dndz = Grad_n[2];
      	break;

      case 7:
      	/* ZEMAX wants the "safe" data. */
         /* this is used by ZEMAX to set the initial values for all parameters and extra data */
         /* when the user first changes to this surface type. */
         /* this is the only time the DLL should modify the data in the FIXED_DATA FD structure */

		// Loading safe values into FD:
		//
		// FD->param[]
		//
		FD->param[1] = 0.2;	// ZEMAX's Delta T parameter

		// Design Wavelength
		//
		FD->param[2] = 0.5875618;	// Design wavelength of the GRIN lens, set to d-line 587.5618 nm

		// Surface parameters, SurfPs
		//
		for (i=0; i<=5; i++)
			FD->param[3+i] = 0.0;	// s_As[0..5], the Even Asphere Terms
		//
		// FD->xdata[]
		//		
		// Surface parameters, SurfPs
		//
		for (i = SURFP_START ; i < ( SURFP_START + SURFP_LEN ) ; i++)
			FD->xdata[i] = 0.0;	// s_As[6..AUMAX-1]

		// Mold parameters, MoldPs
		//
		j = MOLDP_START;
		FD->xdata[j++] = 0.0;		// m_cv

		// Global temperature parameters
		j = TEMPP_START;
		FD->xdata[j++] = 20.0;		// Initial Temp, T0
		FD->xdata[j++] = 20.0;		// Current Temp, T
		FD->xdata[j++] = 0.0;		// TCE Lens x 1E-6, TCE

		// GRIN stack parameters, GrinPs
		//
		j = GRIN_TYPE_START;
		FD->xdata[j] = 1.0;					// nFnType
		j = GRIN_FLAG_START;
		FD->xdata[j] = 0.0;					// n_BFlag
		j = GRIN_START;
		FD->xdata[j++] = 1.5;				// Set nPs[0] to a non-zero value (so index is non-zero)
		for (i = 1 ; i < GRIN_LEN ; i++)	// Set nPs[1..IP_NUM-1] to zero
			FD->xdata[j++] = 0.0;			// nPs[1..IP_NUM-1]

		// Glass material parameters, MatPs
		//
		// Loop to set the MatPs->DispPs[][]
		//
		PosInFD = (MATP_START);			// Initialize position in FD
		// Populate the DispPs[i][] material parameters per material
		i = 0;		// Initialize, which material # we are on
		while ( i < 2 )		// Initialize only two materials, since we're only using two materials
		{
			switch (i)		// Set default material values, using InitializeMaterial() function
			{
			case 0:
				// Set material parameters to PMMA.
				InitializeMaterial(&FD->xdata[PosInFD],PMMA);
				break;
			case 1:
				// Set material parameters to SAN17.
				InitializeMaterial(&FD->xdata[PosInFD],SAN);
				break;
			default:
				// Set material parameters to SAN17.
				InitializeMaterial(&FD->xdata[PosInFD],SAN);
				break;
			}
			PosInFD += MATP_MEM_LEN;		//Increment position counter by one material.
			i++;	// Increment material #
		}
		if (i = 0) return -1;	// If there was no initialization of any materials, something is drastically wrong

		// Finished initializing FD with "safe" data
		break;
/*	case 12:		// For debugging purposes:
		// Testing the MatPs decoder
		Test_Decode_MatPs(FD);
		break;
*/
      }
   return 0;
   }


int InitializeMaterial(double *DispPpos, enum MaterialType Mat)
{
	// PASSED VALUES:
	//
	// *DispPpos	// Points to subsection of an array, where material parameters should be written
	// Mat			// Chooses which specific material's parameters to write into array

	int j;		// Loop counter

	// Note: Case 0 serves as default value
	switch (Mat)
	{
	case PMMA:			// PMMA Coefficients, from Zemax glass catalog MISC.AGF, Version 13 Release 2 SP1, October 8, 2013
		DispPpos[0] = 1;				// Type 1, Schott
		DispPpos[1] = 2.18645820E+000;	// Six coefficients
		DispPpos[2] = -2.44753480E-004;
		DispPpos[3] = 1.41557870E-002;
		DispPpos[4] = -4.43297810E-004;
		DispPpos[5] = 7.76642590E-005;
		DispPpos[6] = -2.99363820E-006;
		for (j=7; j<MATP_MEM_LEN; j++)
			DispPpos[j] = 0.0;			// Set remaining coefficients to zero
		break;
	case SAN:			// SAN Coefficients, from Zemax glass catalog MISC.AGF, Version 13 Release 2 SP1, October 8, 2013
		DispPpos[0] = 1;				// Type 1, Schott
		DispPpos[1] = 2.38687023E+000;	// Six coefficients
		DispPpos[2] = -1.23063994E-003;
		DispPpos[3] = 2.29467817E-002;
		DispPpos[4] = 3.69810122E-004;
		DispPpos[5] = 2.67577106E-005;
		DispPpos[6] = 2.84806182E-006;
		for (j=7; j<MATP_MEM_LEN; j++)
			DispPpos[j] = 0.0;			// Set remaining coefficients to zero
		break;
	default:
		for (j=0; j<MATP_MEM_LEN; j++)
			DispPpos[j] = 0.0;			// Set all coefficients to zero
		break;
	}
	// Set Thermal Material Values
	DispPpos[GP_NUM+2] = 20;			// Char. Temp set to 20 degrees
	//	All other thermal material values default to 0, set above
	return 0;
}
