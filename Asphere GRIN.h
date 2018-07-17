/********************************
* Asphere GRIN.h				*
*  for v4.0.0					*
* Written by Richard Flynn		*
* First Version: Oct  7, 2010	*
* version 2.1 on Nov 23, 2010	*
* Last Revised : May 12, 2014	*
********************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Clearing out obsolete comments left over from legacy versions
	 In some cases, replacing them with current information
	 No code was harmed in the modification of this file.
	 .. ALMOST no code:
	    * Default material variable name changed from SAN17 to SAN,
		  to reflect true source of dispersion equation - it's from
		  Zemax's SAN definition (not PolymerPlus's SAN17 definition)
		* Changed nFnType variable type to "int" from "enum GrinFuncType".
		  This caused insidious bug in another DLL (Ball & Radial), so
		  changed here prophylactically.
*/

/* Enumerates data structures for passing parameters within Asphere GRIN code modules */

// CONSTANTS

#define VERSION_NAME "MGRINv4.0.0"	// Name of version to display in Zemax
#define AUMAX 8		// Max number of terms user is allowed in polynomial expansion
					//   - Determines number of columns labelled & counted in lens table
					//   - Range of indices given by s_As[0]...s_As[AUMAX-1]
#define GP_NUM 10	// Number of glass parameters to be used in dispersion equations
					// Note, each glass actually has GP_NUM+1 parameters.  The first parameter is the type.
					// BEWARE: This is long enough for Dispersion_Eqn() case being used
#define IP_NUM 10	// Number of index parameters used to describe the non-straight GRIN stack
#define MAT_NUM 2	// Default number of different materials in stack
#define XDATA_MAX 200	// Highest row number of FD->xdata[]
						// Used to check for overflow in MatPs->DispPs[][]
#define MAX_POLY 9	// Maximum polynomial degree for a GRIN distribution
					//	Used in polynomial GRIN distribution function
					// BE SURE that MAX_POLY <= IP_NUM, or arrays could overflow

	// FD ARRAY READOUT CONSTANTS
	// note: FD array positions are defined by counting with 1 at first position.
	//		 But, the C arrays FD->param[] and FD->xdata[] have 1 extra member, so addressing can begin at 1 instead of 0.
#define SURFP_START 1		// Start of 7th member of SurfPs->As array in FD->xdata array
										// The first 6 members of SurfPs->As array located in FD->param[3..8]
#define MOLDP_LEN 1			// Length of MoldPs in FD->xdata array
#define TEMPP_LEN 3			// Length of Global thermal parameters in FD->xdata array
#define THERM_NUM 9			// Number of thermal modeling parameters per material

	// Derived readout constants
#define SURFP_LEN (AUMAX - 6)							// Length of portion of SurfPs->As array within FD->xdata
#define MOLDP_START (SURFP_START + SURFP_LEN)			// Start of MoldPs in FD->xdata array
#define TEMPP_START (MOLDP_START + MOLDP_LEN)			// Location of Global thermal parameters in FD->xdata array
#define GRIN_TYPE_START (TEMPP_START + TEMPP_LEN)		// Location of GrinPs->nFnType in FD->xdata array
#define GRIN_FLAG_START (GRIN_TYPE_START + 1)			// Location of GrinPs->n_BFlag in FD->xdata array
#define GRIN_START (GRIN_FLAG_START + 1)				// Start of 1st member of GrinPs->nPs array in FD->xdata array
#define GRIN_LEN (IP_NUM)								// Length of portion of GrinPs->nPs array within FD->xdata
#define MATP_START (GRIN_START + GRIN_LEN)				// Start of MatPs in FD->xdata array
#define MATP_MEM_LEN (GP_NUM + THERM_NUM + 1)			// Length of one material in MatP
#define MAT_LEN (MAT_NUM * MATP_MEM_LEN)				// Length of material data in FD->xdata array

//  ENUMERATED TYPES, for switching readability

enum GrinFuncType			// Enumerates which Grin Index Distribution function to use
{	Flatline = 0,
	Polynomial = 1};

enum MaterialType			// Enumerates Material Types for initializing MatPs->DispPs
{	PMMA	= 1,			// Note: Case 0 serves as default value
	SAN	= 2};

//	STRUCTURE DEFINITIONS, for MoldPs, MatPs, GrinPs, SurfPs

typedef struct
{
	double m_cv;		/* Mold curvature of zeta=0 surface */
}MoldPtype;

typedef struct
{
	double L0;							/* Design wavelength of GRIN optic */
	double L;							/* Current wavelength of ray trace */
	double *DispPs;						/* Dispersion parameters per material */
										/* Pointer to a 1D array masquerading as a 2D array */
										/* DispPs[n][i] gives glass parameter i of material n
											Each glass has 1 type parameter, followed by GP_NUM parameters
											Hence, MATP_MEM_LEN = GP_NUM+1*/
	/* Tip for addressing DispPs:
	DispPs is a 1D array masquerading as a 2D array, here's how to address it:
	Want:	DispPs[i][j]
	Use:	DispPs[i*MATP_MEM_LEN + j]
	*/
}MatPtype;

typedef struct
{
	// GRIN DISTRIBUTION Section
	double n_BFlag;				/* Boundary value flag - handles cases when index function exceeds the material limits
									if 1 then index is clamped to the upper or lower boundary (with a small slop allowance)
									if 0 then index is allowed to exceed the material limits (gives a weight fraction outside 0-1 range) */
	int nFnType;			/* Switch for type of index distribution */
	double *nPs;				/* Parameters to describe the particular distribution 
								   nPs[] is an array of length IP_NUM, but we'll refer to it by a pointer into FD
								   nPs[0] describes which kind of distribution */
}GrinPtype;

typedef struct
{
	double cv;		/* Surface curvature = 1/R */
	double k;		/* Conic constant */
	double thic;	/* Lens thickness */
	int surf;		/* Surface number, used for error reporting */
	double n1;		/* Index of material in front of lens */
	double As[AUMAX];	/* Forbes polynomial coefficient values for Forbes asphere */
}SurfPtype;

typedef struct
{
	// Thermal coefficients
	double TCE;		/* Thermal coefficient of expansion.  Needs 1e-6 scaling when used. */
	double T;		/* Current system temperature */
	double T0;		/* Initial system temperature, for unscaled lens */
	double s;		/* Length scaling factor */
}ThermPtype;
