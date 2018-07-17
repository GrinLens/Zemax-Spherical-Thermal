/************************************************
* Dispersion Functions.c						*
* Written by Guy Beadie & Richard Flynn			*
* version 2.0 on Jun 30, 2008					*
* version 2.1 on Nov 23, 2010					*
* Last Revised:  Jan 16, 2014					*
************************************************/
/* Comments tweaked by Guy Beadie Jun. 1, 2018
     Clearing out obsolete comments left over from legacy versions
	 In some cases, replacing them with current information
	 No code was harmed in the modification of this file.
*/

/*
140116 Matches dispersion codes to Zemax glass catalog types as of version 13 Release 2 SP1
October 8, 2013
Includes a constant dispersion type and a model glass as codes 0 and -1 respectively
*/

/*
Defines the various dispersion equations as outlined in the Zemax user's
guide, version Aug 1, 2006, Chap 18 "Using Glass Catalogs," p. 512-514
*/

/******************************************************************************
* NOTE:  In all cases it is assumed that the wavelength L is given in microns *
******************************************************************************/
#include <math.h>
#define PI 3.1415926535897932

#define LFsq	0.23632500	// Lambda F sq in um sq (blue H)
#define Ldsq	0.34522887	// Lambda d sq in um sq (yellow He)
#define LCsq	0.43069359	// Lambda C sq in um sq (red H)


// ******* FUNCTIONS TO BE EXPORTED ***********
double Dispersion_Eqn(int Type, double L, double Cs[]);

// Function definition:
double Dispersion_Eqn(int Type, double L, double Cs[])
{
	int j;
	double n, nsq, invLsq, Lsq, tmp;
	double nd, Vd, A, B;
	switch (Type)
	{
	case -1:  // Model glass, as calculated from Ditteon "Modern Geometrical Optics," Eqns. (1.37-1.39)
		nd = Cs[0];
		Vd = Cs[1];
		Lsq = L*L;
		B = LFsq*LCsq*(nd-1.0)/Vd/(LCsq-LFsq);
		A = nd - B/Ldsq;
		n = A + B/Lsq;
		break;

	case 0:  // Constant, non-dispersive index
		n = Cs[0];
		break;

	case 1:  // Schott formula
		Lsq = L*L;
		invLsq = 1.0/Lsq;
		nsq = 0.0;
		for (j=5; j>=2; j--)
			nsq = (nsq + Cs[j])*invLsq;
		nsq = nsq + Cs[0] + Cs[1]*Lsq;
		n = sqrt(nsq);
		break;

	case 2: // Sellmeier 1
		Lsq = L*L;
		nsq = 1.0;
		for (j=0; j<=2; j++)
			nsq += Cs[2*j]*Lsq / (Lsq - Cs[2*j+1]);
		n = sqrt(nsq);
		break;

	case 3: // Herzberger
		Lsq = L*L;
		tmp = 1.0 / (Lsq - 0.028);
		n = Lsq*(Cs[3] + Lsq*(Cs[4] + Lsq*Cs[5]));
		n += Cs[0] + tmp*(Cs[1] + tmp*Cs[2]);
		break;

	case 4:  // Sellmeier 2
		Lsq = L*L;
		nsq = 1.0 + Cs[0] + Cs[1]*Lsq/(Lsq - Cs[2]*Cs[2]) + Cs[3]/(Lsq - Cs[4]*Cs[4]);
		n = sqrt(nsq);
		break;

	case 5: // Conrady
		Lsq = L*L;
		tmp = exp(3.5*log(L));  // lambda to the 3.5 power
		n = Cs[0] + Cs[1]/L + Cs[2]/tmp;
		break;

	case 6: // Sellmeier 3
		Lsq = L*L;
		nsq = 1.0;
		for (j=0; j<=3; j++)
			nsq += Cs[2*j]*Lsq / (Lsq - Cs[2*j+1]);
		n = sqrt(nsq);
		break;

	case 7: // Handbook of Optics Formula 1
		Lsq = L*L;
		nsq = Cs[0] + Cs[1]/(Lsq - Cs[2]) - Cs[3]*Lsq;
		n = sqrt(nsq);
		break;

	case 8: // Handbook of Optics Formula 2
		Lsq = L*L;
		nsq = Cs[0] + Cs[1]*Lsq/(Lsq - Cs[2]) - Cs[3]*Lsq;
		n = sqrt(nsq);
		break;

	case 9: // Sellmeier 4
		Lsq = L*L;
		nsq = Cs[0];
		for (j=0; j<=1; j++)
			nsq += Cs[2*j+1]*Lsq / (Lsq - Cs[2*j+2]);
		n = sqrt(nsq);
		break;

	case 10:  // Extended formula
		Lsq = L*L;
		invLsq = 1.0/Lsq;
		nsq = 0.0;
		for (j=7; j>=2; j--)
			nsq = (nsq + Cs[j])*invLsq;
		nsq = nsq + Cs[0] + Cs[1]*Lsq;
		n = sqrt(nsq);
		break;

	case 11: // Sellmeier 5
		Lsq = L*L;
		nsq = 1.0;
		for (j=0; j<=4; j++)
			nsq += Cs[2*j]*Lsq / (Lsq - Cs[2*j+1]);
		n = sqrt(nsq);
		break;

	case 12:  // Extended formula 2
		Lsq = L*L;
		invLsq = 1.0/Lsq;
		nsq = 0.0;
		for (j=5; j>=2; j--)
			nsq = (nsq + Cs[j])*invLsq;
		nsq = nsq + Cs[0] + (Cs[1] + (Cs[6] + Cs[7]*Lsq)*Lsq)*Lsq;
		n = sqrt(nsq);
		break;

	case 13:  // Extended formula 3
		Lsq = L*L;
		invLsq = 1.0/Lsq;
		nsq = 0.0;
		for (j=8; j>=3; j--)
			nsq = (nsq + Cs[j])*invLsq;
		nsq = nsq + Cs[0] + (Cs[1] + Cs[2]*Lsq)*Lsq;
		n = sqrt(nsq);
		break;
		
	default:
		n = -1.0;
		break;
	}
	return n;
}

