#ifndef __EQUATION_DUAL_VEL_CG1_CG1_HH
#define __EQUATION_DUAL_VEL_CG1_CG1_HH

#include <math.h>

#include "ShapeFunctionTetLin.hh"

class EquationDualVel_cG1cG1: public Equation {

private:

  real theta;		// theta=0.5 => cG(1)cG(1) 
  
  real delta1, delta2, r;
  real C1, C2, K;
  
	real Re, dt, t;

	real nu;

	real u1sum, u2sum, u3sum, psum, divU, divUmid, normU;

	real U1U1, U1U2, U1U3, U2U1, U2U2, U2U3, U3U1, U3U2, U3U3;
	real u1U1, u1U2, u1U3, u2U1, u2U2, u2U3, u3U1, u3U2, u3U3;

	real Ii;
	real Iii;
	real Iij;

	real vx, vy, vz;
	real uv1, uv2, uv3, uv4;

	real u1_1, u1_2, u1_3, u1_4, u2_1, u2_2, u2_3, u2_4, u3_1, u3_2,
	    u3_3, u3_4;
	real D1u1, D2u1, D3u1, D1u2, D2u2, D3u2, D1u3, D2u3, D3u3;
	real u1v, u2v, u3v, uv;
	real MASS1, MASS2, MASS3, REA1, REA2, REA3, LAP1, LAP2, LAP3, NL1,
	    NL2, NL3, LSdiv1, LSdiv2, LSdiv3;
	real UDu_UDv_1, UDu_UDv_2, UDu_UDv_3, SD1, SD2, SD3, Dp_DUv_1,
	    Dp_DUv_2, Dp_DUv_3;
	real UDu_DUv_1, UDu_DUv_2, UDu_DUv_3, DUu_DUv_1, DUu_DUv_2,
	    DUu_DUv_3, DUu_UDv_1, DUu_UDv_2, DUu_UDv_3, UDv;
	real d, tf, centerpt_x, centerpt_y, centerpt_z, x_el, y_el, z_el;
	real f_DUv_1, f_DUv_2, f_DUv_3, f_UDv_1, f_UDv_2, f_UDv_3;
	real f1, f2, f3, f_v_1, f_v_2, f_v_3, u1_v_TS, u2_v_TS, u3_v_TS,
	    Dp_v1, Dp_v2, Dp_v3;

      public:

	//-----------------------------------------
	//      User definitions of equations
	//-----------------------------------------
	 EquationDualVel_cG1cG1():Equation(3, 3) {

		pi = KW_PI;

		theta = 0.5;

		C1 = 1.0;
		C2 = 1.0;
		K = 20.0;

		Ii = 0.25;
		Iii = 0.1;
		Iij = 0.05;
	} ~EquationDualVel_cG1cG1() {
	}


	//--------- Equation definition -----------
	inline void jacobian(int el, int gpt, int tst_bf, int trial_bf,
			     int comp) {
	}


	//--------- Load definition --------------
	inline void residual(int el, int gpt, int tst_bf) {
	}




	//--------- Load definition for exact integration on tetrahedrons --------------
	inline
	    void jacobianExactInt(MV_Vector < real > &jac, int &tst_bf,
				  int &trial_bf, int &j_eq, int &el,
				  ShapeFunction *sfcn) {


    int cellsize = sfcn->GetNoNodes();

		if ((tst_bf == 0) && (trial_bf == 0)) {

			h = 2.0 * sfcn->GetCircumRadius();

			//for ( i=0; i < cellsize; i++ ) uprev_TS1(i) = COEFF(0) -> GetNodalValue( sfcn->GetNode(i) );    
			//sfcn->EvalGradient(uprev_TS1,uprev_TS1Grad);
			//for ( i=0; i < cellsize; i++ ) uprev_TS2(i) = COEFF(1) -> GetNodalValue( sfcn->GetNode(i) );    
			//sfcn->EvalGradient(uprev_TS2,uprev_TS2Grad);
			//for ( i=0; i < cellsize; i++ ) uprev_TS3(i) = COEFF(2) -> GetNodalValue( sfcn->GetNode(i) );    
			//sfcn->EvalGradient(uprev_TS3,uprev_TS3Grad);
			//for ( i=0; i < cellsize; i++ ) pprev_TS(i) = COEFF(3) -> GetNodalValue( sfcn->GetNode(i) );    
			//sfcn->EvalGradient(pprev_TS,pprev_TSGrad);

			for (int i = 0; i < cellsize; i++)
				U1(i) =
				    COEFF(4)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(U1, U1Grad);
			for (int i = 0; i < cellsize; i++)
				U2(i) =
				    COEFF(5)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(U2, U2Grad);
			for (int i = 0; i < cellsize; i++)
				U3(i) =
				    COEFF(6)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(U3, U3Grad);

			for (int i = 0; i < cellsize; i++)
				u1(i) =
				    COEFF(8)->GetNodalValue(sfcn->
							    GetNode(i));
			//sfcn->EvalGradient(u1,u1Grad);
			for (int i = 0; i < cellsize; i++)
				u2(i) =
				    COEFF(9)->GetNodalValue(sfcn->
							    GetNode(i));
			//sfcn->EvalGradient(u2,u2Grad);
			for (int i = 0; i < cellsize; i++)
				u3(i) =
				    COEFF(10)->GetNodalValue(sfcn->
							     GetNode(i));
			//sfcn->EvalGradient(u3,u3Grad);
			//for ( i=0; i < cellsize; i++ ) p(i) = COEFF(11) -> GetNodalValue( sfcn->GetNode(i) );    
			//sfcn->EvalGradient(p,pGrad);

			dt = COEFF(12)->Eval();
			//t = COEFF(13) -> Eval();
			Re = COEFF(14)->Eval();

			u1sum = u1(0) + u1(1) + u1(2) + u1(3);
			u2sum = u2(0) + u2(1) + u2(2) + u2(3);
			u3sum = u3(0) + u3(1) + u3(2) + u3(3);
			normU = sqrt(sqr(u1sum) + sqr(u2sum) + sqr(u3sum)) / 4;

			divUmid =
			    theta * (uprev_TS1Grad(0) + uprev_TS2Grad(1) +
				     uprev_TS3Grad(2)) + (1.0 -
							  theta) *
			    (u1Grad(0) + u2Grad(1) + u3Grad(2));
			divU = u1Grad(0) + u2Grad(1) + u3Grad(2);

			if (h * Re > 1.0) {
				//delta1 = 0.5 / sqrt( 1.0/sqr(dt) + sqr(normU/h) );
				delta1 = C1 * h;
			} else {
				delta1 = C1 * sqr(h);
			}
			if (h * Re > 1.0) {
				delta2 = C2 * h;
			} else {
				delta2 = C2 * sqr(h);
			}
		}

		if (tst_bf == 0) {
			uv1 = Iii;
			uv2 = uv3 = uv4 = Iij;
			vx = sfcn->dN1x();
			vy = sfcn->dN1y();
			vz = sfcn->dN1z();
		} else if (tst_bf == 1) {
			uv2 = Iii;
			uv1 = uv3 = uv4 = Iij;
			vx = sfcn->dN2x();
			vy = sfcn->dN2y();
			vz = sfcn->dN2z();
		} else if (tst_bf == 2) {
			uv3 = Iii;
			uv1 = uv2 = uv4 = Iij;
			vx = sfcn->dN3x();
			vy = sfcn->dN3y();
			vz = sfcn->dN3z();
		} else if (tst_bf == 3) {
			uv4 = Iii;
			uv1 = uv2 = uv3 = Iij;
			vx = sfcn->dN4x();
			vy = sfcn->dN4y();
			vz = sfcn->dN4z();
		}

		if (trial_bf == 0) {
			u1_1 = u2_1 = u3_1 = Iii;
			u1_2 = u1_3 = u1_4 = u2_2 = u2_3 = u2_4 = u3_2 =
			    u3_3 = u3_4 = Iij;
			D1u1 = D1u2 = D1u3 = sfcn->dN1x();
			D2u1 = D2u2 = D2u3 = sfcn->dN1y();
			D3u1 = D3u2 = D3u3 = sfcn->dN1z();
		} else if (trial_bf == 1) {
			u1_2 = u2_2 = u3_2 = Iii;
			u1_1 = u1_3 = u1_4 = u2_1 = u2_3 = u2_4 = u3_1 =
			    u3_3 = u3_4 = Iij;
			D1u1 = D1u2 = D1u3 = sfcn->dN2x();
			D2u1 = D2u2 = D2u3 = sfcn->dN2y();
			D3u1 = D3u2 = D3u3 = sfcn->dN2z();
		} else if (trial_bf == 2) {
			u1_3 = u2_3 = u3_3 = Iii;
			u1_1 = u1_2 = u1_4 = u2_1 = u2_2 = u2_4 = u3_1 =
			    u3_2 = u3_4 = Iij;
			D1u1 = D1u2 = D1u3 = sfcn->dN3x();
			D2u1 = D2u2 = D2u3 = sfcn->dN3y();
			D3u1 = D3u2 = D3u3 = sfcn->dN3z();
		} else if (trial_bf == 3) {
			u1_4 = u2_4 = u3_4 = Iii;
			u1_1 = u1_2 = u1_3 = u2_1 = u2_2 = u2_3 = u3_1 =
			    u3_2 = u3_3 = Iij;
			D1u1 = D1u2 = D1u3 = sfcn->dN4x();
			D2u1 = D2u2 = D2u3 = sfcn->dN4y();
			D3u1 = D3u2 = D3u3 = sfcn->dN4z();
		}

		if (tst_bf == trial_bf) {
			u1v = u2v = u3v = uv = Iii;
		} else {
			u1v = u2v = u3v = uv = Iij;
		}

		if (j_eq == 0) {
			u2_1 = u2_2 = u2_3 = u2_4 = u3_1 = u3_2 = u3_3 =
			    u3_4 = 0.0;;
			D1u2 = D2u2 = D3u2 = D1u3 = D2u3 = D3u3 = 0.0;
			u2v = u3v = 0.0;
		} else if (j_eq == 1) {
			u1_1 = u1_2 = u1_3 = u1_4 = u3_1 = u3_2 = u3_3 =
			    u3_4 = 0.0;;
			D1u1 = D2u1 = D3u1 = D1u3 = D2u3 = D3u3 = 0.0;
			u1v = u3v = 0.0;
		} else {
			u1_1 = u1_2 = u1_3 = u1_4 = u2_1 = u2_2 = u2_3 =
			    u2_4 = 0.0;;
			D1u1 = D2u1 = D3u1 = D1u2 = D2u2 = D3u2 = 0.0;
			u1v = u2v = 0.0;
		}

		if ((tst_bf == 0) && (trial_bf == 0)) {
			U1U1 =
			    U1(0) * U1(0) * Iii + U1(0) * U1(1) * Iij +
			    U1(0) * U1(2) * Iij + U1(0) * U1(3) * Iij +
			    U1(1) * U1(0) * Iij + U1(1) * U1(1) * Iii +
			    U1(1) * U1(2) * Iij + U1(1) * U1(3) * Iij +
			    U1(2) * U1(0) * Iij + U1(2) * U1(1) * Iij +
			    U1(2) * U1(2) * Iii + U1(2) * U1(3) * Iij +
			    U1(3) * U1(0) * Iij + U1(3) * U1(1) * Iij +
			    U1(3) * U1(2) * Iij + U1(3) * U1(3) * Iii;
			U1U2 =
			    U1(0) * U2(0) * Iii + U1(0) * U2(1) * Iij +
			    U1(0) * U2(2) * Iij + U1(0) * U2(3) * Iij +
			    U1(1) * U2(0) * Iij + U1(1) * U2(1) * Iii +
			    U1(1) * U2(2) * Iij + U1(1) * U2(3) * Iij +
			    U1(2) * U2(0) * Iij + U1(2) * U2(1) * Iij +
			    U1(2) * U2(2) * Iii + U1(2) * U2(3) * Iij +
			    U1(3) * U2(0) * Iij + U1(3) * U2(1) * Iij +
			    U1(3) * U2(2) * Iij + U1(3) * U2(3) * Iii;
			U1U3 =
			    U1(0) * U3(0) * Iii + U1(0) * U3(1) * Iij +
			    U1(0) * U3(2) * Iij + U1(0) * U3(3) * Iij +
			    U1(1) * U3(0) * Iij + U1(1) * U3(1) * Iii +
			    U1(1) * U3(2) * Iij + U1(1) * U3(3) * Iij +
			    U1(2) * U3(0) * Iij + U1(2) * U3(1) * Iij +
			    U1(2) * U3(2) * Iii + U1(2) * U3(3) * Iij +
			    U1(3) * U3(0) * Iij + U1(3) * U3(1) * Iij +
			    U1(3) * U3(2) * Iij + U1(3) * U3(3) * Iii;
			U2U1 = U1U2;
			U2U2 =
			    U2(0) * U2(0) * Iii + U2(0) * U2(1) * Iij +
			    U2(0) * U2(2) * Iij + U2(0) * U2(3) * Iij +
			    U2(1) * U2(0) * Iij + U2(1) * U2(1) * Iii +
			    U2(1) * U2(2) * Iij + U2(1) * U2(3) * Iij +
			    U2(2) * U2(0) * Iij + U2(2) * U2(1) * Iij +
			    U2(2) * U2(2) * Iii + U2(2) * U2(3) * Iij +
			    U2(3) * U2(0) * Iij + U2(3) * U2(1) * Iij +
			    U2(3) * U2(2) * Iij + U2(3) * U2(3) * Iii;
			U2U3 =
			    U2(0) * U3(0) * Iii + U2(0) * U3(1) * Iij +
			    U2(0) * U3(2) * Iij + U2(0) * U3(3) * Iij +
			    U2(1) * U3(0) * Iij + U2(1) * U3(1) * Iii +
			    U2(1) * U3(2) * Iij + U2(1) * U3(3) * Iij +
			    U2(2) * U3(0) * Iij + U2(2) * U3(1) * Iij +
			    U2(2) * U3(2) * Iii + U2(2) * U3(3) * Iij +
			    U2(3) * U3(0) * Iij + U2(3) * U3(1) * Iij +
			    U2(3) * U3(2) * Iij + U2(3) * U3(3) * Iii;
			U3U1 = U1U3;
			U3U2 = U2U3;
			U3U3 =
			    U3(0) * U3(0) * Iii + U3(0) * U3(1) * Iij +
			    U3(0) * U3(2) * Iij + U3(0) * U3(3) * Iij +
			    U3(1) * U3(0) * Iij + U3(1) * U3(1) * Iii +
			    U3(1) * U3(2) * Iij + U3(1) * U3(3) * Iij +
			    U3(2) * U3(0) * Iij + U3(2) * U3(1) * Iij +
			    U3(2) * U3(2) * Iii + U3(2) * U3(3) * Iij +
			    U3(3) * U3(0) * Iij + U3(3) * U3(1) * Iij +
			    U3(3) * U3(2) * Iij + U3(3) * U3(3) * Iii;
		}

		UDu_UDv_1 =
		    (1.0 - theta) * D1u1 * (vx * U1U1 + vy * U1U2 +
					    vz * U1U3)
		    + (1.0 - theta) * D2u1 * (vx * U2U1 + vy * U2U2 +
					      vz * U2U3)
		    + (1.0 - theta) * D3u1 * (vx * U3U1 + vy * U3U2 +
					      vz * U3U3);

		UDu_UDv_2 =
		    (1.0 - theta) * D1u2 * (vx * U1U1 + vy * U1U2 +
					    vz * U1U3)
		    + (1.0 - theta) * D2u2 * (vx * U2U1 + vy * U2U2 +
					      vz * U2U3)
		    + (1.0 - theta) * D3u2 * (vx * U3U1 + vy * U3U2 +
					      vz * U3U3);

		UDu_UDv_3 =
		    (1.0 - theta) * D1u3 * (vx * U1U1 + vy * U1U2 +
					    vz * U1U3)
		    + (1.0 - theta) * D2u3 * (vx * U2U1 + vy * U2U2 +
					      vz * U2U3)
		    + (1.0 - theta) * D3u3 * (vx * U3U1 + vy * U3U2 +
					      vz * U3U3);


		UDu_DUv_1 =
		    (1.0 -
		     theta) * (D1u1 * U1Grad(0) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u1 * U1Grad(0) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u1 * U1Grad(0) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4) +
			       D1u2 * U2Grad(0) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u2 * U2Grad(0) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u2 * U2Grad(0) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4) +
			       D1u3 * U3Grad(0) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u3 * U3Grad(0) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u3 * U3Grad(0) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4));

		UDu_DUv_2 =
		    (1.0 -
		     theta) * (D1u1 * U1Grad(1) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u1 * U1Grad(1) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u1 * U1Grad(1) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4) +
			       D1u2 * U2Grad(1) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u2 * U2Grad(1) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u2 * U2Grad(1) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4) +
			       D1u3 * U3Grad(1) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u3 * U3Grad(1) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u3 * U3Grad(1) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4));

		UDu_DUv_3 =
		    (1.0 -
		     theta) * (D1u1 * U1Grad(2) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u1 * U1Grad(2) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u1 * U1Grad(2) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4) +
			       D1u2 * U2Grad(2) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u2 * U2Grad(2) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u2 * U2Grad(2) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4) +
			       D1u3 * U3Grad(2) * (U1(0) * uv1 +
						   U1(1) * uv2 +
						   U1(2) * uv3 * U1(3) *
						   uv4) +
			       D2u3 * U3Grad(2) * (U2(0) * uv1 +
						   U2(1) * uv2 +
						   U2(2) * uv3 * U2(3) *
						   uv4) +
			       D3u3 * U3Grad(2) * (U3(0) * uv1 +
						   U3(1) * uv2 +
						   U3(2) * uv3 * U3(3) *
						   uv4));


		u1U1 =
		    u1_1 * U1(0) + u1_1 * U1(1) + u1_1 * U1(2) +
		    u1_1 * U1(3)
		    + u1_2 * U1(0) + u1_2 * U1(1) + u1_2 * U1(2) +
		    u1_2 * U1(3)
		    + u1_3 * U1(0) + u1_3 * U1(1) + u1_3 * U1(2) +
		    u1_3 * U1(3)
		    + u1_4 * U1(0) + u1_4 * U1(1) + u1_4 * U1(2) +
		    u1_4 * U1(3);
		u1U2 =
		    u1_1 * U2(0) + u1_1 * U2(1) + u1_1 * U2(2) +
		    u1_1 * U2(3)
		    + u1_2 * U2(0) + u1_2 * U2(1) + u1_2 * U2(2) +
		    u1_2 * U2(3)
		    + u1_3 * U2(0) + u1_3 * U2(1) + u1_3 * U2(2) +
		    u1_3 * U2(3)
		    + u1_4 * U2(0) + u1_4 * U2(1) + u1_4 * U2(2) +
		    u1_4 * U2(3);
		u1U3 =
		    u1_1 * U3(0) + u1_1 * U3(1) + u1_1 * U3(2) +
		    u1_1 * U3(3)
		    + u1_2 * U3(0) + u1_2 * U3(1) + u1_2 * U3(2) +
		    u1_2 * U3(3)
		    + u1_3 * U3(0) + u1_3 * U3(1) + u1_3 * U3(2) +
		    u1_3 * U3(3)
		    + u1_4 * U3(0) + u1_4 * U3(1) + u1_4 * U3(2) +
		    u1_4 * U3(3);
		u2U1 =
		    u2_1 * U1(0) + u2_1 * U1(1) + u2_1 * U1(2) +
		    u2_1 * U1(3)
		    + u2_2 * U1(0) + u2_2 * U1(1) + u2_2 * U1(2) +
		    u2_2 * U1(3)
		    + u2_3 * U1(0) + u2_3 * U1(1) + u2_3 * U1(2) +
		    u2_3 * U1(3)
		    + u2_4 * U1(0) + u2_4 * U1(1) + u2_4 * U1(2) +
		    u2_4 * U1(3);
		u2U2 =
		    u2_1 * U2(0) + u2_1 * U2(1) + u2_1 * U2(2) +
		    u2_1 * U2(3)
		    + u2_2 * U2(0) + u2_2 * U2(1) + u2_2 * U2(2) +
		    u2_2 * U2(3)
		    + u2_3 * U2(0) + u2_3 * U2(1) + u2_3 * U2(2) +
		    u2_3 * U2(3)
		    + u2_4 * U2(0) + u2_4 * U2(1) + u2_4 * U2(2) +
		    u2_4 * U2(3);
		u2U3 =
		    u2_1 * U3(0) + u2_1 * U3(1) + u2_1 * U3(2) +
		    u2_1 * U3(3)
		    + u2_2 * U3(0) + u2_2 * U3(1) + u2_2 * U3(2) +
		    u2_2 * U3(3)
		    + u2_3 * U3(0) + u2_3 * U3(1) + u2_3 * U3(2) +
		    u2_3 * U3(3)
		    + u2_4 * U3(0) + u2_4 * U3(1) + u2_4 * U3(2) +
		    u2_4 * U3(3);
		u3U1 =
		    u3_1 * U1(0) + u3_1 * U1(1) + u3_1 * U1(2) +
		    u3_1 * U1(3)
		    + u3_2 * U1(0) + u3_2 * U1(1) + u3_2 * U1(2) +
		    u3_2 * U1(3)
		    + u3_3 * U1(0) + u3_3 * U1(1) + u3_3 * U1(2) +
		    u3_3 * U1(3)
		    + u3_4 * U1(0) + u3_4 * U1(1) + u3_4 * U1(2) +
		    u3_4 * U1(3);
		u3U2 =
		    u3_1 * U2(0) + u3_1 * U2(1) + u3_1 * U2(2) +
		    u3_1 * U2(3)
		    + u3_2 * U2(0) + u3_2 * U2(1) + u3_2 * U2(2) +
		    u3_2 * U2(3)
		    + u3_3 * U2(0) + u3_3 * U2(1) + u3_3 * U2(2) +
		    u3_3 * U2(3)
		    + u3_4 * U2(0) + u3_4 * U2(1) + u3_4 * U2(2) +
		    u3_4 * U2(3);
		u3U3 =
		    u3_1 * U3(0) + u3_1 * U3(1) + u3_1 * U3(2) +
		    u3_1 * U3(3)
		    + u3_2 * U3(0) + u3_2 * U3(1) + u3_2 * U3(2) +
		    u3_2 * U3(3)
		    + u3_3 * U3(0) + u3_3 * U3(1) + u3_3 * U3(2) +
		    u3_3 * U3(3)
		    + u3_4 * U3(0) + u3_4 * U3(1) + u3_4 * U3(2) +
		    u3_4 * U3(3);

		DUu_DUv_1 =
		    (1.0 - theta) * (U1Grad(0) * U1Grad(0) * u1v +
				     U1Grad(0) * U2Grad(0) * u2v +
				     U1Grad(0) * U3Grad(0) * u3v +
				     U1Grad(1) * U1Grad(1) * u1v +
				     U1Grad(1) * U2Grad(1) * u2v +
				     U1Grad(1) * U3Grad(1) * u3v +
				     U1Grad(2) * U1Grad(2) * u1v +
				     U1Grad(2) * U2Grad(2) * u2v +
				     U1Grad(2) * U3Grad(2) * u3v);

		DUu_DUv_2 =
		    (1.0 - theta) * (U2Grad(0) * U1Grad(0) * u1v +
				     U2Grad(0) * U2Grad(0) * u2v +
				     U2Grad(0) * U3Grad(0) * u3v +
				     U2Grad(1) * U1Grad(1) * u1v +
				     U2Grad(1) * U2Grad(1) * u2v +
				     U2Grad(1) * U3Grad(1) * u3v +
				     U2Grad(2) * U1Grad(2) * u1v +
				     U2Grad(2) * U2Grad(2) * u2v +
				     U2Grad(2) * U3Grad(2) * u3v);

		DUu_DUv_3 =
		    (1.0 - theta) * (U3Grad(0) * U1Grad(0) * u1v +
				     U3Grad(0) * U2Grad(0) * u2v +
				     U3Grad(0) * U3Grad(0) * u3v +
				     U3Grad(1) * U1Grad(1) * u1v +
				     U3Grad(1) * U2Grad(1) * u2v +
				     U3Grad(1) * U3Grad(1) * u3v +
				     U3Grad(2) * U1Grad(2) * u1v +
				     U3Grad(2) * U2Grad(2) * u2v +
				     U3Grad(2) * U3Grad(2) * u3v);


		DUu_UDv_1 =
		    (1.0 -
		     theta) * (U1Grad(0) * (vx * u1U1 + vy * u1U2 +
					    vz * u1U3) +
			       U2Grad(0) * (vx * u2U1 + vy * u2U2 +
					    vz * u2U3) +
			       U3Grad(0) * (vx * u3U1 + vy * u3U2 +
					    vz * u3U3));

		DUu_UDv_2 =
		    (1.0 -
		     theta) * (U1Grad(1) * (vx * u1U1 + vy * u1U2 +
					    vz * u1U3) +
			       U2Grad(1) * (vx * u2U1 + vy * u2U2 +
					    vz * u2U3) +
			       U3Grad(1) * (vx * u3U1 + vy * u3U2 +
					    vz * u3U3));

		DUu_UDv_3 =
		    (1.0 -
		     theta) * (U1Grad(2) * (vx * u1U1 + vy * u1U2 +
					    vz * u1U3) +
			       U2Grad(2) * (vx * u2U1 + vy * u2U2 +
					    vz * u2U3) +
			       U3Grad(2) * (vx * u3U1 + vy * u3U2 +
					    vz * u3U3));

		MASS1 = u1v;
		MASS2 = u2v;
		MASS3 = u3v;

		REA1 =
		    (1.0 - theta) * (U1Grad(0) * MASS1 +
				     U2Grad(0) * MASS2 +
				     U3Grad(0) * MASS3);
		REA2 =
		    (1.0 - theta) * (U1Grad(1) * MASS1 +
				     U2Grad(1) * MASS2 +
				     U3Grad(1) * MASS3);
		REA3 =
		    (1.0 - theta) * (U1Grad(2) * MASS1 +
				     U2Grad(2) * MASS2 +
				     U3Grad(2) * MASS3);

		nu = (1.0 / Re);
		LAP1 =
		    (1.0 - theta) * nu * (D1u1 * vx + D2u1 * vy +
					  D3u1 * vz);
		LAP2 =
		    (1.0 - theta) * nu * (D1u2 * vx + D2u2 * vy +
					  D3u2 * vz);
		LAP3 =
		    (1.0 - theta) * nu * (D1u3 * vx + D2u3 * vy +
					  D3u3 * vz);

		NL1 =
		    (1.0 -
		     theta) * ((U1(0) * uv1 + U1(1) * uv2 + U1(2) * uv3 +
				U1(3) * uv4) * D1u1 + (U2(0) * uv1 +
						       U2(1) * uv2 +
						       U2(2) * uv3 +
						       U2(3) * uv4) *
			       D2u1 + (U3(0) * uv1 + U3(1) * uv2 +
				       U3(2) * uv3 + U3(3) * uv4) * D3u1);

		NL2 =
		    (1.0 -
		     theta) * ((U1(0) * uv1 + U1(1) * uv2 + U1(2) * uv3 +
				U1(3) * uv4) * D1u2 + (U2(0) * uv1 +
						       U2(1) * uv2 +
						       U2(2) * uv3 +
						       U2(3) * uv4) *
			       D2u2 + (U3(0) * uv1 + U3(1) * uv2 +
				       U3(2) * uv3 + U3(3) * uv4) * D3u2);

		NL3 =
		    (1.0 -
		     theta) * ((U1(0) * uv1 + U1(1) * uv2 + U1(2) * uv3 +
				U1(3) * uv4) * D1u3 + (U2(0) * uv1 +
						       U2(1) * uv2 +
						       U2(2) * uv3 +
						       U2(3) * uv4) *
			       D2u3 + (U3(0) * uv1 + U3(1) * uv2 +
				       U3(2) * uv3 + U3(3) * uv4) * D3u3);

		SD1 = UDu_UDv_1 + DUu_DUv_1 - UDu_DUv_1 - DUu_UDv_1;
		SD2 = UDu_UDv_2 + DUu_DUv_2 - UDu_DUv_2 - DUu_UDv_2;
		SD3 = UDu_UDv_3 + DUu_DUv_3 - UDu_DUv_3 - DUu_UDv_3;

		LSdiv1 = (1.0 - theta) * (D1u1 + D2u2 + D3u3) * vx;
		LSdiv2 = (1.0 - theta) * (D1u1 + D2u2 + D3u3) * vy;
		LSdiv3 = (1.0 - theta) * (D1u1 + D2u2 + D3u3) * vz;

		jac(0) =
		    -MASS1 / dt + LAP1 - NL1 + REA1 + delta1 * SD1 +
		    delta2 * LSdiv1;
		jac(1) =
		    -MASS2 / dt + LAP2 - NL2 + REA2 + delta1 * SD2 +
		    delta2 * LSdiv2;
		jac(2) =
		    -MASS3 / dt + LAP3 - NL3 + REA3 + delta1 * SD3 +
		    delta2 * LSdiv3;

	}

	//--------- Load definition for exact integration on tetrahedrons --------------
	inline
	    void residualExactInt(MV_ColMat < real > &bres, int &bf,
				  int &el, ShapeFunction *sfcn) {
	  
	  int cellsize = sfcn->GetNoNodes();

		if (bf == 0) {

			h = 2.0 * sfcn->GetCircumRadius();

			x(0) = sfcn->x1();
			x(1) = sfcn->x2();
			x(2) = sfcn->x3();
			x(3) = sfcn->x4();

			y(0) = sfcn->y1();
			y(1) = sfcn->y2();
			y(2) = sfcn->y3();
			y(3) = sfcn->y4();

			z(0) = sfcn->z1();
			z(1) = sfcn->z2();
			z(2) = sfcn->z3();
			z(3) = sfcn->z4();

			for (int i = 0; i < cellsize; i++)
				uprev_TS1(i) =
				    COEFF(0)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(uprev_TS1, uprev_TS1Grad);
			for (int i = 0; i < cellsize; i++)
				uprev_TS2(i) =
				    COEFF(1)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(uprev_TS2, uprev_TS2Grad);
			for (int i = 0; i < cellsize; i++)
				uprev_TS3(i) =
				    COEFF(2)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(uprev_TS3, uprev_TS3Grad);
			//for ( i=0; i < cellsize; i++ ) pprev_TS(i) = COEFF(3) -> GetNodalValue( sfcn->GetNode(i) );    
			//sfcn->EvalGradient(pprev_TS,pprev_TSGrad);

			for (int i = 0; i < cellsize; i++)
				U1(i) =
				    COEFF(4)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(U1, U1Grad);
			for (int i = 0; i < cellsize; i++)
				U2(i) =
				    COEFF(5)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(U2, U2Grad);
			for (int i = 0; i < cellsize; i++)
				U3(i) =
				    COEFF(6)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(U3, U3Grad);

			for (int i = 0; i < cellsize; i++)
				u1(i) =
				    COEFF(8)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(u1, u1Grad);
			for (int i = 0; i < cellsize; i++)
				u2(i) =
				    COEFF(9)->GetNodalValue(sfcn->
							    GetNode(i));
			sfcn->EvalGradient(u2, u2Grad);
			for (int i = 0; i < cellsize; i++)
				u3(i) =
				    COEFF(10)->GetNodalValue(sfcn->
							     GetNode(i));
			sfcn->EvalGradient(u3, u3Grad);
			for (int i = 0; i < cellsize; i++)
				p(i) =
				    COEFF(11)->GetNodalValue(sfcn->
							     GetNode(i));
			sfcn->EvalGradient(p, pGrad);

			dt = COEFF(12)->Eval();
			//t = COEFF(13) -> Eval();
			Re = COEFF(14)->Eval();


			u1sum = u1(0) + u1(1) + u1(2) + u1(3);
			u2sum = u2(0) + u2(1) + u2(2) + u2(3);
			u3sum = u3(0) + u3(1) + u3(2) + u3(3);
			normU =
			    sqrt(sqr(u1sum) + sqr(u2sum) + sqr(u3sum)) / 4;

			divUmid =
			    theta * (uprev_TS1Grad(0) + uprev_TS2Grad(1) +
				     uprev_TS3Grad(2)) + (1.0 -
							  theta) *
			    (u1Grad(0) + u2Grad(1) + u3Grad(2));
			divU = u1Grad(0) + u2Grad(1) + u3Grad(2);

			if (h * Re > 1.0) {
				delta1 = C1 * h;
			} else {
				delta1 = C1 * sqr(h);
			}
			if (divU > K) {
				r = 0.5 * divU;
			} else {
				r = 0.0;
			}
			if (h * Re > 1.0) {
				delta2 = C2 * h;
			} else {
				delta2 = C2 * sqr(h);
			}
			r = 0.0;

			psum = p(0) + p(1) + p(2) + p(3);
		}

		vx, vy, vz;
		uv1, uv2, uv3, uv4;
		if (bf == 0) {
			uv1 = Iii;
			uv2 = uv3 = uv4 = Iij;
			vx = sfcn->dN1x();
			vy = sfcn->dN1y();
			vz = sfcn->dN1z();
		} else if (bf == 1) {
			uv2 = Iii;
			uv1 = uv3 = uv4 = Iij;
			vx = sfcn->dN2x();
			vy = sfcn->dN2y();
			vz = sfcn->dN2z();
		} else if (bf == 2) {
			uv3 = Iii;
			uv1 = uv2 = uv4 = Iij;
			vx = sfcn->dN3x();
			vy = sfcn->dN3y();
			vz = sfcn->dN3z();
		} else if (bf == 3) {
			uv4 = Iii;
			uv1 = uv2 = uv3 = Iij;
			vx = sfcn->dN4x();
			vy = sfcn->dN4y();
			vz = sfcn->dN4z();
		}

		if (bf == 0) {
			U1U1 =
			    U1(0) * U1(0) * Iii + U1(0) * U1(1) * Iij +
			    U1(0) * U1(2) * Iij + U1(0) * U1(3) * Iij +
			    U1(1) * U1(0) * Iij + U1(1) * U1(1) * Iii +
			    U1(1) * U1(2) * Iij + U1(1) * U1(3) * Iij +
			    U1(2) * U1(0) * Iij + U1(2) * U1(1) * Iij +
			    U1(2) * U1(2) * Iii + U1(2) * U1(3) * Iij +
			    U1(3) * U1(0) * Iij + U1(3) * U1(1) * Iij +
			    U1(3) * U1(2) * Iij + U1(3) * U1(3) * Iii;
			U1U2 =
			    U1(0) * U2(0) * Iii + U1(0) * U2(1) * Iij +
			    U1(0) * U2(2) * Iij + U1(0) * U2(3) * Iij +
			    U1(1) * U2(0) * Iij + U1(1) * U2(1) * Iii +
			    U1(1) * U2(2) * Iij + U1(1) * U2(3) * Iij +
			    U1(2) * U2(0) * Iij + U1(2) * U2(1) * Iij +
			    U1(2) * U2(2) * Iii + U1(2) * U2(3) * Iij +
			    U1(3) * U2(0) * Iij + U1(3) * U2(1) * Iij +
			    U1(3) * U2(2) * Iij + U1(3) * U2(3) * Iii;
			U1U3 =
			    U1(0) * U3(0) * Iii + U1(0) * U3(1) * Iij +
			    U1(0) * U3(2) * Iij + U1(0) * U3(3) * Iij +
			    U1(1) * U3(0) * Iij + U1(1) * U3(1) * Iii +
			    U1(1) * U3(2) * Iij + U1(1) * U3(3) * Iij +
			    U1(2) * U3(0) * Iij + U1(2) * U3(1) * Iij +
			    U1(2) * U3(2) * Iii + U1(2) * U3(3) * Iij +
			    U1(3) * U3(0) * Iij + U1(3) * U3(1) * Iij +
			    U1(3) * U3(2) * Iij + U1(3) * U3(3) * Iii;
			U2U1 = U1U2;
			U2U2 =
			    U2(0) * U2(0) * Iii + U2(0) * U2(1) * Iij +
			    U2(0) * U2(2) * Iij + U2(0) * U2(3) * Iij +
			    U2(1) * U2(0) * Iij + U2(1) * U2(1) * Iii +
			    U2(1) * U2(2) * Iij + U2(1) * U2(3) * Iij +
			    U2(2) * U2(0) * Iij + U2(2) * U2(1) * Iij +
			    U2(2) * U2(2) * Iii + U2(2) * U2(3) * Iij +
			    U2(3) * U2(0) * Iij + U2(3) * U2(1) * Iij +
			    U2(3) * U2(2) * Iij + U2(3) * U2(3) * Iii;
			U2U3 =
			    U2(0) * U3(0) * Iii + U2(0) * U3(1) * Iij +
			    U2(0) * U3(2) * Iij + U2(0) * U3(3) * Iij +
			    U2(1) * U3(0) * Iij + U2(1) * U3(1) * Iii +
			    U2(1) * U3(2) * Iij + U2(1) * U3(3) * Iij +
			    U2(2) * U3(0) * Iij + U2(2) * U3(1) * Iij +
			    U2(2) * U3(2) * Iii + U2(2) * U3(3) * Iij +
			    U2(3) * U3(0) * Iij + U2(3) * U3(1) * Iij +
			    U2(3) * U3(2) * Iij + U2(3) * U3(3) * Iii;
			U3U1 = U1U3;
			U3U2 = U2U3;
			U3U3 =
			    U3(0) * U3(0) * Iii + U3(0) * U3(1) * Iij +
			    U3(0) * U3(2) * Iij + U3(0) * U3(3) * Iij +
			    U3(1) * U3(0) * Iij + U3(1) * U3(1) * Iii +
			    U3(1) * U3(2) * Iij + U3(1) * U3(3) * Iij +
			    U3(2) * U3(0) * Iij + U3(2) * U3(1) * Iij +
			    U3(2) * U3(2) * Iii + U3(2) * U3(3) * Iij +
			    U3(3) * U3(0) * Iij + U3(3) * U3(1) * Iij +
			    U3(3) * U3(2) * Iij + U3(3) * U3(3) * Iii;
		}

		UDu_UDv_1 =
		    theta * uprev_TS1Grad(0) * (vx * U1U1 + vy * U1U2 +
						vz * U1U3) +
		    theta * uprev_TS1Grad(1) * (vx * U2U1 + vy * U2U2 +
						vz * U2U3) +
		    theta * uprev_TS1Grad(2) * (vx * U3U1 + vy * U3U2 +
						vz * U3U3);

		UDu_UDv_2 =
		    theta * uprev_TS2Grad(0) * (vx * U1U1 + vy * U1U2 +
						vz * U1U3) +
		    theta * uprev_TS2Grad(1) * (vx * U2U1 + vy * U2U2 +
						vz * U2U3) +
		    theta * uprev_TS2Grad(2) * (vx * U3U1 + vy * U3U2 +
						vz * U3U3);

		UDu_UDv_3 =
		    theta * uprev_TS3Grad(0) * (vx * U1U1 + vy * U1U2 +
						vz * U1U3) +
		    theta * uprev_TS3Grad(1) * (vx * U2U1 + vy * U2U2 +
						vz * U2U3) +
		    theta * uprev_TS3Grad(2) * (vx * U3U1 + vy * U3U2 +
						vz * U3U3);

		LSdiv1 =
		    theta * (uprev_TS1Grad(0) + uprev_TS2Grad(1) +
			     uprev_TS3Grad(2)) * vx;
		LSdiv2 =
		    theta * (uprev_TS1Grad(0) + uprev_TS2Grad(1) +
			     uprev_TS3Grad(2)) * vy;
		LSdiv3 =
		    theta * (uprev_TS1Grad(0) + uprev_TS2Grad(1) +
			     uprev_TS3Grad(2)) * vz;

		Dp_DUv_1 =
		    (pGrad(0) * U1Grad(0) + pGrad(1) * U2Grad(0) +
		     pGrad(2) * U3Grad(0)) * Ii;
		Dp_DUv_2 =
		    (pGrad(0) * U1Grad(1) + pGrad(1) * U2Grad(1) +
		     pGrad(2) * U3Grad(1)) * Ii;
		Dp_DUv_3 =
		    (pGrad(0) * U1Grad(2) + pGrad(1) * U2Grad(2) +
		     pGrad(2) * U3Grad(2)) * Ii;

		UDv =
		    ((U1(0) + U1(1) + U1(2) + U1(3)) * vx +
		     (U2(0) + U2(1) + U2(2) + U2(3)) * vy + (U3(0) +
							     U3(1) +
							     U3(2) +
							     U3(3)) * vz) *
		    Ii;

		UDu_DUv_1 =
		    theta * (uprev_TS1Grad(0) * U1Grad(0) *
			     (U1(0) * uv1 + U1(1) * uv2 +
			      U1(2) * uv3 * U1(3) * uv4) +
			     uprev_TS1Grad(1) * U1Grad(0) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS1Grad(2) * U1Grad(0) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4) +
			     uprev_TS2Grad(0) * U2Grad(0) * (U1(0) * uv1 +
							     U1(1) * uv2 +
							     U1(2) * uv3 *
							     U1(3) * uv4) +
			     uprev_TS2Grad(1) * U2Grad(0) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS2Grad(2) * U2Grad(0) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4) +
			     uprev_TS3Grad(0) * U3Grad(0) * (U1(0) * uv1 +
							     U1(1) * uv2 +
							     U1(2) * uv3 *
							     U1(3) * uv4) +
			     uprev_TS3Grad(1) * U3Grad(0) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS3Grad(2) * U3Grad(0) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4));

		UDu_DUv_2 =
		    theta * (uprev_TS1Grad(0) * U1Grad(1) *
			     (U1(0) * uv1 + U1(1) * uv2 +
			      U1(2) * uv3 * U1(3) * uv4) +
			     uprev_TS1Grad(1) * U1Grad(1) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS1Grad(2) * U1Grad(1) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4) +
			     uprev_TS2Grad(0) * U2Grad(1) * (U1(0) * uv1 +
							     U1(1) * uv2 +
							     U1(2) * uv3 *
							     U1(3) * uv4) +
			     uprev_TS2Grad(1) * U2Grad(1) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS2Grad(2) * U2Grad(1) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4) +
			     uprev_TS3Grad(0) * U3Grad(1) * (U1(0) * uv1 +
							     U1(1) * uv2 +
							     U1(2) * uv3 *
							     U1(3) * uv4) +
			     uprev_TS3Grad(1) * U3Grad(1) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS3Grad(2) * U3Grad(1) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4));

		UDu_DUv_3 =
		    theta * (uprev_TS1Grad(0) * U1Grad(2) *
			     (U1(0) * uv1 + U1(1) * uv2 +
			      U1(2) * uv3 * U1(3) * uv4) +
			     uprev_TS1Grad(1) * U1Grad(2) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS1Grad(2) * U1Grad(2) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4) +
			     uprev_TS2Grad(0) * U2Grad(2) * (U1(0) * uv1 +
							     U1(1) * uv2 +
							     U1(2) * uv3 *
							     U1(3) * uv4) +
			     uprev_TS2Grad(1) * U2Grad(2) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS2Grad(2) * U2Grad(2) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4) +
			     uprev_TS3Grad(0) * U3Grad(2) * (U1(0) * uv1 +
							     U1(1) * uv2 +
							     U1(2) * uv3 *
							     U1(3) * uv4) +
			     uprev_TS3Grad(1) * U3Grad(2) * (U2(0) * uv1 +
							     U2(1) * uv2 +
							     U2(2) * uv3 *
							     U2(3) * uv4) +
			     uprev_TS3Grad(2) * U3Grad(2) * (U3(0) * uv1 +
							     U3(1) * uv2 +
							     U3(2) * uv3 *
							     U3(3) * uv4));


		if (bf == 0) {
			u1U1 =
			    uprev_TS1(0) * U1(0) * Iii +
			    uprev_TS1(0) * U1(1) * Iij +
			    uprev_TS1(0) * U1(2) * Iij +
			    uprev_TS1(0) * U1(3) * Iij +
			    uprev_TS1(1) * U1(0) * Iij +
			    uprev_TS1(1) * U1(1) * Iii +
			    uprev_TS1(1) * U1(2) * Iij +
			    uprev_TS1(1) * U1(3) * Iij +
			    uprev_TS1(2) * U1(0) * Iij +
			    uprev_TS1(2) * U1(1) * Iij +
			    uprev_TS1(2) * U1(2) * Iii +
			    uprev_TS1(2) * U1(3) * Iij +
			    uprev_TS1(3) * U1(0) * Iij +
			    uprev_TS1(3) * U1(1) * Iij +
			    uprev_TS1(3) * U1(2) * Iij +
			    uprev_TS1(3) * U1(3) * Iii;
			u1U2 =
			    uprev_TS1(0) * U2(0) * Iii +
			    uprev_TS1(0) * U2(1) * Iij +
			    uprev_TS1(0) * U2(2) * Iij +
			    uprev_TS1(0) * U2(3) * Iij +
			    uprev_TS1(1) * U2(0) * Iij +
			    uprev_TS1(1) * U2(1) * Iii +
			    uprev_TS1(1) * U2(2) * Iij +
			    uprev_TS1(1) * U2(3) * Iij +
			    uprev_TS1(2) * U2(0) * Iij +
			    uprev_TS1(2) * U2(1) * Iij +
			    uprev_TS1(2) * U2(2) * Iii +
			    uprev_TS1(2) * U2(3) * Iij +
			    uprev_TS1(3) * U2(0) * Iij +
			    uprev_TS1(3) * U2(1) * Iij +
			    uprev_TS1(3) * U2(2) * Iij +
			    uprev_TS1(3) * U2(3) * Iii;
			u1U3 =
			    uprev_TS1(0) * U3(0) * Iii +
			    uprev_TS1(0) * U3(1) * Iij +
			    uprev_TS1(0) * U3(2) * Iij +
			    uprev_TS1(0) * U3(3) * Iij +
			    uprev_TS1(1) * U3(0) * Iij +
			    uprev_TS1(1) * U3(1) * Iii +
			    uprev_TS1(1) * U3(2) * Iij +
			    uprev_TS1(1) * U3(3) * Iij +
			    uprev_TS1(2) * U3(0) * Iij +
			    uprev_TS1(2) * U3(1) * Iij +
			    uprev_TS1(2) * U3(2) * Iii +
			    uprev_TS1(2) * U3(3) * Iij +
			    uprev_TS1(3) * U3(0) * Iij +
			    uprev_TS1(3) * U3(1) * Iij +
			    uprev_TS1(3) * U3(2) * Iij +
			    uprev_TS1(3) * U3(3) * Iii;

			u2U1 =
			    uprev_TS2(0) * U1(0) * Iii +
			    uprev_TS2(0) * U1(1) * Iij +
			    uprev_TS2(0) * U1(2) * Iij +
			    uprev_TS2(0) * U1(3) * Iij +
			    uprev_TS2(1) * U1(0) * Iij +
			    uprev_TS2(1) * U1(1) * Iii +
			    uprev_TS2(1) * U1(2) * Iij +
			    uprev_TS2(1) * U1(3) * Iij +
			    uprev_TS2(2) * U1(0) * Iij +
			    uprev_TS2(2) * U1(1) * Iij +
			    uprev_TS2(2) * U1(2) * Iii +
			    uprev_TS2(2) * U1(3) * Iij +
			    uprev_TS2(3) * U1(0) * Iij +
			    uprev_TS2(3) * U1(1) * Iij +
			    uprev_TS2(3) * U1(2) * Iij +
			    uprev_TS2(3) * U1(3) * Iii;
			u2U2 =
			    uprev_TS2(0) * U2(0) * Iii +
			    uprev_TS2(0) * U2(1) * Iij +
			    uprev_TS2(0) * U2(2) * Iij +
			    uprev_TS2(0) * U2(3) * Iij +
			    uprev_TS2(1) * U2(0) * Iij +
			    uprev_TS2(1) * U2(1) * Iii +
			    uprev_TS2(1) * U2(2) * Iij +
			    uprev_TS2(1) * U2(3) * Iij +
			    uprev_TS2(2) * U2(0) * Iij +
			    uprev_TS2(2) * U2(1) * Iij +
			    uprev_TS2(2) * U2(2) * Iii +
			    uprev_TS2(2) * U2(3) * Iij +
			    uprev_TS2(3) * U2(0) * Iij +
			    uprev_TS2(3) * U2(1) * Iij +
			    uprev_TS2(3) * U2(2) * Iij +
			    uprev_TS2(3) * U2(3) * Iii;
			u2U3 =
			    uprev_TS2(0) * U3(0) * Iii +
			    uprev_TS2(0) * U3(1) * Iij +
			    uprev_TS2(0) * U3(2) * Iij +
			    uprev_TS2(0) * U3(3) * Iij +
			    uprev_TS2(1) * U3(0) * Iij +
			    uprev_TS2(1) * U3(1) * Iii +
			    uprev_TS2(1) * U3(2) * Iij +
			    uprev_TS2(1) * U3(3) * Iij +
			    uprev_TS2(2) * U3(0) * Iij +
			    uprev_TS2(2) * U3(1) * Iij +
			    uprev_TS2(2) * U3(2) * Iii +
			    uprev_TS2(2) * U3(3) * Iij +
			    uprev_TS2(3) * U3(0) * Iij +
			    uprev_TS2(3) * U3(1) * Iij +
			    uprev_TS2(3) * U3(2) * Iij +
			    uprev_TS2(3) * U3(3) * Iii;

			u3U1 =
			    uprev_TS3(0) * U1(0) * Iii +
			    uprev_TS3(0) * U1(1) * Iij +
			    uprev_TS3(0) * U1(2) * Iij +
			    uprev_TS3(0) * U1(3) * Iij +
			    uprev_TS3(1) * U1(0) * Iij +
			    uprev_TS3(1) * U1(1) * Iii +
			    uprev_TS3(1) * U1(2) * Iij +
			    uprev_TS3(1) * U1(3) * Iij +
			    uprev_TS3(2) * U1(0) * Iij +
			    uprev_TS3(2) * U1(1) * Iij +
			    uprev_TS3(2) * U1(2) * Iii +
			    uprev_TS3(2) * U1(3) * Iij +
			    uprev_TS3(3) * U1(0) * Iij +
			    uprev_TS3(3) * U1(1) * Iij +
			    uprev_TS3(3) * U1(2) * Iij +
			    uprev_TS3(3) * U1(3) * Iii;
			u3U2 =
			    uprev_TS3(0) * U2(0) * Iii +
			    uprev_TS3(0) * U2(1) * Iij +
			    uprev_TS3(0) * U2(2) * Iij +
			    uprev_TS3(0) * U2(3) * Iij +
			    uprev_TS3(1) * U2(0) * Iij +
			    uprev_TS3(1) * U2(1) * Iii +
			    uprev_TS3(1) * U2(2) * Iij +
			    uprev_TS3(1) * U2(3) * Iij +
			    uprev_TS3(2) * U2(0) * Iij +
			    uprev_TS3(2) * U2(1) * Iij +
			    uprev_TS3(2) * U2(2) * Iii +
			    uprev_TS3(2) * U2(3) * Iij +
			    uprev_TS3(3) * U2(0) * Iij +
			    uprev_TS3(3) * U2(1) * Iij +
			    uprev_TS3(3) * U2(2) * Iij +
			    uprev_TS3(3) * U2(3) * Iii;
			u3U3 =
			    uprev_TS3(0) * U3(0) * Iii +
			    uprev_TS3(0) * U3(1) * Iij +
			    uprev_TS3(0) * U3(2) * Iij +
			    uprev_TS3(0) * U3(3) * Iij +
			    uprev_TS3(1) * U3(0) * Iij +
			    uprev_TS3(1) * U3(1) * Iii +
			    uprev_TS3(1) * U3(2) * Iij +
			    uprev_TS3(1) * U3(3) * Iij +
			    uprev_TS3(2) * U3(0) * Iij +
			    uprev_TS3(2) * U3(1) * Iij +
			    uprev_TS3(2) * U3(2) * Iii +
			    uprev_TS3(2) * U3(3) * Iij +
			    uprev_TS3(3) * U3(0) * Iij +
			    uprev_TS3(3) * U3(1) * Iij +
			    uprev_TS3(3) * U3(2) * Iij +
			    uprev_TS3(3) * U3(3) * Iii;
		}

		DUu_UDv_1 =
		    theta * (U1Grad(0) *
			     (vx * u1U1 + vy * u1U2 + vz * u1U3)
			     + U2Grad(0) * (vx * u2U1 + vy * u2U2 +
					    vz * u2U3)
			     + U3Grad(0) * (vx * u3U1 + vy * u3U2 +
					    vz * u3U3));

		DUu_UDv_2 =
		    theta * (U1Grad(1) *
			     (vx * u1U1 + vy * u1U2 + vz * u1U3)
			     + U2Grad(1) * (vx * u2U1 + vy * u2U2 +
					    vz * u2U3)
			     + U3Grad(1) * (vx * u3U1 + vy * u3U2 +
					    vz * u3U3));

		DUu_UDv_3 =
		    theta * (U1Grad(2) *
			     (vx * u1U1 + vy * u1U2 + vz * u1U3)
			     + U2Grad(2) * (vx * u2U1 + vy * u2U2 +
					    vz * u2U3)
			     + U3Grad(2) * (vx * u3U1 + vy * u3U2 +
					    vz * u3U3));


		DUu_DUv_1 =
		    theta * (U1Grad(0) * U1Grad(0) *
			     (uprev_TS1(0) * uv1 + uprev_TS1(1) * uv2 +
			      uprev_TS1(2) * uv3 * uprev_TS1(3) * uv4) +
			     U1Grad(0) * U2Grad(0) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U1Grad(0) * U3Grad(0) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4) +
			     U1Grad(1) * U1Grad(1) * (uprev_TS1(0) * uv1 +
						      uprev_TS1(1) * uv2 +
						      uprev_TS1(2) * uv3 *
						      uprev_TS1(3) * uv4) +
			     U1Grad(1) * U2Grad(1) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U1Grad(1) * U3Grad(1) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4) +
			     U1Grad(2) * U1Grad(2) * (uprev_TS1(0) * uv1 +
						      uprev_TS1(1) * uv2 +
						      uprev_TS1(2) * uv3 *
						      uprev_TS1(3) * uv4) +
			     U1Grad(2) * U2Grad(2) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U1Grad(2) * U3Grad(2) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4));

		DUu_DUv_2 =
		    theta * (U2Grad(0) * U1Grad(0) *
			     (uprev_TS1(0) * uv1 + uprev_TS1(1) * uv2 +
			      uprev_TS1(2) * uv3 * uprev_TS1(3) * uv4) +
			     U2Grad(0) * U2Grad(0) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U2Grad(0) * U3Grad(0) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4) +
			     U2Grad(1) * U1Grad(1) * (uprev_TS1(0) * uv1 +
						      uprev_TS1(1) * uv2 +
						      uprev_TS1(2) * uv3 *
						      uprev_TS1(3) * uv4) +
			     U2Grad(1) * U2Grad(1) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U2Grad(1) * U3Grad(1) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4) +
			     U2Grad(2) * U1Grad(2) * (uprev_TS1(0) * uv1 +
						      uprev_TS1(1) * uv2 +
						      uprev_TS1(2) * uv3 *
						      uprev_TS1(3) * uv4) +
			     U2Grad(2) * U2Grad(2) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U2Grad(2) * U3Grad(2) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4));

		DUu_DUv_3 =
		    theta * (U3Grad(0) * U1Grad(0) *
			     (uprev_TS1(0) * uv1 + uprev_TS1(1) * uv2 +
			      uprev_TS1(2) * uv3 * uprev_TS1(3) * uv4) +
			     U3Grad(0) * U2Grad(0) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U3Grad(0) * U3Grad(0) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4) +
			     U3Grad(1) * U1Grad(1) * (uprev_TS1(0) * uv1 +
						      uprev_TS1(1) * uv2 +
						      uprev_TS1(2) * uv3 *
						      uprev_TS1(3) * uv4) +
			     U3Grad(1) * U2Grad(1) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U3Grad(1) * U3Grad(1) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4) +
			     U3Grad(2) * U1Grad(2) * (uprev_TS1(0) * uv1 +
						      uprev_TS1(1) * uv2 +
						      uprev_TS1(2) * uv3 *
						      uprev_TS1(3) * uv4) +
			     U3Grad(2) * U2Grad(2) * (uprev_TS2(0) * uv1 +
						      uprev_TS2(1) * uv2 +
						      uprev_TS2(2) * uv3 *
						      uprev_TS2(3) * uv4) +
			     U3Grad(2) * U3Grad(2) * (uprev_TS3(0) * uv1 +
						      uprev_TS3(1) * uv2 +
						      uprev_TS3(2) * uv3 *
						      uprev_TS3(3) * uv4));



		f1 = 0.0;
		f2 = 0.0;
		f3 = 0.0;

		d = 0.125;
		tf = 10.0 - d;
		//tf = 9.0;

		centerpt_x = 1.5;
		centerpt_y = 0.5;
		centerpt_z = 0.5;
		x_el = 0.25 * (x(0) + x(1) + x(2) + x(3));
		y_el = 0.25 * (y(0) + y(1) + y(2) + y(3));
		z_el = 0.25 * (z(0) + z(1) + z(2) + z(3));
		if ((t >= tf) && ((centerpt_x - 0.5 * d) <= x_el)
		    && ((centerpt_x + 0.5 * d) >= x_el)
		    && ((centerpt_y - 0.5 * d) <= y_el)
		    && ((centerpt_y + 0.5 * d) >= y_el)
		    && ((centerpt_z - 0.5 * d) <= z_el)
		    && ((centerpt_z + 0.5 * d) >= z_el)) {
			f1 = 1.0;
			f2 = 0.0;
			f3 = 0.0;
		}

		f_DUv_1 =
		    (f1 * U1Grad(0) + f2 * U2Grad(0) +
		     f3 * U3Grad(0)) * Ii;
		f_DUv_2 =
		    (f1 * U1Grad(1) + f2 * U2Grad(1) +
		     f3 * U3Grad(1)) * Ii;
		f_DUv_3 =
		    (f1 * U1Grad(2) + f2 * U2Grad(2) +
		     f3 * U3Grad(2)) * Ii;

		f_UDv_1 = f1 * vx * (U1(0) + U1(1) + U1(2) + U1(3)) * Ii
		    + f2 * vy * (U2(0) + U2(1) + U2(2) + U2(3)) * Ii
		    + f3 * vz * (U3(0) + U3(1) + U3(2) + U3(3)) * Ii;

		f_UDv_2 = f_UDv_1;
		f_UDv_3 = f_UDv_1;

		f_v_1 = f1 * Ii;
		f_v_2 = f2 * Ii;
		f_v_3 = f3 * Ii;


		NL1 =
		    theta *
		    ((U1(0) * uv1 + U1(1) * uv2 + U1(2) * uv3 +
		      U1(3) * uv4) * uprev_TS1Grad(0) + (U2(0) * uv1 +
							 U2(1) * uv2 +
							 U2(2) * uv3 +
							 U2(3) * uv4) *
		     uprev_TS1Grad(1) + (U3(0) * uv1 + U3(1) * uv2 +
					 U3(2) * uv3 +
					 U3(3) * uv4) * uprev_TS1Grad(2));

		NL2 =
		    theta *
		    ((U1(0) * uv1 + U1(1) * uv2 + U1(2) * uv3 +
		      U1(3) * uv4) * uprev_TS2Grad(0) + (U2(0) * uv1 +
							 U2(1) * uv2 +
							 U2(2) * uv3 +
							 U2(3) * uv4) *
		     uprev_TS2Grad(1) + (U3(0) * uv1 + U3(1) * uv2 +
					 U3(2) * uv3 +
					 U3(3) * uv4) * uprev_TS2Grad(2));

		NL3 =
		    theta *
		    ((U1(0) * uv1 + U1(1) * uv2 + U1(2) * uv3 +
		      U1(3) * uv4) * uprev_TS3Grad(0) + (U2(0) * uv1 +
							 U2(1) * uv2 +
							 U2(2) * uv3 +
							 U2(3) * uv4) *
		     uprev_TS3Grad(1) + (U3(0) * uv1 + U3(1) * uv2 +
					 U3(2) * uv3 +
					 U3(3) * uv4) * uprev_TS3Grad(2));

		nu = (1.0 / Re);
		LAP1 =
		    theta * nu * (uprev_TS1Grad(0) * vx +
				  uprev_TS1Grad(1) * vy +
				  uprev_TS1Grad(2) * vz);
		LAP2 =
		    theta * nu * (uprev_TS2Grad(0) * vx +
				  uprev_TS2Grad(1) * vy +
				  uprev_TS2Grad(2) * vz);
		LAP3 =
		    theta * nu * (uprev_TS3Grad(0) * vx +
				  uprev_TS3Grad(1) * vy +
				  uprev_TS3Grad(2) * vz);

		u1_v_TS =
		    uprev_TS1(0) * uv1 + uprev_TS1(1) * uv2 +
		    uprev_TS1(2) * uv3 + uprev_TS1(3) * uv4;
		u2_v_TS =
		    uprev_TS2(0) * uv1 + uprev_TS2(1) * uv2 +
		    uprev_TS2(2) * uv3 + uprev_TS2(3) * uv4;
		u3_v_TS =
		    uprev_TS3(0) * uv1 + uprev_TS3(1) * uv2 +
		    uprev_TS3(2) * uv3 + uprev_TS3(3) * uv4;

		REA1 =
		    theta * (U1Grad(0) * u1_v_TS + U2Grad(0) * u2_v_TS +
			     U3Grad(0) * u3_v_TS);
		REA2 =
		    theta * (U1Grad(1) * u1_v_TS + U2Grad(1) * u2_v_TS +
			     U3Grad(1) * u3_v_TS);
		REA3 =
		    theta * (U1Grad(2) * u1_v_TS + U2Grad(2) * u2_v_TS +
			     U3Grad(2) * u3_v_TS);

		SD1 =
		    -UDu_UDv_1 - DUu_DUv_1 + UDu_DUv_1 + DUu_UDv_1 -
		    Dp_DUv_1 + pGrad(0) * UDv - f_UDv_1 + f_DUv_1;
		SD2 =
		    -UDu_UDv_2 - DUu_DUv_2 + UDu_DUv_2 + DUu_UDv_2 -
		    Dp_DUv_2 + pGrad(1) * UDv - f_UDv_2 + f_DUv_2;
		SD3 =
		    -UDu_UDv_3 - DUu_DUv_3 + UDu_DUv_3 + DUu_UDv_3 -
		    Dp_DUv_3 + pGrad(2) * UDv - f_UDv_3 + f_DUv_3;

		Dp_v1 = psum * vx * Ii;
		Dp_v2 = psum * vy * Ii;
		Dp_v3 = psum * vz * Ii;

		bres(0, bf) =
		    -u1_v_TS / dt + f_v_1 + NL1 - REA1 - LAP1 + Dp_v1 +
		    delta1 * SD1 - delta2 * LSdiv1;
		bres(1, bf) =
		    -u2_v_TS / dt + f_v_2 + NL2 - REA2 - LAP2 + Dp_v2 +
		    delta1 * SD2 - delta2 * LSdiv2;;
		bres(2, bf) =
		    -u3_v_TS / dt + f_v_3 + NL3 - REA3 - LAP3 + Dp_v3 +
		    delta1 * SD3 - delta2 * LSdiv3;

	}


	/*
	   inline
	   void MatrixVectorProductExactInt( MV_ColMat<real>& bres, int& bf, int& el, FET4n3D& FE ){}

	   inline
	   void NonLinearResidualExactInt( MV_ColMat<real>& bres, int& bf, int& el, FET4n3D& FE ){}
	 */

};

#endif
