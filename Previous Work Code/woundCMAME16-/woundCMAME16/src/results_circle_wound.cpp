/*

RESULTS
circular wound problem.

Read a quad mesh defined by myself.
Then apply boundary conditions.
Solve.

*/

#include "wound.h"
#include "solver.h"
#include "myMeshGenerator.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept> 
#include <math.h> 
#include <string>
#include <time.h>

#include <Eigen/Dense> 
using namespace Eigen;

double frand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


int main(int argc, char *argv[])
{

	std::cout<<"\nResults 3: full domain simulations \n";
	srand (time(NULL));
	
	//---------------------------------//
	// GLOBAL PARAMETERS
	//
	// for normalization
	double rho_phys = 10000.; // [cells/mm^3]
	double c_max = 1.0e-4; // [g/mm3] from tgf beta review, 5e-5g/mm3 was good for tissues
	//
	double k0 = 0.0511; // neo hookean for skin, used previously, in MPa
	double kf = 0.015; // stiffness of collagen in MPa, from previous paper
	double k2 = 0.048; // nonlinear exponential coefficient, non-dimensional
	double t_rho = 0.0045/rho_phys; // force of fibroblasts in MPa, this is per cell. so, in an average sense this is the production by the natural density
	double t_rho_c = (0.045)/rho_phys; // force of myofibroblasts enhanced by chemical, I'm assuming normalized chemical, otherwise I'd have to add a normalizing constant
	double K_t_c = c_max/10.; // saturation of chemical on force. this can be calculated from steady state
	double D_rhorho = 0.0833; // diffusion of cells in [mm^2/hour], not normalized
	double D_rhoc = 1.66e-12/c_max/c_max; // diffusion of chemotactic gradient, an order of magnitude greater than random walk [mm^2/hour], not normalized
	double D_cc = 0.10; // diffusion of chemical TGF, not normalized. 
	double p_rho =0.034; // in 1/hour production of fibroblasts naturally, proliferation rate, not normalized, based on data of doubling rate from commercial use
	double p_rho_c = 0.034; // production enhanced by the chem, if the chemical is normalized, then suggest two fold,
	double p_rho_theta = 0.034; // enhanced production by theta
	double K_rho_c= c_max/10.; // saturation of cell proliferation by chemical, this one is definitely not crucial, just has to be small enough <cmax
	double d_rho = 0.2*p_rho ;// decay of cells, 20 percent of cells die per day
	double K_rho_rho = rho_phys; // saturation of cell by cell, from steady state
	double vartheta_e = 2.0; // physiological state of area stretch
	double gamma_theta = 5.; // sensitivity of heviside function
	double p_c_rho = 90.0e-16/rho_phys;// production of c by cells in g/cells/h
	double p_c_thetaE = 300.0e-16/rho_phys; // coupling of elastic and chemical, three fold
	double K_c_c = 1.;// saturation of chem by chem, from steady state
	double d_c = 0.01; // decay of chemical in 1/hours
	//---------------------------------//
	std::vector<double> global_parameters = {k0,kf,k2,t_rho,t_rho_c,K_t_c,D_rhorho,D_rhoc,D_cc,p_rho,p_rho_c,p_rho_theta,K_rho_c,K_rho_rho,d_rho,vartheta_e,gamma_theta,p_c_rho,p_c_thetaE,K_c_c,d_c};

	//---------------------------------//
	// LOCAL PARAMETERS
	//
	// collagen fraction
	double p_phi = 0.002/rho_phys; // production by fibroblasts, natural rate in percent/hour, 5% per day
	double p_phi_c = p_phi; // production up-regulation, weighted by C and rho
	double p_phi_theta = p_phi; // mechanosensing upregulation. no need to normalize by Hmax since Hmax = 1
	double K_phi_c = 0.0001; // saturation of C effect on deposition. RANDOM?
	double d_phi = 0.000970; // rate of degradation, in the order of the wound process, 100 percent in one year for wound, means 0.000116 effective per hour means degradation = 0.002 - 0.000116
	double d_phi_rho_c = 0.5*0.000970/rho_phys/c_max;//0.000194; // degradation coupled to chemical and cell density to maintain phi equilibrium
	double K_phi_rho = rho_phys*p_phi/d_phi -1; // saturation of collagen fraction itself, from steady state
	//
	//
	// fiber alignment
	double tau_omega = 10./(K_phi_rho+1); // time constant for angular reorientation, think 100 percent in one year
	//
	// dispersion parameter
	double tau_kappa = 1./(K_phi_rho+1); // time constant, on the order of a year
	double gamma_kappa = 5.; // exponent of the principal stretch ratio
	// 
	// permanent contracture/growth
	double tau_lamdaP_a = 1.0/(K_phi_rho+1); // time constant for direction a, on the order of a year
	double tau_lamdaP_s = 1.0/(K_phi_rho+1); // time constant for direction s, on the order of a year
	//
	std::vector<double> local_parameters = {p_phi,p_phi_c,p_phi_theta,K_phi_c,K_phi_rho,d_phi,d_phi_rho_c,tau_omega,tau_kappa,gamma_kappa,tau_lamdaP_a,tau_lamdaP_s,gamma_theta,vartheta_e};
	//
	// solution parameters
	double tol_local = 1e-5; // local tolerance
	int max_local_iter = 35; // max local iter
	//
	// other local stuff
	double PIE = 3.14159;
	//---------------------------------//
	
	
	//---------------------------------//
	// values for the wound
	double rho_wound = 0; // [cells/mm^3]
	double c_wound = 1.0e-4;
	double phif0_wound=0;
	double kappa0_wound = 0.5;
	double a0x = 1.;//frand(0,1.);
	double a0y = 0.;//sqrt(1-a0x*a0x);
	Vector2d a0_wound;a0_wound<<a0x,a0y;
	Vector2d lamda0_wound;lamda0_wound<<1.,1.;
	//---------------------------------//
	
	
	//---------------------------------//
	// values for the healthy
	double rho_healthy = 9000; // [cells/mm^3]
	double c_healthy = 0.0;
	double phif0_healthy=1.;
	double kappa0_healthy = 0.25;
	Vector2d a0_healthy;a0_healthy<<1.0,0.0;
	Vector2d lamda0_healthy;lamda0_healthy<<1.,1.;
	//---------------------------------//
	
	
	//---------------------------------//
	// create mesh (only nodes and elements)
	std::cout<<"Going to create the mesh\n";
	std::vector<double> rectangleDimensions = {0.0,100.0,0.0,100.0};
	std::vector<int> meshResolution = {50,50};
	QuadMesh myMesh = myRectangleMesh(rectangleDimensions, meshResolution);
	std::cout<<"Created the mesh with "<<myMesh.n_nodes<<" nodes and "<<myMesh.n_elements<<" elements\n";
	// print the mesh
	std::cout<<"nodes\n";
	for(int nodei=0;nodei<myMesh.n_nodes;nodei++){
		std::cout<<myMesh.nodes[nodei](0)<<","<<myMesh.nodes[nodei](1)<<"\n";
	}
	std::cout<<"elements\n";
	for(int elemi=0;elemi<myMesh.n_elements;elemi++){
		std::cout<<myMesh.elements[elemi][0]<<","<<myMesh.elements[elemi][1]<<","<<myMesh.elements[elemi][2]<<","<<myMesh.elements[elemi][3]<<"\n";
	}
	std::cout<<"boundary\n";
	for(int nodei=0;nodei<myMesh.n_nodes;nodei++){
		std::cout<<myMesh.boundary_flag[nodei]<<"\n";
	}
	// create the other fields needed in the tissue struct.
	int n_elem = myMesh.n_elements;
	int n_node = myMesh.n_nodes;
	//
	// global fields rho and c initial conditions 
	std::vector<double> node_rho0(n_node,rho_healthy);
	std::vector<double> node_c0 (n_node,c_healthy);
	//
	// values at the integration points
	std::vector<double> ip_phi0(n_elem*4,phif0_healthy);
	std::vector<Vector2d> ip_a00(n_elem*4,a0_healthy);
	std::vector<double> ip_kappa0(n_elem*4,kappa0_healthy);
	std::vector<Vector2d> ip_lamda0(n_elem*4,lamda0_healthy);
	//
	// define ellipse
	double x_center = 50.;
	double y_center = 50.;
	double x_axis = 10.;
	double y_axis = 10.;
	double alpha_ellipse = 0.;
	// boundary conditions and definition of the wound
	double tol_boundary = 1e-5;
	std::map<int,double> eBC_x;
	std::map<int,double> eBC_rho;
	std::map<int,double> eBC_c;
	for(int nodei=0;nodei<n_node;nodei++){
		double x_coord = myMesh.nodes[nodei](0);
		double y_coord = myMesh.nodes[nodei](1);
		// check
		if(myMesh.boundary_flag[nodei]==1){
			// insert the boundary condition for displacement
			std::cout<<"fixing node "<<nodei<<"\n";
			eBC_x.insert ( std::pair<int,double>(nodei*2+0,myMesh.nodes[nodei](0)) ); // x coordinate
			eBC_x.insert ( std::pair<int,double>(nodei*2+1,myMesh.nodes[nodei](1)) ); // y coordinate 
			// insert the boundary condition for rho
			eBC_rho.insert ( std::pair<int,double>(nodei,rho_healthy) ); 
			// insert the boundary condition for c
			eBC_c.insert   ( std::pair<int,double>(nodei,c_healthy) );
		}
		// check if it is in the center of the wound
		double check_ellipse = pow((x_coord-x_center)*cos(alpha_ellipse)+(y_coord-y_center)*sin(alpha_ellipse),2)/(x_axis*x_axis) +\
						pow((x_coord-x_center)*sin(alpha_ellipse)+(y_coord-y_center)*cos(alpha_ellipse),2)/(y_axis*y_axis) ;
		if(check_ellipse<=1){
			// inside
			std::cout<<"wound node "<<nodei<<"\n";
			node_rho0[nodei] = rho_wound;
			node_c0[nodei] = c_wound;
		}
	}
	for(int elemi=0;elemi<n_elem;elemi++){
		std::vector<Vector3d> IP = LineQuadriIP();
		for(int ip=0;ip<IP.size();ip++)
		{
			double xi = IP[ip](0);
			double eta = IP[ip](1);
			// weight of the integration point
			double wip = IP[ip](2);
			std::vector<double> R = evalShapeFunctionsR(xi,eta);
			Vector2d X_IP;X_IP.setZero();
			for(int nodej=0;nodej<4;nodej++){
				X_IP += R[nodej]*myMesh.nodes[myMesh.elements[elemi][nodej]];
			}
			double check_ellipse_ip = pow((X_IP(0)-x_center)*cos(alpha_ellipse)+(X_IP(1)-y_center)*sin(alpha_ellipse),2)/(x_axis*x_axis) +\
						pow((X_IP(0)-x_center)*sin(alpha_ellipse)+(X_IP(1)-y_center)*cos(alpha_ellipse),2)/(y_axis*y_axis) ;
			if(check_ellipse_ip<=1){
				std::cout<<"IP node: "<<4*elemi+ip<<"\n";
				ip_phi0[elemi*4+ip] = phif0_wound;
				ip_a00[elemi*4+ip] = a0_wound;
				ip_kappa0[elemi*4+ip] = kappa0_wound;
				ip_lamda0[elemi*4+ip] = lamda0_wound;
			}
		}
	}
	// no neumann boundary conditions. 
	std::map<int,double> nBC_x;
	std::map<int,double> nBC_rho;
	std::map<int,double> nBC_c;
	
	// initialize my tissue
	tissue myTissue;
	// connectivity
	myTissue.LineQuadri = myMesh.elements;
	// parameters
	myTissue.global_parameters = global_parameters;
	myTissue.local_parameters = local_parameters;
	//
	myTissue.node_X = myMesh.nodes;
	myTissue.node_x = myMesh.nodes;
	myTissue.node_rho_0 = node_rho0;
	myTissue.node_rho = node_rho0;
	myTissue.node_c_0 = node_c0;
	myTissue.node_c = node_c0;
	myTissue.ip_phif_0 = ip_phi0;	
	myTissue.ip_phif = ip_phi0;	
	myTissue.ip_a0_0 = ip_a00;		
	myTissue.ip_a0 = ip_a00;	
	myTissue.ip_kappa_0 = ip_kappa0;	
	myTissue.ip_kappa = ip_kappa0;	
	myTissue.ip_lamdaP_0 = ip_lamda0;	
	myTissue.ip_lamdaP = ip_lamda0;	
	//
	myTissue.eBC_x = eBC_x;
	myTissue.eBC_rho = eBC_rho;
	myTissue.eBC_c = eBC_c;
	myTissue.nBC_x = nBC_x;
	myTissue.nBC_rho = nBC_rho;
	myTissue.nBC_c = nBC_c;
	myTissue.time_final = 24*15;
	myTissue.time_step = 0.1;
	myTissue.tol = 1e-8;
	myTissue.max_iter = 20;
	myTissue.n_node = myMesh.n_nodes;
	myTissue.n_quadri = myMesh.n_elements;
	myTissue.n_IP = 4*myMesh.n_elements;
	//
	std::cout<<"filling dofs...\n";
	fillDOFmap(myTissue);
	std::cout<<"going to eval jacobians...\n";
	evalElemJacobians(myTissue);
	//
	//print out the Jacobians
	std::cout<<"element jacobians\nJacobians= ";
	std::cout<<myTissue.elem_jac_IP.size()<<"\n";
	for(int i=0;i<myTissue.elem_jac_IP.size();i++){
		std::cout<<"element: "<<i<<"\n";
		for(int j=0;j<4;j++){
			std::cout<<"ip; "<<i<<"\n"<<myTissue.elem_jac_IP[i][j]<<"\n";
		}
	}
	// print out the forward dof map
	std::cout<<"Total :"<<myTissue.n_dof<<" dof\n";
	for(int i=0;i<myTissue.dof_fwd_map_x.size();i++){
		std::cout<<"x node*2+coord: "<<i<<", dof: "<<myTissue.dof_fwd_map_x[i]<<"\n";
	}
	for(int i=0;i<myTissue.dof_fwd_map_rho.size();i++){
		std::cout<<"rho node: "<<i<<", dof: "<<myTissue.dof_fwd_map_rho[i]<<"\n";
	}
	for(int i=0;i<myTissue.dof_fwd_map_c.size();i++){
		std::cout<<"c node: "<<i<<", dof: "<<myTissue.dof_fwd_map_c[i]<<"\n";
	}
	//
	// 
	std::cout<<"going to start solver\n";
	// save a node and an integration point to a file
	std::vector<int> save_node;save_node.clear();
	std::vector<int> save_ip;save_ip.clear();
	
		
	std::stringstream ss;
	std::string filename = "data/fullwound06"+ss.str()+"_";
	
	
	
	//----------------------------------------------------------//
	// SOLVE
	sparseWoundSolver(myTissue, filename, 60,save_node,save_ip);
	//----------------------------------------------------------//
	
		
	return 0;	
}