#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "WriteVtuFunction.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
//#include <vtkDataSetMapper.h>
#include <vtkNew.h>
//#include <vtkProperty.h>
//#include <vtkRenderer.h>
#include <vtkTetra.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkHexahedron.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace Eigen;

void WriteVtu(const Eigen::MatrixXd Rho, const MESH myMesh, const FEA myFea) {
	// Access struct parameters
    Eigen::MatrixXd coordinates = myMesh.coordinate;
    Eigen::MatrixXd nodemap = myMesh.nodemap;
    Eigen::MatrixXd mesh_centers = myMesh.mesh_centers;
	Eigen::MatrixXd Vel = myFea.Vel;
	double nel = myMesh.nel;
	int dim = myMesh.dim;
	double nelx = myMesh.nelx;
	double nely = myMesh.nely;
	double nelz = myMesh.nelz;
	Eigen::MatrixXd T = myFea.T;
	Eigen::MatrixXd Kappa_Tau = myFea.Kappa_Tau;
	Eigen::MatrixXd U = myFea.U;
	Eigen::MatrixXd vonMises = myFea.vonMises;
	Eigen::MatrixXd stress_true = myFea.stress_true;


	/*************************************************************/
	std::cout << "Inputing Data into VTU format for visualization" << std::endl;
	vtkNew<vtkPoints> points; // coordinates
  	//int nx = 3;
  	//int ny = 2;
	//std::cout << "Adding Coordinate points" << std::endl;
  	for (int i = 0; i < coordinates.rows(); i++)
  	{
      points->InsertNextPoint(coordinates(i,0), coordinates(i,1), coordinates(i,2));
    }

	vtkNew<vtkUnstructuredGrid> unstructuredGrid; // mesh container
  	unstructuredGrid->SetPoints(points);

	//std::cout << "Adding nodemap" << std::endl;
	for (int i = 0; i < nodemap.rows(); i++)
  	{
		vtkNew<vtkHexahedron> cell;
      	cell->GetPointIds()->SetId(0, nodemap(i,0));
      	cell->GetPointIds()->SetId(1, nodemap(i,1));
      	cell->GetPointIds()->SetId(2, nodemap(i,2));
      	cell->GetPointIds()->SetId(3, nodemap(i,3));
		cell->GetPointIds()->SetId(4, nodemap(i,4));
      	cell->GetPointIds()->SetId(5, nodemap(i,5));
      	cell->GetPointIds()->SetId(6, nodemap(i,6));
      	cell->GetPointIds()->SetId(7, nodemap(i,7));
       	unstructuredGrid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
    }

	//std::cout << "Adding Temperature Data" << std::endl;
	vtkNew<vtkDoubleArray> dataX;
  	dataX->SetNumberOfComponents(1);
  	dataX->SetName("Temperature");
	for (int i = 0; i < coordinates.rows(); i++)
  	{
      dataX->InsertNextTuple1(T(i,0));
    }
	unstructuredGrid->GetPointData()->AddArray(dataX);

	//std::cout << "Adding Displacement Data" << std::endl;
	vtkNew<vtkDoubleArray> displacement;
  	displacement->SetNumberOfComponents(3);
  	displacement->SetName("Displacement");
	for (int i = 0; i < coordinates.rows(); i++)
  	{
		double U0 = U(i*dim); double U1 = U(i*dim + 1); double U2 = U(i*dim + 2);
      	displacement->InsertNextTuple3(U0,U1,U2);
    }
	unstructuredGrid->GetPointData()->AddArray(displacement);

	//std::cout << "Adding Stress Data" << std::endl;
	vtkNew<vtkDoubleArray> stress;
  	stress->SetNumberOfComponents(6);
  	stress->SetName("Stress");
	for (int i = 0; i < nodemap.rows(); i++)
  	{
		Eigen::MatrixXd elem_stress = stress_true.col(i);
		//double S0 = elem_stress(0); double S1 = elem_stress(1); double S2 = elem_stress(2);
		//double S0 = elem_stress(0); double S1 = elem_stress(1); double S2 = elem_stress(2);
      	stress->InsertNextTuple6(elem_stress(0),elem_stress(1),elem_stress(2),elem_stress(3),elem_stress(4),elem_stress(5));
    }
	unstructuredGrid->GetCellData()->AddArray(stress);

	//std::cout << "Adding VonMises Stress Data" << std::endl;
	vtkNew<vtkDoubleArray> vonMisesStress;
  	vonMisesStress->SetNumberOfComponents(1);
  	vonMisesStress->SetName("VonMises Stress");
	for (int i = 0; i < nodemap.rows(); i++)
  	{
		double vonMisesi = sqrt(vonMises(i,0));
      	vonMisesStress->InsertNextTuple1(vonMisesi);
    }
	unstructuredGrid->GetCellData()->AddArray(vonMisesStress);

	//std::cout << "Adding Density Data" << std::endl;
	vtkNew<vtkDoubleArray> rho1;
  	rho1->SetNumberOfComponents(1);
  	rho1->SetName("Rho 1");
	for (int i = 0; i < nodemap.rows(); i++)
  	{
		//double vonMisesi = vonMises(i,0);
      rho1->InsertNextTuple1(Rho(i,0));
    }
	unstructuredGrid->GetCellData()->AddArray(rho1);

	vtkNew<vtkDoubleArray> rho2;
  	rho2->SetNumberOfComponents(1);
  	rho2->SetName("Rho 2");
	for (int i = 0; i < nodemap.rows(); i++)
  	{
		//double vonMisesi = vonMises(i,0);
      rho2->InsertNextTuple1(Rho(i,1));
    }
	unstructuredGrid->GetCellData()->AddArray(rho2);

	//std::cout << "Adding Velocity Data" << std::endl;
	vtkNew<vtkDoubleArray> velocity;
  	velocity->SetNumberOfComponents(3);
  	velocity->SetName("Velocity");
	for (int i = 0; i < nodemap.rows(); i++)
  	{
		double V0 = Vel(i,0); double V1 = Vel(i,1); double V2 = Vel(i,2);
      	velocity->InsertNextTuple3(V0,V1,V2);
    }
	unstructuredGrid->GetCellData()->AddArray(velocity);

	//std::cout << "Adding Kappa_Tau Data" << std::endl;
	vtkNew<vtkDoubleArray> kappa_tau;
  	kappa_tau->SetNumberOfComponents(2);
  	kappa_tau->SetName("Kappa_Tau");
	for (int i = 0; i < nodemap.rows(); i++)
  	{
		double V0 = Kappa_Tau(i,0); double V1 = Kappa_Tau(i,1);
		///(1.204*1000)
      	kappa_tau->InsertNextTuple2(V0,V1);
    }
	unstructuredGrid->GetCellData()->AddArray(kappa_tau);

	// Write file.
  	vtkNew<vtkXMLUnstructuredGridWriter> writer;
	std::string filename = "Results/HeatExchangerResults" + std::to_string(myMesh.opt_iter) + ".vtu";
  	writer->SetFileName(filename.c_str());
  	writer->SetInputData(unstructuredGrid);
  	writer->Write();

}