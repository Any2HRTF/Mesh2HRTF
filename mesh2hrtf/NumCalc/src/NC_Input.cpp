/*===========================================================================*\
 *                                                                            *
 *  File name:      NC_Input.cpp                                              *
 *  Description:    Input data of a job                                       *
 *  Author:         H. Ziegelwanger, W. Kreuzer and Z.S. Chen                 *
 *                                                                            *
 \*===========================================================================*/



//================================= I N C L U D E S ==================================
// local includes                                                                   //
#include"NC_ConstantsVariables.h"                                                   //
#include"NC_TypeDefinition.h"                                                       //
#include"NC_Arrays.h"                                                               //
//                                                                                  //
// system includes                                                                  //
#include<iostream>                                                                  //
#include<iomanip>                                                                   //
#include<string>                                                                    //
#include<fstream> // file class functions                                           //
#include<sstream>                                                                   //
#include<math.h>                                                                    //
using namespace std;                                                                //
//====================================================================================



// local functions
void NC_ReadBasicParametersB(ofstream&,FILE* inputFile_,char*,string[],double*);
void NC_ReadMesh(ofstream&,FILE* inputFile_,char*, string[]);
void NC_ReadBoundaryConditions(ofstream&,FILE* inputFile_,char*,string[],int&,int&,inputLineBoundaryCondition *);
void NC_ReadSoundSources(ofstream&,FILE* inputFile_,char*, string[]);
void NC_ReadFrequencyCurves(ofstream&,FILE* inputFile_,char*, string[]);
void NC_StoreBoundaryConditions(ofstream&,FILE* inputFile_,char*,string[],int&,inputLineBoundaryCondition *);
void NC_WriteParameters(ofstream&,double*);
void NC_ComputeFrequencies(ofstream&,int&, double*,double*,double*);
int NC_GetLine(FILE *,char *,string[]);
int NC_String2Integer(string);
double NC_String2Double(string Ori_Str);



// read the input file
void NC_Read
(
	ofstream& NCout,
	FILE* inputFile_,			// W: input file
	char* chinpline,		// W: input line
	string chterms[NTRM],	// W: entries in the input line
	double* Freqs			// O: array of the frequencies
)
{
	int i, nabspsourc, maxbcline3, nbcline3;
	inputLineBoundaryCondition *bouconlin3;

	// read the basic parameters
	NC_ReadBasicParametersB(NCout, inputFile_, chinpline, chterms, Freqs);

	// alloctate the global arrays
	// these are needed later in the postprocessing
	extNumbersOfNodes = new int[numNodes_];
	isNodeMeshOrEvalNode = new int[numNodes_];
	extNumbersOfElements = new int[numElements_];

	nodesCoordinates = new double*[numNodes_];
	for(i=0; i<numNodes_; i++) nodesCoordinates[i] = new double[NDIM];

	elementsConnectivity = new int*[numElements_];
	for(i=0; i<numElements_; i++) elementsConnectivity[i] = new int[NNODPE];

	listNumberNodesPerElement = new int[numElements_];
	listElementProperty = new int[numElements_];
	listElementsElementGroup = new int[numElements_];

	listSymmetryPlanes = new int[NDIM];
	listSymmetryPlanesCoordinates = new double[NDIM];

	ibval = new int[numElements_];
	zbval0 = new Complex*[numElements_];
	for(i=0; i<numElements_; i++) zbval0[i] = new Complex[NNODPE];
	zbvao0 = new Complex*[numElements_];
	for(i=0; i<numElements_; i++) zbvao0[i] = new Complex[NNODPE];
	zbval1 = new Complex[numElements_];
	zbval2 = new Complex[numElements_];
	zbvao1 = new Complex[numElements_];
	zbvao2 = new Complex[numElements_];

	elenor = new double*[numElements_];
	for(i=0; i<numElements_; i++) elenor[i] = new double[(NNODPE + 1)*NDIM];
	centel = new double*[numElements_];
	for(i=0; i<numElements_; i++) centel[i] = new double[NDIM];
	areael = new double[numElements_];
	jelist = new int*[numElements_];
	for(i=0; i<numElements_; i++) jelist[i] = new int[2];
	for(i=0; i<numElements_; i++) jelist[i][0] = jelist[i][1] = -1;

	indexOfElementGroup = new int[numElementGroups_];
	numberOfElementsInGroup = new int[numElementGroups_];
	propertyOfGroup = new int[numElementGroups_];

	numinw = new int[numIncidentPlaneWaves_];
	nabspsourc = abs(numPointSources_);
	numposo = new int[nabspsourc];
	dirinw = new double*[numIncidentPlaneWaves_];
	for(i=0; i<numIncidentPlaneWaves_; i++) dirinw[i] = new double[NDIM];
	coorps = new double*[nabspsourc];
	for(i=0; i<nabspsourc; i++) coorps[i] = new double[NDIM];
	inwacurv = new int*[numIncidentPlaneWaves_];
	for(i=0; i<numIncidentPlaneWaves_; i++) inwacurv[i] = new int[2];
	posocurv = new int*[nabspsourc];
	for(i=0; i<nabspsourc; i++) posocurv[i] = new int[2];
	zinwav = new Complex[numIncidentPlaneWaves_];
	zposoc = new Complex[nabspsourc];
	zinwao = new Complex[numIncidentPlaneWaves_];
	zposoo = new Complex[nabspsourc];

	// read mesh
	NC_ReadMesh(NCout, inputFile_, chinpline, chterms);

	// read boundary conditions
	maxbcline3 = 2*numElementsOfBoundaryMesh_;
	bouconlin3 = new inputLineBoundaryCondition[maxbcline3];
	NC_ReadBoundaryConditions(NCout, inputFile_, chinpline, chterms, maxbcline3,
		nbcline3, bouconlin3);

	// read the incident waves
	if(numIncidentPlaneWaves_ > 0 || numPointSources_ > 0 || numPointSources_ == -1)
		NC_ReadSoundSources(NCout, inputFile_, chinpline, chterms);

	// read the curves and the data for post process
	NC_ReadFrequencyCurves(NCout, inputFile_, chinpline, chterms);

	// store the boundary conditions
	NC_StoreBoundaryConditions(NCout, inputFile_, chinpline, chterms, nbcline3,
		bouconlin3);
	delete [] bouconlin3;

	// write the input data into the out put file
	NC_WriteParameters(NCout, Freqs);
}

// program to read the basic parameters (part 1)
void NC_ReadBasicParametersA
(
	ofstream& NCout,
	FILE* inputFile_,			// input file
	char* chinpline,		// input line
	string chterms[NTRM]	// entries in the input line
)
{
    int i, nterms, j, j1;

    // version number
    nterms = 0;
    while(!nterms){nterms = NC_GetLine(inputFile_, chinpline, chterms);}
    if(chterms[0].compare("Mesh2HRTF")) NC_Error_Exit_0(NCout, "Key word Mesh2HRTF expected!");

    j1 = 0;
    //    for(j=0; j<(int)chterms[0].length(); j++) versionNumber_[j1++] = chterms[0][j];
    versionNumber_ = chterms[0];
    /* this does not make sense right now, but if in the future we need 
       a version number we have it in one of the chterms
    */
    for(i=1; i<nterms; i++)
    {
      /*        strcat(versionNumber_, " ");
        j1++;
        for(j=0; j<(int)chterms[i].length(); j++) versionNumber_[j1++] = chterms[i][j];
      */
      versionNumber_ = versionNumber_ + " " + chterms[i];
    }

    // job title
    nterms = 0;
    /* it is a bit weird to use NC_GetLine here that splits a string 
       at " ", just to merge it later, but it also gets rid of comments
    */
    while(!nterms){nterms = NC_GetLine(inputFile_, chinpline, chterms);}
    j1 = 0;
    // for(j=0; j<(int)chterms[0].length(); j++) jobTitle_[j1++] = chterms[0][j];
    jobTitle_ = chterms[0];
    for(i=1; i<nterms; i++)
    {
      jobTitle_ = jobTitle_ + " " + chterms[i];
      /*
        strcat(jobTitle_, " ");
        j1++;
        for(j=0; j<(int)chterms[i].length(); j++) jobTitle_[j1++] = chterms[i][j];
      */
    }

    cout << "Job title: " << jobTitle_ << endl;

    // Kontrollparameter I: analysis type, unit system
    while(!NC_GetLine(inputFile_, chinpline, chterms));

    // Kontrollparameter II
    while(!NC_GetLine(inputFile_, chinpline, chterms));
    numFrequencies_ = NC_String2Integer(chterms[1]);
    frequencyStepIncrement_ = NC_String2Double(chterms[2]);
    firstFrequencyStep_ = NC_String2Double(chterms[3]);
}

// program to read the basic parameters (part 2)
void NC_ReadBasicParametersB
(
	ofstream& NCout,
	FILE* inputFile_,			// input file
	char* chinpline,		// input line
	string chterms[NTRM],	// entries in the input line
	double* Freqs			// array of frequencies
)
{

	int i, npoicurv;
	double *p_dcurvfrq, *p_dcurvtime;

	while(!NC_GetLine(inputFile_, chinpline, chterms));
	npoicurv = NC_String2Integer(chterms[1]);

	p_dcurvfrq = new double[npoicurv];
	p_dcurvtime = new double[npoicurv];

	for(i=0; i<npoicurv; i++)
	{
		while(!NC_GetLine(inputFile_, chinpline, chterms));
		p_dcurvtime[i] = NC_String2Double(chterms[0]);
		p_dcurvfrq[i] = NC_String2Double(chterms[1]);
	}

	// compute the frequency array
	NC_ComputeFrequencies(NCout, npoicurv, p_dcurvfrq, p_dcurvtime, Freqs);

	delete [] p_dcurvfrq;
	delete [] p_dcurvtime;

	// BEM Main Parameters I
	while(!NC_GetLine(inputFile_, chinpline, chterms));
	numElementGroups_ = NC_String2Integer(chterms[0]);
	numElements_ = NC_String2Integer(chterms[1]);
	numNodes_ = NC_String2Integer(chterms[2]);
	minExpansionTermsFMM_ = NC_String2Integer(chterms[3]);
	numSymmetricPlanes_ = NC_String2Integer(chterms[4]);
	methodBEM_ = NC_String2Integer(chterms[7]);
	methodSolver_ = NC_String2Integer(chterms[8]);

	if(numElements_ < 1) NC_Error_Exit_1(NCout, "Number of elements must > 0!",
		"Number of elements = ", numElements_);

	if(numNodes_ < 1) NC_Error_Exit_1(NCout, "Number of nodes must > 0!", "Number of nodes = ", numNodes_);

	if(minExpansionTermsFMM_ == 0) minExpansionTermsFMM_ = 8;
	if(minExpansionTermsFMM_ < 5 || minExpansionTermsFMM_ > 10) {
		NC_Error_Exit_1(NCout, "minExpansionTermsFMM_  must be in [5 10]!", "minExpansionTermsFMM_ = ", minExpansionTermsFMM_);
    }

	if(numSymmetricPlanes_ < 0 || numSymmetricPlanes_ > 3)
		NC_Error_Exit_1(NCout, "Number of symmetry planes must be in [0 3]!",
		"Number of symmetry planes = ", numSymmetricPlanes_);

	if(methodSolver_ < 0 || methodSolver_ > 4)
        NC_Error_Exit_1(NCout, "ISOLVER must be 0 or 4!", "ISOLVER = ", methodSolver_);

    if(methodSolver_ > 0 && methodSolver_ < 4)
        NC_Error_Exit_1(NCout, "ISOLVER must be 0 or 4!", "ISOLVER = ", methodSolver_);

	switch(methodBEM_)    // input parameter
	{
	case 0:
		methodFMM_ = 0; // TBEM
		break;
	case 1:
		methodFMM_ = 1; // SLFMBEM
		break;
	case 4:
		methodFMM_ = 3; // DMLFMBEM
		break;
	default:
		NC_Error_Exit_1(NCout, "methodBEM_ must be in {0,1,4}!", "methodBEM_ = ", methodBEM_);
	}

	if(methodFMM_ && methodSolver_ == 4) methodSolver_ = 0;

	// BEM Main Parameters II
	while(!NC_GetLine(inputFile_, chinpline, chterms));
	numIncidentPlaneWaves_ = NC_String2Integer(chterms[0]);
	numPointSources_ = NC_String2Integer(chterms[1]);
	isInternalProblem_ = NC_String2Integer(chterms[2]);
	ClusEdgL0_ = NC_String2Double(chterms[3]);
	methodPreconditioner_ = NC_String2Integer(chterms[4]);

	if(numIncidentPlaneWaves_ < 0) NC_Error_Exit_1(NCout, "Number of the incident plane waves must >= 0!",
		"Number of the incident plane waves = ", numIncidentPlaneWaves_);
	if(numPointSources_ < -1) NC_Error_Exit_1(NCout, "Number of the point sources must >= -1!",
		"Number of the point sources = ", numIncidentPlaneWaves_);
	if(isInternalProblem_ != 0 && isInternalProblem_ != 1)
		NC_Error_Exit_1(NCout, "IFINTE must = 0 or 1!", "IFINTE = ", isInternalProblem_);
	if(methodPreconditioner_ < -3 || methodPreconditioner_ > 2)
		NC_Error_Exit_1(NCout, "IPRECOND must be in [-3 2]!", "IPRECOND = ", methodPreconditioner_);

	if(methodPreconditioner_ <= -3) {
		scanningDegreeLU_ = 3;
		methodPreconditioner_ = 0;
	} else if(methodPreconditioner_ == -2) {
		scanningDegreeLU_ = 2;
		methodPreconditioner_ = 0;
	} else if(methodPreconditioner_ == -1) {
		scanningDegreeLU_ = 1;
		methodPreconditioner_ = 0;
	} else {
		scanningDegreeLU_ = 0;
	}

	// if the direct method is used, no preconditioner
	if(methodSolver_ == 4) methodPreconditioner_ = 2;

	// BEM Main Parameters III
	while(!NC_GetLine(inputFile_, chinpline, chterms));

	// BEM Main Parameters IV
	while(!NC_GetLine(inputFile_, chinpline, chterms));
	speedOfSound_ = NC_String2Double(chterms[0]);
	densityOfMedium_ = NC_String2Double(chterms[1]);
	harmonicTimeFactor_ = NC_String2Double(chterms[2]);

	if(harmonicTimeFactor_ != -1.0) harmonicTimeFactor_ = 1.0;

    if(speedOfSound_ < 300.0) NC_Error_Warning_1(NCout, "speedOfSound_ perhaps too small!",
        "speedOfSound_ = ", speedOfSound_);
    if(densityOfMedium_ < 1.0) NC_Error_Warning_1(NCout, "densityOfMedium_ perhaps too small!",
        "densityOfMedium_ = ", densityOfMedium_);
    epsilon_ = 1.0e-6;

	// compute the number of reflections with respect to the symmetric planes
	numReflectionsOfElements_ = (int)pow(2.0, numSymmetricPlanes_);

	// set the parameter Tao_
	if(isInternalProblem_) Tao_ = -1.0; else Tao_ = 1.0;
}

// read mesh
void NC_ReadMesh
(
	ofstream& NCout,
	FILE* inputFile_,			// input file
	char* chinpline,		// input line
	string chterms[NTRM]	// entries in the input line
)
{
	int i, j = 0, k = 0, nterms;
	double areevl;

	// read the key word NODES
	while(!NC_GetLine(inputFile_, chinpline, chterms));
	if(chterms[0].compare("NODES")) NC_Error_Exit_0(NCout, "Key word NODES expected!");

    i=0;
    while(i<numNodes_)
    {
        while(!NC_GetLine(inputFile_, chinpline, chterms));
        tmpFile_ = fopen(chterms[0].c_str(), "r");
        if(tmpFile_) // If file exists
        {
            while(!NC_GetLine(tmpFile_, chinpline, chterms));
            int tmpFileLength = NC_String2Integer(chterms[0]);
            for(int jj=0; jj<tmpFileLength; jj++)
            {
                while(!NC_GetLine(tmpFile_, chinpline, chterms));
                extNumbersOfNodes[i] = NC_String2Integer(chterms[0]);
                nodesCoordinates[i][0] = NC_String2Double(chterms[1]);
                nodesCoordinates[i][1] = NC_String2Double(chterms[2]);
                nodesCoordinates[i][2] = NC_String2Double(chterms[3]);
                i+=1;
            }
        }
        else // If file does not exist
        {
            cout << "No Nodes file found!" << endl;
            NCout << "No Nodes file found!" << endl;
	    exit(-1);
        }
	fclose(tmpFile_); // everything that is open should be closed
    }

	// read the key word ELEMENTS
	while(!NC_GetLine(inputFile_, chinpline, chterms));
	if(chterms[0].compare("ELEMENTS")) NC_Error_Exit_0(NCout, "Key word ELEMENTS expected!");

    i=0;
    while(i<numElements_)
    {
        while(!NC_GetLine(inputFile_, chinpline, chterms));
        tmpFile_ = fopen(chterms[0].c_str(), "r");
        if(tmpFile_) // If file exists
        {
            while(!NC_GetLine(tmpFile_, chinpline, chterms));
            int tmpFileLength = NC_String2Integer(chterms[0]);

            for(int jj=0; jj<tmpFileLength; jj++)
            {
                nterms = 0;
                while(!nterms){nterms = NC_GetLine(tmpFile_, chinpline, chterms);}

                if(nterms < 7 || nterms > 8) NC_Error_Exit_2(NCout,
                                                        "number of entris of a input line must be 7 or 8!",
                                                        "number of the element under the key word ELEMENTS = ", NC_String2Integer(chterms[0]),
                                                        "number of entrie of the line = ", nterms);

                extNumbersOfElements[i] = NC_String2Integer(chterms[0]);
                listNumberNodesPerElement[i] = nterms - 4;
                for(j=0; j<listNumberNodesPerElement[i]; j++) elementsConnectivity[i][j] = NC_String2Integer(chterms[j + 1]);
                listElementProperty[i] = NC_String2Integer(chterms[nterms - 3]);
                listElementsElementGroup[i] = NC_String2Integer(chterms[nterms - 1]);
                i+=1;
            }
        }
        else // If file does not exist
        {
            cout << "No Elements file found!" << endl;
            NCout << "No Elements file found!" << endl;
	    exit(-1);
        }
	fclose(tmpFile_);
    }

	// compute the number of BEs and the number of evaluation elements
	numElementsOfBoundaryMesh_ = numElementsOfEvaluationMesh_ = 0;
	for(i=0; i<numElements_; i++)
	{
		if(listElementProperty[i] == 0)
		{
			numElementsOfBoundaryMesh_++;
		}
		else if(listElementProperty[i] == 2)
		{
			if(numElementsOfEvaluationMesh_ == 0)
			{
				j = listElementsElementGroup[i];
				k = i;
			}
			else if(listElementsElementGroup[i] != j)
			{
				cout << "\n" << "Element " << extNumbersOfElements[k] << " belongs to the el. group " <<
					listElementsElementGroup[k] << " but" << endl;
				cout << "Element " << extNumbersOfElements[i] << " belongs to the el. group " <<
					listElementsElementGroup[i] << "!" << endl;
				NCout << "\n" << "Element " << extNumbersOfElements[k] << " belongs to the el. group " <<
					listElementsElementGroup[k] << " but" << endl;
				NCout << "Element " << extNumbersOfElements[i] << " belongs to the el. group " <<
					listElementsElementGroup[i] << "!" << endl;
				NC_Error_Exit_0(NCout, "Elements of the evaluation mesh must be of the same el. group!");
			}
			numElementsOfEvaluationMesh_++;
		}
		else NC_Error_Exit_2(NCout, "listElementProperty[i] must be 0 or 2!", "i = ", i, "listElementProperty[i] = ",
			listElementProperty[i]);
	} // end of loop I

	// input the symmetric planes
	if(numSymmetricPlanes_ > 0)
	{
		while(!NC_GetLine(inputFile_, chinpline, chterms));
		if(chterms[0].compare("SYMMETRY")) NC_Error_Exit_0(NCout, "Key word SYMMETRY expected!");

		while(!NC_GetLine(inputFile_, chinpline, chterms));
		for(i=0; i<NDIM; i++) listSymmetryPlanes[i] = NC_String2Integer(chterms[i]);

		while(!NC_GetLine(inputFile_, chinpline, chterms));
		for(i=0; i<NDIM; i++) listSymmetryPlanesCoordinates[i] = NC_String2Double(chterms[i]);

		j = 0;
		for(i=0; i<NDIM; i++)
		{
			if(listSymmetryPlanes[i]) j++;
			switch(listSymmetryPlanes[i])
			{
			case 0:
			case 1:
				break;
			case 2:
				listSymmetryPlanes[i] = -1;
				break;
			case 3:
				listSymmetryPlanes[i] = 2;
				break;
			case 4:
				listSymmetryPlanes[i] = -2;
				break;
			default:
				NC_Error_Exit_2(NCout, "listSymmetryPlanes[i] must be equal to 0, 1, 2, 3, or 4!",
					"i = ", i, "listSymmetryPlanes[i] = ", listSymmetryPlanes[i]);
			}
		}
		if(j != numSymmetricPlanes_) NC_Error_Exit_2(NCout, "j must equal to numSymmetricPlanes_", "j = ", j,
			"numSymmetricPlanes_ = ", numSymmetricPlanes_);
	}

	// change the connectivities: replace the external nodal numbers with the corresponding internal nodal numbers
	for(i=0; i<numElements_; i++)
	{
		for(j=0; j<listNumberNodesPerElement[i]; j++)
		{
			for(k=0; k<numNodes_; k++)
			{
				if(elementsConnectivity[i][j] == extNumbersOfNodes[k])
				{
					elementsConnectivity[i][j] = k;
					goto outofloopk_cc;
				}
			}
			NC_Error_Exit_2(NCout, "Node of element does not exist in the connectivity list!",
				"Number of element = ", extNumbersOfElements[i],
				"Local number of the node in the element = ", j);
outofloopk_cc:
			continue;
		}
	}

	// compute (1) normal vector to each element at its center
	//         (2) area of each element
	//         (3) coordinates of the center point of each element
	Matrix<double>crdel(NNODPE, NDIM);
	Vector<double> crdpoi(NDIM), shfunx(NNODPE), elnorv(NDIM), shnds(NNODPE),
		shndt(NNODPE), dxds(NDIM), dxdt(NDIM);


    averageElementArea_ = areevl = 0.0;
	int nbouel = 0, inode, nevlel = 0;
	double s, t, are0 = 0.0, leng;

	for(i=0; i<numElements_; i++)
	{

		if(listElementProperty[i] < 0 || listElementProperty[i] > 2 || listElementProperty[i] == 1)
			NC_Error_Exit_2(NCout, "listElementProperty[i] must be equal to 0 or 2!",
			"number of the element = ", extNumbersOfElements[i], "listElementProperty[i] = ", listElementProperty[i]);

		inode = listNumberNodesPerElement[i];

		for(j=0; j<inode; j++) for(k=0; k<NDIM; k++) crdel(j, k) = nodesCoordinates[elementsConnectivity[i][j]][k];

		switch(inode)
		{
			case 3:
				s = t = 1.0/3.0;
				are0 = 0.5;
				break;
			case 4:
				s = t = 0.0;
				are0 = 4.0;
				break;
			default: NC_Error_Exit_2(NCout, "Number of nodes of an element must be 3 or 4!",
						 "number of nodes of the element = ", inode,
						 "number of the element = ", extNumbersOfElements[i]);
		}


		NC_ComputeGlobalCoordinates(crdpoi, shfunx, crdel, s, t, inode);
		leng = NC_ComputeUnitVectorOnElement(elnorv, shnds, shndt, dxds, dxdt, crdel, s, t, inode);

		for(k=0; k<NDIM; k++)
		{
			elenor[i][k] = elnorv[k];
			centel[i][k] = crdpoi[k];
		}
		areael[i] = leng*are0;

		if(listElementProperty[i] == 2)
		{
			areevl += areael[i];
			nevlel++;
		}
		else
		{
			averageElementArea_ += areael[i];
			nbouel++;
		}

		for(j=0; j<inode; j++)
		{
			if(inode == NETYP3)
			{
				s = CSI6[j];
				t = ETA6[j];
			}
			else
			{
				s = CSI8[j];
				t = ETA8[j];
			}
		    leng = NC_ComputeUnitVectorOnElement(elnorv, shnds, shndt, dxds, dxdt, crdel,
			                s, t, inode);

			for(k=0; k<NDIM; k++)
			{
				elenor[i][(j + 1)*NDIM + k] = elnorv[k];
			}
		}
	} // end of loop I //

	averageElementArea_ /= (double)nbouel;
	if(nevlel > 0) areevl /= (double)nevlel;

	// compute the array isNodeMeshOrEvalNode and the variables NNODBE and NNODIN
	// isNodeMeshOrEvalNode[i] >=  0: I is node of the BE mesh
	//              = -1: it is a node of the evaluation mesh
	for(i=0; i<numNodes_; i++) isNodeMeshOrEvalNode[i] = -3;

	for(i=0; i<numElements_; i++)
	{
		if(listElementProperty[i] == 0)
		{
			for(j=0; j<listNumberNodesPerElement[i]; j++) isNodeMeshOrEvalNode[elementsConnectivity[i][j]] = 1;
		}
	}

	for(i=0; i<numElements_; i++)
	{
		if(listElementProperty[i] == 2)
		{
			for(j=0; j<listNumberNodesPerElement[i]; j++)
			{
				if(isNodeMeshOrEvalNode[elementsConnectivity[i][j]] != 1)
				{
					isNodeMeshOrEvalNode[elementsConnectivity[i][j]] = -1;
				}
				else NC_Error_Exit_1(NCout, "BE mesh and evaluation mesh have common node!",
						"Number of the common node = ", extNumbersOfNodes[elementsConnectivity[i][j]]);
			}
		}
	}

	numNodesOfBoundaryMesh_ = numNodesOfEvaluationMesh_ = 0;
	for(i=0; i<numNodes_; i++)
	{
		if(isNodeMeshOrEvalNode[i] == 1)
		{
			isNodeMeshOrEvalNode[i] = numNodesOfBoundaryMesh_;
			numNodesOfBoundaryMesh_++;
		}
		else if(isNodeMeshOrEvalNode[i] == -1)
		{
			numNodesOfEvaluationMesh_++;
		}
		else if(isNodeMeshOrEvalNode[i] == -3)
		{
			cout << "\n" << "number of the isolated node: " << extNumbersOfNodes[i] << endl;
			NCout << "\n" << "number of the isolated node: " << extNumbersOfNodes[i] << endl;
			NC_Error_Exit_0(NCout, "Isolated node found!");
		}
	}

	if(numNodesOfBoundaryMesh_ == 0)
		NC_Error_Exit_0(NCout, "Number of boundary elements is zero!");

	// compute the dimension of the box including all boundary elements
	double crdmaxmin[3][2], dtl = sqrt(averageElementArea_)/100.0;

	for(j=0; j<NDIM; j++)
	{
		crdmaxmin[j][0] = -9.9e20;
		crdmaxmin[j][1] = 9.9e20;
	}
	for(i=0; i<numNodes_; i++)
	{
		if(isNodeMeshOrEvalNode[i] < 0) continue;
		for(j=0; j<NDIM; j++)
		{
			if(nodesCoordinates[i][j] > crdmaxmin[j][0]) crdmaxmin[j][0] = nodesCoordinates[i][j];
			if(nodesCoordinates[i][j] < crdmaxmin[j][1]) crdmaxmin[j][1] = nodesCoordinates[i][j];
		}
	}

	for(j=0; j<NDIM; j++) {
		if(numSymmetricPlanes_ == 0) {
			crdmaxmin[j][0] -= crdmaxmin[j][1];
		} else {
			if(listSymmetryPlanes[j] == 0) {
				crdmaxmin[j][0] -= crdmaxmin[j][1];
			} else if(abs(listSymmetryPlanes[j]) <= 2) {
				if(crdmaxmin[j][0] >= listSymmetryPlanesCoordinates[j]-dtl && crdmaxmin[j][1] >= listSymmetryPlanesCoordinates[j]-dtl){
					crdmaxmin[j][0] = 2.0*(crdmaxmin[j][0] - listSymmetryPlanesCoordinates[j]);
				} else if(crdmaxmin[j][0] <= listSymmetryPlanesCoordinates[j]+dtl && crdmaxmin[j][1] <= listSymmetryPlanesCoordinates[j]+dtl) {
					crdmaxmin[j][0] = 2.0*(listSymmetryPlanesCoordinates[j] - crdmaxmin[j][1]);
				} else {
					NC_Error_Warning_1(NCout, "", "Number of the symmetric plane = ", j);
					NC_Error_Exit_2(NCout,
						"All BE nodes must be located on the same side of a symmetric plane!",
						"crdmaxmin[j][0] - listSymmetryPlanesCoordinates[j] = ", crdmaxmin[j][0] - listSymmetryPlanesCoordinates[j],
						"crdmaxmin[j][1] - listSymmetryPlanesCoordinates[j] =  = ", crdmaxmin[j][1] - listSymmetryPlanesCoordinates[j]);
				}
			} else {
				NC_Error_Exit_2(NCout, "listSymmetryPlanes[j] mist be 0, -1, 1, -2, 2!", "listSymmetryPlanes[j] = ", listSymmetryPlanes[j],
					"j = ", j);
			}
		}
	}
	diameterBE_ = 0.0;
	for(j=0; j<NDIM; j++) diameterBE_ += crdmaxmin[j][0]*crdmaxmin[j][0];
	diameterBE_ = sqrt(diameterBE_);

	// compute the number of each element group
	int ielgrp = 0;
	for(i=0; i<numElements_; i++)
	{
		for(j=0; j<ielgrp; j++)
		{
			if(listElementsElementGroup[i] == indexOfElementGroup[j]) goto EndOfLoopI_Elgr;
		}
		indexOfElementGroup[ielgrp] = listElementsElementGroup[i];
		propertyOfGroup[ielgrp] = listElementProperty[i];
		ielgrp++;
EndOfLoopI_Elgr: continue;
	}

	if(ielgrp < numElementGroups_) numElementGroups_ = ielgrp;
	else if(ielgrp > numElementGroups_)
		NC_Error_Exit_2(NCout, "ielgrp must be less than or equal to nelgrp", "ielgrp = ",
		ielgrp, "nelgrp = ", numElementGroups_);

	// compute number of elements of each element group
	for(j=0; j<numElementGroups_; j++) numberOfElementsInGroup[j] = 0;
	for(i=0; i<numElements_; i++)
	{
		for(j=0; j<numElementGroups_; j++)
		{
			if(listElementsElementGroup[i] == indexOfElementGroup[j])
			{
				numberOfElementsInGroup[j]++;
				break;
			}
		}
	}
}

// read the boundary conditions
void NC_ReadBoundaryConditions
(
	ofstream& NCout,
	FILE* inputFile_,				// input file
	char* chinpline,			// input line
	string chterms[NTRM],		// entries in the input line
	int& maxbcline,			//I: maximum number of input lines
	int& nbcline,				//O: actual number of input lines
	inputLineBoundaryCondition *bouconlin	//O: array to store the input lines
)
{
    int i, j, nterms, j_imag_p, j_imag_n, j_nega;

    // read the key word
    while(!NC_GetLine(inputFile_, chinpline, chterms));
    if(chterms[0].compare("BOUNDARY")) NC_Error_Exit_0(NCout, "Key word BOUNDARY expected!");

    // read the input lines
    nbcline = -1;
    for(i=0; i<maxbcline; i++)
    {
        nterms = 0;
        while(!nterms){nterms = NC_GetLine(inputFile_, chinpline, chterms);}

        if(!chterms[0].compare("RETU"))
        {
            nbcline = i;
            break;
        }
        else
        {	// chterms contains:
			// 0: ELEM (Keyword)
			// 1: Index of first element
			// 2: TO (Keyword)
			// 3: Index of last element
			// 4: Boundary condition (ADMI, IMPE, PRES, VELO)
			// 5: real value of boundary condition
			// 6: Curve id to defines the real part or -1 (no curve used)
			// 7: imaginary value of boundary condition
			// 8: Curve id to defines the imaginary part or -1 (no curve used)
			// Indicees of elements to which the boundary condition applies
            if(nterms < 9) NC_Error_Exit_1(NCout, "A input line under BOUNDARY is too short!",
                                           "Number of the input line = ", i);
            if(chterms[0].compare("ELEM")) NC_Error_Exit_0(NCout, "Key word ELEM expected!");
            bouconlin[i].nLow = NC_String2Integer(chterms[1]);
            if(chterms[2].compare("TO")) NC_Error_Exit_0(NCout, "Key word TO expected!");
            bouconlin[i].nHigh = NC_String2Integer(chterms[3]);

			// Type of boundary condition
            bouconlin[i].sKeyword = chterms[4];

			// check if a curve is defined for both parts
			if ( (chterms[6] == "-1" && chterms[8] != "-1") ||
			     (chterms[6] != "-1" && chterms[8] == "-1")){
				NC_Error_Exit_0(NCout,
					"ERROR: Curves have to be specified for all parts of a "
					"Boundary or Source. I.e., ALL curve numbers in one line "
					"must be -1 or ALL curve numbers must not be -1");
			}

			// read data for real part
			bouconlin[i].dReal = NC_String2Double(chterms[5]);
			if (chterms[6] == "-1"){
				// do not use frequency curve
				bouconlin[i].bRef = false;
			} else {
				// use frequency curve
				bouconlin[i].bRef = true;
				bouconlin[i].nRealRef = NC_String2Integer(chterms[6]);
			}
			// read data for negative part
			bouconlin[i].dImag = NC_String2Double(chterms[7]);
			if (chterms[8] == "-1"){
				// do not use frequency curve
				bouconlin[i].bRef = false;
			} else {
				// use frequency curve
				bouconlin[i].bRef = true;
                bouconlin[i].nImagRef = NC_String2Integer(chterms[8]);
			}
			// flag for defining infinitly thin elements that is not contained
			// in Mesh2HRTF
			bouconlin[i].bNega = false;

        } // ELSE
    } // loop I

    if(nbcline == -1) NC_Error_Exit_0(NCout, "Key word RETU expected!");
}
/* // dperecated version that did not read the data correct because it was
// looking for keywords 'IMAG' and 'NEGA'
void NC_ReadBoundaryConditions
(
	ofstream& NCout,
	FILE* inputFile_,				// input file
	char* chinpline,			// input line
	string chterms[NTRM],		// entries in the input line
	int& maxbcline,			//I: maximum number of input lines
	int& nbcline,				//O: actual number of input lines
	inputLineBoundaryCondition *bouconlin	//O: array to store the input lines
)
{
    int i, j, nterms, j_imag_p, j_imag_n, j_nega;

    // read the key word
    while(!NC_GetLine(inputFile_, chinpline, chterms));
    if(chterms[0].compare("BOUNDARY")) NC_Error_Exit_0(NCout, "Key word BOUNDARY expected!");

    // read the input lines
    nbcline = -1;
    for(i=0; i<maxbcline; i++)
    {
        nterms = 0;
        while(!nterms){nterms = NC_GetLine(inputFile_, chinpline, chterms);}

        if(!chterms[0].compare("RETU"))
        {
            nbcline = i;
            break;
        }
        else
        {	// chterms contains:
			// 0: ELEM (Keyword)
			// 1: Index of first element
			// 2: TO (Keyword)
			// 3: Index of last element
			// 4: Boundary condition (ADMI, IMPE, PRES, VELO)
			// 5: real value of boundary condition
			// 6: Curve id to defines the real part or -1 (no curve used)
			// 7: imaginary value of boundary condition
			// 8: Curve id to defines the imaginary part or -1 (no curve used)
			// Indicees of elements to which the boundary condition applies
            if(nterms < 6) NC_Error_Exit_1(NCout, "A input line under BOUNDARY is too short!",
                                           "Number of the input line = ", i);
            if(chterms[0].compare("ELEM")) NC_Error_Exit_0(NCout, "Key word ELEM expected!");
            bouconlin[i].nLow = NC_String2Integer(chterms[1]);
            if(chterms[2].compare("TO")) NC_Error_Exit_0(NCout, "Key word TO expected!");
            bouconlin[i].nHigh = NC_String2Integer(chterms[3]);

			// Type of boundary condition
            bouconlin[i].sKeyword = chterms[4];

            bouconlin[i].bRef = false;
            bouconlin[i].dReal = NC_String2Double(chterms[5]);

            j_imag_p = -1;  // flag to check if imaginary part is specified
			j_nega = -1;    // flag to check if negative side is specified
			j_imag_n = -1;  // flag to check if imaginary part of negative side
			                // is specified
            for(j=6; j<nterms; j++)
            {
                if(!chterms[j].compare("NEGA")) j_nega = j;
            }
            if(j_nega > 0)
            {
                for(j=6; j<j_nega; j++) if(!chterms[j].compare("IMAG")) j_imag_p = j;
                for(j=j_nega+1; j<nterms; j++) if(!chterms[j].compare("IMAG")) j_imag_n = j;
            }
            else
            {
                for(j=6; j<nterms; j++) if(!chterms[j].compare("IMAG")) j_imag_p = j;
            }

			// Use a frequency curve for the real part
            if(j_imag_p == 7 || j_nega == 7 || nterms == 7)
            {
                bouconlin[i].bRef = true;
                bouconlin[i].nRealRef = NC_String2Integer(chterms[6]);
            }

            if(j_imag_p > 0) // there is a imaginary part
            {
                bouconlin[i].dImag = NC_String2Double(chterms[j_imag_p + 1]);
                if(j_nega == j_imag_p + 3 || nterms == j_imag_p + 3)
                {
                    bouconlin[i].bRef = true;
                    bouconlin[i].nImagRef = NC_String2Integer(chterms[j_imag_p + 2]);
                }
            }
            else bouconlin[i].dImag = 0.0;

            if(j_nega > 0) // there is a negative side
            {
                bouconlin[i].bNega = true;
                bouconlin[i].dNegReal = NC_String2Double(chterms[j_nega + 1]);
                if(j_imag_n == j_nega + 3 || nterms == j_nega + 3)
                {
                    bouconlin[i].bRef = true;
                    bouconlin[i].nNegRealRef = NC_String2Integer(chterms[j_nega + 2]);
                }

                if(j_imag_n > 0) // there is a imaginary part
                {
                    bouconlin[i].dNegImag = NC_String2Double(chterms[j_imag_n + 1]);
                    if(nterms == j_imag_n + 3)
                    {
                        bouconlin[i].bRef = true;
                        bouconlin[i].nNegImagRef =  NC_String2Integer(chterms[j_imag_n + 2]);
                    }
                }
                else bouconlin[i].dNegImag = 0.0;
            }
            else bouconlin[i].bNega = false;
        } // ELSE
    } // loop I

    if(nbcline == -1) NC_Error_Exit_0(NCout, "Key word RETU expected!");
}*/

// read the incident waves
void NC_ReadSoundSources
(
	ofstream& NCout,
	FILE* inputFile_,			// input file
	char* chinpline,		// input line
	string chterms[NTRM]	// entries in the input line
)
{
	int i, nterms;
	double dw, re, im;

	// read the plane waves
	if(numIncidentPlaneWaves_ > 0)
	{
		while(!NC_GetLine(inputFile_, chinpline, chterms));
		if(chterms[0].compare("PLANE")) NC_Error_Exit_0(NCout, "Key word PLANE WAVE expected!");

		for(i=0; i<numIncidentPlaneWaves_; i++)
		{
			nterms = 0;
			while(!nterms){nterms = NC_GetLine(inputFile_, chinpline, chterms);}

			if(nterms != 8)
				NC_Error_Exit_2(NCout,
				"Number of terms of an input line for plane wave must be 8!",
				"Number of terms = ", nterms,
				"Number of the input line under the key word PLANE WAVE = ", i+1);

			numinw[i] = NC_String2Integer(chterms[0]);
			dirinw[i][0] = NC_String2Double(chterms[1]);
			dirinw[i][1] = NC_String2Double(chterms[2]);
			dirinw[i][2] = NC_String2Double(chterms[3]);
			Vcnorm_dim3_(dirinw[i], dw);
			re = NC_String2Double(chterms[4]);
			im = NC_String2Double(chterms[6]);
			zinwao[i].set(re, im);
			inwacurv[i][0] = NC_String2Integer(chterms[5]);
			inwacurv[i][1] = NC_String2Integer(chterms[7]);
		}
	}

	if(numPointSources_ > 0)
	{
		while(!NC_GetLine(inputFile_, chinpline, chterms));
		if(chterms[0].compare("POINT")) NC_Error_Exit_0(NCout, "Key word POINT SOURCES expected!");

		for(i=0; i<numPointSources_; i++)
		{
			nterms = 0;
			while(!nterms) {nterms = NC_GetLine(inputFile_, chinpline, chterms);}

			if(nterms != 8)
				NC_Error_Exit_2(NCout,
				"Number of terms of an input line for pont source must be 8!",
				"Number of terms = ", nterms,
				"Number of the input line under the key word POINT SOURCE = ", i+1);

			numposo[i] = NC_String2Integer(chterms[0]);
			coorps[i][0] = NC_String2Double(chterms[1]);
			coorps[i][1] = NC_String2Double(chterms[2]);
			coorps[i][2] = NC_String2Double(chterms[3]);
			re = NC_String2Double(chterms[4]);
			im = NC_String2Double(chterms[6]);
			zposoo[i].set(re, im);
			posocurv[i][0] = NC_String2Integer(chterms[5]);
			posocurv[i][1] = NC_String2Integer(chterms[7]);
		}
	}
}

// read the frequency curves and the data for post process (not be used)
void NC_ReadFrequencyCurves
(
	ofstream& NCout,
	FILE* inputFile_,			// input file
	char* chinpline,		// input line
	string chterms[NTRM]	// entries in the input line
)
{
    int i, j, nterms, mpcurv;
    bool bifrdcurv = false;

    // create and initilize the array IBVCURV
    ibvcurv = new int*[numElements_];
    for(i=0; i<numElements_; i++) ibvcurv[i] = new int[6];
    for(i=0; i<numElements_; i++) for(j=0; j<6; j++) ibvcurv[i][j] = -1;

    // read the key word
    while(!NC_GetLine(inputFile_, chinpline, chterms));
    if(!chterms[0].compare("CURVES"))
    {
        bifrdcurv = true;
    }
    else if(chterms[0].compare("POST"))
    {
        NC_Error_Exit_0(NCout, "Key word CURVES or POST PROCESS expected!");
    }

    // read number of curves and the maximun point number
    numCurvesFrequency_ = mpcurv = 0;
    if(bifrdcurv)
    {
        nterms = 0;
        while(!nterms){nterms = NC_GetLine(inputFile_, chinpline, chterms);}

        if(nterms != 2) NC_Error_Exit_1(NCout,
                                        "Number of terms of the input line direct under the key word CURVE must be 2!",
                                        "Number of terms = ", nterms);
        numCurvesFrequency_ = NC_String2Integer(chterms[0]);
        mpcurv = NC_String2Integer(chterms[1]);
    }
    if(numCurvesFrequency_ < 0) numCurvesFrequency_ = 0;
    if(mpcurv < 0) mpcurv = 0;

    // create the arrays to store the frequency curves
    nucurv = new int[numCurvesFrequency_];
    npcurv = new int[numCurvesFrequency_];
    frqcurv = new double*[numCurvesFrequency_];
    for(i=0; i<numCurvesFrequency_; i++) frqcurv[i] = new double[mpcurv];
    faccurv = new double*[numCurvesFrequency_];
    for(i=0; i<numCurvesFrequency_; i++) faccurv[i] = new double[mpcurv];

    // read the curves
    if(bifrdcurv)
    {
        for(i=0; i<numCurvesFrequency_; i++)
        {
            nterms = 0;
            while(!nterms){nterms = NC_GetLine(inputFile_, chinpline, chterms);}

            if(nterms != 2) NC_Error_Exit_2(NCout,
                                            "Number of terms of the firs input line for a curve must be 2!",
                                            "Number of terms = ", nterms, "number of the curve = ", i);
            nucurv[i] = NC_String2Integer(chterms[0]);
            npcurv[i] = NC_String2Integer(chterms[1]);

            for(j=0; j<npcurv[i]; j++)
            {
                while(!NC_GetLine(inputFile_, chinpline, chterms));
                frqcurv[i][j] = NC_String2Double(chterms[0]);
                faccurv[i][j] = NC_String2Double(chterms[1]);
            }
        }
    }

    // read the key word POST PROCESS
    if(bifrdcurv)
    {
        while(!NC_GetLine(inputFile_, chinpline, chterms));
        if(chterms[0].compare("POST")) NC_Error_Exit_0(NCout, "Key word POST PROCESS expected!");
    }

    // read the remainder input lines up to the key word END
    nterms = 0;
    while(nterms == 0 || nterms > 1 || (nterms == 1 && chterms[0].compare("END")))
    {
        nterms = NC_GetLine(inputFile_, chinpline, chterms);
    }
}

// store the boundary conditions
void NC_StoreBoundaryConditions
(
	ofstream& NCout,
	FILE* inputFile_,				// input file
	char* chinpline,			// input line
	string chterms[NTRM],		// entries in the input line
	int& nbcline,				//I: actual number of input lines
	inputLineBoundaryCondition *bouconlin	//O: array to store the input lines
)
{
    int i, i1, j, j1, k, nbtyp = 0, nrearefp = 0, nimgrefp = 0, nrearefn = 0, nimgrefn = 0,
    maxnuel = 0;
    string keywd;
    bool ifimpe, ifnega, ifref;
    double reap, rean = 0.0, imgp = 0.0, imgn = 0.0, dwk;
    Vector<int> ifimpecurve(numCurvesFrequency_);
    if(numCurvesFrequency_ > 0) for(i=0; i<numCurvesFrequency_; i++) ifimpecurve[i] = 0;
    int *ElGlobToLoca;


    // compute the vector indicating for each global element number the corresponding local number
    for(i=0; i<numElements_; i++) if(extNumbersOfElements[i] > maxnuel) maxnuel = extNumbersOfElements[i];
    if(maxnuel > 10000000)
    {
        NCout << "\nWarning: external element number too large!" << endl;
        NCout << "The maximun value = " << maxnuel << endl;
        cout << "\nWarning: external element number too large!" << endl;
        cout << "The maximun value = " << maxnuel << endl;
    }
    maxnuel++;
    ElGlobToLoca = new int[maxnuel];
    for(i=0; i<maxnuel; i++) ElGlobToLoca[i] = -1;

    for(i=0; i<numElements_; i++)
    {
        if(listElementProperty[i] == 2) continue;

        i1 = extNumbersOfElements[i];
        for(j=0; j<numElements_; j++)
        {
            if(listElementProperty[j] == 2) continue;
            if(extNumbersOfElements[j] == i1)
            {
                ElGlobToLoca[i1] = j;
                break;
            }
        }
    }

    // initialize the arrays
    for(i=0; i<numElements_; i++)
    {
        ibval[i] = 0;
        for(j=0; j<NNODPE; j++) zbvao0[i][j].set(0.0, 0.0);
        zbvao1[i].set(0.0, 0.0);
        zbvao2[i].set(0.0, 0.0);
    }

    // loop over the input lines
    for(i=0; i<nbcline; i++)
    {
        keywd = bouconlin[i].sKeyword;
        ifimpe = false;

        if(keywd == "VELO")
        {
            nbtyp = 0;
        }
        else if(keywd == "PRES")
        {
            nbtyp = 1;
        }
        else if(keywd == "ADMI")
        {
            nbtyp = 2;
        }
        else if(keywd == "IMPE")
        {
            nbtyp = 2;
            ifimpe = true;
        }
        else if(keywd == "TRAD")
        {
            nbtyp = 3;
            NC_Error_Exit_0(NCout, "Transfer admittance boundary condition is not supported!");
        }
        else if(keywd == "TRIM")
        {
            nbtyp = 3;
            ifimpe = true;
            NC_Error_Exit_0(NCout, "Transfer impedance boundary condition is not supported!");
        }
        else
        {
            cout << "The false key word = " << keywd << endl;
            NCout << "The false key word = " << keywd << endl;
            NC_Error_Exit_1(NCout, "False key word in a input line for boundary conditions!",
                            "Input line number = ", i + 1);
        }

        // shorthands of the structure entries
        ifref = bouconlin[i].bRef;
        ifnega = bouconlin[i].bNega;
        reap = bouconlin[i].dReal;
        imgp = bouconlin[i].dImag;

	if(ifimpe) {
	  if(reap == 0 && imgp == 0) {
	    ifimpe = false;
	    nbtyp = 1;  // sound soft PRES condition
	  }
	}

	if(ifref && ifimpe) {
	  cerr << "Sorry Impedance boundary conditions in CURVES is not implemented yet. Please use ADMI conditions. They work with CURVES. Sorry for the inconvinience.\n";
	  exit(-1);
	}

        if(ifimpe)
        {
            dwk = reap*reap + imgp*imgp;
            reap /= dwk;
            imgp /= (-dwk);
        }
        if(ifref)
        {
            nrearefp = bouconlin[i].nRealRef;
            nimgrefp = bouconlin[i].nImagRef;
        }
        if(ifnega)
        {
            rean = bouconlin[i].dNegReal;
            imgn = bouconlin[i].dNegImag;
            if(ifimpe)
            {
                dwk = rean*rean + imgn*imgn;
                rean /= dwk;
                imgn /= (-dwk);
            }
            if(ifref)
            {
                nrearefn = bouconlin[i].nNegRealRef;
                nimgrefn = bouconlin[i].nNegImagRef;
            }
        }

        // loop over the elements
        for(j=bouconlin[i].nLow; j<=bouconlin[i].nHigh; j++)
        {
            j1 = ElGlobToLoca[j]; // local number of the element
            if(j1 == -1) continue; // the element is a evaluation element
            switch(nbtyp)
            {
                case 0:
                    for(k=0; k<listNumberNodesPerElement[j1]; k++) zbvao0[j1][k].set(reap, imgp);
                    break;
                case 1:
                    for(k=0; k<listNumberNodesPerElement[j1]; k++) zbvao0[j1][k].set(reap, imgp);
                    ibval[j1] = nbtyp;
                    break;
                case 2:
                    zbvao1[j1].set(reap, imgp);
                    if(ifnega) zbvao2[j1].set(rean, imgn);
                    if(ifref)
                    {
                        ibvcurv[j1][2] = nrearefp;
                        ibvcurv[j1][3] = nimgrefp;
                        if(ifnega)
                        {
                            ibvcurv[j1][4] = nrearefn;
                            ibvcurv[j1][5] = nimgrefn;
                        }
                    }
                    switch(ibval[j1])
                {
                    case 0:
                        ibval[j1] = nbtyp;
                        break;
                    case 1:
                        cout << "Warning; for el. " << j << " the pressure boundary condition is " <<
                        "prescribed, the additional admittance boundary condition is ignored!" <<
                        endl;
                        break;
                    case 3:
                        ibval[j1] = 4;
                        break;
                    case 2:
                    case 4:
                        break;
                }
                    break;
                case 3:
                    for(k=0; k<listNumberNodesPerElement[j1]; k++) zbvao0[j1][k].set(reap, imgp);
                    if(ifref)
                    {
                        ibvcurv[j1][0] = nrearefp;
                        ibvcurv[j1][1] = nimgrefp;
                    }
                    if(ibval[j1] <= 1)
                    {
                        ibval[j1] = nbtyp;
                    }
                    else if(ibval[j1] == 2)
                    {
                        ibval[j1] = 4;
                    }
                    break;
            }
        } // end of loop J

        // update the array IFIMPECURVE
        if(ifimpe  && ifref)
        {
            ifimpecurve[nrearefp] = 1;
            ifimpecurve[nimgrefp] = 1;
            if(ifnega)
            {
                ifimpecurve[nrearefn] = 1;
                ifimpecurve[nimgrefn] = 1;
            }
        }
    } // end of loop I

    // convert the factors for impedance into that for admittance
    for(i=0; i<numCurvesFrequency_; i++)
    {
        if(ifimpecurve[i])
        {
            for(j=0; j<npcurv[i]; j++)
            {
                if(fabs(faccurv[i][j]) > 1.0e-10)
                {
                    faccurv[i][j] = 1.0/faccurv[i][j];
                }
                else
                {
                    faccurv[i][j] = faccurv[i][j] >= 0.0 ? 1.0e10 : -1.0e10;
                }
            }
        }
    }

    delete [] ElGlobToLoca;
}

// write the job parameters into the output file
void NC_WriteParameters
(
	ofstream& NCout,
	double* Freqs
 )
{
    int i, j, j1, k, k1, numFrequencies_ori = numFrequencies_;
    double freq_max, dmesh_ave= 0.0, dmesh_eqi = 0.0, dmesh_max = 0, dmesh_min = 1.0e20,
    dist, dis_ref;

    // compute the maximum, minimum and average distance between two neighbouring nodes
    k1 = 0;
    for(i=0; i<numElements_; i++)
    {
        j1 = listNumberNodesPerElement[i];

        if(j1 != NETYP3 && j1 != NETYP4) NC_Error_Exit_1(NCout, "An element must have 3 or 4 nodes!",
                                                         "Number of the false element = ", extNumbersOfElements[i]);

        if(listElementProperty[i] == 2) continue;

        if(j1 > NETYP4) j1 /= 2;
        if(areael[i] < 0.0) NC_Error_Exit_2(NCout, "Element with negative area is found!",
                                            "external elemen number = ", extNumbersOfElements[i], "internal element number = ", i);
        dis_ref = sqrt(areael[i]);
        if(j1 == NETYP3) dis_ref *= 1.51967;
        for(j=0; j<j1; j++)
        {
            k = j + 1;
            if(k == j1) k = 0;
            dist = Dispoi_dim3_(nodesCoordinates[elementsConnectivity[i][j]], nodesCoordinates[elementsConnectivity[i][k]]);
            if(dist > dmesh_max) dmesh_max = dist;
            if(dist < dmesh_min) dmesh_min = dist;
            dmesh_ave += dist;
            if(j1 == NETYP3) {dmesh_eqi += dist*0.748;} else {dmesh_eqi += dist;}
            if(dist > dis_ref*10.0 && dist < dis_ref*100.0) NC_Error_Warning_2(NCout,
                                                                               "Element with very long side found!", "external elemen number = ",
                                                                               extNumbersOfElements[i], "internal element number = ", i);
            if(dist > dis_ref*100.0) NC_Error_Exit_2(NCout,
                                                     "Element of very distorted form is found!", "external elemen number = ",
                                                     extNumbersOfElements[i], "internal element number = ", i);
            k1++;
        }
    }
    dmesh_ave /= (double)k1;
    dmesh_eqi /= (double)k1;

    delta_ = dmesh_min/100.0;

    // check the frequncies: the frequecies higher than the highest allowable frequency are ignored

    double ele_leng = dmesh_eqi, elnumperwve;

    elnumperwve = 6;

    freq_max = speedOfSound_/(elnumperwve*ele_leng);

    for(i=0; i<numFrequencies_; i++)
    {
        if(Freqs[i] > freq_max) break;
    }

    //if(i < numFrequencies_) numFrequencies_ = i; Changed by Harald Ziegelwanger to calculate all frequencies

    //if(numFrequencies_ori > numFrequencies_)
    if(numFrequencies_ori > i)
    {
        NCout << endl;
        NC_Error_Warning_2(NCout, "There are some frequencies higher than the highest allowable frequency!",
			   "Number of frequencies outside the range= ", numFrequencies_ori - i,
                           "The highest allowable frequency = ", freq_max);
    }

    // output the results
    NCout << versionNumber_ << endl;
    NCout << jobTitle_ << endl;

    NCout << endl;
    NCout << "Number of frequency steps ................... = " << numFrequencies_ << endl;

    NCout << endl;
    NCout << "First step .................................. = " << firstFrequencyStep_ << endl;
    NCout << "Number of steps ............................. = " << numFrequencies_ << endl;
    NCout << "Step increment .............................. = " << frequencyStepIncrement_ << endl;
    NCout << "Highest allowable frequency .................. = " << freq_max << endl;

    // write the basic BEM parameters
    NCout << "\n\n\n>> B A S I C   P A R A M E T E R S   F O R   B E M <<\n" << endl;

    NCout << "Number of element groups .................... = " << numElementGroups_ << endl;
    NCout << "Total number of elements .................... = " << numElements_ << endl;
    NCout << "Number of boundary elements ................. = " << numElementsOfBoundaryMesh_ << endl;
    NCout << "Number of evaluation elements ............... = " << numElementsOfEvaluationMesh_ << endl;
    NCout << "Maximum edge length of boundary elements .... = " << dmesh_max << endl;
    NCout << "Average edge length of boundary elements .... = " << dmesh_ave << endl;
    NCout << "Minimum edge length of boundary elements .... = " << dmesh_min << endl;

    NCout << endl;
    NCout << "Total number of nodes ....................... = " << numNodes_ << endl;
    NCout << "Number of nodes of the boundary elements .... = " << numNodesOfBoundaryMesh_ << endl;
    NCout << "Number of nodes of the evaluation elements .. = " << numNodesOfEvaluationMesh_ << endl;

    NCout << endl;
    NCout << "Number of symmetric planes .................. = " << numSymmetricPlanes_ << endl;

    NCout << endl;
    NCout << "Number of incident plane waves .............. = " << numIncidentPlaneWaves_ << endl;
    NCout << "Number of incident spharic waves ............ = " << numPointSources_ << endl;

    NCout << endl;
    NCout << "Sound speed ................................. = " << speedOfSound_ << endl;
    NCout << "Density of the medium ....................... = " << densityOfMedium_ <<
    endl;

    NCout << endl;
    NCout << "Diameter of the BE-mesh ..................... = " << diameterBE_ << endl;

    if(methodBEM_ > 0) {
        NCout << endl;
        NCout << "Minimum number of terms in the FMM expansion. = " << minExpansionTermsFMM_
        << endl;
    }

    NCout << endl;

    if(isInternalProblem_ == 0) NCout << "Exterior domain problem" << endl;
    if(isInternalProblem_ == 1) NCout << "Interior domain problem" << endl;

    NCout << endl;
    if(harmonicTimeFactor_ == 1.0)
    {
        NCout << "Harmonic time factor = EXP(-i*omega*t)" << endl;
    }
    else
    {
        NCout << "Harmonic time factor = EXP(i*omega*t)" << endl;
    }

    NCout << endl;

    switch(methodSolver_)
    {
        case 0:
            NCout << "methodSolver_  = 0: CGS iterative solver" << endl;
            break;
        case 4:
            NCout << "methodSolver_  = 4: Gaussean elimination (direct) solver" << endl;
            break;
        default:
            NC_Error_Exit_0(NCout, "methodSolver_ must be equal to 0, 1, 2, 3 or 4!");
    }
    if(methodSolver_ < 4)
    {
        switch(methodPreconditioner_)
        {
            case 0:
                NCout << "methodPreconditioner_ = 0: incomplete LU-decomposition preconditioner"
                << endl;
                break;
            case 1:
                NCout << "methodPreconditioner_ = 1: scaling preconditioner" << endl;
                break;
            case 2:
                NCout << "methodPreconditioner_ = 2: no preconditioner" << endl;
                break;
            default:
                NC_Error_Exit_0(NCout, "methodPreconditioner_ must be equal to 0, 1 or 2!");
        }
    }

    // write the frequency informations
    double ell = sqrt(averageElementArea_);

    NCout << endl;
    NCout << " No.   Freq.   Lambda/L_BE_Mesh  Lambda/L_Element" << endl;

    for(i=0; i<numFrequencies_; i++)
    {
        double wavl = speedOfSound_/Freqs[i];
        NCout << setw(3) << i+1 << " " << setw(7) << Freqs[i]
        << "   " << setw(10) << wavl/diameterBE_
        << "        " << setw(10) << wavl/ell << endl;
    }
    NCout << endl;
}

// get a input line and store the results in a string array
int NC_GetLine
(
	FILE *inputFile_,			// I: input file
	char *chinpline,		// W: input line
	string chterms[NTRM]	// O: array contain the terms in the input line
)
{
	int ii, len, nfi, ss, nn;

	if(fgets(chinpline, SIZE, inputFile_)==NULL) {
	  cout << "\nError parsing the input file\n";
	  exit(-1);
	}
	string inpline(chinpline);

	nn=0,ss=0,ii=0;
	len = (int)inpline.length()-1;
	while (nn < len)
	{
		nn = int(inpline.find(' ',ss));

		if (nn==-1)
		{
			nn = len;
		}
		if (nn-ss > 0)
		{
			chterms[ii++] = inpline.substr(ss, nn-ss);
		}
		ss = nn+1;
	}

	nfi = int(chterms[0].find("#",0));
	if(nfi >= 0) ii = 0;
	return(ii);
}

// change a string into integer
/* This is a bit over the top because c++ offers nicer ways */
int NC_String2Integer
(
	string Ori_Str		// I: string to be coverted
)
{
	int j = 0;
	char chwk[15];
	string chnu = " ";

	if((int)Ori_Str.length() > 15)
	{
		cout << "Length of the string > n_chwk!" << endl;
		cout << "The string: " << Ori_Str << endl;
		exit(1);
	}

	for(j=0; j<15; j++) chwk[j] = chnu[0];
	Ori_Str.copy(chwk, 15, 0);

	return(atoi(chwk));
}

// change a string into double
// kreiza changed everything about because valgrind complained
double NC_String2Double
(
	string Ori_Str		// I: string to be coverted
)
{
  /*	int j = 0;
	double dwk = 0.0;
	char chwk[25] = " ";
	string chnu = " ";*/

	if((int)Ori_Str.length() > 25)
	{
		cout << "Length of the string > n_chwk!" << endl;
		cout << "The string: " << Ori_Str << endl;
		exit(1);
	}

	return atof(Ori_Str.c_str());
	
	/*	for(j=0; j<25; j++) chwk[j] = chnu[0];
	Ori_Str.copy(chwk, 25, 0);
	sscanf(chwk, "%lf", &dwk);

	return(dwk);*/
}

// compute the frequencies
void NC_ComputeFrequencies
(
	ofstream& NCout,
	int& npoicurv,			//I: number of points of the frequency curve
	double* p_dcurvfrq,		//I: frequncy curve
	double* p_dcurvtime,	//I: time curve
	double* Freqs			//O: array of the frequencies
)
{
	int i, j, j1;
	double t;

	for(i=0; i<numFrequencies_; i++)
	{
		t = firstFrequencyStep_ + (double)(i + 1)*frequencyStepIncrement_;
		if(p_dcurvtime[0] >= t)
		{
			Freqs[i] = p_dcurvfrq[0];
			continue;
		}
		j1 = -1;
		for(j=1; j<npoicurv; j++)
		{
			if(p_dcurvtime[j] > t)
			{
				j1 = j;
				break;
			}
			else if(p_dcurvtime[j] == t)
			{
				Freqs[i] = p_dcurvfrq[j];
				j1 = npoicurv;
				break;
			}
		}
		if(j1 == -1)
		{
			Freqs[i] = p_dcurvfrq[npoicurv - 1];
		}
		else if(j1 < npoicurv)
		{
			Freqs[i] = p_dcurvfrq[j1 - 1] +
				(t - p_dcurvtime[j1 - 1])*(p_dcurvfrq[j1] - p_dcurvfrq[j1 - 1])/
				(p_dcurvtime[j1] - p_dcurvtime[j1 - 1]);
		}
	}
}
