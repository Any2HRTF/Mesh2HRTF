/*===========================================================================*\
*                                                                            *
*                                 Mesh2HRTF                                  *
*                Copyright (C) 2015 by Harald Ziegelwanger,                  *
*        Acoustics Research Institute, Austrian Academy of Sciences          *
*                        mesh2hrtf.sourceforge.net                           *
*                                                                            *
*--------------------------------------------------------------------------- *
*                                                                            *
*  Mesh2HRTF is licensed under the GNU Lesser General Public License as      *
*  published by the Free Software Foundation, either version 3 of            *
*  the License, or (at your option) any later version.                       *
*                                                                            *
*  Mesh2HRTF is distributed in the hope that it will be useful,              *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
*  GNU Lesser General Public License for more details.                       *
*                                                                            *
*  You should have received a copy of the GNU LesserGeneral Public           *
*  License along with Mesh2HRTF. If not, see                                 *
*  <http://www.gnu.org/licenses/lgpl.html>.                                  *
*                                                                            *
*  If you use Mesh2HRTF:                                                     *
*  - Provide credits:                                                        *
*    "Mesh2HRTF, H. Ziegelwanger, ARI, OEAW (mesh2hrtf.sourceforge.net)"     *
*  - In your publication, cite both articles:                                *
*    [1] Ziegelwanger, H., Kreuzer, W., and Majdak, P. (2015). "Mesh2HRTF:   *
*        Open-source software package for the numerical calculation of       *
*        head-related transfer functions," in Proceedings of the 22nd        *
*        ICSV, Florence, IT.                                                 *
*    [2] Ziegelwanger, H., Majdak, P., and Kreuzer, W. (2015). "Numerical    *
*        calculation of listener-specific head-related transfer functions    *
*        and sound localization: Microphone model and mesh discretization,"  *
*        The Journal of the Acoustical Society of America, 138, 208-222.     *
*                                                                            *
*  If you use Plugin-MeshGrading:                                            *
*  - In your publication, cite:                                              *
*    [3] Ziegelwanger, H., and Majdak, P., Kreuzer, W. (submitted).          *
*        "A-priori mesh grading for the numerical calculation of the         *
*        head-related transfer functions," Applied Acoustics, , -.           *
*                                                                            *
*============================================================================*
*                                                                            *
*  This file is based on "Plugin-IsotropicRemesher" (OperFlipper)            *
*                                                                            *
*============================================================================*
*                                                                            *
*                              OpenFlipper                                   *
*      Copyright (C) 2001-2014 by Computer Graphics Group, RWTH Aachen       *
*                           www.openflipper.org                              *
*                                                                            *
*--------------------------------------------------------------------------- *
*  OpenFlipper is free software: you can redistribute it and/or modify       *
*  it under the terms of the GNU Lesser General Public License as            *
*  published by the Free Software Foundation, either version 3 of            *
*  the License, or (at your option) any later version with the               *
*  following exceptions:                                                     *
*                                                                            *
*  If other files instantiate templates or use macros                        *
*  or inline functions from this file, or you compile this file and          *
*  link it with other files to produce an executable, this file does         *
*  not by itself cause the resulting executable to be covered by the         *
*  GNU Lesser General Public License. This exception does not however        *
*  invalidate any other reasons why the executable file might be             *
*  covered by the GNU Lesser General Public License.                         *
*                                                                            *
*  OpenFlipper is distributed in the hope that it will be useful,            *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*  GNU Lesser General Public License for more details.                       *
*                                                                            *
\*===========================================================================*/

#include "MeshGradingPlugin.hh"
#include "MeshGradingT.hh"

#include <OpenFlipper/BasePlugin/PluginFunctions.hh>

MeshGradingPlugin::MeshGradingPlugin() :
tool_(0),
toolIcon_(0),
minEdgeLengthGlobal_(1.0),
maxEdgeLengthGlobal_(1.0),
ear_(1),
gradingFunction_(1),
gradingOrder_(1),
thread_(0)
{
}

/// init the Toolbox
void MeshGradingPlugin::initializePlugin() {
  if ( OpenFlipper::Options::gui()) {
    tool_ = new MeshGradingToolBox();

    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->meshGradingButton, SIGNAL(clicked()), this, SLOT(slotMeshGradingButtonClicked()) );

    toolIcon_ = new QIcon(OpenFlipper::Options::iconDirStr()+OpenFlipper::Options::dirSeparator()+"MeshGrading.png");
    emit addToolbox( tr("Mesh Grading") , tool_, toolIcon_ );
  }
}

void MeshGradingPlugin::slotMeshGradingButtonClicked() {
  minEdgeLengthGlobal_ = tool_->minEdgeLength->value();
  maxEdgeLengthGlobal_ = tool_->maxEdgeLength->value();
  gradingOrder_ = tool_->gradingOrder->value();

  if (tool_->ear_left->isChecked()){
      ear_ = 1;
  }
  if (tool_->ear_right->isChecked()){
      ear_ = 2;
  }
  if (tool_->ear_both->isChecked()){
      ear_ = 3;
  }

  if (tool_->POWalpha->isChecked()){
      gradingFunction_ = 1;
  }
  if (tool_->COSalpha->isChecked()){
      gradingFunction_ = 2;
  }

  if ( thread_ == 0){
    thread_ = new OpenFlipperThread(name() + "MeshGrading");                         // Create your thread containing a unique id \n
    connect(thread_,SIGNAL( finished(QString)), this,SIGNAL(finishJob(QString)));                           // connect your threads finish info to the global one ( you can do the same for a local one ) \n
    connect(thread_,SIGNAL( function() ), this,SLOT(slotMeshGrading()),Qt::DirectConnection);           // You can directly give a slot of your app that gets called \n
    connect(this,SIGNAL( finishJob(QString)), this, SLOT(threadFinished(QString)), Qt::QueuedConnection);
  }

  emit startJob( name() + "MeshGrading", "Mesh pre-processing for Mesh2HRTF (mesh grading)." , 0 , 100 , true);  // As long as meshes cannot be locked, this thread has to be blocking. Otherwise, during operation the mesh could be deleted. We don't want that!

  thread_->start();                                                                                       // start thread
  thread_->startProcessing();                                                                             // start processing
}



void MeshGradingPlugin::slotMeshGrading(){

  //read one target objects
  for ( PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS,DataType(DATA_TRIANGLE_MESH | DATA_POLY_MESH)) ;
                                        o_it != PluginFunctions::objectsEnd(); ++o_it)  {

    //check dataType
    if ( o_it->dataType(DATA_TRIANGLE_MESH) ) {
      TriMesh* mesh = PluginFunctions::triMesh(o_it);

      ProgressEmitter prgEmt(name() + "MeshGrading");
      connect (&prgEmt, SIGNAL(changeDescription(QString,QString)), this, SIGNAL(setJobDescription(QString,QString)) );
      connect (&prgEmt, SIGNAL(signalJobState(QString,int)), this, SIGNAL(setJobState(QString,int)) );
      MeshGrading< TriMesh > meshGrader(&prgEmt);

      meshGrader.grading(*mesh, minEdgeLengthGlobal_, maxEdgeLengthGlobal_, ear_, gradingFunction_, gradingOrder_);

      mesh->update_normals();

    }else{
      emit log("Mesh grading currently only implemented for triangle Meshes");
    }
  }
}

void MeshGradingPlugin::threadFinished(QString /*_jobId*/) {

  std::cerr << "threadFinished() called" << std::endl;

  for ( PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS,DataType(DATA_TRIANGLE_MESH | DATA_POLY_MESH)) ;
                                        o_it != PluginFunctions::objectsEnd(); ++o_it)  {

      emit updatedObject( o_it->id(), UPDATE_ALL );

      emit createBackup(o_it->id(),"Mesh Grading");

      emit updateView();

  }
}


/// Initialize the plugin
void MeshGradingPlugin::pluginsInitialized(){

  emit setSlotDescription("MeshGrading(int,double)", "Mesh Grading",
                          QString("object_id,_minEdgeLengthGlobal,_maxEdgeLengthGlobal,_ear,_gradingFunction,_gradingOrder").split(","),
                          QString("id of an object, global min edge length, global max edge length, ear, grading function, grading order").split(","));
}


void MeshGradingPlugin::meshGrading(int _objectID, double _minEdgeLengthGlobal, double _maxEdgeLengthGlobal, unsigned int  _ear, unsigned int  _gradingFunction, double _gradingOrder){
  BaseObjectData* object  = 0;


  if ( PluginFunctions::getObject(_objectID, object) ){

    //check dataType
    if ( object->dataType(DATA_TRIANGLE_MESH)) {

      TriMesh* mesh = PluginFunctions::triMesh(object);

      MeshGrading< TriMesh > meshgrader;

      meshgrader.grading(*mesh, _minEdgeLengthGlobal, _maxEdgeLengthGlobal, _ear, _gradingFunction, _gradingOrder);

      mesh->update_normals();

      emit updatedObject( object->id(), UPDATE_ALL );

      emit scriptInfo("MeshGrading(" + QString::number(_objectID) + ", " + QString::number(_minEdgeLengthGlobal) + ", " + QString::number(_maxEdgeLengthGlobal) + ", " + QString::number(_ear) + ", " + QString::number(_gradingFunction) + ", " + QString::number(_gradingOrder) + ")");

      emit updateView();

      return;

    }else{
      emit log("Mesh grading currently only implemented for triangle Meshes");
      return;
    }
  }else{
    emit log("Unable to get object");
  }
}

#if QT_VERSION < 0x050000
  Q_EXPORT_PLUGIN2( MeshGradingPlugin , MeshGradingPlugin );
#endif
