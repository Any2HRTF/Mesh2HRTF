
/*===========================================================================*\
*                                                                            *
*                              OpenFlipper                                   *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
*                                                                            *
\*===========================================================================*/

#include <pybind11/include/pybind11/pybind11.h>
#include <pybind11/include/pybind11/embed.h>


#include <MeshGradingPlugin.hh>
#include <QString>
#include <QChar>

#include <OpenFlipper/BasePlugin/PythonFunctions.hh>
#include <OpenFlipper/PythonInterpreter/PythonTypeConversions.hh>

namespace py = pybind11;



PYBIND11_EMBEDDED_MODULE(MeshGrading, m) {

  QObject* pluginPointer = getPluginPointer("MeshGrading");

  if (!pluginPointer) {
     std::cerr << "Error Getting plugin pointer for Plugin-MeshGrading" << std::endl;
     return;
   }

  MeshGradingPlugin* plugin = qobject_cast<MeshGradingPlugin*>(pluginPointer);

  if (!plugin) {
    std::cerr << "Error converting plugin pointer for Plugin-MeshGrading" << std::endl;
    return;
  }

  // Export our core. Make sure that the c++ worlds core object is not deleted if
  // the python side gets deleted!!
  py::class_< MeshGradingPlugin,std::unique_ptr<MeshGradingPlugin, py::nodelete> > meshGrading(m, "MeshGrading");

  // On the c++ side we will just return the existing core instance
  // and prevent the system to recreate a new core as we need
  // to work on the existing one.
  meshGrading.def(py::init([plugin]() { return plugin; }));

  meshGrading.def("meshGrading", &MeshGradingPlugin::meshGrading,
                                 QCoreApplication::translate("PythonDocMeshGrading","Pre-processing for Mesh2HRTF (a-priori mesh-grading).").toLatin1().data(),
                                 py::arg(QCoreApplication::translate("PythonDocMeshGrading","ID of the object").toLatin1().data()),
                                 py::arg(QCoreApplication::translate("PythonDocMeshGrading","min edge length").toLatin1().data()),
                                 py::arg(QCoreApplication::translate("PythonDocMeshGrading","max edge length").toLatin1().data()),
                                 py::arg(QCoreApplication::translate("PythonDocMeshGrading","ear (left=1, right=2, both=3)").toLatin1().data()),
                                 py::arg(QCoreApplication::translate("PythonDocMeshGrading","gradingFunction").toLatin1().data()),
                                 py::arg(QCoreApplication::translate("PythonDocMeshGrading","gradingOrder").toLatin1().data())
                  );


}

