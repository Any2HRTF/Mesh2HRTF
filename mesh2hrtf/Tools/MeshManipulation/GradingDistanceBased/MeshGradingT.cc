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


#define MESHGRADING_C

#include "MeshGradingT.hh"

#include <ACG/Geometry/Algorithms.hh>

// -------------------- BSP
#include <ACG/Geometry/bsp/TriangleBSPT.hh>

#define PI 3.14159265

/// do the mesh grading
template< class MeshT >
void MeshGrading< MeshT >::grading( MeshT& _mesh, const double _minEdgeLengthGlobal, const double _maxEdgeLengthGlobal, const unsigned int _ear, const unsigned int _gradingFunction, const double _gradingOrder )
{

  typename MeshT::Point earVertexLeft;
  earVertexLeft[0]=100.0;
  earVertexLeft[1]=0.0;
  earVertexLeft[2]=100.0;
  typename MeshT::Point earVertexRight;
  earVertexRight[0]=100.0;
  earVertexRight[1]=0.0;
  earVertexRight[2]=100.0;
  typename MeshT::VertexIter v_it;
  typename MeshT::VertexIter v_end = _mesh.vertices_end();

  if (_ear==1 || _ear==3)
  {
    for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
    {
      typename MeshT::Point currentVertex = _mesh.point(*v_it);
      if (currentVertex[1]>0)
      {
        if (sqrt(currentVertex[0]*currentVertex[0]+currentVertex[2]*currentVertex[2])<=sqrt(earVertexLeft[0]*earVertexLeft[0]+earVertexLeft[2]*earVertexLeft[2]))
        {
          earVertexLeft=currentVertex;
        }
      }
    }
  }
  if (_ear==2 || _ear==3)
  {
    for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it)
    {
      typename MeshT::Point currentVertex = _mesh.point(*v_it);
      if (currentVertex[1]<0)
      {
        if (sqrt(currentVertex[0]*currentVertex[0]+currentVertex[2]*currentVertex[2])<=sqrt(earVertexRight[0]*earVertexRight[0]+earVertexRight[2]*earVertexRight[2]))
        {
          earVertexRight=currentVertex;
        }
      }
    }
  }

  std::cout << earVertexLeft[0] << ' ' << earVertexLeft[1] << ' ' << earVertexLeft[2] << ' ' << std::endl;
  std::cout << earVertexRight[0] << ' ' << earVertexRight[1] << ' ' << earVertexRight[2] << ' ' << std::endl;


  MeshT meshCopy = _mesh;
  OpenMeshTriangleBSPT< MeshT >* triangleBSP = getTriangleBSP(meshCopy);

  for (int i=0; i < 100; i++){
    if (prgEmt_)
      prgEmt_->sendProgressSignal(i + 0);

    std::cout << "Iteration = " << i << std::endl;

    splitLongEdges(_mesh, _minEdgeLengthGlobal, _maxEdgeLengthGlobal, earVertexLeft, earVertexRight, _ear, _gradingFunction, _gradingOrder);
    if (prgEmt_)
      prgEmt_->sendProgressSignal(i + 0.2);

    collapseShortEdges(_mesh, _minEdgeLengthGlobal, _maxEdgeLengthGlobal, earVertexLeft, earVertexRight, _ear, _gradingFunction, _gradingOrder);
    if (prgEmt_)
      prgEmt_->sendProgressSignal(i + 0.4);

    equalizeValences(_mesh);
    if (prgEmt_)
      prgEmt_->sendProgressSignal(i + 0.6);

    tangentialRelaxation(_mesh);
    if (prgEmt_)
      prgEmt_->sendProgressSignal(i + 0.8);

    projectToSurface(_mesh, meshCopy, triangleBSP);
    if (prgEmt_)
      prgEmt_->sendProgressSignal(i + 1);
  }

  delete triangleBSP;
}

template< class MeshT >
OpenMeshTriangleBSPT< MeshT >* MeshGrading< MeshT >::getTriangleBSP(MeshT& _mesh)
{
  // create Triangle BSP
  OpenMeshTriangleBSPT< MeshT >* triangle_bsp = new OpenMeshTriangleBSPT< MeshT >( _mesh );

  // build Triangle BSP
  triangle_bsp->reserve(_mesh.n_faces());

  typename MeshT::FIter f_it  = _mesh.faces_begin();
  typename MeshT::FIter f_end = _mesh.faces_end();

  for (; f_it!=f_end; ++f_it)
    triangle_bsp->push_back( *f_it );

  triangle_bsp->build(10, 100);

  // return pointer to triangle bsp
  return triangle_bsp;
}
template< class MeshT>
double MeshGrading<MeshT>::compute_distance(const typename MeshT::Point vec_midpoint, const unsigned int _ear,const typename MeshT::Point earVertex, const typename MeshT::Point earVertexLeft, const typename MeshT::Point earVertexRight){
    double distance;

    if(_ear == 1 || _ear == 2){
        typename MeshT::Point vec_distance = vec_midpoint - earVertex;
        distance = sqrt(vec_distance.sqrnorm());
    }else{
        typename MeshT::Point vec_distance_left = vec_midpoint - earVertexLeft;
        typename MeshT::Point vec_distance_right = vec_midpoint - earVertexRight;
        const double distance_left = sqrt(vec_distance_left.sqrnorm());
        const double distance_right = sqrt(vec_distance_right.sqrnorm());
        if(distance_left < distance_right){
            distance = distance_left;
        }else{
            distance = distance_right;
        }
    }
    return distance;
}

/// performs edge splits until all edges are shorter than the threshold
template< class MeshT >
void MeshGrading< MeshT >::splitLongEdges( MeshT& _mesh, const double _minEdgeLengthGlobal, const double _maxEdgeLengthGlobal, const typename MeshT::Point earVertexLeft, const typename MeshT::Point earVertexRight, const unsigned int _ear, const unsigned int _gradingFunction, const double _gradingOrder )
{

  double maxEdgeLength = 0;
  double grading = 0;

  typename MeshT::Point earVertex;
  double headwidth;

  if(_ear == 1){
      headwidth = sqrt(earVertexLeft.sqrnorm())*2;
      earVertex = earVertexLeft;
  }
  if(_ear == 2){
      headwidth = sqrt(earVertexRight.sqrnorm())*2;
      earVertex = earVertexRight;
  }
  if(_ear == 3){
      headwidth = sqrt((earVertexLeft - earVertexRight).sqrnorm());
  }



  typename MeshT::EdgeIter e_it;
  typename MeshT::EdgeIter e_end = _mesh.edges_end();

  //iterate over all edges
  for (e_it = _mesh.edges_begin(); e_it != e_end; ++e_it){
    const typename MeshT::HalfedgeHandle & hh = _mesh.halfedge_handle( *e_it, 0 );

    const typename MeshT::VertexHandle & v0 = _mesh.from_vertex_handle(hh);
    const typename MeshT::VertexHandle & v1 = _mesh.to_vertex_handle(hh);

    const typename MeshT::Point vec_edge = _mesh.point(v1) - _mesh.point(v0);
    const typename MeshT::Point vec_midpoint = (_mesh.point(v1) + _mesh.point(v0))/2;

    const double distance = compute_distance(vec_midpoint, _ear, earVertex, earVertexLeft, earVertexRight);

    grading=distance/headwidth;
    
    switch(_gradingFunction)
    {
        case 1:
            grading=pow(grading,_gradingOrder);
            break;
        case 2:
            grading=1-pow(cos(PI*grading/2),_gradingOrder);
            break;
    }

    if (distance>=0 && distance<=headwidth){
      maxEdgeLength = (_minEdgeLengthGlobal+(_maxEdgeLengthGlobal-_minEdgeLengthGlobal)*grading)*1.05;
    }
    else if (distance>headwidth){
      maxEdgeLength=_maxEdgeLengthGlobal*1.05;
    }

    //edge to long?
    if ( vec_edge.sqrnorm() > maxEdgeLength*maxEdgeLength ){

      const typename MeshT::Point midPoint = _mesh.point(v0) + ( 0.5 * vec_edge );

      //split at midpoint
      typename MeshT::VertexHandle vh = _mesh.add_vertex( midPoint );

      bool hadFeature = _mesh.status(*e_it).feature();

      _mesh.split(*e_it, vh);

      if ( hadFeature ){

        typename MeshT::VOHIter vh_it;
        for (vh_it = _mesh.voh_iter(vh); vh_it.is_valid(); ++vh_it)
          if ( _mesh.to_vertex_handle(*vh_it) == v0 || _mesh.to_vertex_handle(*vh_it) == v1 )
            _mesh.status( _mesh.edge_handle( *vh_it ) ).set_feature( true );
      }
    }
  }
}

/// collapse edges shorter than minEdgeLength if collapsing doesn't result in new edge longer than maxEdgeLength
template< class MeshT >
void MeshGrading< MeshT >::collapseShortEdges( MeshT& _mesh, const double _minEdgeLengthGlobal, const double _maxEdgeLengthGlobal, const typename MeshT::Point earVertexLeft, const typename MeshT::Point earVertexRight, const unsigned int _ear, const unsigned int _gradingFunction, const double _gradingOrder )
{

  double minEdgeLength = 0;
  double maxEdgeLength = 0;
  double grading = 0;

  typename MeshT::Point earVertex;
  double headwidth;

  if(_ear == 1){
      headwidth = sqrt(earVertexLeft.sqrnorm())*2;
      earVertex = earVertexLeft;
  }
  if(_ear == 2){
      headwidth = sqrt(earVertexRight.sqrnorm())*2;
      earVertex = earVertexRight;
  }
  if(_ear == 3){
      headwidth = sqrt((earVertexLeft - earVertexRight).sqrnorm());
  }

  //add checked property
  OpenMesh::EPropHandleT< bool > checked;
  if ( !_mesh.get_property_handle(checked, "Checked Property") )
    _mesh.add_property(checked,"Checked Property" );

  //init property
  typename MeshT::ConstEdgeIter e_it;
  typename MeshT::ConstEdgeIter e_end = _mesh.edges_end();

  for (e_it = _mesh.edges_begin(); e_it != e_end; ++e_it)
    _mesh.property(checked, *e_it) = false;

  bool finished = false;

  while( !finished ){

    finished = true;

    for (e_it = _mesh.edges_begin(); e_it != _mesh.edges_end() ; ++e_it){

        if ( _mesh.property(checked, *e_it) )
          continue;

        _mesh.property(checked, *e_it) = true;

        const typename MeshT::HalfedgeHandle & hh = _mesh.halfedge_handle( *e_it, 0 );

        const typename MeshT::VertexHandle & v0 = _mesh.from_vertex_handle(hh);
        const typename MeshT::VertexHandle & v1 = _mesh.to_vertex_handle(hh);

        const typename MeshT::Point vec_edge = _mesh.point(v1) - _mesh.point(v0);
        const typename MeshT::Point vec_midpoint = (_mesh.point(v1) + _mesh.point(v0))/2;

        const double distance = compute_distance(vec_midpoint, _ear, earVertex, earVertexLeft, earVertexRight);

        const double edgeLength = vec_edge.sqrnorm();

        grading=distance/headwidth;

        switch(_gradingFunction)
        {
            case 1:
                grading=pow(grading,_gradingOrder);
                break;
            case 2:
                grading=1-pow(cos(PI*grading/2),_gradingOrder);
                break;
        }

        if (distance>0 && distance<=headwidth){
          minEdgeLength = (_minEdgeLengthGlobal+(_maxEdgeLengthGlobal-_minEdgeLengthGlobal)*grading)*0.95;
          maxEdgeLength = (_minEdgeLengthGlobal+(_maxEdgeLengthGlobal-_minEdgeLengthGlobal)*grading)*1.05;
        }
        else if (distance>headwidth){
          minEdgeLength = _maxEdgeLengthGlobal*0.95;
          maxEdgeLength = _maxEdgeLengthGlobal*1.05;
        }

        // edge too short but don't try to collapse edges that have length 0
        if ( (edgeLength < minEdgeLength*minEdgeLength) && (edgeLength > DBL_EPSILON) ){

          //check if the collapse is ok
          const typename MeshT::Point & B = _mesh.point(v1);

          bool collapse_ok = true;

          for(typename MeshT::VOHIter vh_it = _mesh.voh_iter(v0); vh_it.is_valid(); ++vh_it)
            if ( (( B - _mesh.point( _mesh.to_vertex_handle(*vh_it) ) ).sqrnorm() > maxEdgeLength*maxEdgeLength )
                 || ( _mesh.status( _mesh.edge_handle( *vh_it ) ).feature())
                 || ( _mesh.is_boundary( _mesh.edge_handle( *vh_it ) ) )){
              collapse_ok = false;
              break;
            }

          if( collapse_ok && _mesh.is_collapse_ok(hh) ) {
            _mesh.collapse( hh );

            finished = false;
          }

        }
    }

  }

  _mesh.remove_property(checked);

  _mesh.garbage_collection();
}

template< class MeshT >
void MeshGrading< MeshT >::equalizeValences( MeshT& _mesh )
{

  typename MeshT::EdgeIter e_it;
  typename MeshT::EdgeIter e_end = _mesh.edges_end();

  for (e_it = _mesh.edges_sbegin(); e_it != e_end; ++e_it){

    if ( !_mesh.is_flip_ok(*e_it) ) continue;
    if ( _mesh.status( *e_it ).feature() ) continue;


    const typename MeshT::HalfedgeHandle & h0 = _mesh.halfedge_handle( *e_it, 0 );
    const typename MeshT::HalfedgeHandle & h1 = _mesh.halfedge_handle( *e_it, 1 );

    if (h0.is_valid() && h1.is_valid())

      if (_mesh.face_handle(h0).is_valid() && _mesh.face_handle(h1).is_valid()){
        //get vertices of corresponding faces
        const typename MeshT::VertexHandle & a = _mesh.to_vertex_handle(h0);
        const typename MeshT::VertexHandle & b = _mesh.to_vertex_handle(h1);
        const typename MeshT::VertexHandle & c = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(h0));
        const typename MeshT::VertexHandle & d = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(h1));

        const int deviation_pre =  abs((int)(_mesh.valence(a) - targetValence(_mesh, a)))
                                  +abs((int)(_mesh.valence(b) - targetValence(_mesh, b)))
                                  +abs((int)(_mesh.valence(c) - targetValence(_mesh, c)))
                                  +abs((int)(_mesh.valence(d) - targetValence(_mesh, d)));
        _mesh.flip(*e_it);

        const int deviation_post = abs((int)(_mesh.valence(a) - targetValence(_mesh, a)))
                                  +abs((int)(_mesh.valence(b) - targetValence(_mesh, b)))
                                  +abs((int)(_mesh.valence(c) - targetValence(_mesh, c)))
                                  +abs((int)(_mesh.valence(d) - targetValence(_mesh, d)));

        if (deviation_pre <= deviation_post)
          _mesh.flip(*e_it);
      }
  }
}

///returns 4 for boundary vertices and 6 otherwise
template< class MeshT >
inline
int MeshGrading< MeshT >::targetValence( MeshT& _mesh, const typename MeshT::VertexHandle& _vh ){

  if (isBoundary(_mesh,_vh))
    return 4;
  else
    return 6;
}

template< class MeshT >
inline
bool MeshGrading< MeshT >::isBoundary( MeshT& _mesh, const typename MeshT::VertexHandle& _vh ){

  typename MeshT::ConstVertexOHalfedgeIter voh_it;

  for (voh_it = _mesh.voh_iter( _vh ); voh_it.is_valid(); ++voh_it )
    if ( _mesh.is_boundary( _mesh.edge_handle( *voh_it ) ) )
      return true;

  return false;
}

template< class MeshT >
inline
bool MeshGrading< MeshT >::isFeature( MeshT& _mesh, const typename MeshT::VertexHandle& _vh ){

  typename MeshT::ConstVertexOHalfedgeIter voh_it;

  for (voh_it = _mesh.voh_iter( _vh ); voh_it.is_valid(); ++voh_it )
    if ( _mesh.status( _mesh.edge_handle( *voh_it ) ).feature() )
      return true;

  return false;
}

template< class MeshT >
void MeshGrading< MeshT >::tangentialRelaxation( MeshT& _mesh )
{
  _mesh.update_normals();

  //add checked property
  OpenMesh::VPropHandleT< typename MeshT::Point > q;
  if ( !_mesh.get_property_handle(q, "q Property") )
    _mesh.add_property(q,"q Property" );

  typename MeshT::ConstVertexIter v_it;
  typename MeshT::ConstVertexIter v_end = _mesh.vertices_end();

  //first compute barycenters
  for (v_it = _mesh.vertices_sbegin(); v_it != v_end; ++v_it){

    typename MeshT::Point tmp(0.0, 0.0, 0.0);
    uint N = 0;

    typename MeshT::VVIter vv_it;
    for (vv_it = _mesh.vv_iter(*v_it); vv_it.is_valid(); ++vv_it){
      tmp += _mesh.point(*vv_it);
      N++;
    }

    if (N > 0)
      tmp /= (double) N;

    _mesh.property(q, *v_it) = tmp;
  }

  //move to new position
  for (v_it = _mesh.vertices_sbegin(); v_it != v_end; ++v_it){
    if ( !isBoundary(_mesh, *v_it) && !isFeature(_mesh, *v_it) )
      _mesh.set_point(*v_it,  _mesh.property(q, *v_it) + (_mesh.normal(*v_it)| (_mesh.point(*v_it) - _mesh.property(q, *v_it) ) ) * _mesh.normal(*v_it));
  }

  _mesh.remove_property(q);
}

template <class MeshT>
template <class SpatialSearchT>
typename MeshT::Point
MeshGrading< MeshT >::findNearestPoint(const MeshT&                   _mesh,
                                          const typename MeshT::Point&   _point,
                                          typename MeshT::FaceHandle&    _fh,
                                          SpatialSearchT*                _ssearch,
                                          double*                        _dbest)
{

  typename MeshT::Point  p_best = _mesh.point(_mesh.vertex_handle(0));
  typename MeshT::Scalar d_best = (_point-p_best).sqrnorm();

  typename MeshT::FaceHandle fh_best;

  if( _ssearch == 0 )
  {
    // exhaustive search
    typename MeshT::ConstFaceIter cf_it  = _mesh.faces_begin();
    typename MeshT::ConstFaceIter cf_end = _mesh.faces_end();

    for(; cf_it != cf_end; ++cf_it)
    {
      typename MeshT::ConstFaceVertexIter cfv_it = _mesh.cfv_iter(*cf_it);

      const typename MeshT::Point& pt0 = _mesh.point(   *cfv_it);
      const typename MeshT::Point& pt1 = _mesh.point( *(++cfv_it));
      const typename MeshT::Point& pt2 = _mesh.point( *(++cfv_it));

      typename MeshT::Point ptn;

      typename MeshT::Scalar d = ACG::Geometry::distPointTriangleSquared( _point,
                     pt0,
                     pt1,
                     pt2,
                     ptn );

      if( d < d_best)
      {
        d_best = d;
        p_best = ptn;

        fh_best = *cf_it;
      }
    }

    // return facehandle
    _fh = fh_best;

    // return distance
    if(_dbest)
      *_dbest = sqrt(d_best);


    return p_best;
  }
  else
  {
    typename MeshT::FaceHandle     fh = _ssearch->nearest(_point).handle;
    typename MeshT::CFVIter        fv_it = _mesh.cfv_iter(fh);

    const typename MeshT::Point&   pt0 = _mesh.point( *(  fv_it));
    const typename MeshT::Point&   pt1 = _mesh.point( *(++fv_it));
    const typename MeshT::Point&   pt2 = _mesh.point( *(++fv_it));

    // project
    d_best = ACG::Geometry::distPointTriangleSquared(_point, pt0, pt1, pt2, p_best);

    // return facehandle
    _fh = fh;

    // return distance
    if(_dbest)
      *_dbest = sqrt(d_best);

    return p_best;
  }
}


template< class MeshT >
template< class SpatialSearchT >
void MeshGrading< MeshT >::projectToSurface( MeshT& _mesh, MeshT& _original, SpatialSearchT* _ssearch )
{

  typename MeshT::VertexIter v_it;
  typename MeshT::VertexIter v_end = _mesh.vertices_end();

  for (v_it = _mesh.vertices_begin(); v_it != v_end; ++v_it){

    if (isBoundary(_mesh, *v_it)) continue;
    if ( isFeature(_mesh, *v_it)) continue;

    typename MeshT::Point p = _mesh.point(*v_it);
    typename MeshT::FaceHandle fhNear;
    double distance;

    typename MeshT::Point pNear = findNearestPoint(_original, p, fhNear, _ssearch, &distance);

    _mesh.set_point(*v_it, pNear);
  }
}
