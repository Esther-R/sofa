/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <SofaGeneralDeformable/TriangularBendingSprings.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/TopologyChange.h>
#include <sofa/type/RGBAColor.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging

#include <SofaBaseTopology/TopologyData.inl>

namespace sofa::component::forcefield
{

typedef core::topology::BaseMeshTopology::EdgesInTriangle EdgesInTriangle;

template< class DataTypes>
void TriangularBendingSprings<DataTypes>::applyEdgeCreation(Index , EdgeInformation &ei, const core::topology::Edge &,
    const sofa::type::vector<Index> &, const sofa::type::vector<double> &)
{
        unsigned int u,v;
        /// set to zero the edge stiffness matrix
        for (u=0; u<N; ++u)
        {
            for (v=0; v<N; ++v)
            {
                ei.DfDx[u][v]=0;
            }
        }

        ei.is_activated=false;
        ei.is_initialized=false;
}



template< class DataTypes>
void TriangularBendingSprings<DataTypes>::applyTriangleCreation(const sofa::type::vector<Index> &triangleAdded, const sofa::type::vector<core::topology::Triangle> &, 
    const sofa::type::vector<sofa::type::vector<Index> > &, const sofa::type::vector<sofa::type::vector<double> > &)
{
    using namespace core::topology;
        double m_ks=getKs();
        double m_kd=getKd();

        unsigned int u,v;

        unsigned int nb_activated = 0;

        const typename DataTypes::VecCoord& restPosition = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();

        type::vector<EdgeInformation>& edgeData = *(edgeInfo.beginEdit());

        for (unsigned int i=0; i<triangleAdded.size(); ++i)
        {

            /// describe the jth edge index of triangle no i
            EdgesInTriangle te2 = m_topology->getEdgesInTriangle(triangleAdded[i]);
            /// describe the jth vertex index of triangle no i
            Triangle t2 = m_topology->getTriangle(triangleAdded[i]);

            for(unsigned int j=0; j<3; ++j)
            {

                EdgeInformation &ei = edgeData[te2[j]]; // edgeInfo
                if(!(ei.is_initialized))
                {

                    unsigned int edgeIndex = te2[j];
                    ei.is_activated=true;

                    /// set to zero the edge stiffness matrix
                    for (u=0; u<N; ++u)
                    {
                        for (v=0; v<N; ++v)
                        {
                            ei.DfDx[u][v]=0;
                        }
                    }

                    const auto& shell = m_topology->getTrianglesAroundEdge(edgeIndex);
                    if (shell.size()==2)
                    {

                        nb_activated+=1;

                        EdgesInTriangle te1;
                        Triangle t1;

                        if(shell[0] == triangleAdded[i])
                        {

                            te1 = m_topology->getEdgesInTriangle(shell[1]);
                            t1 = m_topology->getTriangle(shell[1]);

                        }
                        else   // shell[1] == triangleAdded[i]
                        {

                            te1 = m_topology->getEdgesInTriangle(shell[0]);
                            t1 = m_topology->getTriangle(shell[0]);
                        }

                        int i1 = m_topology->getEdgeIndexInTriangle(te1, edgeIndex); //edgeIndex //te1[j]
                        int i2 = m_topology->getEdgeIndexInTriangle(te2, edgeIndex); // edgeIndex //te2[j]

                        ei.m1 = t1[i1];
                        ei.m2 = t2[i2];

                        //TriangularBendingSprings<DataTypes> *fftest= (TriangularBendingSprings<DataTypes> *)param;
                        ei.ks=m_ks; //(fftest->ks).getValue();
                        ei.kd=m_kd; //(fftest->kd).getValue();

                        Coord u = (restPosition)[ei.m1] - (restPosition)[ei.m2];

                        Real d = u.norm();

                        ei.restlength=(double) d;

                        ei.is_activated=true;

                    }
                    else
                    {

                        ei.is_activated=false;

                    }

                    ei.is_initialized = true;
                }
            }

        }
        edgeInfo.endEdit();    
}


template< class DataTypes>
void TriangularBendingSprings<DataTypes>::applyTriangleDestruction(const sofa::type::vector<Index> &triangleRemoved)
{
    using namespace core::topology;

        double m_ks=getKs(); // typename DataTypes::
        double m_kd=getKd(); // typename DataTypes::

        //unsigned int u,v;

        const typename DataTypes::VecCoord& restPosition = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
        type::vector<EdgeInformation>& edgeData = *(edgeInfo.beginEdit());

        for (unsigned int i=0; i<triangleRemoved.size(); ++i)
        {
            /// describe the jth edge index of triangle no i
            EdgesInTriangle te = m_topology->getEdgesInTriangle(triangleRemoved[i]);
            /// describe the jth vertex index of triangle no i
            //Triangle t = m_topology->getTriangle(triangleRemoved[i]);


            for(unsigned int j=0; j<3; ++j)
            {

                EdgeInformation &ei = edgeData[te[j]]; // edgeInfo
                if(ei.is_initialized)
                {

                    unsigned int edgeIndex = te[j];

                    const auto& shell = m_topology->getTrianglesAroundEdge(edgeIndex);
                    if (shell.size()==3)
                    {

                        EdgesInTriangle te1;
                        Triangle t1;
                        EdgesInTriangle te2;
                        Triangle t2;

                        if(shell[0] == triangleRemoved[i])
                        {
                            te1 = m_topology->getEdgesInTriangle(shell[1]);
                            t1 = m_topology->getTriangle(shell[1]);
                            te2 = m_topology->getEdgesInTriangle(shell[2]);
                            t2 = m_topology->getTriangle(shell[2]);

                        }
                        else
                        {

                            if(shell[1] == triangleRemoved[i])
                            {

                                te1 = m_topology->getEdgesInTriangle(shell[2]);
                                t1 = m_topology->getTriangle(shell[2]);
                                te2 = m_topology->getEdgesInTriangle(shell[0]);
                                t2 = m_topology->getTriangle(shell[0]);

                            }
                            else   // shell[2] == triangleRemoved[i]
                            {

                                te1 = m_topology->getEdgesInTriangle(shell[0]);
                                t1 = m_topology->getTriangle(shell[0]);
                                te2 = m_topology->getEdgesInTriangle(shell[1]);
                                t2 = m_topology->getTriangle(shell[1]);

                            }
                        }

                        int i1 = m_topology->getEdgeIndexInTriangle(te1, edgeIndex);
                        int i2 = m_topology->getEdgeIndexInTriangle(te2, edgeIndex);

                        ei.m1 = t1[i1];
                        ei.m2 = t2[i2];

                        //TriangularBendingSprings<DataTypes> *fftest= (TriangularBendingSprings<DataTypes> *)param;
                        ei.ks=m_ks; //(fftest->ks).getValue();
                        ei.kd=m_kd; //(fftest->kd).getValue();

                        Coord u = (restPosition)[ei.m1] - (restPosition)[ei.m2];
                        Real d = u.norm();

                        ei.restlength=(double) d;

                        ei.is_activated=true;

                    }
                    else
                    {

                        ei.is_activated=false;
                        ei.is_initialized = false;

                    }

                }
                else
                {

                    ei.is_activated=false;
                    ei.is_initialized = false;

                }
            }

        }

        edgeInfo.endEdit();
}


template<class DataTypes>
void TriangularBendingSprings<DataTypes>::applyPointDestruction(const sofa::type::vector<Index> &tab)
{
    using namespace core::topology;
        bool debug_mode = false;

        unsigned int last = m_topology->getNbPoints() -1;
        unsigned int i,j;

        type::vector<EdgeInformation>& edgeInf = *(edgeInfo.beginEdit());

        sofa::type::vector<unsigned int> lastIndexVec;
        for(unsigned int i_init = 0; i_init < tab.size(); ++i_init)
        {

            lastIndexVec.push_back(last - i_init);
        }

        for ( i = 0; i < tab.size(); ++i)
        {

            unsigned int i_next = i;
            bool is_reached = false;
            while( (!is_reached) && (i_next < lastIndexVec.size() - 1))
            {

                i_next += 1 ;
                is_reached = is_reached || (lastIndexVec[i_next] == tab[i]);
            }

            if(is_reached)
            {

                lastIndexVec[i_next] = lastIndexVec[i];

            }

            const auto &shell= m_topology->getTrianglesAroundVertex(lastIndexVec[i]);
            for (j=0; j<shell.size(); ++j)
            {

                Triangle tj = m_topology->getTriangle(shell[j]);

                int vertexIndex = m_topology->getVertexIndexInTriangle(tj, lastIndexVec[i]);

                EdgesInTriangle tej = m_topology->getEdgesInTriangle(shell[j]);

                unsigned int ind_j = tej[vertexIndex];

                if (edgeInf[ind_j].m1 == (int) last)
                {
                    edgeInf[ind_j].m1=(int) tab[i];
                }
                else
                {
                    if (edgeInf[ind_j].m2 == (int) last)
                    {
                        edgeInf[ind_j].m2=(int) tab[i];
                    }
                }
            }

            if(debug_mode)
            {

                for (unsigned int j_loc=0; j_loc<edgeInf.size(); ++j_loc)
                {

                    //bool is_forgotten = false;
                    if (edgeInf[j_loc].m1 == (int) last)
                    {
                        edgeInf[j_loc].m1 =(int) tab[i];
                        //is_forgotten=true;

                    }
                    else
                    {
                        if (edgeInf[j_loc].m2 ==(int) last)
                        {
                            edgeInf[j_loc].m2 =(int) tab[i];
                            //is_forgotten=true;

                        }

                    }

                }
            }

            --last;
        }

        edgeInfo.endEdit();
}


template<class DataTypes>
void TriangularBendingSprings<DataTypes>::applyPointRenumbering(const sofa::type::vector<Index> &tab)
{
        type::vector<EdgeInformation>& edgeInf = *(edgeInfo.beginEdit());
        for (unsigned int i = 0; i < m_topology->getNbEdges(); ++i)
        {
            if(edgeInf[i].is_activated)
            {
                edgeInf[i].m1  = tab[edgeInf[i].m1];
                edgeInf[i].m2  = tab[edgeInf[i].m2];
            }
        }
        edgeInfo.endEdit();
}


template<class DataTypes>
TriangularBendingSprings<DataTypes>::TriangularBendingSprings(/*double _ks, double _kd*/)
    : f_ks(initData(&f_ks,(double) 100000.0,"stiffness","uniform stiffness for the all springs")) //(Real)0.3 ??
    , f_kd(initData(&f_kd,(double) 1.0,"damping","uniform damping for the all springs")) // (Real)1000. ??
    , d_showSprings(initData(&d_showSprings, true, "showSprings", "option to draw springs"))
    , l_topology(initLink("topology", "link to the topology container"))
    , edgeInfo(initData(&edgeInfo, "edgeInfo", "Internal edge data"))
    , updateMatrix(true)
    , edgeHandler(nullptr)
    , m_topology(nullptr)
{
    // Create specific handler for EdgeData
    edgeHandler = new TriangularBSEdgeHandler(&edgeInfo);
    //msg_error()<<"TriangularBendingSprings<DataTypes>::TriangularBendingSprings";
}

template<class DataTypes>
TriangularBendingSprings<DataTypes>::~TriangularBendingSprings()
{
    if(edgeHandler) delete edgeHandler;
}


template<class DataTypes>
void TriangularBendingSprings<DataTypes>::init()
{
    //msg_error() << "initializing TriangularBendingSprings" ;
    this->Inherited::init();

    if (l_topology.empty())
    {
        msg_info() << "link to Topology container should be set to ensure right behavior. First Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    m_topology = l_topology.get();
    msg_info() << "Topology path used: '" << l_topology.getLinkedPath() << "'";

    if (m_topology == nullptr)
    {
        msg_error() << "No topology component found at path: " << l_topology.getLinkedPath() << ", nor in current context: " << this->getContext()->name;
        sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    if (m_topology->d_componentState.getValue() == sofa::core::objectmodel::ComponentState::Invalid || m_topology->getNbTriangles()==0)
    {
        msg_error() << " object must have a Triangular Set Topology.";
        sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }
    edgeInfo.createTopologyHandler(m_topology,edgeHandler);
    edgeInfo.linkToPointDataArray();
    edgeInfo.linkToTriangleDataArray();

    edgeInfo.setCreationCallback([this](Index edgeIndex, EdgeInformation& ei,
        const core::topology::BaseMeshTopology::Edge& edge,
        const sofa::type::vector< Index >& ancestors,
        const sofa::type::vector< double >& coefs)
    {
        applyEdgeCreation(edgeIndex, ei, edge, ancestors, coefs);
    });


    edgeHandler->addCallBack(sofa::core::topology::TopologyChangeType::TRIANGLESADDED, [this](const core::topology::TopologyChange* eventTopo) {
        const core::topology::TrianglesAdded* triAdd = static_cast<const core::topology::TrianglesAdded*>(eventTopo);
        applyTriangleCreation(triAdd->getIndexArray(), triAdd->getElementArray(), triAdd->ancestorsList, triAdd->coefs);
    });

    edgeHandler->addCallBack(sofa::core::topology::TopologyChangeType::TRIANGLESREMOVED, [this](const core::topology::TopologyChange* eventTopo) {
        const core::topology::TrianglesRemoved* triRemove = static_cast<const core::topology::TrianglesRemoved*>(eventTopo);
        applyTriangleDestruction(triRemove->getArray());
    });

    edgeHandler->addCallBack(sofa::core::topology::TopologyChangeType::POINTSREMOVED, [this](const core::topology::TopologyChange* eventTopo) {
        const core::topology::PointsRemoved* pRemove = static_cast<const core::topology::PointsRemoved*>(eventTopo);
        applyPointDestruction(pRemove->getArray());
    });

    edgeHandler->addCallBack(sofa::core::topology::TopologyChangeType::POINTSRENUMBERING, [this](const core::topology::TopologyChange* eventTopo) {
        const core::topology::PointsRenumbering* pRenum = static_cast<const core::topology::PointsRenumbering*>(eventTopo);
        applyPointRenumbering(pRenum->getIndexArray());
    });

    this->reinit();
}


template<class DataTypes>
void TriangularBendingSprings<DataTypes>::reinit()
{
    using namespace core::topology;
    /// prepare to store info in the edge array
    type::vector<EdgeInformation>& edgeInf = *(edgeInfo.beginEdit());
    edgeInf.resize(m_topology->getNbEdges());
    Index i;
    // set edge tensor to 0
    for (i=0; i<m_topology->getNbEdges(); ++i)
    {

        applyEdgeCreation(i, edgeInf[i],
            m_topology->getEdge(i),  (const sofa::type::vector< Index > )0,
            (const sofa::type::vector< double >)0);
    }

    // create edge tensor by calling the triangle creation function
    sofa::type::vector<Index> triangleAdded;
    for (i=0; i<m_topology->getNbTriangles(); ++i)
        triangleAdded.push_back(i);

    applyTriangleCreation(triangleAdded,
        (const sofa::type::vector<Triangle>)0,
        (const sofa::type::vector<sofa::type::vector<Index> >)0,
        (const sofa::type::vector<sofa::type::vector<double> >)0);

    edgeInfo.endEdit();
}

template <class DataTypes>
SReal TriangularBendingSprings<DataTypes>::getPotentialEnergy(const core::MechanicalParams* /* mparams */, const DataVecCoord& /* d_x */) const
{
    msg_error()<<"TriangularBendingSprings::getPotentialEnergy-not-implemented !!!";
    return 0;
}


template<class DataTypes>
void TriangularBendingSprings<DataTypes>::addForce(const core::MechanicalParams* /* mparams */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v)
{
    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    const VecDeriv& v = d_v.getValue();

    size_t nbEdges=m_topology->getNbEdges();
    EdgeInformation *einfo;
    type::vector<EdgeInformation>& edgeInf = *(edgeInfo.beginEdit());

    //const type::vector<Spring>& m_springs= this->springs.getValue();
    //this->dfdx.resize(nbEdges); //m_springs.size()
    f.resize(x.size());
    m_potentialEnergy = 0;
    /*        msg_error()<<"TriangularBendingSprings<DataTypes>::addForce()";*/

#if 0
    const VecCoord& x_rest = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
#endif

    for(unsigned int i=0; i<nbEdges; i++ )
    {
        einfo=&edgeInf[i];

        // safety check
#if 0
        {
            EdgeInformation e2;
            const sofa::type::vector< unsigned int > shell = m_topology->getTrianglesAroundEdge(i);
            if (shell.size() != 2)
                e2.is_activated = false;
            else
            {
                e2.is_activated = true;
                e2.m1 = -1;
                e2.m2 = -1;
                for (int j=0; j<3; j++)
                    if (m_topology->getTriangle(shell[0]][j] != getEdge(i)[0] && m_topology->getTriangle(shell[0])[j] != getEdge(i)[1])
                        e2.m1 = m_topology->getTriangle(shell[0])[j];
                for (int j=0; j<3; j++)
                    if (m_topology->getTriangle(shell[1])[j] != getEdge(i)[0] && m_topology->getTriangle(shell[1])[j] != getEdge(i)[1])
                        e2.m2 = m_topology->getTriangle(shell[1])[j];
                if (e2.m1 >= 0 && e2.m2 >= 0)
                {
                    e2.restlength = (x_rest[e2.m2]-x_rest[e2.m1]).norm();
                }
            }

            if (e2.is_activated != einfo->is_activated) msg_error() << " EdgeInfo["<<i<<"].is_activated = "<<einfo->is_activated<<" while it should be "<<e2.is_activated<<"";
            else if (e2.is_activated)
            {
                if (!((e2.m1 == einfo->m1 && e2.m2 == einfo->m2) || (e2.m1 == einfo->m2 && e2.m2 == einfo->m1)))
                    msg_;() << "EdgeInfo["<<i<<"] points = "<<einfo->m1<<"-"<<einfo->m2<<" while it should be "<<e2.m1<<"-"<<e2.m2<<"";
                if (e2.restlength != einfo->restlength)
                    msg_error() << " EdgeInfo["<<i<<"] length = "<<einfo->restlength<<" while it should be "<<e2.restlength<<"";
            }
        }

#endif

        /*            msg_error()<<"TriangularBendingSprings<DataTypes>::addForce() between "<<springs[i].m1<<" and "<<springs[i].m2;*/

        if(einfo->is_activated)
        {
            //this->addSpringForce(m_potentialEnergy,f,x,v, i, einfo->spring);

            int a = einfo->m1;
            int b = einfo->m2;
            Coord u = x[b]-x[a];
            Real d = u.norm();
            if( d>1.0e-4 )
            {
                Real inverseLength = 1.0f/d;
                u *= inverseLength;
                Real elongation = (Real)(d - einfo->restlength);
                m_potentialEnergy += elongation * elongation * einfo->ks / 2;
                /*      msg_error()<<"TriangularBendingSprings<DataTypes>::addSpringForce, p = "<<p;

                        msg_error()<<"TriangularBendingSprings<DataTypes>::addSpringForce, new potential energy = "<<potentialEnergy;*/
                Deriv relativeVelocity = v[b]-v[a];
                Real elongationVelocity = dot(u,relativeVelocity);
                Real forceIntensity = (Real)(einfo->ks*elongation+einfo->kd*elongationVelocity);
                Deriv force = u*forceIntensity;
                f[a]+=force;
                f[b]-=force;

                updateMatrix=true;

                Mat& m = einfo->DfDx; //Mat& m = this->dfdx[i];
                Real tgt = forceIntensity * inverseLength;
                for( int j=0; j<N; ++j )
                {
                    for( int k=0; k<N; ++k )
                    {
                        m[j][k] = ((Real)einfo->ks-tgt) * u[j] * u[k];
                    }
                    m[j][j] += tgt;
                }
            }
            else // null length, no force and no stiffness
            {
                Mat& m = einfo->DfDx; //Mat& m = this->dfdx[i];
                for( int j=0; j<N; ++j )
                {
                    for( int k=0; k<N; ++k )
                    {
                        m[j][k] = 0;
                    }
                }
            }
        }
    }

    edgeInfo.endEdit();
    d_f.endEdit();
    //for (unsigned int i=0; i<springs.size(); i++)
    //{
    /*            msg_error()<<"TriangularBendingSprings<DataTypes>::addForce() between "<<springs[i].m1<<" and "<<springs[i].m2;*/
    //    this->addSpringForce(m_potentialEnergy,f,x,v, i, springs[i]);
    //}
}

template<class DataTypes>
void TriangularBendingSprings<DataTypes>::addDForce(const core::MechanicalParams* mparams, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
    VecDeriv& df = *d_df.beginEdit();
    const VecDeriv& dx = d_dx.getValue();
    Real kFactor = (Real)sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams, this->rayleighStiffness.getValue());

    size_t nbEdges=m_topology->getNbEdges();
    const EdgeInformation *einfo;
    const type::vector<EdgeInformation>& edgeInf = edgeInfo.getValue();

    df.resize(dx.size());
    //msg_error()<<"TriangularBendingSprings<DataTypes>::addDForce, dx1 = "<<dx1;
    //msg_error()<<"TriangularBendingSprings<DataTypes>::addDForce, df1 before = "<<f1;
    //const type::vector<Spring>& springs = this->springs.getValue();

    for(unsigned int i=0; i<nbEdges; i++ )
    {
        einfo=&edgeInf[i];

        /*            msg_error()<<"TriangularBendingSprings<DataTypes>::addForce() between "<<springs[i].m1<<" and "<<springs[i].m2;*/

        if(einfo->is_activated)
        {
            //this->addSpringDForce(df,dx, i, einfo->spring);

            const int a = einfo->m1;
            const int b = einfo->m2;
            const Coord d = dx[b]-dx[a];
            const Deriv dforce = einfo->DfDx*d; //this->dfdx[i]*d;
            df[a]+= dforce * kFactor;
            df[b]-= dforce * kFactor;
            //msg_error()<<"TriangularBendingSprings<DataTypes>::addSpringDForce, a="<<a<<", b="<<b<<", dforce ="<<dforce;

            //if(updateMatrix){
            //}
            updateMatrix=false;
        }
    }
    d_df.endEdit();
}


template<class DataTypes>
void TriangularBendingSprings<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!d_showSprings.getValue()) { return; }
    if (!vparams->displayFlags().getShowForceFields()) {return;}
    if (!this->mstate) {return;}

    vparams->drawTool()->saveLastState();

    if (vparams->displayFlags().getShowWireFrame()){
        vparams->drawTool()->setPolygonMode(0, true);
    }

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    std::vector<sofa::type::Vector3> vertices;
    std::vector<sofa::type::RGBAColor> colors;

    vparams->drawTool()->disableLighting();
    const type::vector<EdgeInformation>& edgeInfos = edgeInfo.getValue();
    for(auto& edgeInfo : edgeInfos)
    {
        if(edgeInfo.is_activated)
        {
            bool external=true;
            Real d = (x[edgeInfo.m2]-x[edgeInfo.m1]).norm();
            if (external)
            {
                if (d<edgeInfo.restlength*0.9999)
                {
                    colors.push_back(sofa::type::RGBAColor::red());
                }
                else
                {
                    colors.push_back(sofa::type::RGBAColor::green());
                }
            }
            else
            {
                if (d<edgeInfo.restlength*0.9999)
                {
                    colors.push_back(sofa::type::RGBAColor(1,0.5, 0,1));
                }
                else
                {
                    colors.push_back(sofa::type::RGBAColor(0,1,0.5,1));
                }
            }

            vertices.push_back( x[edgeInfo.m1] );
            vertices.push_back( x[edgeInfo.m2] );
        }
    }
    vparams->drawTool()->drawLines(vertices, 1, colors);

    if (vparams->displayFlags().getShowWireFrame()){
        vparams->drawTool()->setPolygonMode(0, false);
    }

    vparams->drawTool()->restoreLastState();
}


} // namespace sofa::component::forcefield