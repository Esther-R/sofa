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
#define SOFA_COMPONENT_ANIMATIONLOOP_MECHANICALMATRIXMAPPER_CPP
#include "MechanicalMatrixMapper.inl"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>
#include <fstream>

namespace sofa::component::interactionforcefield
{

using namespace sofa::defaulttype;

std::vector<int> MatrixProduct::sort_toSort(const std::vector<int> &vector){
    std::vector<int> index(vector.size(), 0);

    for (int i = 0 ; i != (int)index.size() ; i++) {
        index[i] = i;
    }
    sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
        return (vector[a] < vector[b]);
        }
    );
    return index ;
}

void MatrixProduct::computeTransposeB()
{
    const auto writeMatrix = [](Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, std::string path)
    {
        std::ofstream fichier(path.c_str(), std::ios::out | std::ios::trunc);
        if(fichier)  // si l'ouverture a réussi
        {
            fichier << mat << "\n";
            fichier.close();  // on referme le fichier
        }
    };
    writeMatrix(*m_B, "b.txt");

        int row;
        std::vector<int> toSort,columnsBTBeforePerm ;
        int nbrow =m_B->rows(); // number of row
        int nbcols=m_B->cols();
        int nbNonZero=m_B->nonZeros();// number of non zero values

        /////// rowIndexBT /////////

        // filling rowIndexBT with 0
        m_transposeB.rowIndex.clear();
        m_transposeB.rowIndex.resize(nbcols+1,0);

        // calculate how many non-zero per columns there are
        for(int i=0;i<nbNonZero;i++){
            m_transposeB.rowIndex[m_B->innerIndexPtr()[i]+1]++;
        }

        for(int i=1;i<nbcols+1;i++){
            m_transposeB.rowIndex[i]+=m_transposeB.rowIndex[i-1]; // we add up the number of values per row of BT to build rowIndexBT
        }

        /////////// perm and columns ///////////////
        for(int i=0;i<nbrow;i++){
            for(int j=m_B->outerIndexPtr()[i];j<m_B->outerIndexPtr()[i+1];j++){
                row=i;
                toSort.push_back(m_B->innerIndexPtr()[j]*nbrow+row);
                columnsBTBeforePerm.push_back(row); // calculate the row of the non-zero in B
            }
        }

        // calculate perm
        m_transposeB.perm=sort_toSort(toSort);

        // order columnsBT with perm
        for(int i=0;i<nbNonZero;i++){
            m_transposeB.columns.push_back(columnsBTBeforePerm[m_transposeB.perm[i]]);
        }
}


void MatrixProduct::computeIntersection()
{
    int nbNonZeroRow, nbNonZeroCol, offsetA, offsetBT, colA, colBT;
    bool doIntersect, isInfA;
    int index=0; // on l'utilise pour remplir rowC et columnsC
    int index_innerIndexC=0;
    std::vector<int> rowC;
    std::pair<int,int> couple;
    std::vector<std::pair<int,int>> list_couple;

    m_C.resize(m_A->rows(),m_B->cols());


//    writeMatrix(*m_B, "b.txt");

    m_C=(*m_A)*(*m_B);

//    writeMatrix(m_C, "c.txt");

    computeTransposeB();

    for(int r=0;r<m_A->rows();r++){
        nbNonZeroRow=m_A->outerIndexPtr()[r+1]-m_A->outerIndexPtr()[r];
        offsetA=m_A->outerIndexPtr()[r];
        for(int c=0;c<m_B->cols();c++){
            offsetBT=m_transposeB.rowIndex[c];
            nbNonZeroCol=m_transposeB.rowIndex[c+1]-m_transposeB.rowIndex[c];
            int itR=0,itC=0;
            while((itC<nbNonZeroCol) && (itR<nbNonZeroRow)){
                colA=m_A->innerIndexPtr()[offsetA+itR];
                colBT=m_transposeB.columns[offsetBT+itC];
                doIntersect=!(colA-colBT); // true if there is an intersection
                isInfA=colA<colBT;
                if(doIntersect){
                    // construction du couple ([offsetA+itR,perm[ossetBT+itC])
                    couple.first=offsetA+itR;
                    couple.second=m_transposeB.perm[offsetBT+itC];
                    // ajout du couple à la liste de couples
                    list_couple.push_back(couple);

//                    if(index==0){ // si c'est la première fois qu'on passe (pour éviter les doublons)
//                        rowC.push_back(r);
//                        m_C.innerIndexPtr()[index_innerIndexC++]=c;
//                    }
                    index++;
                }
                itC+=1*(doIntersect||!isInfA); // itC incremented if intersection or if colA>colBT
                itR+=1*(doIntersect||isInfA); // itR incremented if intersection or if colA<colBT
            }
            index=0;
            if(list_couple.size()){
                m_intersectionAB.intersection.push_back(list_couple);
            }
                list_couple.clear();
        }
    }

    /////// rowIndexC /////////
    // calculate how many non-zero per columns there are
//    for(int i=0;i<(int)m_intersectionAB.intersection.size();i++){
//        m_C.outerIndexPtr()[rowC[i]+1]++;
//    }

//    for(int i=1;i<m_A->rows();i++){
//        m_C.outerIndexPtr()[i]+=m_C.outerIndexPtr()[i-1]; // we add up the number of values per row of BT to build rowIndexBT
//    }
    m_hasComputedIntersection=true;


}


void MatrixProduct::computeProduct(){
    if (m_hasComputedIntersection == false)
    {
        computeIntersection();
    }
    else
    {
        double value;
        int index=0;
        for(const auto &list_couple:m_intersectionAB.intersection){
            value=0.0;
            for(const auto &couple:list_couple){
                int a=couple.first;
                int b=couple.second;
                value+=m_A->valuePtr()[a]*m_B->valuePtr()[b];
            }
            m_C.valuePtr()[index++]=(value);
        }
    }
}



////////////////////////////////////////////    FACTORY    //////////////////////////////////////////////
int MechanicalMatrixMapperClass = core::RegisterObject("This component allows to map the stiffness (and mass) matrix through a mapping.")

        .add< MechanicalMatrixMapper<Rigid3Types, Rigid3Types> >()
        .add< MechanicalMatrixMapper<Vec3Types, Rigid3Types> >()
        .add< MechanicalMatrixMapper<Vec3Types, Vec3Types> >()
        .add< MechanicalMatrixMapper<Vec1Types, Rigid3Types> >()
        .add< MechanicalMatrixMapper<Vec1Types, Vec3Types> >()
        .add< MechanicalMatrixMapper<Vec1Types, Vec1Types> >()
        .add< MechanicalMatrixMapper<Rigid3Types, Vec1Types> >(true)

        ;
////////////////////////////////////////////////////////////////////////////////////////////////////////

template class SOFA_SOFAGENERALANIMATIONLOOP_API MechanicalMatrixMapper<Rigid3Types, Rigid3Types>;
template class SOFA_SOFAGENERALANIMATIONLOOP_API MechanicalMatrixMapper<Vec3Types, Rigid3Types>;
template class SOFA_SOFAGENERALANIMATIONLOOP_API MechanicalMatrixMapper<Vec3Types, Vec3Types>;
template class SOFA_SOFAGENERALANIMATIONLOOP_API MechanicalMatrixMapper<Vec1Types, Rigid3Types>;
template class SOFA_SOFAGENERALANIMATIONLOOP_API MechanicalMatrixMapper<Vec1Types, Vec3Types>;
template class SOFA_SOFAGENERALANIMATIONLOOP_API MechanicalMatrixMapper<Vec1Types, Vec1Types>;
template class SOFA_SOFAGENERALANIMATIONLOOP_API MechanicalMatrixMapper<Rigid3Types, Vec1Types> ;



} // namespace sofa::component::interactionforcefield
