//
//  vForce.c
//
//  Calculation of forces on different lifting surfaces
//
//

#include "vForce.h"
#include "vLiftsurf.h"

void calcforces(struct liftsurf *plift, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int it, double dt, double UVW[], double rho, int cwing, int cflap, int caileron)
{
    /* Calculate aerodynamic forces acting on all liftsurfs */
    int i,neigh;
    double dGdt, Ux, Uy, Ui, Gx, Gy;
    
    *(plift->aeroforce)=0.0;
    *(plift->aeroforce+1)=0.0;
    *(plift->aeroforce+2)=0.0;
    *(plift->aeroforce+3)=0.0;
    
    for (i=0; i < plift->nface; i++){
        /* Calculate change of vorticity in time */
        if (it == 0){
            dGdt=*(plift->gamma+i+it*plift->nface)/dt;
        }else{
            dGdt=(*(plift->gamma+i+it*plift->nface)-*(plift->gamma+i+(it-1)*plift->nface))/dt;
        }
        
        Ux=(*(plift->uvw+i)+UVW[0])* *(plift->tangx+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangx+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangx+i+2*plift->nface);
        Uy=(*(plift->uvw+i)+UVW[0])* *(plift->tangy+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangy+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangy+i+2*plift->nface);
        Ui= *(plift->wind+i)+ *(plift->uvw+i+2*plift->nface);

        /* Calculate lift contribution of each panel */
        /* Separate calculations for each half-wing */
        
        neigh=*(plift->neighbours+i+2*plift->nface);
        Gx=0;
        Gy=0;
        if (neigh == -1){
            Gx=*(plift->gamma+i+it*plift->nface);
        }else if (neigh < cflap && neigh > cwing-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pwing->gamma+neigh-cwing+it*pwing->nface);
        }else if (neigh < caileron && neigh > cflap-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pflap->gamma+neigh-cflap+it*pflap->nface);
        }else if (neigh > caileron-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(paileron->gamma+neigh-caileron+it*paileron->nface);
        }      
        if (i < plift->nface/2){
            neigh=*(plift->neighbours+i);
            if (neigh == -1){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < cflap && neigh > cwing-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pwing->gamma+neigh-cwing+it*pwing->nface);
            }else if (neigh < caileron && neigh > cflap-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pflap->gamma+neigh-cflap+it*pflap->nface);
            }else if (neigh > caileron-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(paileron->gamma+neigh-caileron+it*paileron->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)+Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }else{
            neigh=*(plift->neighbours+i+plift->nface);
            if (neigh == -1){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < cflap && neigh > cwing-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pwing->gamma+neigh-cwing+it*pwing->nface);
            }else if (neigh < caileron && neigh > cflap-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pflap->gamma+neigh-cflap+it*pflap->nface);
            }else if (neigh > caileron-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(paileron->gamma+neigh-caileron+it*paileron->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)-Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }
        /* Calculate induced drag contribution of each panel */
        *(plift->Deltad+i)=rho*(-Ui*Gx* *(plift->dxy+i+plift->nface)+dGdt* *(plift->nsurf+i));
        
        *(plift->aeroforce) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i);
        *(plift->aeroforce+1) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+1*plift->nface);
        *(plift->aeroforce+2) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+2*plift->nface);
        *(plift->aeroforce+3) += *(plift->Deltad+i);
    }
}

void calcforceshtail(struct liftsurf *plift, struct liftsurf *phtail, struct liftsurf *pelevator, int it, double dt, double UVW[], double rho, int chtail, int celevator)
{
    /* Calculate aerodynamic forces acting on all liftsurfs */
    int i,neigh;
    double dGdt, Ux, Uy, Ui, Gx, Gy;

    *(plift->aeroforce)=0.0;
    *(plift->aeroforce+1)=0.0;
    *(plift->aeroforce+2)=0.0;
    *(plift->aeroforce+3)=0.0;

    for (i=0; i < plift->nface; i++){
        /* Calculate change of vorticity in time */
        if (it == 0){
            dGdt=*(plift->gamma+i+it*plift->nface)/dt;
        }else{
            dGdt=(*(plift->gamma+i+it*plift->nface)-*(plift->gamma+i+(it-1)*plift->nface))/dt;
        }

        Ux=(*(plift->uvw+i)+UVW[0])* *(plift->tangx+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangx+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangx+i+2*plift->nface);
        Uy=(*(plift->uvw+i)+UVW[0])* *(plift->tangy+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangy+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangy+i+2*plift->nface);
        Ui= *(plift->wind+i)+ *(plift->uvw+i+2*plift->nface);

        /* Calculate lift contribution of each panel */
        /* Separate calculations for each half-wing */

        neigh=*(plift->neighbours+i+2*plift->nface);
        Gx=0;
        Gy=0;
        if (neigh < chtail){
            Gx=*(plift->gamma+i+it*plift->nface);
        }else if (neigh < celevator && neigh > chtail-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(phtail->gamma+neigh-chtail+it*phtail->nface);
        }else if (neigh > celevator-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pelevator->gamma+neigh-celevator+it*pelevator->nface);
        }
        if (i < plift->nface/2){
            neigh=*(plift->neighbours+i);
            if (neigh < chtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < celevator && neigh > chtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(phtail->gamma+neigh-chtail+it*phtail->nface);
            }else if (neigh > celevator-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pelevator->gamma+neigh-celevator+it*pelevator->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)+Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }else{
            neigh=*(plift->neighbours+i+plift->nface);
            if (neigh < chtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < celevator && neigh > chtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(phtail->gamma+neigh-chtail+it*phtail->nface);
            }else if (neigh > celevator-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pelevator->gamma+neigh-celevator+it*pelevator->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)-Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }
        /* Calculate induced drag contribution of each panel */
        *(plift->Deltad+i)=rho*(-Ui*Gx* *(plift->dxy+i+plift->nface)+dGdt* *(plift->nsurf+i));

        *(plift->aeroforce) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i);
        *(plift->aeroforce+1) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+1*plift->nface);
        *(plift->aeroforce+2) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+2*plift->nface);
        *(plift->aeroforce+3) += *(plift->Deltad+i);
    }
}

void calcforcesvtail(struct liftsurf *plift, struct liftsurf *pvtail, struct liftsurf *prudder, int it, double dt, double UVW[], double rho, int cvtail, int crudder)
{
    /* Calculate aerodynamic forces acting on all liftsurfs */
    int i,neigh;
    double dGdt, Ux, Uy, Ui, Gx, Gy;

    *(plift->aeroforce)=0.0;
    *(plift->aeroforce+1)=0.0;
    *(plift->aeroforce+2)=0.0;
    *(plift->aeroforce+3)=0.0;

    for (i=0; i < plift->nface; i++){
        /* Calculate change of vorticity in time */
        if (it == 0){
            dGdt=*(plift->gamma+i+it*plift->nface)/dt;
        }else{
            dGdt=(*(plift->gamma+i+it*plift->nface)-*(plift->gamma+i+(it-1)*plift->nface))/dt;
        }

        Ux=(*(plift->uvw+i)+UVW[0])* *(plift->tangx+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangx+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangx+i+2*plift->nface);
        Uy=(*(plift->uvw+i)+UVW[0])* *(plift->tangy+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangy+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangy+i+2*plift->nface);
        Ui= *(plift->wind+i)+ *(plift->uvw+i+2*plift->nface);

        /* Calculate lift contribution of each panel */
        /* Separate calculations for each half-wing */

        neigh=*(plift->neighbours+i+2*plift->nface);
        Gx=0;
        Gy=0;
        if (neigh < cvtail){
            Gx=*(plift->gamma+i+it*plift->nface);
        }else if (neigh < crudder && neigh > cvtail-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pvtail->gamma+neigh-cvtail+it*pvtail->nface);
        }else if (neigh > crudder-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(prudder->gamma+neigh-crudder+it*prudder->nface);
        }
        if (i < plift->nface/2){
            neigh=*(plift->neighbours+i);
            if (neigh < cvtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < crudder && neigh > cvtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pvtail->gamma+neigh-cvtail+it*pvtail->nface);
            }else if (neigh > crudder-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(prudder->gamma+neigh-crudder+it*prudder->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)+Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }else{
            neigh=*(plift->neighbours+i+plift->nface);
            if (neigh < cvtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < crudder && neigh > cvtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pvtail->gamma+neigh-cvtail+it*pvtail->nface);
            }else if (neigh > crudder-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(prudder->gamma+neigh-crudder+it*prudder->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)-Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }
        /* Calculate induced drag contribution of each panel */
        *(plift->Deltad+i)=rho*(-Ui*Gx* *(plift->dxy+i+plift->nface)+dGdt* *(plift->nsurf+i));

        *(plift->aeroforce) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i);
        *(plift->aeroforce+1) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+1*plift->nface);
        *(plift->aeroforce+2) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+2*plift->nface);
        *(plift->aeroforce+3) += *(plift->Deltad+i);
    }
}
