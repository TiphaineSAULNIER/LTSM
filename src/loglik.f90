
module modltsm

  implicit none
  integer,save ::ny,ns,nv,neatot &
       ,maxmes,nobs,ntrtot &
       ,nySPL,ntotvalSPL,npmtot &
       ,nMC,methInteg,nmescur &
       ,nvarxevt,nvarbyevt,nbevt,logspecif,idtrunc,nrisqtot,nxevt,nasso,nasso_s_fp,nasso_s
  integer,dimension(:),allocatable,save::nassoparevt
  integer,save::idshiftasso     
  integer,save::nef_s,nvc_s,nvc_err 
  integer,dimension(:),allocatable,save::nef,nea,nvc
  integer,dimension(:),allocatable,save::nea_fp
  integer,dimension(:),allocatable,save::nef_fp,nef_var,nef_interac
  integer,dimension(:),allocatable,save::typrisq,nz,nprisq,nevtparx
  double precision,dimension(:),allocatable,save::Y,uniqueY,minY,maxY
  double precision,dimension(:,:),allocatable,save ::X,zi
  double precision,dimension(:),allocatable,save::Tsurv0,Tsurv
  integer,dimension(:),allocatable,save::Devt
  integer,dimension(:),allocatable,save::idsurv,idsurvsfp
  integer,save :: idsurvs
  integer,dimension(:,:),allocatable,save::ideafp,idgfp,idg
  integer,dimension(:),allocatable,save::idea,idginterac
  integer,dimension(:),allocatable,save::idgs,indiceY
  integer,dimension(:),allocatable,save::idmodel,idlink,ntr
  integer,dimension(:,:),allocatable,save::nmes
  integer,dimension(:),allocatable,save::nvalSPL
  double precision,dimension(:),allocatable,save :: seqMC
  double precision,dimension(:,:),allocatable,save::zitr
  double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2
  double precision,save::Tmm,Tmm1,&
       Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02,Tmm03,Tim0, &
       Tim01,Tim02,Tim03,Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1,&
       Timt2,Timt3
  double precision,dimension(:),allocatable,save::Tmm_GK,Tmm1_GK,&
       Tmm2_GK,Tmm3_GK,Tmm0_GK,Tmm01_GK,Tmm02_GK,Tmm03_GK     
  integer,dimension(:),allocatable,save::fix
  double precision,dimension(:),allocatable,save::bfix,refix
  double precision,dimension(:),allocatable,save::bp
  double precision,save::nu,rate
  double precision,dimension(:),allocatable,save::asymptL,asymptU
  integer, dimension(:),allocatable,save::npm_t,npm_x,npm_tx,npm_ea_t
  double precision,dimension(:,:),allocatable,save::b0_t,b0_x,b0_tx
  double precision,dimension(:,:),allocatable,save::X0_x
  real,dimension(8),save :: FPpower
  integer,save::nGK
  integer,dimension(:,:),allocatable,save::idst
  double precision,dimension(:),allocatable,save::GKw,GKp
  double precision,dimension(:),allocatable,save::Tsurv_GK,Tsurv0_GK
  double precision,save::centr,hlgth,centr0,hlgth0
  
end module modltsm




subroutine loglik(Y0,X0,Tentr0,Tevt0,Devt0 &
     ,ideafp0,idgfp0,idea0,idg0,idginterac0,idgs0,idsurv0,idsurvs0,idsurvsfp0,sharedtype0,nGK0 &
     ,typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0 &
     ,ny0,ns0,nv0 &
     ,nobs0,nmes0,npm0,b0,nfix0,bfix0 &
     ,idmodel0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0,nvalSPLORD0,fix0 &
     ,methInteg0,nMC0,dimMC0,seqMC0 &
     ,loglik_res)

  use modltsm

  IMPLICIT NONE
  
  !Declaration des variables en entree
  integer,intent(in)::nv0,ny0,nMC0,methInteg0,dimMC0,nfix0
  integer, intent(in)::ns0,nobs0,npm0
  integer,intent(in)::idtrunc0,logspecif0,nbevt0
  integer, dimension(nbevt0,ny0), intent(in)::sharedtype0
  integer, intent(in)::nGK0
  double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
  integer, dimension(ns0),intent(in)::Devt0
  integer, dimension(nv0),intent(in)::idsurv0
  integer, intent(in)::idsurvs0
  integer, dimension(8),intent(in)::idsurvsfp0
  integer,dimension(nbevt0),intent(in)::typrisq0,nz0
  double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0 
  integer, dimension(ny0),intent(in)::idmodel0,idlink0,nbzitr0,nvalSPLORD0
  double precision,dimension(maxval(nbzitr0),ny0),intent(in)::zitr0
  integer,dimension(nobs0),intent(in)::indiceY0
  double precision,dimension(sum(nvalSPLORD0(:))),intent(in)::uniqueY0
  integer,dimension(ny0,8),intent(in)::ideafp0,idgfp0
  integer, dimension(ny0),intent(in)::idea0,idginterac0
  integer, dimension(ny0,nv0),intent(in)::idg0
  integer, dimension(nv0),intent(in)::idgs0
  integer,dimension(ns0,ny0)::nmes0   
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  integer,dimension(npm0+nfix0),intent(in)::fix0
  double precision,dimension(dimMC0*nMC0),intent(in)::seqMC0
  double precision, dimension(npm0), intent(in) :: b0
  double precision, dimension(nfix0), intent(in) :: bfix0
  
  !Declaration des variables en sortie
  double precision,intent(out)::loglik_res
  
  !Variables locales
  integer::jtemp,i,j,ier,k,ktemp,yk,k1,mi,nbfix,l
  integer::npmtot0,ke
  double precision::eps
  double precision,external::vrais
  
  !open(unit=3,file='output.txt') !pour affichages
  
  !print*,"ELEMENTS EN ENTREE FORTRAN"
  !print*,"Y0",Y0
  !print*,"X0",X0
  !print*,"Tentr0",Tentr0
  !print*,"Tevt0",Tevt0
  !print*,"Devt0",Devt0
  !print*,"ideafp0",ideafp0
  !print*,"idgfp0",idgfp0
  !print*,"idea0",idea0
  !print*,"idg0",idg0
  !print*,"idginterac0",idginterac0
  !print*,"idgs0",idgs0
  !print*,"idsurv0",idsurv0
  !print*,"idsurvs0",idsurvs0
  !print*,"idsurvsfp0",idsurvsfp0
  !print*,"sharedtype0",sharedtype0
  !print*,"nGK0",nGK0
  !print*,"typrisq0",typrisq0
  !print*,"nz0",nz0
  !print*,"zi0",zi0
  !print*,"nbevt0",nbevt0
  !print*,"idtrunc0",idtrunc0
  !print*,"logspecif0",logspecif0
  !print*,"ny0",ny0
  !print*,"ns0",ns0
  !print*,"nv0",nv0
  !print*,"nobs0",nobs0
  !print*,"nmes0",nmes0
  !print*,"npm0",npm0
  !print*,"b0",b0
  !print*,"nfix0",nfix0
  !print*,"bfix0",bfix0
  !print*,"idmodel0",idmodel0
  !print*,"idlink0",idlink0
  !print*,"nbzitr0",nbzitr0
  !print*,"zitr0",zitr0
  !print*,"uniqueY0",uniqueY0
  !print*,"indiceY0",indiceY0
  !print*,"nvalSPLORD0",nvalSPLORD0
  !print*,"fix0",fix0
  !print*,"methInteg0",methInteg0
  !print*,"nMC0",nMC0
  !print*,"dimMC0",dimMC0
  !print*,"seqMC0",seqMC0
  !print*,"loglik_res",loglik_res
  !print*,"FIN ELEMENTS EN ENTREE FORTRAN"
  
  maxmes=0
  do i=1,ns0
     mi=sum(nmes0(i,:))
     if (mi.gt.maxmes) then
        maxmes=mi
     end if
  end do
  
  allocate(minY(ny0),maxY(ny0),idmodel(ny0),idlink(ny0),ntr(ny0))
  
  FPpower = [ -2.0 , -1.0 , -0.5 , 0.0 , 0.5 , 1.0 , 2.0 , 3.0 ]
  
  nySPL=0
  do k=1,ny0
     idmodel(k)=idmodel0(k)
     idlink(k)=idlink0(k)
     minY(k)=zitr0(1,k)
     maxY(k)=zitr0(nbzitr0(k),k)
     if (idlink(k).eq.2) then
        nySPL=nySPL+1
     end if
  end do
  
  if(nySPL.gt.0) then 
     allocate(nvalSPL(nySPL))
     nvalSPL=0
  else
     allocate(nvalSPL(1))
     nvalSPL(1) = 0
  end if

  k1=0
  do k=1,ny0
     if(idlink(k).eq.2) then
        k1=k1+1
        nvalSPL(k1)=nvalSPLORD0(k)
     end if
  end do
  ntotvalSPL=sum(nvalSPL(:))

  methInteg = methInteg0
  nMC = nMC0
  
  if(all(idlink.ne.2)) then
     allocate(zitr(1,1))
     allocate(mm(1),mm1(1),mm2(1),im(1),im1(1),im2(1))
     mm(1)=0.d0
     mm1(1)=0.d0
     mm2(1)=0.d0
     im(1)=0.d0
     im1(1)=0.d0
     im2(1)=0.d0
  else
     allocate(zitr(-1:(maxval(nbzitr0)+2),nySPL))
     allocate(mm(ntotvalSPL),mm1(ntotvalSPL),mm2(ntotvalSPL),im(ntotvalSPL),im1(ntotvalSPL),im2(ntotvalSPL))
  end if


  zitr=0.d0  
  k1=0
  do k=1,ny0
     if (idlink(k).eq.0) ntr(k)=0   
     if (idlink(k).eq.1) ntr(k)=2
     if (idlink(k).eq.2) then
        k1=k1+1
        ntr(k)=nbzitr0(k)+2

        zitr(1:nbzitr0(k),k1)=zitr0(1:nbzitr0(k),k)
        zitr(-1,k1)=zitr(1,k1)
        zitr(0,k1)=zitr(1,k1)
        zitr(ntr(k)-1,k1)=zitr(ntr(k)-2,k1)
        zitr(ntr(k),k1)=zitr(ntr(k)-1,k1)
     end if
  end do

  ntrtot = sum(ntr(:))
  
  allocate(Y(nobs0),X(nobs0,nv0),uniqueY(ntotvalSPL) &
       ,ideafp(ny0,8),idgfp(ny0,8),idea(ny0),idg(ny0,nv0),idginterac(ny0),idgs(nv0),nmes(ns0,ny0) &
       ,indiceY(nobs0))
  allocate(asymptL(ny0),asymptU(ny0))
  allocate(npm_t(ny0),npm_x(ny0),npm_tx(ny0),npm_ea_t(ny0))
  allocate(b0_t(ny0,2*8),b0_x(ny0,nv0),b0_tx(ny0,nv0*(2*8)))
  allocate(X0_x(maxmes,nv0))
       
  allocate(nprisq(nbevt0),nevtparx(nv0))     

  if(nbevt0.gt.0) then
  
    allocate(Tsurv0(ns0),Tsurv(ns0),Devt(ns0))
    allocate(typrisq(nbevt0),nz(nbevt0))
    allocate(idsurv(nv0))
    allocate(idsurvsfp(8))
    
    ! zi : contient noeuds pour hazard (ou min et max si Weibull)
    if(any(typrisq0.eq.3)) then !spl
      allocate(zi(-2:maxval(nz0)+3,nbevt0))
    else
      allocate(zi(maxval(nz0),nbevt0))
    end if
  
    ! enregistrement pour les modules
    typrisq=typrisq0
    idtrunc=idtrunc0
    Tsurv0=Tentr0   
    Tsurv=Tevt0    
    devt=devt0    
    logspecif=logspecif0 
    
    !
    nGK=nGK0
    allocate(idst(nbevt0,ny0))
    idst=0
    do ke=1,nbevt0
      do yk=1,ny0
        idst(ke,yk)=sharedtype0(ke,yk)
      end do
    end do
    !print*,"idst",idst
    
    if(any(idst.eq.1)) then   ! any CV in surv
    
      allocate(Tsurv_GK(nGK),Tsurv0_GK(nGK))
      allocate(GKw(nGK),GKp(nGK))
      GKw=0.d0
      GKp=0.d0
      if(nGK.eq.7) then
        ! weights
        GKw(1)=0.1294849661688696932706114326790820
        GKw(2)=0.2797053914892766679014677714237796
        GKw(3)=0.3818300505051189449503697754889751
        GKw(4)=0.4179591836734693877551020408163265
        GKw(5)=0.3818300505051189449503697754889751
        GKw(6)=0.2797053914892766679014677714237796
        GKw(7)=0.1294849661688696932706114326790820
        ! points
        GKp(1)=-0.9491079123427585245261896840478513
        GKp(2)=-0.7415311855993944398638647732807884
        GKp(3)=-0.4058451513773971669066064120769615
        GKp(4)=0.0000000000000000000000000000000000
        GKp(5)=0.4058451513773971669066064120769615
        GKp(6)=0.7415311855993944398638647732807884
        GKp(7)=0.9491079123427585245261896840478513
      end if
      if(nGK.eq.15) then
        ! weights
        GKw(1)=0.022935322010529224963732008058970
        GKw(2)=0.063092092629978553290700663189204
        GKw(3)=0.104790010322250183839876322541518
        GKw(4)=0.140653259715525918745189590510238
        GKw(5)=0.169004726639267902826583426598550
        GKw(6)=0.190350578064785409913256402421014
        GKw(7)=0.204432940075298892414161999234649
        GKw(8)=0.209482141084727828012999174891714
        GKw(9)=0.204432940075298892414161999234649
        GKw(10)=0.190350578064785409913256402421014
        GKw(11)=0.169004726639267902826583426598550
        GKw(12)=0.140653259715525918745189590510238
        GKw(13)=0.104790010322250183839876322541518
        GKw(14)=0.063092092629978553290700663189204
        GKw(15)=0.022935322010529224963732008058970
        ! points
        GKp(1)=-0.9914553711208126392068546975263285
        GKp(2)=-0.9491079123427585245261896840478513
        GKp(3)=-0.8648644233597690727897127886409262
        GKp(4)=-0.7415311855993944398638647732807884
        GKp(5)=-0.5860872354676911302941448382587296
        GKp(6)=-0.4058451513773971669066064120769615
        GKp(7)=-0.2077849550078984676006894037732449	
        GKp(8)=0.0000000000000000000000000000000000
        GKp(9)=0.2077849550078984676006894037732449
        GKp(10)=0.4058451513773971669066064120769615
        GKp(11)=0.5860872354676911302941448382587296
        GKp(12)=0.7415311855993944398638647732807884
        GKp(13)=0.8648644233597690727897127886409262
        GKp(14)=0.9491079123427585245261896840478513
        GKp(15)=0.9914553711208126392068546975263285
      end if
      !print*,"GKw",GKw
      !print*,"GKp",GKp
      
      if(any(typrisq.eq.3)) then
        allocate(Tmm_GK(nGK),Tmm1_GK(nGK),Tmm2_GK(nGK),Tmm3_GK(nGK),&
                 Tmm0_GK(nGK),Tmm01_GK(nGK),Tmm02_GK(nGK),Tmm03_GK(nGK))
      end if
      
    end if
  
  end if
  
  
  eps=1.d-20
  
  nbevt=nbevt0
  ny=ny0
  ns=ns0
  nv=nv0
  nobs=nobs0
  
  
  if (ntotvalSPL.gt.0) uniqueY(1:ntotvalSPL)=uniqueY0(1:ntotvalSPL)
  
  nmes=0
  Y=0.d0
  X=0.d0
  idea=0
  idginterac=0
  idg=0
  idgs=0
  if(nbevt.gt.0) then
     idsurv=0
     idsurvs=idsurvs0
  end if
  ktemp=0
  
  do yk=1,ny
      idea(yk)=idea0(yk)
      idginterac(yk)=idginterac0(yk)
  end do
  
  do k=1,nv
  
     if(nbevt.gt.0) then
       idsurv(k)=idsurv0(k)
     end if
     
     idgs(k)=idgs0(k)
     
     do yk=1,ny
      idg(yk,k)=idg0(yk,k)
     end do
     
     jtemp=0
     do i=1,ns
        do yk=1,ny            
           if (k.eq.1) then
              nmes(i,yk)=nmes0(i,yk)   !dim(nmes)=ns*ny    
              do j=1,nmes(i,yk)
                 jtemp=jtemp+1
                 Y(jtemp)=Y0(jtemp)
                 indiceY(jtemp)=indiceY0(jtemp)
                 ktemp=ktemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           else
              do j=1,nmes(i,yk)
                 ktemp=ktemp+1
                 jtemp=jtemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           end if
        end do
     end do
     
  end do
  
  idgfp=0
  ideafp=0
  do yk=1,ny
    do k=1,8
      idgfp(yk,k)=idgfp0(yk,k)
      ideafp(yk,k)=ideafp0(yk,k)
    end do
  end do
  !print*,'idgfp',idgfp
  !print*,'ideafp',ideafp
  
  if(nbevt0.gt.0) then
     idsurvsfp = 0
     do k=1,8
       idsurvsfp(k)=idsurvsfp0(k)
     end do
  end if
  
  
  ! prm fixes
  npmtot0 = npm0+nfix0
  allocate(fix(npmtot0))
  fix=0
  fix(1:npmtot0)=fix0(1:npmtot0)
  nbfix=sum(fix(:))
  if(nbfix.eq.0) then
     allocate(bfix(1))
  else
     allocate(bfix(nbfix))
  end if
  bfix(1:nbfix)=bfix0(1:nbfix)
  
  
  ! creation des parametres
  
  nprisq=0
  nrisqtot=0
  
  if (nbevt.gt.0) then
  
    do ke=1,nbevt

     nz(ke)=nz0(ke) ! nb de noeuds pour hazard (ou 2 pr Weibull)

     if (typrisq(ke).eq.1) then
        nprisq(ke)=nz(ke)-1
     end if
     if (typrisq(ke).eq.2) then
        nprisq(ke)=2
     end if
     if (typrisq(ke).eq.3) then
        nprisq(ke)=nz(ke)+2
     end if

     nrisqtot = nrisqtot+nprisq(ke)  ! nb total de prm pour hazards
     zi(1:nz(ke),ke)=zi0(1:nz(ke),ke)
     
    end do
    
  end if
  
  ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
  nxevt=0
  nevtparx=0
  if(nbevt.gt.0) then
    do j=1,nv
      if(idsurv(j).eq.1) then
        nevtparx(j) = nbevt
        nxevt = nxevt + 1
      end if
    end do
  end if

  nvarbyevt = nxevt
  nvarxevt = sum(nevtparx(:)) !+ nvdepsurv
  
  nef_s=0
  do k=1,nv
    if (idgs(k).eq.1) then
      nef_s = nef_s + 1
    end if
  end do
  nvc_s = 1
  
  allocate(nea_fp(ny),nea(ny))
  allocate(nef_fp(ny),nef_var(ny),nef_interac(ny),nef(ny))
  nea_fp=0
  nea=0
  nef_fp=0
  nef_var=0
  nef_interac=0
  nef=0
  do yk=1,ny
  
    do k=1,8
      nef_fp(yk)=nef_fp(yk)+idgfp(yk,k)
      nea_fp(yk)=nea_fp(yk)+ideafp(yk,k)
    end do
    do k=1,nv
      if(idg(yk,k).ne.0) nef_var(yk)=nef_var(yk)+1
    end do
    if(idginterac(yk).eq.1) nef_interac(yk)=nef_fp(yk)*nef_var(yk)
    
    nef(yk)=nef_fp(yk)+nef_var(yk)+nef_interac(yk)
    nea(yk)=nea_fp(yk)+idea(yk)
    
    if(idmodel(yk).eq.0) then !logistic
      nef(yk)=nef(yk)+1 ! pr nu
      nea(yk)=nea(yk)+2 ! pr nu et rate
    end if
    
  end do
 !print*,"nef_fp",nef_fp
 !print*,"nef_var",nef_var
 !print*,"nef_interac",nef_interac
 !print*,"nef",nef
  !print*,"nea_fp",nea_fp
  !print*,"nea",nea
  
  allocate(nvc(ny))
  nvc=0
  do yk=1,ny
    if(idmodel(yk).eq.0) then !logistic
      nvc(yk) = nea(yk)
    end if
    if(idmodel(yk).eq.1) then !linear
      nvc(yk) = (nea(yk)+1)*nea(yk)/2 
      if(idlink(yk).ne.0) nvc(yk) = nvc(yk)-1
    end if
  end do
  
  neatot = 0
  neatot = sum(nea) + nvc_s
  !print*,"neatot",neatot
  
  nvc_err = ny 
  
  allocate(nassoparevt(nbevt))
  
  nasso_s = 0
  nasso_s_fp = 0
  nasso = 0
  nassoparevt = 0
  if(nbevt.gt.0) then
  
    nasso_s = idsurvs
    nasso_s_fp = sum(idsurvsfp)
    
    do ke=1,nbevt
      
      !markers' characteristics
      do yk=1,ny
        if(idst(ke,yk).eq.0) then !RE
          nassoparevt(ke) = nassoparevt(ke) + nea(yk)
        end if
        if(idst(ke,yk).eq.1) then !CV
          nassoparevt(ke) = nassoparevt(ke) + 1
        end if
      end do
      
      !shift
      nassoparevt(ke) = nassoparevt(ke) + nasso_s + nasso_s_fp
      
    end do
  end if
  nasso=sum(nassoparevt)
  
   !write(3,*)"nassoparevt =",nassoparevt
   !write(3,*)"nasso =",nasso

  
  
  npmtot = nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sum(nef)+sum(nvc)+ntrtot+nvc_err
  !print*,"npm0",npm0
  !print*,"npmtot",npmtot
  
  ! points qmc
  if(methInteg.ne.3) then 
     allocate(seqMC(1))
  else
     allocate(seqMC(dimMC0*nMC))
     seqMC = seqMC0(1:dimMC0*nMC) 
  end if
  
  
  ! base de splines transfos
  if (any(idlink.eq.2)) then 
     call design_splines(ier)
     if (ier.eq.-1) then
        loglik_res=-1.d9
        go to 1589
     end if
  end if
  
  
  !!!!!!!!!!
  !print*,"deb vrais"

  ! calcul de la vraisemblance
  loglik_res = vrais(b0,npm0)

  !print*,"fin vrais"
  
  !print*,"loglik_res",loglik_res

  1589 continue
  !!!!!!!!!!
  
  
  if(nbevt.gt.0) then
     if(any(idst.eq.1)) then
       deallocate(GKw,GKp)
       deallocate(Tsurv_GK,Tsurv0_GK)
       if(any(typrisq.eq.3)) then
         deallocate(Tmm_GK,Tmm1_GK,Tmm2_GK,Tmm3_GK,&
                    Tmm0_GK,Tmm01_GK,Tmm02_GK,Tmm03_GK)
       end if
     end if
     deallocate(Tsurv0,Tsurv,zi,devt,typrisq,nz,idsurv)
     deallocate(idsurvsfp)
     deallocate(idst)
  end if
  
  deallocate(nprisq,nevtparx)
  deallocate(nassoparevt)
  
  deallocate(nef_fp,nef_var,nef_interac,nef)
  deallocate(nea_fp,nea)
  deallocate(nvc)

  deallocate(Y,X,ideafp,idgfp,idgs,idea,idg,idginterac,nmes,uniqueY,indiceY,ntr)
  deallocate(asymptL,asymptU)
  deallocate(X0_x)
  deallocate(b0_t,b0_x,b0_tx)
  deallocate(npm_t,npm_x,npm_tx,npm_ea_t)

  deallocate(zitr,mm,mm1,mm2,im,im1,im2,minY,maxY,idmodel,idlink,nvalSPL)

  deallocate(fix,bfix,seqMC)
  
  !close(unit=3) !fin affichages

  return
  
end subroutine loglik



!-------------------------------------!
!     LOG-LIKELIHOOD all subjects     !
!-------------------------------------!


double precision function vrais(b,m)

  use modltsm,only:ns,nmes,nmescur

  implicit none

  integer::m,i
  double precision::vrais_i,temp
  double precision,dimension(m)::b
  
  
  nmescur=0
  vrais=0.d0
  do i=1,ns
  
     !print*,"##################"
     !print*,"## new subject ",i 
     
     temp = vrais_i(b,m,i) 
      
     !print*,"vrais_i ",temp
     
     vrais = vrais + temp
     
     !print*,"vrais ",vrais
  
     
     if (temp.eq.-1.d9 .or. temp/temp.ne.1) then 
        !if (temp/temp.ne.1) write(*,*)"i=",i,"vrais= ",temp
        !if (temp.eq.-1.d9) then 
        vrais = -1.d9
        !print*,"dans vrais i=",i," vrais=",vrais," m=",m," b=",b
        !if(verbose==1) write(*,*)"i=",i,"vrais= ",temp
        goto 541
     end if
     
     nmescur = nmescur + sum(nmes(i,:))
  
  end do
  
  541 continue
  return

end function vrais



!-------------------------------------!
!      individual LOG-LIKELIHOOD      !
!-------------------------------------!

double precision function vrais_i(b,npm,i) 

  use modltsm
!  use optim

  IMPLICIT NONE
  integer ::i,j,k,l,m,jj,npm,ll,ii,numSPL,p
  integer ::ier,kk,j1,j2,sumMesYk,yk,ke,sumnrisq
  integer :: sumRe,sumasso
  integer::sumPrmYk
  integer::nxevtcurr
  
  
  double precision,dimension(maxmes,nv) :: X0_x_yk
  double precision,dimension(maxmes,2*8) :: X0_t,Z_t
  double precision,dimension(2*8) :: ui_Z_t
  double precision,dimension(maxmes,nv*(2*8)) :: X0_tx
  
  double precision,dimension(nef_s) :: b_s,X_s
  double precision :: mu_s
  double precision,dimension(maxmes) :: Xtime
  
  double precision,dimension(ny) :: mu_v
  double precision,dimension(nv) :: b_r,X_r
  double precision,dimension(ny) :: mu_r
  
  double precision,dimension(maxmes,maxmes) :: var_Y, inv_var_Y
  double precision,dimension(maxmes*(maxmes+1)/2) :: var_Y_vect
  double precision :: det_var_Y
  
  double precision :: shift
  double precision,dimension(maxmes) :: time_shifted
  double precision,dimension(maxmes) :: denom,mu
  
  double precision,dimension(neatot,neatot) ::Ut
  double precision,dimension(npm) :: b
  double precision,dimension(npmtot)::b1
  double precision,dimension(nvarxevt)::Xevt,bevt
  double precision,dimension(maxval(nprisq))::brisq
  
  double precision :: eps,som,eta0
  double precision ::Y4,jacobien
  double precision,dimension(maxmes) :: Y1,Y2,Y3
  double precision,dimension(-1:maxval(ntr)-3)::splaa
  double precision,dimension(neatot)::ui,usim
  double precision,dimension(nasso)::ui_asso
  double precision,dimension(nbevt)::risq,surv,surv0
  double precision::SX,x22,div,vrais_Y,vrais_surv,varexpsurv
  double precision::surv0_glob,surv_glob,fevt,easurv  
  double precision::som_T0,som_Ti
  
  double precision :: shiftsurv
  
  
  ! definir le nombre total de mesures pour le sujet i : nmestot (valable que pour cette fonction)

  ! if (verbose==1) write(*,*)'i',i 
  b1=0.d0
  eps=1.d-20
  l=0
  m=0
  do k=1,npmtot
     if(fix(k).eq.0) then
        l=l+1
        b1(k)=b(l)
     end if
     if(fix(k).eq.1) then
        m=m+1
        b1(k)=bfix(m)
     end if
  end do
  
  !print*, "b1", b1

  !----------- rappel des parametres utilises ---------

  !        write(*,*)'i',i,nmescur,nmescur + nmes(i)

  l = 0
  ll = 0
  sumPrmYk = 0
  Ut = 0.d0
  Ut(1,1) = abs(b1(nrisqtot+nvarxevt+nasso+nef_s+1)) !shift std
  do yk=1,ny
    if(idmodel(yk).eq.0) then !logistic
      do j=1,nea(yk)
        Ut(1+l+j,1+l+j)=abs(b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+j))
      end do
    end if
    if(idmodel(yk).eq.1) then !linear
      if(idlink(yk).ne.0) then !lin or spl
        ll=0
        Ut(1+l+1,1+l+1)=1
      end if
      if(idlink(yk).eq.0) then !identity
        ll=1
        Ut(1+l+1,1+l+1)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+1)
      end if
      if(nea(yk).gt.1) then
        do j=2,nea(yk)
           do k=1,j
              Ut(1+l+j,1+l+k)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+ll+k-1+j*(j-1)/2)
           end do
        end do
      end if
    end if
    l = l+nea(yk)
    sumPrmYk = sumPrmYk + nef(yk) + nvc(yk) + ntr(yk) + 1
  end do
  ! Ut = matrice de Cholesky, partie triang sup remplie uniquement
 !print*,"Ut",Ut
  
  ! creation de Y1
  Y1=0.d0
  jacobien=0.d0
  splaa=0.d0
  
  sumPrmYk = 0
  sumMesYk = 0
  numSPL=0
  do yk=1,ny

     if (idlink(yk).eq.0) then  ! Identity link : Y1 = Y

        do j=1,nmes(i,yk)
           Y1(sumMesYk+j)=dble(Y(nmescur+sumMesYk+j))
           !jacobien
        end do
        
     else if (idlink(yk).eq.1) then  ! Linear link   
     
        do j=1,nmes(i,yk)
           Y1(sumMesYk+j)=(dble(Y(nmescur+sumMesYk+j))-b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+nvc(yk)+1)) &
                /abs(b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+nvc(yk)+2))
            
           jacobien = jacobien - log(b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+nvc(yk)+2))
        end do

     else if (idlink(yk).eq.2) then  ! Splines link
     
       numSPL=numSPL+1

       splaa=0.d0
       eta0=0.d0
       eta0=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+nvc(yk)+1)
       do kk=2,ntr(yk)
         splaa(kk-3)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+nvc(yk)+kk)**2
       end do
       !if(i==1 .and. id==0 .and. jd==0) print*,"eta0=",eta0,"splaa=",sqrt(splaa)
        
       do j=1,nmes(i,yk)
         ll=0
         !if(i==1 .and. id==0 .and. jd==0) print*,"Y=",Y(nmescur+sumMesYk+j)
         if (Y(nmescur+sumMesYk+j).eq.zitr(ntr(yk)-2,numSPL)) then
           ll=ntr(yk)-3
         end if

         som=0.d0
         do kk = 2,ntr(yk)-2
           if ((Y(nmescur+sumMesYk+j).ge.zitr(kk-1,numSPL)).and. &
                   (Y(nmescur+sumMesYk+j).lt.zitr(kk,numSPL))) then
                 ll=kk-1
           end if
         end do

         if (ll.lt.1.or.ll.gt.ntr(yk)-3) then   
           vrais_i=-1.d9
           !print*,"-1.d9 ll<1 ou ll>ntr-3",ll!," ntr=",ntr(yk)," numSPL=",numSPL," y=",Y(nmescur+sumMesYk+j)
           goto 654
         end if
         if (ll.gt.1) then
           do ii=2,ll
             som=som+splaa(ii-3)
           end do
         end if

         Y1(sumMesYk+j)=eta0 + som & 
            + splaa(ll-2)*im2(indiceY(nmescur+sumMesYk+j))&
            + splaa(ll-1)*im1(indiceY(nmescur+sumMesYk+j))&
            + splaa(ll)*im(indiceY(nmescur+sumMesYk+j))
            
         jacobien = jacobien &
            + log(splaa(ll-2)*mm2(indiceY(nmescur+sumMesYk+j)) &
                +splaa(ll-1)*mm1(indiceY(nmescur+sumMesYk+j))&
                +splaa(ll)*mm(indiceY(nmescur+sumMesYk+j)))

         !print*,"jac=",jacobien
         !print*,"ll =",ll
         !print*,"splaa =",splaa(ll-2)
         !print*,"nmescur+sumMesYk+j =",nmescur+sumMesYk+j
         !print*,"indiceY= ",indiceY(nmescur+sumMesYk+j)
         !print*,"mm =",mm(nmescur+sumMesYk+j)
         !print*,"mm1 =",mm1(nmescur+sumMesYk+j)
         !print*,"mm2 =",mm2(nmescur+sumMesYk+j)
         !print*,"Y =",Y(nmescur+sumMesYk+j)
         !write(*,*)'Y',Y1(sumMesYk+j),sumMesYk,yk,j,jacobien
       end do
       
     end if
     
     sumMesYk=sumMesYk+nmes(i,yk)
     sumPrmYk=sumPrmYk+nef(yk)+nvc(yk)+ntr(yk)+1
     
  end do !fin boucle yk
  !print*,'Y1',Y1
   
    
    
  ! creation des elements specifiques aux sous-modeles / trajectoires des marqueurs   
  
  mu_v=0.d0
  mu_r=0.d0
  asymptL=0.d0
  asymptU=0.d0
  
  b0_t=0.d0
  npm_t=0
  b0_x=0.d0
  X0_x=0.d0
  npm_x=0
  b0_tx=0.d0
  npm_tx=0
  
  sumMesYk=0
  sumPrmYk=0   
  do yk=1,ny
  
    if(idmodel(yk).eq.0) then !logistic
    
      ! location parameter : mean = mu_v
      mu_v(yk)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk))
       
      ! mu_r = taux de progression / rate = beta_r * X_r
      l=0
      b_r=0.d0
      X_r=0.d0
      do k=1,nv
        if(idg(yk,k).eq.1) then
          l=l+1
          X_r(l)=dble(X(nmescur+1,k)) !time-independant so 1st line of patient i
          b_r(l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+l)
        end if
      end do
      mu_r(yk)=DOT_PRODUCT(X_r,b_r)
       
      ! asymptotes
      asymptL(yk)=minY(yk) !lower
      asymptU(yk)=maxY(yk) !upper
    
    end if
    
    if(idmodel(yk).eq.1) then !linear
    
      !!! fixed effects : creation de X0_x et de b0_t, b0_x, b0_tx
       
      !!prm for time effect
      l=0
      do k=1,8
        if(idgfp(yk,k).ne.0) then
          do kk=1,idgfp(yk,k)
            l=l+1
            b0_t(yk,l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+l)
          end do
        end if
      end do
      npm_t(yk)=l   !sum(idgfp(yk,:)) !nb de prms time effect 
       
      
      !!prm and matrix for varexp
      l=0
      do k=1,nv
        if(idg(yk,k).eq.1) then
          !prm
          l=l+1
          b0_x(yk,l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+npm_t(yk)+l)
          !matrix
          do j=1,nmes(i,yk)
            X0_x(sumMesYk+j,l)=dble(X(nmescur+sumMesYk+j,k))
          end do
        end if
      end do
      npm_x(yk)=l   !sum(idg(yk,:)) !nb de prms varexp
      !write(*,*) 'b0_x',b0_x
      !write(*,*) 'X0_x',X0_x
      !write(*,*) 'npm_x',npm_x
 
       
      !!prm for interaction effect
      if(idginterac(yk).eq.1) then
        l=0
        do k=1,npm_t(yk)
          do kk=1,npm_x(yk)
            l=l+1
            b0_tx(yk,l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+npm_t(yk)+npm_x(yk)+l)
          end do
        end do
      end if
      npm_tx(yk)=l  !npm_t(yk)*npm_x(yk) !nb de prms d interaction
      
    end if
    
    !incrementation
    sumPrmYk=sumPrmYk+nef(yk)+nvc(yk)+ntr(yk)+1
    sumMesYk=sumMesYk+nmes(i,yk)
  
  end do
  
  !print*,"mu_v",mu_v
  !print*,"mu_r",mu_r
  !print*,"asymptL",asymptL
  !print*,"asymptU",asymptU
  
  !print*,"b0_t",b0_t
  !print*,"b0_x",b0_x
  !print*,"X0_x",X0_x
  !print*,"b0_tx",b0_tx
  !print*,"npm_t",npm_t
  !print*,"npm_x",npm_x
  !print*,"npm_tx",npm_tx
  
  
  !! creation elements du shift s
  
  ! mu_s = beta_s * X_s
  b_s=0.d0
  X_s=0.d0
  l=0
  do k=1,nv
     if (idgs(k).ne.0) then
        l=l+1
        X_s(l)=dble(X(nmescur+1,k))  !time-independant so 1st line of patient i  
        b_s(l)=b1(nrisqtot+nvarxevt+nasso+l)
     end if
  end do
  mu_s = DOT_PRODUCT(X_s,b_s)
  !print*,"b_s",b_s
  !print*,"X_s",X_s
  !print*,"mu_s",mu_s
  
  ! Xtime = vecteur du temps brut
  Xtime=0.d0
  do j=1,sum(nmes(i,:))
      Xtime(j)=dble(X(nmescur+j,1))  
  end do
  
  
  ! contribution individuelle a la vraisemblance
  ! print*,"i=",i," -ni*log(2pi)=",-sum(nmes(i,:))*dlog(dble(2*3.14159265)), " log(det)=",det
  ! print*,"Vi=",VC
  vrais_Y=0.d0
  vrais_i=0.d0
  
  
  som=0.d0
  som_T0=0.d0 !delayed entry
  som_Ti=0.d0
  
  !print*, "deb MCMC"
  
  do l=1,nMC

     !print*, "nMC", l
     !write(3,*)"nMC =",l
     
     vrais_Y=1.d0
     
     
     !!!!!!!!! MC pour EA !!!!!!!!!

     if(methInteg.eq.1) then 
        ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!

        ! simuler les effets aleatoires
        x22=0.d0
        SX=1.d0
        do j=1,neatot
          call bgos(SX,0,usim(j),x22,0.d0)
          !print*,"usim=",usim(j)
        end do
        ui=0.d0
        ui=matmul(Ut,usim)

     else if(methInteg.eq.2) then 
        ! !!!!!!!!!!!!! MCA !!!!!!!!!!!!!

        if(mod(l,2).eq.0) then
           ! si l est pair on prend l'oppose des precedents
           ui = -ui
        else
           ! sinon on simule des nouveaux
           ! simuler les effets aleatoires
           x22=0.d0
           SX=1.d0
           do j=1,neatot
            call bgos(SX,0,usim(j),x22,0.d0)
           end do
           ui=0.d0
           ui=matmul(Ut,usim)
        end if

     else 
        ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!

        ! simuler les effets aleatoires
        usim=0.d0
        do j=1,neatot
          usim(j)=seqMC(nMC*(j-1)+l)
        end do
        ui=0.d0
        ui=matmul(Ut,usim)

     end if ! fin if methInteg
     
     !print*,"methInteg",methInteg
     !print*,"Ut",Ut
     !print*,"usim",usim
     !print*,"ui", ui

     !write(3,*)"ui =",ui


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!! partie longitudinale
     
     !!!!!!!!!!!!!!!!
     !!!!! time-shift
     
     !! s
     shift = mu_s + ui(1)
     !! t+s
     time_shifted=0.d0
     do j=1,sum(nmes(i,:))
       time_shifted(j) = Xtime(j) + shift
     end do
     !print*, "shift", shift
     !print*, "time_shifted", time_shifted
     
     
     !!!!!!!!!!!!!!!!
     !!!!! outcome
     
     npm_ea_t=0
     Z_t=0.d0
     ui_Z_t=0.d0
     
     sumMesYk=0
     sumPrmYk=0
     sumRe=1 !repere ds vect ui

     do yk =1,ny
     
       !print*,"yk",yk

       if(nmes(i,yk).gt.0) then
          
         if(idmodel(yk).eq.0) then !logistic
         
           ! rate
           rate=0.d0
           rate=mu_r(yk)+ui(sumRe+1)
           !print*, "rate", rate
         
           ! nu = location prm
           nu=0.d0
           nu = mu_v(yk) + ui(sumRe+2)
           !print*, "nu", nu
           
           ! mu denominator
           denom = 0.d0
           do j=1,nmes(i,yk)
             denom(j) = 1 + exp( - rate*time_shifted(sumMesYk+j) )
             denom(j) = denom(j) ** nu
           end do
           !print*, "denom", denom
           
           ! mu
           mu = 0.d0
           do j=1,nmes(i,yk)
             mu(j) = asymptL(yk) + (asymptU(yk) - asymptL(yk)) / denom(j)
           end do
           !print*, "mu", mu
         
           ! + random effects
           ! mu + ea = ~Y = Y - error
           do j=1,nmes(i,yk)
             mu(j) = mu(j) + ui(sumRe+3)
           end do
           !print*, "mu+ea", mu
           
           
           !! var_Y = variance de Y|s,v,u
           ! avec inv_var_Y = matrice inverse de var_Y, comme var_Y diago, alors inverse = coefficient inverse
           ! avec det_var_Y = determinant de var_Y; comme var_Y diago, alors determinant = produit des coefficients diagonaux
           var_Y=0.d0
           inv_var_Y=0.d0
           det_var_Y=1
           do j=1,nmes(i,yk)
             var_Y(j,j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+nvc(yk)+ntr(yk)+1)**2 !std erreur de mesure
             inv_var_Y(j,j) = 1/var_Y(j,j)
             det_var_Y = det_var_Y * var_Y(j,j)
           end do
           !print*,"var_Y",var_Y
           !print*,"inv_var_Y",inv_var_Y
           !print*,"det_var_Y",det_var_Y
              
         end if
          
         if(idmodel(yk).eq.1) then !linear
          
           ! matrice de fct du shifted time X0_t
           X0_t=0.d0
           ll=0
           do k=1,8
             if(idgfp(yk,k).ne.0) then
               ll=ll+1
               do j=1,nmes(i,yk)
                 if(FPpower(k).eq.dble(0)) then
                   X0_t(j,ll)=log(time_shifted(sumMesYk+j)) ! log or ln ?
                 end if
                 if(FPpower(k).ne.dble(0)) then
                   X0_t(j,ll)=time_shifted(sumMesYk+j)**FPpower(k)
                 end if
               end do
               if(idgfp(yk,k).eq.2) then
                 ll=ll+1
                 do j=1,nmes(i,yk)
                   X0_t(j,ll)=X0_t(j,ll-1)*log(time_shifted(sumMesYk+j))
                 end do
               end if
             end if
           end do
           !print*,"X0_t",X0_t
             
           ! matrice des covariables du marqueur en cours X0_x_yk
           X0_x_yk=0.d0
           if(npm_x(yk).ne.0) then
             do k=1,nv
               do j=1,nmes(i,yk)
                 X0_x_yk(j,k)=X0_x(sumMesYk+j,k)
               end do
             end do
           end if
           !print*,"X0_x_yk",X0_x_yk   
             
           ! matrice d interaction entre fct du temps et covariables X0_tx
           X0_tx=0.d0
           if(npm_tx(yk).ne.0) then
            ll=0
            do k=1,npm_t(yk)
              do kk=1,npm_x(yk)
                ll=ll+1
                do j=1,nmes(i,yk)
                  X0_tx(j,ll)=X0_t(j,k)*X0_x(sumMesYk+j,kk)
                end do
              end do
            end do
           end if
           !print*,"X0_tx",X0_tx


           ! fixed effects
           mu = matmul(X0_t,b0_t(yk,:))
           if(npm_x(yk).ne.0) then
             mu = mu + matmul(X0_x_yk,b0_x(yk,:))
           end if
           if(npm_tx(yk).ne.0) then
             mu = mu + matmul(X0_tx,b0_tx(yk,:))
           end if
          !print*,"mu",mu


           ! matrice des EA sur le temps Z
           Z_t=0.d0
           ll=0
           do k=1,8
             if(ideafp(yk,k).ne.0) then
               ll=ll+1
               do j=1,nmes(i,yk)
                 if(FPpower(k).eq.dble(0)) then
                   Z_t(j,ll)=log(time_shifted(sumMesYk+j)) ! log or ln ? !here
                 end if
                 if(FPpower(k).ne.dble(0)) then
                   Z_t(j,ll)=time_shifted(sumMesYk+j)**FPpower(k)
                 end if
               end do
               if(ideafp(yk,k).eq.2) then
                 ll=ll+1
                 do j=1,nmes(i,yk)
                   Z_t(j,ll)=Z_t(j,ll-1)*log(time_shifted(sumMesYk+j))
                 end do
               end if
             end if
           end do
           npm_ea_t(yk)=ll
           !print*,"Z_t",Z_t
           
           
           ! + random effects
           ui_Z_t=0.d0
           ll=1
           if(npm_ea_t(yk).ne.0) then
             do k=1,npm_ea_t(yk)
               ui_Z_t(k)=ui(sumRe+idea(yk)+k)
             end do
           end if
          !print*,'ui_Z_t',ui_Z_t
           
           ! -> mu = esperance conditionnelle
           if(idea(yk).eq.1) mu = mu+ui(sumRe+1)
          !print*,"mu+intercept",mu
           if(npm_ea_t(yk).ne.0) mu = mu+matmul(Z_t,ui_Z_t)
          !print*,"mu+tous_ea",mu

           
           !! variance de Y|s,ui
           var_Y=0.d0
           do j1=1,nmes(i,yk)
             var_Y(j1,j1) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s & 
             +sumPrmYk+nef(yk)+nvc(yk)+ntr(yk)+1)**2 !variance de l'erreur yk
           end do
           !print*,"var_Y",var_Y
              
           ! var_Y en vecteur
           jj=0
           var_Y_vect=0.d0
           do j1=1,nmes(i,yk)
             do j2=j1,nmes(i,yk)
               jj=j1+j2*(j2-1)/2
               var_Y_vect(jj)=var_Y(j1,j2)
             end do
           end do
           !print*,"var_Y_vect",var_Y_vect
              
           ! inversion
           CALL dsinv(var_Y_vect,nmes(i,yk),eps,ier,det_var_Y)
           if (ier.eq.-1) then
             vrais_i=-1.d9
             !print*,"-1.d9 dsinv continu MC"
             !print*,"b=",b
             !print*,"bfix=",bfix
             !print*,"fix=",fix
             goto 654
           end if
           
           ! determinant en exponentiel
           det_var_Y = exp(det_var_Y)
              
           ! retransformation du vecteur var_Y_vect en matrice :
           inv_var_Y=0.d0
           do j1=1,nmes(i,yk)
             do j2=1,nmes(i,yk)
               if (j2.ge.j1) then
                 inv_var_Y(j1,j2)=var_Y_vect(j1+j2*(j2-1)/2)
               else
                 inv_var_Y(j1,j2)=var_Y_vect(j2+j1*(j1-1)/2)
               end if
             end do
           end do
           !print*,"var_Y",var_Y
           !print*,"inv_var_Y",inv_var_Y
           !print*,"det_var_Y",det_var_Y
           
         end if
         
         
         ! calcul de la vrais
         Y2=0.d0
         Y3=0.d0
         Y4=0.d0
         do j=1,nmes(i,yk)
           Y2(j) = Y1(sumMesYk+j) - mu(j)
         end do
        !print*,"Y1",Y1
        !print*,"mu",mu
        !print*,"Y2",Y2
              
         Y3=matmul(inv_var_Y,Y2)
        !print*,"Y3",Y3
              
         Y4=DOT_PRODUCT(Y2,Y3)
        !print*,"Y4",Y4
              
         div = (dble(2*3.14159265)**(dble(nmes(i,yk))/2))*sqrt(det_var_Y)
         !print*,"div",div
              
         vrais_Y = vrais_Y * exp(-Y4/2.d0)/div
         !print*,"density",exp(-Y4/2.d0)/div
         !print*,"vrais_Y",vrais_Y
         !print*,"vrais_Y_k",exp(-Y4/2.d0)/div
          
       end if
       
       !incrementation
       sumMesYk = sumMesYk+nmes(i,yk)
       sumPrmYk=sumPrmYk+nef(yk)+nvc(yk)+ntr(yk)+1
       sumRe=sumRe+nea(yk)
       
       !print*,"vrais_Y",vrais_Y
       
     end do ! fin boucle yk
     !if(i.lt.4) print*,"avant survie, vrais_Y",vrais_Y    
     !print*,"vrais_Y",vrais_Y  
        
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!! partie survie
     
     if (nbevt.ne.0) then
     
        !print*,"partie survie"

        ! calcul de brisq en chaque composante et risq, surv et surv0 pour chaque evt
        risq=0.d0
        surv=0.d0
        surv0=0.d0
        
        sumnrisq=0
        do ke=1,nbevt

           brisq=0.d0
           if (logspecif.eq.1) then
              do k=1,nprisq(ke)
                 brisq(k)=exp(b1(sumnrisq+k))
              end do
           else
              do k=1,nprisq(ke)
                 brisq(k)=b1(sumnrisq+k)*b1(sumnrisq+k)
              end do
           end if
           !print*,"brisq",brisq
           
           !print*,"idst",idst
           if(any(idst(ke,:).eq.1).OR.(nasso_s_fp.ne.0)) then !si au moins une asso par current value ou fct du tps en interaction avec shift ds survie, alors GK quadrature
           
             Tsurv_GK=0.d0
             
             centr=(tsurv(i)-0)/2
             hlgth=(tsurv(i)-0)/2
             !print*,"tsurv(i)",tsurv(i)
             !print*,"centr",centr
             !print*,"hlgth",hlgth
             !print*,"GKp",GKp
             
             do p=1,nGK
               Tsurv_GK(p)=centr+hlgth*GKp(p)
             end do
             !print*,"Tsurv_GK",Tsurv_GK
             
             Tsurv0_GK=0.d0
             if (idtrunc.eq.1) then
               centr0=(tsurv0(i)-0)/2
               hlgth0=(tsurv0(i)-0)/2
               do p=1,nGK
                 Tsurv0_GK(p)=centr0+hlgth0*GKp(p)
               end do
             end if
             
             !print*,"Tsurv0_GK",Tsurv0_GK
           
           
             ! creer base de splines si au moins un hazard splines
             ! specific to subject i
             if(typrisq(ke).eq.3) then
               call splines_i(i,ke) !temps reel
               call splines_GK_i(ke) !temps de quadrature
             end if
           
             call fct_risq_GK(i,shift,ke,brisq,sumnrisq,b1,ui,risq,surv,surv0)
             !risq = instant risq at event time Ti, with risq0(Ti)*exp(alpha*CV(Ti+s)+alpha'*s*delta(Ti))
             !surv = cumul risq computed with GK quadrature for event time Ti
             !surv0 = same for entry time T0i
             
           else
           
             ! creer base de splines si au moins un hazard splines
             ! specific to subject i
             if(typrisq(ke).eq.3) then
               call splines_i(i,ke)
             end if
           
             call fct_risq(i,ke,brisq,risq,surv,surv0)
           
           end if
           
           !print*,"risq",risq
           !print*,"surv",surv
           !print*,"surv0",surv0
           
           sumnrisq = sumnrisq + nprisq(ke) + nvarbyevt + nassoparevt(ke)
        end do

        ! variables explicatives de la survie
        Xevt=0.d0
        bevt=0.d0
        if (nxevt.ne.0) then    ! si varexpsurv

           m=0
           sumnrisq=0
           do ke=1,nbevt
              ll=0
              do k=1,nv 

                    if (idsurv(k).eq.1) then  
                       m=m+1
                       ll=ll+1
                       bevt(m)=b1(sumnrisq+nprisq(ke)+ll)
                       Xevt(m)=X(nmescur+1,k)
                    end if

              end do
              sumnrisq = sumnrisq + nprisq(ke) + nvarbyevt + nassoparevt(ke)
           end do

        end if
        
        
        ! vrais survie
        surv_glob=0.d0
        surv0_glob=0.d0
        varexpsurv=0.d0
        nxevtcurr=0
        fevt=0.d0
        easurv=0.d0
        sumnrisq=0
        do ke=1,nbevt
        
           !write(3,*)"ke =",ke
           !write(3,*)"nxevtcurr =",nxevtcurr
           !write(3,*)"sumnrisq =",sumnrisq
        
           ! calculer Xevt * bevt
           varexpsurv=0.d0
           if (nxevt.ne.0) then   !si varexpsurv
              varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevt))&
                   ,bevt((nxevtcurr+1):(nxevtcurr+nxevt)))
           end if
           !write(3,*)"varexpsurv =",varexpsurv
           
           !! association : shift
           ! partie sans interaction avec temps : intercept (si il y a)
           
           shiftsurv = 0.d0
           
           !intercept
           if(nasso_s.eq.1) shiftsurv = b1(sumnrisq+nprisq(ke)+nvarbyevt+1)
           shiftsurv = shiftsurv * shift
           !print*,"shift",shift
           !print*,"shiftsurv",shiftsurv
           !write(3,*)"shiftsurv =",shiftsurv
           
           
           !! association : effets aleatoires partages
           easurv=0.d0
           if(any(idst(ke,:).eq.0)) then
           
             sumasso = nasso_s + nasso_s_fp
             sumRe=1
             do yk=1,ny
             
                !write(3,*)"yk =",yk
                !write(3,*)"sumRe =",sumRe
                !write(3,*)"sumasso =",sumasso
             
               if(idst(ke,yk).eq.0) then !RE
               
                 ui_asso=0.d0
                 
                 if(idmodel(yk).eq.0) then !logistic
                   !ui_asso(1)= mu_r(yk) + ui(sumRe+1)
                   !ui_asso(2)= mu_v(yk) + ui(sumRe+2)
                   ui_asso(1)= ui(sumRe+1)
                   ui_asso(2)= ui(sumRe+2)
                   ui_asso(3)= ui(sumRe+3)
                 end if 
                 
                 if(idmodel(yk).eq.1) then !linear
                   ui_asso(1:nea(yk))= ui(sumRe+1:nea(yk))
                 end if 
                 
                 !write(3,*)"ui_asso =",ui_asso
                 
                 easurv = easurv + DOT_PRODUCT(ui_asso,b1(sumnrisq+nprisq(ke)+nvarbyevt+sumasso+1:nea(yk)))
                 
                 !write(3,*)"easurv =",easurv
                 
               end if  
                 
               !incrementation
               sumRe=sumRe+nea(yk)
               if(idst(ke,yk).eq.0) sumasso = sumasso + nea(yk)
               if(idst(ke,yk).eq.1) sumasso = sumasso + 1
               
             end do
           
           end if
           !print*,"easurv",easurv
           !write(3,*)"easurv =",easurv
           
           
           ! avoir evt au temps Ti si Devt=1
           if (Devt(i).eq.ke) then     !si sujet i a evt ke
              fevt=risq(ke)*exp(varexpsurv+shiftsurv+easurv)   !fct de risq         
           end if
           !print*,"Devt(i)",Devt(i)
           !print*,"fevt",fevt
           !write(3,*)"Devt(i) =",Devt(i)
           !write(3,*)"risq(ke) =",risq(ke)
           !write(3,*)"fevt =",fevt
           
           !write(3,*)"Surv_glob =",Surv_glob
           
           ! risque cumule jusque Ti
           Surv_glob=surv_glob + &
                  exp(varexpsurv+shiftsurv+easurv)*surv(ke)
           !print*,"Surv_glob",Surv_glob 
           !write(3,*)"surv(ke) =",surv(ke)
           !write(3,*)"Surv_glob =",Surv_glob
           
           ! troncature : risque cumule au temps T0
           if (idtrunc.eq.1) then
            surv0_glob=surv0_glob+surv0(ke)*exp(varexpsurv+shiftsurv+easurv)    
           end if
           
           nxevtcurr=nxevtcurr+nxevt
           sumnrisq = sumnrisq + nprisq(ke) + nvarbyevt + nassoparevt(ke)
           
        end do

        ! vraisemblance de la partie survie
        vrais_surv = exp(-Surv_glob)
        !print*,'vrais_surv',vrais_surv
        !write(3,*)"vrais_surv =",vrais_surv
        
        ! print*,"vrais_surv ok"        
        if(Devt(i).gt.0) vrais_surv = vrais_surv * fevt
        !write(3,*)"vrais_surv =",vrais_surv

        if (idtrunc.eq.1) then
           !vrais_surv = vrais_surv / exp(-surv0_glob)
           som_T0 = som_T0 + exp(-surv0_glob)   !delayed entry
        end if
        
        som_Ti = som_Ti + vrais_surv
        !print*,'som_Ti',som_Ti
        !print*,'vrais_surv',vrais_surv
        
        ! vrais totale 
        som = som + vrais_Y * vrais_surv
        !write(3,*)"som =",som

     
     else !  pas de survie

        som = som + vrais_Y

     end if
    
     !print*, "som", som
     
  end do ! fin boucle nMC
  
  !print*, "fin MCMC"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! contribution individuelle a vraisemblance
  
  vrais_i = vrais_i + log(som) - log(dble(nMC)) + jacobien
  
  if (idtrunc.eq.1) then    !delayedentry
      vrais_i = vrais_i - log(som_T0) + log(dble(nMC))
  end if
  
  !print*, "som", som
  !print*, "nMC", nMC
  !print*, "jacobien", jacobien
  !print*, "vrais_i", vrais_i
  
  
  654 continue

  return

end function vrais_i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      SPLINE BASES      !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------!
! I-SPLINES for link functions !
!------------------------------!

subroutine design_splines(ier)

  use modltsm

  implicit none

  integer ::jj,l,k,ier,yk,q,sumnval,qq,sumqqval
  double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht

  ier=0
  jj=0
  l=0
  q=0
  qq=0
  sumqqval=0
  sumnval=0
  do yk=1,ny
     if (idlink(yk).eq.2) then 
        q=q+1
        do jj=1,nvalSPL(q)      !     ou se trouve la valeur de zi

           do k = 2,ntr(yk)-2
              if ((uniqueY(sumqqval+sumnval+jj).ge.zitr(k-1,q)).and.(uniqueY(sumqqval+sumnval+jj).lt.zitr(k,q))) then
                 l=k-1
              end if
           End do


           if (uniqueY(sumqqval+sumnval+jj).eq.zitr(ntr(yk)-2,q)) then
              l=ntr(yk)-3
           end if

           ht2 = zitr(l+1,q)-uniqueY(sumqqval+sumnval+jj)
           htm= uniqueY(sumqqval+sumnval+jj)-zitr(l-1,q)
           ht = uniqueY(sumqqval+sumnval+jj)-zitr(l,q)
           ht3 = zitr(l+2,q)-uniqueY(sumqqval+sumnval+jj)
           hht = uniqueY(sumqqval+sumnval+jj)-zitr(l-2,q)
           h = zitr(l+1,q)-zitr(l,q)
           hh= zitr(l+1,q)-zitr(l-1,q)
           hn= zitr(l+1,q)-zitr(l-2,q)
           h2n=zitr(l+2,q)-zitr(l-1,q)
           h2= zitr(l+2,q)-zitr(l,q)
           h3= zitr(l+3,q)-zitr(l,q)

           if (uniqueY(sumqqval+sumnval+jj).ne.zitr(ntr(yk)-2,q)) then
              mm2(sumnval+jj) = (3.d0*ht2*ht2)/(hh*h*hn)
              mm1(sumnval+jj) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
              mm(sumnval+jj)  = (3.d0*ht*ht)/(h3*h2*h)

           end if
           if (uniqueY(sumqqval+sumnval+jj).eq.zitr(ntr(yk)-2,q)) then
              mm2(sumnval+jj) = 0.d0
              mm1(sumnval+jj) = 0.d0
              mm(sumnval+jj)  = 3.d0/h
           end if

           if (mm2(sumnval+jj).lt.0.or.mm1(sumnval+jj).lt.0.or.mm(sumnval+jj).lt.0) then
              ier=-1
              goto 765
           end if

           im2(sumnval+jj)=hht*mm2(sumnval+jj)/(3.d0)+ h2n*mm1(sumnval+jj)/(3.d0) &
                +h3*mm(sumnval+jj)/(3.d0)
           im1(sumnval+jj)=htm*mm1(sumnval+jj)/(3.d0)+h3*mm(sumnval+jj)/(3.d0)
           im(sumnval+jj)=ht*mm(sumnval+jj)/(3.d0)

        end do
        sumnval = sumnval + nvalSPL(q)

     end if

  end do


765 continue

end subroutine design_splines



!---------------------------------------!
! M-SPLINES for baseline risk functions !
!     in shared random effects case     !
!---------------------------------------!

subroutine splines_i(i,k)    ! specific to subject i
  use modltsm
  implicit none

  integer::k
  integer::i,kk,n,l
  double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n, &
       hn,hh2,hh3


  l=0
  Tmm=0.d0
  Tmm1=0.d0
  Tmm2=0.d0
  Tmm3=0.d0
  Tim=0.d0
  Tim1=0.d0
  Tim2=0.d0
  Tim3=0.d0
  Tmm0=0.d0
  Tmm01=0.d0
  Tmm02=0.d0
  Tmm03=0.d0
  Tim0=0.d0
  Tim01=0.d0
  Tim02=0.d0
  Tim03=0.d0
  Tmmt=0.d0
  Tmmt1=0.d0
  Tmmt2=0.d0
  Tmmt3=0.d0
  Timt=0.d0
  Timt1=0.d0
  Timt2=0.d0
  Timt3=0.d0

  zi(-2,k)=zi(1,k)
  zi(-1,k)=zi(1,k)
  zi(0,k)=zi(1,k)
  zi(nz(k)+1,k)=zi(nz(k),k)
  zi(nz(k)+2,k)=zi(nz(k),k)
  zi(nz(k)+3,k)=zi(nz(k),k)
  

  n=nz(k)+2
  !------------------- Tsurv ---------------------------
  
     do kk=2,n-2
        if ((Tsurv(i).ge.zi(kk-1,k)).and.  &
             Tsurv(i).lt.zi(kk,k)) then
           l=kk-1
        end if
     end do

     if (Tsurv(i).eq.zi(n-2,k)) then
        l=n-3
     end if

     ht = Tsurv(i)-zi(l,k)
     htm = Tsurv(i)-zi(l-1,k)
     h2t = Tsurv(i)-zi(l+2,k)
     ht2 = zi(l+1,k)-Tsurv(i)
     ht3 = zi(l+3,k)-Tsurv(i)
     hht = Tsurv(i)-zi(l-2,k)
     h = zi(l+1,k)-zi(l,k)
     hh = zi(l+1,k)-zi(l-1,k)
     h2 = zi(l+2,k)-zi(l,k)
     h3 = zi(l+3,k)-zi(l,k)
     h4 = zi(l+4,k)-zi(l,k)
     h3m = zi(l+3,k)-zi(l-1,k)
     h2n = zi(l+2,k)-zi(l-1,k)
     hn = zi(l+1,k)-zi(l-2,k)
     hh3 = zi(l+1,k)-zi(l-3,k)
     hh2 = zi(l+2,k)-zi(l-2,k)

     if (Tsurv(i).ne.zi(n-2,k)) then

        Tmm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        Tmm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
             +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
             +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        Tmm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
             +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
             +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        Tmm = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

     end if

     if (Tsurv(i).eq.zi(n-2,k)) then

        Tmm3 = 0.d0
        Tmm2 = 0.d0
        Tmm1 = 0.d0
        Tmm = 4.d0/h

     end if

     Tim3 = (0.25d0*(Tsurv(i)-zi(l-3,k))*Tmm3) &
          +(0.25d0*hh2*Tmm2)        &
          +(0.25d0*h3m*Tmm1)+(0.25d0*h4*Tmm)
     Tim2 = (0.25d0*hht*Tmm2)  &
          +(h3m*Tmm1*0.25d0)+(h4*Tmm*0.25d0)
     Tim1 = (htm*Tmm1*0.25d0)+(h4*Tmm*0.25d0)
     Tim = ht*Tmm*0.25d0

     !------------------- Tsurv0 --------------------------

     if (idtrunc.eq.1) then

        do kk=2,n-2
           if ((Tsurv0(i).ge.zi(kk-1,k)).and.   &
                Tsurv0(i).lt.zi(kk,k)) then
              l=kk-1
           end if
        end do

        if (Tsurv0(i).eq.zi(n-2,k)) then
           l=n-3
        end if

        ht = Tsurv0(i)-zi(l,k)
        htm = Tsurv0(i)-zi(l-1,k)
        h2t = Tsurv0(i)-zi(l+2,k)
        ht2 = zi(l+1,k)-Tsurv0(i)
        ht3 = zi(l+3,k)-Tsurv0(i)
        hht = Tsurv0(i)-zi(l-2,k)
        h = zi(l+1,k)-zi(l,k)
        hh = zi(l+1,k)-zi(l-1,k)
        h2 = zi(l+2,k)-zi(l,k)
        h3 = zi(l+3,k)-zi(l,k)
        h4 = zi(l+4,k)-zi(l,k)
        h3m = zi(l+3,k)-zi(l-1,k)
        h2n = zi(l+2,k)-zi(l-1,k)
        hn = zi(l+1,k)-zi(l-2,k)
        hh3 = zi(l+1,k)-zi(l-3,k)
        hh2 = zi(l+2,k)-zi(l-2,k)

        if (Tsurv0(i).ne.zi(nz(k)-2,k)) then

           Tmm03 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))

           Tmm02 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))   &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm01 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))    &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm0 = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (Tsurv0(i).eq.zi(n-2,k)) then

           Tmm03 = 0.d0
           Tmm02 = 0.d0
           Tmm01 = 0.d0
           Tmm0 = 4.d0/h

        end if

        Tim03 = (0.25d0*(Tsurv0(i)-zi(l-3,k))*Tmm03)  &
             +(0.25d0*hh2*Tmm02)           &
             +(0.25d0*h3m*Tmm01)+(0.25d0*h4*Tmm0)
        Tim02 = (0.25d0*hht*Tmm02)                  &
             +(h3m*Tmm01*0.25d0)+(h4*Tmm0*0.25d0)
        Tim01 = (htm*Tmm01*0.25d0)+(h4*Tmm0*0.25d0)
        Tim0 = ht*Tmm0*0.25d0

     end if

        Timt3 = Tim3
        Timt2 = Tim2
        Timt1 = Tim1
        Timt = Tim


end subroutine splines_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------!
! M-SPLINES for baseline risk functions !
!     in shared current value case      !
!---------------------------------------!

!bases de splines pr ts les noeuds de la quadrature de Gauss-Kronrod
subroutine splines_GK_i(k)    ! specific to subject i
  use modltsm
  implicit none

  integer::k
  integer::kk,n,l
  integer::p 
  double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n, &
       hn,hh2,hh3


  l=0
  Tmm_GK=0.d0
  Tmm1_GK=0.d0
  Tmm2_GK=0.d0
  Tmm3_GK=0.d0
  Tmm0_GK=0.d0
  Tmm01_GK=0.d0
  Tmm02_GK=0.d0
  Tmm03_GK=0.d0
  
  !zi deja mis en forme a l etape avec splines_i

  n=nz(k)+2
  
  do p=1,nGK
  
     !------------------- Tsurv ---------------------------

     !! HERE

     do kk=2,n-2
        if ((Tsurv_GK(p).ge.zi(kk-1,k)).and.  &
             Tsurv_GK(p).lt.zi(kk,k)) then
           l=kk-1
        end if
     end do

     if (Tsurv_GK(p).eq.zi(n-2,k)) then
        l=n-3
     end if

     ht = Tsurv_GK(p)-zi(l,k)
     htm = Tsurv_GK(p)-zi(l-1,k)
     h2t = Tsurv_GK(p)-zi(l+2,k)
     ht2 = zi(l+1,k)-Tsurv_GK(p)
     ht3 = zi(l+3,k)-Tsurv_GK(p)
     hht = Tsurv_GK(p)-zi(l-2,k)
     h = zi(l+1,k)-zi(l,k)
     hh = zi(l+1,k)-zi(l-1,k)
     h2 = zi(l+2,k)-zi(l,k)
     h3 = zi(l+3,k)-zi(l,k)
     h4 = zi(l+4,k)-zi(l,k)
     h3m = zi(l+3,k)-zi(l-1,k)
     h2n = zi(l+2,k)-zi(l-1,k)
     hn = zi(l+1,k)-zi(l-2,k)
     hh3 = zi(l+1,k)-zi(l-3,k)
     hh2 = zi(l+2,k)-zi(l-2,k)

     if (Tsurv_GK(p).ne.zi(n-2,k)) then

        Tmm3_GK(p) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        Tmm2_GK(p) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
             +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
             +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        Tmm1_GK(p) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
             +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
             +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        Tmm_GK(p) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

     end if

     if (Tsurv_GK(p).eq.zi(n-2,k)) then

        Tmm3_GK(p) = 0.d0
        Tmm2_GK(p) = 0.d0
        Tmm1_GK(p) = 0.d0
        Tmm_GK(p) = 4.d0/h

     end if
     

     !------------------- Tsurv0 --------------------------

     if (idtrunc.eq.1) then

        do kk=2,n-2
           if ((Tsurv0_GK(p).ge.zi(kk-1,k)).and.   &
                Tsurv0_GK(p).lt.zi(kk,k)) then
              l=kk-1
           end if
        end do

        if (Tsurv0_GK(p).eq.zi(n-2,k)) then
           l=n-3
        end if

        ht = Tsurv0_GK(p)-zi(l,k)
        htm = Tsurv0_GK(p)-zi(l-1,k)
        h2t = Tsurv0_GK(p)-zi(l+2,k)
        ht2 = zi(l+1,k)-Tsurv0_GK(p)
        ht3 = zi(l+3,k)-Tsurv0_GK(p)
        hht = Tsurv0_GK(p)-zi(l-2,k)
        h = zi(l+1,k)-zi(l,k)
        hh = zi(l+1,k)-zi(l-1,k)
        h2 = zi(l+2,k)-zi(l,k)
        h3 = zi(l+3,k)-zi(l,k)
        h4 = zi(l+4,k)-zi(l,k)
        h3m = zi(l+3,k)-zi(l-1,k)
        h2n = zi(l+2,k)-zi(l-1,k)
        hn = zi(l+1,k)-zi(l-2,k)
        hh3 = zi(l+1,k)-zi(l-3,k)
        hh2 = zi(l+2,k)-zi(l-2,k)

        if (Tsurv0_GK(p).ne.zi(nz(k)-2,k)) then

           Tmm03_GK(p) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))

           Tmm02_GK(p) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))   &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm01_GK(p) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))    &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm0_GK(p) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (Tsurv0_GK(p).eq.zi(n-2,k)) then

           Tmm03_GK(p) = 0.d0
           Tmm02_GK(p) = 0.d0
           Tmm01_GK(p) = 0.d0
           Tmm0_GK(p) = 4.d0/h

        end if

     end if
     
    end do !end p
  
  
end subroutine splines_GK_i



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      SURVIVAL ELEMENTS      !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------!
!     baseline and cumulative RISK FUNCTIONS    !
!         in shared random effects case         !
!-----------------------------------------------!

subroutine fct_risq(i,k,brisq,risq,surv,surv0)

  use modltsm
        
  implicit none
  
  integer::i,k
  double precision,dimension(nprisq(k))::brisq
  double precision,dimension(nbevt)::risq,surv,surv0
  
  integer::j,l,ll,kk,ii
  double precision::som
  
  if (typrisq(k).eq.2.and.logspecif.eq.1) then
     
     surv(k)=brisq(1)*(tsurv(i)-zi(1,k))**brisq(2)      !zi(1,k)=depart Weibull
     
     risq(k)=brisq(1)*brisq(2)*(tsurv(i)-zi(1,k))**(brisq(2)-1)
     
     if (idtrunc.eq.1) then
        surv0(k)=brisq(1)*(tsurv0(i)-zi(1,k))**brisq(2)
     end if
     
  end if
  
  if (typrisq(k).eq.2.and.logspecif.eq.0) then
     
     surv(k)=(brisq(1)*(tsurv(i)-zi(1,k)))**brisq(2)
     
     risq(k)=brisq(1)*brisq(2)*(brisq(1)*(tsurv(i)-zi(1,k)))**(brisq(2)-1)
     
     if (idtrunc.eq.1) then
        surv0(k)=(brisq(1)*(tsurv0(i)-zi(1,k)))**brisq(2)
     end if
     
  end if
 
  if (typrisq(k).eq.1) then
     do j=1,nz(k)-1
        som=0.d0
        do l=1,j-1
           som=som+brisq(l)*(zi(l+1,k)-zi(l,k))
        end do
        if (idtrunc.eq.1) then
           if (Tsurv0(i).ge.zi(j,k).and.Tsurv0(i).le.zi(j+1,k)) then
              surv0(k)=som+brisq(j)*(Tsurv0(i)-zi(j,k))
           end if
        end if
        if (Tsurv(i).ge.zi(j,k).and.Tsurv(i).le.zi(j+1,k)) then
           surv(k)=som+brisq(j)*(Tsurv(i)-zi(j,k))
           risq(k)=brisq(j)
        end if
     end do
  end if
  
  

  if (typrisq(k).eq.3) then
     !------------ survie et risq pour Tsurv ----------------
     ll=0
     if (Tsurv(i).eq.zi(nz(k),k)) then
        ll=nz(k)-1
     end if
     som=0.d0
     do kk=2,nz(k)
        if ((Tsurv(i).ge.zi(kk-1,k)).and.(Tsurv(i).lt.zi(kk,k))) &
             then
           ll=kk-1
        end if
     end do
     if (ll.gt.1) then
        do ii=1,ll-1
           som=som+brisq(ii)
        end do
     end if
     
     surv(k)=som+brisq(ll)*Tim3+brisq(ll+1)*Tim2 &
          +brisq(ll+2)*Tim1+brisq(ll+3)*Tim
     risq(k)=brisq(ll)*Tmm3+brisq(ll+1)*Tmm2     &
          +brisq(ll+2)*Tmm1+brisq(ll+3)*Tmm
     

     !------------ survie et risq pour Tsurv0 ----------------
     
     if (idtrunc.eq.1) then
        ll=0
        if (Tsurv0(i).eq.zi(nz(k),k)) then
           ll=nz(k)-1
        end if
        som=0.d0
        do kk=2,nz(k)
           if ((Tsurv0(i).ge.zi(kk-1,k)).and.(Tsurv0(i).lt.zi(kk,k))) &
                then
              ll=kk-1
           end if
        end do
        !               if (ll.lt.1.or.ll.gt.nz-1) then
        !                  write(*,*) 'probleme dans fct_risq splines'
        !                  write(*,*) 'll=',ll,'T=',Tsurv0(i)
        !                  stop
        !               end if
        if (ll.gt.1) then
           do ii=1,ll-1
              som=som+brisq(ii)
           end do
        end if
        
        surv0(k)=som+brisq(ll)*Tim03+brisq(ll+1)*Tim02 &
             +brisq(ll+2)*Tim01+brisq(ll+3)*Tim0
        
     end if
     
  end if
  
  
end subroutine fct_risq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------!
!     baseline and cumulative RISK FUNCTIONS    !
!    including time-varying linear predictor    !
!          in shared current level case         !
!-----------------------------------------------!

subroutine fct_risq_GK(i,s,k,brisq,sumnrisq,b1,ui,risq,surv,surv0)


  use modltsm
        
  implicit none
  
  integer::i,k,sumnrisq
  double precision :: s
  double precision,dimension(nprisq(k))::brisq
  double precision,dimension(npmtot)::b1
  double precision,dimension(neatot)::ui
  double precision,dimension(nbevt)::risq,surv,surv0
  
  ! declaration objets
  integer::yk,p,sumasso,sumMesYk
  double precision::predCV_event 
  double precision,dimension(nGK)::risq0GK_event,risq0GK_entry
  double precision,dimension(nGK)::risqGK_event,risqGK_entry
  double precision,dimension(ny,nGK)::predCV_GK_event,predCV_GK_entry
  double precision,dimension(nGK)::predlinCV_GK_event,predlinCV_GK_entry
  
  double precision :: delta_shift
  double precision, dimension(nGK) :: delta_shift_GK,delta_shift0_GK
  
  double precision,external::fct_risq_base,fct_predCV,fct_delta_shift

  
  !!! baseline risq at event time Ti 
  risq(k) = fct_risq_base(Tsurv(i),k,brisq,0,0) 
  !entry = 0 pr event time, GKpoint = 0 pr tps reel
  !print*,"risq(k) de base",risq(k)
  
  !! marker current value prediction at event time Ti + si
  predCV_event = 0.d0
  
  sumasso = nasso_s + nasso_s_fp
  sumMesYk = 0
  do yk=1,ny
  
     if(idst(k,yk).eq.1) then !CV
       
       ! prediction, alpha*CV(t), puis somme
       predCV_event = predCV_event &
         + fct_predCV(i,s,yk,sumMesYk,ui,0,0) * b1(sumnrisq+nprisq(k)+nvarbyevt+sumasso+1) !here
       
       !print*,"yk",yk
       !print*,"sumasso",sumasso
       !print*,"passoCV",b1(sumnrisq+nprisq(k)+nvarbyevt+sumasso+1)
       !print*,"predCV",fct_predCV(i,s,yk,sumMesYk,ui,0,0)
       !print*,"predCV_event",predCV_event
       
     end if
     
     !incrementation
     sumMesYk=sumMesYk+nmes(i,yk)
     if(idst(k,yk).eq.0) sumasso = sumasso + nea(yk)
     if(idst(k,yk).eq.1) sumasso = sumasso + 1
     
  end do
  !print*,"predCV_event",predCV_event
  ! exp(somme)
  predCV_event = EXP(predCV_event)
  !print*,"predCV_event en exp",predCV_event
  
  
  !! shift * delta(T)
  delta_shift = 0.d0
  
  if(nasso_s_fp.ne.0) then
       delta_shift = fct_delta_shift(i,k,sumnrisq,b1,0,0)
  end if
  !print*,"delta_shift",delta_shift
  
  
  !! instant risq
  risq(k) = risq(k) * ( EXP(s*delta_shift) * predCV_event )
  !print*,"risq(k)",risq(k)
  
  
  !!! cumulative risq at GK quadrature time points
  
  ! baseline risq
  do p=1,nGK
     risq0GK_event(p) = fct_risq_base(Tsurv_GK(p),k,brisq,0,p)
     if(idtrunc.eq.1) risq0GK_entry(p) = fct_risq_base(Tsurv0_GK(p),k,brisq,1,p)
  end do
  !print*,'Tsurv_GK',Tsurv_GK
  !print*,'risq0GK_event',risq0GK_event
  
  
  !! marker current value prediction
  predCV_GK_event = 0.d0  !matrix with 1 column per GK point and 1 line per marker
  predCV_GK_entry = 0.d0  
  sumasso = nasso_s + nasso_s_fp
  sumMesYk = 0
  do yk=1,ny
  
     if(idst(k,yk).eq.1) then
       
       ! prediction
       do p=1,nGK
         predCV_GK_event(yk,p) = fct_predCV(i,s,yk,sumMesYk,ui,0,p) 
         if(idtrunc.eq.1) predCV_GK_entry(yk,p) = fct_predCV(i,s,yk,sumMesYk,ui,1,p)
       end do
       
       ! alpha*CV(t)
       predCV_GK_event(yk,:) = predCV_GK_event(yk,:) * b1(sumnrisq+nprisq(k)+nvarxevt+sumasso+1)
       if(idtrunc.eq.1) then
         predCV_GK_entry(yk,:) = predCV_GK_entry(yk,:) * b1(sumnrisq+nprisq(k)+nvarxevt+sumasso+1)
       end if
       
     end if
     
     !incrementation
     sumMesYk=sumMesYk+nmes(i,yk)
     if(idst(k,yk).eq.0) sumasso = sumasso + nea(yk)
     if(idst(k,yk).eq.1) sumasso = sumasso + 1
     
  end do
  
  
  !! shift * delta(T)
  delta_shift_GK = 0.d0
  delta_shift0_GK = 0.d0
  
  if(nasso_s_fp.ne.0) then
     do p=1,nGK
       delta_shift_GK(p) = fct_delta_shift(i,k,sumnrisq,b1,0,p)
       if(idtrunc.eq.1) delta_shift0_GK(p) = fct_delta_shift(i,k,sumnrisq,b1,1,p)
     end do
  end if
  !print*,'delta_shift_GK',delta_shift_GK
  
  ! exp(somme)
  do p=1,nGK
     predlinCV_GK_event(p) = EXP(sum(predCV_GK_event(:,p))) * EXP( s * delta_shift_GK(p))
     if(idtrunc.eq.1) predlinCV_GK_entry(p) = EXP(sum(predCV_GK_entry(:,p)) * EXP( s * delta_shift0_GK(p)))
  end do
  !print*,'predlinCV_GK_event',predlinCV_GK_event
  
  !! instant risq
  do p=1,nGK
    risqGK_event(p) = risq0GK_event(p)*predlinCV_GK_event(p) !DOT_PRODUCT(risq0GK_event,predlinCV_GK_event) !!!!here
    if(idtrunc.eq.1) risqGK_entry(p) = risq0GK_entry(p)*predlinCV_GK_entry(p) !DOT_PRODUCT(risq0GK_entry,predlinCV_GK_entry)
  end do
  !print*,'risq0GK_event',risq0GK_event
  !print*,'predlinCV_GK_event',predlinCV_GK_event
  !print*,'risqGK_event',risqGK_event
  
  !! GK approximation
  ! ponderation
  do p=1,nGK
     risqGK_event(p)= GKw(p) * risqGK_event(p)
     if(idtrunc.eq.1) risqGK_entry(p)= GKw(p) * risqGK_entry(p)
  end do
  !print*,'risqGK_event pond',risqGK_event
  
  ! cumulative risq
  surv(k) = sum(risqGK_event) * hlgth
  if(idtrunc.eq.1) surv0(k) = sum(risqGK_entry) * hlgth0
  !print*,'risqGK_event',risqGK_event
  !print*,'sum(risqGK_event)',sum(risqGK_event)
  !print*,'hlgth',hlgth
  !print*,'surv(k)',surv(k)

end subroutine fct_risq_GK



!-----------------------------------------------!
!             baseline RISK FUNCTION            !
!          in shared current level case         !
!-----------------------------------------------!

!!! baseline risk function
! output : 1 value

double precision function fct_risq_base(t,k,brisq,entry,GKpoint)  
   
   !t time
   !k event
   !brisq prm de base
   !entry = 0 pr event time, 1 pr entry time  !necessary for splines
   !GKpoint = 0 pr tps reel, p=1,...,nGK pr tps de quadrature GK  !necessary for splines
   
   use modltsm
        
   implicit none
  
   double precision::t
   integer::k,entry,GKpoint
   double precision,dimension(nprisq(k))::brisq
   
   integer::j,kk,ll
   double precision::risq0 
   
   
   risq0=0.d0
   
   
   if (typrisq(k).eq.2.and.logspecif.eq.1) then     !Weibull & exponentiel
     
     risq0=brisq(1)*brisq(2)*(t-zi(1,k))**(brisq(2)-1)   !fct risq base au tps t, zi(1,k)=depart Weibull
     
  end if
  
  if (typrisq(k).eq.2.and.logspecif.eq.0) then         !Weibull & quadratique

     risq0=brisq(1)*brisq(2)*(brisq(1)*(t-zi(1,k)))**(brisq(2)-1)
     
  end if
 
  if (typrisq(k).eq.1) then        !piecewise
     do j=1,nz(k)-1
        if (t.ge.zi(j,k).and.t.le.zi(j+1,k)) then
           risq0=brisq(j)
        end if
     end do
  end if
  
  if (typrisq(k).eq.3) then          !splines
     ll=0
     if (t.eq.zi(nz(k),k)) then
        ll=nz(k)-1
     end if
     do kk=2,nz(k)
        if ((t.ge.zi(kk-1,k)).and.(t.lt.zi(kk,k))) &
             then
           ll=kk-1
        end if
     end do

     if(GKpoint.eq.0) then  !real time
        if (entry.eq.0) then  !event time
          risq0=brisq(ll)*Tmm3+brisq(ll+1)*Tmm2     &
            +brisq(ll+2)*Tmm1+brisq(ll+3)*Tmm
        else if (entry.eq.2) then  !entry time
          risq0=brisq(ll)*Tmm03+brisq(ll+1)*Tmm02     &
            +brisq(ll+2)*Tmm01+brisq(ll+3)*Tmm0
        end if
     else   !time of GK quadrature
        if (entry.eq.0) then  !event time
          risq0=brisq(ll)*Tmm3_GK(GKpoint)+brisq(ll+1)*Tmm2_GK(GKpoint)     &
            +brisq(ll+2)*Tmm1_GK(GKpoint)+brisq(ll+3)*Tmm_GK(GKpoint)
        else if (entry.eq.1) then   !entry time
          risq0=brisq(ll)*Tmm03_GK(GKpoint)+brisq(ll+1)*Tmm02_GK(GKpoint)     &
            +brisq(ll+2)*Tmm01_GK(GKpoint)+brisq(ll+3)*Tmm0_GK(GKpoint)
        end if
     end if
     
  end if
  
  
  fct_risq_base = risq0
   

end function fct_risq_base



!-----------------------------------------------!
!    latent process current level predictions   !
!             at each GK time points            !
!          in shared current level case         !
!-----------------------------------------------!

!!! marker(s)' current value prediction
! output : 1 value

double precision function fct_predCV(i,s,yk,sumMesYk,ui,entry,GKpoint)  
   
   !i subject
   !s shift
   !yk marker
   !sumMesYk nb obs before marker yk
   !ui REs
   !entry = 0 pr event time, 1 pr entry time
   !GKpoint = 0 pr tps reel, p=1,...,nGK pr tps de quadrature GK
   
   use modltsm
        
   implicit none
  
   integer::i,yk,sumMesYk,entry,GKpoint
   double precision :: s
   double precision,dimension(neatot)::ui
   
   integer::k,kk,l,ll
   double precision::t
   
   double precision,dimension(2*8)::v_t,vz_t
   double precision,dimension(nv) :: v_x
   double precision,dimension(nv*(2*8)) :: v_tx
   
   double precision::pred 
   
   
   pred=0.d0
   
   !! time of prediction
   t=0.d0
   if(entry.eq.0 .and. GKpoint.eq.0) t=Tsurv(i)
   if(entry.eq.0 .and. GKpoint.ne.0) t=Tsurv_GK(GKpoint)
   if(entry.eq.1 .and. GKpoint.eq.0) t=Tsurv0(i)
   if(entry.eq.1 .and. GKpoint.ne.0) t=Tsurv0_GK(GKpoint)
   !shifted
   t = t+s
   
   !! marker value prediction
   
   if(idmodel(yk).eq.0) then  !logistic
   
     pred = ( 1 + exp( - rate*t ) ) ** nu   ! denominator
     pred = asymptL(yk) + (asymptU(yk) - asymptL(yk)) / pred   ! mean
     pred = pred + ui(1+(yk-1)*sum(nea(1:(yk-1)))+3)   ! + random effects    ! here ok ??

   end if
   
   if(idmodel(yk).eq.1) then  !linear
   
     ! vecteur de fct du temps X0_t
     v_t=0.d0
     ll=0
     do k=1,8
       if(idgfp(yk,k).ne.0) then
         ll=ll+1
         if(FPpower(k).eq.dble(0)) then
           v_t(ll)=log(t) ! log or ln ? !here
         end if
         if(FPpower(k).ne.dble(0)) then
           v_t(ll)=t**FPpower(k)
         end if
         if(idgfp(yk,k).eq.2) then
           ll=ll+1
           v_t(ll)=v_t(ll-1)*log(t)
         end if
       end if
     end do
     !print*,"v_t",v_t
     
     ! vecteur des covariables du marqueur au dernier temps observe
     v_x=0.d0
     if(npm_x(yk).ne.0) then
       do k=1,nv
         v_x(k)=X0_x(sumMesYk+nmes(i,yk),k)
       end do
     end if
     !print*,"v_x",v_x   
     
     ! vecteur d interaction entre fct du temps et covariables
     v_tx=0.d0
     if(npm_tx(yk).ne.0) then
       ll=0
       do k=1,npm_t(yk)
         do kk=1,npm_x(yk)
           ll=ll+1
           v_tx(ll)=v_t(k)*v_x(kk)
         end do
       end do
     end if
     !print*,"v_tx",v_tx


     ! prediction
     pred = DOT_PRODUCT(v_t,b0_t(yk,:))   !here
     if(npm_x(yk).ne.0) then
       pred = pred + DOT_PRODUCT(v_x,b0_x(yk,:))
     end if
     if(npm_tx(yk).ne.0) then
       pred = pred + DOT_PRODUCT(v_tx,b0_tx(yk,:))
     end if
     !print*,"pred fx",pred
     
     
     ! vecteur des EA sur le temps Z
     vz_t=0.d0
     ll=0
     do k=1,8
       if(ideafp(yk,k).ne.0) then
         ll=ll+1
         if(FPpower(k).eq.dble(0)) then
           vz_t(ll)=log(t) ! log or ln ? !here
         end if
         if(FPpower(k).ne.dble(0)) then
           vz_t(ll)=t**FPpower(k)
         end if
         if(ideafp(yk,k).eq.2) then
           ll=ll+1
           vz_t(ll)=vz_t(ll-1)*log(t)
         end if
       end if
     end do
     !print*,"vz_t",vz_t
           
     ! prediction + random effects
     if(idea(yk).ne.0) then
       pred = pred + ui(1+(yk-1)*sum(nea(1:(yk-1)))+1)
     end if
     !print*,"ind_ui_int_ea",1+(yk-1)*sum(nea(1:(yk-1)))+1
     !print*,"ui_int_ea",ui(1+(yk-1)*sum(nea(1:(yk-1)))+1)
     !print*,"pred fx + int_ea",pred
     
     if(npm_ea_t(yk).ne.0) then
       do k=1,npm_ea_t(yk)
         pred = pred + vz_t(k) * ui(1+(yk-1)*sum(nea(1:(yk-1)))+idea(yk)+k)
       end do
     end if

   end if
   !print*,"pred fx + ea",pred
   
   
   fct_predCV = pred
   

end function fct_predCV



!------------------------------------------------------------!
!                       Delta function                       !
!                in interaction with the shift               !
!        as linear predictor in the survival submodel        !
!------------------------------------------------------------!

!!! delta function value
! output : 1 value

double precision function fct_delta_shift(i,k,sumnrisq,b1,entry,GKpoint)  


   !i subject
   !k event
   !sumnrisq nb of survival prms before event k
   !b1 vector of parameters
   !entry = 0 pr event time, 1 pr entry time
   !GKpoint = 0 pr tps reel, p=1,...,nGK pr tps de quadrature GK
   
   use modltsm
        
   implicit none
  
   integer::i,k,sumnrisq,entry,GKpoint
   double precision,dimension(npmtot)::b1
   
   integer::l,ll
   double precision::t
   
   double precision,dimension(2*8)::v_t
   
   double precision::delta 
   
   
   delta=0.d0
   
   
   !! time of prediction
   t=0.d0
   if(entry.eq.0 .and. GKpoint.eq.0) t=Tsurv(i)
   if(entry.eq.0 .and. GKpoint.ne.0) t=Tsurv_GK(GKpoint)
   if(entry.eq.1 .and. GKpoint.eq.0) t=Tsurv0(i)
   if(entry.eq.1 .and. GKpoint.ne.0) t=Tsurv0_GK(GKpoint)
   
   
   !! delta function value
   
   ! vecteur de fct du temps X0_t
   v_t=0.d0
   if(nasso_s_fp.ne.0) then
     ll=0
     do l=1,8
       if(idsurvsfp(l).ne.0) then
         ll=ll+1
         if(FPpower(l).eq.dble(0)) then
           v_t(ll)=log(t) ! log or ln ? !here
         end if
         if(FPpower(l).ne.dble(0)) then
           v_t(ll)=t**FPpower(l)
         end if
         if(idsurvsfp(l).eq.2) then
           ll=ll+1
           v_t(ll)=v_t(ll-1)*log(t)
         end if
       end if
     end do
   end if
   !print*,"v_t",v_t
     
   ! value computation
   if(nasso_s_fp.ne.0) delta = delta + DOT_PRODUCT(v_t,b1(sumnrisq+nprisq(k)+nvarbyevt+nasso_s+1:nasso_s_fp))
   
   fct_delta_shift = delta

end function fct_delta_shift

