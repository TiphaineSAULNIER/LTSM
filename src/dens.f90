
subroutine dens(re0 &
     ,best0 &
     ,Y0,X0,Tentr0,Tevt0,Devt0 &
     ,ideafp0,idgfp0,idea0,idg0,idginterac0,idgs0,idsurv0,idsurvs0,idsurvsfp0,sharedtype0,nGK0 &
     ,typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0 &
     ,ny0,ns0,nv0 &
     ,nobs0,nmes0,npm0,nea0,nfix0,refix0,fix0 &
     ,idmodel0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0,nvalSPLORD0 &
     ,ind0,nmescur0 &
     ,dens_res)

  use modltsm

  IMPLICIT NONE
  
  !Declaration des variables en entree
  integer,intent(in)::nv0,ny0
  integer, intent(in)::ns0,nobs0,npm0,nea0,nfix0
  integer, dimension(nbevt0,ny0), intent(in)::sharedtype0
  integer, intent(in)::nGK0
  integer,intent(in)::idtrunc0,logspecif0,nbevt0
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
  integer, dimension(ny0,8),intent(in)::ideafp0,idgfp0
  integer, dimension(ny0),intent(in)::idea0,idginterac0
  integer, dimension(ny0,nv0),intent(in)::idg0
  integer, dimension(nv0),intent(in)::idgs0
  integer,dimension(ns0,ny0),intent(in)::nmes0   
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  double precision,dimension(nea0),intent(in)::re0
  double precision, dimension(npm0), intent(in) :: best0
  integer,intent(in)::ind0,nmescur0
  integer,dimension(nea0+nfix0),intent(in)::fix0
  double precision, dimension(nfix0), intent(in) :: refix0
  
  !Declaration des variables en sortie
  double precision,intent(out)::dens_res
  
  !Variables locales
  integer::jtemp,i,j,ier,k,ktemp,yk,k1,mi,nbfix,l,m
  integer::ke,neatot0
  double precision, dimension(nea0+nfix0)::re1
  double precision::eps
  double precision,external::pred
  
  
  allocate(bp(npm0)) !vecteur des prms
  do k=1,npm0
     bp(k)=best0(k)
  end do
  !print*,'bp',bp
  
  
  !print*,"ELEMENTS EN ENTREE FORTRAN"
  !print*,"re0",re0
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
  !print*,"nea0",nea0
  !print*,"idmodel0",idmodel0
  !print*,"idlink0",idlink0
  !print*,"nbzitr0",nbzitr0
  !print*,"zitr0",zitr0
  !print*,"uniqueY0",uniqueY0
  !print*,"indiceY0",indiceY0
  !print*,"nvalSPLORD0",nvalSPLORD0
  !print*,"ind0",ind0
  !print*,"nmescur0",nmescur0
  !print*,"dens_res",dens_res
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
    if(any(typrisq0.eq.3)) then
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
    
    if(any(idst.eq.1)) then
    
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
  
  
  ! nb prms
  !npm0
  
  
  ! prm fixes
  neatot0 = nea0+nfix0
  allocate(fix(neatot0))
  fix=0
  fix(1:neatot0)=fix0(1:neatot0)
  nbfix=sum(fix(:))
  if(nbfix.eq.0) then
     allocate(refix(1))
  else
     allocate(refix(nbfix))
  end if
  refix(1:nbfix)=refix0(1:nbfix)
  
  re1=0.d0
  l=0
  m=0
  do k=1,neatot0
     if(fix(k).eq.0) then
        l=l+1
        re1(k)=re0(l)
     end if
     if(fix(k).eq.1) then
        m=m+1
        re1(k)=refix(m)
     end if
  end do
  
  !write(*,*)'re1',re1

  
  
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
      nef(yk)=nef(yk)+1 !+1 pr nu
      nea(yk)=nea(yk)+2 !+2 pr rate et nu
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
  nasso_s_fp =0
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
  
  
  npmtot = nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sum(nef)+sum(nvc)+ntrtot+nvc_err
  
  
  ! base de splines transfos
  if (any(idlink.eq.2)) then 
     call design_splines(ier)
     if (ier.eq.-1) then
        dens_res=-1.d9
        go to 1589
     end if
  end if
  
  
  !!!!!!!!!!
  !print*,"deb pred"

  ! prediction des EAs
  dens_res = pred(re1,neatot0,ind0,nmescur0)

  !print*,"fin pred"
  
  !print*,"dens_res",dens_res
  !!!!!!!!!!


  1589 continue
  
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
  
  deallocate(bp)
  
  deallocate(fix,refix)

  return
  
end subroutine dens



!----------------------------------------------!
!     INDIVIDUAL RANDOM EFFECTS PREDICTION     !
!----------------------------------------------!

double precision function pred(re,m,i,ncur)

  ! re = vecteur des effets aleatoires individuels
  ! m = nb re
  ! i = individu
  ! ncur = nb observations des individus precedents
  
  use modltsm

  implicit none

  integer::m,i,ncur
  double precision,dimension(m)::re
  
  integer ::j,k,l,jj,npm,ll,lll,ii,numSPL,p
  integer ::ier,kk,j1,j2,sumMesYk,yk,ke,sumnrisq
  integer :: sumRe,sumasso
  integer::sumPrmYk
  integer::nxevtcurr
  
  double precision,dimension(maxmes,nv) :: X0_x_yk
  double precision,dimension(maxmes,2*8) :: X0_t,Z_t
  double precision,dimension(2*8) :: re_Z_t
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
  
  double precision,dimension(npmtot)::b1
  double precision,dimension(nxevt)::Xevt,bevt
  double precision,dimension(maxval(nprisq))::brisq
  
  double precision :: eps,som,eta0
  double precision ::Y4,jacobien
  double precision,dimension(maxmes) :: Y1,Y2,Y3
  double precision,dimension(-1:maxval(ntr)-3)::splaa
  double precision,dimension(nasso)::re_asso
  double precision,dimension(nbevt)::risq,surv,surv0
  double precision::SX,x22,div,dens_Y,dens_surv,varexpsurv
  double precision::surv0_glob,surv_glob,fevt,easurv
  
  double precision :: shiftsurv
  
  double precision,dimension(neatot,neatot) :: var_re, inv_var_re
  integer,dimension(neatot)::ind_var_re_0
  double precision,dimension(neatot*(neatot+1)/2) :: var_re_vect
  double precision :: det_var_re
  double precision ::Y4_re
  double precision,dimension(neatot) :: mu_re,Y2_re,Y3_re
  double precision::div_re,dens_re
  
  pred=0.d0
  
  nmescur = ncur
  !print*,'nmescur',nmescur
  
  b1=0.d0
  do k=1,npmtot
     b1(k)=bp(k)
  end do
  !print*, "b1", b1
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! creation de Y1 = vecteur des reponses observees
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
           pred=-1.d9
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
  !print*, "nmes(i,:)", nmes(i,:)
  !print*,'Y1',Y1
  
  
  
  ! creation des elements specifiques aux trajectoires des marqueurs   
  
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
       
      ! mu_r = taux de progression/ rate = beta_r * X_r
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
  end do!! Xtime
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  dens_Y=1.d0
  dens_surv=1.d0
     
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! vecteur des EAs
     
     !print*,"re",re
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! partie longitudinale ! f(Y|re)
     
     
     !print*,"partie longitudinale"
     
     
     !!!!!!!!!!!!!!!!
     !!!!! time-shift
     
     !! s
     shift = mu_s + re(1)
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
     re_Z_t=0.d0
     
     sumMesYk=0
     sumPrmYk=0
     sumRe=1 !repere ds vect re

     do yk =1,ny
     
       !print*,"yk",yk

       if(nmes(i,yk).gt.0) then
          
         if(idmodel(yk).eq.0) then !logistic
         
           ! rate
           rate=0.d0
           rate=mu_r(yk)+re(sumRe+1)
         
           ! nu = location prm
           nu=0.d0
           nu = mu_v(yk) + re(sumRe+2)
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
             mu(j) = mu(j) + re(sumRe+3)
           end do
           !print*, "mu+ea", mu
           
           
           !! var_Y = variance de Y|s,v,u
           ! avec inv_var_Y = matrice inverse de var_Y, comme var_Y diago, alors inverse = coefficient inverse
           ! avec det_var_Y = determinant de var_Y; comme var_Y diago, alors determinant = produit des coefficients diagonaux
           var_Y=0.d0
           inv_var_Y=0.d0
           det_var_Y=1
           do j=1,nmes(i,yk)
             var_Y(j,j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+nvc(yk)+1)**2 !std erreur de mesure
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
                   X0_t(j,ll)=log(time_shifted(sumMesYk+j)) ! log or ln ? here
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
                   Z_t(j,ll)=log(time_shifted(sumMesYk+j)) ! log or ln ? here
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
           re_Z_t=0.d0
           ll=1
           if(npm_ea_t(yk).ne.0) then
             do k=1,npm_ea_t(yk)
               re_Z_t(k)=re(sumRe+idea(yk)+k)
             end do
           end if
           !print*,'re_Z_t',re_Z_t
           
           ! -> mu = esperance conditionnelle
           if(idea(yk).eq.1) mu = mu+re(sumRe+1)
           !print*,"mu+intercept",mu
           if(npm_ea_t(yk).ne.0) mu = mu+matmul(Z_t,re_Z_t)
           !print*,"mu+tous_ea",mu

           
           !! variance de Y|s,re
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
             pred=-1.d9
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
         !print*,"Y2",Y2
              
         Y3=matmul(inv_var_Y,Y2)
         !print*,"Y3",Y3
              
         Y4=DOT_PRODUCT(Y2,Y3)
         !print*,"Y4",Y4
              
         div = (dble(2*3.14159265)**(dble(nmes(i,yk))/2))*sqrt(det_var_Y)
         !print*,"div",div
          
         dens_Y = dens_Y * exp(-Y4/2.d0)/div  
         !print*,"density",exp(-Y4/2.d0)/div
         !print*,"dens_Y",dens_Y
          
       end if
       
       !incrementation
       sumMesYk = sumMesYk+nmes(i,yk)
       sumPrmYk=sumPrmYk+nef(yk)+nvc(yk)+ntr(yk)+1
       sumRe=sumRe+nea(yk)
       
     end do ! fin boucle yk
     !print*,"dens_Y all mqs",dens_Y
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! partie survie ! f(T|re)
  
  
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
           
             call fct_risq_GK(i,shift,ke,brisq,sumnrisq,b1,re,risq,surv,surv0)
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
        m=0
        sumnrisq=0
        do ke=1,nbevt
        
           ! calculer Xevt * bevt
           varexpsurv=0.d0
           if (nxevt.ne.0) then   !si varexpsurv
              varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevt))&
                   ,bevt((nxevtcurr+1):(nxevtcurr+nxevt)))
           end if
           
           
           !! association : shift
           ! partie sans interaction avec temps : intercept (si il y a)
           
           shiftsurv = 0.d0
           
           !intercept
           if(nasso_s.eq.1) shiftsurv = b1(sumnrisq+nprisq(ke)+nvarbyevt+1)
           shiftsurv = shiftsurv * shift
           !print*,"shift",shift
           !print*,"shiftsurv",shiftsurv
           
           
           !! association : effets aleatoires partages
           easurv=0.d0
           if(any(idst(ke,:).eq.0)) then
           
             sumasso = nasso_s + nasso_s_fp
             sumRe=1
             do yk=1,ny
             
               if(idst(ke,yk).eq.0) then !RE
               
                 re_asso=0.d0
                 
                 if(idmodel(yk).eq.0) then !logistic
                   !re_asso(1)= mu_r(yk) + re(sumRe+1)
                   !re_asso(2)= mu_v(yk) + re(sumRe+2)
                   re_asso(1)= re(sumRe+1)
                   re_asso(2)= re(sumRe+2)
                   re_asso(3)= re(sumRe+3)
                 end if 
                 
                 if(idmodel(yk).eq.1) then !linear
                   re_asso(1:nea(yk))= re(sumRe+1:nea(yk))
                 end if 
                 
                 easurv = easurv +DOT_PRODUCT(re_asso,b1(sumnrisq+nprisq(ke)+nvarbyevt+sumasso+1:nea(yk)))
                 
               end if  
                 
               !incrementation
               sumRe=sumRe+nea(yk)
               if(idst(ke,yk).eq.0) sumasso = sumasso + nea(yk)
               if(idst(ke,yk).eq.1) sumasso = sumasso + 1
             
             end do
           
           end if
           !print*,"easurv",easurv
           
            
           ! avoir evt au temps Ti si Devt=1
           if (Devt(i).eq.ke) then     !si sujet i a evt ke
              fevt=risq(ke)*exp(varexpsurv+shiftsurv+easurv)   !fct de risq         
           end if
           
           ! risque cumule jusque Ti
           Surv_glob=surv_glob + &
                  exp(varexpsurv+shiftsurv+easurv)*surv(ke)     
           
           ! troncature : risque cumule au temps T0
           if (idtrunc.eq.1) then
            surv0_glob=surv0_glob+surv0(ke)*exp(varexpsurv+shiftsurv+easurv)    
           end if
           
           nxevtcurr=nxevtcurr+nxevt
           sumnrisq = sumnrisq + nprisq(ke) + nvarbyevt + nassoparevt(ke)
        end do

        ! density de la partie survie
        dens_surv = exp(-Surv_glob)
        
        ! print*,"dens_surv ok"        
        if(Devt(i).gt.0) dens_surv = dens_surv * fevt

        if (idtrunc.eq.1) then
           dens_surv = dens_surv / exp(-surv0_glob) !delayed entry
        end if
        
     else !  pas de survie

        dens_surv = 1.d0

     end if
     
     !print*,"dens_surv",dens_surv
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! partie effets aleatoires ! f(re)
  
     
     !! var_re = variance de s,v,u
     ! avec inv_var_re = matrice inverse de var_re
     ! avec det_var_re = determinant de var_re
     
     l = 0
     ll = 0
     sumPrmYk = 0
     var_re = 0.d0 !diagonal matrix
     var_re(1,1) = b1(nrisqtot+nvarxevt+nasso+nef_s+1) !shift std
     do yk=1,ny
       if(idmodel(yk).eq.0) then !logistic
         do j=1,nea(yk)
           var_re(1+l+j,1+l+j)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+j)
         end do
       end if
       if(idmodel(yk).eq.1) then !linear
         if(idlink(yk).ne.1) then !linear or spl: identifiability contraint
           ll=0
           var_re(1+l+1,1+l+1)=1
         end if
         if(idlink(yk).eq.0) then !identity
           ll=1
           var_re(1+l+1,1+l+1)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+1)
         end if
         if(nea(yk).gt.1) then
           do j=2,nea(yk)
              do k=1,j
                 var_re(1+l+j,1+l+k)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+sumPrmYk+nef(yk)+ll+k-1+j*(j-1)/2)
              end do
           end do
         end if
       end if
       l = l+nea(yk)
       sumPrmYk = sumPrmYk + nef(yk) + nvc(yk) + ntr(yk) + 1
     end do
     
     !0->1 sur diagonale : revient au mm pr calculer inverse et determinant
     ind_var_re_0=0
     do l=1,neatot
      if(var_re(l,l).eq.0) then
        var_re(l,l)=1
        ind_var_re_0(l)=1
      end if
     end do
     !print*,"ind_var_re_0",ind_var_re_0
     
     ! var_re en vecteur
     jj=0
     var_re_vect=0.d0
     do j1=1,neatot
       do j2=j1,neatot
         jj=j1+j2*(j2-1)/2
         var_re_vect(jj)=var_re(j1,j2)
       end do
     end do
     !print*,"var_re_vect",var_re_vect
              
     ! inversion
     CALL dsinv(var_re_vect,neatot,eps,ier,det_var_re)
     if (ier.eq.-1) then
       dens_re=-1.d9
       goto 654
     end if
           
     ! determinant en exponentiel
     det_var_re = exp(det_var_re)
              
     ! retransformation du vecteur var_re_vect en matrice :
     inv_var_re=0.d0
     do j1=1,neatot
       do j2=1,neatot
         if (j2.ge.j1) then
           inv_var_re(j1,j2)=var_re_vect(j1+j2*(j2-1)/2)
         else
           inv_var_re(j1,j2)=var_re_vect(j2+j1*(j1-1)/2)
         end if
       end do
     end do
     !print*,"var_re",var_re
     !print*,"inv_var_re",inv_var_re
     !print*,"det_var_re",det_var_re
     
     !1->0 sur diagonale
     do l=1,neatot
      if(ind_var_re_0(l).eq.1) then
        inv_var_re(l,l)=0
      end if
     end do
     !print*,"inv_var_re",inv_var_re
     
     
     !! mu_re = moyenne de re [s,v,u]
     mu_re=0.d0
     !! moyennes = 0 car moyennes ajoutees au fur et a mesure precedemment
     
     
     ! calcul de la densite
     dens_re=1.d0
     Y2_re=0.d0
     Y3_re=0.d0
     Y4_re=0.d0
     do j=1,neatot
       Y2_re(j) = re(j) - mu_re(j)
     end do
     !print*,"Y2_re",Y2_re
              
     Y3_re=matmul(inv_var_re,Y2_re)
     !print*,"Y3_re",Y3_re
              
     Y4_re=DOT_PRODUCT(Y2_re,Y3_re)
     !print*,"Y4_re",Y4_re
              
     div_re = (dble(2*3.14159265)**(dble(neatot)/2))*sqrt(det_var_re)
     !print*,"div_re",div_re
        
     dens_re = dens_re * exp(-Y4_re/2.d0)/div_re 
     !print*,'density_re',exp(-Y4_re/2.d0)/div_re
     !print*,"dens_re",dens_re
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
  pred = dens_Y * dens_surv * dens_re !density
  pred = log(pred) + jacobien !density in log
     
  !print*, "jacobien", jacobien   
  !print*, "pred", pred
  
  if (pred.eq.-1.d9 .or. pred/pred.ne.1) then 
    pred = -1.d9
    !goto 541
  end if
  
  654 continue

  return
  
end function pred  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      SPLINE BASES      !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------!
! I-SPLINES for link functions !
!------------------------------!

! function design_splines in loglik.f90


!---------------------------------------!
! M-SPLINES for baseline risk functions !
!     in shared random effects case     !
!---------------------------------------!

! functions splines_i and splines_GK_i in loglik.f90



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      SURVIVAL ELEMENTS      !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------!
!     baseline and cumulative RISK FUNCTIONS    !
!         in shared random effects case         !
!-----------------------------------------------!

! function fct_risq in loglik.f90


!-----------------------------------------------!
!     baseline and cumulative RISK FUNCTIONS    !
!    including time-varying linear predictor    !
!          in shared current level case         !
!-----------------------------------------------!

! function fct_risq_GK in loglik.f90


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
