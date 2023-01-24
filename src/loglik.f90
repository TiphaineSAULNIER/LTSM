
module modirtsre

  implicit none
  integer,save ::ny,ns,nv,nea &
       ,maxmes,nobs &
       ,npmtot &
       ,nMC,methInteg,nmescur &
       ,nvarxevt,nbevt,nrisqtot,nasso
  integer,save::nef_s,nvc_s,nef_B,nmu_v,nvc_v,nvc_u,nvc_err     
  double precision,dimension(:),allocatable,save::Y,minY,maxY
  double precision,dimension(:,:),allocatable,save ::X
  integer,dimension(:),allocatable,save::idg
  integer,dimension(:,:),allocatable,save::nmes
  double precision,dimension(:),allocatable,save :: seqMC
  integer,dimension(:),allocatable,save::fix
  double precision,dimension(:),allocatable,save::bfix
  
end module modirtsre




subroutine loglik(Y0,X0 &
     ,idg0 &
     ,ny0,ns0,nv0,nobs0,nmes0 &
     ,npm0,b0,nfix0,bfix0,zitr0 &
     ,fix0,methInteg0,nMC0,dimMC0,seqMC0 &
     ,loglik_res)

  use modirtsre

  IMPLICIT NONE
  
  !Declaration des variables en entree
  integer,intent(in)::nv0,ny0,nMC0,methInteg0,dimMC0,nfix0
  integer, intent(in)::ns0,nobs0,npm0
  !integer,intent(in)::idtrunc0,logspecif0,nbevt0
  !double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
  !integer, dimension(ns0),intent(in)::Devt0
  !integer, dimension(nv0),intent(in)::idsurv0
  !integer,dimension(nbevt0),intent(in)::typrisq0,nz0
  !double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0    
  double precision,dimension(2,ny0),intent(in)::zitr0
  integer, dimension(nv0),intent(in)::idg0
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
  integer::jtemp,i,j,k,ktemp,yk,mi,nbfix
  integer::npmtot0 !ke
  double precision,external::vrais
  
  
  !print*,"b0",b0
  
  !print*,"#"
  !print*,"##"
  !print*,"### IN"
  !print*,"##"
  !print*,"#"
  
  !print*,"ELEMENTS EN ENTREE FORTRAN"
  !print*,"Y0",Y0
  !print*,"X0",X0
  !print*,"Tentr0",Tentr0
  !print*,"Tevt0",Tevt0
  !print*,"Devt0",Devt0
  !print*,"idg0",idg0
  !print*,"idsurv0",idsurv0
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
  !print*,"zitr0",zitr0
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
  
  !rangeY
  allocate(minY(ny0),maxY(ny0))
  do k=1,ny0
     minY(k)=zitr0(1,k)
     maxY(k)=zitr0(2,k)
  end do
  
  methInteg = methInteg0
  nMC = nMC0
  
  allocate(Y(nobs0),X(nobs0,nv0) &
       ,idg(nv0),nmes(ns0,ny0))
       
  nbevt=0
  
  ny=ny0
  ns=ns0
  nv=nv0
  nobs=nobs0
  
  nmes=0
  Y=0.d0
  X=0.d0
  idg=0
  ktemp=0
  
  do k=1,nv
     idg(k)=idg0(k)
     
     jtemp=0
     do i=1,ns
        do yk=1,ny            
           if (k.eq.1) then
              nmes(i,yk)=nmes0(i,yk)   !dim(nmes)=ns*ny    
              do j=1,nmes(i,yk)
                 jtemp=jtemp+1
                 Y(jtemp)=Y0(jtemp)
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
        
  ! prm fixes
  npmtot0 = npm0+nfix0
  allocate(fix(npmtot0))
  fix=0
  fix(1:npmtot0)=fix0(1:npmtot0)
  nbfix=sum(fix)
  if(nbfix.eq.0) then
     allocate(bfix(1))
  else
     allocate(bfix(nbfix))
  end if
  bfix(1:nbfix)=bfix0(1:nbfix)
  
  
  ! creation des parametres
  
  nrisqtot=0
  
  ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
  
  nvarxevt = 0 !+ nvdepsurv

  nef_s=0
  do k=1,nv
     if (idg(k).eq.1) then
        nef_s = nef_s + 1
     end if
  end do
  
  nasso = nbevt*(1+2*ny)
  
  nvc_s = 1
  nef_B = ny
  nmu_v = ny
  nvc_v = ny
  nvc_u = ny 
  nvc_err = ny 
  
  nea = nvc_s + nvc_v + nvc_u ! nb RE
  
  npmtot = nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+nvc_v+nvc_u+nvc_err
  
  ! points qmc
  if(methInteg.ne.3) then 
     allocate(seqMC(1))
  else
     allocate(seqMC(dimMC0*nMC))
     seqMC = seqMC0(1:dimMC0*nMC) 
  end if


  !print*,"deb vrais"

  ! calcul de la vraisemblance
  loglik_res = vrais(b0,npm0)

  !print*,"fin vrais"

  !1589 continue

  deallocate(Y,X,idg,nmes)

  deallocate(minY,maxY)

  deallocate(fix,bfix,seqMC)
  
  !print*,"#"
  !print*,"##"
  !print*,"### OUT"
  !print*,"##"
  !print*,"#"
  
  return
  
end subroutine loglik



!-------------------------------------!
!     LOG-LIKELIHOOD all subjects     !
!-------------------------------------!

double precision function vrais(b,m)

  use modirtsre,only:ns,nmes,nmescur

  implicit none

  integer::m,i
  double precision::vrais_i,temp
  double precision,dimension(m)::b
  
  nmescur=0
  vrais=0.d0
  do i=1,ns
  
     !print*,"## new subject ",i 
  
     temp = vrais_i(b,m,i) 
     
     !print*,"temp ",temp
     
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

  use modirtsre
!  use optim

  IMPLICIT NONE
  integer ::i,j,k,l,m,npm,ll
  integer ::sumMesYk,yk    !,ke  !,sumnrisq
  double precision,dimension(nef_s) :: b_s
  double precision,dimension(maxmes,nef_s) :: X_s
  double precision,dimension(maxmes) :: mu_s
  double precision,dimension(maxmes) :: Xtime
  double precision,dimension(maxmes) :: mu_v
  double precision,dimension(maxmes) :: rate
  double precision,dimension(maxmes) :: asympt_inf,asympt_sup
  double precision,dimension(maxmes,maxmes) :: var_Y, inv_var_Y
  double precision :: det_var_Y
  
  double precision,dimension(maxmes) :: re_z
  double precision,dimension(maxmes) :: re_v,re_u
  
  double precision,dimension(maxmes) :: shift, time_shifted
  double precision,dimension(maxmes) :: nu
  double precision,dimension(maxmes) :: denom,mu,mu_Y
  
  double precision,dimension(nea,nea) ::Ut
  double precision,dimension(npm) :: b
  double precision,dimension(npmtot)::b1
  !double precision,dimension(nxevt)::Xevt,bevt
  !double precision,dimension(maxval(nprisq))::brisq
  !double precision::basso 
  
  double precision :: som
  double precision ::Y4
  double precision,dimension(maxmes) :: Y1,Y2,Y3  !mu
  double precision,dimension(nea)::ui,usim
  !double precision,dimension(nbevt)::risq,surv,surv0
  double precision::SX,x22,div,vrais_Y  !,vrais_surv,varexpsurv
  !double precision::surv0_glob,surv_glob,fevt,easurv  
  !double precision::som_T0,som_Ti
  
  ! definir le nombre total de mesures pour le sujet i : nmestot (valable que pour cette fonction)

  ! if (verbose==1) write(*,*)'i',i 
  b1=0.d0
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

  Ut = 0.d0 !diagonal matrix
  Ut(1,1) = b1(nrisqtot+nvarxevt+nasso+nef_s+1)
  do j=1,nvc_v
    Ut(1+j,1+j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+j)
  end do
  do j=1,nvc_u
    Ut(1+nvc_v+j,1+nvc_v+j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+nvc_v+j)
  end do

  vrais_Y=0.d0
  

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! creation de Y1 = vecteur des reponses observees
  Y1=0.d0
  sumMesYk = 0
  do yk=1,ny
      do j=1,nmes(i,yk)
          Y1(sumMesYk+j)=Y(nmescur+sumMesYk+j)
      end do
      sumMesYk=sumMesYk+nmes(i,yk)
  end do !fin boucle yk
  
  !print*, "nmes(i,:)", nmes(i,:)
  !print*, "Y1", Y1
  
  !! mu_s
  b_s=0.d0
  X_s=0.d0
  l=0
  do k=1,nv
     if (idg(k).ne.0) then
        l=l+1
        do j=1,sum(nmes(i,:))
           X_s(j,l)=dble(X(nmescur+j,k))  
        end do
        b_s(l)=b1(nrisqtot+nvarxevt+nasso+l)
     end if
  end do
  mu_s = matmul(X_s,b_s)
  
  !print*,"b_s",b_s
  !print*,"X_s",X_s
  !print*,"mu_s",mu_s
  
  !! Xtime
  Xtime=0.d0
  do j=1,sum(nmes(i,:))
      Xtime(j)=dble(X(nmescur+j,1))  
  end do
  
  !! mu_v
  mu_v=0.d0
  l=0
  do k=1,ny
     do j=1,nmes(i,k)
        l=l+1
        mu_v(l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+k)
     end do
  end do
  
  !! rate
  rate=0.d0
  l=0
  do k=1,ny
     do j=1,nmes(i,k)
        l=l+1
        rate(l)=b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+k)
     end do
  end do
  
  !! asymptotes
  asympt_inf=0.d0
  asympt_sup=0.d0
  l=0
  do k=1,ny
     do j=1,nmes(i,k)
        l=l+1
        asympt_inf(l)=minY(k)
        asympt_sup(l)=maxY(k)
     end do
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 ! contribution individuelle a la vraisemblance
  vrais_i=0.d0
  
  som=0.d0
  !som_T0=0.d0 !delayed entry
  
  !print*, "deb MCMC"
  
  do l=1,nMC

     !print*, "nMC", l

     vrais_Y=1.d0
     
     
     !!!!!!!!! MC pour EA !!!!!!!!!

     if(methInteg.eq.1) then 
        ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!

        ! simuler les effets aleatoires
        x22=0.d0
        SX=1.d0
        do j=1,nea
          call bgos(SX,0,usim(j),x22,0.d0)
          !print*,"usim=",usim(j)
        end do
        ui=0.d0
        ui=matmul(Ut,usim)
        !print*,"usim=",usim(j)," ui=",ui, " Ut=",Ut, "  nea=",nea

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
           do j=1,nea
            call bgos(SX,0,usim(j),x22,0.d0)
           end do
           ui=0.d0
           ui=matmul(Ut,usim)
        end if

     else 
        ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!

        ! simuler les effets aleatoires
        usim=0.d0
        do j=1,nea
          usim(j)=seqMC(nMC*(j-1)+l)
        end do
        ui=0.d0
        ui=matmul(Ut,usim)

     end if ! fin if methInteg
     
     !print*,"Ut",Ut
     !print*,"usim",usim
     !print*,"ui", ui

     
     !!! decomposition ui
     
     ! zi
     re_z=0.d0
     do j=1,sum(nmes(i,:))
        re_z(j)=ui(1)  
     end do
     
     ! vik and uik
     ll=0
     re_v=0.d0
     re_u=0.d0
     do k=1,ny
        do j=1,nmes(i,k)
            ll=ll+1
            re_v(ll)=ui(1+k)
            re_u(ll)=ui(1+ny+k)
        end do
     end do
     
     !print*,"re_z",re_z
     !print*,"re_v",re_v
     !print*,"re_u",re_u
     
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!! partie longitudinale
     
     
     !!! time-shift
     
     !! s
     shift = mu_s + re_z
     !! t+s
     time_shifted = Xtime + shift
     
     !print*, "mu_s", mu_s
     !print*, "shift", shift
     !print*, "Xtime", Xtime
     !print*, "time_shifted", time_shifted
     
     !! v
     nu = mu_v + re_v
     
     !print*, "mu_v", mu_v
     !print*, "nu", nu
     
     !!! mu
     
     !! denominator
     denom = 0.d0
     do j=1,sum(nmes(i,:))
        denom(j) = 1 + exp( - rate(j)*time_shifted(j) )
        denom(j) = denom(j) ** nu(j)
     end do
     
     !print*, "rate", rate
     !print*, "denom", denom
     
     !! mu
     mu = 0.d0
     do j=1,sum(nmes(i,:))
        mu(j) = asympt_inf(j) + (asympt_sup(j) - asympt_inf(j)) / denom(j)
     end do
     
     !print*, "asympt_inf", asympt_inf
     !print*, "asympt_sup", asympt_sup
     !print*, "mu", mu
     
     !! mu_Y
     mu_Y = mu + re_u
     
     !print*, "mu_Y", mu_Y
     
     
     sumMesYk=0

     do yk =1,ny

          if(nmes(i,yk).gt.0) then
          
              !print*, "yk", yk
          
              !! var_Y = variance de Y|s,v,u
              ! avec inv_var_Y = matrice inverse de var_Y, comme var_Y diago, alors inverse = coefficient inverse
              ! avec det_var_Y = determinant de var_Y, comme var_Y diago, alors determinant = produit des coefficients diagonaux
              var_Y=0.d0
              inv_var_Y=0.d0
              det_var_Y=1
              do j=1,nmes(i,yk)
                  var_Y(j,j) = b1(nrisqtot+nvarxevt+nasso+nef_s+nvc_s+nef_B+nmu_v+nvc_v+nvc_u+yk)**2
                  !print*,"var_Y(j,j)",var_Y(j,j)
                  inv_var_Y(j,j) = 1/var_Y(j,j)
                  det_var_Y = det_var_Y * var_Y(j,j)
              end do
              
              ! calcul de la vrais
              Y2=0.d0
              Y3=0.d0
              Y4=0.d0
              do j=1,nmes(i,yk)
                 Y2(j) = Y1(sumMesYk+j) - mu_Y(sumMesYk+j)
              end do
              
              !print*,"Y2",Y2
              
              Y3=matmul(inv_var_Y,Y2)
              
              !print*,"Y3",Y3
              
              Y4=DOT_PRODUCT(Y2,Y3)
              
              !print*,"Y4",Y4
              
              div = (dble(2*3.14159265)**(dble(nmes(i,yk))/2))*sqrt(det_var_Y) !exp(det_var_Y)
              
              !print*,"div",div
              
              !print*,"density",exp(-Y4/2.d0)/div
              
              vrais_Y = vrais_Y * exp(-Y4/2.d0)/div
              
          end if
        
          sumMesYk = sumMesYk+nmes(i,yk)
     end do ! fin boucle yk
     !if(i.lt.4) print*,"avant survie, vrais_Y",vrais_Y
     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!! partie survie deleted
     
  
     som = som + vrais_Y

     
  end do ! fin boucle nMC
  
  !print*, "fin MCMC"
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! contribution individuelle a vraisemblance
  
  vrais_i = vrais_i + log(som) - log(dble(nMC))
  
  !print*, "som", som
  !print*, "nMC", nMC
  !print*, "vrais_i", vrais_i
  
  !654 continue

  return

end function vrais_i



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!